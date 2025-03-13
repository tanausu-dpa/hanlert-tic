      !> Compute SVD solution
      module svd_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC/NSSCC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     27/02/2023
!  Last version:
!     20/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     20/12/2024:    V4.0.0 - Revised headers (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!#####################################################################
!#####################################################################
!
!  To do:
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  SVD_Solve
!    Solve the system of linear equations for the step in the
!  Levenberg-Marquardt inversion with SVD decomposition
!
!  Check_W
!    Iterate the singular values using the SVD matrix
!
!  Svbksb
!    Calculate the solution of a set of linear equations after SVD
!  decomposition
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use commons_mod
      use parameters_mod, only: TINYSVDS
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the system of linear equations for the step in the
      !! Levenberg-Marquardt inversion with SVD decomposition\n
      !!          A(double(:,:)): System matrix\n
      !!            B(double(:)): Independent term\n
      !!            X(double(:)): Solution\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!       SVD_type(integer): Type of SVD solution
      subroutine SVD_Solve(A,B,X,N,Inf_Nodes,SVD_type)

      ! I/O

      type(Nodes_class), intent(inout):: Inf_Nodes
      integer, intent(in):: N,SVD_type
      double precision, dimension(:), intent(in):: B
      double precision, dimension(:), intent(inout):: X
      double precision, dimension(:,:), intent(in):: A

      ! Local

      character(3):: length
      character(30):: fmt

      integer, parameter:: LWMAX = 2000
      integer:: i,n_in,INFO,LWORK

      double precision:: WMIN,WMAX,TMP_Threshold_Svd
      double precision, dimension(LWMAX):: WORK
      double precision, dimension(N):: W,W0
      double precision, dimension(N,N):: U,V,VT


      ! Master
      if (pid.eq.0) then

        ! Verbose
        umsg = ' - SVD solution'
        call verboseI(3)

        ! Prepare format
        n_in = n
        if (n_in*15.gt.471) then
          do while(n_in*15.gt.467)
            n_in = n_in - 1
          end do
          write(length, "(i3)") n_in
          fmt = '(A,'//trim(adjustl(length))//'es15.4," ...")'
          fmt = trim(adjustl(fmt))
        else
          write(length, "(i3)") n
          fmt = '(A,'//trim(adjustl(length))//'es15.4)'
          fmt = trim(adjustl(fmt))
        end if

        ! Initialize DGESV
        INFO = 0

        ! Get optimal LWORK
        LWORK = -1
        call DGESVD('A','A',N,N,A,N,W0,U,N,VT,N,WORK,LWORK,INFO)
        LWORK = min(LWMAX, int(WORK(1)))

        ! Get SVD
        call DGESVD('A','A',N,N,A,N,W0,U,N,VT,N,WORK,LWORK,INFO)

        ! Error
        if (INFO.gt.0) goto 999

        ! Get transpose
        V = transpose(VT)

        ! Verbose
        write(umsg, fmt=fmt) "   W(original) = ",(W0(i),i=1,N_in)
        call verboseI(3)

        ! Get current threshold
        TMP_Threshold_Svd = Inf_Nodes%Threshold_Svd

        ! Loop until done
        do while (.True.)

          ! copy W
          W = W0

          ! Type of SVD algorithm
          select case(SVD_type)

            ! Traditional
            case(0)

              ! Set limits from singular values of A
              WMAX = W(1)
              WMIN = WMAX*Inf_Nodes%Threshold_Svd

              ! For each singular value
              do i=1,N

                ! If singular value below threshold, do 0
                if (W(i).lt.WMIN) W(i) = 0d0

              end do ! singular values

            ! SIR like
            case(2)

              ! Check the singular values
              call Check_W(Inf_Nodes, N, V, W)

          end select

          ! Solve the system of equations
          call Svbksb(U,W,V,N,N,B,X)

          ! If step value within bounds or too large threshold
          if ((maxval(X).lt.Inf_Nodes%Max_Step.and. &
              minval(X).gt.-Inf_Nodes%Max_Step).or. &
              Inf_Nodes%Threshold_Svd.gt.1d-3) then

            ! Break from the loop
            exit

          ! Too large of a step with room to increase threshold
          else

            ! Increase threshold by five
            Inf_Nodes%Threshold_Svd = Inf_Nodes%Threshold_Svd*5d0

          end if ! Step value within bounds or too large threshold

        end do ! While loop

 

        ! Verbose
        write(umsg,'(A,es15.4)') '   Threshold = ', &
                                 Inf_Nodes%Threshold_Svd
        call verboseI(3)
        write(umsg, FMT=fmt) '   W (modified) = ',(W(i),i=1,N_in)
        call verboseI(3)
        write(umsg, FMT=fmt) '   SVD solution = ',(X(i),i=1,N_in)
        call verboseI(3)

        ! Recover threshold
        Inf_Nodes%Threshold_Svd = TMP_Threshold_Svd

      end if ! Verbose

      ! Share info
999   call MPI_BCAST(INFO,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)

      ! Error
      if (INFO.gt.0) then

        ! Issue error
        umsg = 'The SVD algorithm failed to converge'
        urou = 'SVD_Solve'
        call aborted
        return

      end if ! Error

      ! Share solution
      call MPI_BCAST(X(1),N,MPI_DOUBLE_PRECISION,0,MPI_COMM_RT,ierr)

      return

      end subroutine SVD_Solve

!#####################################################################
!#####################################################################
!#####################################################################

      !> Iterate the singular values using the SVD matrix\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!              N(integer): Size of the system of equations\n
      !!          V(double(N,N)): SVD V matrix\n
      !!             W(double(N): SVD singular values
      subroutine Check_W(Inf_Nodes,N,V,W)

      ! I/O

      type(Nodes_class), intent(in):: Inf_Nodes
      integer, intent(in):: N
      double precision, dimension(N), intent(inout):: W
      double precision, dimension(N,N), intent(in):: V

      ! Local

      integer:: ido,j,i0,i1,i2

      double precision:: wtn,wtx
      double precision, dimension(N):: wt,ww
      double precision, dimension(2,N):: wi


      ! Iterate this process twice
      do ido=1,2

        ! Initialize new values to cero
        ww = 0d0

        ! Initialize index shift
        i2 = 0

        ! For the parameter indexes we have to account for
        do i0=Inf_Nodes%Indx_b,Inf_Nodes%Indx_e

          ! If inverting that parameter
          if (Inf_Nodes%Nodes_Flags(i0)) then

            ! Get index limits for this variable
            i1 = i2 + 1
            i2 = i2 + Inf_Nodes%Num_Vary(i0)

            ! Compute from singular values and the SVD V matrix
            wt = sum(V(i1:i2,:)*V(i1:i2,:), dim=1)*W

            ! Get maximum of wt
            wtx = max(maxval(wt), 0d0)

            ! Get minimum from threshold
            wtn = wtx*Inf_Nodes%Threshold_Svd

            ! For each element
            do j=1,N

              ! If larger or equal than lower limit, add to ww
              if (wt(j).ge.wtn) ww(j) = ww(j) + wt(j)

            end do ! Elements

          end if ! Inverting variable

        end do ! Movel parameters

        ! Update iteration
        wi(ido,:) = ww

      end do ! Iterate wtice

      ! For each element
      do j=1,N

        ! If value in second iteration above a hard-coded
        ! threshold
        if (abs(wi(2,j)).gt.1d-10) then

          ! Get new value from both iterations
          w(j) = wi(1,j)*wi(1,j)/wi(2,j)

        ! Otherwise leave it at zero
        else

          w(j) = 0d0

        end if ! Final value above threshold

      end do ! Elements

      return

      end subroutine Check_W

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the solution of a set of linear equations after SVD
      !! decomposition\n
      !!  U(double(M,N)): SVD U matrix\n
      !!    W(double(N)): Singular values\n
      !!   V(double(N,N): SVD V matrix\n
      !!      M(integer): Dimension of independent term\n
      !!      N(integer): Dimen sion of solution\n
      !!    B(double(M)): Independent term\n
      !!    X(double(N)): Solution
      !! Adapted from numerical recipes
      subroutine Svbksb(U,W,V,M,N,B,X)

      ! I/O

      integer, intent(in):: M, N
      double precision, dimension(N), intent(in):: W
      double precision, dimension(M), intent(in):: B
      double precision, dimension(N), intent(inout):: X
      double precision, dimension(M,N), intent(in):: U
      double precision, dimension(N,N), intent(in):: V

      ! Local

      integer:: j,NMAX

      double precision, dimension(:), allocatable:: tmp


      ! Get maximum dimension
      NMAX = max(M,N)

      ! Allocate auxiliar
      allocate(tmp(NMAX))

      ! For each element singular value
      do j=1,N

        ! If singular value above what we consider 0
        if (abs(W(j)).gt.TINYSVDS) then

          ! Add product U*B
          tmp(j) = sum(U(:,j)*B)/W(j)

        ! If too small
        else

          ! Zero
          tmp(j) = 0d0

        end if ! Singular value's value

      end do ! Elements in solution

      ! For each element in solution
      do j=1,N

        ! Get solution
        X(j) = sum(V(j,:)*tmp)

      end do ! Elements in solution

      ! Free memory
      deallocate(tmp)

      return

      end subroutine Svbksb

!#####################################################################
!#####################################################################
!#####################################################################

      end module svd_mod
