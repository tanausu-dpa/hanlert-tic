      !> Runs the NG acceleration algorithm
      module ng_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     06/29/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     09/20/2018:    V1.0.0 - First Version (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  NG:
!    Runs the NG acceleration algorithm
!
!#####################################################################
!#####################################################################
!#####################################################################

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Processes the NG acceleration algorithm.\n
      !>      Y(dfloat(:)): Vector to accelerate in the current
      !>                    iteration.\n
      !>        M(integer): Dimension of the vector to accelerate.\n
      !>        N(integer): Order of the NG acceleration.\n
      !>   YY(dfloat(:,:)): Matrix with memory of the vector to
      !>                    accelerate.\n
      !>     ntry(integer): Step of the current acceleration.\n
      !>     doit(logical): Indicates if acceleration is ready to be
      !>                    applied.\n
      subroutine NG(M,N,YY,ntry,doit)

      ! I/O
      logical, intent(out):: doit
      integer, intent(in):: M,N
      integer, intent(inout):: ntry
      double precision, dimension(:,:), intent(inout):: YY

      ! Local
      integer:: i,j,k
      integer, dimension(N):: indx
      double precision:: DI, DY, WT
      double precision, dimension(M):: Y
      double precision, dimension(N):: C
      double precision, dimension(N,N):: A

      ! Initialize doit
      doit = .False.

      ! Check order is in the correct range
      if (N.lt.1.or.N.gt.5) return

      ! If we do not have enough information to perform the
      ! acceleration, we are done
      if (ntry.le.N+1) return

      ! Initialize system of equations to 0
      A = 0d0
      C = 0d0

      !
      ! Build system of equations
      !

      ! For each component of the vector to accelerate
      do k=1,M

        WT = 1d0/(1d0 + abs(YY(k,ntry)))

        ! For each row of the system
        do i=1,N

          DY = YY(k,ntry-1) - YY(k,ntry)
          DI = WT*(DY + YY(k,ntry-i) - YY(k,ntry-i-1))
          C(i) = C(i) + DI*DY

          ! For each column
          do j=1,N

            A(i,j) = A(i,j) + DI*(DY + YY(k,ntry-j) - YY(k,ntry-j-1))

          enddo ! Columns
        enddo ! Rows
      enddo ! Component of vector

      ! Solve system of equations A*X=C
      call DGESV(N,1,A,N,indx,C,N,i)

      ! Indicate that we are ready to apply NG
      doit = .True.

      Y = YY(:,ntry)

      do i=1,N
        do k=1,M
          YY(k,ntry) = YY(k,ntry) + C(i)*(YY(k,ntry-i) - Y(k))
        enddo
      enddo

      return

      end subroutine NG

!#####################################################################
!#####################################################################
!#####################################################################

      end module ng_mod
