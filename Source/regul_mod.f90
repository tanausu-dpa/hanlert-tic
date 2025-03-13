      !> Manages the regularization in the inversion
      module regul_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC/NSSCC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     24/02/2023
!  Last version:
!     17/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     17/12/2024:    V4.0.0 - Updated headers (TdPA)
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
!  Init_Regul
!    Initialize the dimensionality of the regularization for each
!  parameter
!
!  Get_Regl_all
!    Calculate the regularization matrix
!
!  Compute_Regl
!    Calculate the regularization matrix for one inversion parameter
!  and add to the total
!
!  Firderv_Regul
!    Evaluate regularization to the first derivative
!
!  Secderv_Regul
!    Evaluate regularization to the second derivative
!
!  Mean_Regul
!    Evaluate regularization to the mean
!
!  Const_Regul
!    Evaluate regularization to a constant
!
!  Constl1_Regul
!    Evaluate regularization to a constant, decreasing penalty if
!  beyond one unit
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the dimensionality of the regularization for each
      !! parameter\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data
      subroutine Init_Regul(Inf_Nodes)

      ! I/O

      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local

      integer:: i


      ! For each variable
      do i=1,11

        ! If regularizing the variable
        if (Inf_Nodes%Nodes_Regul(i)) then

          ! Type of regularization
          select case(Inf_Nodes%Indx_regul(i))

            ! Mean regularization
            case(1)

              ! Size is nodes
              Inf_Nodes%Num_regul(i) = Inf_Nodes%Num_Nodes(i)

            ! Constant regularization
            case(2)

              ! Size is nodes
              Inf_Nodes%Num_regul(i) = Inf_Nodes%Num_Nodes(i)

            ! First derivative regularization
            case(3)

              ! Size is nodes minus 1
              Inf_Nodes%Num_regul(i) = Inf_Nodes%Num_Nodes(i)-1

            ! Second derivative regularization
            case(4)

              ! Size is nodes minus 2
              Inf_Nodes%Num_regul(i) = Inf_Nodes%Num_Nodes(i)-2

            ! Constant regularization with less penalty beyond
            ! one unit
            case(5)

              ! Size is nodes
              Inf_Nodes%Num_regul(i) = Inf_Nodes%Num_Nodes(i)

            ! Nothing
            case default

              ! No regularization, zero size
              Inf_Nodes%Num_regul(i) = 0
              cycle

          end select ! Type of regularization

          ! If nodes to regularize
          if (Inf_Nodes%Num_regul(i).gt.0) then

            ! Scale weight to number of nodes
            Inf_Nodes%Regul_weight(i) = Inf_Nodes%Regul_weight(i)/ &
                                        dble(Inf_Nodes%Num_regul(i))

          ! No nodes
          else

            ! Abort
            write (umsg, "(A,i2,A)") &
                'Node index = ', i,'. '// &
                'The number of the nodes and the regularization '// &
                'function do not match (the chosen regularization'// &
                'likely requires more nodes)'
            urou = 'Init_Regul'
            call aborted
            return

          end if ! Regularization nodes

        ! No regularization
        else

          ! Size zero
          Inf_Nodes%Num_regul(i) = 0

        end if ! If regularizing

      end do ! Variables

      return

      end subroutine Init_Regul

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the regularization matrix\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!    Matrix_Flag(logical): If the matrix is to be saved or just
      !!                          the total penalty\n
      !!        Rgl(Regul_class): Structure with regularization data
      subroutine Get_Regl_all(Inf_Nodes,Matrix_Flag,Rgl)

      ! I/O

      type(Nodes_class), intent(in):: Inf_Nodes
      type(Regul_class), intent(inout):: Rgl
      logical, intent(in):: Matrix_Flag

      ! Local

      integer:: i, indx

      double precision:: Tmp_Penalty
      double precision, dimension(:), allocatable:: Resi
      double precision, dimension(:,:), allocatable:: LL


      ! If outputting matrix
      if (Matrix_Flag) then

        ! Initialize matrix to zero
        Rgl%Regul_F = 0d0
        Rgl%Regul_H = 0d0

      end if ! If producing output matrix

      ! Initialize
      indx = 0
      Rgl%Penalty = 0
      Tmp_Penalty = 0

      ! Run over variables being inverted
      do i=Inf_Nodes%Indx_b,Inf_Nodes%Indx_e

        ! If inverting the variable
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If regularizing the variable
          if (Inf_Nodes%Nodes_Regul(i)) then

            ! Type of regularization
            select case(Inf_Nodes%Indx_regul(i))

              ! Mean regularization
              case(1)
                call Mean_Regul(Inf_Nodes%Node(i)%Var, &
                                Inf_Nodes%Node(i)%H, &
                                Inf_Nodes%Num_Nodes(i), LL, Resi)

              ! Constant regularization
              case(2)
                call Const_Regul(Inf_Nodes%Node(i)%Var, &
                                 Inf_Nodes%Node(i)%H, &
                                 Inf_Nodes%Num_Nodes(i), LL, Resi, &
                                 Inf_Nodes%Const(i))

              ! First derivative regularization
              case(3)
                call Firderv_Regul(Inf_Nodes%Node(i)%Var, &
                                   Inf_Nodes%Node(i)%H, &
                                   Inf_Nodes%Num_Nodes(i), LL, Resi)

              ! Second derivative regularization
              case(4)
                call Secderv_Regul(Inf_Nodes%Node(i)%Var, &
                                   Inf_Nodes%Node(i)%H, &
                                   Inf_Nodes%Num_Nodes(i), LL, Resi)

              ! Constant regularization with less penalty beyond
              ! one unit
              case(5)
                call Constl1_Regul(Inf_Nodes%Node(i)%Var, &
                                   Inf_Nodes%Node(i)%H, &
                                   Inf_Nodes%Num_Nodes(i), LL, Resi, &
                                   Inf_Nodes%Const(i))

              ! Not a valid option
              case default

                ! Move index and skip
                indx = indx+Inf_Nodes%Num_Nodes(i)
                cycle

            end select ! Type of regularization

            ! Scale residual
            Resi = Resi/Inf_Nodes%Scal(i)

            ! Compute regularization
            call Compute_Regl(Inf_Nodes, Matrix_Flag, indx, i, &
                              LL, Resi, Rgl%Regul_H, Rgl%Regul_F, &
                              Rgl%Penalty)

            ! Master
            if (pid.eq.0) then

              ! Verbosity
              write(umsg,'(A,i5,3x,A,es15.4)') &
                ' - Regularization model parameter = ',i, &
                'Penalty = ', Rgl%Penalty-Tmp_Penalty
              call verboseI(3)

            end if ! Master

            ! Get a copy of the penalty right now
            Tmp_Penalty = Rgl%Penalty

            ! Advance index
            indx = indx + Inf_Nodes%Num_Vary(i)

          ! Not regularizing
          else

            ! Advance index
            indx = indx + Inf_Nodes%Num_Vary(i)

          end if ! Regularizing variable
        end if ! Inverting variable

      end do ! Variables

      return

      end subroutine Get_Regl_all

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the regularization matrix for one inversion
      !! parameter and add to the total\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!    Matrix_Flag(logical): If the matrix is to be saved or just
      !!                          the total penalty\n
      !!       Indx_Pos(integer): Position in the regularization
      !!                          matrix\n
      !!      Indx_Para(integer): Index of the parameter\n
      !!         LL(double(:,:)): Derivative of regularization
      !!                          function\n
      !!            Resi(double): Evaluation of the regularization
      !!                          function\n
      !!    Regul_H(double(:,:)): Regularization matrix\n
      !!      Regul_F(double(:)): Regularizaiton vector\n
      !!         Penalty(double): Final penalty
      subroutine Compute_Regl(Inf_Nodes,Matrix_Flag,Indx_Pos, &
                              Indx_Para,LL,Resi,Regul_H,Regul_F, &
                              Penalty)

      ! I/O

      type(Nodes_class), intent(in):: Inf_Nodes
      logical, intent(in):: Matrix_Flag
      integer, intent(in):: Indx_Pos, Indx_Para
      double precision, intent(inout) :: Penalty
      double precision, dimension(:), intent(inout):: Regul_F
      double precision, dimension(:), &
                        allocatable, intent(inout):: Resi
      double precision, dimension(:,:), intent(inout):: Regul_H
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: LL

      ! Local

      integer:: k,l,m,n,NN,shift


      ! Get shift for regularization positions
      shift = 1 - Inf_Nodes%Node_Vary(1,Indx_Para)

      ! Shorter variable
      NN = Inf_Nodes%Num_regul(Indx_Para)

      ! Add penalty
      Penalty = Penalty + Inf_Nodes%Regul_weight(Indx_Para)* &
                         sum(Resi(1:NN)*Resi(1:NN))

      ! If computing matrix
      if (Matrix_Flag) then

        ! For each node that changes (row)
        do k=Inf_Nodes%Node_Vary(1,Indx_Para), &
             Inf_Nodes%Node_Vary(2,Indx_Para)

          ! Get true index
          m = k + shift

          ! Contribution to vector
          Regul_F(Indx_Pos+m) = Regul_F(Indx_Pos+m) - &
                                Inf_Nodes%Regul_weight(Indx_Para)* &
                                sum(Resi(1:NN)*LL(1:NN,k))

          ! For each varying node (column)
          do l=Inf_Nodes%Node_Vary(1,Indx_Para), &
               Inf_Nodes%Node_Vary(2,Indx_Para)

            ! Get true index
            n = l + shift

            ! Contribution to matrix
            Regul_H(Indx_Pos+m,Indx_Pos+n) = &
                                  Regul_H(Indx_Pos+m,Indx_Pos+n) + &
                                  Inf_Nodes%Regul_weight(Indx_Para)* &
                                  sum(LL(1:NN,k)*LL(1:NN,l))

          end do ! Varying nodes (column)
        end do ! Varying nodes (row)

      end if ! Computing matrix

      ! Deallocate regularization data
      deallocate(Resi,LL)

      return

      end subroutine Compute_Regl

!#####################################################################
!#####################################################################
!#####################################################################

      !> Evaluate regularization to the first derivative\n
      !!      Var(double(:)): Node values\n
      !!      Tau(double(:)): Node positions\n
      !!  Num_Nodes(integer): Number of nodes
      !!     LL(double(:,:)): Derivative of regularization function\n
      !!     Resi(double(:)): Regularization function
      subroutine Firderv_Regul(Var,Tau,Num_Nodes,LL,Resi)

      ! I/O

      integer, intent(in):: Num_Nodes
      double precision, dimension(:), intent(in):: Var
      double precision, dimension(:), intent(in):: Tau
      double precision, dimension(:), allocatable, intent(out):: Resi
      double precision, dimension(:,:), allocatable, intent(out):: LL

      ! Local

      integer:: i

      double precision:: dtau


      ! Allocate outputs
      allocate(Resi(Num_Nodes-1),LL(Num_Nodes-1,Num_Nodes))

      ! Initialize
      Resi = 0d0
      LL = 0d0

      ! For each node but the last
      do i=1,Num_Nodes-1

        ! Get delta step with next, scaled to max
        dtau = abs(Tau(i+1) - Tau(i))

        ! Get first derivative
        Resi(i) = (Var(i+1) - Var(i))/dtau

        ! Get derivative of the regularization function
        LL(i,i) = -1d0/dtau
        LL(i,i+1) = 1d0/dtau

      end do ! Nodes

      return

      end subroutine Firderv_Regul

!#####################################################################
!#####################################################################
!#####################################################################

      !> Evaluate regularization to the second derivative\n
      !!      Var(double(:)): Node values\n
      !!      Tau(double(:)): Node positions\n
      !!  Num_Nodes(integer): Number of nodes
      !!     LL(double(:,:)): Derivative of regularization function\n
      !!     Resi(double(:)): Regularization function
      subroutine Secderv_Regul(Var,Tau,Num_Nodes,LL,Resi)

      ! I/O

      integer, intent(in):: Num_Nodes
      double precision, dimension(:), intent(in):: Var
      double precision, dimension(:), intent(in):: Tau
      double precision, dimension(:), allocatable, intent(out):: Resi
      double precision, dimension(:,:), allocatable, intent(out):: LL

      ! Local

      integer:: i

      double precision:: dtau0,dtau1,A,B,C,tmp


      ! Allocate outputs
      allocate(Resi(Num_Nodes-2),LL(Num_Nodes-2,Num_Nodes))

      ! Initialize
      Resi = 0d0
      LL = 0d0

      ! For each node but the two last
      do i=1,Num_Nodes-2

        ! Get two delta steps with two next, scaled to max
        dtau0 = abs(Tau(i+1) - Tau(i))
        dtau1 = abs(Tau(i+2) - Tau(i+1))

        ! Get inverse of average
        tmp = 2d0/(dtau0+dtau1)

        ! Get second derivative coefficients
        A = tmp/dtau1
        C = tmp/dtau0
        B = -2d0/(dtau0*dtau1)

        ! Evaluate regularization
        Resi(i) = A*Var(i+2) + B*Var(i+1) + C*Var(i)

        ! Get derivative of the regularization function
        LL(i,i) = C
        LL(i,i+1) = B
        LL(i,i+2) = A

      end do ! Nodes

      return

      end subroutine Secderv_Regul

!#####################################################################
!#####################################################################
!#####################################################################

      !> Evaluate regularization to the mean\n
      !!      Var(double(:)): Node values\n
      !!      Tau(double(:)): Node positions\n
      !!  Num_Nodes(integer): Number of nodes
      !!     LL(double(:,:)): Derivative of regularization function\n
      !!     Resi(double(:)): Regularization function
      subroutine Mean_Regul(Var,Tau,Num_Nodes,LL,Resi)

      ! I/O

      integer, intent(in):: Num_Nodes
      double precision, dimension(:), intent(in):: Var
      double precision, dimension(:), intent(in):: Tau
      double precision, dimension(:), allocatable, intent(out):: Resi
      double precision, dimension(:,:), allocatable, intent(out):: LL

      ! Local

      integer:: i,j

      double precision:: dtau,avg,numi


      ! Allocate outputs
      allocate(Resi(Num_Nodes),LL(Num_Nodes,Num_Nodes))

      ! Initialize
      Resi = 0d0
      LL = 0d0

      ! Inverse number of nodes
      numi = 1d0/dble(Num_Nodes)

      ! Compute average
      avg = sum(Var(1:Num_Nodes))*numi


      !
      ! First point
      !

      ! Forward step
      dtau = abs(Tau(2) - Tau(1))
      dtau = 1d0/dtau

      ! Evaluate residual
      Resi(1) = (Var(1) - avg)*dtau

      ! Diagonal
      LL(1,1) = (1d0 - numi)*dtau

      ! Rest
      do j=2,Num_Nodes

        ! Contribution to derivative
        LL(1,j) = -dtau*numi

      end do ! Nodes

      !
      ! Last point
      !

      ! Backward step
      dtau = abs(Tau(Num_nodes) - Tau(Num_nodes-1))
      dtau = 1d0/dtau

      ! Evaluate residual
      Resi(Num_nodes) = (Var(Num_nodes) - avg)*dtau

      ! Diagonal
      LL(Num_nodes,Num_nodes) = (1d0 - numi)*dtau

      ! Rest
      do j=1,Num_Nodes-1

        ! Contribution to derivative
        LL(Num_nodes,j) = -dtau*numi

      end do ! Nodes

      !
      ! Rest of nodes
      do i=2,Num_Nodes-1

        ! Centered step
        dtau = 0.5d0*abs(Tau(i+1) - Tau(i-1))
        dtau = 1d0/dtau

        ! Evaluate residual
        Resi(i) = (Var(i) - avg)*dtau

        ! For each node
        do j=1,Num_Nodes

          ! If diagonal
          if (j.eq.i) then

            ! Contribution to derivative
            LL(i,j) = (1d0 - numi)*dtau

          ! If non-diagonal
          else

            ! Contribution to derivative
            LL(i,j) = -dtau*numi

          end if ! Diagonal

        end do ! Nodes
      end do ! Nodes

      return

      end subroutine Mean_Regul

!#####################################################################
!#####################################################################
!#####################################################################

      !> Evaluate regularization to a constant\n
      !!       Var(double(:)): Node values\n
      !!       Tau(double(:)): Node positions\n
      !!   Num_Nodes(integer): Number of nodes
      !!      LL(double(:,:)): Derivative of regularization function\n
      !!      Resi(double(:)): Regularization function\n
      !!  Value_Const(double): Constant value to regularize to
      subroutine Const_Regul(Var,Tau,Num_Nodes,LL,Resi,Value_Const)

      ! I/O

      integer, intent(in):: Num_Nodes
      double precision, intent(in):: Value_Const
      double precision, dimension(:), intent(in):: Var
      double precision, dimension(:), intent(in):: Tau
      double precision, dimension(:), allocatable, intent(out):: Resi
      double precision, dimension(:,:), allocatable, intent(out):: LL

      ! Local

      integer:: i

      double precision:: dtau


      ! Allocate outputs
      allocate(Resi(Num_Nodes),LL(Num_Nodes,Num_Nodes))

      ! Initialize
      Resi = 0d0
      LL = 0d0

      ! If only one node
      if (Num_Nodes.eq.1) then

        ! Just the difference and unit matrix
        Resi(1) = (Var(1) - Value_Const)
        LL(1,1) = 1d0

      ! If more than one node
      else

        !
        ! First point
        !

        ! Forward step
        dtau = abs(Tau(2) - Tau(1))
        dtau = 1d0/dtau

        ! Get scaled difference
        Resi(1) = (Var(1) - Value_Const)*dtau

        ! Only diagonal contribution
        LL(1,1) = dtau

        !
        ! Last point
        !

        ! Backward step
        dtau = abs(Tau(Num_Nodes) - Tau(Num_Nodes-1))
        dtau = 1d0/dtau

        ! Get scaled difference
        Resi(Num_Nodes) = (Var(Num_Nodes) - Value_Const)*dtau

        ! Only diagonal contribution
        LL(Num_Nodes,Num_Nodes) = dtau

        !
        ! Rest
        do i=2,Num_Nodes-1

          ! Centered step
          dtau = 0.5d0*abs(Tau(i+1) - Tau(i-1))
          dtau = 1d0/dtau

          ! Get scaled difference
          Resi(i) = (Var(i) - Value_Const)*dtau

          ! Only diagonal contribution
          LL(i,i) = dtau

        end do ! Nodes

      end if ! Number of nodes

      return

      end subroutine Const_Regul

!#####################################################################
!#####################################################################
!#####################################################################

      !> Evaluate regularization to a constant, decreasing penalty if
      !! beyond one unit\n
      !!       Var(double(:)): Node values\n
      !!       Tau(double(:)): Node positions\n
      !!   Num_Nodes(integer): Number of nodes
      !!      LL(double(:,:)): Derivative of regularization function\n
      !!      Resi(double(:)): Regularization function\n
      !!  Value_Const(double): Constant value to regularize to
      subroutine Constl1_Regul(Var,Tau,Num_Nodes,LL,Resi,Value_Const)

      ! I/O

      integer, intent(in):: Num_Nodes
      double precision, intent(in):: Value_Const
      double precision, dimension(:), intent(in):: Var
      double precision, dimension(:), intent(in):: Tau
      double precision, dimension(:), allocatable, intent(out):: Resi
      double precision, dimension(:,:), allocatable, intent(out):: LL

      ! Local

      integer:: i

      double precision:: dtau,tmp,dif


      ! Allocate outputs
      allocate(Resi(Num_Nodes),LL(Num_Nodes,Num_Nodes))

      ! Initialize
      Resi = 0d0
      LL = 0d0

      ! If only one node
      if (Num_Nodes.eq.1) then

        ! Difference
        dif = abs(Var(1) - Value_Const)

        ! If difference larger than 1
        if (dif.gt.1d0) then

          ! Residual
          Resi(1) = sqrt(dif)

          ! Derivative
          if (Var(1).gt.Value_Const) then
            LL(1,1) =  0.5d0/Resi(1)
          else
            LL(1,1) = -0.5d0/Resi(1)
          end if

        ! Difference within a unit
        else

          ! Residual and derivative
          Resi(1) = Var(1) - Value_Const
          LL(1,1) = 1d0

        end if ! Difference larger than unit

      ! Many nodes
      else

        !
        ! First
        !

        ! Forward step
        dtau = abs(Tau(2) - Tau(1))
        dtau = 1d0/dtau

        ! Get difference
        dif = abs(Var(1) - Value_Const)

        ! If difference larger than 1 unit
        if (dif.gt.1d0) then

          ! Square root
          tmp = sqrt(dif)

          ! Residual
          Resi(1) = tmp*dtau

          ! Derivative
          if (Var(1).gt.Value_Const) then
            LL(1,1) =  0.5d0*dtau/tmp
          else
            LL(1,1) = -0.5d0*dtau/tmp
          end if

        ! Difference within a unit
        else

          ! Resiaudl and derivative
          Resi(1) = (Var(1) - Value_Const)*dtau
          LL(1,1) = dtau

        end if ! Difference size

        !
        ! Last
        !

        ! Backward step
        dtau = abs(Tau(Num_Nodes) - Tau(Num_Nodes-1))
        dtau = 1d0/dtau

        ! Get difference
        dif = abs(Var(Num_Nodes) - Value_Const)

        ! If difference larger than 1 unit
        if (dif.gt.1d0) then

          ! Square root
          tmp = sqrt(dif)

          ! Residual
          Resi(Num_Nodes) = tmp*dtau

          ! Derivative
          if (Var(Num_Nodes).gt.Value_Const) then
            LL(Num_Nodes,Num_Nodes) =  0.5d0*dtau/tmp
          else
            LL(Num_Nodes,Num_Nodes) = -0.5d0*dtau/tmp
          end if

        ! Difference within a unit
        else

          ! Resiaudl and derivative
          Resi(Num_Nodes) = (Var(Num_Nodes) - Value_Const)*dtau
          LL(Num_Nodes,Num_Nodes) = dtau

        end if ! Difference size

        !
        ! Rest
        do i=2,Num_Nodes-1

          ! Centered step
          dtau = 0.5d0*abs(Tau(i+1) - Tau(i-1))
          dtau = 1d0/dtau

          ! Get difference
          dif = abs(Var(1) - Value_Const)

          ! If difference larger than 1 unit
          if (dif.gt.1d0) then

            ! Square root
            tmp = sqrt(dif)

            ! Residual
            Resi(i) = tmp*dtau

            ! Derivative
            if (Var(i).gt.Value_Const) then
              LL(i,i) =  0.5d0*dtau/tmp
            else
              LL(i,i) = -0.5d0*dtau/tmp
            end if

          ! Difference within a unit
          else

            ! Resiaudl and derivative
            Resi(i) = (Var(i) - Value_Const)*dtau
            LL(i,i) = dtau

          end if ! Difference size

        end do ! Nodes

      end if ! Number of nodes

      return

      end subroutine Constl1_Regul

!#####################################################################
!#####################################################################
!#####################################################################

      end module regul_mod
