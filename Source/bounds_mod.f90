      !> Check node bounds
      module bounds_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/23/2023
!  Last version:
!     03/15/2023 V3.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/15/2023:    V3.0.1 - Added the possibility of adding specific
!                             bounds at particular optical depth
!                             ranges (TdPA)
!                           - Removed CheckBoundsturb and CheckBounds3
!                             because their functionality can be now
!                             accounted for in the input (TdPA)
!                           - Changed the algorithm for FoldBounds
!                             for something that just allows the
!                             azimuth to roll over boundaries (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/23/2023:    V0.0.0 - Started from 05/12/2020
!                             TIC@bounds_mod.f90 revision (TdPA)
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
!    CheckBounds:
!      Force parameters within bounds
!
!    FoldBounds:
!      Check azimuths cyclical behavior
!
!#####################################################################
!#####################################################################
!#####################################################################

      use parameters_mod, only: PI , c , TINYA
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Force parameter within bounds\n
      !!  Node(Node_class): Structure with nodes data\n
      !!      num(integer): Number of nodes
      subroutine CheckBounds(Node,Num)

      ! IO
      type(Node_class):: Node
      integer, intent(in):: Num

      ! Local
      integer:: i, j

      ! For each node
      do i=1,Num

        ! For each special limit
        do j=1,Node%nebound

          ! If node within the special limits
          if (Node%H(i).ge.Node%ebound(1,j).and. &
              Node%H(i).le.Node%ebound(2,j)) then

            ! Force lower limit
            if (Node%Var(i).lt.Node%ebound(3,j)) then

              Node%Var(i) = Node%ebound(3,j)

            ! Force upper limit
            else if (Node%Var(i).gt.Node%ebound(4,j)) then

              Node%Var(i) = Node%ebound(4,j)

            end if ! Beyond limits

            ! Continue with next node
            cycle

          end if ! Node within special limits

        end do ! Number of special limits

        ! Force lower limit
        if(Node%Var(i).lt.Node%Bounds(1)) then

          Node%Var(i) = Node%Bounds(1)

        ! Force upper limit
        else if (Node%Var(i).gt.Node%Bounds(2)) then

          Node%Var(i) = Node%Bounds(2)

        end if ! Beyond limits

      end do ! Elements

      return

      end subroutine CheckBounds

!#####################################################################
!#####################################################################
!#####################################################################

      !> Put azimuths between bounds with angle equivalences\n
      !!   Node(Node_class): Azimuth nodes\n
      !!       Num(integer): Number of azimuthal nodes
      subroutine FoldBounds(Node,Num)

      ! IO
      type(Node_class), intent(inout):: Node
      integer, intent(in):: Num

      ! Local
      logical:: skipnormal

      integer:: i, jj

      double precision:: Delta, MDelta, pi2, pih
      double precision, dimension(:), allocatable:: ldelta

      ! Nodes?
      if (Num.le.0) return

      ! Dynamic range default region
      Delta = Node%Bounds(2) - Node%Bounds(1)
      MDelta = Delta

      ! Special ranges
      do jj=1,Node%nebound

        ! Allocate
        if (.not.allocated(ldelta)) allocate(ldelta(Node%nebound))

        ! Get local delta
        ldelta(jj) = Node%ebound(4,jj) - Node%ebound(3,jj)

        ! Update max
        MDelta = max(MDelta,ldelta(jj))

      end do

      ! Get 2pi and pi/2
      pi2 = 2d0*PI - TINYA
      pih = PI*0.5d0

      ! If dynamic range lower than 2pi, then return
      if (MDelta.lt.pi2) return

      ! For each node
      do i=1,Num

        ! Flag
        skipnormal = .False.

        ! Special ranges
        do jj=1,Node%nebound

          ! If node within this range
          if (Node%H(i).ge.Node%ebound(1,jj).and. &
              Node%H(i).le.Node%ebound(2,jj)) then

            ! Skip normal because there was a special limit
            skipnormal = .True.

            ! If full range
            if (ldelta(jj).ge.pi2) then

              ! If closer than pi/2 to lower limit, fold
              if (Node%Var(i) - Node%ebound(1,jj).lt.pih) then

                ! Fold and exit
                Node%Var(i) = Node%Var(i) + 2d0*PI
                exit

              ! If closer than pi/2 to upper limit, fold
              else if (Node%Var(i) - Node%ebound(2,jj).gt.pih) then

                ! Fold and exit
                Node%Var(i) = Node%Var(i) - 2d0*PI
                exit

              end if ! Node within bounds
            end if ! Full range of phi

            ! Exit loop
            exit

          end if ! In this special range

        end do ! Special ranges

        ! If skipping normal, continue
        if (skipnormal) cycle

        ! If full range
        if (Delta.ge.pi2) then

          ! If closer than pi/2 to lower limit, fold
          if (Node%Var(i) - Node%Bounds(1).lt.pih) then

            ! Fold and continue
            Node%Var(i) = Node%Var(i) + 2d0*PI
            cycle

          ! If closer than pi/2 to upper limit, fold
          else if (Node%Var(i) - Node%Bounds(2).gt.pih) then

            ! Fold and continue
            Node%Var(i) = Node%Var(i) - 2d0*PI
            cycle

          end if ! Node within bounds
        end if ! Full range

      end do ! Nodes

      return

      end subroutine FoldBounds

!#####################################################################
!#####################################################################
!#####################################################################

      end module bounds_mod
