      !> Check node bounds
      module bounds_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!  Start:
!     23/02/2023
!  Last version:
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - Revised headers (TdPA)
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
!      Force parameters value to be within the specified bounds
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

      !> Force parameters value to be within the specified bounds\n
      !!  Node(Node_class): Structure with nodes data\n
      !!      num(integer): Number of nodes
      subroutine CheckBounds(Node,Num)

      ! I/O

      type(Node_class), intent(inout):: Node
      integer, intent(in):: Num

      ! Local

      integer:: i, j


      ! For each node
      do i=1,Num

        ! For each especial limit
        do j=1,Node%nebound

          ! If node location within the especial limits
          if (Node%H(i).ge.Node%ebound(1,j).and. &
              Node%H(i).le.Node%ebound(2,j)) then

            ! If below lower limit
            if (Node%Var(i).lt.Node%ebound(3,j)) then

              ! Force limit
              Node%Var(i) = Node%ebound(3,j)

            ! If above upper limit
            else if (Node%Var(i).gt.Node%ebound(4,j)) then

              ! Force limit
              Node%Var(i) = Node%ebound(4,j)

            end if ! Beyond limits

            ! Continue with next node
            cycle

          end if ! Node within especial limits

        end do ! Number of especial limits

        ! If below lower limit
        if(Node%Var(i).lt.Node%Bounds(1)) then

          ! Force limit
          Node%Var(i) = Node%Bounds(1)

        ! If above upper limit
        else if (Node%Var(i).gt.Node%Bounds(2)) then

          ! Force limit
          Node%Var(i) = Node%Bounds(2)

        end if ! Beyond limits

      end do ! Elements

      return

      end subroutine CheckBounds

!#####################################################################
!#####################################################################
!#####################################################################

      !> Fold azimuth value accounting for additional space beyond the
      !! 2pi range\n
      !!   Node(Node_class): Azimuth nodes\n
      !!       Num(integer): Number of azimuthal nodes
      subroutine FoldBounds(Node,Num)

      ! I/O
      type(Node_class), intent(inout):: Node
      integer, intent(in):: Num

      ! Local

      logical:: skipnormal

      integer:: i, jj

      double precision:: Delta, MDelta, pi2, pih
      double precision, dimension(:), allocatable:: ldelta


      ! If there are no azimuth nodes, just return
      if (Num.le.0) return

      ! Dynamic range default region
      Delta = Node%Bounds(2) - Node%Bounds(1)
      MDelta = Delta

      ! Special ranges
      do jj=1,Node%nebound

        ! Allocate ldelta if not already
        if (.not.allocated(ldelta)) allocate(ldelta(Node%nebound))

        ! Get local delta
        ldelta(jj) = Node%ebound(4,jj) - Node%ebound(3,jj)

        ! Update max range
        MDelta = max(MDelta,ldelta(jj))

      end do

      ! Get 2pi and pi/2
      pi2 = 2d0*PI - TINYA
      pih = PI*0.5d0

      ! If dynamic range lower than 2pi, then there is no folding
      ! to do
      if (MDelta.lt.pi2) return

      ! For each node
      do i=1,Num

        ! Flag as no skip normal limit
        skipnormal = .False.

        ! Special ranges
        do jj=1,Node%nebound

          ! If node located within this range
          if (Node%H(i).ge.Node%ebound(1,jj).and. &
              Node%H(i).le.Node%ebound(2,jj)) then

            ! Skip normal because there was a especial limit
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

        ! If skipping normal, continue with next node
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

      ! Free
      if (allocated(ldelta)) deallocate(ldelta)

      return

      end subroutine FoldBounds

!#####################################################################
!#####################################################################
!#####################################################################

      end module bounds_mod
