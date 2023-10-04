      !> Factorials and signs
      module fctsg_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     06/29/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     12/10/2019:    V1.1.2 - Nullifies the jagged arrays (TdPA)
!
!     11/19/2019:    V1.1.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     04/18/2017:    V1.0.0 - First version (TdPA)
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
!    Calculates and stores log-factorials and signum values to be
!  used in 3n-j routines
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Pre-computes factorials and signs\n
      !!    Flgsg(Fctsg_class): Structure with factorials and signs
      subroutine fctsg(Flgsg)

      ! I/O
      type(Fctsg_class), intent(inout):: Flgsg

      ! Local
      integer:: i, i1

      ! Allocations
      allocate(Flgsg%flg(0:nxdim))
      allocate(Flgsg%sg(-nxdim:nxdim))

      ! Initialize first element
      Flgsg%flg(0)=.0D0
      Flgsg%sg(0)=.1D1

      ! Calculate rest of elements
      do i=1,nxdim

        i1 = i - 1

        Flgsg%flg(i) = log(dble(i)) + Flgsg%flg(i1)

        Flgsg%sg(i)  = -Flgsg%sg(i1)
        Flgsg%sg(-i) =  Flgsg%sg(i)

      end do

      ! Nullify jagged arrays
      nullify(Flgsg%J3%d)
      nullify(Flgsg%J6%d)
      nullify(Flgsg%J9%d)

      ! Control
      call control

      return

      end subroutine fctsg

!#####################################################################
!#####################################################################
!#####################################################################

      end module fctsg_mod
