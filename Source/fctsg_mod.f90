      !> Factorials and signs
      module fctsg_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     18/04/2017
!  Last version:
!     03/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  Fctsg
!    Compute factorials in logarithm and the of the powers of -1
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

      !> Compute factorials in logarithm and the of the powers of -1\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      subroutine fctsg(Flgsg)

      ! I/O

      type(Fctsg_class), intent(inout):: Flgsg

      ! Local

      integer:: i, i1


      ! Allocations
      allocate(Flgsg%flg(0:nxdim))
      allocate(Flgsg%sg(-nxdim:nxdim))

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Flgsg%flg)
      MRAMc = MRAMc + 1d-6*sizeof(Flgsg%sg)

      ! Initialize first element
      Flgsg%flg(0)=.0D0
      Flgsg%sg(0)=.1D1

      ! Calculate rest of elements
      do i=1,nxdim

        ! Get index or last element
        i1 = i - 1

        ! Get logarithm of the factorial
        Flgsg%flg(i) = log(dble(i)) + Flgsg%flg(i1)

        ! Get power of -1
        Flgsg%sg(i)  = -Flgsg%sg(i1)
        Flgsg%sg(-i) =  Flgsg%sg(i)

      end do ! Elements to compute

      ! Nullify jagged arrays for memoization
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
