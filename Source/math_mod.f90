      !> Initialization of photoionization quantities
      module math_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     02/15/2019
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
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     03/18/2019:    V1.1.0 - Added ddexp (TdPA)
!
!     02/15/2019:    V1.0.0 - First version (TdPA)
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
!  ddexp:
!    Exponential with argument control for double precision, positive
!    argument
!
!  diexp:
!    Exponential with argument control for double precision, negative
!    argument
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only :  bigexp , smallexp , vbigexp , &
                                   vbigexpv

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Exponmential with argument control\n
      !!    x(dfloat): Argument of exponential (positive)
      double precision function ddexp(x)

      ! I/O

      double precision, intent(in):: x

      ! Control overflows
      if(x.gt.vbigexp)then

        ddexp = vbigexpv

      ! Normal values
      elseif (x.gt.smallexp)then

        ddexp = exp(x)

      ! Control underflow
      else

        ddexp = 1d0 + x + .5d0*x*x

      endif

      return

      end function ddexp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Exponmential with argument control, negative argument\n
      !!    x(dfloat): Argument of exponential (in abs)
      double precision function diexp(x)

      ! I/O

      double precision, intent(in):: x

      ! Control overflows
      if(x.gt.bigexp)then

        diexp = 0d0

      ! Normal values
      elseif (x.gt.smallexp)then

        diexp = exp(-x)

      ! Control underflow
      else

        diexp = 1d0 - x + .5d0*x*x

      endif

      return

      end function diexp

!#####################################################################
!#####################################################################
!#####################################################################

      end module math_mod

