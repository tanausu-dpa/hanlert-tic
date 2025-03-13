      !> Initialization of photoionization quantities
      module math_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     15/02/2019
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Revised header (TdPA)
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
!  ddexp
!    Exponential with argument control for double precision, positive
!  argument
!
!  diexp:
!    Exponential with argument control for double precision, negative
!  argument
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

      !> Exponential with argument control for double precision,
      !! positive argument\n
      !!   x(double): Positive argument of exponential
      double precision function ddexp(x)

      ! I/O

      double precision, intent(in):: x

      ! Control overflows
      if(x.gt.vbigexp)then

        ! Constant big exponential value
        ddexp = vbigexpv

      ! Normal values
      elseif (x.gt.smallexp)then

        ! Compute exponential
        ddexp = exp(x)

      ! Control underflow
      else

        ! Second order Taylor
        ddexp = 1d0 + x + .5d0*x*x

      endif ! Control argument

      return

      end function ddexp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Exponential with argument control for double precision,
      !! negative argument\n
      !!  x(double): Argument of exponential (in abs)
      double precision function diexp(x)

      ! I/O

      double precision, intent(in):: x

      ! Control overflows
      if(x.gt.bigexp)then

        ! Tends to zero
        diexp = 0d0

      ! Normal values
      elseif (x.gt.smallexp)then

        ! Calculate exponential
        diexp = exp(-x)

      ! Control underflow
      else

        ! Taylor series
        diexp = 1d0 - x + .5d0*x*x

      endif

      return

      end function diexp

!#####################################################################
!#####################################################################
!#####################################################################

      end module math_mod
