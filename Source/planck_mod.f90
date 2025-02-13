      !> Planck function
      module planck_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  planck
!    Calculate the Planck function for a given frequency and
!  temperature
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only : c , ConvF , c2 , wien_limit
      use math_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the Planck function for a given frequency and
      !! temperature\n
      !!   freq(double): Frequency\n
      !!      T(double): Temperature
      double precision function planck(freq,T)

      ! I/O

      double precision, intent(in):: freq,T

      ! Local

      double precision:: pnbar_l,En,c2freq


      ! Energy part
      En = ConvF*freq*1d21*(2d0*c)*freq**2d0

      ! Exponential
      c2freq = c2*freq*1d4/T

      !
      ! Denominator

      ! Argument in Wien limit
      if (c2freq.gt.wien_limit) then

        ! Get factor
        pnbar_l = diexp(c2freq)

      ! Argument not in Wien limit
      else

        ! Get factor
        pnbar_l = 1d0/(ddexp(c2freq) - 1d0)

      end if ! Wien limit

      ! Planck function
      planck = pnbar_l*En

      return

      end function planck

!#####################################################################
!#####################################################################
!#####################################################################

      end module planck_mod
