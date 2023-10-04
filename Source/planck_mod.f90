      !> Planck function
      module planck_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
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
!     03/18/2019:    V1.1.0 - Avoids exponential overflow by changing
!                             to the Wien limit (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!  Calculates planck function for freq and T
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

      !> Computes Planck function\n
      !!   freq(dfloat): Frequency\n
      !!      T(dfloat): Temperature
      double precision function planck(freq,T)

      ! I/O
      double precision, intent(in):: freq, T

      ! Local
      double precision:: pnbar_l, En, c2freq

      ! Energy part
      En = ConvF*freq*1d21*(2d0*c)*freq**2d0

      ! Exponential
      c2freq = c2*freq*1d4/T

      !
      ! Denominator

      ! Argument in Wien limit
      if (c2freq.gt.wien_limit) then

        pnbar_l = diexp(c2freq)

      else

        pnbar_l = 1d0/(ddexp(c2freq) - 1d0)

      end if

      ! Planck function
      planck = pnbar_l*En

      return

      end function planck

!#####################################################################
!#####################################################################
!#####################################################################

      end module planck_mod
