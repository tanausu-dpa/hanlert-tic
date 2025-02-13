      !> Rotation matrices
      module rdmat_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     20/04/2017
!  Last version:
!     17/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     17/12/2024:    V4.0.0 - Changed global version (TdPA)
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
!  rdmat
!    Calculate the rotation matrix d[J,M',M](theta) for
!  -180 < theta < 180
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only : RAD
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the rotation matrix d[J,M',M](theta) for
      !! -180 < theta < 180\n
      !!          rJ(dfloat): Angular momentum J\n
      !!         rM1(dfloat): Magnetic quantum number M'\n
      !!          rM(dfloat): Magnetic quantum number M\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!       theta(dfloat): Angle of the rotation
      double precision function rdmat(rJ,rM1,rM,Flgsg,theta)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: rJ,rM,rM1,theta

      ! Local

      integer:: ifs,ifd,ifs1,ifd1,imd,imdf2
      integer:: i,i2,imin,imax,k,k1,k2

      double precision:: thalf,tmp,ss,cc,ess,ecc


      ! Half the theta angle
      thalf = .5d0*theta

      ! Sine and cosine of half the theta angle
      ss = sin(thalf)
      cc = cos(thalf)

      ! Quantum number combinations
      ifs = nint(rJ+rM)
      ifd = nint(rJ-rM)
      ifs1 = nint(rJ+rM1)
      ifd1 = nint(rJ-rM1)
      imd = nint(rM1-rM)
      imdf2 = nint(rJ+rJ) - imd

      ! Index limits in summation
      imax = min(ifs,ifd1)
      imin = max(0,-imd)

      ! Initialize summation
      tmp = 0d0

      ! For each contribution
      do i=imin,imax

        ! Double the index
        i2 = i + i

        ! Compute k1 and k2
        k1 = imdf2 - i2
        k2 = imd + i2

        ! Initialize powers of sine and cosine
        ess = 0d0
        ecc = 0d0

        ! If angle is not 180º
        if (abs(real(theta*RAD)-1.8e2).gt.0e0) then

          ! Power or cosine
          ecc = cc**k1

        ! Otherwise
        else

          ! If k1 is zero, make equal to 1
          if(k1.eq.0) ecc = 1d0

        end if ! Theta is not 180º

        ! If angle is not zero
        if (abs(real(theta*RAD)).gt.0d0) then

          ! Power of sine
          ess = ss**k2

        ! Otherwise
        else

          ! If k2 is zero, make equal to 1
          if(k2.eq.0) ess = 1d0

        end if ! Theta is not 0

        ! Get k index
        k = imd + i

        ! Add to summ
        tmp = ecc*ess*Flgsg%sg(k)* &
              exp(-Flgsg%flg(k) - Flgsg%flg(ifs-i) &
                  -Flgsg%flg(ifd1-i) - Flgsg%flg(i)) + tmp

      end do ! Contributions

      ! Complete rotation matrix
      rdmat = exp(.5d0*(Flgsg%flg(ifs) + Flgsg%flg(ifd) + &
                        Flgsg%flg(ifs1) + Flgsg%flg(ifd1)))*tmp

      return

      end function rdmat

!#####################################################################
!#####################################################################
!#####################################################################

      end module rdmat_mod
