      !> Rotation matrices
      module rdmat_mod
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
!   Calculates the rotation matrix      d[rJ,rM1,rM](theta)
!   for -180 < theta < 180
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

      !> Computes rotation matrix d[J,M',M](theta)\n
      !!          rJ(dfloat): Angular momentum J\n
      !!         rM1(dfloat): Magnetic quantum number M'\n
      !!          rM(dfloat): Magnetic quantum number M\n
      !!  Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!       theta(dfloat): Angle of the rotation
      double precision function rdmat(rJ,rM1,rM,Flgsg,theta)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: rJ,rM,rM1,theta

      ! Local

      integer:: ifs,ifd,ifs1,ifd1,imd,imdf2
      integer:: i,i2,imin,imax,k,k1,k2

      double precision:: thalf,tmp,ss,cc,ess,ecc

      thalf = .5d0*theta

      ss = sin(thalf)
      cc = cos(thalf)

      ifs = nint(rJ+rM)
      ifd = nint(rJ-rM)
      ifs1 = nint(rJ+rM1)
      ifd1 = nint(rJ-rM1)
      imd = nint(rM1-rM)
      imdf2 = nint(rJ+rJ) - imd

      imax = min(ifs,ifd1)
      imin = max(0,-imd)

      tmp = 0d0

      do i=imin,imax

        i2 = i + i

        k1 = imdf2 - i2
        k2 = imd + i2

        ess = 0d0
        ecc = 0d0

        if(abs(real(theta*RAD)-1.8e2).gt.0e0) then
          ecc = cc**k1
        else
          if(k1.eq.0) ecc = 1d0
        end if

        if(abs(real(theta*RAD)).gt.0d0) then
          ess = ss**k2
        else
          if(k2.eq.0) ess = 1d0
        end if

        k = imd + i

        tmp = Flgsg%sg(k)*exp(-Flgsg%flg(k) - Flgsg%flg(ifs-i) &
                              -Flgsg%flg(ifd1-i) - Flgsg%flg(i))* &
              ecc*ess + tmp

      end do

      rdmat = exp(.5d0*(Flgsg%flg(ifs) + Flgsg%flg(ifd) + &
                        Flgsg%flg(ifs1) + Flgsg%flg(ifd1)))*tmp

      return

      end function rdmat

!#####################################################################
!#####################################################################
!#####################################################################

      end module rdmat_mod
