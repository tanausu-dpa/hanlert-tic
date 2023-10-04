      !> Geometrical tensors in the vertical reference frame
      module stens_mod
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
!     04/20/2017:    V1.0.0 - First working version (TdPA)
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
!     This subroutine calculates the (complex) spherical components of
!   the geometric tensor, T_Q^K(i), for K=0,1,2, in the S-frame
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only : cZero, cOne, cImag, sqrt2, sqrt3
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes geometrical tensor in the vectical reference frame\n
      !!       theta(dfloat): Polar angle of direction\n
      !!         phi(dfloat): Azimuthal angle of direction\n
      !!         gam(dfloat): Polarization angle on the normal plane\n
      !!  Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!   TS(dcmplx(:,:,:)): Geometrical tensor in the vertical
      !!                      reference frame
      subroutine Stens(theta,phi,gam,Flgsg,TS)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: theta, phi, gam
      complex(kind=8), intent(out):: TS(0:3,-2:2,0:2)

      ! Local

      integer:: i,K,iQ

      double precision:: sTh,cTh,sTh2,cTh2,s2Gm,c2Gm

      complex(kind=8):: cTT,cexpPh,cexpPh2


      sTh = sin(theta)
      cTh = cos(theta)

      sTh2 = sTh*sTh
      cTh2 = cTh*cTh

      cexpPh = exp(cImag*phi)
      cexpPh2 = cexpPh*cexpPh

      s2Gm = sin(2d0*gam)
      c2Gm = cos(2d0*gam)

      cTT = (.25d0*sqrt3)*cexpPh2

      ! Initialize TS(i,K,Q)
      TS = cZero

      ! K=0
      ! i=0
      TS(0,0,0) = cOne

      ! K=1
      ! i=3
      TS(3,0,1) = (sqrt3/sqrt2)*cTh*cOne
      TS(3,1,1) = -(.5d0*sqrt3)*sTh*cexpPh

      ! K=2
      ! i=0
      TS(0,0,2) = (.5d0/sqrt2)*(3d0*cTh2-1d0)*cOne
      TS(0,1,2) = TS(3,1,1)*cTh
      TS(0,2,2) = cTT*sTh2

      ! i=1
      TS(1,0,2) = -(1.5d0/sqrt2)*sTh2*c2Gm*cOne
      TS(1,1,2) = TS(3,1,1)*dcmplx(cTh*c2Gm,s2Gm)
      TS(1,2,2) = -cTT*dcmplx((1d0+cTh2)*c2Gm,2d0*cTh*s2Gm)

      ! i=2
      TS(2,0,2) = (1.5d0/sqrt2)*sTh2*s2Gm*cOne
      TS(2,1,2) = -TS(3,1,1)*dcmplx(cTh*s2Gm,-c2Gm)
      TS(2,2,2) = cTT*dcmplx((1d0+cTh2)*s2Gm,-2d0*cTh*c2Gm)

      !
      ! Fix the -Q components
      !
      do K=1,2
        do iQ=1,K
          do i=0,3

            TS(i,-iQ,K) = Flgsg%sg(iQ)*conjg(TS(i,iQ,K))

          end do
        end do
      end do

      end subroutine Stens

!#####################################################################
!#####################################################################
!#####################################################################

      end module stens_mod
