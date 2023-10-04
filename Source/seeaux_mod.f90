      !> Transition and relaxation rates
      module seeaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/26/2017
!  Last version:
!     09/29/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.2 - Updated to term- and transition-wise
!                             K cut limits (TdPA)
!
!     02/01/2023:    V3.0.1 - Bugfix: The expected dimension of the
!                             Jrad array in rS and rA was nxtran, but
!                             it should be just the dimension for the
!                             current atom (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added a double precision K in the
!                             forbidden collisional transfer rates
!                             with polarization (TdPA)
!
!     17/02/2021:    V1.1.3 - Added polarization transfer for
!                             forbidden collisions (TdPA)
!
!     08/03/2018:    V1.1.2 - Introduced the Kcut on rA and tA (TdPA)
!
!     10/30/2017:    V1.1.1 - Decided that tAFC and tSFC are zero for
!                             K!=0, as they are not a dipole type
!                             transition (TdPA)
!
!     10/25/2017:    V1.1.0 - In tAFC and tSFC, missing the 6j that
!                             corresponds to K!=0 (TdPA)
!
!     10/24/2017:    V1.0.3 - In rAFC, start loop from the next level,
!                             does not affect results because it was
!                             adding 0 (TdPA)
!
!     09/08/2017:    V1.0.2 - The relaxation rates for forbidden
!                             collisions do not take into account
!                             ionizing collisions if it is not a
!                             rho00 case (TdPA)
!
!     05/05/2017:    V1.0.1 - The factor variable was outside of the
!                             term loop in rAC and rSC, and it should
!                             be inside to be properly reset for each
!                             collisional transition (TdPA)
!                           - The multilevel collisional relaxation
!                             rates only add contributions that are
!                             considered forbidden (TdPA)
!
!     04/26/2017:    V1.0.0 - First version (TdPA)
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
!    This module contains the routines to calculate the rates for the
!  SEE
!
!  rE:
!    Relaxation spontaneous rate.
!  tE:
!    Spontaneous emission transition rate.
!  tS:
!    Stimulated emission transition rate.
!  tSC:
!    Up-down collisional transition rate.
!  tA:
!    Absorption transition rate.
!  tAC:
!    Down-up collisional transition rate.
!  rS:
!    Relaxation stimulated rate.
!  rA:
!    Relaxation absorption rate.
!  rSC:
!    Relaxation rates for superelastic collisions (isotropic)
!  rAC:
!    Relaxation rates for inelastic collisions (isotropic)
!  tSFC:
!    Transition rates for up-down collision (forbidden)
!  tAFC:
!    Transition rates for down-up collision (forbidden)
!  rSFC:
!    Relaxation rates for superelastic forbidden collisions
!  rAFC:
!    Relaxation rates for inelastic forbidden collisions
!  rEP:
!    Relaxation b-f spontaneous rate.
!  tEP:
!    Spontaneous emission b-f transition rate.
!  tSP:
!    Stimulated emission b-f transition rate.
!  tAP:
!    Absorption b-f transition rate.
!  rSP:
!    Relaxation b-f stimulated rate.
!  rAP:
!    Relaxation b-f absorption rate.
!  GammaF:
!    Used in MK
!  MK:
!    Off-diagonal magnetic kernel of the SEE
!
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use funnj_mod
      use parameters_mod , only : TINYJS
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for spontaneous emission\n
      !!        iterm(integer): Term index\n
      !!   Ecoeff(dfloat(:,:)): Einstein coefficient data\n
      !!       rEcoeff(dfloat): Spontaneous emission relaxation rate
      subroutine rE(iterm,Ecoeff,rEcoeff)

      ! I/O

      integer, intent(in):: iterm
      double precision, intent(out):: rEcoeff
      double precision, dimension(:,:), intent(in):: Ecoeff

      ! Local

      integer:: iterml

      rEcoeff = 0d0

      do iterml=1,iterm-1
        rEcoeff = rEcoeff + Ecoeff(iterm,iterml)
      end do

      end subroutine rE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for spontaneous emission\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!           rLL(dfloat): Orbital angular momentum L'\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!             S(dfloat): Spin S\n
      !!        Ecoeff(dfloat): Einstein coefficient data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!       tEcoeff(dfloat): Spontaneous emission transition rate
      subroutine tE(rL,rJ,rJ1,rK,rLL,rJJ,rJJ1,S,Ecoeff,Flgsg,tEcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: rL,rJ,rJ1,rK,S
      double precision, intent(in):: rLL,rJJ,rJJ1
      double precision, intent(in):: Ecoeff
      double precision, intent(out):: tEcoeff

      tEcoeff = Flgsg%sg(1+nint(rK+rJ1+rJJ1))* &
                (2d0*rLL+1d0)*Ecoeff* &
                sqrt((2d0*rJ+1d0)*(2d0*rJ1+1d0)* &
                (2d0*rJJ+1d0)*(2d0*rJJ1+1d0))* &
                fun6j(rJ,rJ1,rK,rJJ1,rJJ,1d0,Flgsg)* &
                fun6j(rLL,rL,1d0,rJ,rJJ,S,Flgsg)* &
                fun6j(rLL,rL,1d0,rJ1,rJJ1,S,Flgsg)

      end subroutine tE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for stimulated emission\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rLL(dfloat): Orbital angular momentum L'\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!             S(dfloat): Spin S\n
      !!        Ecoeff(dfloat): Einstein coefficient data\n
      !!   RadJS(complex(:,:)): Radiation field tensor integrated over
      !!                        emission profile\n
      !!          mKr(integer): Maximum radiative multipole\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!       tScoeff(dfloat): Stimulated emission transfer rate
      subroutine tS(rL,rJ,rJ1,rK,Q, &
                    rLL,rJJ,rJJ1,rKK,QQ,S, &
                    Ecoeff,RadJS,mKr,Flgsg,tScoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: mKr
      double precision, intent(in):: rL,rJ,rJ1,rK,Q,S
      double precision, intent(in):: rLL,rJJ,rJJ1,rKK,QQ
      double precision, intent(in):: Ecoeff
      complex(kind=8), dimension(-2:2,0:2), intent(in):: RadJS
      complex(kind=8), intent(out):: tScoeff

      ! Local

      integer:: Kr,iQr

      double precision:: rKr,Qr

      complex(kind=8):: tmp,tmp1


      tmp = 0d0

      do Kr=0,mKr

        rKr = dble(Kr)

        tmp1 = 0d0

        do iQr=-Kr,Kr

          Qr = dble(iQr)

          tmp1 = fun3j(rK,rKK,rKr,-Q,QQ,-Qr,Flgsg)*RadJS(iQr,Kr) + &
                 tmp1

        end do

        tmp = Flgsg%sg(Kr)*sqrt(2d0*rKr+1d0)* &
              fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,rKr,Flgsg)*tmp1 + &
              tmp

      end do

      tScoeff = Flgsg%sg(nint(rKK+QQ+rJJ1-rJJ))* &
                (2d0*rL+1d0)*Ecoeff* &
                sqrt(3d0*(2d0*rJ+1d0)*(2d0*rJ1+1d0)* &
                (2d0*rJJ+1d0)*(2d0*rJJ1+1d0)* &
                (2d0*rK+1d0)*(2d0*rKK+1d0))* &
                fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)* &
                fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)*tmp

      end subroutine tS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for collision up-down\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rLL(dfloat): Orbital angular momentum L'\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!             S(dfloat): Spin S\n
      !!        Ccoeff(dfloat): Collisional rate data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!      tSCcoeff(dfloat): Collisional transition rate
      subroutine tSC(rL,rJ,rJ1,rK,Q,rLL,rJJ,rJJ1,rKK,QQ,S, &
                     Ccoeff,Flgsg,tSCcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: rL,rJ,rJ1,rK,Q,S
      double precision, intent(in):: rLL,rJJ,rJJ1,rKK,QQ
      double precision, intent(in):: Ccoeff
      double precision, intent(out):: tSCcoeff

      ! Local

      double precision:: tmp

      tSCcoeff = 0d0

      if (abs(rK-rKK).gt..4d0) return

      tmp = fun3j(rK,rKK,0d0,-Q,QQ,0d0,Flgsg)* &
            fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,0d0,Flgsg)

      tSCcoeff = Flgsg%sg(nint(rKK+QQ+rJJ1-rJJ))* &
                 (2d0*rLL+1d0)*Ccoeff* &
                 sqrt(3d0*(2d0*rJ+1d0)*(2d0*rJ1+1d0)* &
                 (2d0*rJJ+1d0)*(2d0*rJJ1+1d0)* &
                 (2d0*rK+1d0)*(2d0*rKK+1d0))* &
                 fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)* &
                 fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)*tmp

      end subroutine tSC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for absorption\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rLL(dfloat): Orbital angular momentum L'\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!             S(dfloat): Spin S\n
      !!   Ecoeff(dfloat(:,:)): Einstein coefficient data\n
      !!    RadJ(complex(:,:)): Radiation field tensor integrated over
      !!                        absorption profile\n
      !!          mKr(integer): Maximum radiative multipole\n
      !!          mKc(integer): Maximum atomic multipole\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!       tAcoeff(dfloat): Absorption transfer rate
      subroutine tA(rL,rJ,rJ1,rK,Q, &
                    rLL,rJJ,rJJ1,rKK,QQ,S, &
                    Ecoeff,RadJ,mKr,mKc,Flgsg,tAcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: mKr,mKc
      double precision, intent(in):: rL,rJ,rJ1,rK,Q,S
      double precision, intent(in):: rLL,rJJ,rJJ1,rKK,QQ
      double precision, intent(in):: Ecoeff
      complex(kind=8), dimension(-2:2,0:2), intent(in):: RadJ
      complex(kind=8), intent(out):: tAcoeff

      ! Local

      integer:: Kr,iQr

      double precision:: rKr,Qr

      complex(kind=8):: tmp,tmp1


      tmp = 0d0

      do Kr=0,mKr

        rKr = dble(Kr)

        if(nint(rKr+rKK).gt.mKc)cycle

        tmp1 = 0d0

        do iQr=-Kr,Kr

          Qr = dble(iQr)

          tmp1 = fun3j(rK,rKK,rKr,-Q,QQ,-Qr,Flgsg)*RadJ(iQr,Kr) + tmp1

        end do

        tmp = sqrt(2d0*rKr+1d0)* &
              fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,rKr,Flgsg)*tmp1 + &
              tmp

      end do

      tAcoeff = Flgsg%sg(nint(rKK+QQ+rJJ1-rJJ))* &
                (2d0*rLL+1d0)*Ecoeff* &
                sqrt(3d0*(2d0*rJ+1d0)*(2d0*rJ1+1d0)* &
                (2d0*rJJ+1d0)*(2d0*rJJ1+1d0)* &
                (2d0*rK+1d0)*(2d0*rKK+1d0))* &
                fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)* &
                fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)*tmp

      end subroutine tA

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for collision down-up\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rLL(dfloat): Orbital angular momentum L'\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!             S(dfloat): Spin S\n
      !!        Ccoeff(dfloat): Collisional rate data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!      tACcoeff(dfloat): Collisional transition rate
      subroutine tAC(rL,rJ,rJ1,rK,Q,rLL,rJJ,rJJ1,rKK,QQ,S, &
                     Ccoeff,Flgsg,tACcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: rL,rJ,rJ1,rK,Q,S
      double precision, intent(in):: rLL,rJJ,rJJ1,rKK,QQ
      double precision, intent(in):: Ccoeff
      double precision, intent(out):: tACcoeff

      ! Local

      double precision tmp

      tACcoeff = 0d0

      if (abs(rK-rKK).gt..4d0) return

      tmp = fun3j(rK,rKK,0d0,-Q,QQ,0d0,Flgsg)* &
            fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,0d0,Flgsg)

      tACcoeff = Flgsg%sg(nint(rKK+QQ+rJJ1-rJJ))* &
                 (2d0*rLL+1d0)*Ccoeff* &
                 sqrt(3d0*(2d0*rJ+1d0)*(2d0*rJ1+1d0)* &
                 (2d0*rJJ+1d0)*(2d0*rJJ1+1d0)* &
                 (2d0*rK+1d0)*(2d0*rKK+1d0))* &
                 fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)* &
                 fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)*tmp

      end subroutine tAC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for stimulated emission\n
      !!        iterm(integer): Term index\n
      !!      irad(integer(:)): Transition indexes\n
      !!        ntran(integer): Number of transitions\n
      !!      rLval(dfloat(:)): Orbital angular momentum values\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!             S(dfloat): Spin S\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!           zJ(logical): Bool with diagonality J,J''\n
      !!          zJ1(logical): Bool with diagonality J',J'''\n
      !!     Ecoeff(dfloat(:)): Einstein coefficient data\n
      !! RadJS(complex(:,:,:)): Radiation field tensor integrated over
      !!                        emission profile\n
      !!       mKr(integer(:)): Maximum radiative multipole\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!       rScoeff(dfloat): Stimulated emission relaxation rate
      subroutine rS(iterm,irad,ntran,rLval,rL,S,rJ,rJ1,rK,Q,rJJ, &
                    rJJ1,rKK,QQ,zJ,zJ1,Ecoeff,RadJS,mKr,Flgsg,rScoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ,zJ1
      integer, intent(in):: iterm,ntran
      integer, dimension(:), intent(in):: irad, mKr
      double precision, intent(in):: rL,S,rJ,rJ1,rK,Q
      double precision, intent(in):: rJJ,rJJ1,rKK,QQ
      double precision, dimension(:), intent(in):: rLval
      double precision, dimension(:), intent(in):: Ecoeff
      complex(kind=8), dimension(-2:2,0:2,ntran), intent(in):: RadJS
      complex(kind=8), intent(out):: rScoeff

      ! Local

      integer:: iterml,itran,Kr,iQr

      double precision:: rKr,Qr,rLl

      complex(kind=8):: tmp,tmp1,tmp2,factor


      tmp = 0d0

      do iterml=1,iterm-1

        itran = irad(iterml)

        if (itran.ne.0) then

          rLl = rLval(iterml)

          tmp1 = 0d0

          do Kr=0,mKr(itran)

            rKr = dble(Kr)

            tmp2 = 0d0

            do iQr=-Kr,Kr

              Qr=dble(iQr)

              tmp2 = fun3j(rK,rKK,rKr,Q,-QQ,Qr,Flgsg)* &
                     RadJS(iQr,Kr,itran) + tmp2

            end do

            factor = 0d0

            if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                             fun6j(rL,rL,rKr,rJJ1,rJ1,S,Flgsg)* &
                             fun6j(rK,rKK,rKr,rJJ1,rJ1,rJ,Flgsg)

            if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK+rKr))* &
                              sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                              fun6j(rL,rL,rKr,rJJ,rJ,S,Flgsg)* &
                              fun6j(rK,rKK,rKr,rJJ,rJ,rJ1,Flgsg) + &
                              factor

            tmp1 = .5d0*factor*Flgsg%sg(Kr)*sqrt(2d0*rKr+1d0)* &
                    fun6j(rL,rL,rKr,1d0,1d0,rLl,Flgsg)*tmp2 + tmp1

          end do

          tmp = Flgsg%sg(nint(rLl+rJ-S))*(2d0*rLl+1d0)* &
                Ecoeff(iterml)*tmp1 + tmp

        end if

      end do

      rScoeff = Flgsg%sg(1+nint(QQ))* &
                sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for absorption\n
      !!        iterm(integer): Term index\n
      !!      irad(integer(:)): Transition indexes\n
      !!        ntran(integer): Number of transitions\n
      !!      rLval(dfloat(:)): Orbital angular momentum values\n
      !!       nMulti(integer): Number of terms in the atom\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!             S(dfloat): Spin S\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!           zJ(logical): Bool with diagonality J,J''\n
      !!          zJ1(logical): Bool with diagonality J',J'''\n
      !!     Ecoeff(dfloat(:)): Einstein coefficient data\n
      !!  RadJ(complex(:,:,:)): Radiation field tensor integrated over
      !!                        absorption profile\n
      !!       mKr(integer(:)): Maximum radiative multipole\n
      !!       mKc(integer(:)): Maximum atomic multipole\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!       rAcoeff(dfloat): Absorption relaxation rate
      subroutine rA(iterm,irad,ntran,rLval,nMulti, &
                    rL,S,rJ,rJ1,rK,Q, &
                    rJJ,rJJ1,rKK,QQ, &
                    zJ,zJ1,Ecoeff,RadJ,mKr,mKc,Flgsg,rAcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ,zJ1
      integer, intent(in):: iterm, nMulti, ntran
      integer, dimension(:), intent(in):: irad, mKr, mKc
      double precision, intent(in):: rL,S,rJ,rJ1,rK,Q
      double precision, intent(in):: rJJ,rJJ1,rKK,QQ
      double precision, dimension(:), intent(in):: rLval
      double precision, dimension(:), intent(in):: Ecoeff
      complex(kind=8), dimension(-2:2,0:2,ntran), intent(in):: RadJ
      complex(kind=8), intent(out):: rAcoeff

      ! Local

      integer:: itermu,itran,Kr,iQr

      double precision:: rKr,Qr,rLu

      complex(kind=8):: tmp,tmp1,tmp2,factor


      tmp = 0d0

      do itermu=iterm+1,nMulti

        itran = irad(itermu)

        if (itran.ne.0) then

          rLu = rLval(itermu)

          tmp1 = 0d0

          do Kr=0,mKr(itran)

            rKr = dble(Kr)

            if(nint(rKr+rKK).gt.max(mKc(iterm),mKc(itermu))) cycle

            tmp2 = 0d0

            do iQr=-Kr,Kr

              Qr = dble(iQr)

              tmp2 = fun3j(rK,rKK,rKr,Q,-QQ,Qr,Flgsg)* &
                     RadJ(iQr,Kr,itran) + tmp2
            end do

            factor = 0d0

            if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                             fun6j(rL,rL,rKr,rJJ1,rJ1,S,Flgsg)* &
                             fun6j(rK,rKK,rKr,rJJ1,rJ1,rJ,Flgsg)

            if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK+rKr))* &
                              sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                              fun6j(rL,rL,rKr,rJJ,rJ,S,Flgsg)* &
                              fun6j(rK,rKK,rKr,rJJ,rJ,rJ1,Flgsg) + &
                              factor

            tmp1 = factor*sqrt(2d0*rKr+1d0)* &
                   fun6j(rL,rL,rKr,1d0,1d0,rLu,Flgsg)*tmp2 + tmp1

          end do

          tmp = Flgsg%sg(nint(rLu+rJ-S))*Ecoeff(itermu)*tmp1 + tmp

        end if

      end do

      rAcoeff = .5d0*Flgsg%sg(1+nint(QQ))*(2d0*rL+1d0)* &
                sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rA

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for superelastic collisions (isotropic)\n
      !!        iterm(integer): Term index\n
      !!      icol(integer(:)): Transition indexes\n
      !!      rLval(dfloat(:)): Orbital angular momentum values\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!             S(dfloat): Spin S\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!           zJ(logical): Bool with diagonality J,J''\n
      !!          zJ1(logical): Bool with diagonality J',J'''\n
      !!     Ccoeff(dfloat(:)): Collisional rates data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!      rSCcoeff(dfloat): Collisional relaxation rate
      subroutine rSC(iterm,icol,rLval,rL,S,rJ,rJ1,rK,Q,rJJ,rJJ1, &
                     rKK,QQ,zJ,zJ1,Ccoeff,Flgsg,rSCcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ,zJ1
      integer, intent(in):: iterm
      integer, dimension(:), intent(in):: icol
      double precision, intent(in):: rL,S,rJ,rJ1,rK,Q
      double precision, intent(in):: rJJ,rJJ1,rKK,QQ
      double precision, intent(out):: rSCcoeff
      double precision, dimension(:), intent(in):: rLval
      double precision, dimension(:), intent(in):: Ccoeff

      ! Local

      integer:: iterml

      double precision:: rLl,tmp,factor


      rSCcoeff = 0d0

      if (abs(rK-rKK).gt..4d0) return

      tmp = 0d0

      do iterml=1,iterm-1

        if (icol(iterml).ne.0) then

          rLl = rLval(iterml)

          factor = 0d0

          if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                           fun6j(rL,rL,0d0,rJJ1,rJ1,S,Flgsg)* &
                           fun6j(rK,rKK,0d0,rJJ1,rJ1,rJ,Flgsg)

          if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK))* &
                            sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                            fun6j(rL,rL,0d0,rJJ,rJ,S,Flgsg)* &
                            fun6j(rK,rKK,0d0,rJJ,rJ,rJ1,Flgsg) + &
                            factor

          if (iterml.eq.iterm) factor = factor*.5d0

          tmp = Flgsg%sg(nint(rLl+rJ-S))* &
                Ccoeff(iterml)*factor* &
                fun3j(rK,rKK,0d0,Q,-QQ,0d0,Flgsg)* &
                fun6j(rL,rL,0d0,1d0,1d0,rLl,Flgsg) + tmp

        end if


      end do

      rSCcoeff = .5d0*Flgsg%sg(1+nint(QQ))*(2d0*rL+1d0)* &
                 sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rSC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for inelastic collisions (isotropic)\n
      !!        iterm(integer): Term index\n
      !!      icol(integer(:)): Transition indexes\n
      !!      rLval(dfloat(:)): Orbital angular momentum values\n
      !!       nMulti(integer): Number of terms in the atom\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!             S(dfloat): Spin S\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!            QQ(dfloat): Multipolar component Q'\n
      !!           zJ(logical): Bool with diagonality J,J''\n
      !!          zJ1(logical): Bool with diagonality J',J'''\n
      !!     Ccoeff(dfloat(:)): Collisional rate data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!      rACcoeff(dfloat): Collisional relaxation rate
      subroutine rAC(iterm,icol,rLval,nMulti,rL,S,rJ,rJ1,rK,Q,rJJ, &
                     rJJ1,rKK,QQ,zJ,zJ1,Ccoeff,Flgsg,rACcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ,zJ1
      integer, intent(in):: iterm,nMulti
      integer, dimension(:), intent(in):: icol
      double precision, intent(in):: rL,S,rJ,rJ1,rK,Q
      double precision, intent(in):: rJJ,rJJ1,rKK,QQ
      double precision, intent(out):: rACcoeff
      double precision, dimension(:), intent(in):: rLval
      double precision, dimension(:), intent(in):: Ccoeff

      ! Local

      integer:: itermu

      double precision:: rLu,tmp,factor


      rACcoeff = 0d0

      if (abs(rK-rKK).gt..4d0) return

      tmp = 0d0

      do itermu=iterm+1,nMulti

        if (icol(itermu).ne.0) then

          rLu = rLval(itermu)

          factor = 0d0

          if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                           fun6j(rL,rL,0d0,rJJ1,rJ1,S,Flgsg)* &
                           fun6j(rK,rKK,0d0,rJJ1,rJ1,rJ,Flgsg)

          if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK))* &
                            sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                            fun6j(rL,rL,0d0,rJJ,rJ,S,Flgsg)* &
                            fun6j(rK,rKK,0d0,rJJ,rJ,rJ1,Flgsg) + &
                            factor

          tmp = Flgsg%sg(nint(rLu+rJ-S))* &
                Ccoeff(itermu)*factor* &
                fun3j(rK,rKK,0d0,Q,-QQ,0d0,Flgsg)* &
                fun6j(rL,rL,0d0,1d0,1d0,rLu,Flgsg) + tmp

        end if

      end do

      rACcoeff = .5d0*Flgsg%sg(1+nint(QQ))*(2d0*rL+1d0)* &
                 sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rAC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for up-down collision (forbidden)\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!       CcoeffJ(dfloat): Collisional rate data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!         fcol(logical): Include transfer rate for K != 0\n
      !!      tSCcoeff(dfloat): Collisional transition rate
      subroutine tSFC(rJ,rJJ,rK,CcoeffJ,Flgsg,fcol,tSCcoeff)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: fcol
      double precision, intent(in):: rJ,rJJ,rK,CcoeffJ
      double precision, intent(out):: tSCcoeff

      ! Local

      integer:: K,Ktilde
      double precision:: rKtilde, f6j

      K = nint(rK)

      ! Population transfer
      if (K.eq.0) then

        tSCcoeff=sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      ! Polarization transfer allowed
      else if (fcol) then

        ! Get K tilde
        Ktilde = nint(abs(rJ - rJJ))

        ! Ktilde >= 1
        if (Ktilde.lt.1) Ktilde = 1

        rKtilde = dble(Ktilde)

        f6j = fun6j(rJ,rJ,0d0,rJJ,rJJ,rKtilde,Flgsg)

        ! If invalid 6J
        if (abs(f6j).lt.TINYJS) then

          tSCcoeff = 0d0

        else

          tSCcoeff=sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))* &
                   Flgsg%sg(K)*CcoeffJ* &
                   fun6j(rJ,rJ,rK,rJJ,rJJ,rKtilde,Flgsg)/f6j

        end if

      ! No polarization transfer allowed
      else

        tSCcoeff = 0d0

      end if

      end subroutine tSFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for down-up collision (forbidden)\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!       CcoeffJ(dfloat): Collisional rate data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!         fcol(logical): Include transfer rate for K != 0\n
      !!      tACcoeff(dfloat): Collisional transition rate
      subroutine tAFC(rJ,rJJ,rK,CcoeffJ,Flgsg,fcol,tACcoeff)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: fcol
      double precision, intent(in):: rJ,rJJ,rK,CcoeffJ
      double precision, intent(out):: tACcoeff

      ! Local

      integer:: K, Ktilde
      double precision:: rKtilde, f6j

      K = nint(rK)

      ! Population transfer
      if (K.eq.0) then

        tACcoeff = sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      ! Polarization transfer allowed
      else if (fcol) then

        ! Get K tilde
        Ktilde = nint(abs(rJ - rJJ))

        ! Ktilde >= 1
        if (Ktilde.lt.1) Ktilde = 1

        rKtilde = dble(Ktilde)

        f6j = fun6j(rJ,rJ,0d0,rJJ,rJJ,rKtilde,Flgsg)

        ! If invalid 6J
        if (abs(f6j).lt.TINYJS) then

          tACcoeff = 0d0

        else

          tACcoeff=sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))* &
                   Flgsg%sg(K)*CcoeffJ* &
                   fun6j(rJ,rJ,rK,rJJ,rJJ,rKtilde,Flgsg)/f6j

        end if

      ! No polarization transfer allowed
      else

        tACcoeff = 0d0

      end if

      end subroutine tAFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for superelastic forbidden collisions\n
      !!       ilevel(integer): Level index\n
      !!            K(integer): Multipolar component K\n
      !!             Q(dfloat): Multipolar component Q\n
      !!    CcoeffJ(dfloat(:)): Collisional rates data\n
      !!     cflag(integer(:)): Forbidden collision flag\n
      !!      rSCcoeff(dfloat): Collisional relaxation rate
      subroutine rSFC(ilevel,K,CcoeffJ,cflag,rSCcoeff)

      ! I/O

      integer, intent(in):: ilevel,K
      integer, dimension(:), intent(in):: cflag
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rSCcoeff

      ! Local

      integer:: ilevell

      rSCcoeff = 0d0

      do ilevell=1,ilevel-1

        if (cflag(ilevell).lt.1.or. &
            (cflag(ilevell).gt.1.and.K.ne.0)) cycle

        rSCcoeff = rSCcoeff + CcoeffJ(ilevell)

      end do

      end subroutine rSFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for inelastic forbidden collisions\n
      !!       ilevel(integer): Level index\n
      !!            K(integer): Multipolar component K\n
      !!       nlevel(integer): Number of levels in atomic model\n
      !!    CcoeffJ(dfloat(:)): Collisional rates data\n
      !!     cflag(integer(:)): Forbidden collision flag\n
      !!      rACcoeff(dfloat): Collisional relaxation rate
      subroutine rAFC(ilevel,K,nlevel,CcoeffJ,cflag,rACcoeff)

      ! I/O

      integer, intent(in):: ilevel,nlevel,K
      integer, dimension(:), intent(in):: cflag
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rACcoeff

      ! Local

      integer:: ilevelu

      rACcoeff = 0d0

      do ilevelu=ilevel+1,nlevel

        if (cflag(ilevelu).lt.1.or. &
            (cflag(ilevelu).gt.1.and.K.ne.0)) cycle

        rACcoeff = rACcoeff + CcoeffJ(ilevelu)

      end do

      end subroutine rAFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for spontaneous recombination\n
      !!       ilevel(integer): Level index\n
      !!     iphot(integer(:)): Transition indexes\n
      !!      phot(Phot_class): Structure with photoionization data\n
      !!           iz(integer): Height index\n
      !!      rEPcoeff(dfloat): Recombination relaxation rate
      subroutine rEP(ilevel,iphot,phot,iz,rEPcoeff)

      ! I/O

      type(Phot_class), dimension(:), intent(in):: phot
      integer, intent(in):: ilevel,iz
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rEPcoeff

      ! Local

      integer:: ilevell, itran

      rEPcoeff = 0d0

      do ilevell=1,ilevel-1

        itran = iphot(ilevell)

        if (itran.lt.1) cycle

        rEPcoeff = rEPcoeff + phot(itran)%TEI(iz)

      end do

      end subroutine rEP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for spontaneous recombination\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!         Intgr(dfloat): Photoionization intensity integral\n
      !!      tEPcoeff(dfloat): Spontaneous recombination transfer
      !!                        rate
      subroutine tEP(rJ,rJJ,Intgr,tEPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tEPcoeff

      tEPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tEP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for stimulated recombination\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!         Intgr(dfloat): Photoionization intensity integral\n
      !!      tSPcoeff(dfloat): Stimulated recombination transfer
      !!                        rate
      subroutine tSP(rJ,rJJ,Intgr,tSPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tSPcoeff

      tSPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tSP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for photoionization\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!         Intgr(dfloat): Photoionization intensity integral\n
      !!      tAPcoeff(dfloat): Photoionization transfer rate
      subroutine tAP(rJ,rJJ,Intgr,tAPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tAPcoeff

      tAPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tAP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for stimulated recombination\n
      !!       ilevel(integer): Level index\n
      !!     iphot(integer(:)): Transition indexes\n
      !!      Intgr(dfloat(:)): Photoionization intensity integral\n
      !!      rSPcoeff(dfloat): Stimulated recombination relaxation
      !!                        rate
      subroutine rSP(ilevel,iphot,Intgr,rSPcoeff)

      ! I/O

      integer, intent(in):: ilevel
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rSPcoeff
      double precision, dimension(:), intent(in):: Intgr

      ! Local

      integer:: ilevell,itran

      rSPcoeff = 0d0

      do ilevell=1,ilevel-1

        itran = iphot(ilevell)

        if (itran.lt.1) cycle

        rSPcoeff = rSPcoeff + Intgr(itran)

      end do

      end subroutine rSP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for photoionization\n
      !!       ilevel(integer): Level index\n
      !!     iphot(integer(:)): Transition indexes\n
      !!       nlevel(integer): Number of levels in atomic model\n
      !!      Intgr(dfloat(:)): Photoionization intensity integral\n
      !!      rAPcoeff(dfloat): Photoionization relaxation rate
      subroutine rAP(ilevel,iphot,nlevel,Intgr,rAPcoeff)

      ! I/O

      integer, intent(in):: ilevel, nlevel
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rAPcoeff
      double precision, dimension(:), intent(in):: Intgr

      ! Local

      integer:: ilevelu,itran

      rAPcoeff = 0d0

      do ilevelu=ilevel+1,nlevel

        itran = iphot(ilevelu)

        if (itran.lt.1) cycle

        rAPcoeff = rAPcoeff + Intgr(itran)

      end do

      end subroutine rAP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Gamma in the magnetic kernel\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!             S(dfloat): Spin S\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!           zJ(logical): Bool with diagonality J,J'\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs
      function GammaF(rL,S,rJ,rJJ,zJ,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ
      double precision, intent(in):: rL,S,rJ,rJJ

      ! Local

      double precision GammaF


      GammaF = Flgsg%sg(nint(rL+S+rJ)+1)* &
               sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0)* &
               S*(S+1d0)*(2d0*S+1d0))* &
               fun6j(rJ,rJJ,1d0,S,S,rL,Flgsg)

      if (zJ) GammaF = sqrt(rJ*(rJ+1d0)*(2d0*rJ+1d0)) + GammaF

      end function GammaF

!#####################################################################
!#####################################################################
!#####################################################################

      !> Off-diagonal magnetic kernel\n
      !!            rL(dfloat): Orbital angular momentum L\n
      !!             S(dfloat): Spin S\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJ1(dfloat): Angular momentum J'\n
      !!            rK(dfloat): Multipolar component K\n
      !!           rJJ(dfloat): Angular momentum J''\n
      !!          rJJ1(dfloat): Angular momentum J'''\n
      !!           rKK(dfloat): Multipolar component K'\n
      !!             Q(dfloat): Multipolar component Q\n
      !!           zJ(logical): Bool with diagonality J,J''\n
      !!          zJ1(logical): Bool with diagonality J',J'''\n
      !!           zK(logical): Bool with diagonality K,K'\n
      !!           dFS(dfloat): Fine structure energy difference\n
      !!        larmor(dfloat): Magnetic field in larmor frequency
      !!                        units\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!      rMKcoeff(dfloat): Magnetic kernel
      subroutine MK(rL,S,rJ,rJ1,rK,rJJ,rJJ1,rKK,Q, &
                    zJ,zJ1,zK,dFS,larmor,Flgsg,rMKcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ,zJ1,zK
      double precision, intent(in):: rL,S,rJ,rJ1,rJJ,rJJ1
      double precision, intent(in):: rK,rKK,Q,dFS,larmor
      double precision, intent(out):: rMKcoeff

      ! Local

      double precision:: tmp

      rMKcoeff = 0d0
      tmp = 0d0

      if (zJ) tmp = Flgsg%sg(nint(rK-rKK))* &
                    GammaF(rL,S,rJJ1,rJ1,zJ1,Flgsg)* &
                    fun6j(rK,rKK,1d0,rJJ1,rJ1,rJ,Flgsg)

      if (zJ1) tmp = GammaF(rL,S,rJ,rJJ,zJ,Flgsg)* &
                     fun6j(rK,rKK,1d0,rJJ,rJ,rJ1,Flgsg) + tmp

      if (zJ.or.zJ1) rMKcoeff = larmor*Flgsg%sg(nint(rJ+rJ1-Q))* &
                                sqrt((2d0*rK+1d0)*(2d0*rKK+1d0))* &
                                fun3j(rK,rKK,1d0,-Q,Q,0d0,Flgsg)*tmp

      if (zJ.and.zJ1.and.zK) rMKcoeff = dFS + rMKcoeff

      end subroutine MK

!#####################################################################
!#####################################################################
!#####################################################################

      end module seeaux_mod
