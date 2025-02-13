      !> Transition and relaxation rates
      module seeaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     26/04/2017
!  Last version:
!     19/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     19/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  rE
!    Relaxation rate for spontaneous emission
!
!  rS
!    Relaxation rate for stimulated emission
!
!  rA
!    Relaxation rate for absorption
!
!  tE
!    Transition rate for spontaneous emission
!
!  tS
!    Transition rate for stimulated emission
!
!  tA
!    Transition rate for absorption
!
!  rSC
!    Relaxation rate for isotropic superelastic collisions
!
!  rAC
!    Relaxation rate for isotropic inelastic collisions
!
!  tSC
!    Transition rate for isotropic superelastic collisions
!
!  tAC
!    Transition rate for isotropic inelastic collisions
!
!  rSFC
!    Relaxation rate for isotropic forbidden superelastic collisions
!
!  rAFC
!    Relaxation rate for isotropic forbidden inelastic collisions
!
!  tSFC
!    Transition rate for isotropic forbidden superelastic collisions
!
!  tAFC
!    Transition rate for isotropic forbidden inelastic collisions
!
!  rEP
!    Relaxation rate for bound-free spontaneous emission
!
!  rSP
!    Relaxation rate for bound-free stimulated emission
!
!  rAP
!    Relaxation rate for bound-free absorption
!
!  tEP
!    Transition rate for bound-free spontaneous emission
!
!  tSP
!    Transition rate for bound-free stimulated emission
!
!  tAP
!    Transition rate for bound-free absorption
!
!  GammaF
!    Gamma in the magnetic kernel for the multi-term atom
!
!  MK
!    Magnetic kernel of a multi-term atom
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
      !!       iterm(integer): Term index\n
      !!  Ecoeff(double(:,:)): Einstein coefficient data\n
      !!      rEcoeff(double): Spontaneous emission relaxation rate
      subroutine rE(iterm,Ecoeff,rEcoeff)

      ! I/O

      integer, intent(in):: iterm
      double precision, intent(out):: rEcoeff
      double precision, dimension(:,:), intent(in):: Ecoeff

      ! Local

      integer:: iterml

      ! Initialize
      rEcoeff = 0d0

      ! For each lower term
      do iterml=1,iterm-1

        ! Add Aul
        rEcoeff = rEcoeff + Ecoeff(iterm,iterml)

      end do ! Lower terms

      end subroutine rE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for stimulated emission\n
      !!         iterm(integer): Term index\n
      !!       irad(integer(:)): Transition indexes\n
      !!         ntran(integer): Number of transitions\n
      !!       rLval(double(:)): Orbital angular momentum values\n
      !!             rL(double): Orbital angular momentum L\n
      !!              S(double): Spin S\n
      !!             rJ(double): Angular momentum J\n
      !!            rJ1(double): Angular momentum J'\n
      !!             rK(double): Multipolar component K\n
      !!              Q(double): Multipolar component Q\n
      !!            rJJ(double): Angular momentum J''\n
      !!           rJJ1(double): Angular momentum J'''\n
      !!            rKK(double): Multipolar component K'\n
      !!             QQ(double): Multipolar component Q'\n
      !!            zJ(logical): Bool with diagonality J,J''\n
      !!           zJ1(logical): Bool with diagonality J',J'''\n
      !!      Ecoeff(double(:)): Einstein coefficient data\n
      !!  RadJS(dcomplx(:,:,:)): Radiation field tensors integrated
      !!                         over the emission profile\n
      !!        mKr(integer(:)): Maximum radiative multipole\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                         J-symbols\n
      !!        rScoeff(double): Stimulated emission relaxation rate
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


      ! Initialize
      tmp = 0d0

      ! For each lower term
      do iterml=1,iterm-1

        ! Get transition index
        itran = irad(iterml)

        ! Skip not valid
        if (itran.le.0) cycle

        ! Get orbital momenta
        rLl = rLval(iterml)

        ! Initialize second sum
        tmp1 = 0d0

        ! For each radiative multipole
        do Kr=0,mKr(itran)

          ! Get real K value
          rKr = dble(Kr)

          ! Initialize third sum
          tmp2 = 0d0

          ! For all possible Q values
          do iQr=-Kr,Kr

            ! Get real Q value
            Qr=dble(iQr)

            ! Add contribution
            tmp2 = fun3j(rK,rKK,rKr,Q,-QQ,Qr,Flgsg)* &
                   RadJS(iQr,Kr,itran) + tmp2

          end do ! Qr values

          ! Initialize factor
          factor = 0d0

          ! If diagonal in J
          if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                           fun6j(rL,rL,rKr,rJJ1,rJ1,S,Flgsg)* &
                           fun6j(rK,rKK,rKr,rJJ1,rJ1,rJ,Flgsg)

          ! If diagonal in J'
          if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK+rKr))* &
                            sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                            fun6j(rL,rL,rKr,rJJ,rJ,S,Flgsg)* &
                            fun6j(rK,rKK,rKr,rJJ,rJ,rJ1,Flgsg) + &
                            factor

           ! Add to second sum
          tmp1 = .5d0*factor*Flgsg%sg(Kr)*sqrt(2d0*rKr+1d0)* &
                  fun6j(rL,rL,rKr,1d0,1d0,rLl,Flgsg)*tmp2 + tmp1

        end do ! Kr values

        ! Add to first sum
        tmp = Flgsg%sg(nint(rLl+rJ-S))*(2d0*rLl+1d0)* &
              Ecoeff(iterml)*tmp1 + tmp

      end do ! Lower terms

      ! Complete rate
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
      !!      rLval(double(:)): Orbital angular momentum values\n
      !!       nMulti(integer): Number of terms in the atom\n
      !!            rL(double): Orbital angular momentum L\n
      !!             S(double): Spin S\n
      !!            rJ(double): Angular momentum J\n
      !!           rJ1(double): Angular momentum J'\n
      !!            rK(double): Multipolar component K\n
      !!             Q(double): Multipolar component Q\n
      !!           rJJ(double): Angular momentum J''\n
      !!          rJJ1(double): Angular momentum J'''\n
      !!           rKK(double): Multipolar component K'\n
      !!            QQ(double): Multipolar component Q'\n
      !!           zJ(logical): Bool with diagonality J,J''\n
      !!          zJ1(logical): Bool with diagonality J',J'''\n
      !!     Ecoeff(double(:)): Einstein coefficient data\n
      !!  RadJ(dcomplx(:,:,:)): Radiation field tensors integrated
      !!                        over the absorption profile\n
      !!       mKr(integer(:)): Maximum radiative multipole\n
      !!       mKc(integer(:)): Maximum atomic multipole\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!       rAcoeff(double): Absorption relaxation rate
      subroutine rA(iterm,irad,ntran,rLval,nMulti, &
                    rL,S,rJ,rJ1,rK,Q,rJJ,rJJ1,rKK,QQ, &
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


      ! Initialize
      tmp = 0d0

      ! For every upper term
      do itermu=iterm+1,nMulti

        ! Get transition index
        itran = irad(itermu)

        ! Skip invalid transition
        if (itran.le.0) cycle

        ! Get orbitan angular momentum
        rLu = rLval(itermu)

        ! Initialize second sum
        tmp1 = 0d0

        ! For each possible radiation multipole
        do Kr=0,mKr(itran)

          ! Get real value of multipole
          rKr = dble(Kr)

          ! Check J symbol rules
          if (nint(rKr+rKK).gt.max(mKc(iterm),mKc(itermu))) cycle

          ! Initialize third sum
          tmp2 = 0d0

          ! For every possible Q
          do iQr=-Kr,Kr

            ! Get real Q value
            Qr = dble(iQr)

            ! Add contribution to sum
            tmp2 = fun3j(rK,rKK,rKr,Q,-QQ,Qr,Flgsg)* &
                   RadJ(iQr,Kr,itran) + tmp2

          end do ! Qr values

          ! Initialize factor
          factor = 0d0

          ! If diagonal in J
          if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                           fun6j(rL,rL,rKr,rJJ1,rJ1,S,Flgsg)* &
                           fun6j(rK,rKK,rKr,rJJ1,rJ1,rJ,Flgsg)

          ! If diagonal in J'
          if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK+rKr))* &
                            sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                            fun6j(rL,rL,rKr,rJJ,rJ,S,Flgsg)* &
                            fun6j(rK,rKK,rKr,rJJ,rJ,rJ1,Flgsg) + &
                            factor

          ! Add contribution to second sum
          tmp1 = factor*sqrt(2d0*rKr+1d0)* &
                 fun6j(rL,rL,rKr,1d0,1d0,rLu,Flgsg)*tmp2 + tmp1

        end do ! Kr values

        ! Add contribution to first sum
        tmp = Flgsg%sg(nint(rLu+rJ-S))*Ecoeff(itermu)*tmp1 + tmp

      end do ! Upper terms

      ! Complete rate
      rAcoeff = .5d0*Flgsg%sg(1+nint(QQ))*(2d0*rL+1d0)* &
                sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rA

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for spontaneous emission\n
      !!          rL(double): Orbital angular momentum L\n
      !!          rJ(double): Angular momentum J\n
      !!         rJ1(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!         rLL(double): Orbital angular momentum L'\n
      !!         rJJ(double): Angular momentum J''\n
      !!        rJJ1(double): Angular momentum J'''\n
      !!           S(double): Spin S\n
      !!      Ecoeff(double): Einstein coefficient data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!     tEcoeff(double): Spontaneous emission transition rate
      subroutine tE(rL,rJ,rJ1,rK,rLL,rJJ,rJJ1,S,Ecoeff,Flgsg,tEcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: rL,rJ,rJ1,rK,S
      double precision, intent(in):: rLL,rJJ,rJJ1
      double precision, intent(in):: Ecoeff
      double precision, intent(out):: tEcoeff

      ! Get rate
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
      !!           rL(double): Orbital angular momentum L\n
      !!           rJ(double): Angular momentum J\n
      !!          rJ1(double): Angular momentum J'\n
      !!           rK(double): Multipolar component K\n
      !!            Q(double): Multipolar component Q\n
      !!          rLL(double): Orbital angular momentum L'\n
      !!          rJJ(double): Angular momentum J''\n
      !!         rJJ1(double): Angular momentum J'''\n
      !!          rKK(double): Multipolar component K'\n
      !!           QQ(double): Multipolar component Q'\n
      !!            S(double): Spin S\n
      !!       Ecoeff(double): Einstein coefficient data\n
      !!  RadJS(dcomplx(:,:)): Radiation field tensor integrated over
      !!                       the emission profile\n
      !!         mKr(integer): Maximum radiative multipole\n
      !!   Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                       J-symbols\n
      !!      tScoeff(double): Stimulated emission transfer rate
      subroutine tS(rL,rJ,rJ1,rK,Q,rLL,rJJ,rJJ1,rKK,QQ,S, &
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


      ! Initialize sum
      tmp = 0d0

      ! For each radiation multipole
      do Kr=0,mKr

        ! Ger real Kr value
        rKr = dble(Kr)

        ! Initialize second sum
        tmp1 = 0d0

        ! For each possible multipole
        do iQr=-Kr,Kr

          ! Get real Qr value
          Qr = dble(iQr)

          ! Add contribution to second sum
          tmp1 = fun3j(rK,rKK,rKr,-Q,QQ,-Qr,Flgsg)*RadJS(iQr,Kr) + &
                 tmp1

        end do ! Qr values

        ! Add contribution to sum
        tmp = Flgsg%sg(Kr)*sqrt(2d0*rKr+1d0)* &
              fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,rKr,Flgsg)*tmp1 + &
              tmp

      end do ! Kr values

      ! Complete rate
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

      !> Transition rate for absorption\n
      !!           rL(double): Orbital angular momentum L\n
      !!           rJ(double): Angular momentum J\n
      !!          rJ1(double): Angular momentum J'\n
      !!           rK(double): Multipolar component K\n
      !!            Q(double): Multipolar component Q\n
      !!          rLL(double): Orbital angular momentum L'\n
      !!          rJJ(double): Angular momentum J''\n
      !!         rJJ1(double): Angular momentum J'''\n
      !!          rKK(double): Multipolar component K'\n
      !!           QQ(double): Multipolar component Q'\n
      !!            S(double): Spin S\n
      !!  Ecoeff(double(:,:)): Einstein coefficient data\n
      !!   RadJ(dcomplx(:,:)): Radiation field tensor integrated over
      !!                       the absorption profile\n
      !!         mKr(integer): Maximum radiative multipole\n
      !!         mKc(integer): Maximum atomic multipole\n
      !!   Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                       J-symbols\n
      !!      tAcoeff(double): Absorption transfer rate
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


      ! Initialize sum
      tmp = 0d0

      ! For each radiation tensor multipole
      do Kr=0,mKr

        ! Get real Kr value
        rKr = dble(Kr)

        ! Check J symbol rules
        if (nint(rKr+rKK).gt.mKc) cycle

        ! Initialize second sum
        tmp1 = 0d0

        ! For possible Q values
        do iQr=-Kr,Kr

          ! Get real Qr value
          Qr = dble(iQr)

          ! Add contribution to second sum
          tmp1 = fun3j(rK,rKK,rKr,-Q,QQ,-Qr,Flgsg)*RadJ(iQr,Kr) + tmp1

        end do ! Qr values

        ! Add contribution to sum
        tmp = sqrt(2d0*rKr+1d0)* &
              fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,rKr,Flgsg)*tmp1 + &
              tmp

      end do ! Kr values

      ! Complete rate
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

      !> Relaxation rate for isotropic superelastic collisions\n
      !!      iterm(integer): Term index\n
      !!    icol(integer(:)): Transition indexes\n
      !!    rLval(double(:)): Orbital angular momentum values\n
      !!          rL(double): Orbital angular momentum L\n
      !!           S(double): Spin S\n
      !!          rJ(double): Angular momentum J\n
      !!         rJ1(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!           Q(double): Multipolar component Q\n
      !!         rJJ(double): Angular momentum J''\n
      !!        rJJ1(double): Angular momentum J'''\n
      !!         rKK(double): Multipolar component K'\n
      !!          QQ(double): Multipolar component Q'\n
      !!         zJ(logical): Bool with diagonality J,J''\n
      !!        zJ1(logical): Bool with diagonality J',J'''\n
      !!   Ccoeff(double(:)): Collisional rates data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!    rSCcoeff(double): Collisional relaxation rate
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


      ! Initialize rate
      rSCcoeff = 0d0

      ! Check J symbol rules
      if (abs(rK-rKK).gt..4d0) return

      ! Initialize sum
      tmp = 0d0

      ! For each lower term
      do iterml=1,iterm-1

        ! Skip if no rate
        if (icol(iterml).le.0) cycle

        ! Get orbital angular momentum
        rLl = rLval(iterml)

        ! Initialize factor
        factor = 0d0

        ! If diagonal in J
        if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                         fun6j(rL,rL,0d0,rJJ1,rJ1,S,Flgsg)* &
                         fun6j(rK,rKK,0d0,rJJ1,rJ1,rJ,Flgsg)

        ! If diagonal in J'
        if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK))* &
                          sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                          fun6j(rL,rL,0d0,rJJ,rJ,S,Flgsg)* &
                          fun6j(rK,rKK,0d0,rJJ,rJ,rJ1,Flgsg) + &
                          factor

        ! If same term, introduce factor
        if (iterml.eq.iterm) factor = factor*.5d0

        ! Add contribution to sum
        tmp = Flgsg%sg(nint(rLl+rJ-S))* &
              Ccoeff(iterml)*factor* &
              fun3j(rK,rKK,0d0,Q,-QQ,0d0,Flgsg)* &
              fun6j(rL,rL,0d0,1d0,1d0,rLl,Flgsg) + tmp

      end do ! Lower terms

      ! Complete rate
      rSCcoeff = .5d0*Flgsg%sg(1+nint(QQ))*(2d0*rL+1d0)* &
                 sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rSC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for isotropic inelastic collisions\n
      !!      iterm(integer): Term index\n
      !!    icol(integer(:)): Transition indexes\n
      !!    rLval(double(:)): Orbital angular momentum values\n
      !!     nMulti(integer): Number of terms in the atom\n
      !!          rL(double): Orbital angular momentum L\n
      !!           S(double): Spin S\n
      !!          rJ(double): Angular momentum J\n
      !!         rJ1(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!           Q(double): Multipolar component Q\n
      !!         rJJ(double): Angular momentum J''\n
      !!        rJJ1(double): Angular momentum J'''\n
      !!         rKK(double): Multipolar component K'\n
      !!          QQ(double): Multipolar component Q'\n
      !!         zJ(logical): Bool with diagonality J,J''\n
      !!        zJ1(logical): Bool with diagonality J',J'''\n
      !!   Ccoeff(double(:)): Collisional rate data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!    rACcoeff(double): Collisional relaxation rate
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


      ! Initialize rate
      rACcoeff = 0d0

      ! Check J symbol rule
      if (abs(rK-rKK).gt..4d0) return

      ! Initialize sum
      tmp = 0d0

      ! For every upper term
      do itermu=iterm+1,nMulti

        ! Skip if no valid collision
        if (icol(itermu).le.0) cycle

        ! Get orbital angular momentum
        rLu = rLval(itermu)

        ! Initialize factor
        factor = 0d0

        ! If diagonal in J
        if (zJ) factor = sqrt((2d0*rJ1+1d0)*(2d0*rJJ1+1d0))* &
                         fun6j(rL,rL,0d0,rJJ1,rJ1,S,Flgsg)* &
                         fun6j(rK,rKK,0d0,rJJ1,rJ1,rJ,Flgsg)

        ! If diagonal in J'
        if (zJ1) factor = Flgsg%sg(nint(rJJ-rJ1+rK+rKK))* &
                          sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0))* &
                          fun6j(rL,rL,0d0,rJJ,rJ,S,Flgsg)* &
                          fun6j(rK,rKK,0d0,rJJ,rJ,rJ1,Flgsg) + &
                          factor

        ! Add to sum
        tmp = Flgsg%sg(nint(rLu+rJ-S))* &
              Ccoeff(itermu)*factor* &
              fun3j(rK,rKK,0d0,Q,-QQ,0d0,Flgsg)* &
              fun6j(rL,rL,0d0,1d0,1d0,rLu,Flgsg) + tmp

      end do ! Upper terms

      ! Complete rate
      rACcoeff = .5d0*Flgsg%sg(1+nint(QQ))*(2d0*rL+1d0)* &
                 sqrt(3d0*(2d0*rK+1d0)*(2d0*rKK+1d0))*tmp

      end subroutine rAC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for isotropic superelastic collisions\n
      !!          rL(double): Orbital angular momentum L\n
      !!          rJ(double): Angular momentum J\n
      !!         rJ1(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!           Q(double): Multipolar component Q\n
      !!         rLL(double): Orbital angular momentum L'\n
      !!         rJJ(double): Angular momentum J''\n
      !!        rJJ1(double): Angular momentum J'''\n
      !!         rKK(double): Multipolar component K'\n
      !!          QQ(double): Multipolar component Q'\n
      !!           S(double): Spin S\n
      !!      Ccoeff(double): Collisional rate data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!    tSCcoeff(double): Collisional transition rate
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


      ! Initialize rate
      tSCcoeff = 0d0

      ! Check J symbol rule
      if (abs(rK-rKK).gt..4d0) return

      ! Get J symbols
      tmp = fun3j(rK,rKK,0d0,-Q,QQ,0d0,Flgsg)* &
            fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,0d0,Flgsg)

      ! Get rate
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

      !> Transition rate for isotropic inelastic collisions\n
      !!          rL(double): Orbital angular momentum L\n
      !!          rJ(double): Angular momentum J\n
      !!         rJ1(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!           Q(double): Multipolar component Q\n
      !!         rLL(double): Orbital angular momentum L'\n
      !!         rJJ(double): Angular momentum J''\n
      !!        rJJ1(double): Angular momentum J'''\n
      !!         rKK(double): Multipolar component K'\n
      !!          QQ(double): Multipolar component Q'\n
      !!           S(double): Spin S\n
      !!      Ccoeff(double): Collisional rate data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!    tACcoeff(double): Collisional transition rate
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


      ! Initialize rate
      tACcoeff = 0d0

      ! Check J symbol rule
      if (abs(rK-rKK).gt..4d0) return

      ! Get J symbols
      tmp = fun3j(rK,rKK,0d0,-Q,QQ,0d0,Flgsg)* &
            fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0,rK,rKK,0d0,Flgsg)

      ! Get rate
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

      !> Relaxation rate for isotropic forbidden superelastic
      !! collisions\n
      !!     ilevel(integer): Level index\n
      !!          K(integer): Multipolar component K\n
      !!           Q(double): Multipolar component Q\n
      !!  CcoeffJ(double(:)): Collisional rates data\n
      !!   cflag(integer(:)): Forbidden collision flag\n
      !!    rSCcoeff(double): Collisional relaxation rate
      subroutine rSFC(ilevel,K,CcoeffJ,cflag,rSCcoeff)

      ! I/O

      integer, intent(in):: ilevel,K
      integer, dimension(:), intent(in):: cflag
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rSCcoeff

      ! Local

      integer:: ilevell


      ! Initialize rate
      rSCcoeff = 0d0

      ! For every lower level
      do ilevell=1,ilevel-1

        ! If the transition is not forbidden or if the multipole is
        ! not zero, skip
        if (cflag(ilevell).lt.1.or. &
            (cflag(ilevell).gt.1.and.K.ne.0)) cycle

        ! Add to rate
        rSCcoeff = rSCcoeff + CcoeffJ(ilevell)

      end do ! Lower levels

      end subroutine rSFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for isotropic forbidden inelastic
      !! collisions\n
      !!     ilevel(integer): Level index\n
      !!          K(integer): Multipolar component K\n
      !!     nlevel(integer): Number of levels in atomic model\n
      !!  CcoeffJ(double(:)): Collisional rates data\n
      !!   cflag(integer(:)): Forbidden collision flag\n
      !!    rACcoeff(double): Collisional relaxation rate
      subroutine rAFC(ilevel,K,nlevel,CcoeffJ,cflag,rACcoeff)

      ! I/O

      integer, intent(in):: ilevel,nlevel,K
      integer, dimension(:), intent(in):: cflag
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rACcoeff

      ! Local

      integer:: ilevelu


      ! Initialize rate
      rACcoeff = 0d0

      ! For every upper level
      do ilevelu=ilevel+1,nlevel

        ! If not forbidden or if multipole is not zero, skip
        if (cflag(ilevelu).lt.1.or. &
            (cflag(ilevelu).gt.1.and.K.ne.0)) cycle

        ! Add contribution
        rACcoeff = rACcoeff + CcoeffJ(ilevelu)

      end do ! Upper levels

      end subroutine rAFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for isotropic forbidden superelastic
      !! collisions\n
      !!          rJ(double): Angular momentum J\n
      !!         rJJ(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!     CcoeffJ(double): Collisional rate data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!       fcol(logical): Include transfer rate for K != 0\n
      !!    tSCcoeff(double): Collisional transition rate
      subroutine tSFC(rJ,rJJ,rK,CcoeffJ,Flgsg,fcol,tSCcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: fcol
      double precision, intent(in):: rJ,rJJ,rK,CcoeffJ
      double precision, intent(out):: tSCcoeff

      ! Local

      integer:: K,Ktilde

      double precision:: rKtilde, f6j


      ! Get integer K
      K = nint(rK)

      ! If population transfer
      if (K.eq.0) then

        ! Get rate
        tSCcoeff=sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      ! If polarization transfer allowed
      else if (fcol) then

        ! Get K tilde
        Ktilde = nint(abs(rJ - rJJ))

        ! Ktilde >= 1
        if (Ktilde.lt.1) Ktilde = 1

        ! Get real K tilde
        rKtilde = dble(Ktilde)

        ! Get 6J
        f6j = fun6j(rJ,rJ,0d0,rJJ,rJJ,rKtilde,Flgsg)

        ! If invalid 6J
        if (abs(f6j).lt.TINYJS) then

          ! No rate
          tSCcoeff = 0d0

        ! Valid 6J
        else

          ! Get rate
          tSCcoeff=sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))* &
                   Flgsg%sg(K)*CcoeffJ* &
                   fun6j(rJ,rJ,rK,rJJ,rJJ,rKtilde,Flgsg)/f6j

        end if ! Valid 6J

      ! No polarization transfer allowed
      else

        ! No rate
        tSCcoeff = 0d0

      end if ! Multipole or allowed polarization transfer

      end subroutine tSFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for isotropic forbidden inelastic
      !! collisions\n
      !!          rJ(double): Angular momentum J\n
      !!         rJJ(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!     CcoeffJ(double): Collisional rate data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!       fcol(logical): Include transfer rate for K != 0\n
      !!    tACcoeff(double): Collisional transition rate
      subroutine tAFC(rJ,rJJ,rK,CcoeffJ,Flgsg,fcol,tACcoeff)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: fcol
      double precision, intent(in):: rJ,rJJ,rK,CcoeffJ
      double precision, intent(out):: tACcoeff

      ! Local

      integer:: K, Ktilde

      double precision:: rKtilde, f6j


      ! Get integer K value
      K = nint(rK)

      ! IF population transfer
      if (K.eq.0) then

        ! Get transfer rate
        tACcoeff = sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      ! If polarization transfer allowed
      else if (fcol) then

        ! Get K tilde
        Ktilde = nint(abs(rJ - rJJ))

        ! Ktilde >= 1
        if (Ktilde.lt.1) Ktilde = 1

        ! Get real value of K tilde
        rKtilde = dble(Ktilde)

        ! Get 6J symbol
        f6j = fun6j(rJ,rJ,0d0,rJJ,rJJ,rKtilde,Flgsg)

        ! If invalid 6J
        if (abs(f6j).lt.TINYJS) then

          ! No rate
          tACcoeff = 0d0

        ! Valid 6J
        else

          ! Get rate
          tACcoeff=sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))* &
                   Flgsg%sg(K)*CcoeffJ* &
                   fun6j(rJ,rJ,rK,rJJ,rJJ,rKtilde,Flgsg)/f6j

        end if ! Valid 6J

      ! No polarization transfer allowed
      else

        ! No rate
        tACcoeff = 0d0

      end if ! Population transfer or allowed polarization transfer

      end subroutine tAFC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for bound-free spontaneous emission\n
      !!    ilevel(integer): Level index\n
      !!  iphot(integer(:)): Transition indexes\n
      !!   phot(Phot_class): Structure with photoionization data\n
      !!        iz(integer): Height index\n
      !!   rEPcoeff(double): Recombination relaxation rate
      subroutine rEP(ilevel,iphot,phot,iz,rEPcoeff)

      ! I/O

      type(Phot_class), dimension(:), intent(in):: phot
      integer, intent(in):: ilevel,iz
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rEPcoeff

      ! Local

      integer:: ilevell, itran


      ! Initialize rate
      rEPcoeff = 0d0

      ! For every lower level
      do ilevell=1,ilevel-1

        ! Get transition index
        itran = iphot(ilevell)

        ! If not valid, skip
        if (itran.lt.1) cycle

        ! Add to rate
        rEPcoeff = rEPcoeff + phot(itran)%TEI(iz)

      end do ! Lower levels

      end subroutine rEP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for bound-free stimulated emission\n
      !!    ilevel(integer): Level index\n
      !!  iphot(integer(:)): Transition indexes\n
      !!   Intgr(double(:)): Photoionization intensity integral\n
      !!   rSPcoeff(double): Stimulated recombination relaxation rate
      subroutine rSP(ilevel,iphot,Intgr,rSPcoeff)

      ! I/O

      integer, intent(in):: ilevel
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rSPcoeff
      double precision, dimension(:), intent(in):: Intgr

      ! Local

      integer:: ilevell,itran


      ! Initialize rate
      rSPcoeff = 0d0

      ! For every lower level
      do ilevell=1,ilevel-1

        ! Get transition index
        itran = iphot(ilevell)

        ! Skip if not valid transition
        if (itran.lt.1) cycle

        ! Add to rate
        rSPcoeff = rSPcoeff + Intgr(itran)

      end do ! Lower levels

      end subroutine rSP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for bound-free absorption\n
      !!    ilevel(integer): Level index\n
      !!  iphot(integer(:)): Transition indexes\n
      !!    nlevel(integer): Number of levels in atomic model\n
      !!   Intgr(double(:)): Photoionization intensity integral\n
      !!   rAPcoeff(double): Photoionization relaxation rate
      subroutine rAP(ilevel,iphot,nlevel,Intgr,rAPcoeff)

      ! I/O

      integer, intent(in):: ilevel, nlevel
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rAPcoeff
      double precision, dimension(:), intent(in):: Intgr

      ! Local

      integer:: ilevelu,itran


      ! Initialize rate
      rAPcoeff = 0d0

      ! For every upper level
      do ilevelu=ilevel+1,nlevel

        ! Get transition index
        itran = iphot(ilevelu)

        ! Skip if no valid transition
        if (itran.lt.1) cycle

        ! Add to rate
        rAPcoeff = rAPcoeff + Intgr(itran)

      end do ! Upper levels

      end subroutine rAP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for bound-free spontaneous emission\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!     Intgr(double): Photoionization intensity integral\n
      !!  tEPcoeff(double): Spontaneous recombination transfer rate
      subroutine tEP(rJ,rJJ,Intgr,tEPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tEPcoeff


      ! Get rate
      tEPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tEP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for bound-free stimulated emission\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!     Intgr(double): Photoionization intensity integral\n
      !!  tSPcoeff(double): Stimulated recombination transfer rate
      subroutine tSP(rJ,rJJ,Intgr,tSPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tSPcoeff


      ! Get rate
      tSPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tSP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for bound-free absorption\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!     Intgr(double): Photoionization intensity integral\n
      !!  tAPcoeff(double): Photoionization transfer rate
      subroutine tAP(rJ,rJJ,Intgr,tAPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tAPcoeff


      ! Get rate
      tAPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tAP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Gamma in the magnetic kernel for the multi-term atom\n
      !!          rL(double): Orbital angular momentum L\n
      !!           S(double): Spin S\n
      !!          rJ(double): Angular momentum J\n
      !!         rJJ(double): Angular momentum J'\n
      !!         zJ(logical): Bool with diagonality J,J'\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      function GammaF(rL,S,rJ,rJJ,zJ,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: zJ
      double precision, intent(in):: rL,S,rJ,rJJ

      ! Local

      double precision GammaF


      ! Get gamma always active contribution
      GammaF = Flgsg%sg(nint(rL+S+rJ)+1)* &
               sqrt((2d0*rJ+1d0)*(2d0*rJJ+1d0)* &
               S*(S+1d0)*(2d0*S+1d0))* &
               fun6j(rJ,rJJ,1d0,S,S,rL,Flgsg)

      ! If diagonal in J, add other term
      if (zJ) GammaF = sqrt(rJ*(rJ+1d0)*(2d0*rJ+1d0)) + GammaF

      end function GammaF

!#####################################################################
!#####################################################################
!#####################################################################

      !> Magnetic kernel of a multi-term atom\n
      !!          rL(double): Orbital angular momentum L\n
      !!           S(double): Spin S\n
      !!          rJ(double): Angular momentum J\n
      !!         rJ1(double): Angular momentum J'\n
      !!          rK(double): Multipolar component K\n
      !!         rJJ(double): Angular momentum J''\n
      !!        rJJ1(double): Angular momentum J'''\n
      !!         rKK(double): Multipolar component K'\n
      !!           Q(double): Multipolar component Q\n
      !!         zJ(logical): Bool with diagonality J,J''\n
      !!        zJ1(logical): Bool with diagonality J',J'''\n
      !!         zK(logical): Bool with diagonality K,K'\n
      !!         dFS(double): Fine structure energy difference\n
      !!      larmor(double): Magnetic field in larmor frequency
      !!                      units\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!    rMKcoeff(double): Magnetic kernel
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


      ! Initialize rate
      rMKcoeff = 0d0

      ! Initialize sum
      tmp = 0d0

      ! If diagonal in J
      if (zJ) tmp = Flgsg%sg(nint(rK-rKK))* &
                    GammaF(rL,S,rJJ1,rJ1,zJ1,Flgsg)* &
                    fun6j(rK,rKK,1d0,rJJ1,rJ1,rJ,Flgsg)

      ! If diagonal in J'
      if (zJ1) tmp = GammaF(rL,S,rJ,rJJ,zJ,Flgsg)* &
                     fun6j(rK,rKK,1d0,rJJ,rJ,rJ1,Flgsg) + tmp

      ! If diagonal in any of the two
      if (zJ.or.zJ1) rMKcoeff = larmor*Flgsg%sg(nint(rJ+rJ1-Q))* &
                                sqrt((2d0*rK+1d0)*(2d0*rKK+1d0))* &
                                fun3j(rK,rKK,1d0,-Q,Q,0d0,Flgsg)*tmp

      ! If diagonal in both and in K
      if (zJ.and.zJ1.and.zK) rMKcoeff = dFS + rMKcoeff

      end subroutine MK

!#####################################################################
!#####################################################################
!#####################################################################

      end module seeaux_mod
