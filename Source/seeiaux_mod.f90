      !> Transfer and relaxation rates for populations
      module seeiaux_mod
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
!    This module contains the routines to calculate the rates for the
!  SEE
!
!  rEI:
!    Relaxation b-b spontaneous rate.
!  rEPI:
!    Relaxation b-f spontaneous rate.
!  tEI:
!    Spontaneous emission b-b transition rate.
!  tEPI:
!    Spontaneous emission b-f transition rate.
!  tSI:
!    Stimulated emission b-b transition rate.
!  tSPI:
!    Stimulated emission b-f transition rate.
!  tAI:
!    Absorption b-b transition rate.
!  tAPI:
!    Absorption b-f transition rate.
!  rSI:
!    Relaxation b-b stimulated rate.
!  rAI:
!    Relaxation b-b absorption rate.
!  rSPI:
!    Relaxation b-f stimulated rate.
!  rAPI:
!    Relaxation b-f absorption rate.
!  tSFCI:
!    Transition rates for up-down collision
!  tAFCI:
!    Transition rates for down-up collision
!  rSFCI:
!    Relaxation rates for superelastic collisions
!  rAFCI:
!    Relaxation rates for inelastic collisions
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use funnj_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for spontaneous emission\n
      !!        iterm(integer): Term index\n
      !!   Ecoeff(dfloat(:,:)): Einstein coefficient data\n
      !!       rEcoeff(dfloat): Spontaneous emission relaxation rate
      subroutine rEI(iterm,Ecoeff,rEcoeff)

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

      end subroutine rEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for spontaneous recombination\n
      !!       ilevel(integer): Level index\n
      !!     iphot(integer(:)): Transition indexes\n
      !!      phot(Phot_class): Structure with photoionization data\n
      !!           iz(integer): Height index\n
      !!      rEPcoeff(dfloat): Recombination relaxation rate
      subroutine rEPI(ilevel,iphot,phot,iz,rEPcoeff)

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

      end subroutine rEPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for spontaneous emission\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!        Ecoeff(dfloat): Einstein coefficient data\n
      !!       tEcoeff(dfloat): Spontaneous emission transition rate
      subroutine tEI(rJ,rJJ,Ecoeff,tEcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ,Ecoeff
      double precision, intent(out):: tEcoeff

      tEcoeff = Ecoeff*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for spontaneous recombination\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!         Intgr(dfloat): Photoionization intensity integral\n
      !!      tEPcoeff(dfloat): Spontaneous recombination transfer
      !!                        rate
      subroutine tEPI(rJ,rJJ,Intgr,tEPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tEPcoeff

      tEPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tEPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for stimulated emission\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!        Ecoeff(dfloat): Einstein coefficient data\n
      !!         RadJS(dfloat): Mean intensity integrated over
      !!                        emission profile\n
      !!       tScoeff(dfloat): Stimulated emission transfer rate
      subroutine tSI(rJ,rJJ,Ecoeff,RadJS,tScoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ,Ecoeff,RadJS
      double precision, intent(out):: tScoeff

      tScoeff = Ecoeff*sqrt((2d0*rJ+1d0)/(2d0*rJJ+1d0))*RadJS

      end subroutine tSI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for stimulated recombination\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!         Intgr(dfloat): Photoionization intensity integral\n
      !!      tSPcoeff(dfloat): Stimulated recombination transfer
      !!                        rate
      subroutine tSPI(rJ,rJJ,Intgr,tSPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tSPcoeff

      tSPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tSPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for absorption\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!   Ecoeff(dfloat(:,:)): Einstein coefficient data\n
      !!          RadJ(dfloat): Mean intensity integrated over
      !!                        absorption profile\n
      !!       tAcoeff(dfloat): Absorption transfer rate
      subroutine tAI(rJ,rJJ,Ecoeff,RadJ,tAcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ,Ecoeff,RadJ
      double precision, intent(out):: tAcoeff

      tAcoeff = Ecoeff*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*RadJ

      end subroutine tAI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for photoionization\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!         Intgr(dfloat): Photoionization intensity integral\n
      !!      tAPcoeff(dfloat): Photoionization transfer rate
      subroutine tAPI(rJ,rJJ,Intgr,tAPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tAPcoeff

      tAPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tAPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for stimulated emission\n
      !!        iterm(integer): Term index\n
      !!           iJ(integer): Level within term index\n
      !!      irad(integer(:)): Transition indexes\n
      !!        fst(FST_class): Structure with fine structure
      !!                        transition data\n
      !!    ifst(integer(:,:)): Indexes for fine structure
      !!                        transition\n
      !!        nJ(integer(:)): Number of J levels\n
      !!      rFval(dfloat(:)): Angular momentum values\n
      !!            rJ(dfloat): Angular momentum J\n
      !!      RadJS(dfloat(:)): Mean intensity integrated over
      !!                        emission profile\n
      !!       rScoeff(dfloat): Stimulated emission relaxation rate
      subroutine rSI(iterm,iJ,irad,fst,ifst,nJ,rJval,rJ, &
                     RadJS,rScoeff)

      ! I/O

      type(FST_class), dimension(:), intent(in):: fst
      integer, intent(in):: iterm, iJ
      integer, dimension(:), intent(in):: irad, nJ
      integer, dimension(:,:), intent(in):: ifst
      double precision, intent(in):: rJ
      double precision, intent(out):: rScoeff
      double precision, dimension(:), intent(in):: RadJS
      double precision, dimension(:,:), intent(in):: rJval

      ! Local

      integer:: iterml,iJl,itran,ftran,fftran

      double precision:: rJl

      rScoeff = 0d0

      do iterml=1,iterm-1

        itran = irad(iterml)

        if (itran.lt.1) cycle

        do iJl=1,nJ(iterml)

          ftran = fst(itran)%irad(iJ,iJl)

          if (ftran.lt.1) cycle

          fftran = ifst(ftran,itran)

          rJl = rJval(iJl,iterml)

          rScoeff = rScoeff + (2d0*rJl+1d0)*RadJS(fftran)* &
                              fst(itran)%Blu(iJl,iJ)/ &
                              (2d0*rJ+1d0)

        end do

      end do

      end subroutine rSI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for absorption\n
      !!        iterm(integer): Term index\n
      !!           iJ(integer): Level within term index\n
      !!      irad(integer(:)): Transition indexes\n
      !!        fst(FST_class): Structure with fine structure
      !!                        transition data\n
      !!    ifst(integer(:,:)): Indexes for fine structure
      !!                        transition\n
      !!        nJ(integer(:)): Number of J levels\n
      !!       nMulti(integer): Number of terms in the atom\n
      !!       RadJ(dfloat(:)): Mean intensity integrated over
      !!                        absorption profile\n
      !!       rAcoeff(dfloat): Stimulated emission relaxation rate
      subroutine rAI(iterm,iJ,irad,fst,ifst,nJ,nMulti, &
                     RadJ,rAcoeff)

      ! I/O

      type(FST_class), dimension(:), intent(in):: fst
      integer, intent(in):: iterm, iJ, nMulti
      integer, dimension(:), intent(in):: irad, nJ
      integer, dimension(:,:), intent(in):: ifst
      double precision, intent(out):: rAcoeff
      double precision, dimension(:), intent(in):: RadJ

      ! Local

      integer:: itermu,iJu,itran,ftran,fftran

      rAcoeff = 0d0

      do itermu=iterm+1,nMulti

        itran = irad(itermu)

        if (itran.lt.1) cycle

        do iJu=1,nJ(itermu)

          ftran = fst(itran)%irad(iJu,iJ)

          if (ftran.lt.1) cycle

          fftran = ifst(ftran,itran)

          rAcoeff = rAcoeff + fst(itran)%Blu(iJ,iJu)* &
                              RadJ(fftran)

        end do

      end do

      end subroutine rAI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for stimulated recombination\n
      !!       ilevel(integer): Level index\n
      !!     iphot(integer(:)): Transition indexes\n
      !!      Intgr(dfloat(:)): Photoionization intensity integral\n
      !!      rSPcoeff(dfloat): Stimulated recombination relaxation
      !!                        rate
      subroutine rSPI(ilevel,iphot,Intgr,rSPcoeff)

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

      end subroutine rSPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for photoionization\n
      !!       ilevel(integer): Level index\n
      !!     iphot(integer(:)): Transition indexes\n
      !!       nlevel(integer): Number of levels in atomic model\n
      !!      Intgr(dfloat(:)): Photoionization intensity integral\n
      !!      rAPcoeff(dfloat): Photoionization relaxation rate
      subroutine rAPI(ilevel,iphot,nlevel,Intgr,rAPcoeff)

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

      end subroutine rAPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for up-down collision (forbidden)\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!       CcoeffJ(dfloat): Collisional rate data\n
      !!      tSCcoeff(dfloat): Collisional transition rate
      subroutine tSFCI(rJ,rJJ,CcoeffJ,tSCcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ, CcoeffJ
      double precision, intent(out):: tSCcoeff

      tSCcoeff = sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      end subroutine tSFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transition rate for down-up collision (forbidden)\n
      !!            rJ(dfloat): Angular momentum J\n
      !!           rJJ(dfloat): Angular momentum J'\n
      !!       CcoeffJ(dfloat): Collisional rate data\n
      !!      tACcoeff(dfloat): Collisional transition rate
      subroutine tAFCI(rJ,rJJ,CcoeffJ,tACcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ, CcoeffJ
      double precision, intent(out):: tACcoeff

      tACcoeff = sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      end subroutine tAFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for superelastic forbidden collisions\n
      !!       ilevel(integer): Level index\n
      !!    CcoeffJ(dfloat(:)): Collisional rates data\n
      !!      rSCcoeff(dfloat): Collisional relaxation rate
      subroutine rSFCI(ilevel,CcoeffJ,rSCcoeff)

      ! I/O

      integer, intent(in):: ilevel
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rSCcoeff

      ! Local

      integer:: ilevell

      rSCcoeff = 0d0

      do ilevell=1,ilevel-1

        rSCcoeff = rSCcoeff + CcoeffJ(ilevell)

      end do

      end subroutine rSFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for inelastic forbidden collisions\n
      !!       ilevel(integer): Level index\n
      !!       nlevel(integer): Number of levels in atomic model\n
      !!    CcoeffJ(dfloat(:)): Collisional rates data\n
      !!      rACcoeff(dfloat): Collisional relaxation rate
      subroutine rAFCI(ilevel,nlevel,CcoeffJ,rACcoeff)

      ! I/O

      integer, intent(in):: ilevel,nlevel
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rACcoeff

      ! Local

      integer:: ilevelu

      rACcoeff = 0d0

      do ilevelu=ilevel,nlevel

        rACcoeff = rACcoeff + CcoeffJ(ilevelu)

      end do

      end subroutine rAFCI

!#####################################################################
!#####################################################################
!#####################################################################

      end module seeiaux_mod
