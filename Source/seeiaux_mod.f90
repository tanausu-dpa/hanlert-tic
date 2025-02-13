      !> Transfer and relaxation rates for populations
      module seeiaux_mod
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
!  rEI
!    Relaxation rate for spontaneous emission
!
!  rSI
!    Relaxation rate for stimulated emission
!
!  rAI
!    Relaxation rate for absorption
!
!  tEI
!    Transfer rate for spontaneous emission
!
!  tSI
!    Transfer rate for stimulated emission
!
!  tAI
!    Transfer rate for absorption
!
!  rSFCI
!    Relaxation rate for superelastic collisions
!
!  rAFCI
!    Relaxation rate for inelastic collisions
!
!  tSFCI
!    Transfer rate for superelastic collisions
!
!  tAFCI
!    Transfer rate for inelastic collisions
!
!  rEPI
!    Relaxation rate for bound-free spontaneous emission
!
!  rSPI
!    Relaxation rate for bound-free stimulated emission
!
!  rAPI
!    Relaxation rate for bound-free absorption
!
!  tEPI
!    Transfer rate for bound-free spontaneous emission
!
!  tSPI
!    Transfer rate for bound-free stimulated emission
!
!  tAPI
!    Transfer rate for bound-free absorption
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
      !!       iterm(integer): Term index\n
      !!  Ecoeff(double(:,:)): Einstein coefficient data\n
      !!      rEcoeff(double): Spontaneous emission relaxation rate
      subroutine rEI(iterm,Ecoeff,rEcoeff)

      ! I/O

      integer, intent(in):: iterm
      double precision, intent(out):: rEcoeff
      double precision, dimension(:,:), intent(in):: Ecoeff

      ! Local

      integer:: iterml


      ! Initialize rate
      rEcoeff = 0d0

      ! For every lower term
      do iterml=1,iterm-1

        ! Add Einstein coefficient
        rEcoeff = rEcoeff + Ecoeff(iterm,iterml)

      end do ! Lower terms

      end subroutine rEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for stimulated emission\n
      !!      iterm(integer): Term index\n
      !!         iJ(integer): Level within term index\n
      !!    irad(integer(:)): Transition indexes\n
      !!      fst(FST_class): Structure with fine structure transition
      !!                      data\n
      !!  ifst(integer(:,:)): Indexes for fine structure transition\n
      !!      nJ(integer(:)): Number of J levels\n
      !!    rFval(double(:)): Angular momentum values\n
      !!          rJ(double): Angular momentum J\n
      !!    RadJS(double(:)): Mean intensity integrated over emission
      !!                      profile\n
      !!     rScoeff(double): Stimulated emission relaxation rate
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


      ! Initialize rate
      rScoeff = 0d0

      ! For each lower term
      do iterml=1,iterm-1

        ! Get transition rate
        itran = irad(iterml)

        ! Skip if not valid transition
        if (itran.lt.1) cycle

        ! For each level in the lower term
        do iJl=1,nJ(iterml)

          ! Get transition component rate
          ftran = fst(itran)%irad(iJ,iJl)

          ! Skip if not valid transition
          if (ftran.lt.1) cycle

          ! Get rolling index
          fftran = ifst(ftran,itran)

          ! Get angular momentum
          rJl = rJval(iJl,iterml)

          ! Add contribution to rate
          rScoeff = rScoeff + (2d0*rJl+1d0)*RadJS(fftran)* &
                              fst(itran)%Blu(iJl,iJ)/ &
                              (2d0*rJ+1d0)

        end do ! Levels in lower term
      end do ! Lower terms

      end subroutine rSI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for absorption\n
      !!      iterm(integer): Term index\n
      !!         iJ(integer): Level within term index\n
      !!    irad(integer(:)): Transition indexes\n
      !!      fst(FST_class): Structure with fine structure transition
      !!                      data\n
      !!  ifst(integer(:,:)): Indexes for fine structure transition\n
      !!      nJ(integer(:)): Number of J levels\n
      !!     nMulti(integer): Number of terms in the atom\n
      !!     RadJ(double(:)): Mean intensity integrated over
      !!                      absorption profile\n
      !!     rAcoeff(double): Stimulated emission relaxation rate
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


      ! Initialize rate
      rAcoeff = 0d0

      ! For every upper level
      do itermu=iterm+1,nMulti

        ! Get transition index
        itran = irad(itermu)

        ! Skip if not valid transition
        if (itran.lt.1) cycle

        ! For every level in the upper level
        do iJu=1,nJ(itermu)

          ! Get transition component index
          ftran = fst(itran)%irad(iJu,iJ)

          ! skip if not valid transition
          if (ftran.lt.1) cycle

          ! Get rolling index
          fftran = ifst(ftran,itran)

          ! Add contribution to rate
          rAcoeff = rAcoeff + fst(itran)%Blu(iJ,iJu)* &
                              RadJ(fftran)

        end do ! Levels in upper term
      end do ! Upper terms

      end subroutine rAI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for spontaneous emission\n
      !!       rJ(double): Angular momentum J\n
      !!      rJJ(double): Angular momentum J'\n
      !!   Ecoeff(double): Einstein coefficient data\n
      !!  tEcoeff(double): Spontaneous emission transition rate
      subroutine tEI(rJ,rJJ,Ecoeff,tEcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ,Ecoeff
      double precision, intent(out):: tEcoeff


      ! Get rate
      tEcoeff = Ecoeff*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for stimulated emission\n
      !!       rJ(double): Angular momentum J\n
      !!      rJJ(double): Angular momentum J'\n
      !!   Ecoeff(double): Einstein coefficient data\n
      !!    RadJS(double): Mean intensity integrated over emission
      !!                   profile\n
      !!  tScoeff(double): Stimulated emission transfer rate
      subroutine tSI(rJ,rJJ,Ecoeff,RadJS,tScoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ,Ecoeff,RadJS
      double precision, intent(out):: tScoeff


      ! Get rate
      tScoeff = Ecoeff*sqrt((2d0*rJ+1d0)/(2d0*rJJ+1d0))*RadJS

      end subroutine tSI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for absorption\n
      !!           rJ(double): Angular momentum J\n
      !!          rJJ(double): Angular momentum J'\n
      !!  Ecoeff(double(:,:)): Einstein coefficient data\n
      !!         RadJ(double): Mean intensity integrated over
      !!                       absorption profile\n
      !!      tAcoeff(double): Absorption transfer rate
      subroutine tAI(rJ,rJJ,Ecoeff,RadJ,tAcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ,Ecoeff,RadJ
      double precision, intent(out):: tAcoeff


      ! Get rate
      tAcoeff = Ecoeff*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*RadJ

      end subroutine tAI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for superelastic collisions\n
      !!     ilevel(integer): Level index\n
      !!  CcoeffJ(double(:)): Collisional rates data\n
      !!    rSCcoeff(double): Collisional relaxation rate
      subroutine rSFCI(ilevel,CcoeffJ,rSCcoeff)

      ! I/O

      integer, intent(in):: ilevel
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rSCcoeff

      ! Local

      integer:: ilevell


      ! Initialize rate
      rSCcoeff = 0d0

      ! For every lower level
      do ilevell=1,ilevel-1

        ! Add contribution
        rSCcoeff = rSCcoeff + CcoeffJ(ilevell)

      end do ! Lower levels

      end subroutine rSFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for inelastic collisions\n
      !!     ilevel(integer): Level index\n
      !!     nlevel(integer): Number of levels in atomic model\n
      !!  CcoeffJ(double(:)): Collisional rates data\n
      !!    rACcoeff(double): Collisional relaxation rate
      subroutine rAFCI(ilevel,nlevel,CcoeffJ,rACcoeff)

      ! I/O

      integer, intent(in):: ilevel,nlevel
      double precision, dimension(:), intent(in):: CcoeffJ
      double precision, intent(out):: rACcoeff

      ! Local

      integer:: ilevelu


      ! Initialize rate
      rACcoeff = 0d0

      ! For every upper level
      do ilevelu=ilevel,nlevel

        ! Add contribution to rate
        rACcoeff = rACcoeff + CcoeffJ(ilevelu)

      end do ! Upper levels

      end subroutine rAFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for superelastic collisions\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!   CcoeffJ(double): Collisional rate data\n
      !!  tSCcoeff(double): Collisional transition rate
      subroutine tSFCI(rJ,rJJ,CcoeffJ,tSCcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ, CcoeffJ
      double precision, intent(out):: tSCcoeff


      ! Get rate
      tSCcoeff = sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      end subroutine tSFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for inelastic collisions\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!   CcoeffJ(double): Collisional rate data\n
      !!  tACcoeff(double): Collisional transition rate
      subroutine tAFCI(rJ,rJJ,CcoeffJ,tACcoeff)

      ! I/O

      double precision, intent(in):: rJ,rJJ, CcoeffJ
      double precision, intent(out):: tACcoeff


      ! Get rate
      tACcoeff = sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))*CcoeffJ

      end subroutine tAFCI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for bound-free spontaneous emission\n
      !!    ilevel(integer): Level index\n
      !!  iphot(integer(:)): Transition indexes\n
      !!   phot(Phot_class): Structure with photoionization data\n
      !!        iz(integer): Height index\n
      !!   rEPcoeff(double): Recombination relaxation rate
      subroutine rEPI(ilevel,iphot,phot,iz,rEPcoeff)

      ! I/O

      type(Phot_class), dimension(:), intent(in):: phot
      integer, intent(in):: ilevel,iz
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rEPcoeff

      ! Local

      integer:: ilevell,itran


      ! Initialize rate
      rEPcoeff = 0d0

      ! For each lower level
      do ilevell=1,ilevel-1

        ! Get transition index
        itran = iphot(ilevell)

        ! Skip if not valid transition
        if (itran.lt.1) cycle

        ! Add contribution to rate
        rEPcoeff = rEPcoeff + phot(itran)%TEI(iz)

      end do ! Lower levels

      end subroutine rEPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for bound-free stimulated emission\n
      !!    ilevel(integer): Level index\n
      !!  iphot(integer(:)): Transition indexes\n
      !!   Intgr(double(:)): Photoionization intensity integral\n
      !!   rSPcoeff(double): Stimulated recombination relaxation rate
      subroutine rSPI(ilevel,iphot,Intgr,rSPcoeff)

      ! I/O

      integer, intent(in):: ilevel
      integer, dimension(:), intent(in):: iphot
      double precision, intent(out):: rSPcoeff
      double precision, dimension(:), intent(in):: Intgr

      ! Local

      integer:: ilevell,itran


      ! Initialize rate
      rSPcoeff = 0d0

      ! For each possible lower level
      do ilevell=1,ilevel-1

        ! Get transition index
        itran = iphot(ilevell)

        ! Skip if no valid transition
        if (itran.lt.1) cycle

        ! Add contribution to rate
        rSPcoeff = rSPcoeff + Intgr(itran)

      end do ! Lower levels

      end subroutine rSPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Relaxation rate for bound-free absorption\n
      !!    ilevel(integer): Level index\n
      !!  iphot(integer(:)): Transition indexes\n
      !!    nlevel(integer): Number of levels in atomic model\n
      !!   Intgr(double(:)): Photoionization intensity integral\n
      !!   rAPcoeff(double): Photoionization relaxation rate
      subroutine rAPI(ilevel,iphot,nlevel,Intgr,rAPcoeff)

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

        ! Skip if no valid index
        if (itran.lt.1) cycle

        ! Add contribution to rate
        rAPcoeff = rAPcoeff + Intgr(itran)

      end do ! Upper levels

      end subroutine rAPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for bound-free spontaneous emission\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!     Intgr(double): Photoionization intensity integral\n
      !!  tEPcoeff(double): Spontaneous recombination transfer rate
      subroutine tEPI(rJ,rJJ,Intgr,tEPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tEPcoeff


      ! Get rate
      tEPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tEPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for bound-free stimulated emission\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!     Intgr(double): Photoionization intensity integral\n
      !!  tSPcoeff(double): Stimulated recombination transfer rate
      subroutine tSPI(rJ,rJJ,Intgr,tSPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tSPcoeff


      ! Get rate
      tSPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tSPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transfer rate for bound-free absorption\n
      !!        rJ(double): Angular momentum J\n
      !!       rJJ(double): Angular momentum J'\n
      !!     Intgr(double): Photoionization intensity integral\n
      !!  tAPcoeff(double): Photoionization transfer rate
      subroutine tAPI(rJ,rJJ,Intgr,tAPcoeff)

      ! I/O

      double precision, intent(in):: Intgr,rJ,rJJ
      double precision, intent(out):: tAPcoeff


      ! Get rate
      tAPcoeff = Intgr*sqrt((2d0*rJJ+1d0)/(2d0*rJ+1d0))

      end subroutine tAPI

!#####################################################################
!#####################################################################
!#####################################################################

      end module seeiaux_mod
