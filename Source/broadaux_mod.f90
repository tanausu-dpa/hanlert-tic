      !> Contributions to the line broadening
      module broadaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/19/2017
!  Last version:
!     08/07/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.2 - Added broad_vdw_LTE, broad_stk_LTE, and
!                             broad_lstk_LTE routines (TdPA)
!
!     07/13/2022:    V3.0.1 - The abundance is now checked from the
!                             Atmo structure, and thus there is no
!                             longer dependence on the chemicaux_mod
!                             module (TdPA)
!                           - The Barklem broadening is pre-processed
!                             elsewhere to avoid reading many times
!                             (potentially) the same files, even if
!                             they are really small. Consequently,
!                             this module is no longer dependent on
!                             the inter_mod module (TdPA)
!                           - The barklem_* routines are no longer
!                             needed and have been removed (TdPA)
!                           - The resouce argument in broad_vdw is
!                             no longer needed (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case return
!                             clauses have been added after every
!                             call to aborted or control (TdPA)
!                           - Initialize the success flag in the
!                             reading of Barklem tables to failure
!                             state (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     07/31/2020:    V1.3.2 - Changed format of exponential when
!                             storing parameters (TdPA)
!                           - Only the master warns about no
!                             continuum (TdPA)
!
!     03/12/2019:    V1.3.1 - Can write what should be the parameters
!                             for VdW  and Stark broadening (TdPA)
!                           - Stark broadening is parametric if there
!                             is no continuum level (TdPA)
!
!     02/20/2019:    V1.3.0 - New verbosity (TdPA)
!                           - Units are now opened in 100 (TdPA)
!                           - Added error handling for failure reading
!                             Barklem data (TdPA)
!
!     11/19/2018:    V1.2.1 - Bugfix: Missing factor 2 in Barklem
!                             broadening. The paper gives half width
!                             half maximum, and the required
!                             quantity was full width half maximum
!                             instead (TdPA)
!
!     06/05/2017:    V1.2.0 - Bugfix: Parametric form in VdW was
!                             multiplying by the amplitude scale
!                             factor twice (TdPA)
!                           - Bugfix: The results of the parametric
!                             option for Van der Waals parametric
!                             broadening were 10^4 larger than the
!                             correspondent RH option, for no reason,
!                             so I put back a 10^-4 factor (TdPA)
!
!     09/15/2017:    V1.1.1 - Variable path for Barklem files (TdPA)
!
!     09/08/2017:    V1.1.0 - Now you can have a single ion in your
!                             model as long as you use parametric
!                             input in Van der Waals and Stark. Stark
!                             linear has a hardcoded ionization energy
!                             if there is no H II (TdPA)
!
!     07/19/2017:    V1.0.3 - Typo in one comment (TdPA)
!
!     05/31/2017:    V1.0.2 - Barklem only valid for neutrals (TdPA)
!                           - The reason for not using Barklem is
!                             clarified (TdPA)
!
!     05/05/2017:    V1.0.1 - Check could be undefined if the
!                             parameters were wrong in subroutine
!                             broad_vdw (TdPA)
!
!     04/19/2017:    V1.0.0 - First version (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    The Helium abundance is taken from the hardwired table, does not
!  depend on the input atoms
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!    This program is the auxiliar that does the real calculations
!    for the collisional broadening of spectral lines
!
!  broad_vdw:
!    Calculates Van der Waals broadening contribution
!
!  broad_stk:
!    Calculates quadratic Stark effect broadening contribution
!
!  broad_lstk:
!    Calculates linear Stark effect broadening contribution
!
!  broad_vdw_LTE:
!    Calculates Van der Waals broadening contribution for LTE lines
!
!  broad_stk_LTE:
!    Calculates quadratic Stark effect broadening contribution for LTE
!  lines
!
!  broad_lstk_LTE:
!    Calculates linear Stark effect broadening contribution for LTE
!  lines
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : ryd , mem , mhm , mhem , kb , pi , &
                                  amu , qel , a4pieps0 , hplanck , &
                                  rb , pi4eps0 , armass
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Van der Waals broadening\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!          iterm(integer): Lower term of transition\n
      !!         iterm1(integer): Upper term of transition\n
      !!          itran(integer): Index of transition\n
      !!            damp(dfloat): Inverse lifetime due to Van der
      !!                          Waals\n
      !!         aparam(logical): Store parameters of damping
      subroutine broad_vdw(Atom,Atmo,iterm,iterm1,itran,damp, &
                           aparam)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: aparam
      integer, intent(in):: iterm,iterm1,itran
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz,itermc,i1

      double precision:: aryd, Z, sigma, alpha, Ahe
      double precision:: v02c, v02che, nhe, d1, C6, sigc
      double precision:: a1,b1,a2,b2
      double precision, dimension(4):: args

      ! Routine name
      urou = 'broad_vdw'

      ! Take the arguments in the input
      args = Atom%broad_args(:,itran)

      ! Get the helium abundance
      Ahe = Atmo%abund(2)

      ! Correct rydberg energy for mass shift and calculate
      ! relative velocities for H and He
      aryd = ryd/(1d0 + mem/Atom%rmass)
      v02c = 8d0*kb*(1d0 + Atom%rmass/mhm)/pi/Atom%rmass/amu
      v02che = 8d0*kb*(1d0 + Atom%rmass/mhem)/pi/Atom%rmass/amu

      !
      ! Find the first continuum
      !

      ! Initialize the index of the continuum so it can be used as
      ! flag
      itermc = -1

      ! Find the next continuum
      do i1=iterm1,Atom%nMulti

        if (Atom%stage(i1).gt.Atom%stage(iterm1)) then
          itermc = i1
          exit
        end if

      end do

      ! If we did not find the continuum
      if (itermc.lt.0.and.Atom%broad_type(itran).ne.2) then
        umsg = 'Could not find continuum in atom '// &
               Atom%Element//' and Van der Waals '// &
               'broadening not set to parametric.'
        call aborted
        return
      end if


      !
      ! Calculate the Van der Waals cross section C6
      !
      if (Atom%broad_type(itran).ne.2) then

        ! Energy part
        d1 = ryd*ryd/(Atom%TRfreq(itermc) - &
                       Atom%TRfreq(iterm1))/ &
                     (Atom%TRfreq(itermc) - &
                       Atom%TRfreq(iterm1))  - &
             ryd*ryd/(Atom%TRfreq(itermc) - &
                       Atom%TRfreq(iterm))/ &
                     (Atom%TRfreq(itermc) - &
                       Atom%TRfreq(iterm))

        ! Next ion charge
        Z = Atom%stage(iterm)

        ! VdW cross section. 1d6 to receive cgs populations
        C6 = 8.08d0*((5d0*qel*qel*a4pieps0*pi* &
             Z*Z*rb*rb*d1/pi4eps0/hplanck)**.4d0)*1d6

      end if ! If not parametric


      !
      ! Calculate the broadening
      !

      ! If the type of broadening if Barklem
      if (Atom%broad_type(itran).eq.0) then

        ! Sigma and alpha calculated in rBarklem
        sigma = Atom%broad_args(1,itran)
        alpha = Atom%broad_args(2,itran)

        ! broadening constant part. 1d6 for cgs populations.
        sigc = 2d0*rb*rb*(4d0/pi)**(0.5d0*alpha)* &
               GAMMA(0.5d0*(4d0 - alpha))*1d4*sigma* &
               (v02c/1d8)**(.5d0*(1d0 - alpha))*1d6

        ! If storing a parameters
        if (aparam) then
          b1 = 0.5*(1d0 - alpha)
          a1 = 1d8*sigc/((1d0 + mhm/Atom%rmass)**b1)
          b2 = .3d0
          a2 = 1d9*C6*(v02che**.3)/((1d0 + mhem/Atom%rmass)**b2)
          write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
            err=1100) itran,iterm1,iterm,'    barklem',a1,b1,a2,b2
        end if

        ! For each height
        do iz=1,NZ

          ! Barklem for H
          damp(iz) = sigc*(Atmo%T(iz)**(.5d0*(1d0 - alpha)))* &
                     Atmo%nh(iz,1) + damp(iz)

          ! Helium population from the atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then
            nhe = Atmo%nhe(iz,1)
          ! Helium population from the abundance
          else
            nhe = Atmo%nh(iz,1)*Ahe
          end if

          ! Unsold for Helium
          damp(iz) = C6*(v02che**.3d0)*nhe*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

      end if ! Barklem inputs


      ! If the type of broadening is Unsold
      if (Atom%broad_type(itran).eq.1) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then
            nhe = Atmo%nhe(iz,1)
          ! Helium ground level pop from abundance
          else
            nhe = Atmo%nh(iz,1)*Ahe
          end if

          ! Unsold broadening
          damp(iz) = C6*((v02c**.3d0)*Atmo%nh(iz,1)*args(1) + &
                    (v02che**.3d0)*nhe*args(3))*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

        ! If storing a parameters
        if (aparam) then
          b1 = .3d0
          a1 = 1d8*C6*(v02c**.3)/((1d0 + mhm/Atom%rmass)**b1)
          b2 = .3d0
          a2 = 1d9*C6*(v02che**.3)/((1d0 + mhem/Atom%rmass)**b2)
          write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
            err=1100) itran,iterm1,iterm,'     unsold',a1,b1,a2,b2
        end if

      end if ! Unsold broadening


      ! If the type of broadening if Parametric
      if (Atom%broad_type(itran).eq.2) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then
            nhe = Atmo%nhe(iz,1)
          ! Helium ground level pop from abundance
          else
            nhe = Atmo%nh(iz,1)*Ahe
          end if

          ! Parametric broadening
          damp(iz) = 1d-8*args(1)*Atmo%nh(iz,1)* &
        (((1d0 + mhm/Atom%rmass)*Atmo%T(iz))**args(2)) + &
                     1d-9*args(3)*nhe* &
        (((1d0 + mhem/Atom%rmass)*Atmo%T(iz))**args(4)) + &
                     damp(iz)

        end do ! Heights

        ! If storing a parameters
        if (aparam) then
          b1 = args(2)
          a1 = args(1)
          b2 = args(4)
          a2 = args(3)
          write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
            err=1100) itran,iterm1,iterm,' parametric',a1,b1,a2,b2
        end if

      end if ! Parametric broadening

      return

1100  umsg = 'Error writing avdwparam file'
      close(600)
      close(700)
      urou = 'broad_vdw'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine broad_vdw

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Quadratic Stark broadening\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!       iterm(integer): Lower term of transition\n
      !!      iterm1(integer): Upper term of transition\n
      !!       itran(integer): Index of transition\n
      !!         damp(dfloat): Inverse lifetime due to quadratic Stark
      !!                       effect\n
      !!      aparam(logical): Store parameters of damping
      subroutine broad_stk(Atom,Atmo,iterm,iterm1,itran,damp,aparam)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: aparam
      integer, intent(in):: iterm,iterm1,itran
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz,itermc

      double precision:: stk,aryd,Z,n1,n2,v02c,C4,sigc

      ! Routine name
      urou = 'broad_stk'

      ! Take the arguments in the input
      stk = Atom%broad_stark(itran)

      ! Constant broadening proportional to density of electrons
      if (stk.lt.0d0) then

        stk = abs(stk)*1d6

        do iz=1,NZ

          damp(iz) = stk*Atmo%ne(iz) + damp(iz)

        end do

        ! Store param
        if (aparam) write(700,'(i10,1x,es15.8)',err=1100) itran,stk

      ! Any positive number means non-constant
      else if (stk.gt.0d0) then

        ! Correct rydberg energy for mass shift and calculate
        ! relative velocities for electrons and ions
        aryd = ryd/(1d0 + mem/Atom%rmass)
        v02c = ((8d0*kb/pi/Atom%rmass/amu)**(1d0/6d0))* &
               (((1d0 + Atom%rmass/mem)**(1d0/6d0)) + &
                ((1d0 + Atom%rmass/armass)**(1d0/6d0)))

        ! Initialize the index of the continuum so it can be used as
        ! flag
        itermc = -1

        ! Find the first continuum
        do iz=iterm1,Atom%nMulti

          if (Atom%stage(iz).gt.Atom%stage(iterm1)) then
            itermc = iz
            exit
          end if

        end do

        ! If we did not find the continuum
        if (itermc.lt.0) then

          if (pid.eq.0) then
            write(umsg,'(A,1x,i4)') &
              ' # WARNING: Could not find continuum in atom '// &
              Atom%Element//'. The Stark broadening will use '// &
              'your parameter in transition',itran
            call verbose
          end if

          ! For each height
          do iz=1,NZ

            ! Damping due to Stark quadratic
            damp(iz) = stk*(Atmo%T(iz)**(1d0/6d0))*Atmo%ne(iz) + &
                       damp(iz)

          end do

        else

          !
          ! Calculate the Stark cross section C4
          !

          ! Charge of next ion
          Z = Atom%stage(iterm)

          ! Effective principal quantum numbers
          n1 = Z*sqrt(aryd/(Atom%TRfreq(itermc) - &
                            Atom%TRfreq(iterm)))
          n2 = Z*sqrt(aryd/(Atom%TRfreq(itermc) - &
                            Atom%TRfreq(iterm1)))

          ! C4 constant. 1d6 to receive cgs populations
          C4 = qel*qel*2d0*pi*rb*rb*rb* &
               ((n2*(5d0*n2*n2 + 1d0))**2d0 - &
                (n1*(5d0*n1*n1 + 1d0))**2d0)*stk/ &
               pi4eps0/18d0/Z/Z/Z/Z/hplanck

          ! Cross section
          sigc = 11.37d0*(C4**(2d0/3d0))*v02c*1d6

          ! Store param
          if (aparam) &
            write(700,'(i10,1x,i10,"-->",2x,i10,1x,es15.8)',err=1100)&
                  itran,iterm1,iterm,sigc

          ! For each height
          do iz=1,NZ

            ! Damping due to Stark quadratic
            damp(iz) = sigc*(Atmo%T(iz)**(1d0/6d0))*Atmo%ne(iz) + &
                       damp(iz)

          end do

        end if ! There is continuum

      end if ! Parameter

      return

1100  umsg = 'Error writing astkparam file'
      close(600)
      close(700)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine broad_stk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes linear Stark broadening\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!       iterm(integer): Lower term of transition\n
      !!      iterm1(integer): Upper term of transition\n
      !!         damp(dfloat): Inverse lifetime due to linear Stark
      !!                       effect
      subroutine broad_lstk(Atom,Atmo,iterm,iterm1,damp)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      integer, intent(in):: iterm,iterm1
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz, itermc

      double precision:: n1, n2, C, sigc, TRfreq

      ! Routine name
      urou = 'broad_lstk'

      ! If not Hydrogen return
      if (Atom%Element.ne.' H') return

      ! Initialize the index of the continuum so it can be used as
      ! flag
      itermc = -1

      ! Find the first continuum
      do iz=iterm1,Atom%nMulti

        if (Atom%stage(iz).gt.Atom%stage(iterm1)) then
          itermc = iz
          exit
        end if

      end do

      ! If we did not find the continuum
      if (itermc.lt.0) then
        if (pid.eq.0) then
          umsg = ' # Could not find continuum '// &
                 'in atom H, using hardcoded '// &
                 'ionization energy.'
          call verbose
        end if
        TRfreq = 1.0967877174307
      else
        TRfreq = Atom%TRfreq(itermc)
      end if

      ! Calculate principal quantum number
      n1 = dble(nint(sqrt(ryd/ &
                (TRfreq - Atom%TRfreq(iterm)))))
      n2 = dble(nint(sqrt(ryd/ &
                (TRfreq - Atom%TRfreq(iterm1)))))

      ! If the principal quantum numbers differ in one
      if (nint(abs(n2 - n1)).eq.1) then

        C = 0.642d0

      ! If differ in something different than 1
      else

        C = 1d0

      end if

      ! Cross section
      sigc = C*0.6d0*(n2*n2 - n1*n1)

      ! For each height
      do iz=1,NZ

        ! Contribution due to linear Stark effect
        damp(iz) = sigc*(Atmo%ne(iz)**(2d0/3d0)) + damp(iz)

      end do

      return

      end subroutine broad_lstk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Van der Waals broadening for a LTE line\n
      !!     line(LTEline_class): Structure with the LTE line data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!            damp(dfloat): Inverse lifetime due to Van der
      !!                          Waals
      subroutine broad_vdw_LTE(line,Atmo,damp)

      ! I/O

      type(LTEline_class):: line
      type(Atmo_class), intent(in):: Atmo
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz

      double precision:: aryd, Z, sigma, alpha, Ahe, Ei
      double precision:: v02c, v02che, nhe, d1, C6, sigc
      double precision, dimension(4):: args

      ! Routine name
      urou = 'broad_vdw_LTE'

      ! Take the arguments in the input
      args = line%broad_args

      ! If gamma, easy
      if (line%broad_type.eq.3) then
        damp = args(1)
        return
      end if

      ! Get the helium abundance
      Ahe = Atmo%abund(2)

      ! Correct rydberg energy for mass shift and calculate
      ! relative velocities for H and He
      aryd = ryd/(1d0 + mem/line%rmass)
      v02c = 8d0*kb*(1d0 + line%rmass/mhm)/pi/line%rmass/amu
      v02che = 8d0*kb*(1d0 + line%rmass/mhem)/pi/line%rmass/amu

      !
      ! Get ionization energy
      !
      Ei = Atmo%ele(line%ele)%Ei(line%stage)

      !
      ! Calculate the Van der Waals cross section C6
      !
      if (line%broad_type.ne.2) then

        ! Energy part
        d1 = ryd*ryd/(Ei - line%Eu)/(Ei - line%Eu)  - &
             ryd*ryd/(Ei - line%El)/(Ei - line%El)

        ! Next ion charge
        Z = line%stage

        ! VdW cross section. 1d6 to receive cgs populations
        C6 = 8.08d0*((5d0*qel*qel*a4pieps0*pi* &
             Z*Z*rb*rb*d1/pi4eps0/hplanck)**.4d0)*1d6

      end if ! If not parametric


      !
      ! Calculate the broadening
      !

      ! If the type of broadening if Barklem
      if (line%broad_type.eq.0) then

        ! Sigma and alpha calculated in rBarklem
        sigma = line%broad_args(1)
        alpha = line%broad_args(2)

        ! broadening constant part. 1d6 for cgs populations.
        sigc = 2d0*rb*rb*(4d0/pi)**(0.5d0*alpha)* &
               GAMMA(0.5d0*(4d0 - alpha))*1d4*sigma* &
               (v02c/1d8)**(.5d0*(1d0 - alpha))*1d6

        ! For each height
        do iz=1,NZ

          ! Barklem for H
          damp(iz) = sigc*(Atmo%T(iz)**(.5d0*(1d0 - alpha)))* &
                     Atmo%nh(iz,1) + damp(iz)

          ! Helium population from the atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then
            nhe = Atmo%nhe(iz,1)
          ! Helium population from the abundance
          else
            nhe = Atmo%nh(iz,1)*Ahe
          end if

          ! Unsold for Helium
          damp(iz) = C6*(v02che**.3d0)*nhe*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

      end if ! Barklem inputs

      ! If the type of broadening is Unsold
      if (line%broad_type.eq.1) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then
            nhe = Atmo%nhe(iz,1)
          ! Helium ground level pop from abundance
          else
            nhe = Atmo%nh(iz,1)*Ahe
          end if

          ! Unsold broadening
          damp(iz) = C6*((v02c**.3d0)*Atmo%nh(iz,1)*args(1) + &
                    (v02che**.3d0)*nhe*args(3))*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

      end if ! Unsold broadening

      ! If the type of broadening if Parametric
      if (line%broad_type.eq.2) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then
            nhe = Atmo%nhe(iz,1)
          ! Helium ground level pop from abundance
          else
            nhe = Atmo%nh(iz,1)*Ahe
          end if

          ! Parametric broadening
          damp(iz) = 1d-8*args(1)*Atmo%nh(iz,1)* &
        (((1d0 + mhm/line%rmass)*Atmo%T(iz))**args(2)) + &
                     1d-9*args(3)*nhe* &
        (((1d0 + mhem/line%rmass)*Atmo%T(iz))**args(4)) + &
                     damp(iz)

        end do ! Heights

      end if ! Parametric broadening

      return

      end subroutine broad_vdw_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Quadratic Stark broadening\n
      !!  line(LTEline_class): Structure with the LTE line data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!         damp(dfloat): Inverse lifetime due to quadratic Stark
      !!                       effect
      subroutine broad_stk_LTE(line,Atmo,damp)

      ! I/O

      type(LTEline_class):: line
      type(Atmo_class), intent(in):: Atmo
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz

      double precision:: stk,aryd,Z,n1,n2,v02c,C4,sigc,Ei

      ! Routine name
      urou = 'broad_stk_LTE'

      ! Take the arguments in the input
      stk = line%broad_stark

      ! Constant broadening proportional to density of electrons
      if (stk.lt.0d0) then

       !stk = abs(stk)*1d6
        stk = abs(stk)

        do iz=1,NZ

         !damp(iz) = stk*Atmo%ne(iz) + damp(iz)
          damp(iz) = stk + damp(iz)

        end do

      ! Any positive number means non-constant
      else if (stk.gt.0d0) then

        ! Correct rydberg energy for mass shift and calculate
        ! relative velocities for electrons and ions
        aryd = ryd/(1d0 + mem/line%rmass)
        v02c = ((8d0*kb/pi/line%rmass/amu)**(1d0/6d0))* &
               (((1d0 + line%rmass/mem)**(1d0/6d0)) + &
                ((1d0 + line%rmass/armass)**(1d0/6d0)))

        !
        ! Get ionization energy
        !
        Ei = Atmo%ele(line%ele)%Ei(line%stage)

        !
        ! Calculate the Stark cross section C4
        !

        ! Charge of next ion
        Z = line%stage

        ! Effective principal quantum numbers
        n1 = Z*sqrt(aryd/(Ei - line%El))
        n2 = Z*sqrt(aryd/(Ei - line%Eu))

        ! C4 constant. 1d6 to receive cgs populations
        C4 = qel*qel*2d0*pi*rb*rb*rb* &
             ((n2*(5d0*n2*n2 + 1d0))**2d0 - &
              (n1*(5d0*n1*n1 + 1d0))**2d0)*stk/ &
             pi4eps0/18d0/Z/Z/Z/Z/hplanck

        ! Cross section
        sigc = 11.37d0*(C4**(2d0/3d0))*v02c*1d6

        ! For each height
        do iz=1,NZ

          ! Damping due to Stark quadratic
          damp(iz) = sigc*(Atmo%T(iz)**(1d0/6d0))*Atmo%ne(iz) + &
                     damp(iz)

        end do

      end if ! Parameter

      return

      end subroutine broad_stk_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes linear Stark broadening\n
      !!  line(LTEline_class): Structure with the LTE line data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!         damp(dfloat): Inverse lifetime due to linear Stark
      !!                       effect
      subroutine broad_lstk_LTE(line,Atmo,damp)

      ! I/O

      type(LTEline_class):: line
      type(Atmo_class), intent(in):: Atmo
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz

      double precision:: n1, n2, C, sigc, TRfreq

      ! Routine name
      urou = 'broad_lstk_LTE'

      ! If not Hydrogen return
      if (line%ele.ne.1) return

      ! Ionization energy
      TRfreq = 1.0967877174307

      ! Calculate principal quantum number
      n1 = dble(nint(sqrt(ryd/(TRfreq - line%El))))
      n2 = dble(nint(sqrt(ryd/(TRfreq - line%Eu))))

      ! If the principal quantum numbers differ in one
      if (nint(abs(n2 - n1)).eq.1) then

        C = 0.642d0

      ! If differ in something different than 1
      else

        C = 1d0

      end if

      ! Cross section
      sigc = C*0.6d0*(n2*n2 - n1*n1)

      ! For each height
      do iz=1,NZ

        ! Contribution due to linear Stark effect
        damp(iz) = sigc*(Atmo%ne(iz)**(2d0/3d0)) + damp(iz)

      end do

      return

      end subroutine broad_lstk_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      end module broadaux_mod
