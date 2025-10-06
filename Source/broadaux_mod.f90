      !> Contributions to the line broadening
      module broadaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     03/10/2025 V4.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/10/2025:    V4.0.4 - Barklem Van der Waals broadenings can
!                             proceed without a continuum level, they
!                             just neglect the helium contribution
!                             with Unsold (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    The Helium abundance is taken from the hard-coded table, does not
!  depend on the input atoms
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
!  broad_vdw
!    Compute Van der Waals broadening
!
!  broad_stk
!    Compute quadratic Stark broadening
!
!  broad_lstk
!    Compute linear Stark broadening
!
!  broad_vdw_LTE
!    Compute Van der Waals broadening for a LTE line
!
!  broad_stk_LTE
!    Compute quadratic Stark broadening for a LTE line
!
!  broad_lstk_LTE
!    Compute linear Stark broadening for a LTE line
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

      !> Compute Van der Waals broadening\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!   iterml(integer): Lower term of transition\n
      !!   itermu(integer): Upper term of transition\n
      !!    itran(integer): Index of transition\n
      !!      damp(double): Inverse lifetime\n
      !!   aparam(logical): If the equivalant parametric parameters
      !!                    need to be stored
      subroutine broad_vdw(Atom,Atmo,iterml,itermu,itran,damp,aparam)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: aparam
      integer, intent(in):: iterml,itermu,itran
      double precision, dimension(:), intent(inout):: damp

      ! Local

      logical:: get_C6

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

      ! Correct Rydberg energy for mass shift and calculate
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

      !
      ! Find the next continuum
      !

      ! For each term above the upper term
      do i1=itermu+1,Atom%nMulti

        ! Check if stage change
        if (Atom%stage(i1).gt.Atom%stage(itermu)) then
          itermc = i1
          exit
        end if

      end do ! Terms above upper term

      ! Do we need to get C6? If not parametric
      get_C6 = Atom%broad_type(itran).ne.2.and. &
               Atom%broad_type(itran).ne.4

      ! Initialize C6
      C6 = 0d0

      ! If we did not find the continuum
      if (itermc.lt.0.and.Atom%broad_type(itran).ne.2.and. &
          Atom%broad_type(itran).ne.4) then

        ! If Unsold
        if (Atom%broad_type(itran).eq.1) then

          ! Notify error and return
          umsg = 'Could not find continuum in atom '// &
                 Atom%Element//' and Van der Waals '// &
                 'broadening is Unsold'
          call aborted
          return

        ! Barklem
        else

          ! Notify issue
          umsg = 'Could not find continuum in atom '// &
                 Atom%Element//' and Van der Waals '// &
                 'broadening is Barklem, neglecting '// &
                 'helium contribution'
          call abortedS(umsg,urou,.False.,.True.)

          ! Do not get C6
          get_C6 = .False.

        end if ! Type of broadening
      end if ! Continuum not found


      !
      ! Calculate the Van der Waals cross section C6
      !

      ! If C6 needed
      if (get_C6) then

        ! Energy part
        d1 = ryd*ryd/(Atom%TRfreq(itermc) - &
                       Atom%TRfreq(itermu))/ &
                     (Atom%TRfreq(itermc) - &
                       Atom%TRfreq(itermu))  - &
             ryd*ryd/(Atom%TRfreq(itermc) - &
                       Atom%TRfreq(iterml))/ &
                     (Atom%TRfreq(itermc) - &
                       Atom%TRfreq(iterml))

        ! Next ion charge
        Z = Atom%stage(iterml)

        ! VdW cross section. 1d6 to receive cgs populations
        C6 = 8.08d0*((5d0*qel*qel*a4pieps0*pi* &
             Z*Z*rb*rb*d1/pi4eps0/hplanck)**.4d0)*1d6

      end if ! C6 needed


      !
      ! Calculate the broadening
      !

      !
      ! If the type of broadening is Barklem
      if (Atom%broad_type(itran).eq.0.or. &
          Atom%broad_type(itran).eq.3) then

        ! Sigma and alpha calculated in rBarklem
        sigma = Atom%broad_args(1,itran)
        alpha = Atom%broad_args(2,itran)

        ! Broadening constant part. 1d6 for cgs populations.
        sigc = 2d0*rb*rb*(4d0/pi)**(0.5d0*alpha)* &
               GAMMA(0.5d0*(4d0 - alpha))*1d4*sigma* &
               (v02c/1d8)**(.5d0*(1d0 - alpha))*1d6

        ! If storing parametric equivalence
        if (aparam) then
          b1 = 0.5*(1d0 - alpha)
          a1 = 1d8*sigc/((1d0 + mhm/Atom%rmass)**b1)
          b2 = .3d0
          a2 = 1d9*C6*(v02che**.3)/((1d0 + mhem/Atom%rmass)**b2)
          if (Atom%broad_type(itran).eq.0) then
            write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
                  err=1100) itran,itermu,iterml, &
                            '    barklem',a1,b1,a2,b2
          else
            write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
                  err=1100) itran,itermu,iterml, &
                            '      cross',a1,b1,a2,b2
          end if
        end if

        ! For each height
        do iz=1,nz

          ! Barklem for H
          damp(iz) = sigc*(Atmo%T(iz)**(.5d0*(1d0 - alpha)))* &
                     Atmo%nh(iz,1) + damp(iz)

          ! Helium population from the atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then

            ! Get population
            nhe = Atmo%nhe(iz,1)

          ! Helium population from the abundance
          else

            ! Get population
            nhe = Atmo%nh(iz,1)*Ahe

          end if

          ! Unsold for Helium
          damp(iz) = C6*(v02che**.3d0)*nhe*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

      end if ! Barklem inputs

      !
      ! If the type of broadening is Unsold
      if (Atom%broad_type(itran).eq.1) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then

            ! Get population
            nhe = Atmo%nhe(iz,1)

          ! Helium ground level pop from abundance
          else

            ! Get population
            nhe = Atmo%nh(iz,1)*Ahe

          end if

          ! Unsold broadening
          damp(iz) = C6*((v02c**.3d0)*Atmo%nh(iz,1)*args(1) + &
                    (v02che**.3d0)*nhe*args(3))*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

        ! If storing parametric equivalence
        if (aparam) then
          b1 = .3d0
          a1 = 1d8*C6*(v02c**.3)/((1d0 + mhm/Atom%rmass)**b1)
          b2 = .3d0
          a2 = 1d9*C6*(v02che**.3)/((1d0 + mhem/Atom%rmass)**b2)
          write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
                err=1100) itran,itermu,iterml, &
                          '     unsold',a1,b1,a2,b2
        end if

      end if ! Unsold broadening

      !
      ! If the type of broadening is Parametric
      if (Atom%broad_type(itran).eq.2) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then

            ! Get population
            nhe = Atmo%nhe(iz,1)

          ! Helium ground level pop from abundance
          else

            ! Get population
            nhe = Atmo%nh(iz,1)*Ahe

          end if

          ! Parametric broadening
          damp(iz) = 1d-8*args(1)*Atmo%nh(iz,1)* &
        (((1d0 + mhm/Atom%rmass)*Atmo%T(iz))**args(2)) + &
                     1d-9*args(3)*nhe* &
        (((1d0 + mhem/Atom%rmass)*Atmo%T(iz))**args(4)) + &
                     damp(iz)

        end do ! Heights

        ! If storing parametric equivalence
        if (aparam) then
          b1 = args(2)
          a1 = args(1)
          b2 = args(4)
          a2 = args(3)
          write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
                err=1100) itran,itermu,iterml, &
                          ' parametric',a1,b1,a2,b2
        end if

      end if ! Parametric broadening


      !
      ! If the type of broadening is from Kurucz
      if (Atom%broad_type(itran).eq.4) then

        ! Constant in input
        sigc = 10d0**Atom%broad_args(1,itran)

        ! If storing parametric equivalence
        if (aparam) then
          b1 = 0d0
          a1 = 1d8*sigc
          b2 = 0d0
          a2 = 0d0
          write(600,'(i10,1x,i10,"-->",2x,i10,1x,A,4(1x,es15.8))', &
                err=1100) itran,itermu,iterml, &
                          '    Kurucz',a1,b1,a2,b2
        end if

        ! Barklem for H
        damp(:) = sigc*Atmo%nh(:,1) + damp(:)

      end if ! Kurucz inputs

      return

1100  umsg = 'Error writing avdwparam file'
      close(600)
      close(700)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine broad_vdw

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute quadratic Stark broadening\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!   iterml(integer): Lower term of transition\n
      !!   itermu(integer): Upper term of transition\n
      !!    itran(integer): Index of transition\n
      !!      damp(dfloat): Inverse lifetime\n
      !!   aparam(logical): If the equivalant parametric parameters
      !!                    need to be stored
      subroutine broad_stk(Atom,Atmo,iterml,itermu,itran,damp,aparam)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: aparam
      integer, intent(in):: iterml,itermu,itran
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz,itermc

      double precision:: stk,aryd,Z,n1,n2,v02c,C4,sigc


      ! Routine name
      urou = 'broad_stk'

      ! Take the arguments in the input
      stk = Atom%broad_stark(itran)

      ! If the parameter is negative
      if (stk.lt.0d0) then

        ! Constant broadening given by the parameter
        stk = abs(stk)

        ! For each height
        do iz=1,NZ

          ! Proportional to electron density
          damp(iz) = stk*Atmo%ne(iz) + damp(iz)

        end do ! Heights

        ! Store parametric equivalence if requested
        if (aparam) write(700,'(i10,1x,es15.8)',err=1100) itran,stk

      ! Any positive number means non-constant
      else if (stk.gt.0d0) then

        ! Correct Rydberg energy for mass shift and calculate
        ! relative velocities for electrons and ions
        aryd = ryd/(1d0 + mem/Atom%rmass)
        v02c = ((8d0*kb/pi/Atom%rmass/amu)**(1d0/6d0))* &
               (((1d0 + Atom%rmass/mem)**(1d0/6d0)) + &
                ((1d0 + Atom%rmass/armass)**(1d0/6d0)))

        ! Initialize the index of the continuum so it can be used as
        ! flag
        itermc = -1

        !
        ! Find the next continuum
        !

        ! For each term above the upper term
        do iz=itermu+1,Atom%nMulti

          ! Check if stage change
          if (Atom%stage(iz).gt.Atom%stage(itermu)) then
            itermc = iz
            exit
          end if

        end do ! Terms above upper term

        ! If we did not find the continuum
        if (itermc.lt.0) then

          ! Master notify the issue
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

        ! We found the continuum level
        else

          !
          ! Calculate the Stark cross section C4
          !

          ! Charge of next ion
          Z = Atom%stage(iterml)

          ! Effective principal quantum numbers
          n1 = Z*sqrt(aryd/(Atom%TRfreq(itermc) - &
                            Atom%TRfreq(iterml)))
          n2 = Z*sqrt(aryd/(Atom%TRfreq(itermc) - &
                            Atom%TRfreq(itermu)))

          ! C4 constant. 1d6 to receive cgs populations
          C4 = qel*qel*2d0*pi*rb*rb*rb* &
               ((n2*(5d0*n2*n2 + 1d0))**2d0 - &
                (n1*(5d0*n1*n1 + 1d0))**2d0)*stk/ &
               pi4eps0/18d0/Z/Z/Z/Z/hplanck

          ! Cross section
          sigc = 11.37d0*(C4**(2d0/3d0))*v02c*1d6

          ! Store parametric equivalence if requested
          if (aparam) &
            write(700,'(i10,1x,i10,"-->",2x,i10,1x,es15.8)', &
                  err=1100) itran,itermu,iterml,sigc

          ! For each height
          do iz=1,NZ

            ! Damping due to Stark quadratic
            damp(iz) = sigc*(Atmo%T(iz)**(1d0/6d0))*Atmo%ne(iz) + &
                       damp(iz)

          end do ! Heights

        end if ! Found a continuum

      end if ! Type of input

      return

1100  umsg = 'Error writing astkparam file'
      close(600)
      close(700)
      call abortedS(umsg,urou,.True.,.True.)
      call control

      end subroutine broad_stk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute linear Stark broadening\n
      !!     Atom(Atom_class): Structure with atomic data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!      iterml(integer): Lower term of transition\n
      !!      itermu(integer): Upper term of transition\n
      !!         damp(dfloat): Inverse lifetime
      subroutine broad_lstk(Atom,Atmo,iterml,itermu,damp)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      integer, intent(in):: iterml,itermu
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz, itermc

      double precision:: n1, n2, C, sigc, TRfreq


      ! If not Hydrogen return
      if (Atom%Element.ne.' H') return

      ! Initialize the index of the continuum so it can be used as
      ! flag
      itermc = -1

      !
      ! Find the next continuum
      !

      ! For each term above the upper term
      do iz=itermu+1,Atom%nMulti

        ! Check if stage change
        if (Atom%stage(iz).gt.Atom%stage(itermu)) then
          itermc = iz
          exit
        end if

      end do ! Terms above upper term

      ! If we did not find the continuum
      if (itermc.lt.0) then

        ! Master notify the issue
        if (pid.eq.0) then
          umsg = ' # Could not find continuum '// &
                 'in atom H, using hardcoded '// &
                 'ionization energy.'
          call verbose
        end if

        ! Use knowable ionization energy
        TRfreq = 1.0967877174307

      ! Found the continuum
      else

        ! Use the model's ionization energy
        TRfreq = Atom%TRfreq(itermc)

      end if ! Found continuum

      ! Calculate principal quantum number
      n1 = dble(nint(sqrt(ryd/ &
                (TRfreq - Atom%TRfreq(iterml)))))
      n2 = dble(nint(sqrt(ryd/ &
                (TRfreq - Atom%TRfreq(itermu)))))

      ! If the principal quantum numbers differ in one
      if (nint(abs(n2 - n1)).eq.1) then

        ! Constant
        C = 0.642d0

      ! If differ in something different than 1
      else

        ! Constant
        C = 1d0

      end if ! Distance between principal quantum numbers

      ! Cross section
      sigc = C*0.6d0*(n2*n2 - n1*n1)

      ! For each height
      do iz=1,NZ

        ! Contribution due to linear Stark effect
        damp(iz) = sigc*(Atmo%ne(iz)**(2d0/3d0)) + damp(iz)

      end do ! Heights

      return

      end subroutine broad_lstk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute Van der Waals broadening for a LTE line\n
      !!  line(LTEline_class): Structure with the LTE line data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!         damp(dfloat): Inverse lifetime
      subroutine broad_vdw_LTE(line,Atmo,damp)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Atmo_class), intent(in):: Atmo
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz

      double precision:: aryd,Z,sigma,alpha,Ahe,Ei
      double precision:: v02c,v02che,nhe,d1,C6,sigc
      double precision, dimension(4):: args


      ! Take the arguments in the input
      args = line%broad_args

      ! If Kurucz, easy
      if (line%broad_type.eq.4) then
        damp = args(1)*Atmo%nh(:,1)
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
      ! Get ionization energy from atmospheric data
      !
      Ei = Atmo%ele(line%ele)%Ei(line%stage)

      ! If not parametric or Kurucz
      if (line%broad_type.ne.2.and. &
          line%broad_type.ne.4) then

        !
        ! Calculate the Van der Waals cross section C6
        !

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

      ! If the type of broadening is Barklem
      if (line%broad_type.eq.0.or. &
          line%broad_type.eq.3) then

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

            ! Get population
            nhe = Atmo%nhe(iz,1)

          ! Helium population from the abundance
          else

            ! Get population
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

            ! Get population
            nhe = Atmo%nhe(iz,1)

          ! Helium ground level pop from abundance
          else

            ! Get population
            nhe = Atmo%nh(iz,1)*Ahe

          end if

          ! Unsold broadening
          damp(iz) = C6*((v02c**.3d0)*Atmo%nh(iz,1)*args(1) + &
                    (v02che**.3d0)*nhe*args(3))*(Atmo%T(iz)**.3d0) + &
                     damp(iz)

        end do ! Heights

      end if ! Unsold broadening

      ! If the type of broadening is Parametric
      if (line%broad_type.eq.2) then

        ! For each height
        do iz=1,NZ

          ! Helium ground level pop from atmosphere
          if (Atmo%nhe(1,1).ge.0d0) then

            ! Get population
            nhe = Atmo%nhe(iz,1)

          ! Helium ground level pop from abundance
          else

            ! Get population
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

      !> Compute quadratic Stark broadening for a LTE line\n
      !!  line(LTEline_class): Structure with the LTE line data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!         damp(dfloat): Inverse lifetime due to quadratic Stark
      !!                       effect
      subroutine broad_stk_LTE(line,Atmo,damp)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Atmo_class), intent(in):: Atmo
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz

      double precision:: stk,aryd,Z,n1,n2,v02c,C4,sigc,Ei


      ! Take the arguments in the input
      stk = line%broad_stark

      ! Constant negative broadening parameter
      if (stk.lt.0d0) then

       !stk = abs(stk)*1d6
        stk = abs(stk)

        ! For each height
        do iz=1,NZ

          ! Parameter in input
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

        end do ! Heights

      end if ! Parameter

      return

      end subroutine broad_stk_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute linear Stark broadening for a LTE line\n
      !!  line(LTEline_class): Structure with the LTE line data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!         damp(dfloat): Inverse lifetime due to linear Stark
      !!                       effect
      subroutine broad_lstk_LTE(line,Atmo,damp)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Atmo_class), intent(in):: Atmo
      double precision, dimension(:), intent(inout):: damp

      ! Local

      integer:: iz

      double precision:: n1, n2, C, sigc, TRfreq


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

      end do ! Heights

      return

      end subroutine broad_lstk_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      end module broadaux_mod
