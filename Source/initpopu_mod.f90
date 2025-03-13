      !> Populations initialization
      module initpopu_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     11/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/12/2024:    V4.0.0 - Removed any reference to threads in the
!                             call to abortedS (TdPA)
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
!    o Inputpopu_LTE has a dummy input, 'mode', which is suppossed to
!      let the user select partition functions to compute the LTE.
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!    Initpopu_file
!      Initialize the atomic populations from a specified file
!
!    Initpopu
!      Initialize the populations in case they are not read from file
!
!    Initpopu_LTE
!      Initialize the total population and the level population in LTE
!
!    Initpopu_LTE_line
!      Initialize the total population and the level population in LTE
!    for a LTE line
!
!    Initcrhox
!      Initialize the density matrix multipoles
!
!    Initcrho_old
!      Initialize the copy of the density matrix to track the values
!    in the last iteration
!
!    Initpopu_CLE
!      Initialize the atomic populations, the density matrix, and
!    calculate the ionization fraction (if data is given) in the CLE
!    case
!
!    ReviseHatmo
!      Update the hydrogen number density in the model atmosphere
!    from the given hydrogen atom
!
!    ReviseHatom
!      Scale the the hydrogen number density in the model atom from
!    the given model atmosphere
!
!    correctpop
!      Normalize or de-normalize the atomic population
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use initpopuaux_mod
      use parameters_mod , only : cZero , TINYR00
      use types_mod

      ! Big relative change for total hydrogen population
      double precision, parameter:: delta = 0.25d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the atomic populations from a specified file\n
      !!        Atom(Atom_class): Structure with atomic data\n
      !!  filename(character(:)): Name of the file to read\n
      !!        Atmo(Atmo_class): Structure with atmospheric data
      subroutine Initpopu_file(Atom,filename,Atmo)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: filename

      ! Local

      logical:: LJstruct


      !
      ! Check files
      !

      ! If there is no input file specified
      if (trim(filename).eq.'N') then

        ! If the atom is Hydrogen
        if (Atom%Element.eq.' H') then

          ! If there are populations in the atmospheric model
          if (Atmo%typo.eq.0) then

            ! Allocate populations
            allocate(Atom%popu(Atom%nlevel,nz))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%popu)

            ! If it is the one hard-coded
            if (Atom%cust.and.Atmo%typo.eq.0) then

              ! Use the information in the atmospheric file as it is
              Atom%popu = transpose(Atmo%nh(:,:))

            ! If not the hard-coded
            else

              ! Populate the hydrogen atomic model with the data in
              ! the model atmosphere
              call InitHpopu(Atom,Atmo,LJstruct)

              ! If master and there is structure
              if (pid.eq.0.and.LJstruct) then

                ! Warning
                umsg = ' - Using custom H atom without '// &
                       'specifying populations. Assuming '// &
                       'LTE for different L and J within an '//&
                       'n and using populations derived from '// &
                       'the atmosphere.'
                call verbose

              end if ! Master and atomic structure
            end if ! Custom hydrogen
          end if ! There are populations in the atmospheric model
        end if ! Not hydrogen

      ! If there is a file specified
      else

        ! Read the population from the file
        call rPopu(filename,Atom)

        ! Protect the atom from molecules
        Atom%mol_protect = .True.

      end if ! File specified

      ! Check if everything is fine
      call control

      return

      end subroutine Initpopu_file

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the populations in case they are not read from
      !! file\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!   active(logical): Bool to specify if this atom is active
      subroutine Initpopu(Atom,Atmo,active)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), intent(inout):: Atom
      logical, intent(in):: active

      ! Local

      logical:: LJstruct

      integer:: iz,il


      !
      ! Determine initial populations
      !

      ! No existing NLTE population
      if (.not.allocated(Atom%popu)) then

        ! Allocate space for populations
        allocate(Atom%popu(Atom%nlevel,nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%popu)

        ! If there are departure coefficients specified
        if (allocated(Atom%depar)) then

          ! Compute populations from departure coefficients
          Atom%popu = Atom%populte*Atom%depar

        ! There are no departure coefficients either
        else

          ! If the atom is Hydrogen
          if (Atom%Element.eq.' H') then

            ! If it is the one hard-coded
            if (Atom%cust) then

              ! Use the information in the atmospheric file
              Atom%popu = transpose(Atmo%nh(:,:))

            ! If normal model atom for hydrogen
            else

              ! Populate the hydrogen atomic model with the data in
              ! the model atmosphere
              call InitHpopu(Atom,Atmo,LJstruct)

              ! If master and there is structure
              if (pid.eq.0.and.LJstruct) then

                ! Issue warning
                umsg = ' - Using custom H atom without '// &
                       'specifying populations. Assuming '// &
                       'LTE for different L and J within an '//&
                       'n and using populations derived from '// &
                       'the atmosphere.'
                call verbose

              end if ! Master and LJ structure
            end if ! Hard-coded H atom

          ! If it is not Hydrogen
          else

            ! Copy from LTE populations
            Atom%popu = Atom%populte

          end if ! Hydrogen atom
        end if ! Departure coefficients supplied
      end if ! NLTE populations supplied

      ! For every active atom
      if (active) then

        !
        ! Check that there are not zero populations
        !

        ! For each height
        do iz=1,nz

          ! For each level running through the term and J indexes
          do il=1,Atom%nlevel

            ! If too small population
            if (Atom%popu(il,iz).lt.TINYR00) then

              ! Add a very small number to the population
              Atom%popu(il,iz) = Atom%popu(il,iz) + TINYR00
              Atom%n(iz) = Atom%n(iz) + TINYR00

            end if ! Too small population

          end do ! Levels
        end do ! Heights

      end if ! Active atom

      ! Check if everything is fine
      call control

      return

      end subroutine Initpopu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the total population and the level population in
      !! LTE\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!   active(logical): Bool to specify if this atom is activ\n
      !!     mode(integer): Type of LTE calculation
      subroutine Initpopu_LTE(Atom,Atmo,active,mode)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: active
      integer, intent(in):: mode

      ! Local

      integer:: iz


      !
      ! Calculate total population of the atom
      !

      ! If NLTE populations exist
      if (allocated(Atom%popu)) then

        ! For each height
        do iz=1,nz

          ! Total population is sum of the NLTE ones
          Atom%n(iz) = sum(Atom%popu(:,iz))

        end do ! Heights

      ! If LTE populations
      else

        ! Get the population from the abundance and the Hydrogen
        ! population

        ! If active
        if (active) then

          ! For each height
          do iz=1,NZ

            ! Get from atomic model
            Atom%n(iz) = Atom%abun_mod*Atom%abun*Atmo%nht(iz)

          end do ! Heights

        ! If passive
        else

          ! For each height
          do iz=1,NZ

            ! Get from atomic model
            Atom%n(iz) = Atom%abun*Atmo%nht(iz)

          end do ! Heights

        end if ! Active atom

      endif ! NLTE/LTE


      !
      ! Determine LTE population
      !

      ! Allocate and initialize LTE population
      allocate(Atom%populte(Atom%nlevel,nz))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%populte)
      Atom%populte = 0d0

      ! Calculate LTE population
      call LTEPopu(Atom,Atmo)

      ! Check if everything is fine
      call control

      return

      ! Deceive compiler
      iz = mode

      end subroutine Initpopu_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the total population and the level population in
      !! LTE for a LTE line\n
      !!  LTEline(LTEline_class): Structure with LTE line data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data
      subroutine Initpopu_LTE_line(LTEline,Atmo)

      ! I/O

      type(LTEline_class), intent(inout):: LTEline
      type(Atmo_class), intent(in):: Atmo

      ! Local

      integer:: nstg,ii,iim,iz,zm

      double precision:: debey,C0,C1
      double precision, dimension(:), allocatable:: Ei
      double precision, dimension(:,:), allocatable:: pf


      !
      ! Calculate total population of the atom
      !

      ! Allocate arrays
      allocate(LTEline%n(nZ),LTEline%nl(nZ),LTEline%nu(nZ))
      MRAMc = MRAMc + 1d-6*sizeof(LTEline%n)
      MRAMc = MRAMc + 1d-6*sizeof(LTEline%nl)
      MRAMc = MRAMc + 1d-6*sizeof(LTEline%nu)

      ! Initialize total population to hydrogen number density
      LTEline%n = Atmo%nht

      ! Get the partition function for the model atmosphere and
      ! this atom
      call getpf(atom_index2char(LTEline%ele),nstg,pf,Ei,Atmo)

      ! If stage above available
      if (LTEline%stage.gt.nstg) then

        ! Make all populations zero
        LTEline%n = 0d0
        LTEline%nl = 0d0
        LTEline%nu = 0d0

        ! And skip
        return

      end if ! Non-available stage

      !
      ! Calculate Debey correction to the ionization potential
      !

      ! Initialize
      debey = 0d0

      ! Determine the charge of the shell nucleus + rest of
      ! electrons and the change of stages between this level
      ! and the ground level of the model
      zm = LTEline%stage - 1
      iim = zm

      ! Add the contribution to the Debey correction
      do ii=1,iim
        debey = debey + zm
        zm = zm + 1
      end do

      ! Constants
      C0 = hplanck*hplanck/2d0/PI/me/kb
      ! 1d3 sqrt(cm^-3 -> m^-3) electron density below
      C1 = sqrt(8d0*PI/kb)*((qel*qel/pi4eps0)**(1.5d0))*1d3

      ! For each height
      do iz=1,nZ

        ! Get LTE population at this height
        call LTEPopu_line(LTEline,Atmo%T(iz),Atmo%ne(iz),C0,C1, &
                          iz,nstg,Ei,pf(:,iz),debey)

      end do ! Heights

      return

      end subroutine Initpopu_LTE_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the density matrix multipoles\n
      !!  Atom(Atom_class): Structure with atomic data
      subroutine Initcrho(Atom)

      ! I/O

      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iz,it,i,iJ,iR

      double precision:: fac


      !
      ! Allocations
      !

      ! Vector with rhoKQ values
      allocate(Atom%crho(Atom%ndim,nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%crho)

      ! Flags to ignore rhoKQ in the RT rates
      allocate(Atom%rhonull(Atom%ndim,nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rhonull)

      ! Initialize the allocated arrays
      Atom%crho = cZero
      Atom%rhonull = .True.

      ! For each height
      do iz=1,nZ

        ! For each term
        do it=1,Atom%nMulti

          ! For each level within the term
          do iJ=1,Atom%nJ(it)

            ! Get the level index
            i = Atom%irho(it)%irho_ij(iJ)

            ! Get the multipolar index for K=Q=0
            iR = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)

            ! Get the normalizing factor for the level
            fac = 1d0/sqrt(2d0*Atom%rJval(ij,it) + 1d0)/Atom%n(iz)

            ! Determine rho00 from the population
            Atom%crho(iR,iz) = dcmplx(Atom%popu(i,iz)*fac,0d0)

            ! If there is population, un-flag the level
            if (Atom%popu(i,iz).gt.0d0) &
              Atom%rhonull(iR,iz) = .False.

          end do ! J levels
        end do ! Term
      end do ! Heights

      ! Check if everything is fine
      call control

      return

      end subroutine Initcrho

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the copy of the density matrix to track the values
      !! in the last iteration\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!  Atom0(Rhoc_class(:)): Structure to store the density matrix
      !!                        of the previous iteration
      subroutine Initcrho_old(Atom,Atom0)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Rhoc_class), intent(inout):: Atom0


      ! Allocation
      allocate(Atom0%crho(Atom%ndim,Rz0:Rz1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom0%crho)

      ! Initialize to initial density matrix
      Atom0%crho = Atom%crho(:,Rz0:Rz1)

      ! Check if everything is fine
      call control

      return

      end subroutine Initcrho_old

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the atomic populations, the density matrix, and
      !! calculate the ionization fraction (if data is given) in the
      !! CLE case\n
      !!            Atom(Atom_class): Structure with atomic data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!             bion(double(:)): Read ionization fraction data\n
      !!     ion_column_ind(integer): Index of column in buffer for
      !!                              ionization data\n
      !!  ion_column_ind(integer(:)): Index of column in buffer for
      !!                              ionization data\n
      !!   ion_value_ind(integer(:)): Index of value in value array
      !!                              for ionization data\n
      !!      chianti(chianti_class): Structure with the CHIANTI
      !!                              data\n
      !!                 ix(integer): Index location along the LOS\n
      !!             active(logical): Bool to specify if this atom is
      !!                              active
      subroutine Initpopu_CLE(Atom,Atmo,bion,ion_column_ind, &
                              ion_value,ion_value_ind,chianti, &
                              ix,active)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(chianti_class), intent(in):: chianti
      type(Atom_class), intent(inout):: Atom
      logical, intent(in):: active
      integer, intent(in):: ion_column_ind,ion_value_ind,ix
      double precision, dimension(:), intent(in):: bion,ion_value

      ! Local

      integer:: ibion

      double precision:: frc


      !
      ! Calculate total population of the atom
      !

      ! If constant value given
      if (ion_value_ind.gt.0) then

        ! Ion fraction from value data
        frc = ion_value(ion_value_ind)

      ! If value is given for the whole LOS
      else if (ion_column_ind.gt.0) then

        ! Location of ion in this point in the LOS
        ibion = (ion_column_ind - 1)*Atmo%nz + ix

        ! Ion fraction from ionization data
        frc = bion(ibion)

      ! CHIANTI data
      else if (chianti%nT.gt.0) then

        ! Get ion fraction from CHIANTI
        call chianti_fraction(Atom,Atmo%T(1),chianti,frc)

      ! No data
      else

        ! Issue warning
        umsg = 'No ionization data for element '// &
               Atom%Element//'. Ionization fraction '// &
               'set to one'
        urou = 'Initpopu_CLE'
        call abortedS(umsg,urou,.False.,.True.)

        ! Ion fraction unknown
        frc = 1d0

      end if ! Type of ionization data, if any

      ! Atomic population
      Atom%n = Atom%abun*Atom%abun_mod*Atmo%nht*frc

      !
      ! Determine LTE population
      !

      ! Allocate and initialize LTE population
      allocate(Atom%populte(Atom%nlevel,nz))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%populte)
      Atom%populte = 0d0

      ! Calculate LTE population
      call LTEPopu(Atom,Atmo)

      ! Allocate and initialize fake NLTE population to use in
      ! background
      allocate(Atom%popu(Atom%nlevel,nz))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%popu)
      Atom%popu = Atom%populte

      ! For every active atom
      if (active) then

        ! Allocate vector with rhoKQ values
        allocate(Atom%crho(Atom%ndim,nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%crho)

        ! Allocate flags to ignore rhoKQ in the RT rates
        allocate(Atom%rhonull(Atom%ndim,nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%rhonull)

        ! Initialize the allocated arrays
        Atom%crho = cZero
        Atom%rhonull = .True.

      end if ! Active atom

      ! Check if everything is fine
      call control

      return

      end subroutine Initpopu_CLE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Update the hydrogen number density in the model atmosphere
      !! from the given hydrogen atom\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data
      subroutine ReviseHatmo(Atom,Atmo)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(inout):: Atmo

      ! Local


      integer, parameter:: nn = 6

      integer:: iz,it,iJ,il,jl

      double precision:: dd
      double precision, dimension(nn):: Emin, Emax, NH


      !
      ! Definition of a H atom energies
      !              0     1     2      3      4        p+
      Emin = (/ -1.0d0,0.5d0,0.9d0,1.00d0,1.04d0,1.09677d0 /)
      Emax = (/  0.5d0,0.9d0,1.0d0,1.04d0,1.06d0,1.00000d9 /)

      ! For each height
      do iz=1,nz

        ! Initialize level index
        il = 0

        ! Initialize populations
        NH = 0d0

        ! For each term
        do it=1,Atom%nmulti

          ! For each level within the term
          do iJ=1,Atom%nJ(it)

            ! Advance level
            il = il + 1

            ! For each range in the tabulation
            do jl=1,nn

              ! Check if this level is in the current range
              if (Atom%FSfreq(iJ,it).gt.Emin(jl).and. &
                  Atom%FSfreq(iJ,it).lt.Emax(jl)) then

                ! Add to population
                NH(jl) = NH(jl) + Atom%popu(il,iz)

              end if ! Check coincidence

            end do ! Levels (atmo)
          end do ! FS level (atom)
        end do ! Term (atom)

        ! At first point
        if (iz.eq.1) then

          ! If no ground level population
          if (NH(1).le.0d0) then

            ! If master send message
            if (pid.eq.0) then
              umsg = ' - You provided H populations, but there '// &
                     'were no populations in the ground level'
              call verbose
            end if

            ! Exit
            return

          end if ! No ground level population

          ! If no proton population
          if (NH(6).le.0d0) then

            ! If master send message
            if (pid.eq.0) then
              umsg = ' - You provided H populations, but there '// &
                     'were no populations for protons'
              call verbose
            end if

            ! Exit
            return

          end if ! No ground level population
        end if ! First point

        ! Update atmosphere densities
        Atmo%nh(iz,:) = NH(:)

        ! Master checks if difference is significative with previous
        ! populations (total)
        if (pid.eq.0) then

          ! Compute relative difference
          dd = abs(Atmo%nha(iz) - sum(NH))/ &
                  (Atmo%nha(iz) + sum(NH))

          ! If relative difference > parameter in header
          if (dd.gt.delta.and.Atmo%nha(iz).gt.0d0) then

            ! If optical depth scale
            if (ztau) then

              ! Issue warning
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,f4.2,1x,A,1x,f4.2,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: After updating the H populations in '// &
                'the atmosphere at tau',iz,Atmo%z(iz), &
                ', the relative change',dd,'is above',delta,'.', &
                Atmo%nha(iz),'==>',sum(NH)
                call verbose

            ! If geometrical height
            else

              ! Issue warning
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,f4.2,1x,A,1x,f4.2,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: After updating the H populations in '// &
                'the atmosphere at height',iz,Atmo%z(iz)*1e-5, &
                'km, the relative change',dd,'is above',delta,'.', &
                Atmo%nha(iz),'==>',sum(NH)
                call verbose

            end if ! Type of height scale
          end if ! Significant change in populations
        end if ! Master

        ! Save total atomic population
        Atmo%nha(iz) = sum(NH)

        ! If less total hydrogen than atomic hydrogen
        if (Atmo%nht(iz).lt.Atmo%nha(iz)) then

          ! Compute relative difference
          dd = abs(Atmo%nha(iz) - Atmo%nht(iz))/ &
                  (Atmo%nha(iz) + Atmo%nht(iz))

          ! Master checks if difference is significative
          if (pid.eq.0.and.Atmo%nht(iz).gt.0d0.and.dd.gt.delta) then

            ! If optical depth scale
            if (ztau) then

              ! Issue warning
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: Input atomic hydrogen larger than '// &
                'previous total hydrogen. Making them equal at '// &
                'tau',iz,Atmo%z(iz),':',Atmo%nht(iz),'==>', &
                Atmo%nha(iz)
              call verbose

            ! If geometrical scale
            else

              ! Issue warning
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: Input atomic hydrogen larger than '// &
                'previous total hydrogen. Making them equal at '// &
                'height',iz,Atmo%z(iz)*1d-5,':',Atmo%nht(iz), &
                '==>',Atmo%nha(iz)
              call verbose

            end if ! Type of scale
          end if ! Master

          ! Make total equal to atomic
          Atmo%nht(iz) = Atmo%nha(iz)

        end if ! If more atomic than total

      end do ! Heights

      ! Check if everything is fine
      call control

      return

      end subroutine ReviseHatmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Scale the the hydrogen number density in the model atom from
      !! the given model atmosphere\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data
      subroutine ReviseHatom(Atom,Atmo)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iz

      double precision:: N0,N1,dd


      ! For each height
      do iz=1,nz

        ! Get total population in atom
        N0 = sum(Atom%popu(:,iz))

        ! Scale population to the total atomic one
        Atom%popu(:,iz) = Atom%popu(:,iz)*Atmo%nHa(iz)/N0

        ! New sum
        N1 = sum(Atom%popu(:,iz))

        ! Master
        if (pid.eq.0) then

          ! Compute relative difference
          dd = abs(N0 - N1)/(N0 + N1)

          ! If relative difference > parameter in header
          if (dd.gt.delta.and.Atmo%nha(iz).gt.0d0) then

            ! Issue warning
            write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                       '1x,f4.2,1x,A,1x,f4.2,A,'// &
                       '1x,es13.6,1x,A,es13.6)') &
              ' # Warning: After updating the H populations in '// &
              'the atom at height',iz,Atmo%z(iz)*1e-5, &
              'km, the relative change',dd,'is above',delta,'.', &
              N0,'==>',N1
            call verbose

          end if ! Big relative difference
        end if ! Master

      end do ! Heights

      ! Check if everything is fine
      call control

      return

      end subroutine ReviseHatom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalize or de-normalize the atomic population\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!      dir(integer): Direction of correction
      subroutine correctpop(Atom,dir)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      integer, intent(in):: dir

      ! Local

      integer:: iz


      ! Normalize
      if (dir.eq.0) then

        ! For each height
        do iz=Rz0,Rz1

          ! Divide by total
          Atom%popu(:,iz) = Atom%popu(:,iz)/Atom%n(iz)

        end do ! Heights

      ! De-normalize
      else if (dir.eq.1) then

        ! For each height
        do iz=Rz0,Rz1

          ! Multiply by total
          Atom%popu(:,iz) = Atom%popu(:,iz)*Atom%n(iz)

        end do ! Heights

      ! Error, wrong argument
      else

        ! Issue error
        umsg = 'Direction not recognized'
        urou = 'correctpop'
        call abortedS(umsg,urou,.True.,.True.)

      end if ! Transformation direction

      ! Check if everything is fine
      call control

      return

      end subroutine correctpop

!#####################################################################
!#####################################################################
!#####################################################################

      end module initpopu_mod
