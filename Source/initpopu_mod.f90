      !> Populations initialization
      module initpopu_mod
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
!     08/07/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.5 - Added Initpopu_LTE_line (TdPA)
!
!     07/03/2023:    V3.0.4 - Split the initialoization of populations
!                             into two suboutines. The new one is
!                             Initpopu_LTE (TdPA)
!
!     11/24/2022:    V3.0.3 - Finished Initpopu_CLE by adding the
!                             use of CHIANTI data (TdPA)
!
!     10/26/2022:    V3.0.2 - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.1 - Moved the initialization of the density
!                             matrix from Initpopu to the new Initcrho
!                             routine, changing the arguments of
!                             Initpopu to adjust for this (TdPA)
!                           - Added Initpopu_CLE routine, to
!                             initialize the population variables in
!                             the CLE case (TdPA)
!                           - Implemented the limitation of the height
!                             axis (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     09/11/2020:    V1.2.5 - Added the option to initialize with
!                             departure coefficients (TdPA)
!
!     03/05/2020:    V1.2.4 - Added Initpopu_file, which loads the
!                             provided populations files (TdPA)
!                           - Added Initpopu does not load files
!                             anymore (TdPA)
!                           - Initpopu calls a routine to initialize
!                             the H populations from the atmosphere
!                             if no file is specified and the model
!                             for H is not the hardcoded (TdPA)
!                           - Added a check for zero populations to
!                             avoid them (TdPA)
!                           - The populations are not normalized in
!                             Initpopu anymore (TdPA)
!                           - Added ReviseHatom routine to scale the
!                             populations to the hydrogen in the
!                             atmosphere (TdPA)
!                           - Added correctpop to normalize or revert
!                             the normalization of atomic populations
!                             of active atoms (TdPA)
!
!     02/10/2020:    V1.2.3 - Added ReviseHatmo routine, to revise the
!                             populations of the atmospheric file and
!                             overwrite them with the ones in the
!                             input file for H (TdPA)
!
!     11/19/2019:    V1.2.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     11/12/2019:    V1.2.1 - The variable Atom%popu now stores the
!                             populations normalized by the total
!                             population of the atom (TdPA)
!
!     09/26/2019:    V1.2.0 - The collisions that were computed here
!                             are now computed elsewhere, both the
!                             ones in the Initpopu routine and the
!                             full setNSCcoeff routine, that does not
!                             exist anymore (TdPA)
!
!     09/24/2019:    V1.1.3 - The total population of the atom, which
!                             used to be computed in abund, is now
!                             computed here (TdPA)
!
!     08/14/2019:    V1.1.2 - Added routine setNSCcoeff (TdPA)
!
!     03/19/2019:    V1.1.1 - Avoid dividing by zero (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     04/19/2017:    V1.0.0 - First version (TdPA)
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
!    This subroutine initializes the populations
!
!    Initpopu_file:
!      Initializes atomic populations from file
!
!    Initpopu:
!      Initializes the populations if there is no file
!
!    Initpopu_LTE:
!      Computes LTE populations
!
!    Initpopu_LTE_line:
!      Computes LTE populations for an LTE line
!
!    Initcrho:
!      Initializes rhoKQ
!
!    Initpopu_CLE:
!      Computes LTE populations, initializes the populations,
!    initializes rhoKQ in the CLE case and compute ionization
!    fraction if data is given
!
!    ReviseHatmo:
!      Reads H population and updates the atmospheric densities
!
!    ReviseHatom:
!      Scales atomic H populations to the ones given in the atmosphere
!
!    correctpop:
!      Normalizes or de-normalizes the atomic populations
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

      double precision, parameter:: delta = 0.25d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initializes the atomic populations from file\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!   filename(character(:)): Name of the file to read, if any\n
      !!         Atmo(Atmo_class): Structure with atmospheric data
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

            ! If it is the one hardwired, use the information in the
            ! atmospheric file as it is
            if (Atom%cust.and.Atmo%typo.eq.0) then

              Atom%popu = transpose(Atmo%nh(:,:))

            ! If not the hardwired, call the routine to do it in a
            ! controlled manner
            else

              call InitHpopu(Atom,Atmo,LJstruct)

              if (pid.eq.0.and.LJstruct) then
                umsg = ' - Using custom H atom without '// &
                       'specifying populations. Assuming '// &
                       'LTE for different L and J within an '//&
                       'n and using populations derived from '// &
                       'the atmosphere.'
                call verbose
              end if

            end if ! Custom hydrogen
          end if ! There are populations in the atmospheric model
        end if ! Not hydrogen

      ! If there is a file specified
      else

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

      !> Initializes the atomic populations\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!          active(logical): Bool to specify if this atom is
      !!                           active or not
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

      ! No NLTE population
      if (.not.allocated(Atom%popu)) then

        allocate(Atom%popu(Atom%nlevel,nz))

        ! If there are departure coefficients
        if (allocated(Atom%depar)) then

          Atom%popu = Atom%populte*Atom%depar

        else

          ! If the atom is Hydrogen
          if (Atom%Element.eq.' H') then

            ! If it is the one hardwired, use the information in the
            ! atmospheric file
            if (Atom%cust) then

              Atom%popu = transpose(Atmo%nh(:,:))

            ! If the user introduced it, assume LTE
            else

              call InitHpopu(Atom,Atmo,LJstruct)

              if (pid.eq.0.and.LJstruct) then
                umsg = ' - Using custom H atom without '// &
                       'specifying populations. Assuming '// &
                       'LTE for different L and J within an '//&
                       'n and using populations derived from '// &
                       'the atmosphere.'
                call verbose
              end if

            end if ! Hardwired H atom

          ! If it is not Hydrogen, use LTE
          else

            Atom%popu = Atom%populte

          end if ! H atom
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
              Atom%popu(il,iz) = Atom%popu(il,iz) + TINYR00
              Atom%n(iz) = Atom%n(iz) + TINYR00
            end if

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

      !> Initializes the atomic populations\n
      !!   Atom(Atom_class): Structure with the atomic data\n
      !!   Atmo(Atmo_class): Structure with atmospheric data\n
      !!    active(logical): Bool to specify if this atom is active
      !!                     or not\n
      !!      mode(integer): Type of LTE
      subroutine Initpopu_LTE(Atom,Atmo,active,mode)

      ! I/O
      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), intent(inout):: Atom
      logical, intent(in):: active
      integer, intent(in):: mode

      ! Local
      integer:: iz


      !
      ! Calculate total population of the atom
      !

      ! If NLTE populations
      if (allocated(Atom%popu)) then

        ! For each height, total population
        do iz=1,nz
          Atom%n(iz) = sum(Atom%popu(:,iz))
        end do

      ! If LTE populations
      else

        ! Get the population from the abundance and the Hydrogen
        ! population

        ! If active
        if (active) then
          do iz=1,NZ
            Atom%n(iz) = Atom%abun_mod*Atom%abun*Atmo%nht(iz)
          end do
        ! If passive
        else
          do iz=1,NZ
            Atom%n(iz) = Atom%abun*Atmo%nht(iz)
          end do
        end if

      endif ! NLTE/LTE


      !
      ! Determine LTE population
      !

      ! LTE population
      allocate(Atom%populte(Atom%nlevel,nz))
      Atom%populte = 0d0

      ! Calculate LTE population
      call LTEPopu(Atom,Atmo)

      ! Check if everything is fine
      call control

      return

      end subroutine Initpopu_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get LTE populations for a LTE spectral line\n
      !!  LTEline(Input_class): Structure with LTE line data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data
      subroutine Initpopu_LTE_line(LTEline,Atmo)

      ! I/O
      type(LTEline_class):: LTEline
      type(Atmo_class):: Atmo

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

      ! Initialize total population to H
      LTEline%n = Atmo%nht

      ! Get the partition function for the model atmosphere and
      ! this atom
      call getpf(atom_index2char(LTEline%ele), &
                 nstg,pf,Ei,Atmo)

      ! If stage above available, skip
      if (LTEline%stage.gt.nstg) then
        LTEline%n = 0d0
        LTEline%nl = 0d0
        LTEline%nu = 0d0
        return
      end if

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

        ! Get LTE populations
        call LTEPopu_line(LTEline,Atmo%T(iz),Atmo%ne(iz),C0,C1, &
                          iz,nstg,Ei,pf(:,iz),debey)

      end do

      return

      end subroutine Initpopu_LTE_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initializes the density matrix\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atom0(Atom_class): A copy of Atom
      subroutine Initcrho(Atom,Atom0)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Rhoc_class), intent(inout):: Atom0

      ! Local

      integer:: iz,it,i,iJ,iR
      double precision:: fac


      !
      ! Determine initial rho00
      !

      ! Allocations
      ! Vector with rhoKQ values
      allocate(Atom%crho(Atom%ndim,Rz0:Rz1))
      ! Vector with rhoKQ values from the last iteration
      allocate(Atom0%crho(Atom%ndim,Rz0:Rz1))
      ! Flags to ignore rhoKQ in the RT rates
      allocate(Atom%rhonull(Atom%ndim,Rz0:Rz1))

      ! Initialize the allocated arrays
      Atom%crho = cZero
      Atom0%crho = cZero
      Atom%rhonull = .True.

      !
      ! Check that there are not zero
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each level running through the term and J indexes
        do it=1,Atom%nMulti
          do iJ=1,Atom%nJ(it)

            ! Get the level index
            i = Atom%irho(it)%irho_ij(iJ)
            ! Get the multipolar index
            iR = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)
            ! Get the normalizing factor for the level
            fac = 1d0/sqrt(2d0*Atom%rJval(ij,it) + 1d0)/Atom%n(iz)

            ! Determine rho00
            Atom%crho(iR,iz) = dcmplx(Atom%popu(i,iz)*fac,0d0)

            ! If there is population, un-flag the level
            if (Atom%popu(i,iz).gt.0d0) &
              Atom%rhonull(iR,iz) = .False.

          end do ! J levels
        end do ! Term
      end do ! Heights

      ! The old populations are the same
      Atom0%crho = Atom%crho

      ! Check if everything is fine
      call control

      return

      end subroutine Initcrho

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initializes the atomic populations in the CLE case\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!         bion(double(:)): Read ionization fraction data\n
      !! ion_column_ind(integer): Index of column in buffer for
      !!                          ionization data\n
      !!    ion_value(double(:)): Numeric constant ionization
      !!                             fraction values\n
      !!  ion_value_ind(integer): Index of value in value array for
      !!                             ionization data\n
      !!  chianti(chianti_class): Structure with the CHIANTI data\n
      !!             ix(integer): Index location along the LOS\n
      !!         active(logical): Bool to specify if this atom is
      !!                          active or not
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

      ! Try value
      if (ion_value_ind.gt.0) then

        ! Ion fraction
        frc = ion_value(ion_value_ind)

      ! Try column
      else if (ion_column_ind.gt.0) then

        ! Location of ion
        ibion = (ion_column_ind - 1)*Atmo%nz + ix

        ! Ion fraction
        frc = bion(ibion)

      ! CHIANTI data
      else if (chianti%nT.gt.0) then

        ! Get ion fraction from CHIANTI
        call chianti_fraction(Atom,Atmo%T(1),chianti,frc)

      ! No data
      else

        umsg = 'No ionization data for element '// &
               Atom%Element//'. Ionization fraction '// &
               'set to one'
        urou = 'Initpopu_CLE'
        call abortedS(umsg,urou,-1,.False.,.True.)

        ! Ion fraction
        frc = 1d0

      end if ! Type of ionization data, if any

      ! Atomic population
      Atom%n = Atom%abun*Atom%abun_mod*Atmo%nht*frc

      !
      ! Determine LTE population
      !

      ! LTE population
      allocate(Atom%populte(Atom%nlevel,nz))
      Atom%populte = 0d0

      ! Calculate LTE population
      call LTEPopu(Atom,Atmo)

      ! Initialize fake NLTE population to use in background
      allocate(Atom%popu(Atom%nlevel,nz))
      Atom%popu = Atom%populte

      ! For every active atom
      if (active) then

        !
        ! Determine initial rho00
        !

        ! Allocations
        ! Vector with rhoKQ values
        allocate(Atom%crho(Atom%ndim,nz))
        ! Flags to ignore rhoKQ in the RT rates
        allocate(Atom%rhonull(Atom%ndim,nz))

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

      !> Manage hydrogen populations. Set Atmos from Atom file\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data
      subroutine ReviseHatmo(Atom,Atmo)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), intent(inout):: Atom

      ! Parameter
      integer, parameter:: nn = 6

      ! Local

      integer:: iz,it,iJ,il,jl
      double precision:: dd
      double precision, dimension(nn):: Emin, Emax, NH


      !
      ! Definitions of an H atom
      !              0     1     2      3      4        p+
      Emin = (/ -1.0d0,0.5d0,0.9d0,1.00d0,1.04d0,1.09677d0 /)
      Emax = (/  0.5d0,0.9d0,1.0d0,1.04d0,1.06d0,1.00000d9 /)


      !
      ! Try to look for the levels of each n level
      !

      ! For each height
      do iz=1,nz

        ! Advance level index
        il = 0

        ! Initialize populations
        NH = 0d0

        ! For each level in the model
        do it=1,Atom%nmulti
          do iJ=1,Atom%nJ(it)

            ! Advance level
            il = il + 1

            ! For each level in the atmosphere
            do jl=1,nn

              ! Check if this term
              if (Atom%FSfreq(iJ,it).gt.Emin(jl).and. &
                  Atom%FSfreq(iJ,it).lt.Emax(jl)) then

                ! Add to population
                NH(jl) = NH(jl) + Atom%popu(il,iz)

              end if ! Check coincidence

            end do ! Levels (atmo)
          end do ! FS level (atom)
        end do ! Term (atom)

        ! At first point, check there is ground term and protons
        if (iz.eq.1) then

          ! If no ground level population
          if (NH(1).le.0d0) then

            ! If master
            if (pid.eq.0) then
              umsg = ' - You provided H populations, but there '// &
                     'were no populations in the ground level'
              call verbose
            end if

            ! Exit
            return

          end if ! No ground level population

          ! If no ground level population
          if (NH(6).le.0d0) then

            ! If master
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

        ! Master checks if difference is significative
        if (pid.eq.0) then

          ! Compute relative difference
          dd = abs(Atmo%nha(iz) - sum(NH))/ &
                  (Atmo%nha(iz) + sum(NH))

          ! If relative difference > parameter in header
          if (dd.gt.delta.and.Atmo%nha(iz).gt.0d0) then
            if (ztau) then
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,f4.2,1x,A,1x,f4.2,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: After updating the H populations in '// &
                'the atmosphere at tau',iz,Atmo%z(iz), &
                ', the relative change',dd,'is above',delta,'.', &
                Atmo%nha(iz),'==>',sum(NH)
                call verbose
            else
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,f4.2,1x,A,1x,f4.2,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: After updating the H populations in '// &
                'the atmosphere at height',iz,Atmo%z(iz)*1e-5, &
                'km, the relative change',dd,'is above',delta,'.', &
                Atmo%nha(iz),'==>',sum(NH)
                call verbose
            end if
          end if
        end if ! Master

        Atmo%nha(iz) = sum(NH)

        ! Less total hydrogen than atomic hydrogen
        if (Atmo%nht(iz).lt.Atmo%nha(iz)) then

          ! Compute relative difference
          dd = abs(Atmo%nha(iz) - Atmo%nht(iz))/ &
                  (Atmo%nha(iz) + Atmo%nht(iz))

          ! Master checks if difference is significative
          if (pid.eq.0.and.Atmo%nht(iz).gt.0d0.and.dd.gt.delta) then
            if (ztau) then
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: Input atomic hydrogen larger than '// &
                'previous total hydrogen. Making them equal at '// &
                'tau',iz,Atmo%z(iz),':',Atmo%nht(iz),'==>', &
                Atmo%nha(iz)
              call verbose
            else
              write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                         '1x,es13.6,1x,A,es13.6)') &
                ' # Warning: Input atomic hydrogen larger than '// &
                'previous total hydrogen. Making them equal at '// &
                'height',iz,Atmo%z(iz)*1d-5,':',Atmo%nht(iz), &
                '==>',Atmo%nha(iz)
              call verbose
            end if
          end if ! Master

          Atmo%nht(iz) = Atmo%nha(iz)

        end if

      end do ! Heights

      ! Check if everything is fine
      call control

      return

      end subroutine ReviseHatmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Scale the total amount of H\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data
      subroutine ReviseHatom(Atom,Atmo)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iz
      double precision:: N0,N1,dd


      ! For each height
      do iz=1,nz

        N0 = sum(Atom%popu(:,iz))
        Atom%popu(:,iz) = Atom%popu(:,iz)*Atmo%nHa(iz)/N0
        N1 = sum(Atom%popu(:,iz))

        ! Master checks if difference is significative
        if (pid.eq.0) then

          ! Compute relative difference
          dd = abs(N0 - N1)/(N0 + N1)

          ! If relative difference > parameter in header
          if (dd.gt.delta.and.Atmo%nha(iz).gt.0d0) then
            write(umsg,'(A,1x,i3,1x,es13.6,1x,A,'// &
                       '1x,f4.2,1x,A,1x,f4.2,A,'// &
                       '1x,es13.6,1x,A,es13.6)') &
              ' # Warning: After updating the H populations in '// &
              'the atom at height',iz,Atmo%z(iz)*1e-5, &
              'km, the relative change',dd,'is above',delta,'.', &
              N0,'==>',N1
              call verbose
          end if
        end if ! Master

      end do ! Heights

      ! Check if everything is fine
      call control

      return

      end subroutine ReviseHatom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Converts into proper populations
      !!  Atom(Atom_class): Structure with the atomic data\n
      !!      dir(integer): Direction of correction
      subroutine correctpop(Atom,dir)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      integer, intent(in):: dir

      ! Local

      integer:: iz

      ! Normalize
      if (dir.eq.0) then
        do iz=Rz0,Rz1
          Atom%popu(:,iz) = Atom%popu(:,iz)/Atom%n(iz)
        end do
      ! De-normalize
      else if (dir.eq.1) then
        do iz=Rz0,Rz1
          Atom%popu(:,iz) = Atom%popu(:,iz)*Atom%n(iz)
        end do
      ! Error
      else
        umsg = 'Direction not recognized'
        urou = 'correctpop'
        call abortedS(umsg,urou,-1,.True.,.True.)
      end if

      ! Check if everything is fine
      call control

      return

      end subroutine correctpop

!#####################################################################
!#####################################################################
!#####################################################################

      end module initpopu_mod
