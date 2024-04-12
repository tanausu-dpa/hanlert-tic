      !> Reading atomic data
      module ratom_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     04/03/2024 V3.0.9
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     04/03/2024:    V3.0.9 - Set to zero the size of unused nblk when
!                             initializing the sizes (TdPA)
!
!     10/04/2023:    V3.0.8 - Bugfix: The number of M levels in LTE
!                             lines has to be defined in here, because
!                             the normalization (where they used to
!                             be) is not mandatory (TdPA)
!
!     09/29/2023:    V3.0.7 - Active atoms can now decide the rank
!                             limits in term-wise way. If they do
!                             this, the rank limit of the transitions
!                             it automatically set from those (TdPA)
!
!     09/25/2023:    V3.0.6 - Added set_atom_label routine (TdPA)
!
!     08/07/2023:    V3.0.5 - Receive frequencies in cm^-1 and
!                             Einstein coeff. in s^-1, and transform
!                             after reading. Added warning for very
!                             small energies that may mean that an
!                             old model format is being used (TdPA)
!                           - Added setup_LTE_transition and
!                             remove_LTE_transition (TdPA)
!
!     07/03/2023:    V3.0.4 - Added the possibility for an atom to
!                             not contribute to the frequency axis
!                             nodes (TdPA)
!
!     10/26/2022:    V3.0.3 - Changed the indexing storage (TdPA)
!                           - Chose a more reasonable size for the
!                             Atom%ifst_ij array (TdPA)
!
!     10/25/2022:    V3.0.2 - Initialize the p_T pointer and clean
!                             it at the end (TdPA)
!                           - Moved the allocation of Atom%eval
!                             and Atom%evec elsewhere (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: The MPI dependent quantities
!                             cannot be initialized in ratom (TdPA)
!                           - Added routine prepareatomMPI (TdPA)
!                           - Added routine initMblock (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Variables with height dimension are
!                                now allocated elsewhere.
!                              o The term damping is initialized
!                                elsewhere.
!                              o Removed the option for explicit
!                                elastic collisions.
!                              o Changed aborted calls to gaborted.
!                              o Added routine prepareatom.
!                             (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     07/31/2020:    V1.5.4 - If there are no photoionizations, a
!                             dummy Atom%phot is allocated (TdPA)
!
!     03/05/2020:    V1.5.3 - Now shiftatoms shift also the population
!                             files (TdPA)
!
!     11/19/2019:    V1.5.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     11/12/2019:    V1.5.1 - Allocates Atom%rif0 and Atom%rif1
!                             variables (TdPA)
!                           - Nullifies Atom%Normp, that now is a
!                             pointer (TdPA)
!                           - Bugfix: Flags the customized H atom as
!                             multilevel (TdPA)
!
!     09/26/2019:    V1.5.0 - All quantities that needed atmospheric
!                             data to be computed have been moved
!                             elsewhere. That means that the
!                             collisional data is stored in a temporal
!                             container to be used later and a
!                             considerable fraction of the source has
!                             been moved (TdPA)
!                           - Now we store the degeneration of the
!                             terms in the atomic structure (TdPA)
!
!     09/24/2019:    V1.4.2 - Moved computation of total atom
!                             population to the initpopu module (TdPA)
!
!     08/16/2019:    V1.4.1 - Bugfix: The level variable was being
!                             used for stage, instead of the term, for
!                             inelastic collisions between atomic
!                             levels (TdPA)
!
!     08/14/2019:    V1.4.0 - Important revamp of how the inelastic
!                             collisions are processed. Added more
!                             options besides electrons (TdPA)
!
!     08/08/2019:    V1.3.8 - Wrong correction when photoionizations
!                             are given in inverse order, was missing
!                             an index (TdPA)
!
!     07/23/2019:    V1.3.7 - Now there are more possible flags for
!                             how the collisions must be
!                             interpolated (TdPA)
!                           - The collisions are not scaled within
!                             intercol anymore, but here (TdPA)
!                           - Bugfix: The indexing of transitions is
!                             the same than in the atomic file, so
!                             every part of the code agrees in the
!                             order (TdPA)
!                           - In multilevel, the term-term collisions
!                             are just copied (TdPA)
!                           - Bugfix: For ML atoms, quantum numbers
!                             are ignored to build the Ecoeff for FS
!                             transitions (TdPA)
!
!     06/04/2019:    V1.3.6 - Added indexing of FS transitions to
!                             avoid searchs in rtcoeffi (TdPA)
!
!     06/03/2019:    V1.3.5 - Bugfix: shiftatoms tried to shift non
!                             exising atoms when no background atom
!                             was specified (TdPA)
!                           - Now splitf is read in the radiative
!                             transition section (TdPA)
!
!     05/08/2019:    V1.3.4 - Implemented needed changes to allow the
!                             new indexing of radiation field
!                             tensors (TdPA)
!
!     03/18/2019:    V1.3.3 - Now saves the terms and levels of
!                             each b-b anf b-f transition (TdPA)
!                           - In multilevel atoms, do not check
!                             if the line is permitted (TdPA)
!
!     03/13/2019:    V1.3.2 - Now gives a more detailed warning for
!                             negative spline interpolation (TdPA)
!
!     03/12/2019:    V1.3.1 - Do not check J relations for transitions
!                             when reading a ML model (TdPA)
!
!     02/20/2019:    V1.3.0 - New verbosity (TdPA)
!                           - Checks for success of python routine
!                             and unit is now 100 (TdPA)
!                           - Added allocateatom and shiftatoms (TdPA)
!
!     09/26/2018:    V1.2.4 - Added automatic correction of order
!                             when photoionizations have the
!                             unexpected order (TdPA)
!
!     11/02/2017:    V1.2.3 - The multilevel version reads Landé
!                             factors (TdPA)
!
!     11/01/2017:    V1.2.2 - Added output of the treatment of
!                             collisions (TdPA)
!
!     10/31/2017:    V1.2.1 - Changed the rules to determine a
!                             forbidden collisional rate, to avoid
!                             transitions within the same term or that
!                             violate parity to be considered allowed
!                             (TdPA)
!
!     10/30/2017:    V1.2.0 - Now does not read landé factors, but
!                             the multilevel case is treated in a
!                             different way (TdPA)
!                           - Stores .25d0 in a parameter (TdPA)
!
!     10/11/2017:    V1.1.1 - Now reads landé factors (TdPA)
!
!     09/22/2017:    V1.1.0 - Added option to limit K values (TdPA)
!
!     09/14/2017:    V1.0.7 - Added a path and ID to the file (TdPA)
!
!     09/08/2017:    V1.0.6 - The forbidden flag for collisional
!                             ionizations is 2 now (TdPA)
!                           - Removed Dfreq2 from colinter call (TdPA)
!
!     07/19/2017:    V1.0.5 - Added message for absence of inelastic
!                             collisions (TdPA)
!                           - Bugfix: Non-active atoms were trying
!                             to deallocate some non-allocated
!                             variables (TdPA)
!                           - Bugfix: Initializing abundances of
!                             passive elements (TdPA)
!                           - Bugfix: Passive atoms should not
!                             deallocate Atom%iphot and Atom%fst,
!                             they are needed for background
!                             calculations (TdPA)
!                           - Changed the symbol of the H default
!                             message from # to - (TdPA)
!
!     06/15/2017:    V1.0.4 - Bugfix: fcflag initialized (TdPA)
!                           - Bugfix: The rule of the total angular
!                             momentum for collisions to determine
!                             fcflag was a singleline if that was
!                             suppossed to contain two, changed to
!                             a then/endif command (TdPA)
!                           - Bugfix: Ionizing collisions must be
!                             always flagged as forbidden,
!                             independently of quantum numbers (TdPA)
!
!     06/12/2017:    V1.0.3 - Allocate if0, if1, W0 and W1 in atom
!                             structure (TdPA)
!
!     05/05/2017:    V1.0.2 - Flagging forbidden collisions
!                             automatically when they do not comply
!                             with the dipole selection rules. Not
!                             sure the previous way was compatible
!                             with the new SEbuild anymore (TdPA)
!
!     04/27/2017:    V1.0.1 - Bugfix: The indexing of the FS
!                             transitions cannot be generated before
!                             finishing reading the transition
!                             information, in rAtom (TdPA)
!                           - Bugfix: Trying to initialize damp to 0,
!                             before allocating it, in AtomH (TdPA)
!
!     04/18/2017:    V1.0.0 - First version (TdPA)
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
!  This subroutine reads the atomic data
!
!  ratom:
!    This is the general atomic file reader
!
!  AtomH:
!    Generates an ad-hoc hydrogen atom
!
!  abund:
!    Renormalize abundance modificators if necessary
!
!  setFScoeff:
!    Calculates Aul and Blu for individual fine-structure transitions
!
!  setup_LTE_transition:
!    Prepare static quantities for LTE lines
!
!  remove_LTE_transition:
!    Check if an LTE transition is within an atomic model, and remove
!  it from the model (TdPA)
!
!  allocateatom:
!    Allocates array of Atom_class
!
!  set_atom_label
!    Ensures that each atom has an unique label for output files
!
!  shiftatoms:
!    Shift to the right the list of atoms
!
!  prepareatom:
!    Allocate some of the variables with height dimension
!
!  prepareatomMPI:
!    Allocate the variables with dependence on the number of RT
!  processes
!
!  initMblock:
!    Initialize the sizes of the M blocks assuming there is a magnetic
!  field
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use funnj_mod
      use parameters_mod , only : c, c2, dopp, PI, ConvF , TINYSP
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Reads a file with atomic data.\n
      !!   filename(character(:)): Name of the file to read\n
      !!     source(character(:)): Path to the source code\n
      !!         ID(character(:)): ID of this run\n
      !!       skip_wave(logical): This atom does not contribute
      !!                           with wavelengths\n
      !! Kcut_input(integer(:,:)): Term wise K cut input data\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!            indx(integer): Index in the list of atoms\n
      !!           isPRD(logical): Bool to store if this atom has
      !!                           lines in PRD\n
      !!          active(logical): Bool to specify if this atom is
      !!                           active or not
      subroutine rAtom(filename,source,ID,skip_wave,Kcut_input, &
                       Atom,indx,isPRD,active)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      logical, intent(in):: active,skip_wave
      logical, intent(inout):: isPRD
      integer, intent(in):: indx
      integer, dimension(:,:), allocatable, intent(in):: Kcut_input

      ! Local

      character(len=1):: cdump

      logical:: flin,changed

      integer:: ios,iterm,iterm1,ilevel,ilevel1
      integer:: ii,jj,i,itran,ftran,icol,stagl,stagu
      integer:: K,iQ,iJ,iJ1,iJl,iJu,up,low,Tindex,minK,maxK,maxfst
      integer:: l1,l2,l3,l4,nTmp,nclin,cind,col_type,nion

      double precision, parameter:: dLSJ = .25d0

      double precision:: d1,d2,d3,rL,S,rJ,rJ1,rJmin,rJmax
      double precision:: freq,ftmp,radAcoeff
      double precision, dimension(:), allocatable:: vaux1,vaux2

      ! Pointer
      type(Tbox_class), pointer:: p_T

      ! Initialize
      nullify(p_T)

      ! Read the atomic data
      if(pid.eq.0) then

        umsg = ' - Read atom '//trim(filename)
        call verbose

        call system('python '//trim(source)//'ratom.py '// &
                    trim(filename)//' '//ID//' '//verbosef)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      open(100,file='tmp_atom_'//ID,status='old',iostat=ios,err=1000)

      ! Success
      read (100,*,err=1100) ios

      ! If no correct file, abort
      if (ios.lt.0) then

        umsg = 'Problem translating the atomic file '//trim(filename)
        goto 1200

      end if

      ! Read atomic ID
      read (100,*,err=1100) Atom%Element
      if (Atom%Element(2:2).eq.' ') then
        Atom%Element(2:2) = Atom%Element(1:1)
        Atom%Element(1:1) = ' '
      end if
      ! Atomic mass
      read (100,*,err=1100) Atom%rmass
      ! Atomic abundance (12 + log(Atom/H))
      read (100,*,err=1100) Atom%abun
      Atom%abun = 1d1**(Atom%abun - 12d0)
      ! Multiplicative modifier to this abundance
      read (100,*,err=1100) Atom%abun_mod
      ! Normalize abundance multiplier to 1 for same ID atoms
      read (100,*,err=1100) ios
      ! Only taken into account if active
      if (active) then
        if (ios.eq.1) then
          Atom%anorm = .True.
        else
          Atom%anorm = .False.
        end if
      ! If passive initialize the population already
      end if

      ! Read type of model
      read (100,*,err=1100) l1
      Atom%ML = l1.eq.1

      ! Read dimensions of the model
      read (100,*,err=1100) Atom%nMulti,Atom%ntran,Atom%nphot,Atom%ngk
      read (100,*,err=1100) Atom%nJmax
      Atom%NNN = Atom%nMulti*Atom%nJmax

      ! Point pointer of extra collisions to null
      nullify(Atom%Ccoeff_special)


      !
      ! Level information
      !

      ! Allocations
      ! Orbital angular momentum
      allocate(Atom%rLval(Atom%nMulti))
      ! Spin angular momentum
      allocate(Atom%Sval(Atom%nMulti))
      ! Number of FS levels
      allocate(Atom%nJ(Atom%nMulti))
      ! Ionization stage
      allocate(Atom%stage(Atom%nMulti))
      ! Angular momentum
      allocate(Atom%rJval(Atom%nJmax,Atom%nMulti))
      Atom%rJval = -1d0
      ! Frequency of FS levels
      allocate(Atom%FSfreq(Atom%nJmax,Atom%nMulti))
      ! Frequency of term
      allocate(Atom%TRfreq(Atom%nMulti))
      ! Term degeneracy
      allocate(Atom%deg(Atom%nMulti))
      ! Lande factor
      if (Atom%ML) allocate(Atom%gL(Atom%nMulti))

      ! Level index
      l3 = 0

      ! For each term
      do iterm=1,Atom%nMulti

        ! Read orbital and spin ang. moments and ionization stage
        read (100,*,err=1100) rL, S, l1

        if (nint(2d0*S).gt.nxS) nxS = nint(2d0*S)
        if (nint(2d0*rL).gt.nxL) nxL = nint(2d0*rL)

        Atom%rLval(iterm) = rL
        Atom%Sval(iterm) = S
        Atom%stage(iterm) = l1

        ! If multilevel
        if (Atom%ML) then

          Atom%nJ(iterm) = 1

          ! Read angular momentum, energy, term index and level index
          read (100,*,err=1100) rJ, freq, d1, l1, l2

          l3 = l3 + 1

          ! Sanity checks
          if (l1.ne.iterm.or.l2.ne.l3) then

            if (l1.ne.iterm) then
              umsg = 'The term index does not '// &
                     'correspond to the term'
            else
              umsg = 'The level index does not '// &
                     'correspond to the level'
            end if

            goto 1200

          end if

          ! Energy units
          freq = freq*1d-5

          Atom%rJval(1,iterm) = rJ
          Atom%deg(iterm) = 2d0*rJ + 1d0
          Atom%FSfreq(1,iterm) = freq
          Atom%TRfreq(iterm) = freq

          if (nint(2d0*rJ).gt.nxJ) nxJ = nint(2d0*rJ)

          ! If input Landé factor non-physical
          if (d1.lt.0d0) then

            ! If J is 0
            if (rJ.lt.dLSJ) then

              Atom%gL(iterm) = 1d0

            ! If J is not 0
            else

              Atom%gL(iterm) = 1d0 + .5d0*(rJ*(rJ+1d0) + S*(S+1d0) - &
                                           rL*(rL+1d0))/rJ/(rJ+1d0)

            end if ! J value

          ! Physical Landé factor
          else

            Atom%gl(iterm) = d1

          end if

        ! If multiterm
        else

          ! Calculate degeneration of term
          Atom%deg(iterm) = (2d0*rL + 1d0)*(2d0*S + 1d0)

          ! Minimum and maximum angular momentums
          rJmin = abs(rL - S)
          rJmax = rL + S

          ! Number of FS levels for this term
          Atom%nJ(iterm) = nint(rJmax - rJmin) + 1

          ! Cumulative sum
          ftmp = 0d0

          ! For each level within the term
          do iJ=1,Atom%nJ(iterm)

            ! Read angular momentum, energy, term index and level
            ! index
            read (100,*,err=1100) rJ, freq, l1, l2

            l3 = l3 + 1

            ! Sanity checks
            if (l1.ne.iterm.or.l2.ne.l3) then

              if (l1.ne.iterm) then
                umsg = 'The term index does not '// &
                       'correspond to the term'
              else
                umsg = 'The level index does not '// &
                       'correspond to the level'
              end if

              goto 1200

            end if

            ! Energy units
            freq = freq*1d-5

            Atom%rJval(iJ,iterm) = rJ
            Atom%FSfreq(iJ,iterm) = freq

            if (nint(2d0*rJ).gt.nxJ) nxJ = nint(2d0*rJ)

            ! Add weighted contribution to the term energy
            ftmp = (2D0*rJ + 1D0)*freq + ftmp

          end do

          ! Calculate the multiplet energy as a weighted average of
          ! the F-levels energies
          Atom%TRfreq(iterm) = ftmp/Atom%deg(iterm)

        end if ! Multilevel or multiterm

      end do

      ! Sanity check
      if (maxval(Atom%TRfreq).lt.1d-4) then

        write(umsg,'(A)') ' # Maximum energy in '// &
             Atom%element//' model atom is smaller than 10 cm^-1, '// &
             'you may be using an old atomic model format'
        call verbose

      end if

      ! Determine the number of FS levels in this atom
      Atom%nlevel = sum(Atom%nJ)

      ! Calculate the maximum K multipole
      Atom%nKmax = nint(maxval(Atom%rJval)*2d0)
      if (Atom%nKmax.gt.nkx.and.active) nkx = Atom%nKmax

      ! Calculate the maximum number of magnetic levels
      Atom%nMmax = nint(maxval(Atom%rJval)*2d0 + 1d0)

      ! Allocate atomic quantities to be used later
      allocate(Atom%nblk(Atom%nMmax,Atom%nMulti))
      allocate(Atom%iJval(Atom%nJmax,Atom%nMmax,Atom%nMulti))

      ! Define "super-indexes" to handle large number of dimensions
      ! Term given the level
      allocate(Atom%term(Atom%nlevel))
      ! Sublevel in its term given the global level
      allocate(Atom%sublevel(Atom%nlevel))
      ! Allocate indexing by term
      allocate(Atom%irho(Atom%nMulti))

      ! Initialize
      ii = 0

      ! For each term
      do iterm=1,Atom%nMulti

        ! Allocate rho indexing
        allocate(Atom%irho(iterm)%irho_ij(Atom%nJ(iterm)))

        ! For each level
        do iJ=1,Atom%nJ(iterm)

          ii = ii + 1
          Atom%irho(iterm)%irho_ij(iJ) = ii
          Atom%term(ii) = iterm
          Atom%sublevel(ii) = iJ

        end do
      end do

      ! Update the maximum number of transitions (b-b and b-f)
      if (active) then
        nxtran = nxtran + Atom%ntran
        nxphot = nxphot + Atom%nphot
      end if

      ! Initialize shifts
      Atom%tshift = 0
      Atom%tfshift = 0
      Atom%pshift = 0


      !
      ! Radiative transitions
      !

      ! Allocations
      ! Structure for FS information
      allocate(Atom%fst(Atom%ntran))
      ! Energy of the term-term transition
      allocate(Atom%Dfreq(Atom%ntran))
      ! Type of elastic broadening
      allocate(Atom%broad_type(Atom%ntran))
      ! Arguments to calculate the broadening
      allocate(Atom%broad_args(4,Atom%ntran))
      ! Argument for Stark broadening
      allocate(Atom%broad_stark(Atom%ntran))
      ! Number of frequencies to use in the line
      allocate(Atom%nfreqt(Atom%ntran))
      ! How many of those frequencies fo to the core region
      allocate(Atom%nfreqtc(Atom%ntran))
      ! Doppler widths that the line spawns
      allocate(Atom%dwvl(Atom%ntran))
      ! Dopple widths that the core of the line spawns
      allocate(Atom%dwvlc(Atom%ntran))
      ! Flag of second order emissivity
      allocate(Atom%lemiss2(Atom%ntran))
      ! Flag of Split components
      allocate(Atom%splitf(Atom%ntran))
      ! Flag of line presence
      allocate(Atom%fflag(Atom%ntran))
      ! Lower frequency limits
      allocate(Atom%if0(Atom%ntran))
      ! Lower absolute frequency limits
      allocate(Atom%rif0(Atom%ntran))
      ! Upper absolute frequency limits
      allocate(Atom%rif1(Atom%ntran))
      ! Upper frequency limits
      allocate(Atom%if1(Atom%ntran))
      ! Lower frequency weight
      allocate(Atom%W0(Atom%ntran))
      ! Upper frequency weight
      allocate(Atom%W1(Atom%ntran))
      ! Indexing of radiative transitions
      allocate(Atom%irad(Atom%nMulti,Atom%nMulti))
      Atom%irad = 0
      ! Einstein coefficients Aul and Blu
      allocate(Atom%Ecoeff(Atom%nMulti,Atom%nMulti))
      Atom%Ecoeff = 0d0
      ! Nullify Norm
      nullify(Atom%Normp)

      ! Reset counters
      Atom%nfreq = 0
      Atom%nftran = 0

      ! For each transition
      do itran=1,Atom%ntran

        ! Read the terms involved and Aul
        read (100,*,err=1100) iterm1,iterm,radAcoeff

        ! Rate unit
        radAcoeff = radAcoeff*1d-8

        ! iterm1 > iterm
        if (iterm1.lt.iterm) then
          l1 = iterm1
          iterm1 = iterm
          iterm = l1
        end if

        ! Enumerate permitted transitions
        Atom%irad(iterm1,iterm) = itran
        Atom%irad(iterm,iterm1) = itran
        Atom%fst(itran)%iterml = iterm
        Atom%fst(itran)%itermu = iterm1

        ! Calculate the frequency of the transition
        Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                                Atom%TRfreq(iterm))

        ! Read collisional broadening parameters
        read (100,*,err=1100) l1
        Atom%broad_type(itran) = l1
        do l1=1,4
          read (100,*,err=1100) d1
          Atom%broad_args(l1,itran) = d1
        end do


        ! Calculate Aul and Blu
        Atom%Ecoeff(iterm1,iterm) = radAcoeff
        Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                    Atom%deg(iterm1)/ &
                                    Atom%deg(iterm)/ &
                                   (ConvF*Atom%Dfreq(itran)* &
                                    1d21*(2d0*c)* &
                                    Atom%Dfreq(itran)**2d0)


        ! Read the Stark broadening parameter, the number of
        ! frequencies for this transition and its core, the
        ! Doppler widths of this transition and its core and
        ! the flag of second order emissivity, split between
        ! components
        read (100,*,err=1100) d3,l1,l2,d1,d2,l3,l4

        ! Store Stark broadening information
        Atom%broad_stark(itran) = d3

        ! Store node information
        if (mod(l1,2).eq.0) l1 = l1 + 1
        if (mod(l2,2).eq.0) l2 = l2 + 1
        Atom%nfreqt(itran) = l1
        Atom%nfreqtc(itran) = l2
        Atom%Dwvl(itran) = d1
        Atom%Dwvlc(itran) = d2

        ! Store the 2nd order emissivity flag
        if (l3.ne.0) then
          Atom%lemiss2(itran) = .True.
          isPRD = .True.
        else
          Atom%lemiss2(itran) = .False.
        end if

        ! Store the split components flag
        if (l4.ne.0) then
          Atom%splitf(itran) = .True.
        else
          Atom%splitf(itran) = .False.
        end if

        ! Allocate FS transition indexing
        allocate(Atom%fst(itran)%irad(Atom%nJ(iterm1), &
                                      Atom%nJ(iterm)))
        Atom%fst(itran)%irad = 0
        ! Allocate Aul for FS transitions
        allocate(Atom%fst(itran)%Aul(Atom%nJ(iterm1),Atom%nJ(iterm)))
        Atom%fst(itran)%Aul = 0d0
        ! Allocate Blu for FS transitions
        allocate(Atom%fst(itran)%Blu(Atom%nJ(iterm),Atom%nJ(iterm1)))
        Atom%fst(itran)%Blu = 0d0

        ! If Multi-level
        if (Atom%ML) then

          Atom%fst(itran)%irad(1,1) = 1
          Atom%fst(itran)%nt = 1
          Atom%nftran = Atom%nftran + 1
          maxfst = 1
          l3 = 1
          allocate(Atom%fst(itran)%ilevell(l3))
          allocate(Atom%fst(itran)%ilevelu(l3))
          Atom%fst(itran)%ilevell = 1
          Atom%fst(itran)%ilevelu = 1

        else

          ! Initialize maximum FS trans
          maxfst = 0

          ! Count FS transitions
          l3 = 0
          do iJu=1,Atom%nJ(iterm1)
            do iJl=1,Atom%nJ(iterm)
              if(abs(Atom%rJval(iJu,iterm1) - &
                     Atom%rJval(iJl,iterm)).gt.1 .or. &
                 Atom%rJval(iJu,iterm1) + &
                 Atom%rJval(iJl,iterm) .lt. .4d0)cycle
              l3 = l3 + 1
              Atom%fst(itran)%irad(iJu,iJl) = l3
            end do
          end do
          Atom%fst(itran)%nt = l3
          if (l3.gt.maxfst) maxfst = l3
          Atom%nftran = Atom%nftran + l3

          ! Allocate
          allocate(Atom%fst(itran)%ilevell(l3))
          Atom%fst(itran)%ilevell = -1
          allocate(Atom%fst(itran)%ilevelu(l3))
          Atom%fst(itran)%ilevelu = -1

          ! Index
          l3 = 0
          do iJu=1,Atom%nJ(iterm1)
            do iJl=1,Atom%nJ(iterm)
              if(abs(Atom%rJval(iJu,iterm1) - &
                     Atom%rJval(iJl,iterm)).gt.1 .or. &
                 Atom%rJval(iJu,iterm1) + &
                 Atom%rJval(iJl,iterm) .lt. .4d0)cycle
              l3 = l3 + 1
              Atom%fst(itran)%ilevell(l3) = iJl
              Atom%fst(itran)%ilevelu(l3) = iJu
            end do
          end do

        end if

        ! Cycle if passive atom
        if (.not.active) cycle
        ! Add the frequency nodes to the total count of this atom
        Atom%nfreq = Atom%nfreq + l1*l3

      end do

      if (active) then

        ! Update the maximum number of FS transitions
        nxt = nxt + Atom%nftran

        ! Index the FS transitions given the term transition and the
        ! sub FS transition
        allocate(Atom%ifst_ij(maxfst,Atom%ntran))
        Atom%ifst_ij = 0
        allocate(Atom%ifst(Atom%nftran))
        Atom%ifst = 0

        ! Reset counter
        ftran = 0

        ! For each transition
        do itran=1,Atom%ntran

          ! For each FS transition
          do l3=1,Atom%fst(itran)%nt

            ftran = ftran + 1
            Atom%ifst_ij(l3,itran) = ftran
            Atom%ifst(ftran) = itran

          end do ! FS transitions
        end do ! Term transition

        ! Add the frequency nodes of this atom to the total count
        if (.not.skip_wave) nfreq = nfreq + Atom%nfreq

        !
        ! Control K cuts
        !

        ! Term-wise K cut and transition-wise Krad
        allocate(Atom%Kcut(Atom%nMulti))
        Atom%Kcut = Kcut
        allocate(Atom%Krad(Atom%ntran))
        Atom%Krad = Krad

        ! If there are inputs, we need to take them into account
        if (allocated(Kcut_input)) then

          ! Initialize flag
          changed = .False.

          ! For each entry
          do ii=1,size(Kcut_input,2)

            ! If atom is index
            if (indx.eq.Kcut_input(1,ii)) then

              ! Sanity check
              if (Kcut_input(2,ii).lt.0.or. &
                  Kcut_input(3,ii).gt.Atom%nMulti) then

                ! Abort
                umsg = 'Wrong index in Kcut for atom '// &
                       trim(Atom%Element)
                urou = 'rAtom'
                call gaborted

              end if

              ! Flag changed
              changed = .True.

              ! For specified range
              do jj=Kcut_input(2,ii),Kcut_input(3,ii)

                ! For this term, apply K cut
                Atom%Kcut(jj) = min(Kcut_input(4,ii), &
                                   nint(2d0*maxval(Atom%rJval(:,jj))))

              end do

            end if

          end do ! Input entries

          ! If there were changes in the Kcut, control Krad as well
          if (changed) then

            ! For each transition
            do itran=1,Atom%ntran

              ! Get terms
              iterm = Atom%fst(itran)%iterml
              iterm1 = Atom%fst(itran)%itermu

              ! Is there PRD?
              if (isPRD) then

                ! Check if PRD transition
                changed = Atom%lemiss2(itran)

                ! If no PRD, check all
                if (.not.changed) then

                  ! Check all terms
                  do icol=1,Atom%nMulti

                    ! Check transition with same upper term
                    ftran = Atom%irad(icol,iterm1)

                    ! Cycle null
                    if (ftran.le.0) cycle

                    ! If PRD
                    if (Atom%lemiss2(ftran)) then
                      changed = .True.
                      exit
                    end if

                  end do ! All terms

                end if ! Current transition no PRD

                ! If PRD in upper term
                if (changed) then

                  ! Check cut
                  if (Atom%Kcut(iterm1).ne.Kcut) then

                    ! Abort
                    umsg = 'Please, do not add specific K '// &
                           'limits to terms involved in '// &
                           'PRD transitions'
                    urou = 'rAtom'
                    call gaborted

                  end if ! Edited K cut

                  ! Check all terms
                  do icol=1,Atom%nMulti

                    ! Check transition with same upper term
                    ftran = Atom%irad(icol,iterm1)

                    ! Cycle null
                    if (ftran.le.0) cycle

                    ! Check cut
                    if (Atom%Kcut(icol).ne.Kcut) then

                      ! Abort
                      umsg = 'Please, do not add specific K '// &
                             'limits to terms involved in '// &
                             'PRD transitions'
                      urou = 'rAtom'
                      call gaborted

                    end if ! Edited K cut

                  end do ! All terms

                end if ! Upper level in PRD transition
              end if ! There is PRD at all

              ! Max K
              Atom%Krad(itran) = min(Atom%Kcut(iterm) + &
                                     Atom%Kcut(iterm1),2)

            end do ! Transitions

          end if ! Changes in term-wise K cut
        end if ! There are custom K cut


        !
        ! Define the index matrix
        !

        ! Initialize running index
        i = 0

        ! For each term
        do iterm=1,Atom%nMulti

          ! Allocate J index
          allocate(Atom%irho(iterm)%Jrho(Atom%nJ(iterm), &
                                         Atom%nJ(iterm)))

          ! For each level
          do iJ=1,Atom%nJ(iterm)

            ! Get J
            rJ = Atom%rJval(iJ,iterm)

            ! For each other level
            do iJ1=1,Atom%nJ(iterm)
           !do iJ1=iJ,iJ ! no J,J'

              ! Get J
              rJ1 = Atom%rJval(iJ1,iterm)

              ! Get K limits
              minK = nint(abs(rJ-rJ1))
              maxK = min(nint(rJ+rJ1),Atom%Kcut(iterm))

              ! Allocate KK tab
              allocate(Atom%irho(iterm)%Jrho(iJ1,iJ)% &
                            kq(-maxK:maxK,minK:maxK))
              Atom%irho(iterm)%Jrho(iJ1,iJ)%kq = 0

              ! For each K
              do K=minK,maxK
                ! For each Q
                do iQ=-K,K

                  i = i + 1
                  Atom%irho(iterm)%Jrho(iJ1,iJ)%kq(iQ,K) = i

                end do
              end do

            end do

          end do
        end do

        ! Update maximum dimensionality of SEE
        Atom%ndim = i

        if (pid.eq.0) then
          write(umsg,'(A,i6)') '   S.E. dimension =',Atom%ndim
          call verbose
        end if

      end if

      ! Sanity check if there are transitions
      if (Atom%ntran.gt.0) then

        ftmp = 0d0
        do iterm=1,Atom%nMulti-1

          ftmp = max(ftmp, &
                     maxval(Atom%Ecoeff(iterm+1:Atom%nMulti,iterm)))
        end do

        if (ftmp.lt.1d-9) then
          write(umsg,'(A)') ' # Maximum Einstein '// &
              'coefficient in '//Atom%element// &
              ' model atom is smaller than 0.1 s^-1, '// &
              'you may be using an old atomic model format'
          call verbose
        end if
      end if


      !
      ! Elastic collisions
      !

      ! If there are inputs
      if (Atom%ngk.ge.1) then

        ! Allocate elastic data
        allocate(Atom%elas(Atom%ngk))

        ! For each entry
        do ii=1,Atom%ngk

          ! Read to which level corresponds and find term and sublevel
          read (100,*,err=1100) ilevel,l1

          ! Store in elastic structure
          Atom%elas(ii)%ilevel = ilevel
          Atom%elas(ii)%nentry = l1

          ! Allocate entries
          allocate(Atom%elas(ii)%datum(l1))

          ! For each line in this entry
          do jj=1,l1

            ! Read the multipole, the type of input and the
            ! dimensionality
            read (100,*,err=1100) K
            read (100,*,err=1100) cdump
            read (100,*,err=1100) l2

            ! Parse into database
            Atom%elas(ii)%datum(jj)%K = K
            Atom%elas(ii)%datum(jj)%nz = l2

            ! If it is a fit input
            if(cdump.eq.'f')then

              ! Read the coefficients
              read (100,*,err=1100) d1,d2,d3

              ! Put in database
              Atom%elas(ii)%datum(jj)%typo = 0
              Atom%elas(ii)%datum(jj)%a = d1
              Atom%elas(ii)%datum(jj)%b = d2
              Atom%elas(ii)%datum(jj)%c = d3

            ! No type of input recognized
            else

              umsg = 'Mode of elastic collisions '// &
                     'not recognized'
              goto 1200

            end if

          end do ! Input sub entry
        end do ! Input entry

        ! If passive, forget everything
        if (.not.active) then

          ! For each entry
          do ii=1,Atom%ngk

            ! For each subentry
            do jj=1,Atom%elas(ii)%nentry

              ! If explicit, deallocate array
              if (Atom%elas(ii)%datum(jj)%typo.eq.1) then
                deallocate(Atom%elas(ii)%datum(jj)%Coeff)
              end if

            end do ! Subentry

            ! Deallocate subentry
            deallocate(Atom%elas(ii)%datum)

          end do ! Entries

          ! Deallocate structure
          deallocate(Atom%elas)

        end if ! Passive

      ! There is no input
      else

        if(pid.eq.0.and.active) then
          write(umsg,'(A)') '   No elastic collisions'
          call verbose
        end if

      end if


      !
      ! Photoionizations
      !

      ! Allocations
      ! Indexing of b-f transitions
      allocate(Atom%iphot(Atom%nlevel,Atom%nlevel))
      Atom%iphot = 0

      ! If there are inputs
      if (Atom%nphot.ge.1) then

        ! Allocate photoionization information
        allocate(Atom%phot(Atom%nphot))

        ! For each input entry
        do ii=1,Atom%nphot

          ! Read the levels involved
          read (100,*,err=1100) ilevel1,ilevel

          ! ilevel1 > ilevel
          if (ilevel1.lt.ilevel) then
            l1 = ilevel1
            ilevel1 = ilevel
            ilevel = l1
          end if

          ! Find the term and level indexes
          iterm = Atom%term(ilevel)
          iJ = Atom%sublevel(ilevel)
          iterm1 = Atom%term(ilevel1)
          iJ1 = Atom%sublevel(ilevel1)

          ! Check that they have actually different ionization stages
          if (Atom%stage(iterm1).eq.Atom%stage(iterm)) then
            umsg = 'Photoionization between terms in '// &
                   'the same stage.'
            goto 1200
          end if

          ! Index the b-f transition
          Atom%iphot(ilevel1,ilevel) = ii
          Atom%iphot(ilevel,ilevel1) = ii
          Atom%phot(ii)%ilevell = ilevel
          Atom%phot(ii)%ilevelu = ilevel1

          ! Calculate the edge of the photoionization
          Atom%phot(ii)%edge = Atom%FSfreq(iJ1,iterm1) - &
                               Atom%FSfreq(iJ,iterm)

          ! Store the ratio between the degeneration of the levels
          Atom%phot(ii)%glu = (2d0*Atom%rJval(iJ,iterm)+1d0)/ &
                              (2d0*Atom%rJval(iJ1,iterm1)+1d0)

          ! Read the type of input and the number of frequencies to
          ! use for this b-f transition
          read (100,*,err=1100) cdump
          read (100,*,err=1100) l1
          Atom%phot(ii)%nfreq = l1


          ! If the input is explicit
          if (cdump.eq.'e') then

            ! Allocate explicit data
            ! Frequencies in the inpu
            allocate(Atom%phot(ii)%infreq(l1))
            ! Cross sections in the input
            allocate(Atom%phot(ii)%inalpha(l1))
            Atom%phot(ii)%mode = 0

            ! Read the input frequencies and cross sections
            do jj=1,l1

              read(100,*,err=1100) d1, d2

              Atom%phot(ii)%infreq(jj) = 1d2/d1
              Atom%phot(ii)%inalpha(jj) = d2

            end do

            ! Check order
            if (l1.gt.1) then

              ! If wrong order, revert it
              if (Atom%phot(ii)%infreq(1).gt. &
                  Atom%phot(ii)%infreq(2)) then

                ! Check not already allocated
                if (.not.allocated(vaux1)) then
                  allocate(vaux1(l1))
                  allocate(vaux2(l1))
                else
                  if (l1.gt.size(vaux1)) then
                    deallocate(vaux1)
                    deallocate(vaux2)
                    allocate(vaux1(l1))
                    allocate(vaux2(l1))
                  end if
                end if

                vaux1 = Atom%phot(ii)%infreq
                vaux2 = Atom%phot(ii)%inalpha

                do jj=1,l1
                  Atom%phot(ii)%infreq(jj) = vaux1(l1 - jj + 1)
                  Atom%phot(ii)%inalpha(jj) = vaux2(l1 - jj + 1)
                end do

              end if
            end if

          ! If the input is hydrogenic
          else if (cdump.eq.'h') then

            ! Allocate hydrogenic data
            ! Maximum frequency to reach
            allocate(Atom%phot(ii)%infreq(1))
            ! Cross section at edge
            allocate(Atom%phot(ii)%inalpha(1))
            Atom%phot(ii)%mode = 1

            ! Read the maximum frequency and cross section at edge
            read(100,*,err=1100) d1, d2

            Atom%phot(ii)%infreq(1) = 1d2/d1
            Atom%phot(ii)%inalpha(1) = d2

          end if

          ! Add the nodes to the total count
          if (active.and..not.skip_wave) &
            nfreq = nfreq + Atom%phot(ii)%nfreq + 1

        end do

      ! If there is no input
      else

        if (pid.eq.0.and.active) then
          write(umsg,'(A)') '   No photoionizations'
          call verbose
        end if

        ! Allocate dummy structure
        allocate(Atom%phot(1))

      end if



      !
      ! Inelastic collisions
      !

      ! Read number of collisional entries
      read (100,*,err=1100) Atom%ncol

      ! If there are inputs
      if(Atom%ncol.gt.0)then

        ! Allocate database
        allocate(Atom%inelas(Atom%ncol))
        ! Flag for forbidden collisions
        allocate(Atom%fcflag(Atom%nlevel,Atom%nlevel))
        Atom%fcflag = 0
        ! Nullify Tbox pointer
        nullify(Atom%Tbox)

        ! Read number of entries (not the same that number of
        ! collisions)
        read (100,*,err=1100) nclin

        ! Reset index
        icol = 0

        ! Index of temperature block
        Tindex = 0

        ! For each entry
        do ii=1,nclin

          ! Check type of entry
          read (100,*,err=1100) cind

          ! It is a temperature entry
          if (cind.lt.0) then

            ! Read the type of the collisions to come, the type of
            ! interpolation and the number of temperatures in the
            ! input table
            read (100,*,err=1100) col_type
            read (100,*,err=1100) nion
            read (100,*,err=1100) nTmp

            ! Advance the temperature index
            Tindex = Tindex + 1

            ! Change the box
            if (associated(Atom%Tbox)) then
              p_T => Atom%Tbox
              do while (associated(p_T%next))
                p_T => p_T%next
              end do
              allocate(p_T%next)
              p_T => p_T%next
              nullify(p_T%next)
            else
              allocate(Atom%Tbox)
              p_T => Atom%Tbox
              nullify(p_T%next)
            end if

            ! Put data in the box
            p_T%nTmp = nTmp
            p_T%ind = Tindex
            p_T%col_type = col_type
            p_T%nion = nion
            allocate(p_T%temp(nTmp))

            ! Read temperatures
            read (100,*,err=1100) p_T%temp

            ! Determine type of interpolation
            if (nion.eq.2.or.nion.eq.3.or.nion.eq.-2) then
              flin = .True.
            else
              flin = .False.
            endif

            ! Complete the box
            p_T%flin = flin

            ! Nullify
            nullify(p_T)

          ! If it is a collisional entry
          else if (cind.ge.0) then

            ! Run up the index
            icol = icol + 1

            ! Read if forbidden
            read (100,*,err=1100) l3

            ! Read the terms involced in the collision
            read (100,*,err=1100) up,low

            ! up > low
            if (low.gt.up) then
              l2 = up
              up = low
              low = l2
            end if

            ! Store in database
            Atom%inelas(icol)%ind = Tindex
            Atom%inelas(icol)%col_type = cind
            Atom%inelas(icol)%up = up
            Atom%inelas(icol)%low = low
            Atom%inelas(icol)%forbid = l3

            ! Allocate space
            allocate(Atom%inelas(icol)%Cul(nTmp))

            ! Read the data in the table
            read (100,*,err=1100) Atom%inelas(icol)%Cul

            ! If between terms, get stages
            if (col_type.eq.0) then

              ! Stages
              stagl = Atom%stage(low)
              stagu = Atom%stage(up)

            ! If between levels, get degeneracy and flags
            else if (col_type.eq.1) then

              ! Stages
              stagl = Atom%stage(Atom%term(low))
              stagu = Atom%stage(Atom%term(up))

            end if ! If between levels

            !
            ! b-b symmetric
            !

            ! If it is a symmetric b-b excitation
            if (cind.eq.0.or.cind.eq.2.or.cind.eq.3) then

              if (stagu.ne.stagl) then
                umsg = 'Exciting collisions cannot '// &
                       'change ion'
                goto 1200
              end if


            !
            ! b-f
            !
            else if (cind.eq.1.or.cind.eq.4) then

              ! Ionizing collisions can only be level wise
              if (col_type.eq.0) then
                umsg = 'Ionizing collisions only admit '// &
                       'level wise rates'
                goto 1200
              end if

              if (stagu.eq.stagl) then
                umsg = 'Ionizing collisions must '// &
                       'change ion'
                goto 1200
              end if


            !
            ! Charge transfer
            !
            else if (cind.eq.5.or.cind.eq.6) then

              ! Ionizing collisions can only be level wise
              if (col_type.eq.0) then
                umsg = 'Charge transfer collisions only admit '// &
                       'level wise rates'
                goto 1200
              end if

              if (stagu.eq.stagl) then
                umsg = 'Charge transfer collisions must '// &
                       'change ion'
                goto 1200
              end if

            end if ! Type of collision
          end if ! Collisional rate input

        end do ! Collision input lines in text

      ! No collisions
      else

        if (active.and.pid.eq.0) then
          write(umsg,'(A)') '   No collisions'
          call verbose
        end if

      end if

      ! Calculate the atomic part of the Doppler width
      Atom%cDopp = dopp/sqrt(Atom%rmass)

      close (100)

      ! Control that everything went fine
      call control

      ! Remove temporal file
      if (pid.eq.0) then
        call system('rm tmp_atom_'//ID)
      end if

      ! If multilevel, trick the quantum numbers
      if (Atom%ML) then

        Atom%Sval = 0d0

        do iterm=1,Atom%nMulti
          Atom%rLval(iterm) = Atom%rJval(1,iterm)
        end do

      end if

      ! Deallocate unecessary data if not active atom
      if (.not.active) then
        deallocate(Atom%nfreqt)
        deallocate(Atom%nfreqtc)
        deallocate(Atom%lemiss2)
        deallocate(Atom%fflag)
      end if

      ! Control that everything went fine
      call control

      return

1000  umsg = 'Error opening atom file'
      urou = 'rAtom'
      call gaborted
1100  umsg = 'Error reading atom file'
1200  close(100)
      urou = 'rAtom'
      call gaborted

      end subroutine rAtom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generates a 6 levels hydrogen model\n
      !!    Atom(Atom_class): Structure with the atomic data
      subroutine AtomH(Atom)

      ! I/O

      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iterm,iterm1,ii,itran,iJ,nTmp
      double precision:: radAcoeff
      double precision, dimension(:), allocatable:: deg

      ! Pointer
      type(Tbox_class), pointer:: p_T

      ! Initialize
      nullify(p_T)

      ! Informative message
      if (pid.eq.0) then
        umsg = ' - Could not find hydrogen atom, loading '// &
               'default with 6-n levels for the background'
        call verbose
      end if

      ! Atomic ID
      Atom%Element = ' H'
      ! Identify this as a custom atom
      Atom%cust = .True.
      ! Atomic mass
      Atom%rmass = 1.00794d0
      ! Abundance modifier
      Atom%abun_mod = 1d0
      ! Abundance with respect to Hydrogen
      Atom%abun = 1d0

      ! Number of 'terms'
      Atom%nMulti = 6
      ! Number of b-b transitions
      Atom%ntran = 10
      ! Number of b-f transitions
      Atom%nphot = 5
      ! Number of elastic collisions
      Atom%ngk = 0
      ! Maximum of sublevels
      Atom%nJmax = 1
      ! Maximum of levels
      Atom%NNN = Atom%nMulti*Atom%nJmax

      ! Allocations
      ! Orbital angular momentum
      allocate(Atom%rLval(Atom%nMulti))
      ! Spin angular momentum
      allocate(Atom%Sval(Atom%nMulti))
      ! Degeneration (2L+1)*(2S+1)
      allocate(deg(Atom%nMulti))
      ! Number of FS levels
      allocate(Atom%nJ(Atom%nMulti))
      ! Ionization stage
      allocate(Atom%stage(Atom%nMulti))
      ! Angular momentum
      allocate(Atom%rJval(Atom%nJmax,Atom%nMulti))
      Atom%rJval = -1d0
      ! Frequency of FS levels
      allocate(Atom%FSfreq(Atom%nJmax,Atom%nMulti))
      ! Frequency of term
      allocate(Atom%TRfreq(Atom%nMulti))

      ! (g-1)/2, with g statistical weight of n level
      Atom%rLval = (/ .5d0, 3.5d0, 8.5d0, 15.5d0, 24.5d0, 0d0 /)
      ! The spin is artificially put to 0, this atom is given in n
      ! levels
      Atom%Sval = 0d0
      ! Ionization state
      Atom%stage = 1
      Atom%stage(Atom%nMulti) = 2
      ! Degeneration of each level
      deg = (/ 2d0, 8d0, 18d0, 32d0, 50d0, 1d0 /)
      ! Number of sublevels for each level
      Atom%nJ = 1
      ! J = L because S = 0
      Atom%rJval(1,:) = Atom%rLval
      ! Energies
      Atom%TRfreq = (/ 0d0, .82258211d0, .97491219d0, &
                       1.02822766d0, 1.05290508d0, &
                       1.09677617d0 /)
      ! Multilevel
      Atom%ML = .True.
      ! nJ = 1 ==> FSfreq = TRfreq
      Atom%FSfreq(1,:) = Atom%TRfreq

      ! Determine the number of FS levels in this atom
      Atom%nlevel = sum(Atom%nJ)

      ! Define "super-indexes" to handle large number of dimensions
      ! Term given the level
      allocate(Atom%term(Atom%nlevel))
      ! Sublevel in its term given the global level
      allocate(Atom%sublevel(Atom%nlevel))
      ! Degeneration of each term
      allocate(Atom%deg(Atom%nlevel))
      Atom%deg = deg
      ! Allocate indexing by term
      allocate(Atom%irho(Atom%nMulti))

      ! Initialize
      ii = 0

      ! For each term
      do iterm=1,Atom%nMulti

        ! Allocate rho indexing
        allocate(Atom%irho(iterm)%irho_ij(Atom%nJ(iterm)))

        ! For each level
        do iJ=1,Atom%nJ(iterm)

          ii = ii + 1
          Atom%irho(iterm)%irho_ij(iJ) = ii
          Atom%term(ii) = iterm
          Atom%sublevel(ii) = iJ

        end do
      end do



      !
      ! Radiative transitions
      !

      ! Allocations
      ! Energy of the term-term transition
      allocate(Atom%Dfreq(Atom%ntran))
      ! Type of elastic broadening
      allocate(Atom%broad_type(Atom%ntran))
      ! Arguments to calculate the broadening
      allocate(Atom%broad_args(4,Atom%ntran))
      ! Argument for Stark broadening
      allocate(Atom%broad_stark(Atom%ntran))
      ! Doppler widths that the line spawns
      allocate(Atom%dwvl(Atom%ntran))
      ! Dopple widths that the core of the line spawns
      allocate(Atom%dwvlc(Atom%ntran))
      ! Indexing of radiative transitions
      allocate(Atom%irad(Atom%nMulti,Atom%nMulti))
      Atom%irad = 0
      ! Einstein coefficients Aul and Blu
      allocate(Atom%Ecoeff(Atom%nMulti,Atom%nMulti))
      Atom%Ecoeff = 0d0
      ! Nullify Norm
      nullify(Atom%Normp)

      ! Each transition hardwired
      itran = 1
      iterm1 = 2
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = 4.6986d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 600
      Atom%Dwvlc(itran) = 15

      itran = 2
      iterm1 = 3
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .55751d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 10

      itran = 3
      iterm1 = 4
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .12785d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 100
      Atom%Dwvlc(itran) = 3

      itran = 4
      iterm1 = 5
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .04125d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 100
      Atom%Dwvlc(itran) = 3

      itran = 5
      iterm1 = 3
      iterm = 2
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .44101d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 3

      itran = 6
      iterm1 = 4
      iterm = 2
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .084193d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 3

      itran = 7
      iterm1 = 5
      iterm = 2
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .025304d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 3

      itran = 8
      iterm1 = 4
      iterm = 3
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .089860d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 30
      Atom%Dwvlc(itran) = 2

      itran = 9
      iterm1 = 5
      iterm = 3
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .022008d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 30
      Atom%Dwvlc(itran) = 2

      itran = 10
      iterm1 = 5
      iterm = 4
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .026993d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 30
      Atom%Dwvlc(itran) = 1



      !
      ! Photoionizations
      !

      ! Allocations
      ! Indexing of b-f transitions
      allocate(Atom%iphot(Atom%nMulti,Atom%nMulti))
      Atom%iphot = 0
      ! Photoionization information
      allocate(Atom%phot(Atom%nphot))

      ! Each transition hardwired
      ii = 1
      iterm1 = 6
      iterm = 1
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/22.794d0
      Atom%phot(ii)%inalpha(1) = 6.152d-22

      ii = 2
      iterm1 = 6
      iterm = 2
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/91.176d0
      Atom%phot(ii)%inalpha(1) = 1.379d-21

      ii = 3
      iterm1 = 6
      iterm = 3
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/205.147d0
      Atom%phot(ii)%inalpha(1) = 2.149d-21

      ii = 4
      iterm1 = 6
      iterm = 4
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/364.705d0
      Atom%phot(ii)%inalpha(1) = 2.923d-21

      ii = 5
      iterm1 = 6
      iterm = 5
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/569.852d0
      Atom%phot(ii)%inalpha(1) = 3.699d-21


      !
      ! Inelastic collisions
      !

      ! Number of colisional transitions
      Atom%ncol = 10

      ! Allocate inelastic database
      allocate(Atom%inelas(Atom%ncol))

      ! Temperature data points
      nTmp = 6

      ! Allocate Tbox
      allocate(Atom%Tbox)
      p_T => Atom%Tbox
      nullify(p_T%next)

      ! Put data in the box
      p_T%nTmp = nTmp
      p_T%ind = 1
      p_T%col_type = 0
      p_T%nion = 0
      p_T%flin = .False.

      ! Allocate temperature
      allocate(p_T%temp(nTmp))
      p_T%temp = (/ 3d3, 5d3, 7d3, 1d4, 2d4, 3d4 /)

      ii = 1
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 2
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 1.3351d-14,1.0780d-14,9.4856d-15,8.4125d-15, &
               7.0994d-15,6.7550d-15 /)

      ii = 2
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 3
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 8.7453d-16,7.1253d-16,6.3196d-16,5.6633d-16, &
               4.8995d-16,4.7362d-16 /)

      ii = 3
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 4
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 1.6240d-16,1.3263d-16,1.1792d-16,1.0600d-16, &
               9.2277d-17,8.9644d-17 /)

      ii = 4
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 4.7192d-17,3.8580d-17,3.4337d-17,3.0892d-17, &
               2.6995d-17,2.6265d-17 /)

      ii = 5
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 2
      Atom%inelas(ii)%up = 3
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 2.7435d-13,2.5384d-13,2.4973d-13,2.5293d-13, &
               2.7775d-13,2.9945d-13 /)

      ii = 6
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 2
      Atom%inelas(ii)%up = 4
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 1.8623d-14,1.7872d-14,1.8024d-14,1.8705d-14, &
               2.1454d-14,2.3746d-14 /)

      ii = 7
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 2
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 3.5405d-15,3.4405d-15,3.4966d-15,3.6592d-15, &
               4.2698d-15,4.7832d-15 /)

      ii = 8
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 3
      Atom%inelas(ii)%up = 4
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 9.5940d-13,1.0457d-12,1.1455d-12,1.2881d-12, &
               1.6451d-12,1.8677d-12 /)

      ii = 9
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 3
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 6.15d-14,6.8731d-14,7.6113d-14,8.64d-14, &
               1.1348d-13,1.3281d-13 /)

      ii = 10
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 4
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      Atom%inelas(ii)%Cul = &
            (/ 2.7090d-12,3.3113d-12,3.8548d-12,4.5498d-12, &
               6.1112d-12,6.9947d-12 /)

      ! Calculate the atomic part of the Doppler width
      Atom%cDopp = dopp/sqrt(Atom%rmass)

      ! Nullify local pointer
      if (associated(p_T)) nullify(p_T)

      ! Control that everything went fine
      call control

      return

      end subroutine AtomH

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalizes abundances (of the species)\n
      !!    Atom(Atom_class): Structure with the atomic data
      subroutine abund(Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom

      ! Local

      logical, dimension(NA):: renorm, checked

      integer:: i1,i2,nat,natr
      integer, dimension(NA):: list

      double precision:: abmod


      ! Control logical variables
      checked = .False.
      renorm = .False.

      ! Collect the flags into a vector
      do i1=1,NA
        renorm(i1) = Atom(i1)%anorm
      end do

      ! If any of them ask for normalization, normalize
      if (any(renorm)) then

        do i1=1,NA

          ! If we already checked this atom, skip
          if (checked(i1)) cycle

          nat = 1
          natr = 0
          list = 0
          list(1) = i1
          abmod = Atom(i1)%abun_mod

          if (renorm(i1)) natr = natr + 1

          do i2=1,NA

            ! If the atom is checked or if it is the one already
            ! selected, skip
            if (checked(i2).or.i1.eq.i2) cycle

            ! If they are the same element, add to the counter
            if (Atom(i1)%Element.eq.Atom(i2)%Element) then

              ! Add to the list of checked atoms
              nat = nat + 1
              list(nat) = i2

              ! Add this atom to the ones to modify
              if (renorm(i2)) natr = natr + 1

              ! Accumulate the modification factors
              abmod = abmod + Atom(i2)%abun_mod

            end if

          end do

          do i2=1,nat
            checked(list(i2)) = .True.
          end do

          ! If there was an atom asking for renorm
          if (natr.gt.0) then

            if (pid.eq.0) then
              umsg = ' - Atom '//trim(Atom(i1)%Element)// &
                     ' abundance modifier normalized'
              call verbose
            end if

            ! Normalize the atoms on the list
            do i2=1,nat

              Atom(list(i2))%abun_mod = Atom(list(i2))%abun_mod/abmod

            end do

          end if

        end do

      end if


      end subroutine abund

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes fine structure Einstein coefficients from the term
      !! wise one.\n
      !!    Atom(Atom_class): Structure with the atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials and signs
      subroutine setFScoeff(Atom,Flgsg)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg

      ! Local

      integer:: iterm,iterm1,iJ,iJ1,ilevel,ilevel1,itran,ftran

      double precision:: W6,S,rL,rL1,rJ,rJ1,Aul,AJul,BJlu,Dfreq

      ! If the atom was ML, ignore quantum numbers
      if (Atom%ML) then

        ! For each term pair
        do ilevel=1,Atom%nMulti-1
          do ilevel1=ilevel+1,Atom%nMulti

            ! Check for a transition
            itran = Atom%irad(ilevel,ilevel1)

            if (itran.lt.1) cycle

            rJ = Atom%rJval(1,ilevel)
            rJ1 = Atom%rJval(1,ilevel1)

            Dfreq = Atom%FSfreq(1,ilevel1) - &
                    Atom%FSfreq(1,ilevel)

            ! Aul
            AJul = Atom%Ecoeff(ilevel1,ilevel)

            ! Calculate the corresponding FS-Blu
            BJlu = AJul*(2d0*rJ1+1d0)/(2d0*rJ+1d0)/ &
                   (ConvF*1d21*(2d0*c)*Dfreq**3d0)

            Atom%fst(itran)%Aul(1,1) = AJul
            Atom%fst(itran)%Blu(1,1) = BJlu

          end do ! Upper level
        end do ! Lower level

      ! If multiterm, do the branching
      else

        ! For each term pair
        do iterm=1,Atom%nMulti-1
          do iterm1=iterm+1,Atom%nMulti

            ! Check for a transition
            itran = Atom%irad(iterm,iterm1)

            if (itran.lt.1) cycle

            rL = Atom%rLval(iterm)
            rL1 = Atom%rLval(iterm1)
            S = Atom%Sval(iterm)
            Aul = Atom%Ecoeff(iterm1,iterm)

            ! For each pair of FS levels
            do iJ=1,Atom%nJ(iterm)
              do iJ1=1,Atom%nJ(iterm1)

                ! Check for a transition
                ftran = Atom%fst(itran)%irad(iJ1,iJ)

                if (ftran.lt.1) cycle

                rJ = Atom%rJval(iJ,iterm)
                rJ1 = Atom%rJval(iJ1,iterm1)

                Dfreq = Atom%FSfreq(iJ1,iterm1) - &
                        Atom%FSfreq(iJ,iterm)

                ! Branch the Aul
                W6 = fun6j(rL1,rL,1d0,rJ,rJ1,S,Flgsg)
                AJul = (2d0*rL1+1d0)*(2d0*rJ+1d0)*W6*W6*Aul

                ! Calculate the corresponding FS-Blu
                BJlu = AJul*(2d0*rJ1+1d0)/(2d0*rJ+1d0)/ &
                       (ConvF*1d21*(2d0*c)*Dfreq**3d0)

                Atom%fst(itran)%Aul(iJ1,iJ) = AJul
                Atom%fst(itran)%Blu(iJ,iJ1) = BJlu

              end do ! iJ1
            end do ! iJ
          end do ! iterm1
        end do ! iterm

      end if ! ML or MT

      return

      end subroutine setFScoeff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute some derived quantities for the LTE line\n
      !!    Atom(Atom_class): Structure with the atomic data for
      !!                      passive atoms
      !! line(LTEline_class): Structure with the LTE line data
      subroutine setup_LTE_transition(Atom,line)

      ! I/O
      type(Atom_class), dimension(:), intent(inout):: Atom
      type(LTEline_class), intent(inout):: line


      ! If passive atom model
      if (line%is_passive) then

        ! Get from model
        line%rmass = Atom(line%ia)%rmass
        line%abund = Atom(line%ia)%abun
        line%cDopp = Atom(line%ia)%cDopp

      ! No model
      else

        ! Get tabulated
        line%rmass = recallmass_ind(line%ele)
        line%abund = recallabund_ind(line%ele)
        line%cDopp = dopp/sqrt(line%rmass)

      end if

      ! Transition freq. and Einstein c.
      line%Dfreq = line%Eu - line%El
      line%Blu = line%Aul*(2d0*line%Ju + 1d0)/ &
                          (2d0*line%Jl + 1d0)/ &
                 (ConvF*line%Dfreq*1d21*(2d0*c)* &
                  line%Dfreq**2d0)

      ! Number of M levels indexes
      line%nMu = nint(2d0*line%Ju+1d0)
      line%nMl = nint(2d0*line%Jl+1d0)

      ! Check frequencies
      if (line%nfreq.gt.0) then
        if (mod(line%nfreq,2).eq.0) line%nfreq = line%nfreq + 1
        if (mod(line%nfreqc,2).eq.0) line%nfreqc = line%nfreqc + 1
      end if

      ! Add frequencies to global count
      nfreq = nfreq + line%nfreq

      end subroutine setup_LTE_transition

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks if a LTE transition is in the model and removes it\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!  line(LTEline_class): Structure with the LTE line data
      subroutine remove_LTE_transition(Atom,line)

      ! I/O
      type(Atom_class), intent(inout):: Atom
      type(LTEline_class), intent(in):: line

      ! Local
      logical:: no_found_l, no_found_u

      integer:: iterm,iJ,Jl2,Ju2,iterml,iJl,iJu,itermu

      double precision:: El,Eu,factl,factu

      !
      ! Find levels
      !

      ! Initialize
      no_found_l = .True.
      no_found_u = .True.
      Jl2 = nint(line%Jl*2d0)
      Ju2 = nint(line%Ju*2d0)
      El = line%El
      Eu = line%Eu

      ! Relative factors
      if (abs(El).gt.0d0) then
        factl = 1d0/abs(El)
      else
        factl = 1d0
      end if
      if (abs(Eu).gt.0d0) then
        factu = 1d0/abs(Eu)
      else
        factu = 1d0
      end if

      ! For each term
      do iterm=1,Atom%nMulti

        ! If not even same stage
        if (Atom%stage(iterm).ne.line%stage) cycle

        ! For each level
        do iJ=1,Atom%nJ(iterm)

          ! If not found lower
          if (no_found_l) then

            ! Check J
            if (nint(Atom%rJval(iJ,iterm)*2d0).eq.Jl2) then

              ! Compare
              if (abs(El - Atom%FSfreq(iJ,iterm))*factl.le. &
                  TINYSP) then

                  ! Found
                  no_found_l = .False.
                  iterml = iterm
                  iJl = iJ

              end if ! Same energy
            end if ! Same J
          end if ! Not found l

          ! If not found upper
          if (no_found_u) then

            ! Check J
            if (nint(Atom%rJval(iJ,iterm)*2d0).eq.Ju2) then

              ! Compare
              if (abs(Eu - Atom%FSfreq(iJ,iterm))*factu.le. &
                  TINYSP) then

                  ! Found
                  no_found_u = .False.
                  itermu = iterm
                  iJu = iJ

              end if ! Same energy
            end if ! Same J
          end if ! Not found u

        end do ! Levels

        ! If found, stop
        if (.not.no_found_l.and..not.no_found_u) exit

      end do ! Terms

      ! If not found, the line is not here
      if (no_found_l.or.no_found_u) return

      ! Multi-level
      if (Atom%ML) then

        ! Just remove the transition (from the calculations)
        Atom%irad(itermu,iterml) = 0
        Atom%irad(iterml,itermu) = 0

      ! Multi-term
      else

        ! Remove if not forbidden
        if (abs(Jl2-Ju2).le.1.and.Jl2+Ju2.gt.0) then
          Atom%irad(itermu,iterml) = 0
          Atom%irad(iterml,itermu) = 0
        end if

      end if

      return

      end subroutine remove_LTE_transition

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates array of Atom_class.\n
      !!  Atom(Atom_class): Structure to allocate\n
      !!       nn(integer): Size to allocate
      subroutine allocateatom(Atom,nn)

      ! I/O
      type(Atom_class), dimension(:), allocatable:: Atom
      integer, intent(in):: nn

      allocate(Atom(nn))

      ! Control that everything went fine
      call control

      return

      end subroutine allocateatom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up labels for atomic files.\n
      !!  Atom(Atom_class): Structure to allocate\n
      !!       nn(integer): Size to allocate
      subroutine set_atom_label(Atom,nn)

      ! I/O
      type(Atom_class), dimension(:), allocatable:: Atom
      integer, intent(in):: nn

      ! Local
      character(len=10):: label

      logical:: repeated
      logical, dimension(nn):: check

      integer:: ia,ja,ka,irep


      ! Each atom, create default label
      do ia=1,nn

        ! Generate label
        Atom(ia)%file_label = '          '
        if (Atom(ia)%Element(1:1).eq.' ') then
          Atom(ia)%file_label(1:1) = Atom(ia)%Element(2:2)
        else
          Atom(ia)%file_label = Atom(ia)%Element
        end if
      end do

      ! Initialize
      check = .False.

      ! Now check repeated
      do ia=1,nn-1

        ! Checked?
        if (check(ia)) cycle

        ! Current label
        label = Atom(ia)%file_label

        ! Initialized
        repeated = .False.
        irep = 0

        ! Check repeated
        do ja=ia+1,nn

          ! Repeated?
          if (trim(Atom(ja)%file_label).eq.trim(label)) then
            repeated = .True.
            exit
          end if

        end do

        ! If repeated
        if (repeated) then

          ! Change current and initialize index
          Atom(ia)%file_label = trim(label)//'-1'
          ka = 1

          ! Check repeated
          do ja=ia+1,nn

            ! Repeated?
            if (trim(Atom(ja)%file_label).eq.trim(label)) then

              ! Advance
              ka = ka + 1

              ! Get complement
              if (ka.lt.10) then
                write(Atom(ja)%file_label,'("-",i1)') ka
              else if (ka.lt.100) then
                write(Atom(ja)%file_label,'("-",i2)') ka
              else if (ka.lt.1000) then
                write(Atom(ja)%file_label,'("-",i3)') ka
              end if

              ! Get label
              Atom(ja)%file_label = trim(label)// &
                                    trim(Atom(ja)%file_label)
              ! Flag
              check(ja) = .True.

            end if

          end do

        end if ! Repeated

        ! Flag current
        check(ia) = .True.

      end do ! Atoms

      return

      end subroutine set_atom_label

!#####################################################################
!#####################################################################
!#####################################################################

      !> Shift to the right the list of atoms and file names.\n
      !!   Atom(Atom_class): Array to shift\n
      !!  files(Atom_class): Array to shift
      subroutine shiftatoms(Atom,files)

      ! I/O
      type(Atom_class), dimension(:), allocatable:: Atom
      type(strarr_class), dimension(:), allocatable:: files

      ! Local
      type(Atom_class), dimension(:), allocatable:: Atomaux
      type(strarr_class), dimension(:), allocatable:: filesaux


      ! If there are no background atoms, no need to shift
      if (nAb.lt.1) then

          nAb = 1
          allocate(Atom(nAb))
          allocate(files(nAb))
          files(1)%str = 'N'

      ! If we have to shift the loaded atoms
      else

        ! Copy data into auxiliar
        allocate(Atomaux(nAb))
        allocate(filesaux(nAb))
        Atomaux = Atom
        filesaux = files
        deallocate(Atom)
        deallocate(files)

        nAb = nAb + 1

        ! Create a larger one and copy back
        allocate(Atom(nAb))
        allocate(files(nAb))
        Atom(2:nAb) = Atomaux
        files(2:nAb) = filesaux
        files(1)%str = 'N'

        ! Deallocate auxiliar
        deallocate(Atomaux)
        deallocate(filesaux)

      end if ! Previous inputs

      ! Check if everything is fine
      call control

      return

      end subroutine shiftatoms

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates height dependent population\n
      !!  Atom(Atom_class): Atom structures\n
      !!       nn(integer): Size of the structure
      subroutine prepareatom(Atom,nn)

      ! I/O
      type(Atom_class), dimension(:), intent(inout):: Atom
      integer, intent(in):: nn

      ! Local
      integer:: ia

      ! For each atom
      do ia=1,nn

        ! Allocate total population
        allocate(Atom(ia)%n(nz))

      end do

      end subroutine prepareatom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates variables that depend on the number of processes\n
      !!  Atom(Atom_class): Atom structures
      subroutine prepareatomMPI(Atom)

      ! I/O
      type(Atom_class), dimension(:), intent(inout):: Atom

      ! Local
      integer:: ia,ii

      ! If master of MPI
      if (pid.eq.0.and.nproc.gt.1) then

        ! For each atom
        do ia=1,na

          ! Lower frequency limits for Master
          allocate(Atom(ia)%Mif0(Atom(ia)%ntran,0:nproc-1))
          ! Upper frequency limits for Master
          allocate(Atom(ia)%Mif1(Atom(ia)%ntran,0:nproc-1))
          ! Lower frequency weight for Master
          allocate(Atom(ia)%MW0(Atom(ia)%ntran,0:nproc-1))
          ! Upper frequency weight for Master
          allocate(Atom(ia)%MW1(Atom(ia)%ntran,0:nproc-1))

          ! For each b-f transition
          do ii=1,Atom(ia)%nphot
            ! Allocate limits for master
            allocate(Atom(ia)%phot(ii)%Mif0(0:nproc-1))
            allocate(Atom(ia)%phot(ii)%Mif1(0:nproc-1))
            ! Allocate weights for master
            allocate(Atom(ia)%phot(ii)%MW0(0:nproc-1))
            allocate(Atom(ia)%phot(ii)%MW1(0:nproc-1))
          end do
        end do

      end if ! Master of MPI

      end subroutine prepareatomMPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the size of the M blocks\n
      !!    Atom(Atom_class): Atom structures
      subroutine initMblock(Atom)

      ! I/O
      type(Atom_class), dimension(:), intent(inout):: Atom

      ! Local
      integer:: ia,nM,iM,iJ,i,iterm
      double precision:: rL,S,rJ,rJmin,rJmax,pS,rM,rJm


      ! For each atom
      do ia=1,nA

        ! For each term
        do iterm=1,Atom(ia)%nMulti

          ! Get term quantities
          rL = Atom(ia)%rLval(iterm)
          S = Atom(ia)%Sval(iterm)

          ! Multilevel
          if (Atom(ia)%ML) then

            rJ = Atom(ia)%rJval(1,iterm)
            nM = nint(2d0*rJ + 1d0)

            ! For each M
            do iM=1,nM

              ! Size of the block of this M (just one)
              Atom(ia)%nblk(iM,iterm) = 1

            end do

            ! Rest of empty space
            do iM=nM+1,Atom(ia)%nMmax

              ! Size of the block of this M (just one)
              Atom(ia)%nblk(iM,iterm) = 0

            end do

          ! Multiterm
          else

            rJmin = abs(rL - S)
            rJmax = rL + S

            pS = S*(S+1d0)*(2d0*S + 1d0)

            nM = nint(2d0*rJmax + 1d0)

            ! Run over the magnetic block
            do iM=1,nM

              rM = -rJmax + dble(iM-1)

              rJm = max(abs(rM),rJmin)

              ! initialize the column index
              i=0

              ! Column loop
              do iJ = 1,Atom(ia)%nJ(iterm)

                rJ = Atom(ia)%rJval(iJ,iterm)

                if (rJ.ge.rJm) then

                  ! Increase the column index
                  i = i + 1

                end if ! Test J compatibility with M

              end do ! iJ

              ! Block dimension
              Atom(ia)%nblk(iM,iterm) = i

            end do ! iM

            ! Rest of empty space
            do iM=nM+1,Atom(ia)%nMmax

              ! Size of the block of this M (just one)
              Atom(ia)%nblk(iM,iterm) = 0

            end do

          end if ! Multilevel or multiterm

        end do ! Terms
      end do ! Atoms

      end subroutine initMblock

!#####################################################################
!#####################################################################
!#####################################################################

      end module ratom_mod
