      !> Reading settings
      module rinput_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!  Start:
!     17/04/2017
!  Last version:
!     24/09/2025 V4.0.11
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     24/09/2025:   V4.0.11 - Read Input%Const_regul (TdPA)
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
!  rInput
!    Read the input file with the configuration of the run
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the input file with the configuration of the run\n
      !!  Input(Input_class): Structure with configuration data
      subroutine rInput(Input)

      ! I/O

      type(Input_class), intent(inout):: Input

      ! Local

      character(LEN=1):: cdump
      character(LEN=2):: cdump2
      character(LEN=5):: CPUC

      integer:: ios,ia,i1,i2,i3,i4,i5

      double precision:: ddump,ddump2

      ! Routine name
      urou = 'rInput'


      !
      ! Translate input with python
      !

      if(pid.eq.0) call system('python '//trim(Input%source)// &
                               'rinput.py '//trim(Input%input)// &
                               ' '//Input%ID)

      ! Wait until master is finished with the input file
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !
      ! Read the translated input file
      !

      ! Open temporal file
      open(100, file='tmp_input_'//Input%ID, status='old', &
           iostat=ios, err=1000)

      ! Success
      read(100,*,err=1100) ios

      ! If no correct file, abort
      if (ios.lt.0) then

        umsg = 'Problem translating the input file'
        goto 1200

      end if

      ! Type of run
      read(100,*,err=1100) Input%run_mode
      run_mode = Input%run_mode

      ! Atmospheric file
      read(100,'(A)',err=1100) Input%atmo

      ! Atmospheric scale
      read(100,'(A)',err=1100) Input%atm_scale
      read(100,*,err=1100) Input%omega_ref

      ! Atmosphere pressure/density scale
      read(100,*,err=1100) Input%atmo_char

      ! Atom fix population
      read(100,'(A)',err=1100) cdump
      Input%respect_zalt = cdump.eq.'Y'

      ! Number of atoms
      read(100,*,err=1100) Input%nA
      nA = Input%nA

      ! If there are atoms to allocate
      if (nA.gt.0) then

        ! Allocations
        ! Filename of atomic models
        allocate(Input%atom(nA))
        ! Filename of atomic populations
        allocate(Input%popu(nA))
        ! Keep populations fixed
        allocate(Input%fixp(nA))
        ! Keep populations fixed for lower terms
        allocate(Input%fixplt(nA))
        ! Zero ion
        allocate(Input%zero_ion(nA))
        ! No wavelengths
        allocate(Input%skip_wave(nA))
        ! Ionization files for CLE
        if (Input%run_mode.eq.2) allocate(Input%ionf(nA))

        ! Read atomic file, population file names, and if fixing
        ! populations
        do ia=1,nA

          ! Add to RAM (fixp, zero_ion, atom, and fixplt are
          ! deallocated before any relevant call)
          MRAMc = MRAMc + 1d-6*(sizeof(Input%popu(ia)) + &
                                sizeof(Input%skip_wave(ia)))
         if (Input%run_mode.eq.2) &
             MRAMc = MRAMc + 1d-6*sizeof(Input%ionf)

          ! Atom file model
          read(100,'(A)',err=1100) Input%atom(ia)
          ! Atom population file
          read(100,'(A)',err=1100) Input%popu(ia)

          ! Atom fix population
          read(100,'(A)',err=1100) cdump
          Input%fixp(ia) = cdump.eq.'F'

          ! Atom zero ion
          read(100,'(A)',err=1100) cdump
          Input%zero_ion(ia) = cdump.eq.'F'

          ! Ionization files for CLE
          if (Input%run_mode.eq.2) then

            ! Read type of ionization specification
            read(100,*,err=1100) Input%ionf(ia)%typ

            ! If file
            if (Input%ionf(ia)%typ.eq.0) then

              ! Read filename
              read(100,'(A)',err=1100) Input%ionf(ia)%str

            ! If number
            else if (Input%ionf(ia)%typ.eq.1) then

              ! Read value
              read(100,*,err=1100) Input%ionf(ia)%val

            end if ! Type of input

          ! No CLE
          else

            ! Read dummy -1
            read(100,*,err=1100) ios

          end if ! Value

          ! Atom skip wavelengths
          read(100,'(A)',err=1100) cdump
          Input%skip_wave(ia) = .not.(cdump.eq.'N')

          ! Atom fix terms populations
          read(100,'(A)',err=1100) cdump
          Input%fixplt(ia) = cdump.eq.'F'

          ! Check cascades of conditions for fixing populations
          if (Input%fixp(ia)) Input%fixplt(ia) = .False.

        end do

      end if ! Atoms to read

      ! Number of background atoms
      read(100,*,err=1100) Input%nAb
      nAb = Input%nAb

      ! If larger than 0
      if (Input%nAb.gt.0) then

        ! Allocate

        ! Filename of atomic models for background
        allocate(Input%atomback(Input%nAb))
        ! Filename of population file for background
        allocate(Input%popuback(Input%nAb))

        ! Read atomic file and population file names
        do ia=1,Input%nAb

          ! Add to RAM counter (atomback is deallocated before
          ! any relevant call)
          MRAMc = MRAMc + 1d-6*sizeof(Input%popuback(ia))
          read(100,'(A)',err=1100) Input%atomback(ia)
          read(100,'(A)',err=1100) Input%popuback(ia)

        end do ! Passive atoms

      end if ! There are background atoms

      ! Number of molecules
      read(100,*,err=1100) Input%nM
      nM = Input%nM

      ! If larger than 0
      if (Input%nM.gt.0) then

        ! Allocate molecule file name (it is deallocated before
        ! any relevant call and thus not added to RAM counter)
        allocate(Input%mol(Input%nM))

        ! Read molecule file names
        do ia=1,Input%nM
          read(100,'(A)',err=1100) Input%mol(ia)
        end do

      end if ! There are molecules

      ! Type of Magnetic field input
      read(100,'(A)',err=1100) cdump

      ! If file
      if (cdump.eq.'F') then

        ! Flag no number and get filename
        Input%bfieldn = .False.
        read(100,'(A)',err=1100) Input%bfield

      ! If number
      else

        ! Flag numeric and get numbers
        Input%bfieldn = .True.
        read(100,*,err=1100) Input%bfieldv(1)
        read(100,*,err=1100) Input%bfieldv(2)
        read(100,*,err=1100) Input%bfieldv(3)

      end if ! Type of magnetic field input

      ! File with fudge factors
      read(100,'(A)',err=1100) Input%fudge

      ! File with input spectra
      read(100,'(A)',err=1100) Input%spect_input

      ! Use Allen tabulation for intensity
      read(100,*,err=1100) cdump
      Input%use_allen = cdump.eq.'Y'

      ! Type of CLV
      read(100,*,err=1100) cdump
      if (cdump.eq.'A') then
        Input%clv_type = 0
      else
        Input%clv_type = -1
      end if

      ! Assume fully flat spectrum at input (if no input spectra)
      read(100,*,err=1100) cdump
      Input%flat_cle_in = cdump.eq.'Y'

      ! Path to CHIANTI database
      read(100,'(A)',err=1100) Input%chianti_path

      ! T_rad
      read(100,*,err=1100) Input%T_rad

      ! R_star
      read(100,*,err=1100) Input%R_star

      ! Neglect continuum in CLE
      read(100,*,err=1100) cdump
      Input%add_cont_cle = .not.(cdump.eq.'Y')

      ! File with partition functions
      read(100,'(A)',err=1100) Input%pf

      ! File with abundances
      read(100,'(A)',err=1100) Input%abund

      ! File with Barklem SP data
      read(100,'(A)',err=1100) Input%bark_sp

      ! File with Barklem PD data
      read(100,'(A)',err=1100) Input%bark_pd

      ! File with Barklem DF data
      read(100,'(A)',err=1100) Input%bark_df

      ! Ignore b-b background transitions
      read(100,*,err=1100) cdump
      Input%addbb = .not.(cdump.eq.'Y')

      ! Read number of Kurucz files
      read(100,*,err=1100) Input%NK

      ! If Kurucz files
      if (Input%NK.ge.1) then

        ! Allocate Kurucz file names
        allocate(Input%kurucz(Input%NK))

        ! Read Kurucz file names
        do ia=1,Input%NK
          MRAMc = MRAMc + 1d-6*sizeof(Input%kurucz(ia))
          read(100,'(A)',err=1100) Input%kurucz(ia)
        end do

      end if ! Kurucz files

      ! Read number of LTE lines
      read(100,*,err=1100) Input%nLTE
      nLTEl = Input%nLTE

      ! If LTE lines
      if (Input%nLTE.ge.1) then

        ! Allocate structure
        allocate(Input%LTEline(Input%nLTE))

        ! For each entry
        do ia=1,Input%nLTE

          ! Initialize
          Input%LTEline(ia)%is_passive = .False.
          Input%LTEline(ia)%taulim_l = .False.
          Input%LTEline(ia)%Tlim_l = .False.
          Input%LTEline(ia)%ia = -1
          Input%LTEline(ia)%ia = -1

          ! Memory count
          MRAMc = MRAMc + 1d-6*sizeof(Input%LTEline(ia))

          ! Read type of atom label
          read(100,*,err=1100) i1

          ! If label
          if (i1.eq.0) then

            ! Read label
            read(100,*,err=1100) cdump2

            ! Get index
            Input%LTEline(ia)%ele = atom_char2index(cdump2)

          ! If number
          else if (i1.eq.1) then

            ! Read index
            read(100,*,err=1100) Input%LTEline(ia)%ele

          end if ! Type of atomic label

          ! Stage
          read(100,*,err=1100) Input%LTEline(ia)%stage

          ! Lower level
          read(100,*,err=1100) Input%LTEline(ia)%el
          read(100,*,err=1100) Input%LTEline(ia)%Jl
          read(100,*,err=1100) Input%LTEline(ia)%gl

          ! Upper level
          read(100,*,err=1100) Input%LTEline(ia)%eu
          read(100,*,err=1100) Input%LTEline(ia)%Ju
          read(100,*,err=1100) Input%LTEline(ia)%gu

          ! Transition
          read(100,*,err=1100) Input%LTEline(ia)%Aul
          read(100,*,err=1100) Input%LTEline(ia)%broad_type
          allocate(Input%LTEline(ia)%broad_args(4))
          MRAMc = MRAMc + 1d-6*sizeof(Input%LTEline(ia)%broad_args)
          read(100,*,err=1100) Input%LTEline(ia)%broad_args
          read(100,*,err=1100) Input%LTEline(ia)%broad_stark
          read(100,*,err=1100) Input%LTEline(ia)%f_c
          read(100,*,err=1100) Input%LTEline(ia)%broad_rad
          read(100,*,err=1100) Input%LTEline(ia)%nfreq
          read(100,*,err=1100) Input%LTEline(ia)%nfreqc
          read(100,*,err=1100) Input%LTEline(ia)%Dwvl
          read(100,*,err=1100) Input%LTEline(ia)%Dwvlc
          read(100,*,err=1100) Input%LTEline(ia)%ltype

          ! Limits
          read(100,*,err=1100) Input%LTEline(ia)%taulim
          read(100,*,err=1100) Input%LTEline(ia)%Tlim

          ! No wavelengths
          read(100,*,err=1100) cdump
          Input%LTEline(ia)%nowave = cdump.eq.'T'

        end do ! Entries

      end if ! LTE lines

      ! Read number of wavelength files
      read(100,*,err=1100) Input%NW

      ! If Wavelength files
      if (Input%NW.ge.1) then

        ! Allocate Wavelength file names
        allocate(Input%waves(Input%NW))

        ! Read Wavelength file names
        do ia=1,Input%NW
          MRAMc = MRAMc + 1d-6*sizeof(Input%waves(ia))
          read(100,'(A)',err=1100) Input%waves(ia)
        end do

      end if ! Wavelength files

      !
      ! Sanity check wavelengths
      !

      ! If there are atoms
      if (nA.gt.0) then

        ! If no wavelengths
        if (all(Input%skip_wave).and.Input%NW.eq.0) then

          ! Error
          umsg = 'Cannot neglect automatically generated '// &
                 'wavelengths without including LTE lines or '// &
                 'wavelength files'
          call gaborted

        end if ! No wavelengths
      end if ! Active atoms

      ! Asymmetry inputs
      read(100,*,err=1100) Input%nasym

      ! If there are inputs
      if (Input%nasym.gt.0) then

        ! Number of asymmetry inputs as numbers
        read(100,*,err=1100) Input%nasym_num

        ! If numerical values, allocate space
        if (Input%nasym_num.gt.0) then
          allocate(Input%asym_num(2,Input%nasym_num))
          MRAMc = MRAMc + 1d-6*sizeof(Input%asym_num)
        end if

        ! Number of asymmetry inputs as files
        read(100,*,err=1100) Input%nasym_fil

        ! If files, allocate space
        if (Input%nasym_fil.gt.0) then
          allocate(Input%asym_fil(Input%nasym_fil))
          MRAMc = MRAMc + 1d-6*sizeof(Input%asym_fil(ia))
        end if

        ! Initialize indexes for numerical and files
        i2 = 0
        i3 = 0

        ! For each entry
        do i1=1,Input%nasym

          ! Read character
          read(100,'(A)',err=1100) cdump

          ! If value
          if (cdump.eq.'V') then

            ! If value, read and store
            i2 = i2 + 1
            read(100,*,err=1100) i4, i5, ddump, ddump2
            Input%asym_num(1,i2) = dcmplx(i4,i5)
            Input%asym_num(2,i2) = dcmplx(ddump,ddump2)

          end if ! Value

          ! If file
          if (cdump.eq.'F') then

            ! Read
            i3 = i3 + 1
            read(100,'(A)',err=1100) Input%asym_fil(i3)

          end if ! File

        end do ! Entries

      end if ! Asymmetry inputs

      ! Force asymm
      read(100,*,err=1100) cdump
      force_asym = cdump.eq.'Y'

      ! Stimulated emission
      read(100,*,err=1100) cdump
      stm = cdump.eq.'Y'

      ! Output folder
      read(100,'(A)',err=1100) Input%folder

      ! Cache file
      Input%cache = trim(Input%folder)//'/cache'

      ! Directional nodes for formal solver
      read(100,*,err=1100) Input%nTh
      read(100,*,err=1100) Input%nPh
      read(100,*,err=1100) Input%nThI
      read(100,*,err=1100) Input%nPhI

      ! Directions for emergent solution
      ! Polar
      read(100,*,err=1100) Input%nThLOS
      ! If polar directions
      if (Input%nThLOS.gt.0) then
        allocate(Input%L_mu(Input%nThLOS))
        MRAMc = MRAMc + 1d-6*sizeof(Input%L_mu)
        read(100,*,err=1100) Input%L_mu(1:Input%nThLOS)
      end if
      ! If azimuthal directions
      read(100,*,err=1100) Input%nPhLOS
      if (Input%nPhLOS.gt.0) then
        allocate(Input%L_phi(Input%nPhLOS))
        MRAMc = MRAMc + 1d-6*sizeof(Input%L_phi)
        read(100,*,err=1100) Input%L_phi(1:Input%nPhLOS)
      end if

      ! Sanity check
      if (Input%nThLOS.le.0.or.Input%nPhLOS.le.0) then
        Input%nThLOS = 0
        Input%nPhLOS = 0
      end if

      ! Mode of solution
      read(100,*,err=1100) Input%mode

      ! Force type of calculation
      read(100,*,err=1100) Input%force

      ! Two-step velocity for intensity
      read(100,*,err=1100) cdump
      Input%two_step_I_v = cdump.eq.'Y'

      ! Restrict problem in T bottom strictly
      read(100,*,err=1100) cdump
      Input%rest_Tlo_strc = cdump.eq.'Y'

      ! Restrict problem in T bottom
      read(100,*,err=1100) cdump
      Input%rest_Tlo = cdump.eq.'Y'

      ! Get range if they exist
      if (Input%rest_Tlo) then
        read(100,*,err=1100) Input%r1tt
      end if

      ! Restrict problem in T up strictly
      read(100,*,err=1100) cdump
      Input%rest_Tup_strc = cdump.eq.'Y'

      ! Restrict problem in T up
      read(100,*,err=1100) cdump
      Input%rest_Tup = cdump.eq.'Y'

      ! Get range if they exist
      if (Input%rest_Tup) then
        read(100,*,err=1100) Input%r0tt
      end if

      ! Restrict problem in tau_c strictly
      read(100,*,err=1100) cdump
      Input%rest_tau_strc = cdump.eq.'Y'

      ! Restrict problem in tau_c
      read(100,*,err=1100) cdump
      Input%rest_tau = cdump.eq.'Y'

      ! Get ranges if they exist
      if (Input%rest_tau) then
        read(100,*,err=1100) Input%r0tc
        read(100,*,err=1100) Input%r1tc
      end if

      ! Restrict problem in height strictly
      read(100,*,err=1100) cdump
      Input%rest_z_strc = cdump.eq.'Y'

      ! Restrict problem in height
      read(100,*,err=1100) cdump
      Input%rest_z = cdump.eq.'Y'

      ! Get ranges if they exist
      if (Input%rest_z) then
        read(100,*,err=1100) Input%r0z
        read(100,*,err=1100) Input%r1z
      end if

      ! Mode of Zeeman effect
      read(100,*,err=1100) Input%zeeman_mode

      ! Correction of JKQ between multi-level and multi-term
      read(100,*,err=1100) cdump
      Input%Pcorr = cdump.eq.'Y'

      ! Forbidden collisions can transfer alignment
      read(100,*,err=1100) cdump
      fcol_transfer = cdump.eq.'Y'

      ! MIT transitions
      read(100,*,err=1100) Input%MIT_input

      ! MIT node/factor
      read(100,*,err=1100) Input%MIT_node
      Input%MIT_node = 1d0/Input%MIT_node

      ! Skip saving solution
      read(100,*,err=1100) cdump
      Input%keep_sol = .not.(cdump.eq.'Y')

      ! Input solution file
      read(100,'(A)',err=1100) Input%solution

      ! Dummy line
      read(100,*,err=1100) cdump

      ! Keep intensity solution
      read(100,*,err=1100) cdump
      Input%keepIsol = cdump.eq.'Y'

      ! Keep Stokes in solution
      read(100,*,err=1100) cdump
      KSTK = cdump.eq.'Y'

      ! Keep Stokes in solution
      read(100,*,err=1100) cdump
      Input%anisotropy_only = cdump.eq.'Y'

      ! Cut rhoKQ multipole orders
      read(100,*,err=1100) Kcut

      ! Sanity check
      if (Kcut.lt.0) Kcut = 1000

      ! Cut rhoKQ multipole orders in absorb1
      read(100,*,err=1100) cdump
      KcutAB = cdump.eq.'Y'

      ! K cut term-wise entries
      read(100,*,err=1100) ios

      ! If there are entries
      if (ios.gt.0) then

        ! Allocate entries
        allocate(Input%Kcut_input(4,ios))
        MRAMc = MRAMc + 1d-6*sizeof(Input%Kcut_input)

        ! And read each entry
        do ia=1,ios
          read(100,*,err=1100) Input%Kcut_input(:,ia)
        end do

      end if ! K cut term-wise entries

      ! Cut JKQ multipole orders
      read(100,*,err=1100) Krad

      ! Sanity check
      if (Krad.lt.0.or.Krad.gt.2) Krad = 2

      ! Initialize K limit for transitions
      Kradl = Krad

      ! Store J symbols in memoization
      read(100,*,err=1100) cdump
      Input%memo = cdump.eq.'Y'

      ! Store photoionization quantities in RAM
      read(100,*,err=1100) cdump
      PIRAM = cdump.eq.'Y'

      ! Type of intensity Voigt profile
      read(100,*,err=1100) VOITY

      ! Store Voigt profiles for intensity in RAM
      read(100,*,err=1100) cdump
      VIRAM = cdump.eq.'Y'

      ! Store Voigt profiles for intensity in RAM for LTE lines
      read(100,*,err=1100) cdump
      LVIRAM = cdump.eq.'Y'

      ! Do not store for LTE lines if not storing for active lines
      if (.not.VIRAM.and.LVIRAM) LVIRAM = .False.

      ! Store Voigt profiles for polarization in RAM
      read(100,*,err=1100) cdump
      VPRAM = cdump.eq.'Y'

      ! Store Voigt profiles for polarization in RAM for LTE lines
      read(100,*,err=1100) cdump
      LVPRAM = cdump.eq.'Y'

      ! Do not store for LTE lines if not storing for active lines
      if (.not.VPRAM.and.LVPRAM) LVPRAM = .False.

      ! Limit in RAM for profiles
      read(100,*,err=1100) RLIM

      ! Arbitrary limit
      if (RLIM.lt.0) RLIM = 1000000000

      ! If we need to report on RAM use
      read(100,*,err=1100) cdump
      Input%RAMreport = cdump.eq.'Y'

      ! Raman scattering
      read(100,*,err=1100) cdump
      Input%Raman = cdump.eq.'Y'

      ! Non coherent lower term
      read(100,*,err=1100) cdump
      NCHLT = cdump.eq.'Y'

      ! Restrict redistribution in height
      read(100,*,err=1100) cdump
      Input%rest_z_red = cdump.eq.'Y'

      ! Minimum z for PRD
      if (Input%rest_z_red) &
        read(100,*,err=1100) Input%r1z_prd

      ! Restrict redistribution in tau_c
      read(100,*,err=1100) cdump
      Input%rest_tau_red = cdump.eq.'Y'

      ! Maximum tau for PRD
      if (Input%rest_tau_red) &
        read(100,*,err=1100) Input%r1tc_prd

      ! Coherent wings in the observers frame
      read(100,*,err=1100) Input%dcohw
      Input%cohw = Input%dcohw.gt.0d0

      ! Coherent wings in the observers frame
      read(100,*,err=1100) Input%dcohwi
      Input%cohwi = Input%dcohwi.gt.0d0

      ! Type of interpolation to convert to the observer's
      ! frame
      read(100,*,err=1100) Input%PRD_int_mode

      ! Type of redistribution (AA or not)
      read(100,*,err=1100) cdump
      Input%AV = cdump.eq.'A'
      AV = Input%AV

      ! Force angle-averaged intensity redistribution
      read(100,*,err=1100) cdump
      Input%AVI = cdump.eq.'Y'

      ! Either by force or because the problem is that way
      AVI = Input%AVI.or.AV

      ! Two-step for angle-dependent
      read(100,*,err=1100) cdump
      Input%two_step_AD = cdump.eq.'Y'

      ! Check necessity of two-step
      if (Input%force.eq.'I'.and.AVI) then
        Input%two_step_AD = .False.
      else if (AV) then
        Input%two_step_AD = .False.
      end if

      ! Store Wfunc for intensity in RAM
      read(100,*,err=1100) cdump
      IRAM = cdump.eq.'Y'

      ! Store Wfunc for polarization in RAM
      read(100,*,err=1100) cdump
      PRAM = cdump.eq.'Y'

      ! Directional nodes for angle average redistribution
      read(100,*,err=1100) Input%nThAA
      read(100,*,err=1100) Input%nThAAI

      ! Parameters for the selection of input frequencies
      read(100,*,err=1100) ddump
      Input%red_pars(1) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(2) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(3) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(4) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(5) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(6) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(7) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(8) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(9) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(10) = ddump
      read(100,*,err=1100) ddump
      Input%red_pars(11) = ddump

      ! Parameters for the selection of input frequencies (intensity)
      read(100,*,err=1100) ddump
      Input%redi_pars(1) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(2) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(3) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(4) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(5) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(6) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(7) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(8) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(9) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(10) = ddump
      read(100,*,err=1100) ddump
      Input%redi_pars(11) = ddump

      ! Type of Doppler width to build output frequency axis
      read(100,*,err=1100) Input%dws

      ! If the input is flagged as numeric, read the number
      if (Input%dws.eq.'NUM') read(100,*) Input%dw

      ! Force microturbulence
      read(100,*,err=1100) Input%fvmicro

      ! Minimum expected temperature
      read(100,*,err=1100) Input%minT

      ! Maximum expected temperature
      read(100,*,err=1100) Input%maxT

      ! Maximum expected velocity
      read(100,*,err=1100) Input%maxV

      ! Number of CPU per 1.5D column
      read(100,*,err=1100) Input%rt_group_n

      ! Force unmagnetized
      read(100,*,err=1100) cdump
      Input%unmagnetized = cdump.eq.'Y'

      ! Force static
      read(100,*,err=1100) cdump
      Input%static = cdump.eq.'Y'

      ! Force static intensity
      read(100,*,err=1100) cdump
      Input%static_int = cdump.eq.'Y'

      ! Skip star disk
      read(100,*,err=1100) cdump
      Input%skip_disk = cdump.eq.'Y'

      ! Initialize radiation field with b-b transitions
      read(100,*,err=1100) cdump
      Input%init_J_bb = cdump.eq.'Y'

      ! Index of first iteration
      read(100,*,err=1100) Input%iter_min
      read(100,*,err=1100) Input%iteri_min

      ! Maximum number of iterations
      read(100,*,err=1100) Input%iter_max
      read(100,*,err=1100) Input%iteri_max

      ! Maximum number of AD iterations
      read(100,*,err=1100) Input%iterad_max
      read(100,*,err=1100) Input%iteriad_max

      ! Order of emissivity
      read(100,*,err=1100) cdump
      PRD = cdump.eq.'Y'
      if (PRD) then
        Input%iter_ord = 2
      else
        Input%iter_ord = 1
      end if

      ! If solving in two steps
      read(100,*,err=1100) cdump
      Input%two_step_pol = cdump.eq.'Y'

      ! Maximum relative change to consider convergence
      read(100,*,err=1100) Input%mrcj
      read(100,*,err=1100) Input%mrc_i
      read(100,*,err=1100) Input%mrci_i
      read(100,*,err=1100) Input%mrc_p
      read(100,*,err=1100) Input%mrc_adi
      read(100,*,err=1100) Input%mrci_adi
      read(100,*,err=1100) Input%mrc_adp

      ! Iterations of radiation field
      read(100,*,err=1100) Input%iter_j

      ! PRD internal iterations (intensity)
      read(100,*,err=1100) Input%iteri_prd

      ! Sanity check
      if (Input%iteri_prd.lt.1) Input%iteri_prd = 1
      if (Input%iteri_prd.gt.99) Input%iteri_prd = 99

      ! PRD internal iterations MRC
      read(100,*,err=1100) Input%mrci_r

      ! PRD internal iterations
      read(100,*,err=1100) Input%iter_prd

      ! Sanity check
      if (Input%iter_prd.lt.1) Input%iter_prd = 1
      if (Input%iter_prd.gt.99) Input%iter_prd = 99

      ! PRD internal iterations MRC
      read(100,*,err=1100) Input%mrc_r
      read(100,*,err=1100) Input%mrc_p_r

      ! NG acceleration
      read(100,*,err=1100) cdump
      Input%NG = cdump.eq.'Y'

      ! NG acceleration order
      read(100,*,err=1100) Input%NG_ord

      ! Sanity check order
      if (Input%NG_ord.lt.1.or.Input%NG_ord.gt.5) Input%NG = .False.

      ! NG acceleration delay
      read(100,*,err=1100) Input%NG_delay

      ! NG acceleration for intensity
      read(100,*,err=1100) cdump
      Input%NGI = cdump.eq.'Y'

      ! NG acceleration order for intensity
      read(100,*,err=1100) Input%NGI_ord

      ! Sanity check order
      if (Input%NGI_ord.lt.1.or.Input%NGI_ord.gt.5) &
        Input%NGI = .False.

      ! NG acceleration delay for intensity
      read(100,*,err=1100) Input%NGI_delay

      ! PRD in intensity delay
      read(100,*,err=1100) Input%PRD_delay

      ! ALI photoionizations
      read(100,*,err=1100) cdump
      Input%ALI_photo = cdump.eq.'Y'

      ! ALI method delay
      read(100,*,err=1100) Input%ALI_delay

      ! ALI force
      read(100,*,err=1100) cdump
      Input%ALI_force = cdump.eq.'Y'

      ! ALI allow off
      read(100,*,err=1100) cdump
      Input%ALI_allow_off = cdump.eq.'Y'

      ! Append into the MRC file
      read(100,*,err=1100) cdump
      Input%appendMRC = cdump.eq.'Y'

      ! Append into the MRC file for intensity
      read(100,*,err=1100) cdump
      Input%appendMRCI = cdump.eq.'Y'

      ! Allowed iterations with non-physical Stokes, rho, and populations
      read(100,*,err=1100) Input%allownphys_stk
      read(100,*,err=1100) Input%allownphys_rho
      read(100,*,err=1100) Input%allownphys_pop

      ! Solution Box in 1.5D or inversion modes
      if(run_mode.eq.1.or.run_mode.eq.-1) then

        ! Allocate and read
        allocate(Input%sol_box(4))
        MRAMc = MRAMc + 1d-6*sizeof(Input%sol_box(4))
        read(100,*,err=1100) Input%sol_box

      ! Otherwise
      else

        ! Read dummy values
        read(100,*,err=1100) ios
        read(100,*,err=1100) ios
        read(100,*,err=1100) ios
        read(100,*,err=1100) ios

      end if ! Solution box

      ! Excluded pixels
      read(100,*,err=1100) Input%nexcl
      Input%lexcl = Input%nexcl.gt.0

      ! If there are excluded pixels
      if (Input%lexcl) then

        ! Allocate space
        allocate(Input%excl(2,Input%nexcl))
        MRAMC = MRAMc + 1d-6*sizeof(Input%excl)

        ! Read excluded pixels
        do i1=1,Input%nexcl
          read(100,*,err=1100) Input%excl(:,i1)
        end do

      end if ! There are excluded pixels

      ! Steps between storing solution
      read(100,*,err=1100) Input%store_step
      Input%store = Input%store_step.ge.1

      ! Correct dummy number
      if (.not.Input%store) Input%store_step = 1

      ! Steps between storing solution in intensity problem
      read(100,*,err=1100) Input%storei_step
      Input%storeI = Input%storei_step.ge.1

      ! Correct dummy number
      if (.not.Input%storeI) Input%storei_step = 1

      ! Output contribution function
      read(100,*,err=1100) cdump
      Input%out_contr = cdump.eq.'Y'

      ! Output tau=1 height
      read(100,*,err=1100) cdump
      Input%out_tau1 = cdump.eq.'Y'

      ! Keep background
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_back = .True.
      else
        Input%keep_back = .False.
      end if

      ! Keep damp
      read(100,*,err=1100) cdump
      Input%keep_damp = cdump.eq.'Y'

      ! Keep qel
      read(100,*,err=1100) cdump
      Input%keep_qel = cdump.eq.'Y'

      ! Keep equivalent-parametric damping parameters
      read(100,*,err=1100) cdump
      Input%keep_aparam = cdump.eq.'Y'

      ! Keep cols
      read(100,*,err=1100) cdump
      Input%keep_cols = cdump.eq.'Y'

      ! Keep atmo
      read(100,*,err=1100) cdump
      Input%keep_atmo = cdump.eq.'Y'

      ! Keep pop
      read(100,*,err=1100) cdump
      Input%keep_pop = cdump.eq.'Y'

      ! Keep dep
      read(100,*,err=1100) cdump
      Input%keep_dep = cdump.eq.'Y'

      ! Keep rhoKQ
      read(100,*,err=1100) cdump
      Input%keep_rhoKQ = cdump.eq.'Y'

      ! Keep JKQ
      read(100,*,err=1100) cdump
      Input%keep_JKQ = cdump.eq.'Y'

      ! Keep stokes in quadrature
      read(100,*,err=1100) cdump
      Input%keep_stokesQ = cdump.eq.'Y'

      ! Keep JKQnu
      read(100,*,err=1100) cdump
      Input%keep_jkqnu = cdump.eq.'Y'

      ! Keep MRC
      read(100,*,err=1100) cdump
      Input%keep_MRC = cdump.eq.'Y'

      ! Keep COL log
      read(100,*,err=1100) cdump
      Input%keep_coll = cdump.eq.'Y'

      ! Keep MPI log
      read(100,*,err=1100) cdump
      Input%keep_mpil = cdump.eq.'Y'

      ! Keep MPI detail log
      read(100,*,err=1100) cdump
      Input%keep_mpidl = cdump.eq.'Y'

      ! Limit Stokes output
      read(100,*,err=1100) Input%lim_stk%nran

      ! If ranges to limit Stokes output
      if (Input%lim_stk%nran.gt.0) then

        ! Allocate space
        allocate(Input%lim_stk%doub(2,Input%lim_stk%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_stk%doub)

        ! And read data
        do i1=1,Input%lim_stk%nran
          read(100,*,err=1100) Input%lim_stk%doub(1,i1)
          read(100,*,err=1100) Input%lim_stk%doub(2,i1)
        end do

      end if ! Ranges to limit Stokes output

      ! Limit contribution output
      read(100,*,err=1100) Input%lim_ctr%nran

      ! If ranges to limit contribution output
      if (Input%lim_ctr%nran.gt.0) then

        ! Allocate space
        allocate(Input%lim_ctr%doub(2,Input%lim_ctr%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_ctr%doub)

        ! And read data
        do i1=1,Input%lim_ctr%nran
          read(100,*,err=1100) Input%lim_ctr%doub(1,i1)
          read(100,*,err=1100) Input%lim_ctr%doub(2,i1)
        end do

      end if ! Ranges to limit contribution output

      ! Limit tau output
      read(100,*,err=1100) Input%lim_tau%nran

      ! If ranges to limit tau output
      if (Input%lim_tau%nran.gt.0) then

        ! Allocate space
        allocate(Input%lim_tau%doub(2,Input%lim_tau%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_tau%doub)

        ! And read data
        do i1=1,Input%lim_tau%nran
          read(100,*,err=1100) Input%lim_tau%doub(1,i1)
          read(100,*,err=1100) Input%lim_tau%doub(2,i1)
        end do

      end if ! Ranges to limit tau output

      ! Limit cols-TT output
      read(100,*,err=1100) Input%lim_cols_tt%nran

      ! If ranges to limit cols-TT output
      if (Input%lim_cols_tt%nran.gt.0) then

        ! Allocate space
        allocate(Input%lim_cols_tt%indx(3,Input%lim_cols_tt%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_cols_tt%indx)

        ! And read data
        do i1=1,Input%lim_cols_tt%nran
          read(100,*,err=1100) Input%lim_cols_tt%indx(1,i1)
          read(100,*,err=1100) Input%lim_cols_tt%indx(2,i1)
          read(100,*,err=1100) Input%lim_cols_tt%indx(3,i1)
        end do

      end if ! Ranges to limit cols-TT output

      ! Limit cols-LL output
      read(100,*,err=1100) Input%lim_cols_ll%nran

      ! If ranges to limit cols-LL output
      if (Input%lim_cols_ll%nran.gt.0) then

        ! Allocate space
        allocate(Input%lim_cols_ll%indx(3,Input%lim_cols_ll%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_cols_ll%indx)

        ! And read data
        do i1=1,Input%lim_cols_ll%nran
          read(100,*,err=1100) Input%lim_cols_ll%indx(1,i1)
          read(100,*,err=1100) Input%lim_cols_ll%indx(2,i1)
          read(100,*,err=1100) Input%lim_cols_ll%indx(3,i1)
        end do

      end if ! Ranges to limit cols-LL output

      ! Limit damping parameter output
      read(100,*,err=1100) Input%lim_damp%nran

      ! If ranges to limit damping parameter output
      if (Input%lim_damp%nran.gt.0) then

        ! Allocate
        allocate(Input%lim_damp%indx(2,Input%lim_damp%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_damp%indx)

        ! And read data
        do i1=1,Input%lim_damp%nran
          read(100,*,err=1100) Input%lim_damp%indx(1,i1)
          read(100,*,err=1100) Input%lim_damp%indx(2,i1)
        end do

      end if ! Ranges to limit damping parameter output

      ! Limit elastic collisional rates output
      read(100,*,err=1100) Input%lim_qel%nran

      ! If ranges to limit elastic collisional rates output
      if (Input%lim_qel%nran.gt.0) then

        ! Allocate
        allocate(Input%lim_qel%indx(2,Input%lim_qel%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_qel%indx)

        ! And read data
        do i1=1,Input%lim_qel%nran
          read(100,*,err=1100) Input%lim_qel%indx(1,i1)
          read(100,*,err=1100) Input%lim_qel%indx(2,i1)
        end do

      end if ! Ranges to limit elastic collisional rates output

      ! Limit background output
      read(100,*,err=1100) Input%lim_back%nran

      ! If ranges to limit background output
      if (Input%lim_back%nran.gt.0) then

        ! Allocate
        allocate(Input%lim_back%doub(2,Input%lim_back%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_back%doub)

        ! And read data
        do i1=1,Input%lim_back%nran
          read(100,*,err=1100) Input%lim_back%doub(1,i1)
          read(100,*,err=1100) Input%lim_back%doub(2,i1)
        end do

      end if ! Ranges to limit background output

      ! Limit population output
      read(100,*,err=1100) Input%lim_pop%nran

      ! If ranges to limit population output
      if (Input%lim_pop%nran.gt.0) then

        ! Allocate
        allocate(Input%lim_pop%indx(2,Input%lim_pop%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_pop%indx)

        ! And read data
        do i1=1,Input%lim_pop%nran
          read(100,*,err=1100) Input%lim_pop%indx(1,i1)
          read(100,*,err=1100) Input%lim_pop%indx(2,i1)
        end do

      end if ! Ranges to limit population output

      ! Recalculate electron density
      read(100,*,err=1100) Input%redo_ne

      ! Update atmospheric model after completing
      ! intensity calculations, and what to store
      read(100,*,err=1100) Input%update_atmos

      ! Keep input Hydrogen densities (atmos or atmo) after
      ! equation of state
      read(100,*,err=1100) cdump
      Input%protect_H = cdump.eq.'Y'

      ! Keep input Hydrogen minus densities (atmos) after
      ! equation of state
      read(100,*,err=1100) cdump
      Input%protect_Hm = cdump.eq.'Y'

      ! Do not change the atomic populations in the
      ! chemical equilibrium
      read(100,*,err=1100) cdump
      Input%chem_protect_all = cdump.eq.'Y'

      ! Write general performance
      read(100,*,err=1100) cdump
      Input%g_perf = cdump.eq.'Y'

      ! Write MPI performance
      read(100,*,err=1100) cdump
      Input%mpi_perf = cdump.eq.'Y'

      ! Verbosity
      read(100,*,err=1100) cdump
      verbosity = cdump.eq.'Y'

      ! Verbosity file
      verbosef = trim(Input%folder)//'/verbose'

      ! This MUST be false if not inverting
      Input%force_inv_freq = .False.

      !
      ! Inversion mode
      !
      if (Input%run_mode.eq.-1) then

        ! Allocation of inversion inputs
        allocate(Input%Node_type(Input%nvar))
        allocate(Input%Node(Input%nvar))
        allocate(Input%Num_nodes(Input%nvar))
        allocate(Input%extrapolation(Input%nvar))
        allocate(Input%Nodes_flags(Input%nvar))
        allocate(Input%Nodes_Regul(Input%nvar))
        allocate(Input%Indx_regul(Input%nvar))
        allocate(Input%Const_regul(Input%nvar))
        allocate(Input%Regul_weight(Input%nvar))
        allocate(Input%Scal(Input%nvar))
        allocate(Input%Perturb(Input%nvar))
        allocate(Input%min_rel_Pert(Input%nvar))
        MRAMc = MRAMc + 1d-6*(sizeof(Input%Node_type) + &
                              sizeof(Input%Node) + &
                              sizeof(Input%Num_nodes) + &
                              sizeof(Input%extrapolation) + &
                              sizeof(Input%Nodes_flags) + &
                              sizeof(Input%Nodes_Regul) + &
                              sizeof(Input%Indx_regul) + &
                              sizeof(Input%Const_regul) + &
                              sizeof(Input%Regul_weight) + &
                              sizeof(Input%Scal) + &
                              sizeof(Input%Perturb) + &
                              sizeof(Input%min_rel_Pert))

        ! Verbosity level
        read(100,*,err=1100) vlevel
        read(100,*,err=1100) slevel

        ! Data file
        read(100,'(A)',err=1100) Input%Filename_ob

        ! Type of inversion
        read(100,*,err=1100) Input%Type_Inversion

        ! Automatic Stokes weights
        read(100,'(A)',err=1100) cdump
        Input%auto_weight = cdump.eq.'Y'

        ! If weights automatic
        if (Input%auto_weight) then

          ! No read
          Input%Num_weight = 0
          Input%linv_weight = .False.

        ! No automatic weights
        else

          ! Read the number of the weights
          read(100,*,err=1100) Input%Num_Weight

          ! If larger than 0
          if (Input%Num_Weight.gt.0) then

            ! Allocate
            allocate(Input%Weight(0:5,Input%Num_Weight))
            MRAMc = MRAMc + 1d-6*sizeof(Input%Weight)

            ! For each set of weights
            do i1=1,Input%Num_Weight
              read(100,*,err=1100) Input%weight(:,i1)
            end do ! Sets of weights

          end if ! Numeric input

          ! Check if file
          read(100,*,err=1100) i1
          Input%linv_weight = i1.ne.0

          ! If file, read
          if (Input%linv_weight) &
            read(100,'(A)',err=1100) Input%inv_weight

          ! Weight factors
          read(100,*,err=1100) ios

          ! If weights factor to read
          if (ios.gt.0) then

            ! Allocate
            allocate(Input%Weight_Factor(4,ios))
            MRAMc = MRAMc + 1d-6*sizeof(Input%Weight_Factor)

            ! Read entries
            do i1=1,ios
              read(100,*,err=1100) Input%Weight_Factor(:,i1)
            end do

          end if ! Additional factors to read
        end if ! No automatic weights

        ! Sigma factors
        read(100,*,err=1100) ios

        ! If sigma factor to read
        if (ios.gt.0) then

          ! Allocate
          allocate(Input%Sigma_Factor(4,ios))
          MRAMc = MRAMc + 1d-6*sizeof(Input%Sigma_Factor)

          ! Read entries
          do i1=1,ios
            read(100,*,err=1100) Input%Sigma_Factor(:,i1)
          end do

        end if ! Additional factors to read

        ! Inversion initialization
        read(100,'(A)',err=1100) Input%Inv_init

        ! Inversion mask
        read(100,'(A)',err=1100) Input%Inv_mask

        ! Centered derivative
        read(100,'(A)',err=1100) cdump
        Input%centered = cdump.eq.'Y'

        ! Maximum iterations
        read(100,*,err=1100) Input%Num_Iter

        ! Read type/method of node
        do i1=1,Input%nvar
          read(100,*,err=1100) Input%Node_type(i1)
        end do

        ! Read interpolation
        read(100,*,err=1100) Input%Interpolation

        ! Read location/number for each variable
        do i1=1,Input%nvar

          ! Read flag
          read(100,*,err=1100) i2

          ! Read number of nodes
          read(100,*,err=1100) Input%Num_nodes(i1)

          ! If flag is 1, locations are expected now
          if (i2.eq.1) then

            ! Allocate space for nodes
            allocate(Input%Node(i1)%H(Input%Num_nodes(i1)))
            MRAMc = MRAMc + 1d-6*sizeof(Input%Node(i1)%H)

            ! And read positions
            read(100,*,err=1100) Input%Node(i1)%H

          end if ! Explicit locations

        end do ! For each variable

        ! Read extrapolation for each variable
        read(100,*,err=1100) Input%extrapolation

        ! Type of magnetic field vector
        read(100,*,err=1100) Input%btype

        ! Type of velocity vector
        read(100,*,err=1100) Input%vtype

        ! Flag each variable as changing or not
        do i1=1,Input%nvar

          ! Read string and changes if not fixed and non-zero nodes
          read(100,'(A)',err=1100) cdump
          Input%Nodes_Flags(i1) = Input%Num_nodes(i1).gt.0.and. &
                                  .not.(cdump.eq.'Y')

        end do ! For each variable

        ! Flag to correct positions
        read(100,'(A)',err=1100) cdump
        Input%Pos_Correction = cdump.eq.'Y'

        ! Regularization type
        do i1=1,Input%nvar

          ! Read index
          read(100,*,err=1100) Input%Indx_regul(i1)

          ! If not none
          if (Input%Indx_regul(i1).gt.0) then

            ! Read
            read(100,*,err=1100) Input%Regul_weight(i1)

            ! If constant, real value
            if (Input%Indx_regul(i1).eq.2) &
              read(100,*,err=1100) Input%Const_regul(i1)

          ! No regularization
          else

            ! No weight
            Input%Regul_weight(i1) = 0d0

          end if ! Non-none regularization

        end do ! Variables

        ! Regularization limit
        read(100,*,err=1100) Input%Regul_Limit

        ! Regularization factor
        read(100,*,err=1100) Input%Regul_factor

        ! Chi2 Threshold
        read(100,*,err=1100) Input%Threshold_chisq

        ! Inversion MRC
        read(100,*,err=1100) Input%Chisq_fraction

        ! Read SVD type
        read(100,*,err=1100) Input%SVD_type

        ! Threshold for SVD
        read(100,*,err=1100) Input%Threshold_svd

        ! Type of calculation for pressure
        read(100,*,err=1100) i1
        Input%hydroeq = i1.eq.1

        ! Boundary value for the pressure
        read(100,*,err=1100) Input%Pg_bound

        ! Diffuse light factor
        read(100,*,err=1100) Input%f_diff

        ! Read path to initial atomic model
        read(100,'(A)',err=1100) Input%atmo

        ! Hard-coded initializations
        if (trim(Input%atmo).eq.'NONE') then

          ! FALC hardcoded
          Input%atmo = 'NONE'
          Input%Init_Thermal = 0

        else if (trim(Input%atmo).eq.'$$C$$') then

          ! FALC hardcoded
          Input%atmo = 'NONE'
          Input%Init_Thermal = 0

        else if (trim(Input%atmo).eq.'$$P$$') then

          ! FALP hardcoded
          Input%atmo = 'NONE'
          Input%Init_Thermal = 1

        end if ! Hard-coded initializations

        ! Number of atmospheric nodes in synthesis
        read(100,*,err=1100) Input%Atmo_Input

        ! Maximum step for SVD
        read(100,*,err=1100) Input%Max_Step

        ! Type of error
        read(100,*,err=1100) Input%Err_Type

        ! PSF FWHM
        read(100,*,err=1100) ios

        ! If reading PSF
        if (ios.gt.0) then

          ! Allocate structures and names
          allocate(Input%lim_fwhm(ios))
          allocate(Input%fwhm_fil(ios))
          MRAMc = MRAMc + 1d-6*(sizeof(Input%lim_fwhm) + &
                                sizeof(Input%fwhm_fil))

          ! For each range
          do i1=1,ios

            ! Save size data (redundant)
            Input%lim_fwhm(i1)%nn = ios

            ! Read type
            read(100,*,err=1100) i2

            ! Check if Gaussian
            Input%lim_fwhm(i1)%gaussian = i2.eq.0

            ! If gaussian
            if (Input%lim_fwhm(i1)%gaussian) then

              ! Allocate doubles
              allocate(Input%lim_fwhm(i1)%doub(3))
              MRAMc = MRAMc + 1d-6*sizeof(Input%lim_fwhm(i1)%doub)

              ! Read entry
              read(100,*,err=1100) Input%lim_fwhm(i1)%doub

            ! If NOT gaussian
            else

              ! Allocate doubles
              allocate(Input%lim_fwhm(i1)%doub(2))
              MRAMc = MRAMc + 1d-6*sizeof(Input%lim_fwhm(i1)%doub)

              ! Read entry
              read(100,*,err=1100) Input%lim_fwhm(i1)%doub

              ! Get file
              read(100,'(A)',err=1100) Input%fwhm_fil(i1)

            end if ! Gaussian

          end do ! For each range

        end if ! There is PSF to read

        ! Bounds for parameters
        do i1=1,Input%nvar
          read(100,*,err=1100) Input%Node(i1)%Bounds
        end do

        ! Check if special limits
        read(100,*,err=1100) ia

        ! If there were special limits
        if (ia.gt.0) then

          ! For each variable
          do i1=1,Input%nvar

            ! Read number of inputs
            read(100,*,err=1100) ios

            ! If there are inputs
            if (ios.gt.0) then

              ! Allocate special limits
              allocate(Input%Node(i1)%ebound(4,ios))
              MRAMc = MRAMc + 1d-6*sizeof(Input%Node(i1)%ebound)

              ! Save number
              Input%Node(i1)%nebound = ios

              ! For each entry
              do i2=1,ios

                ! Read limits and boundaries
                read(100,*,err=1100) Input%Node(i1)%ebound(:,i2)

              end do ! Number of special boundaries

            end if ! Input for this variable

          end do ! Variables

        end if ! There were special limits

        ! Scale for parameters
        do i1=1,Input%nvar
          read(100,*,err=1100) Input%Scal(i1)
        end do

        ! Perturbation for parameters
        do i1=1,Input%nvar
          read(100,*,err=1100) Input%Perturb(i1)
        end do

        ! Minimum relative perturbation
        do i1=1,Input%nvar
          read(100,*,err=1100) Input%min_rel_Pert(i1)
        end do

        ! Value when too small Bpos initialization
        read(100,*,err=1100) Input%ini_Bpos

        ! Value when too small Bazi initialization
        read(100,*,err=1100) Input%ini_Bazi

        ! Value when too small Vpos initialization
        read(100,*,err=1100) Input%ini_vpos

        ! Value when too small Vazi initialization
        read(100,*,err=1100) Input%ini_vazi

        ! Guess the polarity of the field
        read(100,'(A)',err=1100) cdump
        Input%guess_polarity_l = cdump.eq.'Y'

        ! If guessing polarity
        if (Input%guess_polarity_l) then

          ! Wavelength limits
          read(100,*,err=1100) Input%gp_l
          read(100,*,err=1100) Input%gp_r
          read(100,*,err=1100) Input%gp_g
          read(100,*,err=1100) Input%gp_w

        end if

        ! Invert fractional Stokes
        read(100,'(A)',err=1100) cdump
        Input%Fractional = cdump.eq.'Y'

        ! Range of tau for atmosphere
        read(100,*,err=1100) Input%Tau_Range

        ! Modifications on the stratification of the
        ! model atmosphere
        read(100,*,err=1100) i1

        ! If input data for modification of stratification
        if (i1.gt.0) then

          ! Allocate
          allocate(Input%Atmo_strat(3,i1))
          MRAMc = MRAMc + 1d-6*sizeof(Input%Atmo_strat)

          ! For each input
          do i2=1,i1

            ! Read parameters
            read(100,*,err=1100) Input%Atmo_strat(:,i2)

          end do ! For each input

        end if ! Input data for stratification

        ! Broyden in LM
        read(100,'(A)',err=1100) cdump
        Input%Broyden = cdump.eq.'Y'

        ! Index of LM method
        read(100,*,err=1100) Input%LM_Method

        ! Index of LM mode for backtracking
        read(100,*,err=1100) Input%LM_Back_Mode

        ! Value of lambda to check big
        read(100,*,err=1100) Input%LM_lam_big_test

        ! Value of lambda to check small
        read(100,*,err=1100) Input%LM_lam_small_test

        ! Value of lambda to prove big
        read(100,*,err=1100) Input%LM_lam_big_prove

        ! Value of lambda to prove small
        read(100,*,err=1100) Input%LM_lam_small_prove

        ! Order of tracking
        read(100,*,err=1100) Input%Lam_track
        ! Zero or negative not tracking
        if (Input%Lam_track.lt.1) then
            Input%Lam_track = 0
            Input%l_Lam_track = .False.
        else
            Input%l_Lam_track = .True.
            ! Limit is second order
            if (Input%Lam_track.gt.3) Input%Lam_track = 3
        end if

        ! Range of lambda for LM
        read(100,*,err=1100) Input%Lam_Range

        ! Accepted lambda factor
        read(100,*,err=1100) Input%factoraccept

        ! Rejected lambda factor
        read(100,*,err=1100) Input%factorreject

        ! Project the B vector (negative strength values become
        ! negative polarity)
        read(100,'(A)',err=1100) cdump
        Input%Projection = cdump.eq.'Y'

        ! Initialize from previous solution when computing
        ! response functions
        read(100,'(A)',err=1100) cdump
        Input%Popuinit = cdump.eq.'Y'

        ! Initialize from previous solution when computing
        ! trials if thermo is fixed
        read(100,'(A)',err=1100) cdump
        Input%trialinit = cdump.eq.'Y'

        ! Initialize from previous solution when computing
        ! trials if T and Pg are fixed
        read(100,'(A)',err=1100) cdump
        Input%trialtpinit = cdump.eq.'Y'

        ! Neglect sigma
        read(100,'(A)',err=1100) cdump
        Input%Sigma_neglect = cdump.eq.'Y'

        ! Keep RF
        read(100,'(A)',err=1100) cdump
        Input%Keep_RF = cdump.eq.'Y'

        ! Steps between storing solution
        read(100,*,err=1100) Input%storeinv_step
        Input%storeinv = Input%storeinv_step.ge.1

        ! Correct dummy number
        if (.not.Input%storeinv) Input%storeinv_step = 1

        ! Force inversion frequencies
        read(100,'(A)',err=1100) cdump
        Input%force_inv_freq = cdump.eq.'Y'

      end if ! Inversion mode

      ! Close translated file
      close(100)

      ! If doing PRD
      if (PRD) then

        ! Correct NGI_delay and iterations if doing PRD
        Input%NGI_delay = Input%NGI_delay + Input%PRD_delay - 1
        Input%iteri_max = Input%iteri_max + Input%PRD_delay - 1

      ! CRD
      else

        ! No delay in PRD
        Input%PRD_delay = 0

      end if ! PRD/CRD

      ! Get CPU ID in character
      write(CPUC,'(I0.5)') pid

      ! Store name for error files
      errorf = trim(Input%folder)//'/ERROR'//CPUC

      ! Control that everything went fine
      call control

      ! Master
      if(pid.eq.0) then

        ! Delete the temporal input file and communicate success
        call system('rm tmp_input'//'_'//Input%ID)
        umsg = ' - Input read'
        call verbose

      end if

      return

1000  umsg = 'Error opening input file'
      call gaborted
1100  umsg = 'Error reading input file'
1200  close(100)
      call gaborted

      end subroutine rInput

!#####################################################################
!#####################################################################
!#####################################################################

      end module rinput_mod
