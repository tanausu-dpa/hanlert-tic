      !> Reading settings
      module rinput_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Hao Li (IAC)
!     Roberto Casini (HAO)
!  Start:
!     04/17/2017
!  Last version:
!     09/29/2023 V3.0.16
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:   V3.0.16 - Read Input%Kcut_input, Input%keep_coll,
!                             Input%keep_mpil, and variable
!                             Input%keep_mpidl (TdPA)
!                           - Initialize Kradl (TdPA)
!                           - Removed verbosity for when trying to
!                             store LTE profiles, but not active
!                             atoms (TdPA)
!
!     09/08/2023:   V3.0.15 - Change destiny of verbosity level to
!                             vlevel common variable (TdPA)
!                           - Read slevel (TdPA)
!
!     08/24/2023:   V3.0.14 - Added the possibility of input files
!                             for the inversion FWHM (TdPA)
!                           - Read (or initialize if not inversion)
!                             Input%force_inv_freq (TdPA)
!
!     08/11/2023:   V3.0.13 - Update the read of Input%weight (TdPA)
!                           - Read Input%inv_weight and
!                             Input%linv_weight (TdPA)
!
!     08/07/2023:   V3.0.12 - Read LTE lines data into Input%LTEline()
!                             structure array (TdPA)
!                           - LVIRAM and LVPRAM (TdPA)
!
!     07/03/2023:   V3.0.11 - Read Input%skip_wave, Input%fixplt,
!                             Input%vtype, Input%f_diff,
!                             Input%min_rel_pert, Input%ini_vpos,
!                             and Input%ini_vazi (TdPA)
!                           - Removed Input%initpixel and
!                             Input%Restartpixel (TdPA)
!                           - Put back projection (TdPA)
!                           - Removed region (TdPA)
!                           - Renamed Input%Restore_File to
!                             Input%Inv_init (TdPA)
!                           - Renamed Input%Pg_flag to
!                             Input%hydroeq (TdPA)
!                           - Renamed Input%FWHM to Input%lim_fwhm
!                             and updated to the new format (TdPA)
!                           - Inversion reads Input%sol_box (TdPA)
!                           - Check that nodes are larger than 0
!                             before confirming that a variable is
!                             to be inverted (TdPA)
!                           - Added sanity check in case you skip
!                             wavelengths for every atom without
!                             specifying a wavelength file (TdPA)
!
!     06/12/2023:   V3.0.10 - Removed Projection, NUM_FILE and
!                             INIT_INV_PIXEL (HL)
!                           - Added REGION (HL)
!                           - Update RESTART_INV_PIXEL[2] for 2
!                             dimensions (HL)
!                           - Update read FWHM for multi wavelength
!                             ranges (HL)
!
!     04/25/2023:    V3.0.9 - Bugfix: Added the NONE case for the
!                             initial atmospheric model in the
!                             inversion (TdPA)
!                           - Read Input%Keep_RF (TdPA)
!
!     04/11/2023:    V3.0.8 - Read weights for different ranges of
!                             wavelengths (HL)
!                           - Removed unused inversion inputs that
!                             were eliminated (see rinput.py; HL)
!
!     03/15/2023:    V3.0.7 - Removed unused inversion inputs that
!                             were eliminated (see rinput.py; TdPA)
!                           - Read new input quantities to configure
!                             the inversion (see rinput.py; TdPA)
!
!     03/08/2023:    V3.0.6 - Read new input quantities to configure
!                             the inversion (see rinput.py; TdPA)
!
!     02/14/2023:    V3.0.5 - Read new input quantities to split
!                             the geometry and PRD integral between
!                             intensity and polarization (TdPA)
!                           - Read Input%AVI (TdPA)
!                           - Read Input%static_int (TdPA)
!                           - Geom is not an argument anymore (TdPA)
!
!     11/10/2022:    V3.0.4 - Read Input%zero_ion and force_asym
!                             variables (TdPA)
!                           - Changed where Input%ionf is read (TdPA)
!                           - Bugfix: Wrong format when reading the
!                             Input%asym_fil input (TdPA)
!                           - Now stm is set directly to the global
!                             instance (TdPA)
!
!     10/25/2022:    V3.0.3 - Read Input%ionf, Input%spect_input,
!                             Input%chianti_path, Input%T_rad,
!                             Input%R_star, Input%r0tc, Input%r1tc,
!                             Input%r0z, Input%r1z, Input%skip_disk,
!                             and Input%sol_box (TdPA)
!
!     07/18/2022:    V3.0.2 - Read Input%MPIdetail, Input%operform,
!                             and Input%IWskip (TdPA)
!
!     07/13/2022:    V3.0.1 - Read Input%pf, Input%abund,
!                             Input%bark_sp, Input%bark_pd, and
!                             Input%bark_df (TdPA)
!
!     07/08/2022:    V3.0.0 - Fixed a comment (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case now rInput
!                             reads Input%run_mode, Input%atm_scale,
!                             Input%omega_ref, Input%atmo_char,
!                             Input%cache, Input%keep_sol, Input%minT,
!                             Input%maxV, Input%rt_group_n,
!                             Input%unmagnetized, Input%static,
!                             Input%keep_pop, Input%keep_dep,
!                             Input%keep_rhoKQ, Input%keep_JKQ,
!                             Input%keep_stokesQ, Input%keep_MRC,
!                             Input%lim_stk, Input%lim_ctr,
!                             Input%lim_tau, Input%lim_cols_tt,
!                             Input%lim_cols_ll, Input%lim_damp,
!                             Input%lim_back, Input%lim_pop (TdPA)
!                           - Changed calls to aborted into calls
!                             to gaborted (TdPA)
!
!     06/21/2022:    V2.0.2 - Read Input%dcohw (TdPA)
!
!     04/07/2022:    V2.0.1 - Read Input%keep_jkqnu (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed read of Input%ND (TdPA)
!
!     02/17/2021:   V1.2.24 - Read fcol_transfer (TdPA)
!                           - Put back error checker when reading
!                             asymmetry variables (TdPA)
!
!     02/12/2021:   V1.2.23 - Read Input%g_perf (TdPA)
!                           - Read Input%mpi_perf (TdPA)
!
!     02/04/2021:   V1.2.22 - Read Input%MPI_input (TdPA)
!                           - Read Input%MIT_node (TdPA)
!
!     01/13/2021:   V1.2.21 - Read asymmetry inputs (TdPA)
!                           - Read KSTK value (TdPA)
!
!     11/12/2020:   V1.2.20 - Read Input%PRD_delay and correct the
!                             Input%NGI_delay and Input%iteri_max
!                             if doing PRD. Read the input variable
!                             Input%chem_protect_all (TdPA)
!
!     09/20/2020:   V1.2.19 - Define the name of the error file (TdPA)
!
!     09/11/2020:   V1.2.18 - Read the variables Input%RAMreport,
!                             TIRAM, and TPRAM (TdPA)
!
!     07/22/2020:   V1.2.17 - Read the variables Input%NGI_delay,
!                             Input%NGI, and Input%NGI_ord (TdPA)
!
!     06/26/2020:   V1.2.16 - Read the variable KcutAB (TdPA)
!
!     06/01/2020:   V1.2.15 - Read the variable NCHLT (TdPA)
!
!     03/05/2020:   V1.2.14 - Read the variable Input%protect_H (TdPA)
!
!     02/10/2020:   V1.2.13 - Read the variables Input%redo_ne and
!                             Input%update_atmos (TdPA)
!
!     12/17/2019:   V1.2.12 - Read Input%zeeman_mode (TdPA)
!
!     12/10/2019:   V1.2.11 - Now it can receive no LOS angles (TdPA)
!                           - Reading VOITY and Input%memo (TdPA)
!
!     11/19/2019:   V1.2.10 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - Reading WIFIL and WPFIL (TdPA)
!
!     11/13/2019:    V1.2.9 - Reading PIRAM, VIFIL, and VPFIL (TdPA)
!
!     09/26/2019:    V1.2.8 - Added keep_atmo read (TdPA)
!
!     09/13/2019:    V1.2.7 - Added waves read (TdPA)
!
!     08/09/2019:    V1.2.6 - Added ignore_bb read (TdPA)
!
!     06/11/2019:    V1.2.5 - The magnetic field in numeric form was
!                             expecting to write into characters, not
!                             numbers (TdPA)
!
!     05/08/2019:    V1.2.4 - Reads allownphys_stk and
!                             allownphys_rho (TdPA)
!
!     04/15/2019:    V1.2.3 - Reads alt_bcast (TdPA)
!
!     03/22/2019:    V1.2.2 - Reads NK and kurucz (TdPA)
!
!     03/13/2019:    V1.2.1 - Reads ALI_delay and KEEP_APARAM (TdPA)
!
!     02/20/2019:    V1.2.0 - New verbosity (TdPA)
!                           - Checks for success of python routine
!                             and unit is now 100 (TdPA)
!                           - Now admits numerical inputs for
!                             magnetic fields (TdPA)
!
!     11/06/2018     V1.1.3 - Reading keep_back, keep_damp, and
!                             keep_cols (TdPA)
!
!     09/20/2018     V1.1.2 - Reading NG, NG_ord, and NG_delay (TdPA)
!
!     09/06/2018     V1.1.1 - Reading SOLUTION_KEEPI (TdPA)
!
!     08/09/2018     V1.1.0 - Allows to use a custom filename (TdPA)
!
!     08/04/2018    V1.0.12 - Reading VOI_IRAM, VOI_PRAM, and RAM_LIM
!                             (TdPA)
!
!     08/03/2018    V1.0.11 - Reading Krad (TdPA)
!
!     12/18/2017:   V1.0.10 - Bugfix: Was defining AV after reading
!                             Raman, no impact in the run becase was
!                             redefined later (TdPA)
!
!     12/05/2017:    V1.0.9 - Reading Raman (TdPA)
!
!     11/27/2017:    V1.0.8 - Reading Pcorr (TdPA)
!
!     09/22/2017:    V1.0.7 - Reading Kcut (TdPA)
!
!     09/14/2017:    V1.0.6 - Added a path and ID to the file (TdPA)
!
!     08/24/2017:    V1.0.5 - Reading five more red_params (TdPA)
!
!     06/16/2017:    V1.0.4 - Added read of PRAM and changed the read
!                             of RAM to the read of IRAM (TdPA)
!
!     06/12/2017:    V1.0.3 - Added read of RAM (TdPA)
!
!     06/09/2017:    V1.0.2 - Added input iteri_prd, iter_j, and
!                             mrci_r (TdPA)
!
!     06/08/2017:    V1.0.1 - Added input Geom%nThAA (TdPA)
!
!     04/17/2017:    V1.0.0 - First version (TdPA)
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
!    This program reads an input file (INPUT file where it is
!    executed, this is hardcoded) with information about the
!    execution that is not contained in the other input files
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

      !> Reads a file with the settings.\n
      !!     Input(Input_class): Structure with settings data
      subroutine rInput(Input)

      ! I/O

      type(Input_class), intent(inout):: Input

      ! Local

      character(LEN=1):: cdump
      character(LEN=2):: cdump2
      character(LEN=5):: CPUC

      integer:: ios, ia, i1, i2, i3, i4, i5

      double precision:: ddump, ddump2


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

      ! Number of atoms
      read(100,*,err=1100) Input%nA
      nA = Input%nA

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

        ! Atom file model
        read(100,'(A)',err=1100) Input%atom(ia)
        ! Atom population file
        read(100,'(A)',err=1100) Input%popu(ia)

        ! Atom fix population
        read(100,'(A)',err=1100) cdump

        ! Translate input
        if (cdump.eq.'F') then
          Input%fixp(ia) = .True.
        else
          Input%fixp(ia) = .False.
        end if

        ! Atom zero ion
        read(100,'(A)',err=1100) cdump

        ! Translate input
        if (cdump.eq.'F') then
          Input%zero_ion(ia) = .True.
        else
          Input%zero_ion(ia) = .False.
        end if

        ! Ionization files for CLE
        if (Input%run_mode.eq.2) then

          read(100,*,err=1100) Input%ionf(ia)%typ
          if (Input%ionf(ia)%typ.eq.0) then
            read(100,'(A)',err=1100) Input%ionf(ia)%str
          else if (Input%ionf(ia)%typ.eq.1) then
            read(100,*,err=1100) Input%ionf(ia)%val
          end if

        ! No CLE
        else

          ! Read dummy -1
          read(100,*,err=1100) ios

        end if

        ! Atom skip wavelengths
        read(100,'(A)',err=1100) cdump

        ! Translate input
        if (cdump.eq.'N') then
          Input%skip_wave(ia) = .False.
        else
          Input%skip_wave(ia) = .True.
        end if

        ! Atom fix terms populations
        read(100,'(A)',err=1100) cdump

        ! Translate input
        if (cdump.eq.'F') then
          Input%fixplt(ia) = .True.
        else
          Input%fixplt(ia) = .False.
        end if

        ! Check cascades of conditions for fixing populations
        if (Input%fixp(ia)) Input%fixplt(ia) = .False.

      end do

      ! Number of background atoms
      read(100,*,err=1100) Input%nAb
      nAb = Input%nAb

      ! If greater than 0
      if (Input%nAb.gt.0) then

        ! Allocate
        ! Filename of atomic models for background
        allocate(Input%atomback(Input%nAb))
        ! Filename of population file for background
        allocate(Input%popuback(Input%nAb))

        ! Read atomic file and population file names
        do ia=1,Input%nAb
          read(100,'(A)',err=1100) Input%atomback(ia)
          read(100,'(A)',err=1100) Input%popuback(ia)
        end do

      end if ! There are background atoms

      ! Number of molecules
      read(100,*,err=1100) Input%nM
      nM = Input%nM

      ! If greater than 0
      if (Input%nM.gt.0) then

        ! Allocate molecule file name
        allocate(Input%mol(Input%nM))

        ! Read molecule file names
        do ia=1,Input%nM
          read(100,'(A)',err=1100) Input%mol(ia)
        end do

      end if ! There are molecules

      ! Magnetic field
      read(100,'(A)',err=1100) cdump
      if (cdump.eq.'F') then
        Input%bfieldn = .False.
        read(100,'(A)',err=1100) Input%bfield
      else
        Input%bfieldn = .True.
        read(100,*,err=1100) Input%bfieldv(1)
        read(100,*,err=1100) Input%bfieldv(2)
        read(100,*,err=1100) Input%bfieldv(3)
      end if

      ! File with fudge factors
      read(100,'(A)',err=1100) Input%fudge

      ! File with input spectra
      read(100,'(A)',err=1100) Input%spect_input

      ! Path to CHIANTI database
      read(100,'(A)',err=1100) Input%chianti_path

      ! T_rad
      read(100,*,err=1100) Input%T_rad

      ! R_star
      read(100,*,err=1100) Input%R_star

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
      if(cdump.eq.'Y')then
        Input%addbb = .False.
      else
        Input%addbb = .True.
      endif

      ! Read number of Kurucz files
      read(100,*,err=1100) Input%NK

      ! If Kurucz files
      if (Input%NK.ge.1) then

        ! Allocate Kurucz file names
        allocate(Input%kurucz(Input%NK))

        ! Read Kurucz file names
        do ia=1,Input%NK
          read(100,'(A)',err=1100) Input%kurucz(ia)
        end do

      end if

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
          nullify(Input%LTEline(ia)%prof)

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
          read(100,'(A)',err=1100) Input%waves(ia)
        end do

      end if

      ! Sanity check
      if (all(Input%skip_wave).and.Input%NW.eq.0) then
        umsg = 'Cannot neglect automatically generated '// &
               'wavelengths without including wavelength files'
        call gaborted
      end if

      ! Asymmetry inputs
      read(100,*,err=1100) Input%nasym

      ! If there are inputs
      if (Input%nasym.gt.0) then

        ! Asymmetry inputs as numbers
        read(100,*,err=1100) Input%nasym_num

        ! If numerical values
        if (Input%nasym_num.gt.0) &
          allocate(Input%asym_num(2,Input%nasym_num))

        ! Asymmetry inputs as files
        read(100,*,err=1100) Input%nasym_fil

        ! If files
        if (Input%nasym_fil.gt.0) &
          allocate(Input%asym_fil(Input%nasym_fil))

        ! Initialize indexes for numerical and files
        i2 = 0
        i3 = 0

        ! For each entry
        do i1=1,Input%nasym

          ! Read character
          read(100,'(A)',err=1100) cdump

          ! If value
          if (cdump.eq.'V') then
            i2 = i2 + 1
            read(100,*,err=1100) i4, i5, ddump, ddump2
            Input%asym_num(1,i2) = dcmplx(i4,i5)
            Input%asym_num(2,i2) = dcmplx(ddump,ddump2)
          end if

          ! If file
          if (cdump.eq.'F') then
            i3 = i3 + 1
            read(100,'(A)',err=1100) Input%asym_fil(i3)
          end if

        end do ! Entries

      end if ! Asymmetry inputs

      ! Force asymm
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        force_asym = .True.
      else
        force_asym = .False.
      endif

      ! Stimulated emission
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        stm = .True.
      else
        stm = .False.
      endif

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
      read(100,*,err=1100) Input%nThLOS
      if (Input%nThLOS.gt.0) then
        allocate(Input%L_mu(Input%nThLOS))
        read(100,*,err=1100) Input%L_mu(1:Input%nThLOS)
      end if
      read(100,*,err=1100) Input%nPhLOS
      if (Input%nPhLOS.gt.0) then
        allocate(Input%L_phi(Input%nPhLOS))
        read(100,*,err=1100) Input%L_phi(1:Input%nPhLOS)
      end if

      ! Mode of solution
      read(100,*,err=1100) Input%mode

      ! Force type of calculation
      read(100,*,err=1100) Input%force

      ! Restrict problem in tau_c
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%rest_tau = .True.
        read(100,*,err=1100) Input%r0tc
        read(100,*,err=1100) Input%r1tc
      else
        Input%rest_tau = .False.
      end if

      ! Restrict problem in height
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%rest_z = .True.
        read(100,*,err=1100) Input%r0z
        read(100,*,err=1100) Input%r1z
      else
        Input%rest_z = .False.
      end if

      ! Mode of Zeeman effect
      read(100,*,err=1100) Input%zeeman_mode

      ! Correction of rho
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%Pcorr = .True.
      else
        Input%Pcorr = .False.
      end if

      ! Forbidden collisions can transfer alignment
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        fcol_transfer = .True.
      else
        fcol_transfer = .False.
      end if

      ! MIT transitions
      read(100,*,err=1100) Input%MIT_input

      ! MIT factor
      read(100,*,err=1100) Input%MIT_node
      Input%MIT_node = 1d0/Input%MIT_node

      ! Skip saving solution
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_sol = .False.
      else
        Input%keep_sol = .True.
      end if

      ! Solution file
      read(100,'(A)',err=1100) Input%solution

      ! Dummy line
      read(100,*,err=1100) cdump

      ! Keep intensity solution
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keepIsol = .True.
      else
        Input%keepIsol = .False.
      end if

      ! Keep Stokes in solution
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        KSTK = .True.
      else
        KSTK = .False.
      end if

      ! Cut rhoKQ multipole orders
      read(100,*,err=1100) Kcut

      if (Kcut.lt.0) Kcut = 1000

      ! Cut rhoKQ multipole orders in absorb1
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        KcutAB = .True.
      else
        KcutAB = .False.
      end if

      ! K cut term-wise
      read(100,*,err=1100) ios

      ! If there are entries
      if (ios.gt.0) then

        ! Allocate
        allocate(Input%Kcut_input(4,ios))

        ! And read
        do ia=1,ios
          read(100,*,err=1100) Input%Kcut_input(:,ia)
        end do

      end if

      ! Cut JKQ multipole orders
      read(100,*,err=1100) Krad

      if (Krad.lt.0.or.Krad.gt.2) Krad = 2
      Kradl = Krad

      ! Store J symbols
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%memo = .True.
      else
        Input%memo = .False.
      end if

      ! Store photoionization quantities
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        PIRAM = .True.
      else
        PIRAM = .False.
      end if

      ! Type of intensity Voigt profile
      read(100,*,err=1100) VOITY

      ! Store Voigt profiles for intensity in RAM
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        VIRAM = .True.
      else
        VIRAM = .False.
      end if

      ! Store Voigt profiles for intensity in file
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        VIFIL = .True.
      else
        VIFIL = .False.
      end if

      ! Change RAM if yes file
      if(VIFIL.and.VIRAM)then
        VIRAM = .False.
        if (pid.eq.0) then
          umsg = ' - You have chosen to store the intensity '// &
                 'Voigt profiles in a file, so storing them '// &
                 'in RAM has been disabled'
          call verbose
        end if
      end if

      ! Store Voigt profiles for intensity in RAM for LTE
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        LVIRAM = .True.
      else
        LVIRAM = .False.
      end if

      ! Change RAM if no profiles
      if(.not.(VIRAM.or.VIFIL).and.LVIRAM) LVIRAM = .False.

      ! Store Voigt profiles for polarization in RAM
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        VPRAM = .True.
      else
        VPRAM = .False.
      end if

      ! Store Voigt profiles for polarization in file
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        VPFIL = .True.
      else
        VPFIL = .False.
      end if

      ! Change RAM if yes file
      if(VPFIL.and.VPRAM)then
        VPRAM = .False.
        if (pid.eq.0) then
          umsg = ' - You have chosen to store the Voigt '// &
                 'profiles in a file, so storing them in RAM has '// &
                 'been disabled'
          call verbose
        end if
      end if

      ! Store Voigt profiles for intensity in RAM for LTE
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        LVPRAM = .True.
      else
        LVPRAM = .False.
      end if

      ! Change RAM if no profiles
      if (.not.(VPRAM.or.VPFIL).and.LVPRAM) LVPRAM = .False.

      ! Limit in RAM for profiles
      read(100,*,err=1100) RLIM

      if (RLIM.lt.0) RLIM = 1000000000

      ! If we need to report on RAM use
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%RAMreport = .True.
      else
        Input%RAMreport = .False.
      endif

      ! Raman scattering
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%Raman = .True.
      else
        Input%Raman = .False.
      endif

      ! Non coherent lower term
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        NCHLT = .True.
      else
        NCHLT = .False.
      endif

      ! Coherent wings in the observers frame
      read(100,*,err=1100) Input%dcohw
      Input%cohw = Input%dcohw.gt.0d0

      ! Type of redistribution (AA or not)
      read(100,*,err=1100) cdump
      if(cdump.eq.'A')then
        Input%AV = .True.
      else
        Input%AV = .False.
      endif
      AV = Input%AV

      ! Force angle-averaged intensity redistribution
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%AVI = .True.
      else
        Input%AVI = .False.
      endif
      ! Either by force or because the problem is
      AVI = Input%AVI.or.AV

      ! Store Wfunc for intensity in RAM
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        IRAM = .True.
      else
        IRAM = .False.
      end if

      ! Store Wfunc for intensity in file
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        WIFIL = .True.
      else
        WIFIL = .False.
      end if

      ! Change RAM if yes file
      if(WIFIL.and.IRAM)then
        IRAM = .False.
        if (pid.eq.0) then
          umsg = ' - You have chosen to store the intensity '// &
                 'redistribution in a file, so storing them '// &
                 'in RAM has been disabled'
          call verbose
        end if
      end if

      ! Store Wfunc for polarization in RAM
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        PRAM = .True.
      else
        PRAM = .False.
      end if

      ! Store Wfunc in file
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        WPFIL = .True.
      else
        WPFIL = .False.
      end if

      ! Store interpolation data for intensity
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        TIRAM = .True.
      else
        TIRAM = .False.
      end if

      ! Store interpolation data for polarization
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        TPRAM = .True.
      else
        TPRAM = .False.
      end if

      ! Change RAM if yes file
      if(WPFIL.and.PRAM)then
        PRAM = .False.
        if (pid.eq.0) then
          umsg = ' - You have chosen to store the '// &
                 'redistribution in a file, so storing them '// &
                 'in RAM has been disabled'
          call verbose
        end if
      end if

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
      if(Input%dws.eq.'NUM')read(100,*) Input%dw

      ! Minimum expected temperature
      read(100,*,err=1100) Input%minT

      ! Maximum expected velocity
      read(100,*,err=1100) Input%maxV

      ! Number of CPU per 1.5D column
      read(100,*,err=1100) Input%rt_group_n

      ! Force unmagnetized
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%unmagnetized = .True.
      else
        Input%unmagnetized = .False.
      end if

      ! Force static
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%static = .True.
      else
        Input%static = .False.
      end if

      ! Force static intensity
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%static_int = .True.
      else
        Input%static_int = .False.
      end if

      ! Skip star disk
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%skip_disk = .True.
      else
        Input%skip_disk = .False.
      end if

      ! Index of first iteration
      read(100,*,err=1100) Input%iter_min
      read(100,*,err=1100) Input%iteri_min

      ! Maximum number of iterations
      read(100,*,err=1100) Input%iter_max
      read(100,*,err=1100) Input%iteri_max

      ! Order of emissivity
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%iter_ord = 2
        PRD = .True.
      else
        Input%iter_ord = 1
        PRD = .False.
      end if

      ! Maximum relative change to consider convergence
      read(100,*,err=1100) Input%mrc_i
      read(100,*,err=1100) Input%mrci_i
      read(100,*,err=1100) Input%mrc_p

      ! Iterations of radiation field
      read(100,*,err=1100) Input%iter_j

      ! PRD internal iterations
      read(100,*,err=1100) Input%iteri_prd
      if (Input%iteri_prd.lt.1) Input%iteri_prd = 1
      if (Input%iteri_prd.gt.9) Input%iteri_prd = 9
      read(100,*,err=1100) Input%mrci_r

      ! NG acceleration
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%NG = .True.
      else
        Input%NG = .False.
      endif

      ! NG acceleration order
      read(100,*,err=1100) Input%NG_ord
      if (Input%NG_ord.lt.1.or.Input%NG_ord.gt.5) then
        Input%NG = .False.
      end if

      ! NG acceleration delay
      read(100,*,err=1100) Input%NG_delay

      ! NG acceleration for intensity
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%NGI = .True.
      else
        Input%NGI = .False.
      endif

      ! NG acceleration order for intensity
      read(100,*,err=1100) Input%NGI_ord
      if (Input%NGI_ord.lt.1.or.Input%NGI_ord.gt.5) then
        Input%NGI = .False.
      end if

      ! NG acceleration delay for intensity
      read(100,*,err=1100) Input%NGI_delay

      ! PRD in intensity delay
      read(100,*,err=1100) Input%PRD_delay

      ! ALI method delay
      read(100,*,err=1100) Input%ALI_delay

      ! Append into the MRC files
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%appendMRC = .True.
      else
        Input%appendMRC = .False.
      endif
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%appendMRCI = .True.
      else
        Input%appendMRCI = .False.
      endif

      ! Type of broadcasting
      read(100,*,err=1100) ia
      if(ia.eq.0)then
        Input%altbcast = .False.
      else
        Input%altbcast = .True.
      end if

      ! Non physical Stokes and rho
      read(100,*,err=1100) Input%allownphys_stk
      read(100,*,err=1100) Input%allownphys_rho

      ! Solution Box
      if(run_mode.eq.1.or.run_mode.eq.-1) then
        allocate(Input%sol_box(4))
        read(100,*,err=1100) Input%sol_box
      else
        read(100,*,err=1100) ios
        read(100,*,err=1100) ios
        read(100,*,err=1100) ios
        read(100,*,err=1100) ios
      end if

      ! Steps between storing solution
      read(100,*,err=1100) Input%store_step
      if(Input%store_step.lt.1)then
        Input%store = .False.
        Input%store_step = 1
      else
        Input%store = .True.
      end if
      read(100,*,err=1100) Input%storei_step
      if(Input%storei_step.lt.1)then
        Input%storeI = .False.
        Input%storei_step = 1
      else
        Input%storeI = .True.
      end if

      ! Output contribution function
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%out_contr = .True.
      else
        Input%out_contr = .False.
      end if

      ! Output tau=1 height
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%out_tau1 = .True.
      else
        Input%out_tau1 = .False.
      end if

      ! Keep background
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_back = .True.
      else
        Input%keep_back = .False.
      end if

      ! Keep damp
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_damp = .True.
      else
        Input%keep_damp = .False.
      end if

      ! Keep a parameters
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_aparam = .True.
      else
        Input%keep_aparam = .False.
      end if

      ! Keep cols
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_cols = .True.
      else
        Input%keep_cols = .False.
      end if

      ! Keep atmo
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_atmo = .True.
      else
        Input%keep_atmo = .False.
      end if

      ! Keep pop
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_pop = .True.
      else
        Input%keep_pop = .False.
      end if

      ! Keep dep
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_dep = .True.
      else
        Input%keep_dep = .False.
      end if

      ! Keep rhoKQ
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_rhoKQ = .True.
      else
        Input%keep_rhoKQ = .False.
      end if

      ! Keep JKQ
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_JKQ = .True.
      else
        Input%keep_JKQ = .False.
      end if

      ! Keep stokes in quadrature
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_stokesQ = .True.
      else
        Input%keep_stokesQ = .False.
      end if

      ! Keep JKQnu
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_jkqnu = .True.
      else
        Input%keep_jkqnu = .False.
      end if

      ! Keep MRC
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_MRC = .True.
      else
        Input%keep_MRC = .False.
      end if

      ! Keep COL log
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_coll = .True.
      else
        Input%keep_coll = .False.
      end if

      ! Keep MPI log
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_mpil = .True.
      else
        Input%keep_mpil = .False.
      end if

      ! Keep MPI detail log
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%keep_mpidl = .True.
      else
        Input%keep_mpidl = .False.
      end if

      ! Limit Stokes output
      read(100,*,err=1100) Input%lim_stk%nran
      if (Input%lim_stk%nran.gt.0) then
        allocate(Input%lim_stk%doub(2,Input%lim_stk%nran))
        do i1=1,Input%lim_stk%nran
          read(100,*,err=1100) Input%lim_stk%doub(1,i1)
          read(100,*,err=1100) Input%lim_stk%doub(2,i1)
        end do
      end if

      ! Limit contribution output
      read(100,*,err=1100) Input%lim_ctr%nran
      if (Input%lim_ctr%nran.gt.0) then
        allocate(Input%lim_ctr%doub(2,Input%lim_ctr%nran))
        do i1=1,Input%lim_ctr%nran
          read(100,*,err=1100) Input%lim_ctr%doub(1,i1)
          read(100,*,err=1100) Input%lim_ctr%doub(2,i1)
        end do
      end if

      ! Limit tau output
      read(100,*,err=1100) Input%lim_tau%nran
      if (Input%lim_tau%nran.gt.0) then
          allocate(Input%lim_tau%doub(2,Input%lim_tau%nran))
        do i1=1,Input%lim_tau%nran
          read(100,*,err=1100) Input%lim_tau%doub(1,i1)
          read(100,*,err=1100) Input%lim_tau%doub(2,i1)
        end do
      end if

      ! Limit cols-TT output
      read(100,*,err=1100) Input%lim_cols_tt%nran
      if (Input%lim_cols_tt%nran.gt.0) then
        allocate(Input%lim_cols_tt%indx(3,Input%lim_cols_tt%nran))
        do i1=1,Input%lim_cols_tt%nran
          read(100,*,err=1100) Input%lim_cols_tt%indx(1,i1)
          read(100,*,err=1100) Input%lim_cols_tt%indx(2,i1)
          read(100,*,err=1100) Input%lim_cols_tt%indx(3,i1)
        end do
      end if

      ! Limit cols-ll output
      read(100,*,err=1100) Input%lim_cols_ll%nran
      if (Input%lim_cols_ll%nran.gt.0) then
        allocate(Input%lim_cols_ll%indx(3,Input%lim_cols_ll%nran))
        do i1=1,Input%lim_cols_ll%nran
          read(100,*,err=1100) Input%lim_cols_ll%indx(1,i1)
          read(100,*,err=1100) Input%lim_cols_ll%indx(2,i1)
          read(100,*,err=1100) Input%lim_cols_ll%indx(3,i1)
        end do
      end if

      ! Limit damp output
      read(100,*,err=1100) Input%lim_damp%nran
      if (Input%lim_damp%nran.gt.0) then
        allocate(Input%lim_damp%indx(2,Input%lim_damp%nran))
        do i1=1,Input%lim_damp%nran
          read(100,*,err=1100) Input%lim_damp%indx(1,i1)
          read(100,*,err=1100) Input%lim_damp%indx(2,i1)
        end do
      end if

      ! Limit back output
      read(100,*,err=1100) Input%lim_back%nran
      if (Input%lim_back%nran.gt.0) then
          allocate(Input%lim_back%doub(2,Input%lim_back%nran))
        do i1=1,Input%lim_back%nran
          read(100,*,err=1100) Input%lim_back%doub(1,i1)
          read(100,*,err=1100) Input%lim_back%doub(2,i1)
        end do
      end if

      ! Limit pop output
      read(100,*,err=1100) Input%lim_pop%nran
      if (Input%lim_pop%nran.gt.0) then
          allocate(Input%lim_pop%indx(2,Input%lim_pop%nran))
        do i1=1,Input%lim_pop%nran
          read(100,*,err=1100) Input%lim_pop%indx(1,i1)
          read(100,*,err=1100) Input%lim_pop%indx(2,i1)
        end do
      end if

      ! Recalculate electron density
      read(100,*,err=1100) Input%redo_ne

      ! Update atmospheric model after completing
      ! intensity calculations, and what to store
      read(100,*,err=1100) Input%update_atmos

      ! Keep input Hydrogen densities (atmos or atmo) after
      ! equation of state
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%protect_H = .True.
      else
        Input%protect_H = .False.
      end if

      ! Do not change the atomic populations in the
      ! chemical equilibrium
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%chem_protect_all = .True.
      else
        Input%chem_protect_all = .False.
      end if

      ! Write general performance
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%g_perf = .True.
      else
        Input%g_perf = .False.
      endif

      ! Write MPI performance
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%mpi_perf = .True.
      else
        Input%mpi_perf = .False.
      endif

      ! Read MPI detailed file
      read(100,*,err=1100) Input%MPIdetail

      ! Read previous performance file
      read(100,*,err=1100) Input%operform

      ! Read if skipping first iteration in performance
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        Input%IWskip = .True.
      else
        Input%IWskip = .False.
      endif

      ! Verbosity
      read(100,*,err=1100) cdump
      if(cdump.eq.'Y')then
        verbosity = .True.
      else
        verbosity = .False.
      endif
      verbosef = trim(Input%folder)//'/verbose'

      ! This MUST be false if not inverting
      Input%force_inv_freq = .False.

      !
      ! Inversion mode
      !
      if (Input%run_mode.eq.-1) then

        ! Allocations
        allocate(Input%Node_type(Input%nvar))
        allocate(Input%Node(Input%nvar))
        allocate(Input%Num_nodes(Input%nvar))
        allocate(Input%Nodes_flags(Input%nvar))
        allocate(Input%Nodes_Regul(Input%nvar))
        allocate(Input%Indx_regul(Input%nvar))
        allocate(Input%Regul_weight(Input%nvar))
        allocate(Input%Scal(Input%nvar))
        allocate(Input%Perturb(Input%nvar))
        allocate(Input%min_rel_Pert(Input%nvar))

        ! Verbosity level
        read(100,*,err=1100) vlevel
        read(100,*,err=1100) slevel

        ! Data file
        read(100,'(A)',err=1100) Input%Filename_ob

        ! Type of inversion
        read(100,*,err=1100) Input%Type_Inversion

        ! Automatic Stokes weights
        read(100,'(A)',err=1100) cdump
        if(cdump.eq.'Y')then
          Input%auto_weight = .True.
        else
          Input%auto_weight = .False.
        endif

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

            ! For each set of weights
            do i1=1,Input%Num_Weight
              read(100,*,err=1100) Input%weight(:,i1)
            end do ! Sets of weights

          end if ! Numeric input

          ! Check if file
          read(100,*,err=1100) i1
          Input%linv_weight = i1.ne.0

          ! If file
          if (Input%linv_weight) &
            read(100,'(A)',err=1100) Input%inv_weight

        end if ! No automatic weights

        ! Inversion initialization
        read(100,'(A)',err=1100) Input%Inv_init

        ! Centered derivative
        read(100,'(A)',err=1100) cdump
        if(cdump.eq.'Y')then
          Input%centered = .True.
        else
          Input%centered = .False.
        endif

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

            ! And read positions
            read(100,*,err=1100) Input%Node(i1)%H

          end if

        end do ! For each variable

        ! Type of magnetic field vector
        read(100,*,err=1100) Input%btype

        ! Type of velocity vector
        read(100,*,err=1100) Input%vtype

        ! Flag each variable as changing or not
        do i1=1,Input%nvar

          ! Read string and assign the negative
          read(100,'(A)',err=1100) cdump
          if (cdump.eq.'Y') then
            Input%Nodes_Flags(i1) = .False.
          else
            ! True only if there are nodes
            Input%Nodes_Flags(i1) = Input%Num_nodes(i1).gt.0
          end if

        end do ! For each variable

        ! Flag to correct positions
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Pos_Correction = .True.
        else
          Input%Pos_Correction = .False.
        end if

        ! Regularization type
        do i1=1,Input%nvar

          ! Read index
          read(100,*,err=1100) Input%Indx_regul(i1)
          ! If not none
          if (Input%Indx_regul(i1).gt.0) then
            read(100,*,err=1100) Input%Regul_weight(i1)
          else
            Input%Regul_weight(i1) = 0d0
          end if

        end do ! Variables

        ! Regularization limit
        read(100,*,err=1100) Input%Regul_Limit

        ! CHI2 Threshold
        read(100,*,err=1100) Input%Threshold_chisq

        ! Inversion MRC
        read(100,*,err=1100) Input%Chisq_fraction

        ! Read SVD type
        read(100,*,err=1100) Input%SVD_type

        ! Threshold for SVD
        read(100,*,err=1100) Input%Threshold_svd

        ! Type of calculation for pressure
        read(100,*,err=1100) i1
        if (i1.eq.1) Input%hydroeq = .True.

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

          ! For each range
          do i1=1,ios

            ! Save size data (redundant)
            Input%lim_fwhm(i1)%nn = ios

            ! Read type
            read(100,*,err=1100) i2

            ! Gaussian?
            Input%lim_fwhm(i1)%gaussian = i2.eq.0

            ! If gaussian
            if (Input%lim_fwhm(i1)%gaussian) then

              ! Allocate doubles
              allocate(Input%lim_fwhm(i1)%doub(3))

              ! Read entry
              read(100,*,err=1100) Input%lim_fwhm(i1)%doub

            ! If NOT gaussian
            else

              ! Allocate doubles
              allocate(Input%lim_fwhm(i1)%doub(2))

              ! Read entry
              read(100,*,err=1100) Input%lim_fwhm(i1)%doub

              ! Get file
              read(100,'(A)',err=1100) Input%fwhm_fil(i1)

            end if ! Gaussian

          end do

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

              ! Save number
              Input%Node(i1)%nebound = ios

              ! Read boundaries
              do i2=1,ios

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

        ! Invert fractional Stokes
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Fractional = .True.
        else
          Input%Fractional = .False.
        end if

        ! Range of tau for atmosphere
        read(100,*,err=1100) Input%Tau_Range

        ! Modifications on the stratification of the
        ! model atmosphere
        read(100,*,err=1100) i1

        ! If input data for modification of stratification
        if (i1.gt.0) then

          ! Allocate
          allocate(Input%Atmo_strat(3,i1))

          ! For each input
          do i2=1,i1

            ! Read parameters
            read(100,*,err=1100) Input%Atmo_strat(:,i2)

          end do ! For each input

        end if ! Input data for stratification

        ! Broyden in LM
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Broyden = .True.
        else
          Input%Broyden = .False.
        end if

        ! Index of LM method
        read(100,*,err=1100) Input%LM_Method

        ! Range of lambda for LM
        read(100,*,err=1100) Input%Lam_Range

        ! Accepted lambda factor
        read(100,*,err=1100) Input%factoraccept

        ! Rejected lambda factor
        read(100,*,err=1100) Input%factorreject

        ! Project the B vector (negative strength values become
        ! negative polarity)
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Projection = .True.
        else
          Input%Projection = .False.
        end if

        ! Initialize from previous solution when computing
        ! response functions
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Popuinit = .True.
        else
          Input%Popuinit = .False.
        end if

        ! Neglect sigma
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Sigma_neglect = .True.
        else
          Input%Sigma_neglect = .False.
        end if

        ! Keep RF
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%Keep_RF = .True.
        else
          Input%Keep_RF = .False.
        end if

        ! Force inversion frequencies
        read(100,'(A)',err=1100) cdump
        if (cdump.eq.'Y') then
          Input%force_inv_freq = .True.
        else
          Input%force_inv_freq = .False.
        end if

      end if ! Inversion mode

      close(100)

      ! Correct NGI_delay and iterations if doing PRD
      if (PRD) then
        Input%NGI_delay = Input%NGI_delay + Input%PRD_delay - 1
        Input%iteri_max = Input%iteri_max + Input%PRD_delay - 1
      else
        Input%PRD_delay = 0
      end if

      ! Get ID in character
      write(CPUC,'(I0.5)') pid

      ! Store name for error files
      errorf = trim(Input%folder)//'/ERROR'//CPUC

      ! Control that everything went fine
      call control

      ! Delete the temporal input file and communicate success
      if(pid.eq.0) then
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
