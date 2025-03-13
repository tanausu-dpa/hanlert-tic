      !> Share variables of common use
      module commons_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     17/04/2017
!  Last version:
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - Added nJs (TdPA)
!                           - Removed variables related to store
!                             profiles in files, interpolation
!                             data storage, and threads (TdPA)
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
!    MPI library common access and global variables
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Include MPI header
#ifdef oldmpi
      use mpi
#else
      use mpi_f08
#endif

      ! File for verbosity, for error messages, and for inversion
      ! verbosity
      character(len=500), save:: verbosef, errorf, verbosefv

      ! Universal message
      character(len=650), save:: umsg

      ! Universal routine name
      character(len=20), save:: urou

      ! Stimulated emission, angle-averaged redistribution, angle-
      ! dependent distribution with angle-averaged input, axial
      ! symmetry, dynamics, PRD computation, use RAM to store Wfunc
      ! for intensity, use RAM to store Wfunc for polarization, use
      ! RAM to store Voigt profiles for intensity, use RAM to store
      ! Voigt profiles for polarization, if radiation transfer is
      ! axial (I+Q/I-Q), status of each CPU, verbosity, non physical
      ! Stokes allowed, non physical rhoKQ allowed, if aborted has
      ! sent a message already, if vertical scale is in optical depth,
      ! storing photoionization quantities, non-coherent lower term
      ! approx., apply Kcut in absorb1, store the Stokes parameters
      ! in the Solution file, if forbidden collisions can transfer
      ! K!=0, force asymmetry input, axial symmetry in intensity
      ! problem, angle averaged force in intensity problem, if
      ! synthesis not in inversion mode (important for verbosity
      ! output), use RAM to store Voigt profiles for LTE lines in
      ! intensity, use RAM to store Voigt profiles for LTE lines in
      ! polarization
      logical, save:: stm, AV, tbAD, axial, dyn, PRD, IRAM, &
                      PRAM, VIRAM, VPRAM , RTaxial, laborted, &
                      verbosity, nphysS, nphysR, vaborted, ztau, &
                      PIRAM, NCHLT, KcutAB, KSTK, fcol_transfer, &
                      force_asym, axiali, AVI, ninv_mode, LVIRAM, &
                      LVPRAM

      ! Global integer to use as error in MPI routines
      integer:: ierr

      ! Number of processors and ID of processor for RT groups
      integer, save:: nproc, pid

      ! Number of processors and ID of processor for WORLD
      integer, save:: gnproc, gpid

      ! MPI communicators and status
#ifdef oldmpi
      integer, save:: MPI_COMM_RT, MPI_COMM_CTRL, MPI_STAT
#else
      type(MPI_comm), save:: MPI_COMM_RT, MPI_COMM_CTRL
      type(MPI_Status):: MPI_STAT
#endif

      ! MPI variables for unit handling (no need to keep them)
#ifdef oldmpi
      integer:: funit,info
#else
      type(MPI_Info):: info
      type(MPI_File):: funit
#endif

      ! Number of active atoms, number of background atoms, number
      ! of molecules, number of heights, number of frequencies,
      ! maximum multipole K, maximum dimension of SEE, total number
      ! of transitions, total number of FS transitions, total
      ! number of photoionization transitions, maximum number of
      ! blends, global K cut limit, global K rad limit, global K rad
      ! limit for lines, limit on RAM to allocate in MB, maximum J
      ! number, maximum L number, maximum S number, type of intensity
      ! Voigt profile, type of run, number of LTE lines, verbosity
      ! level, verbosity shut-up level, number of stored J symbols
      integer, save:: nA, nAb, nM, nZ, nfreq, nkx, nxdim, nxtran, &
                      nxt, nxphot, nxb, Kcut, Krad, Kradl, RLIM, &
                      nxJ, nxS, nxL, VOITY, run_mode, nLTEl, &
                      vlevel, slevel, nJs

      ! Height limits due to restrictions, lower, higher, and total
      ! Height limit restriction for PRD
      integer, save:: Rz0, Rz1, Rnz, Rz1_PRD

      ! Global boundary indexes for Stokes
      integer, save:: giz0, giz1

      ! Current coordinates for 1.5D
      integer, dimension(3), save:: icoords

      ! Maximum number for MPI_OFFSET_KIND
      double precision, save:: offlimit

      !
      ! Allocated RAM
      !

      ! Pre-calculated quantities for photoionizations
      double precision, save:: PRAMc

      ! Voigt profiles
      double precision, save:: VRAMc

      ! Redistribution function
      double precision, save:: WRAMc

      ! Memoization
      double precision, save:: ERAMc

      ! Radiation field
      double precision, save:: RRAMc

      ! Self-consistent solution
      double precision, save:: SRAMc

      ! Background quantities
      double precision, save:: BRAMc

      ! Predicted for the RT problem
      double precision, save:: TRAMc

      ! 1st order part of the PRD emissivity
      double precision, save:: ORAMc

      ! Frequency data for PRD emissivity
      double precision, save:: FRAMc

      ! Expected RAM for normalization data
      double precision, save:: DRAMc,DRAM2c

      ! Rest of variables (miscellaneous)
      double precision, save:: MRAMc

      end module commons_mod
