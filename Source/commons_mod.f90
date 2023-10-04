      !> Share variables of common use
      module commons_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/17/2017
!  Last version:
!     09/29/2023 V3.0.8
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.8 - Added Kradl (TdPA)
!
!     09/08/2023:    V3.0.7 - Added vlevel and slevel (TdPA)
!
!     08/07/2023:    V3.0.6 - Added LVIRAM, LVPRAM, and nLTEl (TdPA)
!
!     03/08/2023:    V3.0.5 - Added verbosefv and ninv_mode (TdPA)
!
!     02/14/2023:    V3.0.4 - Added aixali and AVI (TdPA)
!
!     11/10/2022:    V3.0.3 - Added force_asym (TdPA)
!
!     10/25/2022:    V3.0.2 - Added Rz0, Rz1, and Rnz (TdPA)
!
!     07/27/2022:    V3.0.1 - Now we use 'use mpi' or 'use mpi_f08'
!                             (configurable in configure script)
!                             instead of 'import mpif.h', which is
!                             stricter and standard (TdPA)
!                           - Common blocks have been removed in
!                             favor of module variables (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case added the
!                             variables gnproc, gpid, MPI_COMM_RT,
!                             MPI_COMM_CTRL, run_mode, nhori, and
!                             icoords (TdPA)
!
!     03/23/2021:    V2.0.1 - Increased size of umsg (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added nthread (TdPA)
!
!     02/17/2021:   V1.1.12 - Added fcol_transfer (TdPA)
!
!     01/12/2021:   V1.1.11 - Added KSTK (TdPA)
!
!     09/20/2020:   V1.1.10 - Added errorf (TdPA)
!
!     09/11/2020:    V1.1.9 - Added TIRAM and TPRAM (TdPA)
!
!     07/10/2020:    V1.1.8 - Added KcutAB (TdPA)
!
!     06/01/2020:    V1.1.7 - Added NCHLT (TdPA)
!
!     11/20/2019:    V1.1.6 - Added VOITY (TdPA)
!
!     11/19/2019:    V1.1.5 - Added WIFIL and WPFIL (TdPA)
!
!     11/12/2019:    V1.1.4 - Added VIFIL, VPFIL, PIRAM, and
!                             offlimit (TdPA)
!
!     09/13/2019:    V1.1.3 - Added ztau logical variable, that
!                             indicates if the height axis is in tau
!                             units (TdPA)
!
!     05/31/2019:    V1.1.2 - Removed emerging variable (TdPA)
!
!     05/08/2019:    V1.1.1 - Added nphysS, nphysR, vaborted, nxJ,
!                             nxS, and nxL (TdPA)
!
!     02/20/2019:    V1.1.0 - Removed UIO, added laborted, verbosity,
!                             verbosef, umsg, and urou (TdPA)
!                           - Bugfix: RTaxial was not in the common
!                             blocks (TdPA)
!
!     30/01/2019:    V1.0.8 - Added RTaxial (TdPA)
!
!     08/04/2018:    V1.0.7 - Added RLIM integer variable and VIRAM
!                             and VPRAM logical variables (TdPA)
!
!     08/03/2018:    V1.0.6 - Added Krad integer variable (TdPA)
!
!     09/22/2017:    V1.0.5 - Added Kcut integer variable (TdPA)
!
!     06/16/2017:    V1.0.4 - Added PRAM logical variable and changed
!                             RAM to IRAM (TdPA)
!
!     06/12/2017:    V1.0.3 - Added RAM logical variable (TdPA)
!
!     06/08/2017:    V1.0.2 - Removed ALI logical variable (TdPA)
!
!     05/08/2017:    V1.0.1 - Added ALI logical variable (TdPA)
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
!    Common of frequently used numbers
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

      ! File for verbosity, universal message, file for inversion
      ! verbosity
      character(len=500), save:: verbosef, errorf, verbosefv
      character(len=650), save:: umsg

      ! Universal routine
      character(len=20), save:: urou

      ! Stimulated emission, angle-averaged redistribution, angle-
      ! dependent distribution with angle-averaged input, axial
      ! symmetry, dynamics, PRD computation, use RAM to store Wfunc
      ! for intensity, use RAM to store Wfunc for polarization, use
      ! RAM to store Voigt profiles for intensity, use RAM to store
      ! Voigt profiles for polarization, if radiation transfer is
      ! axial (I+Q/I-Q), status of each CPU, verbosity, non physical
      ! Stokes allowed, non physical rhoKQ allowed, if aborted has
      ! sent a message already, using Voigt intensity file, using
      ! Voigt polarization file, storing photoionization quantities
      ! using redistribution intensity file, using redistribution
      ! polarization file, non-coherent lower term approx, apply
      ! Kcut in absorb1, save interpolation data for intensity and
      ! for polarization, store the Stokes parameters in the
      ! Solution file, if forbidden collisions can transfer K!=0,
      ! force asymmetry input, axial symmetry in intensity problem,
      ! angle averaged force in intensity problem, if synthesis not
      ! in inversion mode important for verbosity output, use RAM
      ! to store Voigt profiles for LTE lines in intensity, use
      ! RAM to store Voigt profiles for LTE lines in polarization
      logical, save:: stm, AV, tbAD, axial, dyn, PRD, IRAM, &
                      PRAM, VIRAM, VPRAM , RTaxial, laborted, &
                      verbosity, nphysS, nphysR, vaborted, ztau, &
                      VIFIL, VPFIL, PIRAM, WIFIL, WPFIL, NCHLT, &
                      KcutAB, TIRAM, TPRAM, KSTK, fcol_transfer, &
                      force_asym, axiali, AVI, ninv_mode, LVIRAM, &
                      LVPRAM

      ! Number of processors, ID of processor, integer for error
      ! output, number of threads
      integer, save:: nproc, pid, ierr, nthread

      ! Global number of preocessors, global ID of processor
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
      ! maximum multipole K, maximum dimension of SEE, maximum number
      ! of transitions, maximum number of FS transitions, maximum
      ! number of photoionization transitions, maximum number of
      ! blends, global K cut limit, global K rad limit, global K rad
      ! limit for lines, limit on RAM for profiles in MB, maximum J
      ! number maximum L number, maximum S number, type of intensity
      ! Voigt profile, type of run, number of LTE lines, verbosity
      ! level and shut-up level
      integer, save:: nA, nAb, nM, nZ, nfreq, nkx, nxdim, nxtran, &
                      nxt, nxphot, nxb, Kcut, Krad, Kradl, RLIM, &
                      nxJ, nxS, nxL, VOITY, run_mode, nLTEl, &
                      vlevel, slevel

      ! Height limits due to restrictions, lower and higher and
      ! total
      integer, save:: Rz0, Rz1, Rnz

      ! Global boundary indexes for Stokes
      integer, save:: giz0, giz1

      ! Current coordinates for 1.5D
      integer, dimension(3), save:: icoords

      ! Maximum number for MPI_OFFSET_KIND
      double precision, save:: offlimit

      end module commons_mod
