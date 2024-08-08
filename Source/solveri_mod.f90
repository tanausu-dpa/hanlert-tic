      !> Solution of NLTE problem of the first kind
      module solveri_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Contributors:
!     Hao Li (IAC)
!  Start:
!     04/20/2017
!  Last version:
!     08/08/2024 V3.0.18
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/08/2024:   V3.0.17 - Force at least two iterations if doing
!                             AD starting from AA (TdPA)
!
!     02/23/2024:   V3.0.17 - Manage the force ALI option (TdPA)
!                           - Change the format and header for the
!                             MRC messages for the new limit of
!                             possible PRD iterations (TdPA)
!                           - Now Input%allownphys_pop decides if
!                             negative populations are allowed (TdPA)
!                           - Added bool argument to RTCoeffi to
!                             specify if data2*, rLine*, and
!                             rPhoto* need to be returned (TdPA)
!                           - Now the J iterations can include the
!                             atomic transitions (TdPA)
!                           - Now the J iterations can stop if a
!                             convergence has been achieved (TdPA)
!
!     02/20/2024:   V3.0.16 - Bugfix: Wrong initial index in the
!                             tau1 calculation when restricting
!                             heights (TdPA)
!
!     11/14/2023:   V3.0.15 - Added termination paths in case of
!                             error (TdPA)
!
!     10/31/2023:   V3.0.14 - Call scattering_manage in angular
!                             loops in solvers (TdPA)
!                           - Call get_scattering_los when computing
!                             emergent LOS (TdPA)
!
!     10/16/2023:   V3.0.13 - Made LTElines allocatable to satisfy
!                             memory warnings (TdPA)
!                           - Added number of azimuths variable to
!                             the JKQgen* calls. This is used to
!                             allocate Stokes with its real size,
!                             which can be different to Geom%nPh
!                             if running a two-step solution (TdPA)
!                           - Bugfix: Wrong order of indexes when
!                             copying axial intensity in the
!                             JKQgen_alt routine (TdPA)
!
!     08/28/2023:   V3.0.12 - Wrong check in emergenceI_serial
!                             to free the RT memory (TdPA)
!
!     08/17/2023:   V3.0.11 - When aborting, go to the memory
!                             cleaning instead of returning (TdPA)
!                           - Moved control call before memory
!                             cleaning (TdPA)
!
!     08/07/2023:   V3.0.10 - Added arguments for LTE lines (TdPA)
!                           - Added solveI, emergentI, solveJ, and
!                             JKQgenerate (TdPA)
!                           - Regarding the issue in V3.0.9, I
!                             removed the fix. The reason is that
!                             there is already a MPI_WAIT for
!                             request5 before contr_s is going to
!                             be overwritten. Adding an MPI_WAIT
!                             after the call to MPI_ISEND makes it
!                             fully equivalent to a MPI_SEND, losing
!                             any advantage of asynchronous
!                             communication. My naive guess about
!                             the issue would be that the laptop
!                             does not have the resources to handle
!                             the amount of communication for the
!                             size of the problem. If the issue
!                             persist, it requires further
!                             investigation, but the patch is
!                             performance-affecting (TdPA)
!
!     07/31/2023:    V3.0.9 - It works well on diablos/dracarys.
!                             However, in my laptop (gfortran
!                             13.1.0 openmpi 4.1.5 macos 12.4),it
!                             reports an error. A mpi_wait can solve
!                             the problem, although I do not
!                             understand why (HL)
!
!     07/03/2023:    V3.0.8 - If not doing ALI, point the relevant
!                             pointers to the profile so they are not
!                             undefined (TdPA)
!                           - chi500 is now it the Atmo structure and
!                             not in Cont (TdPA)
!                           - In the inversion, the result of
!                             the emergence is stored in the SolF
!                             structure, setup by calling settau,
!                             setstkI, and setctrI (TdPA)
!
!     04/11/2023:    V3.0.7 - Fix a typo (HL)
!
!     03/21/2023:    V3.0.6 - Bugfix: Wrong specification of CPU rank
!                             when running in 1.5D to output the final
!                             MRC (TdPA)
!                           - Bugfix: In the serial emergence, the
!                             sizes of J00C and J00 were interchanged,
!                             leading to obvious memory issues (TdPA)
!                           - Bugfix: Ensured the OpenMP version
!                             compiles after some years of changes,
!                             albeit did not test it works (TdPA)
!
!     02/14/2023:    V3.0.5 - Also write Stokes file in inversion
!                             mode (HL)
!                           - Now AVI and axiali exist, which do not
!                             imply AV or axial (TdPA)
!                           - In the JKQgen routines, account for
!                             axial and axiali being different in
!                             general (TdPA)
!                           - Bugfix: One height loop was not changed
!                             to account for the axis restriction.
!                             Only relevant with NG, restricting the
!                             axis, and with the need to store the
!                             full stokes array (TdPA)
!
!     10/26/2022:    V3.0.4 - Bugfix: Wrong height axis size when
!                             unrolling indexes in master process to
!                             call the JKQ integrator (TdPA)
!                           - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.3 - Implemented the restriction of the
!                             height axis (TdPA)
!                           - Made sure to clean every pointer (TdPA)
!                           - Added the option to retrieve emergint
!                             profiles in emergence and
!                             emergence_serial routines (TdPA)
!                           - New way of dealing with Frec%exu in
!                             the serial case (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!                           - funit is now global (TdPA)
!
!     06/30/2022:    V3.0.1 - Bugfixes: Quick changes while preparing
!                             the last commit led to mistakes. Solved
!                             two typos (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Changed the MPI communicator from
!                                MPI_COMM_WORLD to MPI_COMM_RT.
!                              o Added returns for when we need to
!                                abort the run.
!                              o Atmo%v has changed to Atmo%vx,%vy,
!                                and %vz.
!                              o Changed the arguments to writesol,
!                                writetau, writectr, and writestk.
!                              o Distinguished how MRC is written
!                                depending on the type of run. It
!                                is also optional.
!                             (TdPA)
!
!     04/07/2022:    V2.0.3 - Added an argument to writesolI (TdPA)
!
!     09/30/2021:    V2.0.2 - Ensure no zero size allocation (TdPA)
!
!     03/23/2021:    V2.0.1 - Bugfix: Was not possible to append in
!                             the MRCI file (TdPA)
!                           - Made some changes to improve error
!                             handling (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Completely changed the Master section
!                             and introduced OpenMP (TdPA)
!
!     02/12/2021:   V1.8.12 - Call report_mpi_timeI in the solver
!                             subroutines (TdPA)
!
!     02/09/2021:   V1.8.11 - Removed some remaning checks when
!                             allocating (TdPA)
!                           - Bugfix: Add checks for presence
!                             of spectral lines when building the
!                             Lambda operator in slave processes. It
!                             could result in out of bounds in some
!                             cases (TdPA)
!
!     01/12/2021:   V1.8.10 - Using KSTK to determine how to store
!                             the Stokes parameters (TdPA)
!
!     11/12/2020:    V1.8.9 - Added option to force a number of CRD
!                             iterations when doing PRD (TdPA)
!
!     09/11/2020:    V1.8.8 - Changes to avoid unallocated arrays
!                             when the photoionization exponential
!                             is not stored, either because it is
!                             specified or because there are no
!                             photoionizations (TdPA)
!                           - Remove some remains of status check
!                             in some deallocations (TdPA)
!
!     07/31/2020:    V1.8.7 - Change the p_exu pointer so it only
!                             points to the stored data if there
!                             is something to point to (TdPA)
!                           - Now initializing J00S (TdPA)
!
!     07/22/2020:    V1.8.6 - Introduced flags for the NG acceleration
!                             independent of the ones in the solver
!                             module for polarization (TdPA)
!                           - Bugfix: When an intensity solution with
!                             AA was loaded to do AD polarization, the
!                             JKQgen* routines could not correctly
!                             deal with it (TdPA)
!
!     07/10/2020:    V1.8.5 - tau1 needs to be allocated in the slaves
!                             also when computing only the
!                             contribution function (TdPA)
!
!     03/05/2020:    V1.8.4 - If doing PRD and the first iteration
!                             in CRD gives a very small MRC, do
!                             at least a second one (TdPA)
!
!     01/16/2020:    V1.8.3 - Bugfix: The wrong variable was checked
!                             to decide if exponentials were stored
!                             in RAM (TdPA)
!
!     11/19/2019:    V1.8.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     11/13/2019:    V1.8.1 - Atom%popu needs update when NG is
!                             applied (TdPA)
!
!     09/13/2019:    V1.8.0 - Now Stokes needs to be stored if the
!                             calculation if dynamic and angle-
!                             averaged too (TdPA)
!                           - The vertical scale can now be optical
!                             depth (TdPA)
!
!     08/19/2019:    V1.7.1 - Ignoring J00S. Commented the definition,
!                             communications, and changes J00S to
!                             J00 in the calls of writesol and seei.
!                             It is still used in JKQgen (TdPA)
!
!     05/31/2019:    V1.7.0 - Changed the dimensionality of the
!                             profile and ratio variables. Now it runs
!                             sequentially on atoms, transitions and
!                             frequencies to save memory and reduce
!                             the size of data shared through MPI
!                             messages (TdPA)
!
!     05/08/2019:    V1.6.6 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!
!     04/16/2019:    V1.6.5 - Implemented back the broadcasting from
!                             the mpi libraries (TdPA)
!
!     03/20/2019:    V1.6.4 - Now resets the MRCI file (TdPA)
!
!     03/18/2019:    V1.6.2 - Added logic to allow for precomputed
!                             exponentials for J integrals (TdPA)
!
!     03/13/2019:    V1.6.2 - Now can skip ALI for Input%ALI_delay
!                             iterations (TdPA)
!
!     03/06/2019:    V1.6.1 - Only the master says that the emergent
!                             direction is finished (TdPA)
!
!     02/20/2019:    V1.6.0 - New verbosity (TdPA)
!                           - Using unit 800 for MRC file (TdPA)
!                           - MRC file is not kept open (TdPA)
!                           - When possible, check communications by
!                             MPI (TdPA)
!
!     11/28/2018:    V1.5.1 - Removed duplicated writing block in
!                             solverI_serial (TdPA)
!
!     09/20/2018:    V1.5.0 - Added NG acceleration (TdPA)
!
!     09/06/2018:    V1.4.1 - Added argument to writesolI (TdPA)
!
!     09/04/2018:    V1.4.0 - Added solverI_alt, solverJ_alt, and
!                             JKQgen_alt (TdPA)
!
!     05/16/2018:    V1.3.5 - Changed azimuthal dimension of Stokes
!                             from nPh2 to nPh (TdPA)
!                           - Bugfix: Stokes was not being updated
!                             when going from the CRD first step
!                             to the angle-dependent case when not
!                             reading a previous solution, in the
!                             serial mode (TdPA)
!                           - Removed unnecessary variable op (TdPA)
!
!     11/28/2017:    V1.3.4 - Added option to skip rho00 multiterm
!                             correction (TdPA)
!
!     11/01/2017:    V1.3.3 - Changed MRC given information (TdPA)
!
!     10/03/2017:    V1.3.2 - Changed inputs of Termprof for to be
!                             consistent with its changes (TdPA)
!                           - Moved the exit of the iteration loop
!                             after the request checking (TdPA)
!
!     09/28/2017:    V1.3.1 - Bugfix: Forgot to add request checks
!                             for isend in solveri in PRD iterations,
!                             the master was destroying the data
!                             before some slaves received it (TdPA)
!
!     09/27/2017:    V1.3.0 - Implemented send_tree algorithm to
!                             avoid the terrible scaling of the native
!                             mpi_bcast (TdPA)
!
!     09/08/2017:    V1.2.1 - Avoid divisions by 0 (TdPA)
!
!     08/30/2017:    V1.2.0 - Change how the received data is
!                             reshaped, now it used pointers and there
!                             is no copy (TdPA)
!
!     08/28/2017:    V1.1.0 - Increased the buffer size, so each
!                             process does never wait the master to
!                             receive, and changes related with this
!                             modification (TdPA)
!
!     08/09/2017:   V1.0.17 - Changed loop order in flat part of
!                             JKQgen and JKQgen_serial (TdPA)
!                           - Do not do a CRD iterations if loading
!                             a file (TdPA)
!
!     07/21/2017:   V1.0.16 - Bugfix: Everyone have to wait until
!                             all communications have been finished
!                             before starting another PRD iteration
!                             in solver (TdPA)
!
!     07/06/2017:   V1.0.15 - Bugfix: In emergencwI_serial, the tau
!                             was calculated in the inverse direction,
!                             that is, making the lower boundary 0 and
!                             counting outwards (TdPA)
!
!     07/05/2017:   V1.0.14 - tauM needs a shift in index because
!                             it is a pointer (TdPA)
!
!     06/28/2017:   V1.0.13 - Receiving and passing Red (TdPA)
!
!     06/22/2017:   V1.0.12 - Added JKQgen and JKQgen_serial (TdPA)
!                           - Added request11 for dataO (TdPA)
!                           - Added J00P to writesolI (TdPA)
!
!     06/20/2017:   V1.0.11 - Changed etaIM, etaIO, and tauM in
!                             emergence to pointers too (TdPA)
!                           - Changed many MPI%ifx to ifx (TdPA)
!                           - Corrected bug in the indexing of
!                             transitions when the Lambda operator
!                             if multiplied by rLine (TdPA)
!
!     06/19/2017:   V1.0.10 - Changed data1, data2, rLine and rPhot
!                             arrays into two pointers each (three
!                             for data1), this avoids copies in
!                             exchange of more allocation/
!                             deallocation (TdPA)
!
!     06/16/2017:    V1.0.9 - Made the modifications consistent with
!                             the change in the Frec limits (TdPA)
!
!     06/13/2017:    V1.0.8 - Initialize rLine and rPhot (TdPA)
!                           - Added indexes for Fint and Jcalc calls,
!                             consistent with the changes in those
!                             routines (TdPA)
!
!     06/12/2017:    V1.0.7 - The frequency limits are passed to the
!                             RTcoeffi routines (TdPA)
!                           - Removed cont%l
!
!     06/09/2017:    V1.0.6 - Added solverJ and solverJ_serial (TdPA)
!                           - First iteration always CRD (TdPA)
!                           - Introduced partial redistribution
!                             radiation field iterations (TdPA)
!
!     05/13/2017:    V1.0.5 - Bugfix: Emergence, there was a 'o',
!                             instead of a 'p' in one of the calls of
!                             RTcoeffIe (TdPA)
!
!     05/12/2017:    V1.0.4 - If doing PRD, first CRD converges and
!                             then it switches off ALI (TdPA)
!                           - Bugfix: The iteration limits now use
!                             the correct variable from the Input
!                             structure (TdPA)
!
!     05/05/2017:    V1.0.3 - RTStep now knows which point is doing
!                             during the call (TdPA)
!
!     05/04/2017:    V1.0.2 - There was no wait after sending the
!                             radiative data to the last CPU, the
!                             master had chances of zeroing the
!                             buffer (TdPA)
!                           - Avoiding the use of status arguments
!                             in MPI calls (TdPA)
!
!     05/01/2017:    V1.0.1 - Revamp of the tau calculation for the
!                             height of tau=1 and the contribution
!                             function so not to rely on the master
!                             knowing about the continuum (TdPA)
!                           - Bugfix: Wrong receiving and sending
!                             index for data1 for both solver and
!                             emergence (TdPA)
!
!     04/20/2017:    V1.0.0 - Started coding (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    The dimensions are ready for blends, but they are not taken
!  into account.
!
!    Emergence now waits for every process to finish each direction
!  before letting any of them to start the next one. This is very
!  inefficient for the domain decomposition (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  solveI:
!    Manage what solver must be called depending on the MPI
!  configuration
!
!  solverI:
!    This subroutine creates and solves the radiation transfer problem
!  of the first kind
!
!  solverI_alt:
!    This subroutine creates and solves the radiation transfer problem
!  of the first kind with the alternative and slower MPI scheme
!
!  solverI_serial:
!    This subroutine creates and solves the radiation transfer problem
!  of the first kind, with just one processor
!
!  emergentI:
!    Manage what emergence must be called depending on the MPI
!  configuration
!
!  emergenceI:
!    This subroutine solves the radiation transfer equation for the
!  specified directions given the rhoKQ and JKQ
!
!  emergenceI_serial:
!    This subroutine solves the radiation transfer equation for the
!  specified directions given the rhoKQ and JKQ, with just one
!  processor
!
!  solveJ:
!    Manage what solver must be called depending on the MPI
!  configuration
!
!  solverJ:
!    This subroutine creates and solves the radiation transfer problem
!  of the first kind with only continuum
!
!  solverJ_alt:
!    This subroutine creates and solves the radiation transfer problem
!  of the first kind with only continuum with the alternative and
!  slower MPI scheme
!
!  solverJ_serial:
!    This subroutine creates and solves the radiation transfer problem
!  of the first kind with only continuum, with just one processor
!
!  JKQgenerate:
!    Manage what solver must be called depending on the MPI
!  configuration
!
!  JKQgen:
!    Corrects the JKQ for the change from multi-level to multi-term
!
!  JKQgen_alt:
!    Corrects the JKQ for the change from multi-level to multi-term
!  with the alternative and slower MPI scheme
!
!  JKQgen_serial:
!    Corrects the JKQ for the change from multi-level to multi-term
!  with just one processor
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use boundary_mod
      use commons_mod
      use fieldb_mod
      use iosolution_mod
      use jcalci_mod
      use mrc_mod
      use ng_mod
      use omp_mod
      use parameters_mod , only : kb, cSaha, fktoJ
      use rtcoeffi_mod
      use rtstepi_mod
      use seei_mod
      use setmpi_mod
      use types_mod

      ! Maximum buffer for NG_int
      double precision, parameter:: maxbuffer_NG = 500d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the solver for the NLTE problem for intensity\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!            lload(logical): Bool that says if previous
      !!                            solution was read\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine solveI(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                        Geom,MPID,Input,lload,Stokes, &
                        J00,J00S,J00C,J00P)
      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      logical:: lload
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00
      double precision, dimension(nxt,Rz0:Rz1):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P

#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,-2)
#endif
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P,Input%folder,-2)
#endif
      ! If MPI
      if (MPID%mpi) then

        if (MPID%alternI) then
          call solverI_alt(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                           Red,Geom,MPID,Input,lload,Stokes, &
                           J00,J00S,J00C,J00P)
        else
          call solverI(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                       Red,Geom,MPID,Input,lload,Stokes, &
                       J00,J00S,J00C,J00P)
        end if

      ! Serial
      else

        call solverI_serial(Atom,LTElines,Rho_old,Atmo,Cont, &
                            Frec,Red,Geom,MPID,Input,lload, &
                            Stokes,J00,J00S,J00C,J00P)
      end if ! MPI
#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,-1)
#endif
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P,Input%folder,-1)
#endif

      end subroutine solveI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem for intensity using the ALI method,
      !! with several CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!            lload(logical): Bool that says if previous
      !!                            solution was read\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine solverI(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                         Geom,MPID,Input,lload,Stokes, &
                         J00,J00S,J00C,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      logical:: lload
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00
      double precision, dimension(nxt,Rz0:Rz1):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P

      ! Local


      type(MRC_class):: MRC

      logical:: goout,gooutprd,AD,ADT,ADD,PRDl,doNG,laux,lALI
      logical:: RIRAM,deal,force_ALI

      character(LEN=20):: iterS

      integer:: iaux,ios,iz0,iz1,diz,m,o,p,op
      integer:: ith,iph,ia,iz,ifreq,if0,if1,iil,iip
      integer:: iproc,iter,itran,jtran,ftran,jftran,fftran,id,iterr
      integer:: if0l,if1l,ip0l,ip1l,nfl,nftl,nfpl,istep
      integer:: NG_dim,NG_entry,iterm,iJ,ing
      integer:: ntpz,npz,ntp,itpz,tid
#ifdef _OPENMP
      integer:: itz
#endif

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset

      double precision:: daux,WA,mu_inv,dsm,dsp
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, &
             dimension(2,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: Norm
      double precision, &
             dimension(2,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: BStk
      double precision, &
             dimension(nxb,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: BLam
      double precision, dimension(nxb,nxt,Rz0:Rz1):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1):: LambdaP
      double precision, dimension(nfreq,Rz0:Rz1):: J00Cold

      double precision, dimension(:), allocatable, target:: LO

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Stokes_r
      double precision, dimension(:), allocatable, target:: LambdaL_r
      double precision, dimension(:), allocatable, target:: LambdaP_r
      double precision, dimension(:), allocatable, target:: Prof_r
      ! Senders
      double precision, dimension(:,:,:,:), allocatable:: Stokes_s
      double precision, dimension(:,:,:,:,:), allocatable:: Prof_s
      double precision, dimension(:,:,:,:), allocatable:: rLine_s
      double precision, dimension(:,:,:,:), allocatable:: rPhot_s
      ! Dual
      integer:: info_b


      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:,:), pointer:: data2O,data2P
      double precision, dimension(:), pointer:: rLineO,rLineP
      double precision, dimension(:), pointer:: rPhotO,rPhotP
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP, p_LO
      double precision, dimension(:), pointer:: p_exu
      double precision, dimension(:,:), pointer:: p_MStk
      double precision, dimension(:,:), pointer:: p_MrLine,p_MrPhot
      double precision, dimension(:,:,:), pointer:: p_MProf


      ! Initialize converged flag
      goout = .False.

      ! Initialize force ALI
      force_ALI = .False.

      ! Initialize angle depended flag
      AD = .not.AVI

      ! Angle dependent or dynamic
      ADD = AD.or.dyn

      ! Really go back to angle dependent
      ADT = .not.AV.and.AD

      ! If storing redistribution
      RIRAM = IRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AVI = .True.
        if (ADT) IRAM = .False.
      end if

      ! Store if it was PRD (first iteration always CRD unless we are
      ! reading)
      PRDl = PRD
      if (.not.lload) PRD = .False.


      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.

      ! If NG acceleration of rho00
      if (Input%NGI) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + RnZ*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRDl) then

          ! If we need Stokes
          if (ADD) then

            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            NG_dim = NG_dim + nfreq*RnZ

          end if ! ADD
        end if ! PRD

        ! If it requires too much buffer
        if (dble(NG_dim)*8d-6.gt.maxbuffer_NG) then

          Input%NGI = .False.

          if (pid.eq.0) then
            umsg = ' # The buffer for NG acceleration '// &
                   'is too big. Not doing NG.'
            call verbose
          end if

        end if ! Buffer

        ! If finally doing it, allocate
        if (Input%NGI) then

          ! Master
          if (pid.eq.0) then

            allocate(NG_scratch(NG_dim, Input%NGI_ord+2))

          ! Slaves
          else

            allocate(NG_scratch(NG_dim, 1))

          end if

        end if ! NG and master
      end if ! NG

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)

      ! Master
      if(pid.eq.0)then

        ! Dimensions
        npz = Geom%nph*Rnz
        ntp = Geom%nth*Geom%nph
        ntpz = Geom%nth*npz

        ! To receive Intensity chunks
        iaux = MPID%nxfreq*ntpz
        allocate(Stokes_r(iaux))

        ! To receive Lambda operator for b-b transitions
        iaux = MPID%nxtfreqi*nxb*ntpz
        if (iaux.lt.1) iaux = 1
        allocate(LambdaL_r(iaux))

        ! To receive Lambda operator for b-f transitions
        iaux = MPID%nxpfreq*nxb*ntpz
        if (iaux.lt.1) iaux = 1
        allocate(LambdaP_r(iaux))

        ! To receive profile information
        iaux = MPID%nxtfreqi*2*ntpz
        if (iaux.lt.1) iaux = 1
        allocate(Prof_r(iaux))

      ! Slave
      else


        !
        ! Allocation of buffers
        !

        ! Common (Master and slave)
        ! Allocate O pointers
        allocate(data2O(Frec%ntfreqi,2))
        allocate(rLineO(Frec%ntfreqi))
        allocate(rPhotO(Frec%npfreq))

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(if0:if1,0:2))
        allocate(data1O(if0:if1,0:2))

        ! Allocate vector for Lambda operator
        allocate(LO(if0:if1))

        ! To send Intensity chunks
        allocate(Stokes_s(if0:if1,Rz0:Rz1,Geom%nPh,Geom%nTh))

        ! To send profile information
        allocate(Prof_s(Frec%ntfreqi,2,Rz0:Rz1,Geom%nPh,Geom%nTh))

        ! To send Lambda operator for b-b transitions
        allocate(rLine_s(Frec%ntfreqi,Rz0:Rz1,Geom%nPh,Geom%nTh))

        ! To send Lambda operator for b-f transitions
        allocate(rPhot_s(Frec%npfreq,Rz0:Rz1,Geom%nPh,Geom%nTh))

      end if ! Master or Slave


      !
      ! Initialization messages
      !

      ! Global Master
      if (gpid.eq.0) then

        ! Announce we are starting
        umsg = '    Iteration          MRC(rho^0_0) Atom_index '// &
               'Level_index Height_index Height(km)'
        call verbose
      end if

      ! Master storing MRC
      if (gpid.eq.0.and.Input%keep_MRC) then

        ! Open the file to store MRC

        ! If appending
        if (Input%appendMRC) then
          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRCI', exist=laux)
          if (.not.laux) then
            open(800, file=trim(Input%folder)//'/MRCI', &
                 action='write',iostat=ios,err=1000)
            write(800,'(A)') '!   Iteration          MRC(rho^0_0) '// &
                         'Atom_index Level_index Height_index '// &
                         'Height(km)'
            close(800)
          end if
        ! If not appending, create
        else
          open(800, file=trim(Input%folder)//'/MRCI', &
               action='write',iostat=ios,err=1000)
          write(800,'(A)') '!   Iteration          MRC(rho^0_0) '// &
                       'Atom_index Level_index Height_index '// &
                       'Height(km)'
          close(800)
        end if
      end if

      ! Control
      call control
      if (laborted) goto 2000

      ! Initialize stimulated radiation tensor, becuase it is
      ! commented in the main loop
      J00S = 0d0

      ! If measuring performance
      if (Input%mpi_perf.and.pid.eq.0) &
        call report_mpi_timeI(Input%folder,Input%ID, &
                              0,0,0,.False.)


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=Input%iteri_min,Input%iteri_max

        ! Flags for physics
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics
        if (iter.le.Input%allownphys_pop) then
          if (.not.nphysR) nphysR = .True.
        else
          if (nphysR) nphysR = .False.
        end if


        !
        ! The old population will be the current one
        !

        ! For each atom
        do ia=1,nA
          Rho_old(ia)%crho = Atom(ia)%crho
        end do

        ! Flag for ALI
        lALI = iter.gt.Input%ALI_delay.or.force_ALI

        ! Internal PRD iterations
        do iterr=1,Input%iteri_prd


          !
          ! Master
          !
          if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(id,deal,WA,itpz,itz,ith,iph,iz,p_exu,jtran,jftran,ia) &
!$omp private(daux,dsm) &
!$omp shared(J00,J00P,J00C,BStk,BLam,lALI,LambdaL,LambdaP,Norm) &
!$omp shared(info_b,Stokes_r,laborted,Prof_r,LambdaL_r,LambdaP_r) &
!$omp shared(iter,iterr,p_MStk,p_MProf,p_MrLine,p_MrPhot,PIRAM) &
!$omp shared(ntpz,npz,na,nz,nxb,if0l,if1l,ip0l,ip1l,nfl,nftl,nfpl) &
!$omp shared(MPID,Input,Frec,Geom,Atom,Atmo,Stokes,KSTK,J00Cold,PRD) &
!$omp shared(Rz0,Rz1,Rnz,ierr) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT)

            ! Save old radiation field
            if (PRD) then
!$omp workshare
              J00Cold = J00C
!$omp end workshare
            end if
            ! Reset radiation field variables
!$omp workshare
            J00 = 0d0
           !J00S = 0d0
            J00P = 0d0
            J00C = 0d0
            BStk = 0d0
            BLam = 0d0
            Norm = 0d0
!$omp end workshare
            if (lALI) then
!$omp workshare
              LambdaL = 0d0
              LambdaP = 0d0
!$omp end workshare
            end if

            ! Each frequency cut
            do id=1,MPID%nnd

!$omp single
              !
              ! Receive data from a slave
              !

              ! Receive indexing data
              do while (.True.)
                call MPI_recv(info_b,1,MPI_INTEGER, &
                              MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do

              ! Flag error
              if (info_b.lt.0) laborted = .True.
!$omp end single

              ! Continue?
              if (info_b.lt.0) cycle

!$omp single
              ! Receive intensity
              do while (.True.)
                call MPI_recv(Stokes_r(1), MPID%sizei4(info_b), &
                              MPI_DOUBLE_PRECISION, info_b, &
                              info_b, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do

              ! Receive profile
              do while (.True.)
                call MPI_recv(Prof_r(1), MPID%sizei5(info_b), &
                              MPI_DOUBLE_PRECISION, info_b, &
                              1000000+info_b, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do

              ! Only if ALI iteration
              if (lALI) then

                ! Receive Lambda operator for b-b transition
                do while (.True.)
                  call MPI_recv(LambdaL_r(1), MPID%sizei9(info_b), &
                                MPI_DOUBLE_PRECISION, info_b, &
                                2000000+info_b, MPI_COMM_RT, &
                                MPI_STATUS_IGNORE, ierr)
                  if (ierr.eq.0) exit
                end do

                ! Receive Lambda operator for b-f transition
                do while (.True.)
                  call MPI_recv(LambdaP_r(1), MPID%sizei0(info_b), &
                                MPI_DOUBLE_PRECISION, info_b, &
                                3000000+info_b, MPI_COMM_RT, &
                                MPI_STATUS_IGNORE, ierr)
                  if (ierr.eq.0) exit
                end do

              end if ! ALI iteration

              ! If measuring performance
              if (Input%mpi_perf) &
                call report_mpi_timeI(Input%folder,Input%ID, &
                                      info_b,iter,iterr,.True.)

              ! Shorter variables
              if0l = MPID%if0(info_b)
              if1l = MPID%if1(info_b)
              ip0l = Frec%Mpif0(info_b)
              ip1l = Frec%Mpif1(info_b)
              nfl = MPID%nf(info_b)
              nftl = Frec%Mntfreqi(info_b)
              nfpl = Frec%Mnpfreq(info_b)

              ! Pointers
              p_MStk(if0l:if1l,1:ntpz) => &
                                      Stokes_r(1:MPID%sizei4(info_b))
              p_MProf(1:nftl,1:2,1:ntpz) => &
                                        Prof_r(1:MPID%sizei5(info_b))
              ! Point to actual Lambda data
              if (lALI) then
                p_MrLine(1:nftl,1:ntpz) => &
                                     LambdaL_r(1:MPID%sizei9(info_b))
                p_MrPhot(1:nfpl,1:ntpz) => &
                                     LambdaP_r(1:MPID%sizei0(info_b))
              else
                ! Point to whatever
                p_MrLine(1:1,1:ntpz) => Prof_r(1:ntpz)
                p_MrPhot(1:1,1:ntpz) => Prof_r(1:ntpz)
              end if

!$omp end single

              deal = .False.

              ! If we are aborting, just skip this
              if (laborted) cycle

              ! Compute line quantities
!$omp do
              do itpz=1,ntpz

                ! Get indexes
                ith = (itpz-1)/npz
                iph = (itpz - npz*ith - 1)/Rnz
                iz = itpz - Rnz*iph - npz*ith + Rz0 - 1
                ith = ith + 1
                iph = iph + 1

                ! Determine where to store intensity
                if (KSTK.or.iz.eq.Rz0) &
                  Stokes(if0l:if1l,iph,ith,iz) = &
                                                p_MStk(if0l:if1l,itpz)

                ! Calculate frequency integral for b-b quantities
                call FIntI_line(Atom,MPID,Frec%W_freq,info_b, &
                                p_MStk(:,itpz),p_MrLine(:,itpz), &
                                p_MProf(:,:,itpz), &
                                Norm(:,:,iph,ith,iz), &
                                BStk(:,:,iph,ith,iz), &
                                BLam(:,:,iph,ith,iz), &
                                lALI)
#ifdef _OPENMP
              end do
!$omp end do

              !
              ! Compute continuum and photoionizations
              !

              ! For each height
!$omp do
              do iz=Rz0,Rz1
#endif
                ! Point to exu values
                if (PIRAM.and.ip1l.ge.ip0l) then
                  p_exu => Frec%exu(ip0l:ip1l,iz)
                else
                  allocate(p_exu(1))
                  deal = .True.
                end if
#ifdef _OPENMP
                ! For each polar direction
                do ith=1,Geom%nth

                  ! Partial indexing
                  itz = iz + (ith-1)*npz - Rz0 + 1

                  ! For each azimuth
                  do iph=1,Geom%nph

                    ! Get running index
                    itpz = itz + Rnz*(iph-1)
#endif
                    ! Get angular weight
                    WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                    !
                    ! Calculate rest of integrals
                    !
                    call FIntI_rest(Atom,MPID,Frec%omega, &
                                    Frec%W_freq,ip0l,ip1l, &
                                    Atmo%T(iz),info_b,WA, &
                                    p_MStk(:,itpz),p_MrPhot(:,itpz), &
                                    J00P(:,:,iz),J00C(:,iz), &
                                    LambdaP(:,:,:,iz),lALI, &
                                    p_exu)

#ifdef _OPENMP
                  end do
                end do
#endif
                ! Nullify pointer
                if (deal) deallocate(p_exu)
                nullify(p_exu)

              end do
!$omp end do
            end do ! frequency domains

!$omp single
            ! Nullify pointers
            nullify(p_MStk,p_MProf,p_MrLine,p_MrPhot)
!$omp end single

            !
            ! Apply weights to J00, J00S, and Lambda operator
            ! and normalize
            !


            ! For each height
!$omp do
            do iz=Rz0,Rz1

              ! For each polar direction
              do ith=1,Geom%nTh

                ! For each azimuthal direction
                do iph=1,Geom%nph

                  ! Get the angular integral weight
                  WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                  ! For each atom
                  do ia=1,nA

                    ! For each FS transition
                    do ftran=1,Atom(ia)%nftran

                      ! Apply shift
                      jftran = ftran + Atom(ia)%tfshift

                      ! Get the weight
                      if (Norm(1,jftran,iph,ith,iz).gt.0d0) then

                        daux = WA/Norm(1,jftran,iph,ith,iz)

                        ! Integrate angle
                        J00(jftran,iz) = J00(jftran,iz) + &
                                         BStk(1,jftran,iph,ith,iz)* &
                                         daux

                        ! For each transition blended with ftran
                        if (lALI) then

                          do fftran=1,nxb
                            LambdaL(fftran,jftran,iz) = &
                                     LambdaL(fftran,jftran,iz) + &
                                     BLam(fftran,jftran,iph,ith,iz)* &
                                     daux
                          end do

                        end if !LI

                      end if

                     !! If there is stimulated emission
                     !if (stm) then

                     !  ! Get the weight
                     !  if (Norm(2,jftran,iph,ith,iz).gt.0d0) then

                     !    daux = WA/Norm(2,jftran,iph,ith,iz)

                     !    ! Integrate angle
                     !    J00S(jftran,iz) = J00S(jftran,iz) + &
                     !                    BStk(2,jftran,iph,ith,iz)* &
                     !                    daux

                     !  end if

                     !end if ! Stimulated emission

                    end do ! FS transition
                  end do ! Atoms
                end do ! azimuthal directions
              end do ! polar directions

              !
              ! Add the Saha factor (ne*Zeta) to J00P
              !

              ! Argument of the exponential
              WA = fktoJ/kb/Atmo%T(iz)

              ! Part that does not depend on the line
              daux = cSaha*Atmo%ne(iz)/(Atmo%T(iz)**(1.5d0))

              ! For each atom
              do ia=1,nA

                ! For each b-f transition
                do itran=1,Atom(ia)%nphot

                  ! Apply shift
                  jtran = itran + Atom(ia)%pshift

                  ! Calculate the multiplicative factor
                  dsm = daux*exp(Atom(ia)%phot(itran)%edge*WA)* &
                        Atom(ia)%phot(itran)%glu

                  ! Apply it to the emission integral
                  J00P(jtran,2,iz) = J00P(jtran,2,iz)*dsm

                  ! Apply it to the emission Lambda operator
                  if (lALI) &
                    LambdaP(:,jtran,2,iz) = LambdaP(:,jtran,2,iz)*dsm

                end do ! b-f transitions
              end do ! atoms
            end do ! heights
!$omp end do
!$omp end parallel

            ! Calculate MRC for J if PRD
            if (PRD.and.Input%iteri_prd.gt.1) then

              ! Call the routine
              call MRCJ_sb(J00C,J00Cold,MRC)

              ! Convert cm into km
              MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

              ! Only global Máster writes
              if (gpid.eq.0) then

                if (iterr.eq.1.and.iter.eq.2) then
                  umsg = '         PRD            MRC(J^0_0) '// &
                         'Freq_index  Wavelength Height_index '// &
                         'Height(km)'
                  call verbose
                end if

                ! Write in stdout
                write(umsg,'(2x,"PRD it:",1x,i3,2x,es20.12,'// &
                           '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                           iterr,MRC%values(2,1),MRC%indexes(1,1), &
                           1d2/Frec%omega(MRC%indexes(1,1)), &
                           MRC%indexes(2,1),MRC%values(1,1)
                call verbose

              end if

              ! Check exit criteria
              if (MRC%values(2,1).le.Input%mrci_r.or. &
                  iterr.eq.Input%iteri_prd) then
                gooutprd = .True.
              else
                gooutprd = .False.
              end if

            ! Always go out if no PRD
            else

              gooutprd = .True.

            end if

            ! If aborting, obviously go out
            if (laborted) gooutprd = .True.

          !
          ! Slave
          !
          else

            !
            ! Ratiation Transfer
            !

            !  For each polar direction
            do ith=1,Geom%nTh

              ! Calculate inverse of cosine of polar direction
              mu_inv = 1d0/Geom%V_mu(ith)

              ! Determine the direction of propagation for indexes
              diz = -int(sign(1d0, Geom%V_mu(ith)))

              !
              ! Determine the boundaries and the domain decomposition
              ! partnerts
              !

              ! Determine the first and last height indexes to run
              ! over
              iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
              iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

              ! For each azimuthal direction
              do iph=1,Geom%nPh

                !
                ! If angle-dependent, manage scattering angles
                if (.not.AVI.and.PRD) &
                  call scattering_manage(Geom,ith,iph)


                !
                ! First height
                !

                ! If going down, get top boundary
                if(diz.eq.1)then

                  ! Call top boundary
                  call topI(MPID,data1M(:,2))

                ! If going up, get bottom boundary
                else

                  ! Call bottom boundary
                  call bottomI(Frec%omega,Atmo%T(iz0), &
                               Atmo%vx(iz0),Atmo%vy(iz0), &
                               Atmo%vz(iz0),Geom%V_mu(ith), &
                               Geom%V_mux(iph),Geom%V_muy(iph), &
                               MPID,data1M(:,2))

                endif ! propagation direction

                ! Identify current height
                o = iz0

                ! Index for Stokes
                if (PRDl.and.ADD) op = o

                ! Calculate radiative coefficients
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              o,ith,iph,if0,if1,J00(:,o), &
                              J00C(:,o),Cont%ndir, &
                              Cont%c(:,:,:,o),Stokes(:,:,:,op), &
                              rLineO(:),rPhotO(:), &
                              data1M(:,0:1),data2O(:,:),.True.)
                if (laborted) goto 3000

                !
                ! Store in buffer
                !

                ! Intensity
                Stokes_s(:,o,iph,ith) = data1M(:,2)

                ! Profiles
                Prof_s(:,:,o,iph,ith) = data2O(:,:)

                ! If ALI iteration
                if (lALI) then

                  ! b-b Lambda operator (bottom boundary does not
                  ! contribute)
                  rLine_s(:,o,iph,ith) = 0d0

                  ! b-f Lambda operator (bottom boundary does not
                  ! contribute)
                  rPhot_s(:,o,iph,ith) = 0d0

                end if

                ! Identify next height
                p = iz0 + diz

                ! Index for Stokes
                if (PRDl.and.ADD) op = p

                ! Calculate radiative coefficients
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              p,ith,iph,if0,if1,J00(:,p), &
                              J00C(:,p),Cont%ndir, &
                              Cont%c(:,:,:,p),Stokes(:,:,:,op), &
                              rLineO(:),rPhotO(:), &
                              data1O(:,0:1),data2O(:,:),.True.)
                if (laborted) goto 3000


                !
                ! Intermedium heights
                !

                ! For each height this CPU has assigned
                do iz=iz0,iz1,diz

                  ! We treat the boundaries outside
                  if(iz.eq.iz0.or.iz.eq.iz1)cycle

                  ! Allocate P pointers
                  allocate(data1P(if0:if1,0:2))
                  allocate(data2P(Frec%ntfreqi,2))
                  allocate(rLineP(Frec%ntfreqi))
                  allocate(rPhotP(Frec%npfreq))

                  ! Identify heights
                  m = iz - diz
                  o = iz
                  p = iz + diz

                  ! Calculate distance to previous point
                  dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                  ! Calculate distance to the next point
                  dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

                  ! If tau scale
                  if (ztau) then
                    dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                        Atmo%chi500(m))
                    dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                        Atmo%chi500(p))
                  end if

                  ! Index for Stokes
                  if (PRDl.and.ADD) op = p

                  ! Calculate radiative coefficients
                  call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID, &
                                Geom,p,ith,iph,if0,if1,J00(:,p), &
                                J00C(:,p),Cont%ndir, &
                                Cont%c(:,:,:,p),Stokes(:,:,:,op), &
                                rLineP(:),rPhotP(:), &
                                data1P(:,0:1),data2P(:,:),.True.)

                  ! Point to the data
                  p_K0M  => data1M(:,0)
                  p_SM   => data1M(:,1)
                  p_StkM => data1M(:,2)
                  p_K0O  => data1O(:,0)
                  p_SO   => data1O(:,1)
                  p_StkO => data1O(:,2)
                  p_K0P  => data1P(:,0)
                  p_SP   => data1P(:,1)
                  p_LO => LO(:)

                  ! Apply short characteristics BESSER
                  call RTStepI(o,ith,iph,MPID%nf(pid), &
                               dsm,dsp,p_K0M,p_SM,p_K0O, &
                               p_SO,p_K0P,p_SP,p_StkM, &
                               p_StkO,p_LO,lALI,.True.)


                  !
                  ! Combine the value of lambda operator with the
                  ! transition strength
                  !
                  if (lALI) then

                    ! Initialize indexes for rLine and rPhot
                    iil = 0
                    iip = 0

                    ! For each atom
                    do ia=1,nA

                      ! For each b-b trantision
                      do itran=1,Atom(ia)%ntran

                        ! If this CPU does not have frequencies in
                        ! this line, skip
                        if (Atom(ia)%fflag(itran)%absent) cycle

                        ! For each FS transition
                        do ftran=1,Atom(ia)%fst(itran)%nt

                          ! Get the sequential index of FS transition
                          fftran = Atom(ia)%ifst_ij(ftran,itran)

                          ! Apply shift
                          jftran = fftran + Atom(ia)%tfshift

                          ! For each frequency
                          do ifreq=Atom(ia)%if0(itran), &
                                   Atom(ia)%if1(itran)

                            iil = iil + 1
                            rLineO(iil) = LO(ifreq)*rLineO(iil)

                          end do ! frequency
                        end do ! FS transition
                      end do ! b-b transition

                      ! For each b-f transition
                      do itran=1,Atom(ia)%nphot

                        ! If this CPU does not have frequencies in
                        ! this transition, skip
                        if (Atom(ia)%phot(itran)%absent) cycle

                        ! Apply shift
                        jtran = itran + Atom(ia)%pshift

                        ! For each frequency
                        do ifreq=Atom(ia)%phot(itran)%if0, &
                                 Atom(ia)%phot(itran)%if1

                          iip = iip + 1
                          rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                        end do ! frequency
                      end do ! b-f transition
                    end do ! atom

                    ! b-b Lambda operator
                    rLine_s(:,o,iph,ith) = rLineO(:)

                    ! Send b-f Lambda operator
                    rPhot_s(:,o,iph,ith) = rPhotO(:)

                  end if ! ALI


                  !
                  ! Store in buffer
                  !

                  ! Intensity
                  Stokes_s(:,o,iph,ith) = data1O(:,2)

                  ! Profiles
                  Prof_s(:,:,o,iph,ith) = data2O(:,:)


                  ! Shift data (O->M, P->O)
                  deallocate(data1M,data2O)
                  data1M => data1O
                  data1O => data1P
                  data2O => data2P
                  nullify(data1P,data2P)
                  if (lALI) then
                    deallocate(rLineO,rPhotO)
                    rLineO => rLineP
                    rPhotO => rPhotP
                  else
                    deallocate(rLineP,rPhotP)
                  end if
                  nullify(rLineP,rPhotP)

                  ! Error
                  if (laborted) goto 3000

                end do ! Intermedium heights


                !
                ! Last height
                !

                ! Identify heights
                m = iz1 - diz
                o = iz1

                ! Calculate distance to previous point
                dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                ! If tau scale
                if (ztau) &
                  dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(m))

                ! Point to the data
                p_K0M  => data1M(:,0)
                p_SM   => data1M(:,1)
                p_StkM => data1M(:,2)
                p_K0O  => data1O(:,0)
                p_SO   => data1O(:,1)
                p_StkO => data1O(:,2)
                p_LO => LO(:)

                ! Apply short characteristics LINEAR
                call RTStepI(o,ith,iph,MPID%nf(pid), &
                             dsm,dsp,p_K0M,p_SM,p_K0O, &
                             p_SO,p_K0P,p_SP,p_StkM, &
                             p_StkO,p_LO,lALI,.False.)
                if (laborted) goto 3000


                !
                ! Combine the value of lambda operator with the
                ! transition strength
                !

                if (lALI) then

                  ! Initialize indexes for rLine and rPhot
                  iil = 0
                  iip = 0

                  ! For each atom
                  do ia=1,nA

                    ! For each b-b trantision
                    do itran=1,Atom(ia)%ntran

                      ! If this CPU does not have frequencies in
                      ! this line, skip
                      if (Atom(ia)%fflag(itran)%absent) cycle

                      ! For each FS transition
                      do ftran=1,Atom(ia)%fst(itran)%nt

                        ! Get the sequential index of FS transition
                        fftran = Atom(ia)%ifst_ij(ftran,itran)

                        ! Apply atomic shift
                        jftran = fftran + Atom(ia)%tfshift

                        ! For each frequency
                        do ifreq=Atom(ia)%if0(itran), &
                                 Atom(ia)%if1(itran)

                          iil = iil + 1
                          rLineO(iil) = LO(ifreq)*rLineO(iil)

                        end do ! frequency
                      end do ! FS transition
                    end do ! b-b transition

                    ! For each b-f transition
                    do itran=1,Atom(ia)%nphot

                      ! If this CPU does not have frequencies in
                      ! this transition, skip
                      if (Atom(ia)%phot(itran)%absent) cycle

                      ! Apply atomic shift
                      jtran = itran + Atom(ia)%pshift

                      ! For each frequency
                      do ifreq=Atom(ia)%phot(itran)%if0, &
                               Atom(ia)%phot(itran)%if1

                        iip = iip + 1
                        rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                      end do ! frequency
                    end do ! b-f transition
                  end do ! atom

                  ! b-b Lambda operator
                  rLine_s(:,o,iph,ith) = rLineO(:)

                  ! b-f Lambda operator
                  rPhot_s(:,o,iph,ith) = rPhotO(:)

                end if ! ALI


                !
                ! Store in buffer
                !

                ! Intensity
                Stokes_s(:,o,iph,ith) = data1O(:,2)

                ! Profiles
                Prof_s(:,:,o,iph,ith) = data2O(:,:)

              enddo ! azimuthal angles
            enddo ! polar angles


            !
            ! Send to master
            !

            ! If had an error
3000        if (laborted) then

              ! Send error
              do while (.True.)
                call MPI_SEND(-pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                              ierr)
                if (ierr.eq.0) exit
              end do

            ! No problem
            else

              ! Send indexes
              do while (.True.)
                call MPI_SEND(pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                            ierr)
                if (ierr.eq.0) exit
              end do

              ! Send Stokes
              do while (.True.)
                call MPI_SEND(Stokes_s(if0,Rz0,1,1), &
                              MPID%sizei4(pid), &
                              MPI_DOUBLE_PRECISION, 0, pid, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

              ! Send profiles
              do while (.True.)
                call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                              MPID%sizei5(pid), &
                              MPI_DOUBLE_PRECISION, &
                              0, 1000000+pid, MPI_COMM_RT, &
                              ierr)
                if (ierr.eq.0) exit
              end do

              ! If ALI
              if (lALI) then

                ! Send b-b Lambda operator
                do while (.True.)
                  call MPI_SEND(rLine_s(1,Rz0,1,1), &
                                MPID%sizei9(pid), &
                                MPI_DOUBLE_PRECISION, &
                                0, 2000000+pid, MPI_COMM_RT, &
                                ierr)
                  if (ierr.eq.0) exit
                end do

                ! Send b-f Lambda operator
                do while (.True.)
                  call MPI_SEND(rPhot_s(1,Rz0,1,1), &
                                MPID%sizei0(pid), &
                                MPI_DOUBLE_PRECISION, &
                                0, 3000000+pid, MPI_COMM_RT, &
                                ierr)
                  if (ierr.eq.0) exit
                end do

              end if ! ALI
            end if ! Had a problem
          end if ! Master or Slave

          !
          ! Share if we are finished or not (in PRD)
          !

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive goout
              call MPI_RECV(gooutprd, 1, MPI_LOGICAL, &
                            MPID%recv, 7000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No Master

            ! For each send
            do istep=1,MPID%nsend

              ! Send goout
              call MPI_ISEND(gooutprd, 1, MPI_LOGICAL, &
                             MPID%lsend(istep), &
                             7000000+MPID%lsend(istep), &
                             MPI_COMM_RT, MPID%requestA(istep,1), &
                             ierr)
            end do ! sends

          ! Normal bcast
          else

            call MPI_BCAST(gooutprd, 1, MPI_LOGICAL, 0, &
                           MPI_COMM_RT, ierr)

          end if

          ! Control
          call control
          if (laborted) goto 2000


          !
          ! Share the radiation information
          !

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive J00
              call MPI_RECV(J00(1,Rz0), MPID%sizei6(0), &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

             !! Receive J00S if stimulated emission
             !if (stm) &
             !  call MPI_RECV(J00S(1,Rz0), MPID%sizei6(0), &
             !                MPI_DOUBLE_PRECISION,  &
             !                MPID%recv, 1000000+pid, &
             !                MPI_COMM_RT, MPI_STATUS_IGNORE, &
             !                ierr)

              ! Receive J00C
              call MPI_RECV(J00C(1,Rz0), MPID%sizei7(0), &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, 2000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

              ! Receive intensity if doing A-D PRD
              if (PRDl.and.ADD) &
                call MPI_RECV(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                              MPI_DOUBLE_PRECISION,  &
                              MPID%recv, 3000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              ! If we are going to do SEE
              if (gooutprd) then

                ! Receive J00 for b-f transitions
                call MPI_RECV(J00P(1,1,Rz0), MPID%sizei3(0), &
                              MPI_DOUBLE_PRECISION,  &
                              MPID%recv, 4000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

                ! ALI
                if (lALI) then

                  ! Receive Lambda operator for b-b transitions
                  call MPI_RECV(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, 5000000+pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                  ! Receive Lambda operator for b-f transitions
                  call MPI_RECV(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, 6000000+pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                end if ! ALI
              end if ! To do SEE
            end if ! No master

            ! For each send
            do istep=1,MPID%nsend

              ! Send J00
              call MPI_ISEND(J00(1,Rz0), MPID%sizei6(0), &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,1), ierr)

              ! Send J00 if stimulated emission
             !if (stm) &
             !  call MPI_ISEND(J00S(1,Rz0), MPID%sizei6(0), &
             !                 MPI_DOUBLE_PRECISION, &
             !                 MPID%lsend(istep), &
             !                 1000000+MPID%lsend(istep), &
             !                 MPI_COMM_RT, &
             !                 MPID%requestA(istep,2), ierr)

              ! Send J00C
              call MPI_ISEND(J00C(1,Rz0), MPID%sizei7(0), &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             2000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,3), ierr)

              ! Send intensity if doing A-D PRD
              if (PRDl.and.ADD) &
                call MPI_ISEND(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               3000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,4), ierr)

              ! If we are going to do SEE
              if (gooutprd) then

                ! Send J00P
                call MPI_ISEND(J00P(1,1,Rz0), MPID%sizei3(0), &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               4000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,5), ierr)

                ! ALI
                if (lALI) then

                  ! Send Lambda operator for b-b transitions
                  call MPI_ISEND(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 5000000+MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,6), ierr)

                  ! Send Lambda operator for b-f transitions
                  call MPI_ISEND(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 6000000+MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,7), ierr)
                end if ! ALI
              end if ! To do SEE

            end do ! Sends

          ! Normal bcast
          else

            ! Share J00
            call MPI_BCAST(J00(1,Rz0), MPID%sizei6(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

           !! Share J00S
           !if (stm) &
           !  call MPI_BCAST(J00S(1,Rz0), MPID%sizei6(0), &
           !                 MPI_DOUBLE_PRECISION, 0, &
           !                 MPI_COMM_RT, ierr)

            ! Share J00C
            call MPI_BCAST(J00C(1,Rz0), MPID%sizei7(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Share intensity if doing A-D PRD
            if (PRDl.and.ADD) &
              call MPI_BCAST(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                            MPI_DOUBLE_PRECISION, 0, &
                            MPI_COMM_RT, ierr)

            ! If we are going to do SEE
            if (gooutprd) then

              ! Share J00 for b-f transitions
              call MPI_BCAST(J00P(1,1,Rz0), MPID%sizei3(0), &
                            MPI_DOUBLE_PRECISION, 0, &
                            MPI_COMM_RT, ierr)

              ! ALI
              if (lALI) then

                ! Share Lambda operator for b-b transitions
                call MPI_BCAST(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                              MPI_DOUBLE_PRECISION, 0, &
                              MPI_COMM_RT, ierr)

                ! Share Lambda operator for b-f transitions
                call MPI_BCAST(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                              MPI_DOUBLE_PRECISION, 0, &
                              MPI_COMM_RT, ierr)

              end if ! ALI
            end if ! To do SEE

          end if ! Type of bcast

          ! If we are finished, everyone exits
          if (.not.gooutprd.and.MPID%altbcast) then

            ! For each slave
            do iproc=1,MPID%nsend

              ! Wait for everyone to receive the radiation data
              ! before continuing and reset the buffers
              call MPI_WAIT(MPID%requestA(iproc,1), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,2), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,3), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,4), &
                            MPI_STATUS_IGNORE,ierr)

              end do ! Processors

          end if ! If not going out

          ! If we are finished, everyone exits
          if (gooutprd) exit

        end do ! PRD iteration

        ! Control
        if (laborted) goto 2000


        !
        ! Solve SEE
        !
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P, &
                                  Input%folder,iter)
#endif
#ifdef DEBUGLAMBDA
      if (pid.eq.0) call dump_lambda(Atom,LambdaL,LambdaP, &
                                     Input%folder,iter)
#endif

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,fftran,jftran,urou,umsg,tid) &
!$omp shared(nA,nz,Atom,Rho_old,J00,J00P,LambdaL,LambdaP,lALI) &
!$omp shared(laborted,vaborted,Rz0,Rz1)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1


          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL(:,itran:jtran,iz), &
                      LambdaP(:,fftran:jftran,:,iz),iz,lALI,tid,INPUT)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL(:,itran:jtran,iz), &
                      LambdaP(:,fftran:jftran,:,iz),iz,lALI,tid)
#endif

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel
#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,iter)
#endif

        ! Control
        call control
        if (laborted) goto 2000


        !
        ! NG acceleration
        !

        ! Check if doing NG acceleration
        if(Input%NGI.and.iter.gt.Input%NGI_delay)then

          ! Advance ntry. The master does not need it if
          ! if accelerated the intensity
          NG_entry = NG_entry + 1

          ! If Master
          if (pid.eq.0) then

            ! Initialize index
            o = 0

            ! For each atom
            do ia=1,NA

              ! For each term
              do iterm=1,Atom(ia)%nMulti

                ! For each level
                do iJ=1,Atom(ia)%nJ(iterm)

                  ! Get level index
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    o = o + 1

                    ! Store rho
                    NG_scratch(o,NG_entry) = dble(Atom(ia)%crho(p,iz))

                  enddo ! Heights
                enddo ! Levels
              enddo ! Term
            enddo ! Atoms

            ! If PRD
            if (PRDl) then

              ! If ADD
              if (ADD) then

                ! For each height
                do iz=Rz0,Rz1

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Advance index
                        o = o + 1

                        ! Store Stokes
                        NG_scratch(o,NG_entry) = &
                                            Stokes(ifreq,iph,ith,iz)

                      end do ! Frequencies
                    end do ! Azimuthal directions
                  end do ! Polar directions
                end do ! Heights

              ! If AV
              else

                ! For each height
                do iz=Rz0,Rz1

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! Advance index
                    o = o + 1

                    ! Store Stokes
                    NG_scratch(o,NG_entry) = J00C(ifreq,iz)

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD or dyn
            end if ! PRD

            ! Call NG and check if it should be processed
            call NG(NG_dim,Input%NGI_ord,NG_scratch,NG_entry,doNG)

          ! Slave
          else

            ! If wrong order
            if (Input%NGI_ord.lt.1.or.Input%NGI_ord.gt.5) then
              doNG = .False.
            ! Valid order
            else
              ! Check if Master is in NG step
              if (NG_entry.gt.(Input%NGI_ord+1)) then
                doNG = .True.
              ! Not a NG step
              else
                doNG = .False.
              end if ! NG step
            end if ! order

          end if ! Master of slave

          ! If communication is needed
          if (doNG) then

            ! If Master, send last
            if (pid.eq.0) then
              ing = NG_entry
            ! Slaves send the one they have
            else
              ing = 1
            end if

            ! Alternative bcast
            if (MPID%altbcast) then

              ! If not master, receive first
              if (pid.ne.0) then

                ! Receive NG iteration
                call MPI_RECV(NG_scratch(1,1), NG_dim, &
                              MPI_DOUBLE_PRECISION,  &
                              MPID%recv, 7000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No master

              ! For each send
              do istep=1,MPID%nsend

                ! Send NG iteration
                call MPI_ISEND(NG_scratch(1,ing), NG_dim, &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               7000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,8), ierr)

              end do ! Sends

            ! Normal bcast
            else

              ! Share NG iteration
              call MPI_BCAST(NG_scratch(1,ing), NG_dim, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

            end if ! Type of bcast

            ! Reconstruct NG data

            ! Initialize index
            o = 0

            ! For each atom
            do ia=1,NA

              ! For each term
              do iterm=1,Atom(ia)%nMulti

                ! For each level
                do iJ=1,Atom(ia)%nJ(iterm)

                  ! Get level index
                  daux = sqrt(2d0*Atom(ia)%rJval(iJ,iterm) + 1d0)
                  m = Atom(ia)%irho(iterm)%irho_ij(iJ)
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Advance index
                    o = o + 1

                    ! Accelerate rho
                    Atom(ia)%crho(p,iz) = &
                                        dcmplx(NG_scratch(o,ing), 0d0)
                    Atom(ia)%popu(m,iz) = NG_scratch(o,ing)*daux

                  enddo ! Heights
                enddo ! Levels
              enddo ! Terms
            enddo ! Atoms

            ! If PRD
            if (PRDl) then

              ! If ADD
              if (ADD) then

                ! For each height
                do iz=Rz0,Rz1

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Advance index
                        o = o + 1

                        ! Store Stokes
                        Stokes(ifreq,iph,ith,iz) = NG_scratch(o,ing)

                      end do ! Frequencies
                    end do ! Azimuthal directions
                  end do ! Polar directions
                end do ! Heights

              ! If AV
              else

                ! For each height
                do iz=Rz0,Rz1

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! Advance index
                    o = o + 1

                    ! Store Stokes
                    J00C(ifreq,iz) = NG_scratch(o,ing)

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD or dyn
            end if ! PRD

            if(gpid.eq.0) then
              umsg = 'NG acceleration'
              call verbose
            end if

            NG_entry = 0

            ! Master wait to free the buffer if alternative bcast
            if (pid.eq.0.and.MPID%altbcast) then

              ! For each slave
              do iproc=1,MPID%nsend

                ! Wait for everyone to receive the radiation data
                ! before continuing and reset the buffers
                call MPI_WAIT(MPID%requestA(iproc,8), &
                              MPI_STATUS_IGNORE,ierr)

              end do ! Processors

            end if ! Master
          end if ! Communication
        endif ! NG acceleration


        !
        !  Calculate MRC
        !

        ! Only the master does
        if (pid.eq.0) then

          ! Calculate MRC
          call MRCI_sb(Atom,Rho_old,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Check exit criteria
          if (MRC%values(2,1).le.Input%mrci_i) goout = .True.

          ! Check in first iteration if was a fix population case
          ! but we want PRD
          if (iter.eq.Input%iteri_min.and.PRDl.and. &
            MRC%values(2,1).lt.1d-16) goout = .False.
          ! If going to AD from AA
          if (iter.eq.Input%iteri_min.and.tbAD.and.ADT) &
            goout = .False.

        end if

        ! We can swith now to AD if we had AV input
        if (tbAD.and.ADT) then
          AVI = .False.
          tbAD = .False.
          IRAM = RIRAM
        end if

        ! Only the global Master does
        if (gpid.eq.0) then

          ! Write in stdout
          write(umsg,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                     '2x,f9.3)') &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(2,1),MRC%values(1,1)
          call verbose

          ! File
          if (Input%keep_MRC) then

            ! Write in MRC file
            open(800, file=trim(Input%folder)//'/MRCI', &
                 iostat=ios,err=1000,position='append')
            write(800,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                      '2x,f9.3)', err=1100) &
            iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
            MRC%indexes(2,1),MRC%values(1,1)
            close(800)

          end if
        end if


        !
        ! Save partial solution
        !

        ! Check if we want to store partial results
        if(Input%storei.and.mod(iter,Input%storei_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesolI(Input,iterS,Frec%omega,Geom, &
                        !Atom,Atmo%z,Stokes,J00,J00S,J00C,J00P, &
                         Atom,Atmo%z,Stokes,J00,J00,J00C,J00P, &
                         .False.)
        ! Or have a control check
        else

          call control
          if (laborted) goto 2000

        end if

        ! If in the mandatory non-PRD
        if (iter.le.Input%PRD_delay) goout = .False.

        ! Recover the PRD variable and avoid to go out if
        ! started with CRD
        if (iter.eq.Input%PRD_delay) PRD = PRDl

        !
        ! Share if we are finished or not

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive goout
            call MPI_RECV(goout, 1, MPI_LOGICAL, &
                          MPID%recv, 7000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, ierr)

          end if ! No Master

          ! For each send
          do istep=1,MPID%nsend

            ! Send goout
            call MPI_ISEND(goout, 1, MPI_LOGICAL, &
                           MPID%lsend(istep), &
                           7000000+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,1), &
                           ierr)
          end do ! sends

        ! Normal bcast
        else

          ! Share goout
          call MPI_BCAST(goout, 1, MPI_LOGICAL, 0, &
                         MPI_COMM_RT, ierr)

        end if

        ! If alternative bcast
        if (MPID%nsend.gt.0.and.MPID%altbcast) then

          ! For each slave
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing and reset the buffers
            call MPI_WAIT(MPID%requestA(iproc,1), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,2), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,3), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,4), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,5), &
                          MPI_STATUS_IGNORE,ierr)
            if (lALI) then
              call MPI_WAIT(MPID%requestA(iproc,6), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,7), &
                            MPI_STATUS_IGNORE,ierr)
            end if ! ALI

          end do ! Processors

        end if ! Domain decomposition

        ! If going out, but with no ALI
        if (goout.and..not.lALI.and.INPUT%ALI_force) then

          ! Don't leave and activate ALI
          goout = .False.
          force_ALI = .True.

        end if

        ! If we are finished, everyone exits
        if (goout) exit

      end do ! Iterations

      ! If 1.5D, Master save MRC
      if (pid.eq.0.and.run_mode.eq.1.and.Input%keep_MRC) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(Input%folder)//'/MRC', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1001

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1011
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1011

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = 0e0

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1101

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      ! Control
      call control

      !
      ! Clean slave pointers
      !
2000  if (pid.ne.0) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
        deallocate(rLineO,rPhotO)
        nullify(rLineO,rPhotO)
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO,p_LO)
      end if

      return

1000  umsg = 'Error opening MRCI file'
      urou = 'solverI'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing MRCI file'
      urou = 'solverI'
      close(800)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening MRC file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1011  umsg = 'Error seeking MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine solverI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem for intensity using the ALI method,
      !! with several CPU. Alternative communication scheme.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!            lload(logical): Bool that says if previous
      !!                            solution was read\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine solverI_alt(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                             Red,Geom,MPID,Input,lload,Stokes, &
                             J00,J00S,J00C,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      logical:: lload
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00
      double precision, dimension(nxt,Rz0:Rz1):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P

      ! Local


      type(MRC_class):: MRC

      logical:: goout,gooutprd,AD,ADT,ADD,PRDl,doNG,laux,lALI
      logical:: RIRAM,deal,force_ALI

      character(LEN=20):: iterS

      integer:: iaux,ios,iz0,iz1,diz,m,o,p,op
      integer:: ith,iph,iph1,ia,iz,ifreq,if0,if1,iil,iip
      integer:: iproc,iter,itran,jtran,ftran,jftran,fftran,id,iterr
      integer:: if0l,if1l,ip0l,ip1l,nfl,nftl,nfpl,istep,ith1
      integer:: NG_dim,NG_entry,iterm,iJ,ing,tid

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset
      double precision:: daux,WA,mu_inv,dsm,dsp
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, &
             dimension(2,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: Norm
      double precision, &
             dimension(2,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: BStk
      double precision, &
             dimension(nxb,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: BLam
      double precision, dimension(nxb,nxt,Rz0:Rz1):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1):: LambdaP
      double precision, dimension(nfreq,Rz0:Rz1):: J00Cold

      double precision, dimension(:), allocatable, target:: LO

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Stokes_r
      double precision, dimension(:), allocatable, target:: LambdaL_r
      double precision, dimension(:), allocatable, target:: LambdaP_r
      double precision, dimension(:), allocatable, target:: Prof_r
      ! Senders
      double precision, dimension(:,:), allocatable:: Stokes_s
      double precision, dimension(:,:,:), allocatable:: Prof_s
      double precision, dimension(:,:), allocatable:: rLine_s
      double precision, dimension(:,:), allocatable:: rPhot_s
      ! Dual
      integer:: info_b
      integer, dimension(3):: info_c


      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:,:), pointer:: data2O,data2P
      double precision, dimension(:), pointer:: rLineO,rLineP
      double precision, dimension(:), pointer:: rPhotO,rPhotP
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP, p_LO
      double precision, dimension(:), pointer:: p_exu
      double precision, dimension(:,:), pointer:: p_MStk
      double precision, dimension(:,:), pointer:: p_MrLine,p_MrPhot
      double precision, dimension(:,:,:), pointer:: p_MProf


      ! Initialize converged flag
      goout = .False.

      ! Initialize force ALI
      force_ALI = .False.

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! Really go back to angle dependent
      ADT = .not.AV.and.AD

      ! If storing redistribution
      RIRAM = IRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AVI = .True.
        if (ADT) IRAM = .False.
      end if

      ! Store if it was PRD (first iteration always CRD unless we are
      ! reading)
      PRDl = PRD
      if (.not.lload) PRD = .False.


      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.

      ! If NG acceleration of rho00
      if (Input%NGI) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + RnZ*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRDl) then

          ! If we need Stokes
          if (ADD) then

            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            NG_dim = NG_dim + nfreq*RnZ

          end if ! ADD
        end if ! PRD

        ! If it requires too much buffer
        if (dble(NG_dim)*8d-6.gt.maxbuffer_NG) then

          Input%NGI = .False.

          if (pid.eq.0) then
            umsg = ' # The buffer for NG acceleration '// &
                   'is too big. Not doing NG.'
            call verbose
          end if ! Master

        end if ! Buffer

        ! If finally doing it, allocate
        if (Input%NGI) then

          ! Master
          if (pid.eq.0) then

            allocate(NG_scratch(NG_dim, Input%NGI_ord+2))

          ! Slaves
          else

            allocate(NG_scratch(NG_dim, 1))

          end if

        end if ! NG and master
      end if ! NG

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)

      ! Master
      if(pid.eq.0)then

        ! To receive Intensity chunks
        iaux = MPID%nxfreq*Rnz
        allocate(Stokes_r(iaux))

        ! To receive Lambda operator for b-b transitions
        iaux = MPID%nxtfreqi*nxb*Rnz
        allocate(LambdaL_r(iaux))

        ! To receive Lambda operator for b-f transitions
        iaux = MPID%nxpfreq*nxb*Rnz
        allocate(LambdaP_r(iaux))

        ! To receive profile information
        iaux = MPID%nxtfreqi*2*Rnz
        allocate(Prof_r(iaux))

      ! Slave
      else


        !
        ! Allocation of buffers
        !

        ! Common (Master and slave)
        ! Allocate O pointers
        allocate(data2O(Frec%ntfreqi,2))
        allocate(rLineO(Frec%ntfreqi))
        allocate(rPhotO(Frec%npfreq))

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(if0:if1,0:2))
        allocate(data1O(if0:if1,0:2))

        ! Allocate vector for Lambda operator
        allocate(LO(if0:if1))

        ! To send Intensity chunks
        allocate(Stokes_s(if0:if1,Rz0:Rz1))

        ! To send profile information
        allocate(Prof_s(Frec%ntfreqi,2,Rz0:Rz1))

        ! To send Lambda operator for b-b transitions
        allocate(rLine_s(Frec%ntfreqi,Rz0:Rz1))

        ! To send Lambda operator for b-f transitions
        allocate(rPhot_s(Frec%npfreq,Rz0:Rz1))

      end if ! Master or Slave


      !
      ! Initialization messages
      !

      ! Only global Máster
      if (gpid.eq.0) then

        ! Announce we are starting
        umsg = '    Iteration          MRC(rho^0_0) Atom_index '// &
               'Level_index Height_index Height(km)'
        call verbose

      end if

      ! Only global Máster
      if (gpid.eq.0.and.Input%keep_MRC) then

        ! Open the file to store MRC

        ! If appending
        if (Input%appendMRC) then
          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRCI', exist=laux)
          if (.not.laux) then
            open(800, file=trim(Input%folder)//'/MRCI', &
                 action='write',iostat=ios,err=1000)
            write(800,'(A)') '!   Iteration          MRC(rho^0_0) '//&
                         'Atom_index Level_index Height_index '// &
                         'Height(km)'
            close(800)
          end if
        ! If not appending, create
        else
          open(800, file=trim(Input%folder)//'/MRCI', &
               action='write',iostat=ios,err=1000)
          write(800,'(A)') '!   Iteration          MRC(rho^0_0) '// &
                       'Atom_index Level_index Height_index '// &
                       'Height(km)'
          close(800)
        end if
      end if

      ! Control
      call control
      if (laborted) goto 2000

      ! Initialize stimulated radiation tensor, becuase it is
      ! commented in the main loop
      J00S = 0d0

      ! If measuring performance
      if (Input%mpi_perf.and.pid.eq.0) &
        call report_mpi_timeI(Input%folder,Input%ID, &
                             0,0,0,.False.)


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=Input%iteri_min,Input%iteri_max

        ! Flags for physics
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics
        if (iter.le.Input%allownphys_pop) then
          if (.not.nphysR) nphysR = .True.
        else
          if (nphysR) nphysR = .False.
        end if

        !
        ! The old population will be the current one
        !

        ! For each atom
        do ia=1,nA
          Rho_old(ia)%crho = Atom(ia)%crho
        end do

        ! Flag for ALI
        lALI = iter.gt.Input%ALI_delay.or.force_ALI

        ! Internal PRD iterations
        do iterr=1,Input%iteri_prd

          !
          ! Master
          !
          if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(ith1,iph1,iz,id,info_c,WA,deal,ia,ftran,jftran,fftran) &
!$omp private(daux,dsm,p_exu,jtran) &
!$omp shared(J00,J00P,J00C,BStk,BLam,Norm,lALI,LambdaL,LambdaP) &
!$omp shared(info_b,Stokes_r,Prof_r,LambdaL_r,LambdaP_r,J00Cold) &
!$omp shared(ith,iph,laborted,if0l,if1l,ip0l,ip1l,nfl,nftl,nfpl) &
!$omp shared(iter,iterr,p_MStk,p_MProf,p_MrLine,p_MrPhot,Stokes) &
!$omp shared(nz,nxb,na,KSTK,PIRAM,PRD,Rz0,Rz1,Rnz,giz0) &
!$omp shared(MPID,Frec,Atom,Geom,Atmo,Input,ierr) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT)

            ! Save old radiation field
            if (PRD) then
!$omp workshare
              J00Cold = J00C
!$omp end workshare
            end if
            ! Reset radiation field variables
!$omp workshare
            J00 = 0d0
           !J00S = 0d0
            J00P = 0d0
            J00C = 0d0
            BStk = 0d0
            BLam = 0d0
            Norm = 0d0
!$omp end workshare
            if (lALI) then
!$omp workshare
              LambdaL = 0d0
              LambdaP = 0d0
!$omp end workshare
            end if

            ! For each dimension in the problem
            ! Each polar direction
            do ith1=1,Geom%nTh

              ! Each azimuthal direction
              do iph1=1,Geom%nPh

                ! Each frequency cut
                do id=1,MPID%nnd

!$omp single
                  !
                  ! Receive data from a slave
                  !

                  ! Receive indexing data
                  do while (.True.)
                    call MPI_recv(info_c(1),3,MPI_INTEGER, &
                                  MPI_ANY_SOURCE, 0, &
                                  MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)
                    if (ierr.eq.0) exit
                  end do

                  info_b = info_c(1)
                  ith = info_c(2)
                  iph = info_c(3)

                  ! Flag error
                  if (info_b.lt.0) laborted = .True.
!$omp end single
                  ! Continue?
                  if (info_b.lt.0) cycle
!$omp single
                  ! Receive intensity
                  do while (.True.)
                    call MPI_recv(Stokes_r(1), MPID%sizei4(info_b), &
                                  MPI_DOUBLE_PRECISION, info_b, &
                                  info_b, MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)
                    if (ierr.eq.0) exit
                  end do

                  ! Receive profile
                  do while (.True.)
                    call MPI_recv(Prof_r(1), MPID%sizei5(info_b), &
                                  MPI_DOUBLE_PRECISION, info_b, &
                                  1000000+info_b, MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)
                    if (ierr.eq.0) exit
                  end do

                  ! If ALI
                  if (lALI) then

                    ! Receive Lambda operator for b-b transition
                    do while (.True.)
                      call MPI_recv(LambdaL_r(1), &
                                    MPID%sizei9(info_b), &
                                    MPI_DOUBLE_PRECISION, info_b, &
                                    2000000+info_b, &
                                    MPI_COMM_RT, &
                                    MPI_STATUS_IGNORE, ierr)
                      if (ierr.eq.0) exit
                    end do

                    ! Receive Lambda operator for b-f transition
                    do while (.True.)
                      call MPI_recv(LambdaP_r(1), &
                                    MPID%sizei0(info_b), &
                                    MPI_DOUBLE_PRECISION, info_b, &
                                    3000000+info_b, &
                                    MPI_COMM_RT, &
                                    MPI_STATUS_IGNORE, ierr)
                      if (ierr.eq.0) exit
                    end do

                  end if ! ALI

                  ! If measuring performance
                  if (Input%mpi_perf) &
                    call report_mpi_timeI(Input%folder,Input%ID, &
                                          info_b,iter,iterr,.True.)

                  ! Shorter variables
                  if0l = MPID%if0(info_b)
                  if1l = MPID%if1(info_b)
                  ip0l = Frec%Mpif0(info_b)
                  ip1l = Frec%Mpif1(info_b)
                  nfl = MPID%nf(info_b)
                  nftl = Frec%Mntfreqi(info_b)
                  nfpl = Frec%Mnpfreq(info_b)

                  ! Pointers
                  p_MStk(if0l:if1l,1:Rnz) => &
                                      Stokes_r(1:MPID%sizei4(info_b))
                  p_MProf(1:nftl,1:2,1:Rnz) => &
                                        Prof_r(1:MPID%sizei5(info_b))
                  if (lALI) then
                    p_MrLine(1:nftl,1:Rnz) => &
                                     LambdaL_r(1:MPID%sizei9(info_b))
                    p_MrPhot(1:nfpl,1:Rnz) => &
                                     LambdaP_r(1:MPID%sizei0(info_b))
                  end if

!$omp end single

                  ! Get angular weight
                  WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                  deal = .False.

                  ! Each height
!$omp do
                  do iz=Rz0,Rz1

                    ! Determine where to store intensity
                    if (KSTK.or.iz.eq.giz0) &
                      Stokes(if0l:if1l,iph,ith,iz) = &
                                                  p_MStk(if0l:if1l,iz)
                    ! Point to exu values
                    if (PIRAM.and.ip1l.ge.ip0l) then
                      p_exu => Frec%exu(ip0l:ip1l,iz)
                    else
                      allocate(p_exu(1))
                      deal = .True.
                    end if

                    ! Calculate frequency integral for b-b quantities
                    call FIntI_line(Atom,MPID,Frec%W_freq,info_b, &
                                    p_MStk(:,iz),p_MrLine(:,iz), &
                                    p_MProf(:,:,iz), &
                                    Norm(:,:,iph,ith,iz), &
                                    BStk(:,:,iph,ith,iz), &
                                    BLam(:,:,iph,ith,iz), &
                                    lALI)

                    ! Calculate rest of integrals
                    call FIntI_rest(Atom,MPID,Frec%omega, &
                                    Frec%W_freq,ip0l,ip1l, &
                                    Atmo%T(iz),info_b,WA, &
                                    p_MStk(:,iz),p_MrPhot(:,iz), &
                                    J00P(:,:,iz),J00C(:,iz), &
                                    LambdaP(:,:,:,iz),lALI, &
                                    p_exu)

                    ! Nullify exponential pointer
                    if (deal) deallocate(p_exu)
                    nullify(p_exu)

                  end do ! heights
!$omp end do
                end do ! frequency domains
              enddo ! azimuthal directions
            enddo ! polar directions
!$omp single
            ! Nullify pointers
            nullify(p_MStk,p_MProf,p_MrLine,p_MrPhot)
!$omp end single


            !
            ! Apply weights to J00, J00S, and Lambda operator
            ! and normalize
            !

            ! For each height
!$omp do
            do iz=Rz0,Rz1

              ! For each polar direction
              do ith1=1,Geom%nTh

                ! For each azimuthal direction
                do iph1=1,Geom%nph

                  ! Get the angular integral weight
                  WA = Geom%W_mu(ith1)*Geom%W_mux(iph1)

                  ! For each atom
                  do ia=1,nA

                    ! For each FS transition
                    do ftran=1,Atom(ia)%nftran

                      ! Apply atomic shift
                      jftran = ftran + Atom(ia)%tfshift

                      ! Get the weight
                      if (Norm(1,jftran,iph1,ith1,iz).gt.0d0) then

                        daux = WA/Norm(1,jftran,iph1,ith1,iz)

                        ! Integrate angle
                        J00(jftran,iz) = J00(jftran,iz) + &
                                        BStk(1,jftran,iph1,ith1,iz)* &
                                        daux

                        ! For each transition blended with ftran
                        if (lALI) then

                          do fftran=1,nxb
                            LambdaL(fftran,jftran,iz) = &
                                   LambdaL(fftran,jftran,iz) + &
                                   BLam(fftran,jftran,iph1,ith1,iz)* &
                                   daux
                          end do

                        end if ! ALI
                      end if ! Norma

                      ! If there is stimulated emission
                     !if (stm) then

                     !  ! Get the weight
                     !  if (Norm(2,jftran,iph1,ith1,iz).gt.0d0) then

                     !    daux = WA/Norm(2,jftran,iph1,ith1,iz)

                     !    ! Integrate angle
                     !    J00S(jftran,iz) = J00S(jftran,iz) + &
                     !                  BStk(2,jftran,iph1,ith1,iz)* &
                     !                  daux

                     !  end if

                     !end if ! Stimulated emission

                    end do ! FS transition
                  end do ! Atoms
                end do ! azimuthal directions
              end do ! polar directions

              !
              ! Add the Saha factor (ne*Zeta) to J00P
              !

              ! Argument of the exponential
              WA = fktoJ/kb/Atmo%T(iz)

              ! Part that does not depend on the line
              daux = cSaha*Atmo%ne(iz)/(Atmo%T(iz)**(1.5d0))

              ! For each atom
              do ia=1,nA

                ! For each b-f transition
                do itran=1,Atom(ia)%nphot

                  ! Apply atomic shift
                  jtran = itran + Atom(ia)%pshift

                  ! Calculate the multiplicative factor
                  dsm = daux*exp(Atom(ia)%phot(itran)%edge*WA)* &
                        Atom(ia)%phot(itran)%glu

                  ! Apply it to the emission integral
                  J00P(jtran,2,iz) = J00P(jtran,2,iz)*dsm

                  ! Apply it to the emission Lambda operator
                  if (lALI) &
                    LambdaP(:,jtran,2,iz) = LambdaP(:,jtran,2,iz)*dsm

                end do ! b-f transitions
              end do ! atoms
            end do ! heights
!$omp end do
!$omp end parallel

            ! Calculate MRC for J if PRD
            if (PRD.and.Input%iteri_prd.gt.1) then

              ! Call the routine
              call MRCJ_sb(J00C,J00Cold,MRC)

              ! Convert cm into km
              MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

              ! Global Máster talks
              if (gpid.eq.0) then

                if (iterr.eq.1.and.iter.eq.2) then
                  umsg = '         PRD            MRC(J^0_0)'// &
                         ' Freq_index  Wavelength '// &
                         'Height_index Height(km)'
                  call verbose
                end if

                ! Write in stdout
                write(umsg,'(2x,"PRD it:",1x,i3,2x,es20.12,'// &
                           '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                           iterr,MRC%values(2,1),MRC%indexes(1,1), &
                           1d2/Frec%omega(MRC%indexes(1,1)), &
                           MRC%indexes(2,1),MRC%values(1,1)
                call verbose

              end if

              ! Check exit criteria
              if (MRC%values(2,1).le.Input%mrci_r.or. &
                  iterr.eq.Input%iteri_prd) then
                gooutprd = .True.
              else
                gooutprd = .False.
              end if

            ! Always go out if no PRD
            else

              gooutprd = .True.

            end if

            ! If aborting, obviously go out
            if (laborted) gooutprd = .True.


          !
          ! Slave
          !
          else

            !
            ! Ratiation Transfer
            !

            !  For each polar direction
            do ith=1,Geom%nTh

              ! Calculate inverse of cosine of polar direction
              mu_inv = 1d0/Geom%V_mu(ith)

              ! Determine the direction of propagation for indexes
              diz = -int(sign(1d0, Geom%V_mu(ith)))

              ! Determine the first and last height indexes to run
              ! over

              iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
              iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

              ! For each azimuthal direction
              do iph=1,Geom%nPh

                !
                ! If angle-dependent, manage scattering angles
                if (.not.AVI.and.PRD) &
                  call scattering_manage(Geom,ith,iph)

                !
                ! First height
                !

                ! If going down, get top boundary
                if(diz.eq.1)then

                  ! Call top boundary
                  call topI(MPID,data1M(:,2))

                ! If going up, get bottom boundary
                else

                  ! Call bottom boundary
                  call bottomI(Frec%omega,Atmo%T(iz0), &
                               Atmo%vx(iz0),Atmo%vy(iz0), &
                               Atmo%vz(iz0),Geom%V_mu(ith), &
                               Geom%V_mux(iph),Geom%V_muy(iph), &
                               MPID,data1M(:,2))

                endif ! propagation direction

                ! Identify current height
                o = iz0

                ! Index for Stokes
                if (PRDl.and.ADD) op = o

                ! Calculate radiative coefficients
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              o,ith,iph,if0,if1,J00(:,o), &
                              J00C(:,o),Cont%ndir, &
                              Cont%c(:,:,:,o),Stokes(:,:,:,op), &
                              rLineO(:),rPhotO(:), &
                              data1M(:,0:1),data2O(:,:),.True.)
                if (laborted) goto 3000


                !
                ! Store in buffer
                !

                ! Intensity
                Stokes_s(:,o) = data1M(:,2)

                ! Profiles
                Prof_s(:,:,o) = data2O(:,:)

                ! ALI
                if (lALI) then

                  ! b-b Lambda operator (bottom boundary does not
                  ! contribute)
                  rLine_s(:,o) = 0d0

                  ! b-f Lambda operator (bottom boundary does not
                  ! contribute)
                  rPhot_s(:,o) = 0d0

                end if ! ALI

                ! Identify next height
                p = iz0 + diz

                ! Index for Stokes
                if (PRDl.and.ADD) op = p

                ! Calculate radiative coefficients
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              p,ith,iph,if0,if1,J00(:,p), &
                              J00C(:,p),Cont%ndir, &
                              Cont%c(:,:,:,p),Stokes(:,:,:,op), &
                              rLineO(:),rPhotO(:), &
                              data1O(:,0:1),data2O(:,:),.True.)
                if (laborted) goto 3000

                !
                ! Intermedium heights
                !

                ! For each height this CPU has assigned
                do iz=iz0,iz1,diz

                  ! We treat the boundaries outside
                  if(iz.eq.iz0.or.iz.eq.iz1)cycle

                  ! Allocate P pointers
                  allocate(data1P(if0:if1,0:2))
                  allocate(data2P(Frec%ntfreqi,2))
                  allocate(rLineP(Frec%ntfreqi))
                  allocate(rPhotP(Frec%npfreq))

                  ! Identify heights
                  m = iz - diz
                  o = iz
                  p = iz + diz

                  ! Calculate distance to previous point
                  dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                  ! Calculate distance to the next point
                  dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

                  ! If tau scale
                  if (ztau) then
                    dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                        Atmo%chi500(m))
                    dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                        Atmo%chi500(p))
                  end if

                  ! Index for Stokes
                  if (PRDl.and.ADD) op = p

                  ! Calculate radiative coefficients
                  call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID, &
                                Geom,p,ith,iph,if0,if1,J00(:,p), &
                                J00C(:,p),Cont%ndir, &
                                Cont%c(:,:,:,p),Stokes(:,:,:,op), &
                                rLineP(:),rPhotP(:), &
                                data1P(:,0:1),data2P(:,:),.True.)

                  ! Point to the data
                  p_K0M  => data1M(:,0)
                  p_SM   => data1M(:,1)
                  p_StkM => data1M(:,2)
                  p_K0O  => data1O(:,0)
                  p_SO   => data1O(:,1)
                  p_StkO => data1O(:,2)
                  p_K0P  => data1P(:,0)
                  p_SP   => data1P(:,1)
                  p_LO => LO(:)

                  ! Apply short characteristics BESSER
                  call RTStepI(o,ith,iph,MPID%nf(pid), &
                               dsm,dsp,p_K0M,p_SM,p_K0O, &
                               p_SO,p_K0P,p_SP,p_StkM, &
                               p_StkO,p_LO,lALI,.True.)

                  !
                  ! Combine the value of lambda operator with the
                  ! transition strength
                  !
                  if (lALI) then

                    ! Initialize indexes for rLine and rPhot
                    iil = 0
                    iip = 0

                    ! For each atom
                    do ia=1,nA

                      ! For each b-b trantision
                      do itran=1,Atom(ia)%ntran

                        ! If this CPU does not have frequencies in
                        ! this line, skip
                        if (Atom(ia)%fflag(itran)%absent) cycle

                        ! For each FS transition
                        do ftran=1,Atom(ia)%fst(itran)%nt

                          ! Get the sequential index of FS transition
                          fftran = Atom(ia)%ifst_ij(ftran,itran)

                          ! Apply atomic shift
                          jftran = fftran + Atom(ia)%tfshift

                          ! For each frequency
                          do ifreq=Atom(ia)%if0(itran), &
                                   Atom(ia)%if1(itran)

                            iil = iil + 1
                            rLineO(iil) = LO(ifreq)*rLineO(iil)

                          end do ! frequency
                        end do ! FS transition
                      end do ! b-b transition

                      ! For each b-f transition
                      do itran=1,Atom(ia)%nphot

                        ! If this CPU does not have frequencies in
                        ! this transition, skip
                        if (Atom(ia)%phot(itran)%absent) cycle

                        ! Apply atomic shift
                        jtran = itran + Atom(ia)%pshift

                        ! For each frequency
                        do ifreq=Atom(ia)%phot(itran)%if0, &
                                 Atom(ia)%phot(itran)%if1

                          iip = iip + 1
                          rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                        end do ! frequency
                      end do ! b-f transition
                    end do ! atom

                    ! b-b Lambda operator
                    rLine_s(:,o) = rLineO(:)

                    ! Send b-f Lambda operator
                    rPhot_s(:,o) = rPhotO(:)

                  end if ! ALI


                  !
                  ! Store in buffer
                  !

                  ! Intensity
                  Stokes_s(:,o) = data1O(:,2)

                  ! Profiles
                  Prof_s(:,:,o) = data2O(:,:)

                  ! Shift data (O->M, P->O)
                  deallocate(data1M,data2O)
                  data1M => data1O
                  data1O => data1P
                  data2O => data2P
                  nullify(data1P,data2P)
                  if (lALI) then
                    deallocate(rLineO,rPhotO)
                    rLineO => rLineP
                    rPhotO => rPhotP
                  else
                    deallocate(rLineP,rPhotP)
                  end if
                  nullify(rLineP,rPhotP)

                  ! Error
                  if (laborted) goto 3000

                end do ! Intermedium heights


                !
                ! Last height
                !

                ! Identify heights
                m = iz1 - diz
                o = iz1

                ! Calculate distance to previous point
                dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                ! If tau scale
                if (ztau) &
                  dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(m))

                ! Point to the data
                p_K0M  => data1M(:,0)
                p_SM   => data1M(:,1)
                p_StkM => data1M(:,2)
                p_K0O  => data1O(:,0)
                p_SO   => data1O(:,1)
                p_StkO => data1O(:,2)
                p_LO => LO(:)

                ! Apply short characteristics LINEAR
                call RTStepI(o,ith,iph,MPID%nf(pid), &
                             dsm,dsp,p_K0M,p_SM,p_K0O, &
                             p_SO,p_K0P,p_SP,p_StkM, &
                             p_StkO,p_LO,lALI,.False.)
                if (laborted) goto 3000

                !
                ! Combine the value of lambda operator with the
                ! transition strength
                !
                if (lALI) then

                  ! Initialize indexes for rLine and rPhot
                  iil = 0
                  iip = 0

                  ! For each atom
                  do ia=1,nA

                    ! For each b-b trantision
                    do itran=1,Atom(ia)%ntran

                      ! If this CPU does not have frequencies in
                      ! this line, skip
                      if (Atom(ia)%fflag(itran)%absent) cycle

                      ! For each FS transition
                      do ftran=1,Atom(ia)%fst(itran)%nt

                        ! Get the sequential index of FS transition
                        fftran = Atom(ia)%ifst_ij(ftran,itran)

                        ! Apply atomic shift
                        jftran = fftran + Atom(ia)%tfshift

                        ! For each frequency
                        do ifreq=Atom(ia)%if0(itran), &
                                 Atom(ia)%if1(itran)

                          iil = iil + 1
                          rLineO(iil) = LO(ifreq)*rLineO(iil)

                        end do ! frequency
                      end do ! FS transition
                    end do ! b-b transition

                    ! For each b-f transition
                    do itran=1,Atom(ia)%nphot

                      ! If this CPU does not have frequencies in
                      ! this transition, skip
                      if (Atom(ia)%phot(itran)%absent) cycle

                      ! Apply atomic shift
                      jtran = itran + Atom(ia)%pshift

                      ! For each frequency
                      do ifreq=Atom(ia)%phot(itran)%if0, &
                               Atom(ia)%phot(itran)%if1

                        iip = iip + 1
                        rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                      end do ! frequency
                    end do ! b-f transition
                  end do ! atom

                  ! b-b Lambda operator
                  rLine_s(:,o) = rLineO(:)

                  ! b-f Lambda operator
                  rPhot_s(:,o) = rPhotO(:)

                end if ! ALI

                !
                ! Send to master if had problem
                !

                ! If had an error
3000            if (laborted) then

                  ! Send error
                  info_c = (/ -pid, ith, iph /)
                  do while (.True.)
                    call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                                  MPI_COMM_RT,ierr)
                    if (ierr.eq.0) exit
                  end do

                  cycle

                end if


                !
                ! Store in buffer
                !

                ! Intensity
                Stokes_s(:,o) = data1O(:,2)

                ! Profiles
                Prof_s(:,:,o) = data2O(:,:)

                !
                ! Send to master
                !

                ! Send indexes
                info_c = (/ pid, ith, iph /)
                do while (.True.)
                  call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                                MPI_COMM_RT, ierr)
                  if (ierr.eq.0) exit
                end do

                ! Send Stokes
                do while (.True.)
                  call MPI_SEND(Stokes_s(if0,Rz0), MPID%sizei4(pid), &
                                MPI_DOUBLE_PRECISION, 0, pid, &
                                MPI_COMM_RT, ierr)
                  if (ierr.eq.0) exit
                end do

                ! Send profiles
                do while (.True.)
                  call MPI_SEND(Prof_s(1,1,Rz0), &
                                MPID%sizei5(pid), &
                                MPI_DOUBLE_PRECISION, &
                                0, 1000000+pid, MPI_COMM_RT, &
                                ierr)
                  if (ierr.eq.0) exit
                end do

                ! ALI
                if (lALI) then

                  ! Send b-b Lambda operator
                  do while (.True.)
                    call MPI_SEND(rLine_s(1,Rz0), &
                                  MPID%sizei9(pid), &
                                  MPI_DOUBLE_PRECISION, &
                                  0, 2000000+pid, MPI_COMM_RT, &
                                  ierr)
                    if (ierr.eq.0) exit
                  end do

                  ! Send b-f Lambda operator
                  do while (.True.)
                    call MPI_SEND(rPhot_s(1,Rz0), &
                                  MPID%sizei0(pid), &
                                  MPI_DOUBLE_PRECISION, &
                                  0, 3000000+pid, MPI_COMM_RT, &
                                  ierr)
                    if (ierr.eq.0) exit
                  end do

                end if ! ALI

              enddo ! azimuthal angles
            enddo ! polar angles

          end if ! Master or Slave

          !
          ! Share if we are finished or not (in PRD)
          !

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive goout
              call MPI_RECV(gooutprd, 1, MPI_LOGICAL, &
                            MPID%recv, 7000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No Master

            ! For each send
            do istep=1,MPID%nsend

              ! Send goout
              call MPI_ISEND(gooutprd, 1, MPI_LOGICAL, &
                             MPID%lsend(istep), &
                             7000000+MPID%lsend(istep), &
                             MPI_COMM_RT, MPID%requestA(istep,1), &
                             ierr)
            end do ! sends

          ! Normal bcast
          else

            ! Share goout
            call MPI_BCAST(gooutprd, 1, MPI_LOGICAL, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast

          ! Control
          call control
          if (laborted) goto 2000


          !
          ! Share the radiation information
          !

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive J00
              call MPI_RECV(J00(1,Rz0), MPID%sizei6(0), &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

             !! Receive J00S if stimulated emission
             !if (stm) &
             !  call MPI_RECV(J00S(1,Rz0), MPID%sizei6(0), &
             !                MPI_DOUBLE_PRECISION,  &
             !                MPID%recv, 1000000+pid, &
             !                MPI_COMM_RT, MPI_STATUS_IGNORE, &
             !                ierr)

              ! Receive J00C
              call MPI_RECV(J00C(1,Rz0), MPID%sizei7(0), &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, 2000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

              ! Receive intensity if doing A-D PRD
              if (PRDl.and.ADD) &
                call MPI_RECV(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                              MPI_DOUBLE_PRECISION,  &
                              MPID%recv, 3000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              ! If we are going to do SEE
              if (gooutprd) then

                ! Receive J00 for b-f transitions
                call MPI_RECV(J00P(1,1,Rz0), MPID%sizei3(0), &
                              MPI_DOUBLE_PRECISION,  &
                              MPID%recv, 4000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

                ! ALI
                if (lALI) then

                  ! Receive Lambda operator for b-b transitions
                  call MPI_RECV(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, 5000000+pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                  ! Receive Lambda operator for b-f transitions
                  call MPI_RECV(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, 6000000+pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                end if ! ALI
              end if ! To do SEE

            end if ! No master

            ! For each send
            do istep=1,MPID%nsend

              ! Send J00
              call MPI_ISEND(J00(1,Rz0), MPID%sizei6(0), &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,1), ierr)

             !! Send J00 if stimulated emission
             !if (stm) &
             !  call MPI_ISEND(J00S(1,Rz0), MPID%sizei6(0), &
             !                 MPI_DOUBLE_PRECISION, &
             !                 MPID%lsend(istep), &
             !                 1000000+MPID%lsend(istep), &
             !                 MPI_COMM_RT, &
             !                 MPID%requestA(istep,2), ierr)

              ! Send J00C
              call MPI_ISEND(J00C(1,Rz0), MPID%sizei7(0), &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             2000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,3), ierr)

              ! Send intensity if doing A-D PRD
              if (PRDl.and.ADD) &
                call MPI_ISEND(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               3000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,4), ierr)

              ! If we are going to do SEE
              if (gooutprd) then

                ! Send J00P
                call MPI_ISEND(J00P(1,1,Rz0), MPID%sizei3(0), &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               4000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,5), ierr)

                ! ALI
                if (lALI) then

                  ! Send Lambda operator for b-b transitions
                  call MPI_ISEND(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 5000000+MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,6), ierr)

                  ! Send Lambda operator for b-f transitions
                  call MPI_ISEND(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 6000000+MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,7), ierr)

                end if ! ALI
              end if ! To do SEE

            end do ! Sends

          ! Normal bcast
          else

            ! Share J00
            call MPI_BCAST(J00(1,Rz0), MPID%sizei6(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

           !! Share J00S
           !if (stm) &
           !  call MPI_BCAST(J00S(1,Rz0), MPID%sizei6(0), &
           !                 MPI_DOUBLE_PRECISION, 0, &
           !                 MPI_COMM_RT, ierr)

            ! Share J00C
            call MPI_BCAST(J00C(1,Rz0), MPID%sizei7(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Share intensity if doing A-D PRD
            if (PRDl.and.ADD) &
              call MPI_BCAST(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                            MPI_DOUBLE_PRECISION, 0, &
                            MPI_COMM_RT, ierr)

            ! If we are going to do SEE
            if (gooutprd) then

              ! Share J00 for b-f transitions
              call MPI_BCAST(J00P(1,1,Rz0), MPID%sizei3(0), &
                            MPI_DOUBLE_PRECISION, 0, &
                            MPI_COMM_RT, ierr)

              ! ALI
              if (lALI) then

                ! Share Lambda operator for b-b transitions
                call MPI_BCAST(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                              MPI_DOUBLE_PRECISION, 0, &
                              MPI_COMM_RT, ierr)

                ! Share Lambda operator for b-f transitions
                call MPI_BCAST(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                              MPI_DOUBLE_PRECISION, 0, &
                              MPI_COMM_RT, ierr)

              end if ! ALI
            end if ! To do SEE

          end if ! Type of bcast

          ! If we are finished, everyone exits
          if (.not.gooutprd.and.MPID%altbcast) then

            ! For each slave
            do iproc=1,MPID%nsend

              ! Wait for everyone to receive the radiation data
              ! before continuing and reset the buffers
              call MPI_WAIT(MPID%requestA(iproc,1), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,2), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,3), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,4), &
                            MPI_STATUS_IGNORE,ierr)

              end do ! Processors

          end if ! If not going out

          ! If we are finished, everyone exits
          if (gooutprd) exit

        end do ! PRD iteration


        !
        ! Solve SEE
        !
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P, &
                                  Input%folder,iter)
#endif
#ifdef DEBUGLAMBDA
      if (pid.eq.0) call dump_lambda(Atom,LambdaL,LambdaP, &
                                     Input%folder,iter)
#endif

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,fftran,jftran,umsg,urou,tid) &
!$omp shared(nA,nz,Atom,Rho_old,J00,J00P,LambdaL,LambdaP,lALI) &
!$omp shared(laborted,vaborted,Rz0,Rz1)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
           !          J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL(:,itran:jtran,iz), &
                      LambdaP(:,fftran:jftran,:,iz),iz,lALI,tid,INPUT)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
           !          J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL(:,itran:jtran,iz), &
                      LambdaP(:,fftran:jftran,:,iz),iz,lALI,tid)
#endif

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel
#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,iter)
#endif

        ! Control
        call control
        if (laborted) goto 2000


        !
        ! NG acceleration
        !

        ! Check if doing NG acceleration
        if(Input%NGI.and.iter.gt.Input%NGI_delay)then

          ! Advance ntry. The master does not need it if
          ! if accelerated the intensity
          NG_entry = NG_entry + 1

          ! If Master
          if (pid.eq.0) then

            ! Initialize index
            o = 0

            ! For each atom
            do ia=1,NA

              ! For each term
              do iterm=1,Atom(ia)%nMulti

                ! For each level
                do iJ=1,Atom(ia)%nJ(iterm)

                  ! Get level index
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    o = o + 1

                    ! Store rho
                    NG_scratch(o,NG_entry) = dble(Atom(ia)%crho(p,iz))

                  enddo ! Heights
                enddo ! Levels
              enddo ! Term
            enddo ! Atoms

            ! If PRD
            if (PRDl) then

              ! If ADD
              if (ADD) then

                ! For each height
                do iz=Rz0,Rz1

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Advance index
                        o = o + 1

                        ! Store Stokes
                        NG_scratch(o,NG_entry) = &
                                            Stokes(ifreq,iph,ith,iz)

                      end do ! Frequencies
                    end do ! Azimuthal directions
                  end do ! Polar directions
                end do ! Heights

              ! If AV
              else

                ! For each height
                do iz=Rz0,Rz1

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! Advance index
                    o = o + 1

                    ! Store Stokes
                    NG_scratch(o,NG_entry) = J00C(ifreq,iz)

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD or dyn
            end if ! PRD

            ! Call NG and check if it should be processed
            call NG(NG_dim,Input%NGI_ord,NG_scratch,NG_entry,doNG)

          ! Slave
          else

            ! If wrong order
            if (Input%NGI_ord.lt.1.or.Input%NGI_ord.gt.5) then
              doNG = .False.
            ! Valid order
            else
              ! Check if Master is in NG step
              if (NG_entry.gt.(Input%NGI_ord+1)) then
                doNG = .True.
              ! Not a NG step
              else
                doNG = .False.
              end if ! NG step
            end if ! order

          end if ! Master of slave

          ! If communication is needed
          if (doNG) then

            ! If Master, send last
            if (pid.eq.0) then
              ing = NG_entry
            ! Slaves send the one they have
            else
              ing = 1
            end if

            ! Alternative bcast
            if (MPID%altbcast) then

              ! If not master, receive first
              if (pid.ne.0) then

                ! Receive NG iteration
                call MPI_RECV(NG_scratch(1,1), NG_dim, &
                              MPI_DOUBLE_PRECISION,  &
                              MPID%recv, 7000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No master

              ! For each send
              do istep=1,MPID%nsend

                ! Send NG iteration
                call MPI_ISEND(NG_scratch(1,ing), NG_dim, &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               7000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,8), ierr)

              end do ! Sends

            ! Normal bcast
            else

              ! Share NG iteration
              call MPI_BCAST(NG_scratch(1,ing), NG_dim, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

            end if ! Type of bcast

            ! Reconstruct NG data

            ! Initialize index
            o = 0

            ! For each atom
            do ia=1,NA

              ! For each term
              do iterm=1,Atom(ia)%nMulti

                ! For each level
                do iJ=1,Atom(ia)%nJ(iterm)

                  ! Get level index
                  daux = sqrt(2d0*Atom(ia)%rJval(iJ,iterm) + 1d0)
                  m = Atom(ia)%irho(iterm)%irho_ij(iJ)
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Advance index
                    o = o + 1

                    ! Accelerate rho
                    Atom(ia)%crho(p,iz) = &
                                        dcmplx(NG_scratch(o,ing), 0d0)
                    Atom(ia)%popu(m,iz) = NG_scratch(o,ing)*daux

                  enddo ! Heights
                enddo ! Levels
              enddo ! Terms
            enddo ! Atoms

            ! If PRD
            if (PRDl) then

              ! If ADD
              if (ADD) then

                ! For each height
                do iz=Rz0,Rz1

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Advance index
                        o = o + 1

                        ! Store Stokes
                        Stokes(ifreq,iph,ith,iz) = NG_scratch(o,ing)

                      end do ! Frequencies
                    end do ! Azimuthal directions
                  end do ! Polar directions
                end do ! Heights

              ! If AV
              else

                ! For each height
                do iz=Rz0,Rz1

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! Advance index
                    o = o + 1

                    ! Store Stokes
                    J00C(ifreq,iz) = NG_scratch(o,ing)

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD or dyn
            end if ! PRD

            if (gpid.eq.0) then
              umsg = 'NG acceleration'
              call verbose
            end if

            NG_entry = 0

            ! Master wait to free the buffer if alternative bcast
            if (pid.eq.0.and.MPID%altbcast) then

              ! For each slave
              do iproc=1,MPID%nsend

                ! Wait for everyone to receive the radiation data
                ! before continuing and reset the buffers
                call MPI_WAIT(MPID%requestA(iproc,8), &
                              MPI_STATUS_IGNORE,ierr)

              end do ! Processors

            end if ! Master
          end if ! Communication
        endif ! NG acceleration


        !
        !  Calculate MRC
        !

        ! Only the master does
        if (pid.eq.0) then

          ! Calculate MRC
          call MRCI_sb(Atom,Rho_old,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Check exit criteria
          if (MRC%values(2,1).le.Input%mrci_i) goout = .True.

          ! Check in first iteration if was a fix population case
          ! but we want PRD
          if (iter.eq.Input%iteri_min.and.PRDl.and. &
            MRC%values(2,1).lt.1d-16) goout = .False.
          ! If going to AD from AA
          if (iter.eq.Input%iteri_min.and.tbAD.and.ADT) &
            goout = .False.

        end if

        ! We can swith now to AD if we had AV input
        if (tbAD.and.ADT) then
          AVI = .False.
          tbAD = .False.
          IRAM = RIRAM
        end if

        ! Global Máster
        if (gpid.eq.0) then

          ! Write in stdout
          write(umsg,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                     '2x,f9.3)') &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(2,1),MRC%values(1,1)
          call verbose

          ! File
          if (Input%keep_MRC) then

            ! Write in MRC file
            open(800, file=trim(Input%folder)//'/MRCI', &
                 iostat=ios,err=1000,position='append')
            write(800,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                      '2x,f9.3)', err=1100) &
            iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
            MRC%indexes(2,1),MRC%values(1,1)
            close(800)

          end if
        end if


        !
        ! Save partial solution
        !

        ! Check if we want to store partial results
        if(Input%storei.and.mod(iter,Input%storei_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesolI(Input,iterS,Frec%omega,Geom, &
         !               Atom,Atmo%z,Stokes,J00,J00S,J00C,J00P, &
                         Atom,Atmo%z,Stokes,J00,J00,J00C,J00P, &
                         .False.)

        ! Or have a control check
        else

          call control
          if (laborted) goto 2000

        end if

        ! If in the mandatory non-PRD
        if (iter.le.Input%PRD_delay) goout = .False.

        ! Recover the PRD variable
        if (iter.eq.Input%PRD_delay) PRD = PRDl

        !
        ! Share if we are finished or not

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive goout
            call MPI_RECV(goout, 1, MPI_LOGICAL, &
                          MPID%recv, 7000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, ierr)

          end if ! No Master

          ! For each send
          do istep=1,MPID%nsend

            ! Send goout
            call MPI_ISEND(goout, 1, MPI_LOGICAL, &
                           MPID%lsend(istep), &
                           7000000+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,1), &
                           ierr)
          end do ! sends

        ! Normal bcast
        else

          ! Share goout
          call MPI_BCAST(goout, 1, MPI_LOGICAL, 0, &
                         MPI_COMM_RT, ierr)

        end if

        ! If alternative bcast
        if (MPID%nsend.gt.0.and.MPID%altbcast) then

          ! For each slave
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing and reset the buffers
            call MPI_WAIT(MPID%requestA(iproc,1), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,2), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,3), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,4), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,5), &
                          MPI_STATUS_IGNORE,ierr)
            ! ALI
            if (lALI) then
              call MPI_WAIT(MPID%requestA(iproc,6), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,7), &
                            MPI_STATUS_IGNORE,ierr)
            end if ! ALI

          end do ! Processors

        end if ! Domain decomposition

        ! If going out, but with no ALI
        if (goout.and..not.lALI.and.INPUT%ALI_FORCE) then

          ! Don't leave and activate ALI
          goout = .False.
          force_ALI = .True.

        end if

        ! If we are finished, everyone exits
        if (goout) exit

      end do ! Iterations

      ! If 1.5D, Master save MRC
      if (pid.eq.0.and.run_mode.eq.1.and.Input%keep_MRC) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(Input%folder)//'/MRC', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1001

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1011
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1011

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = 0e0

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1101

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      ! Control
      call control

      !
      ! Clean slave pointers
      !
2000  if (pid.ne.0) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
        deallocate(rLineO,rPhotO)
        nullify(rLineO,rPhotO)
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO,p_LO)
      end if

      return

1000  umsg = 'Error opening MRCI file'
      urou = 'solverI_alt'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing MRCI file'
      urou = 'solverI_alt'
      close(800)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening MRC file'
      urou = 'solverI_alt'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1011  umsg = 'Error seeking MRC file'
      urou = 'solverI_alt'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing MRC file'
      urou = 'solverI_alt'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine solverI_alt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem for intensity using the ALI method,
      !! with one CPU (serial).\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!            lload(logical): Bool that says if previous
      !!                            solution was read\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine solverI_serial(Atom,LTElines,Rho_old,Atmo,Cont, &
                                Frec,Red,Geom,MPID,Input,lload, &
                                Stokes,J00,J00S,J00C,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      logical:: lload
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00
      double precision, dimension(nxt,Rz0:Rz1):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P

      ! Local

      type(MRC_class):: MRC

      logical:: goout,gooutprd,AD,ADT,ADD,PRDl,doNG,laux,lALI,RIRAM
      logical:: lp_exu,force_ALI

      character(LEN=20):: iterS

      integer:: ios,iz0,iz1,diz,m,o,p,op,iil,iip
      integer:: itran,jtran,ftran,jftran,fftran
      integer:: iter,iterr,ith,iph,ia,iz,ifreq,if0,if1
      integer:: NG_dim,NG_entry,iterm,iJ,tid

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset
      double precision:: daux,mu_inv,dsm,dsp
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh, &
                                  giz0:giz1), target:: Stokes_n
      double precision, dimension(nxt,Rz0:Rz1):: J00_n
      double precision, dimension(nxt,Rz0:Rz1):: J00S_n
      double precision, dimension(nfreq,Rz0:Rz1):: J00C_n
      double precision, dimension(nfreq), target:: LO
      double precision, dimension(nxb,nxt,Rz0:Rz1):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1):: LambdaP

      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:), pointer:: rLineO,rLineP
      double precision, dimension(:), pointer:: rPhotO,rPhotP
      double precision, dimension(:,:), pointer:: data2O,data2P
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP, p_LO
      double precision, dimension(:), pointer:: p_exu


      ! Initialize converged flag
      goout = .False.

      ! Initialize force ALI
      force_ALI = .False.

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! Really go back to angle dependent
      ADT = .not.AV.and.AD

      ! If storing redistribution
      RIRAM = IRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AVI = .True.
        if (ADT) IRAM = .False.
      end if

      ! Store if it was PRD (first iteration always CRD)
      PRDl = PRD
      if(.not.lload) PRD = .False.


      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.

      ! If NG acceleration
      if (Input%NGI) then

        ! Initialize NG dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + Rnz*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRDl) then

          ! If we need Stokes
          if (ADD) then

            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            NG_dim = NG_dim + nfreq*RnZ

          end if ! ADD
        end if ! PRD

        ! Check if it requires too much buffer
        if (dble(NG_dim)*8d-6.gt.maxbuffer_NG) then

          Input%NGI = .False.

          umsg = ' # The buffer for NG acceleration '// &
                 'is too big. Not doing NG.'
          call verbose

        end if

        ! If finally doing it, allocate
        if (Input%NGI) then

          allocate(NG_scratch(NG_dim, Input%NGI_ord+2))

        end if ! NG
      end if ! NG

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = 1
      if1 = nfreq

      !
      ! Allocations
      !

      ! Allocate O pointers
      allocate(data2O(Frec%ntfreqi,2))
      allocate(rLineO(Frec%ntfreqi))
      allocate(rPhotO(Frec%npfreq))

      ! Allocate M and O pointers for RT coeff
      allocate(data1M(nfreq,0:2))
      allocate(data1O(nfreq,0:2))

      ! Initialize exponential is not in memory
      if (PIRAM.and.Frec%pif1.ge.Frec%pif0) then
        lp_exu = .True.
      ! No exu allocated or needed
      else
        lp_exu = .False.
        allocate(p_exu(1))
      end if


      !
      ! Initialization messages
      !

      ! Only global Máster
      if (gpid.eq.0) then

        ! Announce we are starting
        umsg = '    Iteration          MRC(rho^0_0) Atom_index '// &
               'Level_index Height_index Height(km)'
        call verbose

      end if

      ! Only global Máster
      if (gpid.eq.0.and.Input%keep_MRC) then

        ! Open the file to store MRC
        ! If appending
        if (Input%appendMRC) then
          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRCI', exist=laux)
          if (.not.laux) then
            open(800, file=trim(Input%folder)//'/MRCI', &
                 action='write',iostat=ios,err=1000)
            write(800,'(A)') '!   Iteration          MRC(rho^0_0) '//&
                         'Atom_index Level_index Height_index '// &
                         'Height(km)'
            close(800)
          end if
        ! If not appending, create
        else
          open(800, file=trim(Input%folder)//'/MRCI', &
               action='write',iostat=ios,err=1000)
          write(800,'(A)') '!   Iteration          MRC(rho^0_0) '// &
                       'Atom_index Level_index Height_index '// &
                       'Height(km)'
          close(800)
        end if
      end if

      ! Initialize stimulated radiation tensor, becuase it is
      ! commented in the main loop
      J00S_n = 0d0

      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=Input%iteri_min,Input%iteri_max

        ! Flags for physics
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics
        if (iter.le.Input%allownphys_pop) then
          if (.not.nphysR) nphysR = .True.
        else
          if (nphysR) nphysR = .False.
        end if

        !
        ! The old population will be the current one
        !

        ! For each atom
        do ia=1,nA
          Rho_old(ia)%crho = Atom(ia)%crho
        end do

        ! Flag for ALI
        lALI = iter.gt.Input%ALI_delay.or.force_ALI

        ! Internal PRD iterations
        do iterr=1,Input%iteri_prd

          ! Reset radiation field variables
!$omp parallel default(none) &
!$omp shared(J00_n,J00P,J00C_n,LambdaL,LambdaP,lALI)
!$omp workshare
          J00_n = 0d0
         !J00S_n = 0d0
          J00P = 0d0
          J00C_n = 0d0
!$omp end workshare
          if (lALI) then
!$omp workshare
            LambdaL = 0d0
            LambdaP = 0d0
!$omp end workshare
          end if
!$omp end parallel


          !
          ! Ratiation Transfer
          !

          !  For each polar direction
          do ith=1,Geom%nTh

            ! Calculate inverse of cosine of polar direction
            mu_inv = 1d0/Geom%V_mu(ith)

            ! Determine the direction of propagation for indexes
            diz = -int(sign(1d0, Geom%V_mu(ith)))

            ! Determine the first and last height indexes to run over
            iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
            iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

            ! For each azimuthal direction
            do iph=1,Geom%nPh

              !
              ! If angle-dependent, manage scattering angles
              if (.not.AVI.and.PRD) &
                call scattering_manage(Geom,ith,iph)


              !
              ! First height
              !

              ! If going down, get top boundary
              if(diz.eq.1)then

                ! Call top boundary
                call topI(MPID,data1M(:,2))

              ! If going up, get bottom boundary
              else

                ! Call bottom boundary
                call bottomI(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                             Atmo%vy(iz0),Atmo%vz(iz0), &
                             Geom%V_mu(ith),Geom%V_mux(iph), &
                             Geom%V_muy(iph),MPID,data1M(:,2))

              endif

              ! Identify current height
              o = iz0

              ! Index for Stokes
              if (PRDl.and.ADD) op = o

              ! Calculate radiative coefficients
              call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                            o,ith,iph,if0,if1,J00(:,o), &
                            J00C(:,o),Cont%ndir, &
                            Cont%c(:,:,:,o),Stokes(:,:,:,op), &
                            rLineO,rPhotO,data1M(:,0:1),data2O,.True.)
              if (laborted) goto 2000

              if (KSTK) Stokes_n(:,iph,ith,o) = data1M(:,2)

              ! The initial point does not contribute to Lambda
              ! operator
              if (lALI) then
                rLineO = 0d0
                rPhotO = 0d0
              end if

              ! Point to exu values
              if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

              !
              ! Calculate integrals
              !
              call JcalcI(Atom,Geom,Frec%omega,Frec%W_freq, &
                          Frec%pif0,Frec%pif1, &
                          Atmo%T(o),Atmo%ne(o),iph,ith, &
                          data1M(:,2),rLineO,rPhotO,data2O, &
                          J00_n(:,o),J00S_n(:,o),J00P(:,:,o), &
                          J00C_n(:,o),LambdaL(:,:,o), &
                          LambdaP(:,:,:,o),lALI,p_exu)

              ! Identify next height
              p = iz0 + diz

              ! Index for Stokes
              if (PRDl.and.ADD) op = p

              ! Calculate radiative coefficients
              call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                            p,ith,iph,if0,if1,J00(:,p), &
                            J00C(:,p),Cont%ndir, &
                            Cont%c(:,:,:,p),Stokes(:,:,:,op), &
                            rLineO,rPhotO,data1O(:,0:1),data2O,.True.)
              if (laborted) goto 2000

              !
              ! Intermedium heights
              !

              ! For each height this CPU has assigned
              do iz=iz0,iz1,diz

                ! We treat the boundaries outside
                if(iz.eq.iz0.or.iz.eq.iz1)cycle

                ! Allocate P pointers
                allocate(data1P(nfreq,0:2))
                allocate(data2P(Frec%ntfreqi,2))
                allocate(rLineP(Frec%ntfreqi))
                allocate(rPhotP(Frec%npfreq))

                ! Identify heights
                m = iz - diz
                o = iz
                p = iz + diz

                ! Calculate distance to previous point
                dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                ! Caculate quantities of the next point
                dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

                ! If tau scale
                if (ztau) then
                  dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(m))
                  dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(p))
                end if

                ! Index for Stokes
                if (PRDl.and.ADD) op = p

                ! RT coefficients
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              p,ith,iph,if0,if1,J00(:,p), &
                              J00C(:,p),Cont%ndir, &
                              Cont%c(:,:,:,p),Stokes(:,:,:,op), &
                              rLineP,rPhotP, &
                              data1P(:,0:1),data2P,.True.)

                ! Point to the data
                p_K0M  => data1M(:,0)
                p_SM   => data1M(:,1)
                p_StkM => data1M(:,2)
                p_K0O  => data1O(:,0)
                p_SO   => data1O(:,1)
                p_StkO => data1O(:,2)
                p_K0P  => data1P(:,0)
                p_SP   => data1P(:,1)
                p_LO => LO(:)

                ! Apply short characteristics BESSER
                call RTStepI(o,ith,iph,nfreq, &
                             dsm,dsp,p_K0M,p_SM,p_K0O, &
                             p_SO,p_K0P,p_SP,p_StkM, &
                             p_StkO,p_LO,lALI,.True.)

                !
                ! Combine the value of lambda operator with the
                ! transition strength
                !
                if (lALI) then

                  ! Initialize indexes for rLine and rPhot
                  iil = 0
                  iip = 0

                  ! For each atom
                  do ia=1,nA

                    ! For each b-b trantision
                    do itran=1,Atom(ia)%ntran

                      ! For each FS transition
                      do ftran=1,Atom(ia)%fst(itran)%nt

                        ! Get the sequential index of FS transition
                        fftran = Atom(ia)%ifst_ij(ftran,itran)

                        ! Apply atomic shift
                        jftran = fftran + Atom(ia)%tfshift

                        ! For each frequency
                        do ifreq=Atom(ia)%if0(itran), &
                                 Atom(ia)%if1(itran)

                          iil = iil + 1
                          rLineO(iil) = LO(ifreq)*rLineO(iil)

                        end do ! frequency
                      end do ! FS transition
                    end do ! b-b transition

                    ! For each b-f transition
                    do itran=1,Atom(ia)%nphot

                      ! Apply atomic shift
                      jtran = itran + Atom(ia)%pshift

                      ! For each frequency
                      do ifreq=Atom(ia)%phot(itran)%if0, &
                               Atom(ia)%phot(itran)%if1

                        iip = iip + 1
                        rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                      end do ! frequency
                    end do ! b-f transition
                  end do ! atom

                end if ! ALI


                if (KSTK) Stokes_n(:,iph,ith,o) = data1O(:,2)

                ! Point to exu values
                if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

                !
                ! Calculate integrals
                !
                call JcalcI(Atom,Geom,Frec%omega,Frec%W_freq, &
                            Frec%pif0,Frec%pif1, &
                            Atmo%T(o),Atmo%ne(o),iph,ith, &
                            data1O(:,2),rLineO,rPhotO,data2O, &
                            J00_n(:,o),J00S_n(:,o), &
                            J00P(:,:,o),J00C_n(:,o), &
                            LambdaL(:,:,o),LambdaP(:,:,:,o), &
                            lALI,p_exu)

                ! Shift data (O->M, P->O)
                deallocate(data1M,data2O)
                data1M => data1O
                data1O => data1P
                data2O => data2P
                nullify(data1P,data2P)
                if (lALI) then
                  deallocate(rLineO,rPhotO)
                  rLineO => rLineP
                  rPhotO => rPhotP
                else
                  deallocate(rLineP,rPhotP)
                end if
                nullify(rLineP,rPhotP)

                ! Error
                if (laborted) goto 2000

              end do

              !
              ! Last height
              !

              ! Identify heights
              m = iz1 - diz
              o = iz1

              ! Calculate distance to previous point
              dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

              ! If tau scale
              if (ztau) &
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))

              ! Point to the data
              p_K0M  => data1M(:,0)
              p_SM   => data1M(:,1)
              p_StkM => data1M(:,2)
              p_K0O  => data1O(:,0)
              p_SO   => data1O(:,1)
              p_StkO => data1O(:,2)
              p_LO => LO(:)

              ! Apply short characteristics LINEAR
              call RTStepI(o,ith,iph,MPID%nf(pid), &
                           dsm,dsp,p_K0M,p_SM,p_K0O, &
                           p_SO,p_K0P,p_SP,p_StkM, &
                           p_StkO,p_LO,lALI,.False.)
              if (laborted) goto 2000

              !
              ! Combine the value of lambda operator with the
              ! transition strength
              !
              if (lALI) then

                ! Initialize indexes for rLine and rPhot
                iil = 0
                iip = 0

                ! For each atom
                do ia=1,nA

                  ! For each b-b trantision
                  do itran=1,Atom(ia)%ntran

                    ! For each FS transition
                    do ftran=1,Atom(ia)%fst(itran)%nt

                      ! Get the sequential index of FS transition
                      fftran = Atom(ia)%ifst_ij(ftran,itran)

                      ! Apply atomic shift
                      jftran = fftran + Atom(ia)%tfshift

                      ! For each frequency
                      do ifreq=Atom(ia)%if0(itran), &
                               Atom(ia)%if1(itran)

                        iil = iil + 1
                        rLineO(iil) = LO(ifreq)*rLineO(iil)

                      end do ! frequency
                    end do ! FS transition
                  end do ! b-b transition

                  ! For each b-f transition
                  do itran=1,Atom(ia)%nphot

                    ! Apply atomic shift
                    jtran = itran + Atom(ia)%pshift

                    ! For each frequency
                    do ifreq=Atom(ia)%phot(itran)%if0, &
                             Atom(ia)%phot(itran)%if1

                      iip = iip + 1
                      rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                    end do ! frequency
                  end do ! b-f transition
                end do ! atom

              end if ! ALI

              if (KSTK.or.o.eq.1) &
                Stokes_n(:,iph,ith,o) = data1O(:,2)

              ! Point to exu values
              if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)


              !
              ! Calculate integrals
              !
              call JcalcI(Atom,Geom,Frec%omega,Frec%W_freq, &
                          Frec%pif0,Frec%pif1, &
                          Atmo%T(o),Atmo%ne(o),iph,ith, &
                          data1O(:,2),rLineO,rPhotO,data2O, &
                          J00_n(:,o),J00S_n(:,o),J00P(:,:,o), &
                          J00C_n(:,o),LambdaL(:,:,o), &
                          LambdaP(:,:,:,o),lALI,p_exu)

            enddo ! azimuthal directions
          enddo ! polar directions


          ! Calculate MRC for J if PRD
          if (PRD.and.Input%iteri_prd.gt.1) then

            ! Call the routine
            call MRCJ_sb(J00C_n,J00C,MRC)

            ! Convert cm into km
            MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

            ! Only global Máster
            if (gpid.eq.0) then

              if (iterr.eq.1.and.iter.eq.2) then
                umsg = '         PRD            MRC(J^0_0) '// &
                       'Freq_index  Wavelength '// &
                       'Height_index Height(km)'
                call verbose
              end if

              ! Write in stdout
              write(umsg,'(2x,"PRD it:",1x,i3,2x,es20.12,'// &
                         '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                         iterr,MRC%values(2,1),MRC%indexes(1,1), &
                         1d2/Frec%omega(MRC%indexes(1,1)), &
                         MRC%indexes(2,1),MRC%values(1,1)
              call verbose

            end if

            ! Check exit criteria
            if (MRC%values(2,1).le.Input%mrci_r.or. &
               iterr.eq.Input%iteri_prd) then
              gooutprd = .True.
            else
              gooutprd = .False.
            end if

          ! Always go out if no PRD
          else

            gooutprd = .True.

          end if

          ! Shift the new values into the proper variables
!$omp parallel default(none) &
!$omp shared(KSTK,Stokes,Stokes_n,J00,J00_n,J00C,J00C_n)
          if (KSTK) then
!$omp workshare
            Stokes = Stokes_n
!$omp end workshare
          end if
!$omp workshare
          J00 = J00_n
         !J00S = J00S_n
          J00C = J00C_n
!$omp end workshare
!$omp end parallel

          ! If we are finished, everyone exits
          if (gooutprd) exit

        end do ! PRD iterations

        ! Control
        if (laborted) goto 2000


        !
        ! Solve SEE
        !
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P, &
                                  Input%folder,iter)
#endif
#ifdef DEBUGLAMBDA
      if (pid.eq.0) call dump_lambda(Atom,LambdaL,LambdaP, &
                                     Input%folder,iter)
#endif

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,fftran,jftran,urou,umsg,tid) &
!$omp shared(nA,nz,Atom,Rho_old,J00_n,J00P,LambdaL,LambdaP,lALI) &
!$omp shared(laborted,vaborted,Rz0,Rz1)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00_n(itran:jtran,iz), &
           !          J00S_n(itran:jtran,iz), &
                      J00_n(itran:jtran,iz), &
                      J00P(fftran:jftran,:,iz), &
                      LambdaL(:,itran:jtran,iz), &
                      LambdaP(:,fftran:jftran,:,iz),iz,lALI,tid,INPUT)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00_n(itran:jtran,iz), &
           !          J00S_n(itran:jtran,iz), &
                      J00_n(itran:jtran,iz), &
                      J00P(fftran:jftran,:,iz), &
                      LambdaL(:,itran:jtran,iz), &
                      LambdaP(:,fftran:jftran,:,iz),iz,lALI,tid)
#endif

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel
#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,iter)
#endif


        !
        ! NG acceleration
        !

        ! Check if doing NG acceleration
        if(Input%NGI.and.iter.gt.Input%NGI_delay)then

          ! Advance ntry
          NG_entry = NG_entry + 1

          ! Initialize index
          o = 0

          ! For each atom
          do ia=1,NA

            ! For each term
            do iterm=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(iterm)

                ! Get level index
                p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                ! For each height
                do iz=Rz0,Rz1

                  o = o + 1

                  ! Store rho
                  NG_scratch(o,NG_entry) = dble(Atom(ia)%crho(p,iz))

                enddo ! Heights
              enddo ! Levels
            enddo ! Term
          enddo ! Atoms

          ! If PRD
          if (PRDl) then

            ! If ADD
            if (ADD) then

              ! For each height
              do iz=giz0,giz1

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! For each azimuthal direction
                  do iph=1,Geom%nPh

                    ! For each frequency
                    do ifreq=1,nfreq

                      ! Advance index
                      o = o + 1

                      ! Store Stokes
                      NG_scratch(o,NG_entry) = &
                                            Stokes(ifreq,iph,ith,iz)

                    end do ! Frequencies
                  end do ! Azimuthal directions
                end do ! Polar directions
              end do ! Heights

            ! If AV
            else

              ! For each height
              do iz=Rz0,Rz1

                ! For each frequency
                do ifreq=1,nfreq

                  ! Advance index
                  o = o + 1

                  ! Store Stokes
                  NG_scratch(o,NG_entry) = J00C(ifreq,iz)

                end do ! Frequencies
              end do ! Heights

            end if ! AV or AD or dyn
          end if ! PRD

          ! Call NG and check if it should be processed
          call NG(NG_dim,Input%NGI_ord,NG_scratch,NG_entry,doNG)

          ! If applying NG
          if(doNG)then

            ! Initialize index
            o = 0

            ! For each atom
            do ia=1,NA

              ! For each term
              do iterm=1,Atom(ia)%nMulti

                ! For each level
                do iJ=1,Atom(ia)%nJ(iterm)

                  ! Get level index
                  daux = sqrt(2d0*Atom(ia)%rJval(iJ,iterm) + 1d0)
                  m = Atom(ia)%irho(iterm)%irho_ij(iJ)
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Advance index
                    o = o + 1

                    ! Accelerate rho
                    Atom(ia)%crho(p,iz) = &
                                   dcmplx(NG_scratch(o,NG_entry), 0d0)
                    Atom(ia)%popu(m,iz) = NG_scratch(o,NG_entry)*daux

                  enddo ! Heights
                enddo ! Levels
              enddo ! Terms
            enddo ! Atoms

            ! If PRD
            if (PRDl) then

              ! If ADD
              if (ADD) then

                ! For each height
                do iz=giz0,giz1

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Advance index
                        o = o + 1

                        ! Store Stokes
                        Stokes(ifreq,iph,ith,iz) = &
                                               NG_scratch(o,NG_entry)

                      end do ! Frequencies
                    end do ! Azimuthal directions
                  end do ! Polar directions
                end do ! Heights

              ! If AV
              else

                ! For each height
                do iz=Rz0,Rz1

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! Advance index
                    o = o + 1

                    ! Store Stokes
                    J00C(ifreq,iz) = NG_scratch(o,NG_entry)

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD or dyn
            end if ! PRD

            ! Send message if global Máster
            if (gpid.eq.0) then
              umsg = 'NG acceleration'
              call verbose
            end if
            ! Reset entry
            NG_entry = 0

          end if ! Applying NG in this iteration
        endif ! Doing NG acceleration


        !
        !  Calculate MRC
        !

        call MRCI_sb(Atom,Rho_old,MRC)

        ! Convert cm into km
        MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

        ! Global Máster outputs
        if (gpid.eq.0) then

          ! Write in stdout
          write(umsg,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                     '2x,f9.3)') &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(2,1),MRC%values(1,1)
          call verbose

          ! File
          if (Input%keep_MRC) then

            ! Write in MRC file
            open(800, file=trim(Input%folder)//'/MRCI', &
                 iostat=ios,err=1000,position='append')
            write(800,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                       '2x,f9.3)', err=1100) &
            iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
            MRC%indexes(2,1),MRC%values(1,1)
            close(800)

          end if
        end if

        ! Check exit criteria
        if (MRC%values(2,1).le.Input%mrci_i) goout = .True.

        ! Check in first iteration if was a fix population case
        ! but we want PRD
        if (iter.eq.Input%iteri_min.and.PRDl.and. &
          MRC%values(2,1).lt.1d-16) goout = .False.

        ! If in the mandatory non-PRD
        if (iter.le.Input%PRD_delay) goout = .False.

        ! If going to AD from AA
        if (iter.eq.Input%iteri_min.and.tbAD.and.ADT) goout = .False.

        ! We can swith now to AD if we had AV input
        if (tbAD.and.ADT) then
          AVI = .False.
          tbAD = .False.
          IRAM = RIRAM
        end if

        ! Recover PRD variable
        if (iter.eq.Input%PRD_delay) PRD = PRDl

        ! If going out, but with no ALI
        if (goout.and..not.lALI.and.INPUT%ALI_FORCE) then

          ! Don't leave and activate ALI
          goout = .False.
          force_ALI = .True.

        end if

        ! Save data of this steps if proceeds
        if(Input%storei.and.mod(iter,Input%storei_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesolI(Input,iterS,Frec%omega,Geom, &
         !               Atom,Atmo%z,Stokes,J00,J00S,J00C,J00P, &
                         Atom,Atmo%z,Stokes,J00,J00,J00C,J00P, &
                         .False.)

        end if

        ! And update Stokes if necessary
        if (KSTK) Stokes = Stokes_n

        ! Control
        if (laborted) goto 2000

        ! If we are finished, everyone exits
        if (goout) exit

      end do ! Iterations

      ! If 1.5D, Master save MRC
      if (run_mode.eq.1.and.Input%keep_MRC) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(Input%folder)//'/MRC', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1001

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1011
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1011

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = 0e0

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1101

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      !
      ! Clean pointers
      !
2000  deallocate(data1M,data1O,data2O)
      nullify(data1M,data1O,data2O)
      deallocate(rLineO,rPhotO)
      nullify(rLineO,rPhotO)
      nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
      nullify(p_StkM,p_StkO,p_LO)
      if (lp_exu) then
        nullify(p_exu)
      else
        deallocate(p_exu)
        nullify(p_exu)
      end if

      return

1000  umsg = 'Error opening MRCI file'
      urou = 'solverI_serial'
      call aborted
      return
1100  umsg = 'Error writing MRCI file'
      urou = 'solverI_serial'
      close(800)
      call aborted
      return
1001  umsg = 'Error opening MRC file'
      urou = 'solverI_serial'
      call aborted
      call control
      return
1011  umsg = 'Error seeking MRC file'
      urou = 'solverI_serial'
      call MPI_FILE_CLOSE(funit, ierr)
      call aborted
      call control
      return
1101  umsg = 'Error writing MRC file'
      urou = 'solverI_serial'
      call MPI_FILE_CLOSE(funit, ierr)
      call aborted
      call control
      return

      end subroutine solverI_serial

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the emergent intensity for the specified LOS\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      subroutine emergentI(Atom,LTElines,Atmo,Cont,Frec,Red,Geom, &
                           MPID,Input,Stokes,J00,J00C,SolF)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Solution_F_class):: SolF
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00

      if (MPID%mpi) then
        call emergenceI(Atom,LTElines,Atmo,Cont,Frec,Red,Geom, &
                        MPID,Input,Stokes,J00,J00C,SolF)
      else
        call emergenceI_serial(Atom,LTElines,Atmo,Cont,Frec,Red, &
                               Geom,MPID,Input,Stokes, &
                               J00,J00C,SolF)
      end if

      end subroutine emergentI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emergent intensity for the specified LOS, with
      !! several CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      subroutine emergenceI(Atom,LTElines,Atmo,Cont,Frec,Red,Geom, &
                            MPID,Input,Stokes,J00,J00C,SolF)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Solution_F_class):: SolF
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00

      ! Local

      logical:: AD, ADD

      integer:: iaux,icount,ncount
      integer:: ith,iph,iz,if0,if1,id
      integer:: iz0,iz1,diz,m,o,op,p
      integer:: tau1size,sshift

      double precision:: mu_inv,dsm,dsp,dzm,dzp
      double precision, dimension(1):: daux1
      double precision, dimension(:), allocatable:: Contr
      double precision, dimension(:,:), allocatable, target:: tau
      double precision, dimension(:,:), allocatable, target:: tau1
      double precision, dimension(:,:), allocatable:: ContrG

      ! Buffers

      ! Receivers
      integer:: rpid
      double precision, dimension(:), allocatable:: Stokes_r
      double precision, dimension(:), allocatable:: Contr_r
      double precision, dimension(nfreq):: Stokes_out
      ! Senders
      double precision, dimension(:), allocatable:: Stokes_s
      double precision, dimension(:,:), allocatable:: Contr_s
      double precision, dimension(:,:), allocatable:: tau1_s
      ! Dual
      integer, dimension(2):: info_b

      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:), pointer:: etaIM,etaIO
      double precision, dimension(:), pointer:: tauM
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP


      !
      ! Initializations
      !


      ! Reset progress counter
      icount = 0

      ! Determine number of directions to do
      ncount = Geom%nThLOS*Geom%nPhLOS

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)


      !
      ! Allocation of buffers
      !

      ! Master
      if(pid.eq.0)then

        ! To receive Intensity chunks
        iaux = MPID%nxfreq
        allocate(Stokes_r(iaux))

        ! If calculating height of tau=1, allocate
        if (Input%out_tau1) allocate(tau1(2,nfreq))

        ! If calculating contribution function, allocate
        if (Input%out_contr) then
          iaux = MPID%nxfreq*Rnz
          allocate(Contr_r(iaux))
          allocate(ContrG(nfreq,Rz0:Rz1))
        end if

        ! If inverting, need to return the output
        if (run_mode.eq.-1) then
          allocate(SolF%e_Stk(0:0,nfreq,Geom%nPhLOS,Geom%nThLOS))
          if (Input%out_tau1) &
          allocate(SolF%e_tau1(Input%lim_tau%nn,Geom%nPhLOS, &
                               Geom%nThLOS))
          if (Input%out_contr) &
          allocate(SolF%e_Ctr(0:0,Input%lim_ctr%nn,Rz0:Rz1, &
                              Geom%nPhLOS,Geom%nThLOS))
        end if ! Inversion

      ! Slave
      else

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(if0:if1,0:2))
        allocate(data1O(if0:if1,0:2))

        ! To send Intensity chunks
        allocate(Stokes_s(MPID%if0(pid):MPID%if1(pid)))

        ! If calculating tau 1 or contribution function, allocate
        if (input%out_tau1.or.input%out_contr) then

          allocate(tau(if0:if1,Rz0:Rz1))
          allocate(tau1(2,if0:if1))
          allocate(tau1_s(2,if0:if1))
          tau1size = MPID%nf(pid)*2
          allocate(etaIM(if0:if1))

          ! If calculating contribution function, allocate
          if (input%out_contr) then

            allocate(Contr(1:MPID%nf(pid)))
            allocate(Contr_s(1:MPID%nf(pid),Rz0:Rz1))

          end if ! contribution output

        else
          nullify(etaIM,etaIO,tauM)
        end if ! tau1 output

      end if ! master or slave

      ! Control
      call control
      if (laborted) goto 2000


      !
      ! Formal solutions
      !

      ! For each polar LOS direction
      do ith=1,Geom%nThLOS

        ! Calculate inverse of cosine of polar direction
        mu_inv = 1d0/Geom%L_mu(ith)

        ! Determine the direction of propagation for indexes
        diz = -int(sign(1d0, Geom%L_mu(ith)))

        !
        ! Determine the boundaries and the domain decomposition
        ! partnerts
        !

        ! If slave
        if (pid.gt.0) then

          ! Determine the first and last height indexes to run over
          iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
          iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

        end if

        ! For each azimuthal LOS direction
        do iph=1,Geom%nPhLOS

          !
          ! Master
          !
          if (pid.eq.0) then

            ! Advance the counter
            icount = icount + 1

            ! Write progress into stdout if global master
            if (gpid.eq.0) then
              write(umsg,'(A,i4,A,i4)') &
                         '   Doing direction ',icount,' of ',ncount
              call verbose
            end if


            !
            ! If calculating height of tau=1
            !
            if (Input%out_tau1) then

              do id=1,MPID%nnd

                call MPI_recv(rpid,1,MPI_INTEGER, &
                              MPI_ANY_SOURCE, 2000000, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

                tau1size = MPID%nf(rpid)*2

                call MPI_recv(tau1(1,MPID%if0(rpid)), &
                              tau1size, &
                              MPI_DOUBLE_PRECISION, rpid, &
                              3000000+rpid, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)

              end do

              ! Inversion
              if (run_mode.eq.-1) then

                ! Keep tau1
                call settau(SolF%e_tau1(:,iph,ith),tau1(2,:), &
                            Input%lim_tau)

              ! Synthesis
              else

                ! Store the height where tau=1
                call writetau(Input%folder,iph,ith,Frec%omega,Geom, &
                              tau1(2,:),Input%lim_tau)
                call control
                if (laborted) goto 2000

              end if ! Inversion
            end if ! calculate tau=1

            ! For each frequency domain
            do id=1,MPID%nnd

              !
              ! Receive data from a slave
              !

              ! Receive indexing data
              call MPI_recv(info_b,2,MPI_INTEGER, &
                            MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)

              ! If calculating contribution function
              if (Input%out_contr) then

                ! Receive contribution function data
                call MPI_recv(Contr_r(1), MPID%sizei14(info_b(1)), &
                              MPI_DOUBLE_PRECISION, info_b(1), &
                              1000000+info_b(1), MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)

                ! Reset shift in index
                sshift = 0

                ! For each height
                do iz=Rz0,Rz1

                  ! Rearrange the contribution function
                  ContrG(MPID%if0(info_b(1)):MPID%if1(info_b(1)), &
                         iz) = Contr_r(sshift+1: &
                                       sshift+MPID%nf(info_b(1)))

                  ! Update the shift in the buffer
                  sshift = sshift + MPID%nf(info_b(1))

                end do ! heights

              end if ! If contribution function

              ! If it is not a boundary, we are not receiving
              ! anything else
              if (info_b(2).lt.0) cycle

              ! Receive intensity
              call MPI_recv(Stokes_r(1), MPID%nf(info_b(1)), &
                            MPI_DOUBLE_PRECISION, info_b(1), &
                            info_b(1), MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)

              ! Rearrange the intensity
              Stokes_out(MPID%if0(info_b(1)):MPID%if1(info_b(1))) = &
                                       Stokes_r(1:MPID%nf(info_b(1)))

            end do ! Frequency domains


          !
          ! Slave
          !
          else

            !
            ! Get geometry if PRD AD
            !
            if (PRD.and..not.AVI) &
              call get_scattering_los(Geom,ith,iph)


            !
            ! If calculating height of tau=1
            !
            if (Input%out_tau1.or.Input%out_contr) then

              ! Reset cummulative quantities
              tau = 0d0
              tau1(1,:) = 0
              tau1(2,:) = Atmo%z(Rz0)
              tauM => tau1(1,if0:if1)

              ! First height

              ! Top boundary
              o = Rz0

              ! Calculate absorptivity
              call RTAbsI(Frec,Atom,LTElines,Atmo,Geom, &
                          o,ith,iph,if0,if1,Cont%ndir, &
                          Cont%c(:,:,:,o),etaIM)

              ! Middle heights
              do iz=Rz0,Rz1

                if (iz.eq.Rz0.or.iz.eq.Rz1) cycle

                ! Allocate O pointer
                allocate(etaIO(if0:if1))

                ! Determine indexes
                m = iz - 1
                o = iz

                ! Calculate distance to previous point
                dsm = abs(Atmo%z(o) - Atmo%z(m))*mu_inv

                ! If tau scale
                if (ztau) &
                  dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(m))

                ! Calculate opacity current point
                call RTAbsI(Frec,Atom,LTElines,Atmo,Geom, &
                            o,ith,iph,if0,if1,Cont%ndir, &
                            Cont%c(:,:,:,o),etaIO)

                ! Accumulate tau
                call RTtauI(dsm,MPID%nf(pid),Atmo%z(m),Atmo%z(o), &
                            etaIM(if0:if1),etaIO(if0:if1), &
                            tauM(1:MPID%nf(pid)),tau(if0:if1,o), &
                            tau1(:,if0:if1))


                ! Shift the opacity value and previous tau
                deallocate(etaIM)
                etaIM => etaIO
                nullify(etaIO)
                tauM => tau(if0:if1,o)

              end do

              ! Allocate O pointer
              allocate(etaIO(if0:if1))

              ! Identify heights
              m = Rz1 - 1
              o = Rz1

              ! Calculate distance to previous point
              dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

              ! If tau scale
              if (ztau) &
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))

              ! Calculate opacity current point
              call RTAbsI(Frec,Atom,LTElines,Atmo,Geom, &
                          o,ith,iph,if0,if1,Cont%ndir, &
                          Cont%c(:,:,:,o),etaIO)

              ! Accumulate tau
              call RTtauI(dsm,MPID%nf(pid),Atmo%z(m),Atmo%z(o), &
                          etaIM(if0:if1),etaIO(if0:if1), &
                          tauM(1:MPID%nf(pid)), &
                          tau(if0:if1,o),tau1(:,if0:if1))

              ! Deallocate etaI pointers
              deallocate(etaIO)
              nullify(etaIO)

              !
              ! Send to master the tau=1
              !

              ! If sending tau
              if (Input%out_tau1) then

                ! Wait for last send to finish
                call MPI_WAIT(MPID%request7,MPI_STATUS_IGNORE, &
                              ierr)
                call MPI_WAIT(MPID%request8,MPI_STATUS_IGNORE, &
                              ierr)

                ! Send indexes
                call MPI_ISEND(pid,1,MPI_INTEGER, &
                               0,2000000,MPI_COMM_RT, &
                               MPID%request7,ierr)


                ! Send Intensity
                tau1_s = tau1
                call MPI_ISEND(tau1_s(1,if0), tau1size, &
                               MPI_DOUBLE_PRECISION, &
                               0,3000000+pid,MPI_COMM_RT, &
                               MPID%request8,ierr)

                ! If synthesis, control
                if (run_mode.ne.-1) then
                  call control
                  if (laborted) goto 2000
                end if

              end if
            end if ! calculate tau=1

            !
            ! Actual intensity emergence
            !

            !
            ! First height
            !

            ! If going down, get top boundary
            if(diz.eq.1)then

              ! Call top boundary
              call topI(MPID,data1M(:,2))

            ! If going up, get bottom boundary
            else

              ! Call bottom boundary
              call bottomI(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                           Atmo%vy(iz0),Atmo%vz(iz0), &
                           Geom%L_mu(ith),cos(Geom%L_phi(iph)), &
                           sin(Geom%L_phi(iph)),MPID,data1M(:,2))

            endif ! propagation direction

            ! Identify current height
            o = iz0

            ! Index for Stokes
            if (PRD.and.ADD) op = o

            ! Calculate radiative coefficients
            call RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom, &
                           o,ith,iph,if0,if1,J00(:,o),J00C(:,o), &
                           Cont%ndir,Cont%c(:,:,:,o), &
                           Stokes(:,:,:,op),data1M(:,0:1))

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Wait till last communication was received
              call MPI_WAIT(MPID%request5,MPI_STATUS_IGNORE,ierr)

              ! The first does not contribute
              contr_s(:,o) = 0d0

            end if ! computing contribution function

            ! Identify next height
            p = iz0 + diz

            ! Index for Stokes
            if (PRD.and.ADD) op = p

            ! Calculate radiative coefficients
            call RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom, &
                           p,ith,iph,if0,if1,J00(:,p),J00C(:,p), &
                           Cont%ndir,Cont%c(:,:,:,p), &
                           Stokes(:,:,:,op),data1O(:,0:1))


            !
            ! Intermedium heights
            !

            ! For each height this CPU has assigned
            do iz=iz0,iz1,diz

              ! We treat the boundaries outside
              if(iz.eq.iz0.or.iz.eq.iz1)cycle

              ! Allocate P pointers
              allocate(data1P(if0:if1,0:2))

              ! Identify heights
              m = iz - diz
              o = iz
              p = iz + diz

              ! Calculate distance to previous point
              dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

              ! Calculate distance to the next point
              dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

              ! If tau scale
              if (ztau) then
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))
                dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(p))
              end if

              ! Index for Stokes
              if (PRD.and.ADD) op = p

              ! Calculate radiative coefficients
              call RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom, &
                             p,ith,iph,if0,if1,J00(:,p),J00C(:,p), &
                             Cont%ndir,Cont%c(:,:,:,p), &
                             Stokes(:,:,:,op),data1P(:,0:1))

              ! Point to the data
              p_K0M  => data1M(:,0)
              p_SM   => data1M(:,1)
              p_StkM => data1M(:,2)
              p_K0O  => data1O(:,0)
              p_SO   => data1O(:,1)
              p_StkO => data1O(:,2)
              p_K0P  => data1P(:,0)
              p_SP   => data1P(:,1)

              ! Apply short characteristics BESSER
              call RTStepI(o,ith,iph,MPID%nf(pid), &
                           dsm,dsp,p_K0M,p_SM,p_K0O, &
                           p_SO,p_K0P,p_SP,p_StkM, &
                           p_StkO,daux1,.False.,.True.)

              ! If calculating contribution function
              if (Input%out_contr) then

                ! Calculate vertical distance between points
                if (ztau) then
                  dzm = dsm/mu_inv
                  dzp = dsp/mu_inv
                else
                  dzm = Atmo%z(o) - Atmo%z(m)
                  dzp = Atmo%z(p) - Atmo%z(o)
                end if

                ! Compute contribution
                call RTContrI(dsm,dsp,MPID%nf(pid),dzm,dzp, &
                              p_K0M,p_K0O,p_SO,p_K0P, &
                              tau(MPID%if0(pid):MPID%if1(pid),o), &
                              Contr,.True.)

                ! Store in buffer
                contr_s(:,o) = contr

              end if ! calculating contribution

              ! Shift data (O->M, P->O)
              deallocate(data1M)
              data1M => data1O
              data1O => data1P
              nullify(data1P)

            end do ! Intermedium heights


            !
            ! Last height
            !

            ! Identify heights
            m = iz1 - diz
            o = iz1

            ! Calculate distance to previous point
            dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

            ! If tau scale
            if (ztau) &
              dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(m))

            ! Point to the data
            p_K0M  => data1M(:,0)
            p_SM   => data1M(:,1)
            p_StkM => data1M(:,2)
            p_K0O  => data1O(:,0)
            p_SO   => data1O(:,1)
            p_StkO => data1O(:,2)

            ! Apply short characteristics LINEAR
            call RTStepI(o,ith,iph,MPID%nf(pid), &
                         dsm,dsp,p_K0M,p_SM,p_K0O, &
                         p_SO,p_K0P,p_SP,p_StkM, &
                         p_StkO,daux1,.False.,.False.)

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Calculate vertical distance between points
              if (ztau) then
                dzm = dsm/mu_inv
              else
                dzm = Atmo%z(o) - Atmo%z(m)
              end if

              ! Calculate contribution function
              call RTContrI(dsm,dsp,MPID%nf(pid),dzm,dzp,p_K0M, &
                            p_K0O, p_SO, p_K0P, &
                            tau(MPID%if0(pid):MPID%if1(pid),o), &
                            Contr,.False.)

            end if ! calculating contribution function

            ! Wait till last communication was received
            call MPI_WAIT(MPID%request3,MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%request4,MPI_STATUS_IGNORE,ierr)

            !
            ! Send to master
            !

            ! Send indexes
            info_b = (/ pid, 1 /)
            call MPI_ISEND(info_b(1),2,MPI_INTEGER, &
                           0,0, MPI_COMM_RT,MPID%request3, &
                           ierr)

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Send contribution function
              contr_s(:,o) = contr
              call MPI_ISEND(contr_s(1,Rz0), &
                             MPID%sizei14(pid), &
                             MPI_DOUBLE_PRECISION, &
                             0, 1000000+pid, MPI_COMM_RT, &
                             MPID%request5, ierr)
            end if

            ! Send intensity
            Stokes_s = data1O(:,2)
            call MPI_ISEND(Stokes_s(if0), MPID%nf(pid), &
                           MPI_DOUBLE_PRECISION, 0, pid, &
                           MPI_COMM_RT, MPID%request4, ierr)

          end if ! Master or slave

          ! Master
          if (pid.eq.0) then

            ! If inverting
            if (run_mode.eq.-1) then

              ! Keep Stokes
              call setstkI(SolF%e_Stk(:,:,iph,ith),Stokes_out, &
                           Input%lim_stk,.False.)

              ! Keep contribution function
              if (Input%out_contr) &
                call setctrI(SolF%e_Ctr(:,:,:,iph,ith),ContrG, &
                             Input%lim_ctr)

            ! Synthesis, to file
            else

              ! Write stokes
              call writestkI(Input%folder,iph,ith,Frec%omega,Geom, &
                             Stokes_out,Input%lim_stk)
              if (laborted) goto 2000

              ! Write contribution function
              if (Input%out_contr) then
                call writectrI(Input%folder,iph,ith,Frec%omega,Geom, &
                               Atmo%z,ContrG,Input%lim_ctr)
                if (laborted) goto 2000
              end if

              ! Call control
              call control
              if (laborted) goto 2000

            end if

            ! If stdout terminal, remove bar and say completed
            if (gpid.eq.0) then
              umsg = '   Completed'
              call verbose
            end if

          ! Slave in synthesis
          else if (run_mode.ne.-1) then

            ! Control
            call control
            if (laborted) goto 2000

          end if

        enddo ! azimuthal LOS directions
      enddo ! polar LOS directions

      ! Control
      call control

      !
      ! Clean
      !
2000  if (pid.ne.0) then
        deallocate(data1M,data1O)
        nullify(data1M,data1O)
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO)
        if (associated(tauM)) nullify(tauM)
        if (associated(etaIM)) then
          deallocate(etaIM)
          nullify(etaIM)
        end if
      end if

      return

      end subroutine emergenceI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emergent intensity for the specified LOS, with
      !! one CPU (serial).\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      subroutine emergenceI_serial(Atom,LTElines,Atmo,Cont,Frec,Red, &
                                   Geom,MPID,Input, &
                                   Stokes,J00,J00C,SolF)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Solution_F_class):: SolF
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C
      double precision, dimension(nxt,Rz0:Rz1):: J00

      ! Local

      logical:: AD,ADD

      integer:: icount,ncount,ith,iph,iz
      integer:: iz0,iz1,diz,m,o,op,p,if0,if1

      double precision:: mu_inv,dsm,dsp,dzm,dzp
      double precision, dimension(1):: daux1
      double precision, dimension(:,:), allocatable:: tau, tau1
      double precision, dimension(:,:), allocatable:: Contr

      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:), pointer:: etaIM,etaIO
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP


      !
      ! Initializations
      !


      ! Reset progress counter
      icount = 0

      ! Determine number of directions to do
      ncount = Geom%nThLOS*Geom%nPhLOS

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = 1
      if1 = nfreq

      ! Allocate M and O pointers for RT coeff
      allocate(data1M(nfreq,0:2))
      allocate(data1O(nfreq,0:2))

      ! If calculating height of tau=1, allocate
      if (Input%out_tau1.or.Input%out_contr) then
        allocate(tau(nfreq,Rz0:Rz1))
        allocate(etaIM(nfreq))
        allocate(tau1(2,nfreq))
        allocate(Contr(nfreq,Rz0:Rz1))
      else
        nullify(etaIM,etaIO)
      end if

      ! If inverting, need to return the output
      if (run_mode.eq.-1) then
        allocate(SolF%e_Stk(0:0,nfreq,Geom%nPhLOS,Geom%nThLOS))
        if (Input%out_tau1) &
        allocate(SolF%e_tau1(Input%lim_tau%nn, &
                             Geom%nPhLOS,Geom%nThLOS))
        if (Input%out_contr) &
        allocate(SolF%e_Ctr(0:0,Input%lim_ctr%nn,Rz0:Rz1, &
                            Geom%nPhLOS,Geom%nThLOS))
      end if ! Inversion


      !
      ! Formal solutions
      !

      ! For each polar LOS direction
      do ith=1,Geom%nThLOS

        ! Calculate inverse of cosine of polar direction
        mu_inv = 1d0/Geom%L_mu(ith)

        ! Determine the direction of propagation for indexes
        diz = -int(sign(1d0, Geom%L_mu(ith)))

        ! Determine the first and last height indexes to run over
        iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
        iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

        ! For each azimuthal LOS direction
        do iph=1,Geom%nPhLOS

          ! Advance the index
          icount = icount + 1

          ! Communicate which direction we are doing if global Máster
          if (gpid.eq.0) then
            write(umsg,'(A,i4,A,i4)') &
                        '   Doing direction ',icount,' of ',ncount
            call verbose
          end if

          !
          ! Get geometry if PRD AD
          !
          if (PRD.and..not.AVI) &
            call get_scattering_los(Geom,ith,iph)

          ! If calculating contribution function
          if (Input%out_tau1.or.Input%out_contr) then

            ! Reset tau
            tau = 0
            tau1(1,:) = 0
            tau1(2,:) = Atmo%z(Rz0)

            ! Initial point
            o = Rz0

            ! Calculate absorptivity
            call RTAbsI(Frec,Atom,LTElines,Atmo,Geom, &
                        o,ith,iph,if0,if1,Cont%ndir, &
                        Cont%c(:,:,:,o),etaIM)

            ! For each intermediate point
            do iz=Rz0+1,Rz1

              ! Allocate O pointer
              allocate(etaIO(nfreq))

              ! Determine indexes
              m = iz - 1
              o = iz

              ! Calculate distance to previous point
              dsm = abs(Atmo%z(o) - Atmo%z(m))*mu_inv

              ! If tau scale
              if (ztau) &
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))

              ! Calculate opacity current point
              call RTAbsI(Frec,Atom,LTElines,Atmo,Geom, &
                          o,ith,iph,if0,if1,Cont%ndir, &
                          Cont%c(:,:,:,o),etaIO)

              ! Calculate run of tau
              call RTtauI(dsm,nfreq,Atmo%z(m),Atmo%z(o),etaIM, &
                          etaIO,tau(:,m),tau(:,o),tau1)

              ! Shift the opacity variable
              deallocate(etaIM)
              etaIM => etaIO
              nullify(etaIO)

            end do ! heights

            ! If calculating height where tau=1
            if (Input%out_tau1) then

              ! Inversion
              if (run_mode.eq.-1) then

                ! Keep tau1
                call settau(SolF%e_tau1(:,iph,ith),tau1(2,:), &
                            Input%lim_tau)

              ! Synthesis
              else

                ! Write to file
                call writetau(Input%folder,iph,ith,Frec%omega,Geom, &
                              tau1(2,:),Input%lim_tau)

              end if ! Inversion
            end if ! Tau output
          end if ! If calculating contribution function or tau=1


          !
          ! First height
          !

          ! If going down, get top boundary
          if(diz.eq.1)then

            ! Call top boundary
            call topI(MPID,data1M(:,2))

          ! If going up, get bottom boundary
          else

            ! Call bottom boundary
            call bottomI(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                         Atmo%vy(iz0),Atmo%vz(iz0), &
                         Geom%L_mu(ith),cos(Geom%L_phi(iph)), &
                         sin(Geom%L_phi(iph)),MPID,data1M(:,2))

          endif ! propagation direction

          ! Identify current height
          o = iz0

          ! Index for Stokes
          if (PRD.and.ADD) op = o

          ! Calculate radiative coefficients
          call RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom, &
                         o,ith,iph,if0,if1,J00(:,o),J00C(:,o), &
                         Cont%ndir,Cont%c(:,:,:,o), &
                         Stokes(:,:,:,op),data1M(:,0:1))

          ! Identify next height
          p = iz0 + diz

          ! Index for Stokes
          if (PRD.and.ADD) op = p

          ! Calculate radiative coefficients
          call RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom, &
                         p,ith,iph,if0,if1,J00(:,p),J00C(:,p), &
                         Cont%ndir,Cont%c(:,:,:,p), &
                         Stokes(:,:,:,op),data1O(:,0:1))

          ! Contribution function at first point is 0
          if (Input%out_contr) Contr(:,o) = 0d0


          !
          ! Intermedium heights
          !

          ! Medium heights
          do iz=iz0,iz1,diz

            ! We treat the boundaries outside
            if(iz.eq.iz0.or.iz.eq.iz1)cycle

            ! Allocate P pointers
            allocate(data1P(nfreq,0:2))

            ! Identify heights
            m = iz - diz
            o = iz
            p = iz + diz

            ! Calculate distance to previous point
            dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

            ! Calculate distance to the next point
            dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

            ! If tau scale
            if (ztau) then
              dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(m))
              dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(p))
            end if

            ! Index for Stokes
            if (PRD.and.ADD) op = p

            ! Calculate radiative coefficients
            call RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom, &
                           p,ith,iph,if0,if1,J00(:,p),J00C(:,p), &
                           Cont%ndir,Cont%c(:,:,:,p), &
                           Stokes(:,:,:,op),data1P(:,0:1))


            ! Point to the data
            p_K0M  => data1M(:,0)
            p_SM   => data1M(:,1)
            p_StkM => data1M(:,2)
            p_K0O  => data1O(:,0)
            p_SO   => data1O(:,1)
            p_StkO => data1O(:,2)
            p_K0P  => data1P(:,0)
            p_SP   => data1P(:,1)

            ! Apply short characteristics BESSER
            call RTStepI(o,ith,iph,nfreq, &
                         dsm,dsp,p_K0M,p_SM,p_K0O, &
                         p_SO,p_K0P,p_SP,p_StkM, &
                         p_StkO,daux1,.False.,.True.)

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Calculate vertical distance between points
              if (ztau) then
                dzm = dsm/mu_inv
                dzp = dsp/mu_inv
              else
                dzm = Atmo%z(o) - Atmo%z(m)
                dzp = Atmo%z(p) - Atmo%z(o)
              end if

              ! Compute contribution
              call RTContrI(dsm,dsp,nfreq,dzm,dzp,p_K0M, &
                            p_K0O,p_SO,p_K0P,tau(:,o), &
                            Contr(:,o),.True.)

            end if ! calculating contribution

            ! Shift data (O->M, P->O)
            deallocate(data1M)
            data1M => data1O
            data1O => data1P
            nullify(data1P)

          end do ! Intermedium heights


          !
          ! Last height
          !

          ! Identify heights
          m = iz1 - diz
          o = iz1

          ! Calculate distance to previous point
          dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

          ! If tau scale
          if (ztau) &
            dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                Atmo%chi500(m))

          ! Point to the data
          p_K0M  => data1M(:,0)
          p_SM   => data1M(:,1)
          p_StkM => data1M(:,2)
          p_K0O  => data1O(:,0)
          p_SO   => data1O(:,1)
          p_StkO => data1O(:,2)

          ! Apply short characteristics LINEAR
          call RTStepI(o,ith,iph,nfreq, &
                       dsm,dsp,p_K0M,p_SM,p_K0O, &
                       p_SO,p_K0P,p_SP,p_StkM, &
                       p_StkO,daux1,.False.,.False.)

          ! If calculating contribution function
          if (Input%out_contr) then

            ! Calculate vertical distance between points
            dzm = Atmo%z(o) - Atmo%z(m)

            ! Calculate vertical distance between points
            if (ztau) then
              dzm = dsm/mu_inv
            else
              dzm = Atmo%z(o) - Atmo%z(m)
            end if

            ! Calculate contribution function
            call RTContrI(dsm,dsp,nfreq,dzm,dzp,p_K0M, &
                          p_K0O,p_SO,p_K0P,tau(:,o), &
                          Contr(:,o),.False.)

            ! If inverting
            if (run_mode.eq.-1) then

              ! Keep contribution function
              call setctrI(SolF%e_Ctr(:,:,:,iph,ith),Contr, &
                           Input%lim_ctr)

            ! Synthesis
            else

              ! Write contribution function
              call writectrI(Input%folder,iph,ith,Frec%omega,Geom, &
                             Atmo%z,Contr,Input%lim_ctr)
            end if

          end if ! calculating contribution function

          ! If inverting
          if (run_mode.eq.-1) then

            ! Keep Stokes
            call setstkI(SolF%e_Stk(:,:,iph,ith),p_StkO, &
                         Input%lim_stk,.False.)

          ! Synthesis
          else

            ! Write stokes
            call writestkI(Input%folder,iph,ith,Frec%omega,Geom, &
                           p_StkO,Input%lim_stk)

          end if

          ! Communicate we finished this direction
          if (gpid.eq.0) then
            umsg = '   Completed'
            call verbose
          end if

        enddo ! azimuthal LOS directions
      enddo ! polar LOS directions

      !
      ! Clean slave pointers
      !
      deallocate(data1M,data1O)
      nullify(data1M,data1O)
      nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
      nullify(p_StkM,p_StkO)
      if (associated(etaIM)) then
        deallocate(etaIM)
        nullify(etaIM)
      end if

      return

      end subroutine emergenceI_serial

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the solver for the NLTE problem for intensity with just
      !! the continuum opacities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence
      subroutine solveJ(Atmo,Atom,LTElines,Cont,Frec,Geom, &
                        MPID,Input,Stokes,J00C)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C

      ! MPI
      if (MPID%mpi) then
        if (MPID%alternJ) then
          call solverJ_alt(Atmo,Atom,LTElines,Cont,Frec,Geom, &
                           MPID,Input,Stokes,J00C)
        else
          call solverJ(Atmo,Atom,LTElines,Cont,Frec,Geom, &
                       MPID,Input,Stokes,J00C)
        end if
      ! Serial
      else
        call solverJ_serial(Atmo,Atom,LTElines,Cont,Frec,Geom, &
                            MPID,Input,Stokes,J00C)
      end if

      end subroutine solveJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem for intensity with just the continuum
      !! opacities, using lambda iteration, with several CPU.\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence
      subroutine solverJ(Atmo,Atom,LTElines,Cont,Frec, &
                         Geom,MPID,Input,Stokes,J00C)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C

      ! Local

      type(MRC_class):: MRC

      logical:: AD,ADD,PRDl,goout

      integer:: iaux,iz0,iz1,diz,m,o,p,op,ith,iph,iz,if0,if1
      integer:: iproc,iter,id,if0l,if1l,nfl
      integer:: ntpz,npz,ntp,itpz,itz
      integer:: istep

      double precision:: mu_inv,dsm,dsp,WA
      double precision, dimension(1):: daux1
      double precision, dimension(nfreq,Rz0:Rz1):: J00Cold

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Stokes_r
      ! Senders
      double precision, dimension(:,:,:,:), allocatable:: Stokes_s
      ! Dual
      integer:: info_b

      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP
      double precision, dimension(:,:), pointer:: p_MStk

      ! Dummy
      type(Red_class):: Red
      double precision, dimension(:), allocatable:: ad1
      double precision, dimension(:,:), allocatable:: ad2
      double precision, dimension(:,:,:), allocatable:: ad3


      ! Announce we are starting
      if (gpid.eq.0) then
        umsg = '   Iteration            MRC(J^0_0) Freq_index  '// &
               'Wavelength Height_index Height(km)'
        call verbose
      end if

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! B-B
      if (Input%init_J_bb) then
        ! Store if it was PRD
        PRDl = PRD
        PRD = .False.
      end if

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)


      !
      ! Allocation of buffers
      !

      ! Master
      if(pid.eq.0)then

        ! Dimensions
        npz = Geom%nph*Rnz
        ntp = Geom%nth*Geom%nph
        ntpz = Geom%nth*npz

        ! To receive Intensity chunks
        iaux = MPID%nxfreq*Geom%nTh*Geom%nPh*Rnz
        allocate(Stokes_r(iaux))

      ! Slave
      else

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(if0:if1,0:2))
        allocate(data1O(if0:if1,0:2))

        ! To send Intensity chunks
        allocate(Stokes_s(if0:if1,Rz0:Rz1,Geom%nPh,Geom%nTh))

      end if ! Master or Slave

      ! Control
      call control
      if (laborted) goto 2000


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=1,Input%iter_j

        !
        ! Master
        !
        if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(id,Stokes_r,iz,ith,iph,itz,itpz,WA) &
!$omp shared(J00Cold,J00C,MPID,if0l,if1l,nfl,npz,ntp,ntpz) &
!$omp shared(KSTK,p_MStk,Frec,nz,Geom,Stokes,info_b,laborted) &
!$omp shared(Rz0,Rz1,Rnz,ierr) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT)

          ! Save old radiation field
!$omp workshare
          J00Cold = J00C
!$omp end workshare

          ! Reset radiation field variables
!$omp workshare
          J00C = 0d0
!$omp end workshare

          ! Each frequency cut
          do id=1,MPID%nnd

!$omp single
            !
            ! Receive data from a slave
            !

            ! Receive indexing data
            do while (.True.)
              call MPI_recv(info_b,1, &
                            MPI_INTEGER, &
                            MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! Flag error
            if (info_b.lt.0) laborted = .True.
!$omp end single
            ! Continue?
            if (info_b.lt.0) cycle
!$omp single

            ! Receive intensity
            do while (.True.)
              call MPI_recv(Stokes_r(1), &
                            MPID%sizei4b(info_b), &
                            MPI_DOUBLE_PRECISION, info_b, &
                            info_b, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! Shorter variables
            if0l = MPID%if0(info_b)
            if1l = MPID%if1(info_b)
            nfl = MPID%nf(info_b)

            ! Pointer
            p_MStk(if0l:if1l,1:ntpz) => &
                                     Stokes_r(1:MPID%sizei4b(info_b))
!$omp end single
!$omp do
            ! For each height
            do iz=Rz0,Rz1

              ! For each polar direction
              do ith=1,Geom%nth

                ! Partial indexing
                itz = iz + (ith-1)*npz - Rz0 + 1

                ! For each azimuth
                do iph=1,Geom%nph

                  ! Get frequency weight
                  WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                  ! Get running index
                  itpz = itz + Rnz*(iph-1)

                  ! Determine where to store intensity
                  if (KSTK.or.iz.eq.Rz0) &
                    Stokes(if0l:if1l,iph,ith,iz) = &
                                                p_MStk(if0l:if1l,itpz)

                  !
                  ! Calculate frequency integral
                  !
                  call FIntJ(MPID,WA,info_b,p_MStk(:,itpz),J00C(:,iz))


                enddo ! azimuthal directions
              enddo ! polar directions
            enddo ! Heights
!$omp end do
          end do ! frequency domains
!$omp end parallel

          ! Nullify pointer
          nullify(p_MStk)

          !
          ! Calculate MRC for J
          !

          ! Call the routine
          call MRCJ_sb(J00C,J00Cold,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Write in stdout
          if (gpid.eq.0) then
            write(umsg,'(3x,"J it:",1x,i3,2x,es20.12,'// &
                       '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                       iter,MRC%values(2,1),MRC%indexes(1,1), &
                       1d2/Frec%omega(MRC%indexes(1,1)), &
                       MRC%indexes(2,1),MRC%values(1,1)
            call verbose
          end if

          ! Check convergence
          goout = MRC%values(2,1).lt.INPUT%mrcj


        !
        ! Slave
        !
        else


          !
          ! Ratiation Transfer
          !

          !  For each polar direction
          do ith=1,Geom%nTh

            ! Calculate inverse of cosine of polar direction
            mu_inv = 1d0/Geom%V_mu(ith)

            ! Determine the direction of propagation for indexes
            diz = -int(sign(1d0, Geom%V_mu(ith)))


            ! Determine the first and last height indexes to run
            ! over

            iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
            iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

            ! For each azimuthal direction
            do iph=1,Geom%nPh

              !
              ! First height
              !

              ! If going down, get top boundary
              if(diz.eq.1)then

                ! Call top boundary
                call topI(MPID,data1M(:,2))

              ! If going up, get bottom boundary
              else

                ! Call bottom boundary
                call bottomI(Frec%omega,Atmo%T(iz0), &
                             Atmo%vx(iz0),Atmo%vy(iz0), &
                             Atmo%vz(iz0),Geom%V_mu(ith), &
                             Geom%V_mux(iph),Geom%V_muy(iph), &
                             MPID,data1M(:,2))

              endif ! propagation direction

              ! Identify current height
              o = iz0

              ! Calculate radiative coefficients
              if (Input%init_J_bb) then
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              o,ith,iph,if0,if1,ad1, &
                              J00C(:,o),Cont%ndir, &
                              Cont%c(:,:,:,o),ad3, &
                              ad1,ad1,data1M(:,0:1),ad2,.False.)
              else
                call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,o), &
                              Cont%ndir,Cont%c(:,:,:,o), &
                              data1M(:,0:1))
              end if


              !
              ! Store in buffer
              !

              ! Intensity
              Stokes_s(:,o,iph,ith) = data1M(:,2)

              ! Identify next height
              p = iz0 + diz

              ! Calculate radiative coefficients
              if (Input%init_J_bb) then
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              p,ith,iph,if0,if1,ad1, &
                              J00C(:,p),Cont%ndir, &
                              Cont%c(:,:,:,p),ad3, &
                              ad1,ad1,data1O(:,0:1),ad2,.False.)
              else
                call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,p), &
                              Cont%ndir,Cont%c(:,:,:,p), &
                              data1O(:,0:1))
              end if

              !
              ! Intermedium heights
              !

              ! For each height this CPU has assigned
              do iz=iz0,iz1,diz

                ! We treat the boundaries outside
                if(iz.eq.iz0.or.iz.eq.iz1)cycle

                ! Allocate P pointers
                allocate(data1P(if0:if1,0:2))

                ! Identify heights
                m = iz - diz
                o = iz
                p = iz + diz

                ! Calculate distance to previous point
                dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                ! Calculate distance to the next point
                dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

                ! If tau scale
                if (ztau) then
                  dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(m))
                  dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(p))
                end if

                ! Calculate radiative coefficients
                if (Input%init_J_bb) then
                  call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID, &
                                Geom,p,ith,iph,if0,if1,ad1, &
                                J00C(:,p),Cont%ndir, &
                                Cont%c(:,:,:,p),ad3, &
                                ad1,ad1,data1P(:,0:1),ad2,.False.)
                else
                  call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,p), &
                                Cont%ndir,Cont%c(:,:,:,p), &
                                data1P(:,0:1))
                end if

                ! Point to the data
                p_K0M  => data1M(:,0)
                p_SM   => data1M(:,1)
                p_StkM => data1M(:,2)
                p_K0O  => data1O(:,0)
                p_SO   => data1O(:,1)
                p_StkO => data1O(:,2)
                p_K0P  => data1P(:,0)
                p_SP   => data1P(:,1)

                ! Apply short characteristics BESSER
                call RTStepI(o,ith,iph,MPID%nf(pid), &
                             dsm,dsp,p_K0M,p_SM,p_K0O, &
                             p_SO,p_K0P,p_SP,p_StkM, &
                             p_StkO,daux1,.False.,.True.)

                !
                ! Store in buffer
                !

                ! Send Intensity
                Stokes_s(:,o,iph,ith) = data1O(:,2)

                ! Shift data (O->M, P->O)
                deallocate(data1M)
                data1M => data1O
                data1O => data1P
                nullify(data1P)

              end do ! Intermedium heights


              !
              ! Last height
              !

              ! Identify heights
              m = iz1 - diz
              o = iz1

              ! Calculate distance to previous point
              dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

              ! If tau scale
              if (ztau) &
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))

              ! Point to the data
              p_K0M  => data1M(:,0)
              p_SM   => data1M(:,1)
              p_StkM => data1M(:,2)
              p_K0O  => data1O(:,0)
              p_SO   => data1O(:,1)
              p_StkO => data1O(:,2)

              ! Apply short characteristics LINEAR
              call RTStepI(o,ith,iph,MPID%nf(pid), &
                           dsm,dsp,p_K0M,p_SM,p_K0O, &
                           p_SO,p_K0P,p_SP,p_StkM, &
                           p_StkO,daux1,.False.,.False.)

              !
              ! Store in buffer
              !

              ! Intensity
              Stokes_s(:,o,iph,ith) = data1O(:,2)

            enddo ! azimuthal angles
          enddo ! polar angles


          !
          ! Send to master
          !

          ! If had an error
          if (laborted) then

            ! Send error
            do while (.True.)
              call MPI_SEND(-pid,1,MPI_INTEGER,0,0, &
                            MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do

          ! No problems
          else

            ! Send indexes
            do while (.True.)
              call MPI_SEND(pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                            ierr)
              if (ierr.eq.0) exit
            end do

            ! Send Stokes
            do while (.True.)
              call MPI_SEND(Stokes_s(if0,Rz0,1,1), &
                            MPID%sizei4b(pid), &
                            MPI_DOUBLE_PRECISION, 0, pid, &
                            MPI_COMM_RT, ierr)
              if (ierr.eq.0) exit
            end do

          end if ! Errors
        end if ! Master or Slave

        ! Control
        call control
        if (laborted) goto 2000

        !
        ! Share the radiation information
        !

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive J00C
            call MPI_RECV(J00C(1,Rz0), MPID%sizei7(0), &
                          MPI_DOUBLE_PRECISION,  &
                          MPID%recv, 2000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Receive intensity if doing A-D PRD
            if (PRD.and.ADD) &
              call MPI_RECV(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, 3000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            ! Finished?
            call MPI_RECV(goout,1,MPI_LOGICAL,MPID%recv, &
                          4000000+pid,MPI_COMM_RT, &
                          MPI_STATUS_IGNORE,ierr)

          end if ! No master

          ! For each send
          do istep=1,MPID%nsend

            ! Send J00C
            call MPI_ISEND(J00C(1,Rz0), MPID%sizei7(0), &
                           MPI_DOUBLE_PRECISION, &
                           MPID%lsend(istep), &
                           2000000+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,3), ierr)

            ! Send intensity if doing A-D PRD
            if (PRD.and.ADD) &
              call MPI_ISEND(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             3000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,4), ierr)

            ! Finished
            call MPI_ISEND(goout,1,MPI_LOGICAL,MPID%lsend(istep), &
                           4000000+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,5), ierr)

          end do ! Sends

        ! Normal bcast
        else

          ! Share J00C
          call MPI_BCAST(J00C(1,Rz0), MPID%sizei7(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Share intensity if doing A-D PRD
          if (PRD.and.ADD) &
            call MPI_BCAST(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Finished
          call MPI_BCAST(goout,1,MPI_LOGICAL,0, &
                         MPI_COMM_RT, ierr)

        end if ! Type of bcast

        ! If alternative bcast
        if (MPID%nsend.gt.0.and.MPID%altbcast) then

          ! For each slave
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing and reset the buffers
            call MPI_WAIT(MPID%requestA(iproc,3), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,4), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,5), &
                          MPI_STATUS_IGNORE,ierr)

          end do ! Processors

        end if ! Domain decomposition

        ! Finish?
        if (goout) exit

      end do ! Iterations

      ! B-B
      if (Input%init_J_bb) then
        ! Restore PRD
        PRD = PRDl
      end if

      !
      ! Clean slave pointers
      !
2000  if (pid.ne.0) then
        deallocate(data1M,data1O)
        nullify(data1M,data1O)
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO)
      end if

      return

      end subroutine solverJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem for intensity with just the continuum
      !! opacities, using lambda iteration, with several CPU.
      !! Alternative communication scheme.\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence
      subroutine solverJ_alt(Atmo,Atom,LTElines,Cont,Frec, &
                             Geom,MPID,Input,Stokes,J00C)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C

      ! Local

      type(MRC_class):: MRC

      logical:: AD,ADD,PRDl,goout

      integer:: iaux,iz0,iz1,diz,m,o,p,op,ith,iph,ith1,iph1,iz
      integer:: if0,if1,iproc,iter,id
      integer:: if0l,if1l,nfl
      integer:: istep

      double precision:: mu_inv,dsm,dsp,WA
      double precision, dimension(1):: daux1
      double precision, dimension(nfreq,Rz0:Rz1):: J00Cold

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Stokes_r
      ! Senders
      double precision, dimension(:,:), allocatable:: Stokes_s
      ! Dual
      integer:: info_b
      integer,dimension(3):: info_c


      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP
      double precision, dimension(:,:), pointer:: p_MStk

      ! Dummy
      type(Red_class):: Red
      double precision, dimension(:), allocatable:: ad1
      double precision, dimension(:,:), allocatable:: ad2
      double precision, dimension(:,:,:), allocatable:: ad3


      ! Announce we are starting
      if (gpid.eq.0) then
        umsg = '   Iteration            MRC(J^0_0) Freq_index  '// &
               'Wavelength Height_index Height(km)'
        call verbose
      end if

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! B-B
      if (Input%init_J_bb) then
        ! Store if it was PRD
        PRDl = PRD
        PRD = .False.
      end if

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)


      !
      ! Allocation of buffers
      !

      ! Master
      if(pid.eq.0)then

        ! To receive Intensity chunks
        iaux = MPID%nxfreq*Rnz
        allocate(Stokes_r(iaux))

      ! Slave
      else

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(if0:if1,0:2))
        allocate(data1O(if0:if1,0:2))

        ! To send Intensity chunks
        allocate(Stokes_s(if0:if1,Rz0:Rz1))

      end if ! Master or Slave

      ! Control
      call control
      if (laborted) goto 2000


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=1,Input%iter_j

        !
        ! Master
        !
        if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(ith1,iph1,id,info_c,Stokes_r,WA) &
!$omp shared(J00C,J00Cold,MPID,ith,iph,info_b,if0l,if1l,nfl) &
!$omp shared(Geom,p_MStk,KSTK,nz,Stokes,laborted,Rz0,Rz1,Rnz) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT,ierr)

          ! Save old radiation field
!$omp workshare
          J00Cold = J00C
!$omp end workshare

          ! Reset radiation field variables
!$omp workshare
          J00C = 0d0
!$omp end workshare

          ! For each dimension in the problem
          ! Each polar direction
          do ith1=1,Geom%nTh

            ! Each azimuthal direction
            do iph1=1,Geom%nPh

              ! Each frequency cut
              do id=1,MPID%nnd

!$omp single
                !
                ! Receive data from a slave
                !

                ! Receive indexing data
                do while (.True.)
                  call MPI_recv(info_c(1),3, &
                                MPI_INTEGER, &
                                MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                                MPI_STATUS_IGNORE, ierr)
                  if (ierr.eq.0) exit
                end do

                info_b = info_c(1)
                ith = info_c(2)
                iph = info_c(3)

                ! Flag error
                if (info_b.lt.0) laborted = .True.
!$omp end single
                ! Continue?
                if (info_b.lt.0) cycle
!$omp single
                ! Receive intensity
                do while (.True.)
                  call MPI_recv(Stokes_r(1), &
                                MPID%sizei4b(info_b), &
                                MPI_DOUBLE_PRECISION, info_b, &
                                info_b, MPI_COMM_RT, &
                                MPI_STATUS_IGNORE, ierr)
                  if (ierr.eq.0) exit
                end do

                ! Shorter variables
                if0l = MPID%if0(info_b)
                if1l = MPID%if1(info_b)
                nfl = MPID%nf(info_b)

                ! Pointers
                p_MStk(if0l:if1l,Rz0:Rz1) => &
                                     Stokes_r(1:MPID%sizei4b(info_b))
!$omp end single

                ! Get angular weight
                WA = Geom%W_mu(ith)*Geom%W_mux(iph)

!$omp do
                ! Each height
                do iz=Rz0,Rz1

                  ! Determine where to store intensity
                  if (KSTK.or.iz.eq.Rz0) &
                    Stokes(if0l:if1l,iph,ith,iz) = &
                                                  p_MStk(if0l:if1l,iz)

                  !
                  ! Calculate frequency integral
                  !
                  call FIntJ(MPID,WA,info_b,p_MStk(:,iz),J00C(:,iz))

                end do ! heights
!$omp end do
              end do ! frequency domains
            enddo ! azimuthal directions
          enddo ! polar directions
!$omp end parallel

          ! Nullify pointer
          nullify(p_MStk)

          !
          ! Calculate MRC for J
          !

          ! Call the routine
          call MRCJ_sb(J00C,J00Cold,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Write in stdout
          if (gpid.eq.0) then
            write(umsg,'(3x,"J it:",1x,i3,2x,es20.12,'// &
                       '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                       iter,MRC%values(2,1),MRC%indexes(1,1), &
                       1d2/Frec%omega(MRC%indexes(1,1)), &
                       MRC%indexes(2,1),MRC%values(1,1)
            call verbose
          end if

          ! Check convergence
          goout = MRC%values(2,1).lt.INPUT%mrcj


        !
        ! Slave
        !
        else


          !
          ! Ratiation Transfer
          !

          !  For each polar direction
          do ith=1,Geom%nTh

            ! Calculate inverse of cosine of polar direction
            mu_inv = 1d0/Geom%V_mu(ith)

            ! Determine the direction of propagation for indexes
            diz = -int(sign(1d0, Geom%V_mu(ith)))

            ! Determine the first and last height indexes to run
            ! over

            iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
            iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

            ! For each azimuthal direction
            do iph=1,Geom%nPh

              !
              ! First height
              !

              ! If going down, get top boundary
              if(diz.eq.1)then

                ! Call top boundary
                call topI(MPID,data1M(:,2))

              ! If going up, get bottom boundary
              else

                ! Call bottom boundary
                call bottomI(Frec%omega,Atmo%T(iz0), &
                             Atmo%vx(iz0),Atmo%vy(iz0), &
                             Atmo%vz(iz0),Geom%V_mu(ith), &
                             Geom%V_mux(iph),Geom%V_muy(iph), &
                             MPID,data1M(:,2))

              endif ! propagation direction

              ! Identify current height
              o = iz0

              ! Calculate radiative coefficients
              if (Input%init_J_bb) then
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              o,ith,iph,if0,if1,ad1, &
                              J00C(:,o),Cont%ndir, &
                              Cont%c(:,:,:,o),ad3, &
                              ad1,ad1,data1M(:,0:1),ad2,.False.)
              else
                call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,o), &
                              Cont%ndir,Cont%c(:,:,:,o), &
                              data1M(:,0:1))
              end if


              !
              ! Store in buffer
              !

              ! Intensity
              Stokes_s(:,o) = data1M(:,2)

              ! Identify next height
              p = iz0 + diz

              ! Calculate radiative coefficients
              if (Input%init_J_bb) then
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                              p,ith,iph,if0,if1,ad1, &
                              J00C(:,p),Cont%ndir, &
                              Cont%c(:,:,:,p),ad3, &
                              ad1,ad1,data1O(:,0:1),ad2,.False.)
              else
                call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,p), &
                              Cont%ndir,Cont%c(:,:,:,p), &
                              data1O(:,0:1))
              end if

              !
              ! Intermedium heights
              !

              ! For each height this CPU has assigned
              do iz=iz0,iz1,diz

                ! We treat the boundaries outside
                if(iz.eq.iz0.or.iz.eq.iz1)cycle

                ! Allocate P pointers
                allocate(data1P(if0:if1,0:2))

                ! Identify heights
                m = iz - diz
                o = iz
                p = iz + diz

                ! Calculate distance to previous point
                dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

                ! Calculate distance to the next point
                dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

                ! If tau scale
                if (ztau) then
                  dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(m))
                  dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                      Atmo%chi500(p))
                end if

                ! Calculate radiative coefficients
                if (Input%init_J_bb) then
                  call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID, &
                                Geom,p,ith,iph,if0,if1,ad1, &
                                J00C(:,p),Cont%ndir, &
                                Cont%c(:,:,:,p),ad3, &
                                ad1,ad1,data1P(:,0:1),ad2,.False.)
                else
                  call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,p), &
                                Cont%ndir,Cont%c(:,:,:,p), &
                                data1P(:,0:1))
                end if

                ! Point to the data
                p_K0M  => data1M(:,0)
                p_SM   => data1M(:,1)
                p_StkM => data1M(:,2)
                p_K0O  => data1O(:,0)
                p_SO   => data1O(:,1)
                p_StkO => data1O(:,2)
                p_K0P  => data1P(:,0)
                p_SP   => data1P(:,1)

                ! Apply short characteristics BESSER
                call RTStepI(o,ith,iph,MPID%nf(pid), &
                             dsm,dsp,p_K0M,p_SM,p_K0O, &
                             p_SO,p_K0P,p_SP,p_StkM, &
                             p_StkO,daux1,.False.,.True.)


                !
                ! Store in buffer
                !

                ! Send Intensity
                Stokes_s(:,o) = data1O(:,2)

                ! Shift data (O->M, P->O)
                deallocate(data1M)
                data1M => data1O
                data1O => data1P
                nullify(data1P)

              end do ! Intermedium heights


              !
              ! Last height
              !

              ! Identify heights
              m = iz1 - diz
              o = iz1

              ! Calculate distance to previous point
              dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

              ! If tau scale
              if (ztau) &
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))

              ! Point to the data
              p_K0M  => data1M(:,0)
              p_SM   => data1M(:,1)
              p_StkM => data1M(:,2)
              p_K0O  => data1O(:,0)
              p_SO   => data1O(:,1)
              p_StkO => data1O(:,2)

              ! Apply short characteristics LINEAR
              call RTStepI(o,ith,iph,MPID%nf(pid), &
                           dsm,dsp,p_K0M,p_SM,p_K0O, &
                           p_SO,p_K0P,p_SP,p_StkM, &
                           p_StkO,daux1,.False.,.False.)

              !
              ! Send to master if error
              !

              ! If had an error
              if (laborted) then

                ! Send error
                info_c = (/ -pid, ith, iph /)
                do while (.True.)
                  call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                                MPI_COMM_RT,ierr)
                  if (ierr.eq.0) exit
                end do

                cycle

              end if

              !
              ! Store in buffer
              !

              ! Intensity
              Stokes_s(:,o) = data1O(:,2)

              !
              ! Send to master
              !

              ! Send indexes
              info_c = (/ pid, ith, iph /)
              do while (.True.)
                call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                              MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do

              ! Send Stokes
              do while (.True.)
                call MPI_SEND(Stokes_s(if0,Rz0), &
                              MPID%sizei4b(pid), &
                              MPI_DOUBLE_PRECISION, 0, pid, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

            enddo ! azimuthal angles
          enddo ! polar angles

        end if ! Master or Slave

        !
        ! Share the radiation information
        !

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive J00C
            call MPI_RECV(J00C(1,Rz0), MPID%sizei7(0), &
                          MPI_DOUBLE_PRECISION,  &
                          MPID%recv, 2000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Receive intensity if doing A-D PRD
            if (PRD.and.ADD) &
              call MPI_RECV(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, 3000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

          end if ! No master

          ! For each send
          do istep=1,MPID%nsend

            ! Send J00C
            call MPI_ISEND(J00C(1,Rz0), MPID%sizei7(0), &
                           MPI_DOUBLE_PRECISION, &
                           MPID%lsend(istep), &
                           2000000+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,3), ierr)

            ! Send intensity if doing A-D PRD
            if (PRD.and.ADD) &
              call MPI_ISEND(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             3000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,4), ierr)

            ! Finished?
            call MPI_RECV(goout,1,MPI_LOGICAL,MPID%recv, &
                          4000000+pid,MPI_COMM_RT, &
                          MPI_STATUS_IGNORE,ierr)

            ! Finished
            call MPI_ISEND(goout,1,MPI_LOGICAL,MPID%lsend(istep), &
                           4000000+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,5), ierr)

          end do ! Sends

        ! Normal bcast
        else

          ! Share J00C
          call MPI_BCAST(J00C(1,Rz0), MPID%sizei7(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Share intensity if doing A-D PRD
          if (PRD.and.ADD) &
            call MPI_BCAST(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Finished
          call MPI_BCAST(goout,1,MPI_LOGICAL,0, &
                         MPI_COMM_RT, ierr)

        end if ! Type of bcast

        ! If alternative bcast
        if (MPID%nsend.gt.0.and.MPID%altbcast) then

          ! For each slave
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing and reset the buffers
            call MPI_WAIT(MPID%requestA(iproc,3), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,4), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,5), &
                          MPI_STATUS_IGNORE,ierr)

          end do ! Processors

        end if ! Domain decomposition

        ! Control
        call control
        if (laborted) exit

        ! Finish?
        if (goout) exit

      end do ! Iterations

      ! B-B
      if (Input%init_J_bb) then
        ! Restore PRD
        PRD = PRDl
      end if

      !
      ! Clean slave pointers
      !
2000  if (pid.ne.0) then
        deallocate(data1M,data1O)
        nullify(data1M,data1O)
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO)
      end if

      return

      end subroutine solverJ_alt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem for intensity with just the continuum
      !! opacities, using lambda iteration, with one CPU (serial).\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence
      subroutine solverJ_serial(Atmo,Atom,LTElines,Cont,Frec, &
                                Geom,MPID,Input,Stokes,J00C)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target:: Stokes
      double precision, dimension(nfreq,Rz0:Rz1):: J00C

      ! Local

      type(MRC_class):: MRC

      logical:: AD,ADD,PRDl,goout

      integer:: iz0,iz1,diz,m,o,p
      integer:: iter,ith,iph,iz,if0,if1

      double precision:: mu_inv,dsm,dsp
      double precision, dimension(1):: daux1
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh, &
                                  giz0:giz1), target:: Stokes_n
      double precision, dimension(nfreq,Rz0:Rz1):: J00C_n

      ! Pointers
      double precision, dimension(:,:), pointer:: data1M,data1O,data1P
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP

      ! Dummy
      type(Red_class):: Red
      double precision, dimension(:), allocatable:: ad1
      double precision, dimension(:,:), allocatable:: ad2
      double precision, dimension(:,:,:), allocatable:: ad3


      ! Announce we are starting
      if (gpid.eq.0) then
        umsg = '   Iteration            MRC(J^0_0) Freq_index  '// &
               'Wavelength Height_index Height(km)'
        call verbose
      end if

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! B-B
      if (Input%init_J_bb) then
        ! Store if it was PRD
        PRDl = PRD
        PRD = .False.
      end if

      ! CPU limits
      if0 = 1
      if1 = nfreq

      !
      ! Allocations
      !

      ! Allocate M and O pointers for RT coeff
      allocate(data1M(nfreq,0:2))
      allocate(data1O(nfreq,0:2))


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=1,Input%iter_j

        ! Reset radiation field variables
        J00C_n = 0d0


        !
        ! Ratiation Transfer
        !

        !  For each polar direction
        do ith=1,Geom%nTh

          ! Calculate inverse of cosine of polar direction
          mu_inv = 1d0/Geom%V_mu(ith)

          ! Determine the direction of propagation for indexes
          diz = -int(sign(1d0, Geom%V_mu(ith)))

          ! Determine the first and last height indexes to run over
          iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
          iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

          ! For each azimuthal direction
          do iph=1,Geom%nPh


            !
            ! First height
            !

            ! If going down, get top boundary
            if(diz.eq.1)then

              ! Call top boundary
              call topI(MPID,data1M(:,2))

            ! If going up, get bottom boundary
            else

              ! Call bottom boundary
              call bottomI(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                           Atmo%vy(iz0),Atmo%vz(iz0), &
                           Geom%V_mu(ith),Geom%V_mux(iph), &
                           Geom%V_muy(iph),MPID,data1M(:,2))

            endif

            ! Identify current height
            o = iz0

            ! Calculate radiative coefficients
            if (Input%init_J_bb) then
              call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                            o,ith,iph,if0,if1,ad1, &
                            J00C(:,o),Cont%ndir, &
                            Cont%c(:,:,:,o),ad3, &
                            ad1,ad1,data1M(:,0:1),ad2,.False.)
            else
              call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,o), &
                            Cont%ndir,Cont%c(:,:,:,o),data1M(:,0:1))
            end if

            if (KSTK) Stokes_n(:,iph,ith,o) = data1M(:,2)

            !
            ! Calculate integrals
            !
            call JcalcJ(Geom,iph,ith,data1M(:,2),J00C_n(:,o))

            ! Identify next height
            p = iz0 + diz

            ! Calculate radiative coefficients
            if (Input%init_J_bb) then
              call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom, &
                            p,ith,iph,if0,if1,ad1, &
                            J00C(:,p),Cont%ndir, &
                            Cont%c(:,:,:,p),ad3, &
                            ad1,ad1,data1O(:,0:1),ad2,.False.)
            else
              call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,p), &
                            Cont%ndir,Cont%c(:,:,:,p),data1O(:,0:1))
            end if


            !
            ! Intermedium heights
            !

            ! For each height this CPU has assigned
            do iz=iz0,iz1,diz

              ! We treat the boundaries outside
              if(iz.eq.iz0.or.iz.eq.iz1)cycle

              ! Allocate P pointers
              allocate(data1P(nfreq,0:2))

              ! Identify heights
              m = iz - diz
              o = iz
              p = iz + diz

              ! Calculate distance to previous point
              dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

              ! Caculate quantities of the next point
              dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

              ! If tau scale
              if (ztau) then
                dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(m))
                dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                    Atmo%chi500(p))
              end if

              ! RT coefficients
              if (Input%init_J_bb) then
                call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID, &
                              Geom,p,ith,iph,if0,if1,ad1, &
                              J00C(:,p),Cont%ndir, &
                              Cont%c(:,:,:,p),ad3, &
                              ad1,ad1,data1P(:,0:1),ad2,.False.)
              else
                call RTCoeffJ(Geom,ith,iph,if0,if1,J00C(:,p), &
                              Cont%ndir,Cont%c(:,:,:,p),data1P(:,0:1))
              end if

              ! Point to the data
              p_K0M  => data1M(:,0)
              p_SM   => data1M(:,1)
              p_StkM => data1M(:,2)
              p_K0O  => data1O(:,0)
              p_SO   => data1O(:,1)
              p_StkO => data1O(:,2)
              p_K0P  => data1P(:,0)
              p_SP   => data1P(:,1)

              ! Apply short characteristics BESSER
              call RTStepI(o,ith,iph,nfreq, &
                           dsm,dsp,p_K0M,p_SM,p_K0O, &
                           p_SO,p_K0P,p_SP,p_StkM, &
                           p_StkO,daux1,.False.,.True.)

              if (KSTK) Stokes_n(:,iph,ith,o) = data1O(:,2)


              !
              ! Calculate integrals
              !
              call JcalcJ(Geom,iph,ith,data1O(:,2),J00C_n(:,o))

              ! Shift data (O->M, P->O)
              deallocate(data1M)
              data1M => data1O
              data1O => data1P
              nullify(data1P)

            end do

            !
            ! Last height
            !

            ! Identify heights
            m = iz1 - diz
            o = iz1

            ! Calculate distance to previous point
            dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

            ! If tau scale
            if (ztau) &
              dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(m))

            ! Point to the data
            p_K0M  => data1M(:,0)
            p_SM   => data1M(:,1)
            p_StkM => data1M(:,2)
            p_K0O  => data1O(:,0)
            p_SO   => data1O(:,1)
            p_StkO => data1O(:,2)

            ! Apply short characteristics LINEAR
            call RTStepI(o,ith,iph,MPID%nf(pid), &
                         dsm,dsp,p_K0M,p_SM,p_K0O, &
                         p_SO,p_K0P,p_SP,p_StkM, &
                         p_StkO,daux1,.False.,.False.)

            if (KSTK.or.o.eq.1) &
              Stokes_n(:,iph,ith,o) = data1O(:,2)


            !
            ! Calculate integrals
            !
            call JcalcJ(Geom,iph,ith,data1O(:,2),J00C_n(:,o))

          enddo ! azimuthal directions
        enddo ! polar directions

        !
        ! Calculate MRC for J
        !

        ! Call the routine
        call MRCJ_sb(J00C_n,J00C,MRC)

        ! Convert cm into km
        MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

        ! Write in stdout
        if (gpid.eq.0) then
          write(umsg,'(3x,"J it:",1x,i3,2x,es20.12,'// &
                     '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                     iter,MRC%values(2,1),MRC%indexes(1,1), &
                     1d2/Frec%omega(MRC%indexes(1,1)), &
                     MRC%indexes(2,1),MRC%values(1,1)
          call verbose
        end if

        ! Check convergence
        goout = MRC%values(2,1).lt.INPUT%mrcj

        ! Shift the new values into the proper variables
        Stokes = Stokes_n
        J00C = J00C_n

        ! Control
        if (laborted) exit

        ! Finish?
        if (goout) exit

      end do ! Iterations

      ! B-B
      if (Input%init_J_bb) then
        ! Restore PRD
        PRD = PRDl
      end if

      !
      ! Clean pointers
      !
      deallocate(data1M,data1O)
      nullify(data1M,data1O)
      nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
      nullify(p_StkM,p_StkO)

      return

      end subroutine solverJ_serial

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the generator of JKQ for multi-term\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            Pcorr(logical): Bool that says if the mean
      !!                            intensities are to be corrected\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!             rnPh(integer): Allocation size for Stokes\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine JKQgenerate(Atom,Rho_old,Atmo,Frec,Geom,MPID,Input, &
                             Flgsg,Pcorr,Bfield,rnPh, &
                             Stokes0,J00,J00S,J00C, &
                             Stokes,JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Frequency_class):: Frec
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: Pcorr
      integer, intent(in):: rnPh
      double precision, dimension(:,:,:,:), allocatable:: Stokes0
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! MPI
      if (MPID%mpi) then
        if (MPID%alternJgen) then
          call JKQgen_alt(Atom,Rho_old,Atmo,Frec,Geom, &
                          MPID,Input,Flgsg,Pcorr,Bfield,rnPh, &
                          Stokes0,J00,J00S,J00C, &
                          Stokes,JKQ,JKQS,JKQC,J00P)
        else
          call JKQgen(Atom,Rho_old,Atmo,Frec,Geom,MPID, &
                      Input,Flgsg,Pcorr,Bfield,rnPh, &
                      Stokes0,J00,J00S,J00C, &
                      Stokes,JKQ,JKQS,JKQC,J00P)
        end if
      ! Serial
      else
        call JKQgen_serial(Atom,Rho_old,Atmo,Frec,Geom, &
                           MPID,Input,Flgsg,Pcorr,Bfield,rnPh, &
                           Stokes0,J00,J00S,J00C, &
                           Stokes,JKQ,JKQS,JKQC,J00P)
      end if

      end subroutine JKQgenerate

!#####################################################################
!#####################################################################
!#####################################################################

      !> Does one NLTE iteration initializing with multilevel and only
      !! intensity data, in order to generate a consistent set of
      !! populations and mean intensities for the correspondent
      !! multiterm model atom, with several CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            Pcorr(logical): Bool that says if the mean
      !!                            intensities are to be corrected\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!             rnPh(integer): Allocation size for Stokes\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine JKQgen(Atom,Rho_old,Atmo,Frec,Geom,MPID,Input, &
                        Flgsg,Pcorr,Bfield,rnPh, &
                        Stokes0,J00,J00S,J00C, &
                        Stokes,JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Frequency_class):: Frec
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: Pcorr
      integer, intent(in):: rnPh
      double precision, dimension(:,:,:,:), allocatable:: Stokes0
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      type(MRC_class):: MRC

      logical:: AD,ADD

      integer:: iaux,if0,if1,if0l,if1l,nfl,nftl,istep
      integer:: ith,iph,ia,iz,ifreq,itran,jtran,jftran,fftran,ftran
      integer:: nth,nph,id,iproc,tid
      integer:: ntpz,npz,ntp,itpz
      integer, dimension(0:nproc-1):: psize

      double precision:: WA,daux
      double precision, dimension(:,:), allocatable:: LambdaL
      double precision, dimension(:,:,:), allocatable:: LambdaP
      double precision, dimension(:,:,:,:,:), allocatable:: Norm
      double precision, dimension(:,:,:,:,:), allocatable:: BStk

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Prof_r
      ! Senders
      double precision, dimension(:,:,:,:,:), allocatable:: Prof_s
      ! Duals
      integer:: info_b

      ! Pointers
      double precision, dimension(:,:), pointer:: data2O
      double precision, dimension(:,:,:), pointer:: p_MProf

      ! Trick to have AV input for AD calculation
      if (tbAD) AV = .True.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! Select the angular limits
      if (PRD.and.ADD) then
        nth = Geom%nTh
        nph = Geom%nPh
      else
        nth = 1
        nph = 1
      end if

      ! Calculate size of transmition package
      do iproc=0,nproc-1
        psize(iproc) = 2*Frec%Mntfreq(iproc)*nph*nth*Rnz
      end do

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)

      ! Initialize lambda operator to 0
      iaux = 1
      istep = 1
      do ia=1,nA
        if (Atom(ia)%nftran.gt.iaux) iaux = Atom(ia)%nftran
        if (Atom(ia)%nphot.gt.istep) istep = Atom(ia)%nphot
      end do
      allocate(LambdaL(nxb,iaux))
      LambdaL = 0d0
      allocate(LambdaP(nxb,istep,2))
      LambdaP = 0d0


      !
      ! Allocation of buffers
      !

      ! Master
      if(pid.eq.0)then

        ! Dimensions
        npz = nPh*Rnz
        ntp = nTh*nPh
        ntpz = nTh*npz

        ! Norm
        allocate(Norm(2,nxtran,nph,nth,Rz0:Rz1))

        ! Cumulative JKQ
        allocate(BStk(2,nxtran,nph,nth,Rz0:Rz1))

        ! To receive profile information
        iaux = MPID%nxtfreq*2*nTh*nPh*Rnz
        allocate(Prof_r(iaux))

      ! Slave
      else

        ! Allocate O pointers
        allocate(data2O(Frec%ntfreq,2))

        ! To send profile information
        allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1,nph,nth))

      end if ! Master or Slave


      ! For each atom
      do ia=1,nA
        Rho_old(ia)%crho = Atom(ia)%crho
      end do

      !
      ! Allocate Radiation field
      !
      ! Allocate JKQ and JKQS and initialize to 0
      allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
      allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))


      ! Control
      call control
      if (laborted) goto 2000


      !
      ! Master
      !
      if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(id,itpz,ith,iph,iz,ia,WA,itran,jtran,daux,Prof_r) &
!$omp shared(JKQ,JKQs,Norm,Bstk,info_b,if0l,if1l,nfl,nftl) &
!$omp shared(p_MProf,ntpz,npz,nTh,nPh,nz,Atom,MPID,Frec) &
!$omp shared(Stokes0,J00C,NA,stm,psize,PRD,ADD,Geom,laborted) &
!$omp shared(Rz0,Rz1,Rnz,axial,axiali,ierr) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT)

        ! Reset radiation field variables
!$omp workshare
        JKQ = cZero
        JKQS = cZero
        Norm = 0d0
        BStk = 0d0
!$omp end workshare

        ! Each frequency cut
        do id=1,MPID%nnd

!$omp single
          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_b,1, MPI_INTEGER, MPI_ANY_SOURCE, &
                          0, MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)
            if (ierr.eq.0) exit
          end do

          ! Flag error
          if (info_b.lt.0) laborted = .True.
!$omp end single

          ! Continue?
          if (info_b.lt.0) cycle

!$omp single

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), psize(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do


          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreq(info_b)

          ! Pointers
          p_MProf(1:nftl,1:2,1:ntpz) => Prof_r(1:psize(info_b))

!$omp end single

          ! Compute line quantities
!$omp do
          do itpz=1,ntpz

            ! Get indexes
            ith = (itpz-1)/npz
            iph = (itpz - npz*ith - 1)/Rnz
            iz = itpz - Rnz*iph - npz*ith + Rz0 - 1
            ith = ith + 1
            iph = iph + 1

            !
            ! Calculate frequency integral
            !
            if (PRD.and.ADD) then
              if (axiali) then
                call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                            Stokes0(:,1,ith,iz), &
                            p_MProf(:,:,itpz),Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz))
              else
                call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                            Stokes0(:,iph,ith,iz), &
                            p_MProf(:,:,itpz),Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz))
              end if
            else
              call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                          J00C(:,iz),p_MProf(:,:,itpz), &
                          Norm(:,:,iph,ith,iz), &
                          BStk(:,:,iph,ith,iz))

            end if

          end do ! heights/directions
!$omp end do
        end do ! frequency domains

!$omp single
        ! Nullify pointers
        nullify(p_MProf)
!$omp end single


        !
        ! Apply weights to JKQ, JKQS, and normalize if
        !

        ! For each height
!$omp do
        do iz=Rz0,Rz1

          ! For each polar direction
          do ith=1,nTh

            ! For each azimuthal direction
            do iph=1,nph

              ! Get the angular integral weight
              if (PRD.and.ADD) then
                WA = Geom%W_mu(ith)*Geom%W_mux(iph)
              else
                WA = 1d0
              end if

              ! For each atom
              do ia=1,nA

                ! For each FS transition
                do itran=1,Atom(ia)%ntran

                  ! Apply shift
                  jtran = itran + Atom(ia)%tshift

                  ! Get the weight
                  if (Norm(1,jtran,iph,ith,iz).gt.0d0) then

                    daux = WA/Norm(1,jtran,iph,ith,iz)

                    ! Integrate angle
                    JKQ(0,0,jtran,iz) = JKQ(0,0,jtran,iz) + &
                          dcmplx(BStk(1,jtran,iph,ith,iz)*daux, 0d0)

                  end if

                  ! If there is stimulated emission
                  if (stm) then

                    ! Get the weight
                    if (Norm(2,jtran,iph,ith,iz).gt.0d0) then

                      daux = WA/Norm(2,jtran,iph,ith,iz)

                      ! Integrate angle
                      JKQS(0,0,jtran,iz) = JKQS(0,0,jtran,iz) + &
                          dcmplx(BStk(2,jtran,iph,ith,iz)*daux, 0d0)

                    end if

                  end if ! Stimulated emission
                end do ! transitions
              end do ! atoms
            end do ! azimuthal directions
          end do ! polar directions
        end do ! heights
!$omp end do
!$omp end parallel

      !
      ! Slave
      !
      else


        !  For each polar direction
        do ith=1,nTh

          ! For each azimuthal direction
          do iph=1,nPh

            ! For each height this CPU has assigned
            do iz=Rz0,Rz1

              ! RT coefficients
              call Termprof(Frec,Atom,Atmo,MPID,Flgsg,Geom, &
                            Bfield,iz,ith,iph,if0,if1,data2O)


              !
              ! Store in buffer
              !

              ! Send profiles
              Prof_s(:,:,iz,iph,ith) = data2O(:,:)

            end do ! Intermedium heights
          enddo ! azimuthal angles
        enddo ! polar angles


        !
        ! Send to master
        !

        ! If had an error
        if (laborted) then

          ! Send indexes
          do while (.True.)
            call MPI_SEND(-pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

        ! No problems
        else

          ! Send indexes
          do while (.True.)
            call MPI_SEND(pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

          ! Send profiles
          do while (.True.)
            call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                          psize(pid), MPI_DOUBLE_PRECISION, &
                          0, pid, MPI_COMM_RT, ierr)
            if (ierr.eq.0) exit
          end do
        end if ! Errors
      end if ! Master or Slave


      !
      ! Share the radiation information
      !

      ! Control
      call control
      if (laborted) goto 2000

      ! If alternative bcast
      if (MPID%altbcast) then

        ! If not master, receive first
        if (pid.ne.0) then

          ! Receive JKQ
          call MPI_RECV(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                        MPI_DOUBLE_COMPLEX,  &
                        MPID%recv, pid, &
                        MPI_COMM_RT, MPI_STATUS_IGNORE, &
                        ierr)

          ! Receive JKQS if stimulated emission
          if (stm) &
            call MPI_RECV(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                          MPI_DOUBLE_COMPLEX, &
                          MPID%recv, 1000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

        end if ! No Master

        ! For each send
        do istep=1,MPID%nsend

          ! Send JKQ
          call MPI_ISEND(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                         MPI_DOUBLE_COMPLEX, &
                         MPID%lsend(istep), &
                         MPID%lsend(istep), &
                         MPI_COMM_RT, &
                         MPID%requestA(istep,1), ierr)

          ! Send JKQS if stimulated emission
          if (stm) &
            call MPI_ISEND(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                           MPI_DOUBLE_COMPLEX, &
                           MPID%lsend(istep), &
                           1000000+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,2), ierr)

        end do ! Sends

      ! Normal bcast
      else

        ! Share JKQ
        call MPI_BCAST(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                       MPI_DOUBLE_COMPLEX, 0, &
                       MPI_COMM_RT, ierr)

        ! Share JKQS if stimulated emission
        if (stm) &
          call MPI_BCAST(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

      end if ! Type of bcast

      ! If Polarization initial correction
      if (Pcorr) then

        !
        ! Make the J00 and J00S flat
        !

        ! For each height
        do iz=Rz0,Rz1

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Apply shift
              jtran = itran + Atom(ia)%tshift

              ! For each FS transition within this transition
              do ftran=1,Atom(ia)%fst(itran)%nt

                ! Get the global transition index
                fftran = Atom(ia)%ifst_ij(ftran,itran)

                ! Apply shift
                jftran = fftran + Atom(ia)%tfshift

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00(jftran,iz) = dble(JKQ(0,0,jtran,iz))

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00S(jftran,iz) = dble(JKQS(0,0,jtran,iz))

              end do ! fs transition
            end do ! term transition
          end do ! atom
        end do ! height


        !
        ! Solve SEE
        !

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,fftran,jftran,umsg,urou,tid) &
!$omp shared(nA,nz,Atom,Rho_old,J00,J00S,J00P,LambdaL,LambdaP) &
!$omp shared(Rz0,Rz1,laborted,vaborted)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,tid,INPUT)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,tid)
#endif

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel

        !
        !  Calculate MRC
        !

        ! Only the master does
        if (pid.eq.0) then

          ! Calculate MRC
          call MRC_sb(Atom,Rho_old,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Write in stdout
          if (gpid.eq.0) then
            write(umsg,'(A,1x,es22.12)') ' - The conversion to '// &
                       'terms have changed the populations a '// &
                       'maximum of ',MRC%values(2,1)
            call verbose
          end if

        end if

        ! If alternative bcast
        if (MPID%nsend.gt.0.and.MPID%altbcast) then

          ! For each slave
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing and reset the buffers
            call MPI_WAIT(MPID%requestA(iproc,1), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,2), &
                          MPI_STATUS_IGNORE,ierr)

          end do ! Processors

        end if ! Alternative bcast
      end if ! P correction


      !
      ! Convert to polarization quantities
      !

      !
      ! Stokes
      !

!$omp parallel default(none) &
!$omp private(iz,ifreq) &
!$omp shared(KSTK,Stokes,Stokes0,J00,J00S,JKQC,J00C,nz,nfreq,Geom) &
!$omp shared(Rz0,Rz1,axial,axiali)

      ! Allocate the extra Stokes parameters
!$omp single
      if (KSTK) then
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz1))
      else
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz0+1))
      end if
!$omp end single

      ! Initialize to 0
!$omp workshare
      Stokes = 0d0
!$omp end workshare

      ! Get the intensity from the intensity array
!$omp single
      if (KSTK) then
        if (axiali.and..not.axial) then
          do iph=1,Geom%nPh
            Stokes(0,:,iph,:,:) = Stokes0(:,1,:,:)
          end do
        else
          Stokes(0,:,:,:,:) = Stokes0
        end if
      end if
!$omp end single

!$omp single
      ! Deallocate the intensity array
      deallocate(Stokes0)

      ! Deallocate J00 and J00S
      deallocate(J00)
      deallocate(J00S)


      !
      ! JKQ frequency dependent
      !

      ! Allocate JKQC
      allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
!$omp end single
!$omp workshare
      JKQC = cZero
!$omp end workshare

      ! Get J00C into the complex array
!$omp do
      do iz=Rz0,Rz1
        do ifreq=1,nfreq

          JKQC(0,0,ifreq,iz) = dcmplx(J00C(ifreq,iz),0d0)

        end do ! frequencies
      end do ! heights
!$omp end do
!$omp end parallel

      ! Deallocate J00C
      deallocate(J00C)

      ! Control
      call control

      ! Put back the AD
2000  if (tbAD) AV = .False.

      !
      ! Clean slave pointers
      !
      if (pid.ne.0) then
        deallocate(data2O)
        nullify(data2O)
      end if

      return

      end subroutine JKQgen

!#####################################################################
!#####################################################################
!#####################################################################

      !> Does one NLTE iteration initializing with multilevel and only
      !! intensity data, in order to generate a consistent set of
      !! populations and mean intensities for the correspondent
      !! multiterm model atom, with several CPU. Alternative
      !! communication.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            Pcorr(logical): Bool that says if the mean
      !!                            intensities are to be corrected\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!             rnPh(integer): Allocation size for Stokes\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine JKQgen_alt(Atom,Rho_old,Atmo,Frec,Geom, &
                            MPID,Input,Flgsg,Pcorr,Bfield,rnPh, &
                            Stokes0,J00,J00S,J00C, &
                            Stokes,JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Frequency_class):: Frec
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: Pcorr
      integer, intent(in):: rnPh
      double precision, dimension(:,:,:,:), allocatable:: Stokes0
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      type(MRC_class):: MRC

      logical:: AD,ADD

      integer:: iaux,if0,if1,if0l,if1l,nfl,nftl,istep
      integer:: iproc,ith,iph,ith1,iph1,ia,iz,ifreq,itran,jtran
      integer:: id,fftran,ftran,jftran,nth,nph,tid
      integer, dimension(0:nproc-1):: psize

      double precision:: WA,daux
      double precision, dimension(:,:), allocatable:: LambdaL
      double precision, dimension(:,:,:), allocatable:: LambdaP
      double precision, dimension(:,:,:,:,:), allocatable:: Norm
      double precision, dimension(:,:,:,:,:), allocatable:: BStk

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Prof_r
      ! Senders
      double precision, dimension(:,:,:), allocatable:: Prof_s
      ! Duals
      integer:: info_b
      integer, dimension(3):: info_c

      ! Pointers
      double precision, dimension(:,:), pointer:: data2O
      double precision, dimension(:,:,:), pointer:: p_MProf

      ! Trick to have AV input for AD calculation
      if (tbAD) AV = .True.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! Select the angular limits
      if (PRD.and.ADD) then
        nth = Geom%nTh
        nph = Geom%nPh
      else
        nth = 1
        nph = 1
      end if

      ! Calculate size of transmition package
      do iproc=0,nproc-1
        psize(iproc) = 2*Frec%ntfreq*Rnz
      end do

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)


      ! Initialize lambda operator to 0
      iaux = 0
      istep = 0
      do ia=1,nA
        if (Atom(ia)%nftran.gt.iaux) iaux = Atom(ia)%nftran
        if (Atom(ia)%nphot.gt.istep) istep = Atom(ia)%nphot
      end do
      allocate(LambdaL(nxb,iaux))
      LambdaL = 0d0
      allocate(LambdaP(nxb,istep,2))
      LambdaP = 0d0


      !
      ! Allocation of buffers
      !

      ! Master
      if(pid.eq.0)then

        ! Norm
        allocate(Norm(2,nxtran,nph,nth,Rz0:Rz1))

        ! Cumulative JKQ
        allocate(BStk(2,nxtran,nph,nth,Rz0:Rz1))

        ! To receive profile information
        iaux = MPID%nxtfreq*2*Rnz
        allocate(Prof_r(iaux))

      ! Slave
      else

        ! Allocate O pointers
        allocate(data2O(Frec%ntfreq,2))

        ! To send profile information
        allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1))

      end if ! Master or Slave


      ! For each atom
      do ia=1,nA
        Rho_old(ia)%crho = Atom(ia)%crho
      end do

      !
      ! Allocate Radiation field
      !
      ! Allocate JKQ and JKQS and initialize to 0
      allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
      allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))

      ! Control
      call control
      if (laborted) goto 2000


      !
      ! Master
      !
      if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(ith1,iph1,id,info_c,iz,WA,itran,ia,jtran,daux,Prof_r) &
!$omp shared(JKQ,JKQS,Norm,BStk,info_b,ith,iph,psize,if0l,if1l,MPID) &
!$omp shared(nfl,nftl,p_MProf,nth,nph,Stokes0,J00C,stm,na,Frec,nz) &
!$omp shared(Rz0,Rz1,Rnz,PRD,ADD,Atom,Geom,laborted,axial,axiali) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT,ierr)

        ! Reset radiation field variables
!$omp workshare
        JKQ = cZero
        JKQS = cZero
        Norm = 0d0
        BStk = 0d0
!$omp end workshare

        ! For each dimension in the problem
        ! Each polar direction
        do ith1=1,nTh

          ! Each azimuthal direction
          do iph1=1,nPh

            ! Each frequency cut
            do id=1,MPID%nnd

!$omp single
              !
              ! Receive data from a slave
              !

              ! Receive indexing data
              do while (.True.)
                call MPI_recv(info_c(1), 3, MPI_INTEGER, &
                              MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do

              info_b = info_c(1)
              ith = info_c(2)
              iph = info_c(3)

              ! Flag error
              if (info_b.lt.0) laborted = .True.
!$omp end single
              ! Continue?
              if (info_b.lt.0) cycle
!$omp single
              ! Receive profile
              do while (.True.)
                call MPI_recv(Prof_r(1), psize(info_b), &
                              MPI_DOUBLE_PRECISION, info_b, &
                              info_b, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do

              ! Shorter variables
              if0l = MPID%if0(info_b)
              if1l = MPID%if1(info_b)
              nfl = MPID%nf(info_b)
              nftl = Frec%Mntfreq(info_b)

              ! Pointers
              p_MProf(1:nftl,1:2,Rz0:Rz1) => Prof_r(1:psize(info_b))
!$omp end single

              ! Each height
!$omp do
              do iz=Rz0,Rz1

                !
                ! Calculate frequency integral
                !
                if (PRD.and.ADD) then
                  if (axiali) then
                    call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                                Stokes0(:,1,ith,iz), &
                                p_MProf(:,:,iz),Norm(:,:,iph,ith,iz), &
                                BStk(:,:,iph,ith,iz))
                  else
                    call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                                Stokes0(:,iph,ith,iz), &
                                p_MProf(:,:,iz),Norm(:,:,iph,ith,iz), &
                                BStk(:,:,iph,ith,iz))
                  end if
                else
                  call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                              J00C(:,iz),p_MProf(:,:,iz), &
                              Norm(:,:,iph,ith,iz), &
                              BStk(:,:,iph,ith,iz))

                end if

              end do ! heights
!$omp end do

            end do ! frequency domains
          enddo ! azimuthal directions
        enddo ! polar directions
!$omp single
        ! Nullify pointers
        nullify(p_MProf)
!$omp end single


        !
        ! Apply weights to JKQ, JKQS, and normalize if
        !

        ! For each height
!$omp do
        do iz=Rz0,Rz1

          ! For each polar direction
          do ith1=1,nTh

            ! For each azimuthal direction
            do iph1=1,nph

              ! Get the angular integral weight
              if (PRD.and.ADD) then
                WA = Geom%W_mu(ith1)*Geom%W_mux(iph1)
              else
                WA = 1d0
              end if

              ! For each atom
              do ia=1,nA

                ! For each FS transition
                do itran=1,Atom(ia)%ntran

                  ! Apply atomic shift
                  jtran = itran + Atom(ia)%tshift

                  ! Get the weight
                  if (Norm(1,jtran,iph1,ith1,iz).gt.0d0) then

                    daux = WA/Norm(1,jtran,iph1,ith1,iz)

                    ! Integrate angle
                    JKQ(0,0,jtran,iz) = JKQ(0,0,jtran,iz) + &
                          dcmplx(BStk(1,jtran,iph1,ith1,iz)*daux, 0d0)

                  end if

                  ! If there is stimulated emission
                  if (stm) then

                    ! Get the weight
                    if (Norm(2,jtran,iph1,ith1,iz).gt.0d0) then

                      daux = WA/Norm(2,jtran,iph1,ith1,iz)

                      ! Integrate angle
                      JKQS(0,0,jtran,iz) = JKQS(0,0,jtran,iz) + &
                          dcmplx(BStk(2,jtran,iph1,ith1,iz)*daux, 0d0)

                    end if

                  end if ! Stimulated emission
                end do ! transitions
              end do ! atoms
            end do ! azimuthal directions
          end do ! polar directions
        end do ! heights
!$omp end do
!$omp end parallel

      !
      ! Slave
      !
      else

        !  For each polar direction
        do ith=1,nTh

          ! For each azimuthal direction
          do iph=1,nPh

            ! For each height this CPU has assigned
            do iz=Rz0,Rz1

              ! RT coefficients
              call Termprof(Frec,Atom,Atmo,MPID,Flgsg,Geom, &
                            Bfield,iz,ith,iph,if0,if1,data2O)


              !
              ! Store in buffer
              !

              ! Send profiles
              Prof_s(:,:,iz) = data2O(:,:)

            end do ! Intermedium heights

            !
            ! Send to master if error
            !

            ! If had an error
            if (laborted) then

              ! Send error
              info_c = (/ -pid, ith, iph /)
              do while (.True.)
                call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                              MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do

              cycle

            end if

            !
            ! Send to master
            !

            ! Send indexes
            info_c = (/ pid, ith, iph /)
            do while (.True.)
              call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                            MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do

            ! Send profiles
            do while (.True.)
              call MPI_SEND(Prof_s(1,1,Rz0), &
                            psize(pid), MPI_DOUBLE_PRECISION, &
                            0, pid, MPI_COMM_RT, ierr)
              if (ierr.eq.0) exit
            end do

          enddo ! azimuthal angles
        enddo ! polar angles

      end if ! Master or Slave


      !
      ! Share the radiation information
      !

      ! Control
      call control
      if (laborted) goto 2000

      ! Alternative bcast
      if (MPID%altbcast) then

        ! If not master, receive first
        if (pid.ne.0) then

          ! Receive JKQ
          call MPI_RECV(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                        MPI_DOUBLE_COMPLEX,  &
                        MPID%recv, pid, &
                        MPI_COMM_RT, MPI_STATUS_IGNORE, &
                        ierr)

          ! Receive JKQS if stimulated emission
          if (stm) &
            call MPI_RECV(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                          MPI_DOUBLE_COMPLEX, &
                          MPID%recv, 1000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

        end if ! No Master

        ! For each send
        do istep=1,MPID%nsend

          ! Send JKQ
          call MPI_ISEND(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                         MPI_DOUBLE_COMPLEX, &
                         MPID%lsend(istep), &
                         MPID%lsend(istep), &
                         MPI_COMM_RT, &
                         MPID%requestA(istep,1), ierr)

          ! Send JKQS if stimulated emission
          if (stm) &
            call MPI_ISEND(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                           MPI_DOUBLE_COMPLEX, &
                           MPID%lsend(istep), &
                           1000000+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,2), ierr)

        end do ! Sends

      ! Normal bcast
      else

        ! Share JKQ
        call MPI_BCAST(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                       MPI_DOUBLE_COMPLEX, 0, &
                       MPI_COMM_RT, ierr)

        ! Share JKQS if stimulated emission
        if (stm) &
          call MPI_BCAST(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

      end if ! Type of bcast

      ! If Polarization initial correction
      if (Pcorr) then

        !
        ! Make the J00 and J00S flat
        !

        ! For each height
        do iz=Rz0,Rz1

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Apply atomic shift
              jtran = itran + Atom(ia)%tshift

              ! For each FS transition within this transition
              do ftran=1,Atom(ia)%fst(itran)%nt

                ! Get the global transition index
                fftran = Atom(ia)%ifst_ij(ftran,itran)

                ! Apply atomic shift
                jftran = fftran + Atom(ia)%tfshift

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00(jftran,iz) = dble(JKQ(0,0,jtran,iz))

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00S(jftran,iz) = dble(JKQS(0,0,jtran,iz))

              end do ! fs transition
            end do ! term transition
          end do ! atom
        end do ! height


        !
        ! Solve SEE
        !

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,fftran,jftran,urou,umsg,tid) &
!$omp shared(nA,nz,Atom,Rho_old,J00,J00S,J00P,LambdaL,LambdaP) &
!$omp shared(Rz0,Rz1) &
!$omp shared(laborted,vaborted)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,tid,INPUT)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,tid)
#endif

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel

        !
        !  Calculate MRC
        !

        ! Only the master does
        if (pid.eq.0) then

          ! Calculate MRC
          call MRC_sb(Atom,Rho_old,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Write in stdout
          if (gpid.eq.0) then
            write(umsg,'(A,1x,es22.12)') ' - The conversion to '// &
                       'terms have changed the populations a '// &
                       'maximum of ',MRC%values(2,1)
            call verbose
          end if

        end if

        ! If Alternative broadcast
        if (MPID%nsend.gt.0) then

          ! For each slave
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing and reset the buffers
            call MPI_WAIT(MPID%requestA(iproc,1), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,2), &
                          MPI_STATUS_IGNORE,ierr)

          end do ! Processors

        end if ! Domain decomposition
      end if ! P correction


      !
      ! Convert to polarization quantities
      !

      !
      ! Stokes
      !

!$omp parallel default(none) &
!$omp private(iz,ifreq) &
!$omp shared(KSTK,Stokes,Stokes0,J00,J00S,JKQC,J00C,nz,nfreq,Geom) &
!$omp shared(Rz0,Rz1,axial,axiali)

      ! Allocate the extra Stokes parameters
!$omp single
      if (KSTK) then
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz1))
      else
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz0+1))
      end if
!$omp end single

      ! Initialize to 0
!$omp single
      if (KSTK) then
        if (axiali.and..not.axial) then
          do iph=1,Geom%nPh
            Stokes(0,:,iph,:,:) = Stokes0(:,1,:,:)
          end do
        else
          Stokes(0,:,:,:,:) = Stokes0
        end if
      end if
!$omp end single

      ! Get the intensity from the intensity array
!$omp workshare
      Stokes(0,:,:,:,:) = Stokes0
!$omp end workshare

!$omp single
      ! Deallocate the intensity array
      deallocate(Stokes0)

      ! Deallocate J00 and J00S
      deallocate(J00)
      deallocate(J00S)


      !
      ! JKQ frequency dependent
      !

      ! Allocate JKQC
      allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
!$omp end single
!$omp workshare
      JKQC = cZero
!$omp end workshare

      ! Get J00C into the complex array
!$omp do
      do iz=Rz0,Rz1
        do ifreq=1,nfreq

          JKQC(0,0,ifreq,iz) = dcmplx(J00C(ifreq,iz),0d0)

        end do ! frequencies
      end do ! heights
!$omp end do
!$omp end parallel

      ! Deallocate J00C
      deallocate(J00C)

      ! Control
      call control

      !
      ! Clean slave pointers
      !
2000  if (pid.ne.0) then
        deallocate(data2O)
        nullify(data2O)
      end if

      ! Put back the AD
      if (tbAD) AV = .False.

      return

      end subroutine JKQgen_alt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Does one NLTE iteration initializing with multilevel and only
      !! intensity data, in order to generate a consistent set of
      !! populations and mean intensities for the correspondent
      !! multiterm model atom, with one CPU (serial).\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            Pcorr(logical): Bool that says if the mean
      !!                            intensities are to be corrected\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!             rnPh(integer): Allocation size for Stokes\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine JKQgen_serial(Atom,Rho_old,Atmo,Frec,Geom, &
                               MPID,Input,Flgsg,Pcorr,Bfield,rnPh, &
                               Stokes0,J00,J00S,J00C, &
                               Stokes,JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Frequency_class):: Frec
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: Pcorr
      integer, intent(in):: rnPh
      ! Apparently, deferring allocatables keep the custom limits
      double precision, dimension(:,:,:,:), allocatable:: Stokes0
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      type(MRC_class):: MRC

      logical:: AD,ADD

      integer:: ith,iph,ia,iz,ifreq,itran,jtran,fftran,jftran,ftran
      integer:: if0,if1,nth,nph,tid

      double precision:: WA
      double precision, dimension(:,:), allocatable:: LambdaL
      double precision, dimension(:,:,:), allocatable:: LambdaP

      ! Pointers
      double precision, dimension(:,:), pointer:: data2O

      ! Trick to have AV input for AD calculation
      if (tbAD) AV = .True.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! Select the angular limits
      if (PRD.and.ADD) then
        nth = Geom%nTh
        nph = Geom%nPh
      else
        nth = 1
        nph = 1
      end if

      ! Initialize lambda operator to 0
      if0 = 1
      if1 = 1
      do ia=1,nA
        if (Atom(ia)%nftran.gt.if0) if0 = Atom(ia)%nftran
        if (Atom(ia)%nphot.gt.if1) if1 = Atom(ia)%nphot
      end do
      allocate(LambdaL(nxb,if0))
      LambdaL = 0d0
      allocate(LambdaP(nxb,if1,2))
      LambdaP = 0d0

      ! CPU limits
      if0 = 1
      if1 = nfreq


      !
      ! Allocation of buffers
      !

      ! Common (Master and slave)
      ! Allocate O pointers
      allocate(data2O(Frec%ntfreq,2))


      ! For each atom
      do ia=1,nA
        Rho_old(ia)%crho = Atom(ia)%crho
      end do

      !
      ! Allocate Radiation field
      !
      ! Allocate JKQ and JKQS and initialize to 0
      allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
      allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))


      ! Reset radiation field variables
      JKQ = cZero
      JKQS = cZero


      !  For each polar direction
      do ith=1,nTh

        ! For each azimuthal direction
        do iph=1,nPh

          ! Select the angular limits
          if (PRD.and.ADD) then
            WA = Geom%W_mu(ith)*Geom%W_mux(iph)
          else
            WA = 1d0
          end if

          ! For each height this CPU has assigned
          do iz=Rz0,Rz1

            ! RT coefficients
            call Termprof(Frec,Atom,Atmo,MPID,Flgsg,Geom, &
                          Bfield,iz,ith,iph,if0,if1,data2O)

            if (PRD.and.ADD) then
              if (axiali) then
                call Jgen(Atom,Frec%W_freq,WA, &
                          Stokes0(:,1,ith,iz),data2O, &
                          JKQ(:,:,:,iz),JKQS(:,:,:,iz))
              else
                call Jgen(Atom,Frec%W_freq,WA, &
                          Stokes0(:,iph,ith,iz),data2O, &
                          JKQ(:,:,:,iz),JKQS(:,:,:,iz))
              end if
            else
              call Jgen(Atom,Frec%W_freq,WA,J00C(:,iz),data2O, &
                        JKQ(:,:,:,iz),JKQS(:,:,:,iz))
            end if

          end do ! Intermedium heights
        enddo ! azimuthal angles
      enddo ! polar angles


      ! If Polarization initial correction
      if (Pcorr) then

        !
        ! Make the J00 and J00S flat
        !

        ! For each height
        do iz=Rz0,Rz1

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Apply atomic shift
              jtran = itran + Atom(ia)%tshift

              ! For each FS transition within this transition
              do ftran=1,Atom(ia)%fst(itran)%nt

                ! Get the global transition index
                fftran = Atom(ia)%ifst_ij(ftran,itran)

                ! Apply atomic shift
                jftran = fftran + Atom(ia)%tfshift

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00(jftran,iz) = dble(JKQ(0,0,jtran,iz))

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00S(jftran,iz) = dble(JKQS(0,0,jtran,iz))

              end do ! fs transition
            end do ! term transition
          end do ! atom
        end do ! height


        !
        ! Solve SEE
        !

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,fftran,jftran,umsg,urou,tid) &
!$omp shared(nA,nz,Atom,Rho_old,J00,J00S,J00P,LambdaL,LambdaP) &
!$omp shared(Rz0,Rz1) &
!$omp shared(laborted,vaborted)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,tid,INPUT)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,tid)
#endif

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel

        !
        !  Calculate MRC
        !

        ! Calculate MRC
        call MRC_sb(Atom,Rho_old,MRC)

        ! Convert cm into km
        MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

        ! Write in stdout
        if (gpid.eq.0) then
          write(umsg,'(A,1x,es22.12)') ' - The conversion to '//&
                     'terms have changed the populations a '// &
                     'maximum of ',MRC%values(2,1)
          call verbose
        end if

      end if ! P correction


      !
      ! Convert to polarization quantities
      !

      !
      ! Stokes
      !

!$omp parallel default(none) &
!$omp private(iz,ifreq) &
!$omp shared(KSTK,Stokes,Stokes0,J00,J00S,JKQC,J00C,nz,nfreq,Geom) &
!$omp shared(Rz0,Rz1,axial,axiali,data2O)

      ! Allocate the extra Stokes parameters
!$omp single
      if (KSTK) then
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz1))
      else
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz0+1))
      end if
!$omp end single

      ! Initialize to 0
!$omp workshare
      Stokes = 0d0
!$omp end workshare

      ! Get the intensity from the intensity array
!$omp single
      if (KSTK) then
        if (axiali.and..not.axial) then
          do iph=1,Geom%nPh
            Stokes(0,:,iph,:,:) = Stokes0(:,1,:,:)
          end do
        else
          Stokes(0,:,:,:,:) = Stokes0
        end if
      end if
!$omp end single

!$omp single
      ! Deallocate the intensity array
      deallocate(Stokes0)

      ! Deallocate J00 and J00S
      deallocate(J00)
      deallocate(J00S)

      ! Clean slave pointers
      deallocate(data2O)
      nullify(data2O)


      !
      ! JKQ frequency dependent
      !

      ! Allocate JKQC
      allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
!$omp end single
!$omp workshare
      JKQC = cZero
!$omp end workshare

      ! Get J00C into the complex array
!$omp do
      do iz=Rz0,Rz1
        do ifreq=1,nfreq

          JKQC(0,0,ifreq,iz) = dcmplx(J00C(ifreq,iz),0d0)

        end do ! frequencies
      end do ! heights
!$omp end do
!$omp end parallel

      ! Deallocate J00C
      deallocate(J00C)

      ! Put back the AD
      if (tbAD) AV = .False.

      ! Control
      call control

      return

      end subroutine JKQgen_serial

!#####################################################################
!#####################################################################
!#####################################################################

      end module solveri_mod
