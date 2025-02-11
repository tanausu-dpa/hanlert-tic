      !> Solves the NLTE problem of the second kind
      module solver_mod
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
!     09/23/2024 V3.0.17
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/23/2024:   V3.0.17 - Changed MPI tags to comply with the
!                             standard (TdPA)
!
!     08/08/2024:   V3.0.16 - Force at least two iterations if doing
!                             AD starting from AA (TdPA)
!
!     02/20/2024:   V3.0.15 - Bugfix: Wrong initial index in the
!                             tau1 calculation when restricting
!                             heights (TdPA)
!
!     11/14/2023:   V3.0.14 - Call scattering_manage in angular
!                             loops in solvers (TdPA)
!                           - Call get_scattering_los when computing
!                             emergent LOS (TdPA)
!                           - Added termination paths in case of
!                             error (TdPA)
!
!     10/16/2023:   V3.0.13 - Made LTElines allocatable to satisfy
!                             memory warnings (TdPA)
!
!     09/29/2023:   V3.0.12 - Added calculation of radiation field
!                             tensors with Q<0 from conjugation rules
!                             to avoid computing them (TdPA)
!
!     08/17/2023:   V3.0.11 - When aborting, go to the memory
!                             cleaning instead of returning (TdPA)
!                           - Moved control call before memory
!                             cleaning (TdPA)
!
!     08/07/2023:   V3.0.10 - Added arguments for LTE lines (TdPA)
!                           - Added solve and emergent (TdPA)
!
!     07/03/2023:    V3.0.9 - chi500 is now it the Atmo structure and
!                             not in Cont (TdPA)
!                           - In the inversion, the result of
!                             the emergence is stored in the SolF
!                             structure, setup by calling settau,
!                             setstk, and setctr (TdPA)
!
!     05/24/2023:    V3.0.8 - Bugfix: Wrong argument in RTContr lead
!                             to gibberish in the output contribution
!                             function (TdPA)
!
!     03/23/2023:    V3.0.7 - Bugfix: Ensured the OpenMP version
!                             compiles after some years of changes,
!                             albeit did not test it works (TdPA)
!
!     02/14/2023:    V3.0.6 - Also write Stokes file in inversion
!                             mode (HL)
!
!     11/10/2022:    V3.0.5 - Added JKQ_asymm as argument in the calls
!                             to RTCoeff and RTCoeffe (TdPA)
!                           - Added JKQ_asymm is an input to emergence
!                             and emergence_serial to pass them to
!                             RTCoeffe (TdPA)
!                           - Bugfix: In solver_serial, the call to
!                             addJKQasym the JKQ in the arguments
!                             must be the new ones, not the old ones,
!                             JKQ -> JKQ_n (TdPA)
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
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!                           - funit is now global (TdPA)
!
!     06/30/2022:    V3.0.1 - Bugfixes: Quick changes while preparing
!                             the last commit led to mistakes. Solved
!                             a typo result of copy and paste (TdPA)
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
!                              o The geometrical tensors for the
!                                emergent solution are defined just
!                                before starting the formal
!                                solution.
!                             (TdPA)
!
!     04/07/2022:    V2.0.2 - Bugfix: There was a frequency loop
!                             in the serial version of the solver
!                             with polarization. This means that
!                             during the formal solution the Stokes
!                             parameters at the last point in each
!                             direction were wrong (TdPA)
!                           - Added an argument to writesol (TdPA)
!
!     03/23/2021:    V2.0.1 - Bugfix: Was not possible to append in
!                             the MRC file (TdPA)
!                           - Made some changes to improve error
!                             handling (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Completely changed the Master section
!                             and introduced OpenMP (TdPA)
!
!     02/12/2021:    V1.8.9 - Call report_mpi_time in the solver
!                             subroutines (TdPA)
!
!     01/13/2021:    V1.8.8 - Call the subroutine to add ad-hoc values
!                             to the JKQ tensors from the master in
!                             the solvers (TdPA)
!
!     01/12/2021:    V1.8.7 - Using KSTK to determine how to store
!                             the Stokes parameters (TdPA)
!
!     09/11/2020:    V1.8.6 - Changes to avoid unallocated arrays
!                             when the photoionization exponential
!                             is not stored, either because it is
!                             specified or because there are no
!                             photoionizations (TdPA)
!
!     07/31/2020:    V1.8.5 - Change the p_exu pointer so it only
!                             points to the stored data if there
!                             is something to point to (TdPA)
!
!     07/10/2020:    V1.8.4 - tau1 needs to be allocated in the slaves
!                             also when computing only the
!                             contribution function (TdPA)
!
!     05/11/2020:    V1.8.3 - Bugfix: When loading a solution file
!                             with AA redistribution in order to
!                             compute with AD redistribution, during
!                             the first iteration (which assumes AA)
!                             the RAM storage of the redistribution
!                             must be switched off because dimensions
!                             do not coincide (TdPA)
!
!     01/16/2020:    V1.8.2 - Bugfix: The wrong variable was checked
!                             to decide if exponentials were stored
!                             in RAM (TdPA)
!
!     11/19/2019:    V1.8.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     09/13/2019:    V1.8.0 - Now Stokes needs to be stored if the
!                             calculation if dynamic and angle-
!                             averaged too (TdPA)
!                           - The vertical scale can now be optical
!                             depth (TdPA)
!
!     05/31/2019:    V1.7.0 - Changed the dimensionality of the
!                             profile variable. Now it runs
!                             sequentially on atoms, transitions and
!                             frequencies to save memory and reduce
!                             the size of data shared through MPI
!                             messages (TdPA)
!                           - Now the emerging routines use RTCoeffe
!                             instead of RTCoeff. It is not necessary
!                             to define emerging anymore (TdPA)
!
!     05/08/2019:    V1.6.5 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!                           - Deals with parameters to allow for bad
!                             physics (TdPA)
!
!     04/17/2019:    V1.6.4 - Implemented back the broadcasting from
!                             the mpi libraries (TdPA)
!
!     03/20/2019:    V1.6.3 - Now resets the MRC file (TdPA)
!
!     03/18/2019:    V1.6.2 - Added logic to allow for precomputed
!                             exponentials for J integrals (TdPA)
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
!     09/25/2018:    V1.5.0 - Added NG acceleration (TdPA)
!
!     09/04/2018:    V1.4.0 - Added solver_alt (TdPA)
!
!     08/06/2018:    V1.3.6 - Removed J' from the rho00 index of the
!                             MRC in the output, because rho00 is
!                             always diagonal (TdPA)
!
!     05/16/2018:    V1.3.5 - Changed azimuthal dimension of Stokes
!                             from nPh2 to nPh (TdPA)
!
!     11/01/2017:    V1.3.2 - Changed legend message for iteration
!                             information (TdPA)
!
!     10/03/2017:    V1.3.1 - Added Bfield to RTAbs (TdPA)
!                           - Moved the exit of the iteration loop
!                             after the request checking (TdPA)
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
!     08/28/2017:    V1.1.3 - Emergence now waits for every process
!                             to finish each direction before letting
!                             any of them to start the next one. This
!                             is very inefficient for the domain
!                             decomposition (TdPA)
!
!     08/22/2017:    V1.1.2 - Some improvements in emergence and
!                             non-code-breaking bugs were fixes (TdPA)
!
!     08/21/2017:    V1.1.1 - Erased temporal lines from V1.1.0 (TdPA)
!                           - Avoiding making slaves to wait to
!                             communicate between them (TdPA)
!                           - Each processor starts with the direction
!                             that makes them closer to the boundary
!                             (TdPA)
!                           - Applied the logic in change of version
!                             1.1.0 to the contribution function in
!                             emergence (TdPA)
!
!     08/17/2017:    V1.1.0 - Increased the buffer size, so each
!                             process does never wait the master to
!                             receive (TdPA)
!
!     07/10/2017:   V1.0.10 - Bugfix: In emergence, tauM needed a
!                             shift in index because it is a pointer,
!                             also, there was a sender of contribution
!                             function that was sending the variable
!                             instead of the sender, what could result
!                             in changes of the buffer before is
!                             received (TdPA)
!
!     07/06/2017:    V1.0.9 - Bugfix: In emergence_serial, the tau
!                             was calculated in the inverse direction,
!                             that is, making the lower boundary 0 and
!                             counting outwards (TdPA)
!
!     06/28/2017:    V1.0.8 - Receiving and passing Red (TdPA)
!
!     06/20/2017:    V1.0.7 - Changed data1, data2, rLine and rPhot
!                             arrays into two pointers each (three
!                             for data1), this avoids copies in
!                             exchange of more allocation/
!                             deallocation (TdPA)
!                           - Changed etaIM, etaIO, and tauM in
!                             emergence to pointers too (TdPA)
!                           - Changed many MPI%ifx to ifx (TdPA)
!
!     06/16/2017:    V1.0.6 - Made the modifications consistent with
!                             the change in the Frec limits (TdPA)
!
!     06/15/2017:    V1.0.5 - Bugfix: Missing AD condition when
!                             determining the Stokes storage for the
!                             master in solver (TdPA)
!
!     06/13/2017:    V1.0.4 - Added indexes for Fint and Jcalc calls,
!                             consistent with the changes in those
!                             routines (TdPA)
!
!     06/12/2017:    V1.0.3 - The frequency limits are passed to the
!                             RTcoeffi routines (TdPA)
!                           - Removed cont%l
!
!     05/05/2017:    V1.0.2 - RTStep now knows which point is doing
!                             during the call (TdPA)
!
!     05/04/2017:    V1.0.1 - There was no wait after sending the
!                             radiative data to the last CPU, the
!                             master had chances of zeroing the
!                             buffer (TdPA)
!                           - Avoiding the use of status arguments
!                             in MPI calls (TdPA)
!
!     04/20/2017:    V1.0.0 - First Version (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
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
!  solve:
!    Manage what solver must be called depending on the MPI
!  configuration
!
!  solver:
!    This subroutine creates and solves the radiation transfer problem
!  of the second kind
!
!  solver_alt:
!    This subroutine creates and solves the radiation transfer problem
!  of the second kind with the alternative and slower MPI
!  implementation
!
!  solver_serial:
!    This subroutine creates and solves the radiation transfer problem
!  of the second kind, with just one processor
!
!  emergent:
!    Manage what emergence must be called depending on the MPI
!  configuration
!
!  emergence:
!    This subroutine solves the radiation transfer equation for the
!  specified directions given the rhoKQ and JKQ
!
!  emergence_serial:
!    This subroutine solves the radiation transfer equation for the
!  specified directions given the rhoKQ and JKQ, with just one
!  processor
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use boundary_mod
      use commons_mod
      use gauss_mod , only : setTKQLOS
      use iosolution_mod
      use jcalc_mod
      use mrc_mod
      use ng_mod
      use omp_mod
      use parameters_mod , only : B2L , cZero , kb, cSaha , fktoJ
      use rtcoeff_mod
      use rtstep_mod
      use see_mod
      use setmpi_mod
      use types_mod

      ! Maximum buffer for NG_int
      double precision, parameter:: maxbuffer_NG = 500d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the solver for the NLTE problem of the second kind.\n
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
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      subroutine solve(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                       Bfield,Geom,MPID,Input,Flgsg,JKQ_asym, &
                       Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                      giz0:giz1), target:: Stokes
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,-2)
#endif
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,-2)
#endif
      ! MPI
      if (MPID%mpi) then

        if(MPID%alternP)then
          call solver_alt(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                          Red,Bfield,Geom,MPID,Input,Flgsg, &
                          JKQ_asym,Stokes,JKQ,JKQS,JKQC)
        else
          call solver(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                      Bfield,Geom,MPID,Input,Flgsg,JKQ_asym, &
                      Stokes,JKQ,JKQS,JKQC)
        end if

      ! Serial
      else

        call solver_serial(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                           Red,Bfield,Geom,MPID,Input,Flgsg, &
                           JKQ_asym,Stokes,JKQ,JKQS,JKQC)
      end if
#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,-1)
#endif
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,-1)
#endif

      end subroutine solve

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem of the second kind with lambda
      !! iteration with several CPU.\n
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
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      subroutine solver(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                        Bfield,Geom,MPID,Input,Flgsg,JKQ_asym, &
                        Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                      giz0:giz1), target:: Stokes
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

      ! Local

      type(MRC_class):: MRC

      character(LEN=20):: iterS

      logical:: goout, AD, ADD, doNG, laux, RPRAM, deal

      integer:: iaux,ios,iz0,iz1,diz,m,o,p,op,if0,if1,iphot
      integer:: iproc,iter,ith,iph,ia,iz,ifreq,itran,id,if0l,if1l,nfl
      integer:: ntpz,npz,ntp,itpz,tid,K,iQ
#ifdef _OPENMP
      integer:: itz
#endif
      integer:: istep,ip0l,ip1l,jtran,jphot
      integer:: NG_dim,NG_entry,iterm,iJ,ing,nftl

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset
      double precision:: daux,WA,mu_inv,dsm,dsp,larmor
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, &
             dimension(2,nxtran,Geom%nph,Geom%nth,Rz0:Rz1):: Norm
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P

      complex(kind=8), &
         dimension(0:2,0:2,2,nxtran,Geom%nph,Geom%nth,Rz0:Rz1):: BStk

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Stokes_r
      double precision, dimension(:), allocatable, target:: Prof_r
      ! Senders
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes_s
      double precision, dimension(:,:,:,:,:), allocatable:: Prof_s
      ! Duals
      integer:: info_b

      ! Pointers
      double precision, dimension(:,:,:), pointer:: data1M,data1O, &
                                                    data1P
      double precision, dimension(:,:), pointer:: data2O,data2P
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP
      double precision, dimension(:), pointer:: p_exu
      double precision, dimension(:,:,:), pointer:: p_MStk
      double precision, dimension(:,:,:), pointer:: p_MProf


      ! Routine name
      urou = 'solver'

      ! Initialize converged flag
      goout = .False.

      ! Initialize angle depended flag
      AD = .not.AV

      ! Angle dependent or dynamic
      ADD = AD.or.dyn

      ! If storing redistribution
      RPRAM = PRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AV = .True.
        PRAM = .False.
      end if

      ! Initialize index of Stokes
      op = 1


      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.

      ! If NG acceleration of rho00
      if (Input%NG) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + RnZ*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRD) then

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

          Input%NG = .False.

          if (pid.eq.0) then
            umsg = ' # The buffer for NG acceleration '// &
                   'is too big. Not doing NG.'
            call verbose
          end if

        end if ! Buffer

        ! If finally doing it, allocate
        if (Input%NG) then

          ! Master
          if (pid.eq.0) then

            allocate(NG_scratch(NG_dim, Input%NG_ord+2))

          ! Slaves
          else

            allocate(NG_scratch(NG_dim, 1))

          end if

        end if ! NG and master
      end if ! NG

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
        iaux = 4*MPID%nxfreq*Geom%nTh*Geom%nPh*Rnz
        allocate(Stokes_r(iaux))

        ! To receive profile information
        iaux = MPID%nxtfreq*2*Geom%nTh*Geom%nPh*Rnz
        allocate(Prof_r(iaux))

      ! Slave
      else

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(0:3,if0:if1,0:5))
        allocate(data1O(0:3,if0:if1,0:5))

        ! Allocate O pointers
        allocate(data2O(Frec%ntfreq,2))

        ! To send Intensity chunks
        allocate(Stokes_s(0:3,if0:if1,Rz0:Rz1,Geom%nPh,Geom%nTh))

        ! To send profile information
        allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1,Geom%nPh,Geom%nTh))

      end if ! Master or Slave


      !
      ! Initialization messages
      !

      ! Announce that we are starting
      if (gpid.eq.0) then
        umsg = '   Iteration          MRC(rho^0_0) Atom_index '// &
               'Term_index    2J  Height_index Height(km) |   '// &
               '       MRC(rho^K_Q) Atom_index Term_index '// &
               "   2J   2J' Height_index Height(km)  K  Q"
        call verbose
      end if


      ! Open the file to store MRC
      if (gpid.eq.0.and.Input%keep_MRC) then
        ! If appending
        if (Input%appendMRC) then
          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRC', exist=laux)
          if (.not.laux) then
            open(800, file=trim(Input%folder)//'/MRC', &
                 action='write',iostat=ios,err=1000)
            write(800,'(A)') '!  Iteration          MRC(rho^0_0) '// &
                         'Atom_index Term_index    2J  '// &
                         'Height_index Height(km) |   '// &
                         '       MRC(rho^K_Q) Atom_index '// &
                         'Term_index    2J   2J'// &
                         "' Height_index Height(km)  K  Q"
            close(800)
          end if
        ! If not appending, create
        else
          open(800, file=trim(Input%folder)//'/MRC', &
               action='write',iostat=ios,err=1000)
          write(800,'(A)') '!  Iteration          MRC(rho^0_0) '// &
                       'Atom_index Term_index    2J  '// &
                       'Height_index Height(km) |   '// &
                       '       MRC(rho^K_Q) Atom_index '// &
                       'Term_index    2J   2J'// &
                       "' Height_index Height(km)  K  Q"
          close(800)
        end if
      end if

      ! Control
      call control
      if (laborted) goto 2000

      ! If measuring performance
      if (Input%mpi_perf.and.pid.eq.0) &
        call report_mpi_time(Input%folder,Input%ID, &
                             0,0,.False.)


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=Input%iter_min,Input%iter_max

        ! Flags for physics
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics
        if (iter.le.Input%allownphys_rho) then
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


        !
        ! Master
        !
        if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(id,Stokes_r,Prof_r,itpz,ith,iph,iz,deal,WA) &
!$omp private(itran,jtran,ia,daux,dsm,p_exu,itz) &
!$omp shared(JKQ,JKQS,J00P,Norm,BStk,JKQC,iter,p_MStk,p_MProf) &
!$omp shared(if0l,if1l,ip0l,ip1l,nfl,nftl,ntpz,npz,nz,na,Stokes) &
!$omp shared(Rz0,Rz1,Rnz,info_b,KSTK,stm,PIRAM,laborted) &
!$omp shared(MPID,Input,Frec,Atom,Geom,Atmo,ierr) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT)

          ! Reset radiation field variables
!$omp workshare
          JKQ = cZero
          JKQS = cZero
          J00P = 0d0
          Norm = 0d0
          BStk = cZero
          JKQC = cZero
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
              call MPI_recv(Stokes_r(1), MPID%size4(info_b), &
                            MPI_DOUBLE_PRECISION, info_b, &
                            info_b, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! Receive profile
            do while (.True.)
              call MPI_recv(Prof_r(1), MPID%size5(info_b), &
                            MPI_DOUBLE_PRECISION, info_b, &
                            1+info_b, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! If measuring performance
            if (Input%mpi_perf) &
              call report_mpi_time(Input%folder,Input%ID, &
                                   info_b,iter,.True.)

            ! Shorter variables
            if0l = MPID%if0(info_b)
            if1l = MPID%if1(info_b)
            ip0l = Frec%Mpif0(info_b)
            ip1l = Frec%Mpif1(info_b)
            nfl = MPID%nf(info_b)
            nftl = Frec%Mntfreq(info_b)

            ! Pointers
            p_MStk(0:3,if0l:if1l,1:ntpz) => &
                                       Stokes_r(1:MPID%size4(info_b))
            p_MProf(1:nftl,1:2,1:ntpz) => Prof_r(1:MPID%size5(info_b))
!$omp end single

            ! Do not process if leaving
            if (laborted) cycle

            deal = .False.

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
                Stokes(0:3,if0l:if1l,iph,ith,iz) = &
                                            p_MStk(0:3,if0l:if1l,itpz)

              ! Calculate frequency integral for b-b quantities
              call FInt_line(Atom,MPID,Geom,Frec%W_freq, &
                             Frec%Mlif0(info_b),Frec%Mlif1(info_b), &
                             info_b,iph,ith,iz,p_MStk(:,:,itpz), &
                             p_MProf(:,:,itpz), &
                             Norm(:,:,iph,ith,iz), &
                             BStk(:,:,:,:,iph,ith,iz))

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
                  call FInt_rest(Atom,MPID,Geom,Frec%omega, &
                                 Frec%W_freq,ip0l,ip1l, &
                                 Atmo%T(iz),info_b,iph,ith,WA, &
                                 p_MStk(:,:,itpz),J00P(:,:,iz), &
                                 JKQC(:,:,:,iz),p_exu)
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
            nullify(p_MStk,p_MProf)
!$omp end single

          !
          ! Apply weights to JKQ, JKQS, and normalize
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
                  do itran=1,Atom(ia)%ntran

                    ! Apply atomic shift
                    jtran = itran + Atom(ia)%tshift

                    ! Get the weight
                    if (Norm(1,jtran,iph,ith,iz).gt.0d0) then

                      daux = WA/Norm(1,jtran,iph,ith,iz)

                      ! Integrate angle
                      JKQ(0:2,:,jtran,iz) = JKQ(0:2,:,jtran,iz) + &
                                     BStk(:,:,1,jtran,iph,ith,iz)*daux

                    end if

                    ! If there is stimulated emission
                    if (stm) then

                      ! Get the weight
                      if (Norm(2,jtran,iph,ith,iz).gt.0d0) then

                        daux = WA/Norm(2,jtran,iph,ith,iz)

                        ! Integrate angle
                        JKQS(0:2,:,jtran,iz) = &
                                     JKQS(0:2,:,jtran,iz) + &
                                     BStk(:,:,2,jtran,iph,ith,iz)*daux

                      end if

                    end if ! Stimulated emission
                  end do ! transitions
                end do ! atoms
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

              end do ! b-f transitions
            end do ! atoms
          end do ! heights
!$omp end do
!$omp end parallel

          ! If not axial, need to complete JKQ
          if (.not.axial) then

            ! Continuum
            do K=1,Krad
              do iQ=1,K
                JKQC(-iQ,K,:,:) = Flgsg%sg(iQ)*conjg(JKQC(iQ,K,:,:))
              end do
            end do

            ! For each atom
            do ia=1,nA

              ! For each FS transition
              do itran=1,Atom(ia)%ntran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tshift

                ! K
                do K=1,Atom(ia)%Krad(itran)

                  ! Q
                  do iQ=1,K

                    ! Conjugate
                    JKQ(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                         conjg(JKQ(iQ,K,jtran,:))

                    ! Stimulated, conjugate
                    if (stm) &
                      JKQS(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                            conjg(JKQS(iQ,K,jtran,:))

                  end do ! Q
                end do ! K
              end do ! Transition
            end do ! Atom
          end if ! Not axial

          !
          ! Add the Ad-hoc asymmetries
          !
          if (Input%nasym.gt.0) &
            call addJKQasym(Bfield,Flgsg,JKQ_asym,JKQ,JKQS,JKQC)

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

            ! Determine the first and last height indexes to run over
            iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
            iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

            ! For each azimuthal direction
            do iph=1,Geom%nPh

              !
              ! If angle-dependent, manage scattering angles
              if (.not.AV.and.PRD) &
                call scattering_manage(Geom,ith,iph)

              !
              ! First height
              !

              ! If going down, get top boundary
              if(diz.eq.1)then

                ! Call top boundary
                call top(MPID,data1M(:,:,5))

              ! If going up, get bottom boundary
              else

                ! Call bottom boundary
                call bottom(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                            Atmo%vy(iz0),Atmo%vz(iz0), &
                            Geom%V_mu(ith),Geom%V_mux(iph), &
                            Geom%V_muy(iph),MPID,data1M(:,:,5))

              endif ! propagation direction

              ! Identify current height
              o = iz0

              ! Index for Stokes
              if (PRD.and.ADD) op = o

              ! Calculate radiative coefficients
              call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                           Geom,o,ith,iph,if0,if1,JKQ_asym, &
                           JKQ(:,:,:,o),JKQC(:,:,:,o),Cont%ndir, &
                           Cont%c(:,:,:,o),Bfield, &
                           Stokes(:,:,:,:,op), data1M(:,:,0:4), &
                           data2O)
              if (laborted) goto 3000

              !
              ! Store in buffer
              !

              ! Stokes
              Stokes_s(:,:,o,iph,ith) = data1M(:,:,5)

              ! Profiles
              Prof_s(:,:,o,iph,ith) = data2O

              ! Identify next height
              p = iz0 + diz

              ! Index for Stokes
              if (PRD.and.ADD) op = p

              ! Calculate radiative coefficients
              call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                           Geom,p,ith,iph,if0,if1,JKQ_asym, &
                           JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                           Cont%c(:,:,:,p),Bfield, &
                           Stokes(:,:,:,:,op),data1O(:,:,0:4), &
                           data2O)
              if (laborted) goto 3000


              !
              ! Intermedium heights
              !

              ! For each height this CPU has assigned
              do iz=iz0,iz1,diz

                ! We treat the boundaries outside
                if(iz.eq.iz0.or.iz.eq.iz1)cycle

                ! Allocate P pointers
                allocate(data1P(0:3,if0:if1,0:5))
                allocate(data2P(Frec%ntfreq,2))

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
                if (PRD.and.ADD) op = p

                ! RT coefficients
                call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                             Geom,p,ith,iph,if0,if1,JKQ_asym, &
                             JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                             Cont%c(:,:,:,p),Bfield, &
                             Stokes(:,:,:,:,p),data1P(:,:,0:4), &
                             data2P)

                ! Point to the data
                p_K0M  => data1M(:,:,0)
                p_K1M  => data1M(:,:,1)
                p_K2M  => data1M(:,:,2)
                p_SM   => data1M(:,:,4)
                p_StkM => data1M(:,:,5)
                p_K0O  => data1O(:,:,0)
                p_K1O  => data1O(:,:,1)
                p_K2O  => data1O(:,:,2)
                p_SO   => data1O(:,:,4)
                p_StkO => data1O(:,:,5)
                p_K0P  => data1P(:,:,0)
                p_SP   => data1P(:,:,4)

                ! Apply short characteristics BESSER
                call RTStep(o,ith,iph,MPID%nf(pid), &
                            dsm,dsp,p_K0M,p_K1M,p_K2M, &
                            p_SM,p_K0O,p_K1O,p_K2O, &
                            p_SO,p_K0P,p_SP,p_StkM, &
                            p_StkO,.True.)

                !
                ! Store in buffer
                !

                ! Stokes
                Stokes_s(:,:,o,iph,ith) = data1O(:,:,5)

                ! Profiles
                Prof_s(:,:,o,iph,ith) = data2O

                ! Shift data (O->M, P->O)
                deallocate(data1M,data2O)
                data1M => data1O
                data1O => data1P
                data2O => data2P
                nullify(data1P,data2P)

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
              p_K0M  => data1M(:,:,0)
              p_K1M  => data1M(:,:,1)
              p_K2M  => data1M(:,:,2)
              p_SM   => data1M(:,:,4)
              p_StkM => data1M(:,:,5)
              p_K0O  => data1O(:,:,0)
              p_K1O  => data1O(:,:,1)
              p_K2O  => data1O(:,:,2)
              p_SO   => data1O(:,:,4)
              p_StkO => data1O(:,:,5)

              ! Apply short characteristics LINEAR
              call RTStep(o,ith,iph,MPID%nf(pid), &
                          dsm,dsp,p_K0M,p_K1M,p_K2M, &
                          p_SM,p_K0O,p_K1O,p_K2O, &
                          p_SO,p_K0P,p_SP,p_StkM, &
                          p_StkO,.False.)
              if (laborted) goto 3000

              !
              ! Store in buffer
              !

              ! Stokes
              Stokes_s(:,:,o,iph,ith) = data1O(:,:,5)

              ! Profiles
              Prof_s(:,:,o,iph,ith) = data2O

            enddo ! azimuthal angles
          enddo ! polar angles


          !
          ! Send to master
          !

          ! If had an error
3000      if (laborted) then

            ! Send error
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

            ! Send Stokes
            do while (.True.)
              call MPI_SEND(Stokes_s(0,if0,Rz0,1,1), &
                            MPID%size4(pid), &
                            MPI_DOUBLE_PRECISION, 0, pid, &
                            MPI_COMM_RT, ierr)
              if (ierr.eq.0) exit
            end do

            ! Send profiles
            do while (.True.)
              call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                            MPID%size5(pid), MPI_DOUBLE_PRECISION, &
                            0, 1+pid, MPI_COMM_RT, ierr)
              if (ierr.eq.0) exit
            end do

          end if ! Error

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
                          MPID%recv, 1+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Receive JKQC
            call MPI_RECV(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                          MPI_DOUBLE_COMPLEX, &
                          MPID%recv, 2+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Receive Stokes if doing A-D PRD
            if (PRD.and.ADD) &
            call MPI_RECV(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                          MPI_DOUBLE_PRECISION, &
                          MPID%recv, 3+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Bradcast J00 for b-f transitions
            call MPI_RECV(J00P(1,1,Rz0), MPID%sizei3(0), &
                          MPI_DOUBLE_PRECISION, &
                          MPID%recv, 4+pid, &
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
                           1+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,2), ierr)

            ! Send JKQC
            call MPI_ISEND(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                           MPI_DOUBLE_COMPLEX, &
                           MPID%lsend(istep), &
                           2+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,3), ierr)

            ! Send Stokes if doing A-D PRD
            if (PRD.and.ADD) &
            call MPI_ISEND(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                           MPI_DOUBLE_PRECISION, &
                           MPID%lsend(istep), &
                           3+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,4), ierr)

            ! Send J00 for b-f transitions
            call MPI_ISEND(J00P(1,1,Rz0), MPID%sizei3(0), &
                           MPI_DOUBLE_PRECISION, &
                           MPID%lsend(istep), &
                           4+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,5), ierr)

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

          ! Share JKQC
          call MPI_BCAST(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

          ! Share Stokes if doing A-D PRD
          if (PRD.and.ADD) &
          call MPI_BCAST(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Share J00 for b-f transitions
          call MPI_BCAST(J00P(1,1,Rz0), MPID%sizei3(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Type of bcast


        !
        ! Solve SEE
        !
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,iter)
#endif

        ! For each atom
!$omp parallel default(none) if (iter > Input%iter_min) &
!$omp private(ia,iz,itran,jtran,iphot,jphot,larmor,tid) &
!$omp private(urou,umsg) &
!$omp shared(nA,nz,Atom,JKQ,JKQS,J00P,Flgsg,Bfield,Rz0,Rz1) &
!$omp shared(laborted,vaborted)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tshift + 1
          jtran = itran + Atom(ia)%ntran - 1
          iphot = Atom(ia)%pshift + 1
          jphot = iphot + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Set magnetic data
            larmor = B2L*Bfield%Bstrength(iz)

            ! Solve the SEE
            call SEE(Atom(ia),JKQ(:,:,itran:jtran,iz), &
                     JKQS(:,:,itran:jtran,iz), &
                     J00P(iphot:jphot,:,iz),larmor,Flgsg,iz,tid)

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel
#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,iter)
#endif

        ! Control
        call control
        if (laborted) goto 2000


        !
        ! NG acceleration
        !

        ! Check if doing NG acceleration
        if(Input%NG.and.iter.gt.Input%NG_delay)then

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
            if (PRD) then

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
                                            Stokes(0,ifreq,iph,ith,iz)

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
                    NG_scratch(o,NG_entry) = dble(JKQC(0,0,ifreq,iz))

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD
            end if ! PRD

            ! Call NG and check if it should be processed
            call NG(NG_dim,Input%NG_ord,NG_scratch,NG_entry,doNG)

          ! Slave
          else

            ! If wrong order
            if (Input%NG_ord.lt.1.or.Input%NG_ord.gt.5) then
              doNG = .False.
            ! Valid order
            else
              ! Check if Master is in NG step
              if (NG_entry.gt.(Input%NG_ord+1)) then
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
                              MPID%recv, 7+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No master

              ! For each send
              do istep=1,MPID%nsend

                ! Send NG iteration
                call MPI_ISEND(NG_scratch(1,ing), NG_dim, &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               7+MPID%lsend(istep), &
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
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Advance index
                    o = o + 1

                    ! Accelerate rho
                    Atom(ia)%crho(p,iz) = &
                                        dcmplx(NG_scratch(o,ing), 0d0)
                  enddo ! Heights
                enddo ! Levels
              enddo ! Terms
            enddo ! Atoms

            ! If PRD
            if (PRD) then

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
                        Stokes(0,ifreq,iph,ith,iz) = NG_scratch(o,ing)

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
                    JKQC(0,0,ifreq,iz) = &
                                        dcmplx(NG_scratch(o,ing), 0d0)
                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD
            end if ! PRD

            if (gpid.eq.0) then
              umsg = 'NG acceleration'
              call verbose
            end if

            NG_entry = 0

            ! Master wait to free the buffer and alternative bcast
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
          call MRC_sb(Atom,Rho_old,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5
          MRC%values(1,2) = Atmo%z(MRC%indexes(2,2))*1d-5

          ! Check exit criteria
          if (MRC%values(2,1).le.Input%mrc_i.and. &
              MRC%values(2,2).le.Input%mrc_p) &
            goout = .True.

        end if

        ! Only the global Master
        if (gpid.eq.0) then

          ! Write in stdout
          write(umsg, &
          '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
          '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,1x,i4,3x,i11,'// &
          '2x,f9.3,1x,i2,1x,i2)') &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(5,1),MRC%indexes(2,1), &
          MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
          MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
          MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
          MRC%indexes(6,2)
          call verbose

          ! File
          if (Input%keep_MRC) then

            ! Write in MRC file
            open(800, file=trim(Input%folder)//'/MRC', &
                 iostat=ios,err=1000,position='append')
            write(800, &
            '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
            '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,2x,i4,2x,i11,'// &
            '2x,f9.3,1x,i2,1x,i2)', err=1100) &
            iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
            MRC%indexes(5,1),MRC%indexes(2,1), &
            MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
            MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
            MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
            MRC%indexes(6,2)
            close(800)

          end if
        end if

        ! We can swith now to AD if we had AV input
        if (tbAD) then
          AV = .False.
          tbAD = .False.
          PRAM = RPRAM
          ! Force at least two iterations
          goout = .False.
        end if


        !
        ! Save partial solution
        !

        ! Check if we want to store partial results
        if(Input%store.and.mod(iter,Input%store_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesol(Input,iterS,Frec%omega,Geom,Flgsg, &
                        Bfield,Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)

        ! Or have a control check
        else

          call control
          if (laborted) goto 2000

        end if

        !
        ! Share if we are finished or not

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive goout
            call MPI_RECV(goout, 1, MPI_LOGICAL, &
                          MPID%recv, 5+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, ierr)

          end if ! No Master

          ! For each send
          do istep=1,MPID%nsend

            ! Send goout
            call MPI_ISEND(goout, 1, MPI_LOGICAL, &
                           MPID%lsend(istep), &
                           5+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,6), &
                           ierr)

          end do ! sends

        ! Normal bcast
        else

          ! Share goout
          call MPI_BCAST(goout, 1, MPI_LOGICAL, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Type of bcast

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

          end do ! Processors

        end if ! Domain decomposition

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
        if (ierr.ne.0) goto 1000

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = real(MRC%values(2,2))

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

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
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO)
      end if

      return

1000  umsg = 'Error opening MRC file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing MRC file'
      close(800)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine solver

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem of the second kind with lambda
      !! iteration with several CPU. Alternative communication.\n
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
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      subroutine solver_alt(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                            Red,Bfield,Geom,MPID,Input,Flgsg, &
                            JKQ_asym,Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                      giz0:giz1), target:: Stokes
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

      ! Local

      type(MRC_class):: MRC

      character(LEN=20):: iterS

      logical:: goout, AD, ADD, doNG, laux, RPRAM, deal

      integer:: iaux,ios,iz0,iz1,diz,m,o,p,op,if0,if1,iphot,K,iQ
      integer:: iproc,iter,ith,iph,ia,iz,ifreq,itran,id,if0l,if1l,nfl
      integer:: ith1,iph1,istep,ip0l,ip1l,jtran,jphot
      integer:: NG_dim,NG_entry,iterm,iJ,ing,nftl,tid

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset
      double precision:: daux,WA,mu_inv,dsm,dsp,larmor
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, &
             dimension(2,nxtran,Geom%nph,Geom%nth,Rz0:Rz1):: Norm
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P

      complex(kind=8), &
         dimension(0:2,0:2,2,nxtran,Geom%nph,Geom%nth,Rz0:Rz1):: BStk

      ! Buffers

      ! Receivers
      double precision, dimension(:), allocatable, target:: Stokes_r
      double precision, dimension(:), allocatable, target:: Prof_r
      ! Senders
      double precision, dimension(:,:,:), allocatable:: Stokes_s
      double precision, dimension(:,:,:), allocatable:: Prof_s
      ! Duals
      integer:: info_b
      integer, dimension(3):: info_c

      ! Pointers
      double precision, dimension(:,:,:), pointer:: data1M,data1O, &
                                                    data1P
      double precision, dimension(:,:), pointer:: data2O,data2P
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP
      double precision, dimension(:,:,:), pointer:: p_MStk
      double precision, dimension(:,:,:), pointer:: p_MProf
      double precision, dimension(:), pointer:: p_exu


      ! Routine name
      urou = 'solver_alt'

      ! Initialize converged flag
      goout = .False.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! If storing redistribution
      RPRAM = PRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AV = .True.
        PRAM = .False.
      end if

      ! Initialize index of Stokes
      op = 1


      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.

      ! If NG acceleration of rho00
      if (Input%NG) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + RnZ*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRD) then

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

          Input%NG = .False.

          if (pid.eq.0) then
            umsg = ' # The buffer for NG acceleration '// &
                   'is too big. Not doing NG.'
            call verbose
          end if

        end if ! Buffer

        ! If finally doing it, allocate
        if (Input%NG) then

          ! Master
          if (pid.eq.0) then

            allocate(NG_scratch(NG_dim, Input%NG_ord+2))

          ! Slaves
          else

            allocate(NG_scratch(NG_dim, 1))

          end if

        end if ! NG and master
      end if ! NG

      ! CPU limits
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)


      !
      ! Allocation of buffers
      !

      ! Master
      if(pid.eq.0)then

        ! To receive Intensity chunks
        iaux = 4*MPID%nxfreq*Rnz
        allocate(Stokes_r(iaux))

        ! To receive profile information
        iaux = MPID%nxtfreq*2*Rnz
        allocate(Prof_r(iaux))

      ! Slave
      else

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(0:3,if0:if1,0:5))
        allocate(data1O(0:3,if0:if1,0:5))

        ! Allocate O pointers
        allocate(data2O(Frec%ntfreq,2))

        ! To send Intensity chunks
        allocate(Stokes_s(0:3,if0:if1,Rz0:Rz1))

        ! To send profile information
        allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1))

      end if ! Master or Slave


      !
      ! Initialization messages
      !

      ! Announce that we are starting
      if (gpid.eq.0) then
        umsg = '   Iteration          MRC(rho^0_0) Atom_index '// &
               'Term_index    2J  Height_index Height(km) |   '// &
               '       MRC(rho^K_Q) Atom_index Term_index '// &
               "   2J   2J' Height_index Height(km)  K  Q"
        call verbose
      end if


      ! Open the file to store MRC
      if (gpid.eq.0.and.Input%keep_MRC) then
        ! If appending
        if (Input%appendMRC) then
          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRC', exist=laux)
          if (.not.laux) then
            open(800, file=trim(Input%folder)//'/MRC', &
                 action='write',iostat=ios,err=1000)
            write(800,'(A)') '!  Iteration          MRC(rho^0_0) '// &
                         'Atom_index Term_index    2J  '// &
                         'Height_index Height(km) |   '// &
                         '       MRC(rho^K_Q) Atom_index '// &
                         'Term_index    2J   2J'// &
                         "' Height_index Height(km)  K  Q"
            close(800)
          end if
        ! If not appending, create
        else
          open(800, file=trim(Input%folder)//'/MRC', &
               action='write',iostat=ios,err=1000)
          write(800,'(A)') '!  Iteration          MRC(rho^0_0) '// &
                       'Atom_index Term_index    2J  '// &
                       'Height_index Height(km) |   '// &
                       '       MRC(rho^K_Q) Atom_index '// &
                       'Term_index    2J   2J'// &
                       "' Height_index Height(km)  K  Q"
          close(800)
        end if
      end if

      ! Control
      call control
      if (laborted) goto 2000

      ! If measuring performance
      if (Input%mpi_perf.and.pid.eq.0) &
        call report_mpi_time(Input%folder,Input%ID, &
                             0,0,.False.)


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=Input%iter_min,Input%iter_max

        ! Flags for physics
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics
        if (iter.le.Input%allownphys_rho) then
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


        !
        ! Master
        !
        if (pid.eq.0) then

!$omp parallel default(none) &
!$omp private(ith1,iph1,id,info_c,Stokes_r,Prof_r,WA,deal,iz,p_exu) &
!$omp private(ia,itran,jtran,daux,dsm) &
!$omp shared(JKQ,JKQS,J00P,Norm,BStk,JKQC,info_b,ith,iph,if0l,iter) &
!$omp shared(if1l,ip0l,ip1l,nfl,nftl,p_MProf,p_MStk,Stokes,stm) &
!$omp shared(nz,na,KSTK,PIRAM,Atom,Atmo,MPID,Geom,Input,Frec) &
!$omp shared(Rz0,Rz1,ierr) &
!$omp shared(MPI_STATUS_IGNORE,MPI_COMM_RT,laborted)

          ! Reset radiation field variables
!$omp workshare
          JKQ = cZero
          JKQS = cZero
          J00P = 0d0
          Norm = 0d0
          BStk = cZero
          JKQC = cZero
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
                  call MPI_recv(Stokes_r(1), MPID%size4(info_b), &
                                MPI_DOUBLE_PRECISION, info_b, &
                                info_b, MPI_COMM_RT, &
                                MPI_STATUS_IGNORE, ierr)
                  if (ierr.eq.0) exit
                end do

                ! Receive profile
                do while (.True.)
                  call MPI_recv(Prof_r(1), MPID%size5(info_b), &
                                MPI_DOUBLE_PRECISION, info_b, &
                                1+info_b, MPI_COMM_RT, &
                                MPI_STATUS_IGNORE, ierr)
                  if (ierr.eq.0) exit
                end do

                ! If measuring performance
                if (Input%mpi_perf) &
                  call report_mpi_time(Input%folder,Input%ID, &
                                       info_b,iter,.True.)

                ! Shorter variables
                if0l = MPID%if0(info_b)
                if1l = MPID%if1(info_b)
                ip0l = Frec%Mpif0(info_b)
                ip1l = Frec%Mpif1(info_b)
                nfl = MPID%nf(info_b)
                nftl = Frec%Mntfreq(info_b)

                ! Pointers
                p_MStk(0:3,if0l:if1l,Rz0:Rz1) => &
                                       Stokes_r(1:MPID%size4(info_b))
                p_MProf(1:nftl,1:2,Rz0:Rz1) => &
                                         Prof_r(1:MPID%size5(info_b))

!$omp end single

                ! Do not process if leaving
                if (laborted) cycle

                ! Get angular weight
                WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                deal = .False.

                ! Each height
!$omp do
                do iz=Rz0,Rz1

                  ! Determine where to store intensity
                  if (KSTK.or.iz.eq.Rz0) &
                    Stokes(0:3,if0l:if1l,iph,ith,iz) = &
                                              p_MStk(0:3,if0l:if1l,iz)
                  ! Point to exu values
                  if (PIRAM.and.ip1l.ge.ip0l) then
                    p_exu => Frec%exu(ip0l:ip1l,iz)
                  else
                    allocate(p_exu(1))
                    deal = .True.
                  end if

                  ! Calculate frequency integral for b-b quantities
                  call FInt_line(Atom,MPID,Geom,Frec%W_freq, &
                                 Frec%Mlif0(info_b), &
                                 Frec%Mlif1(info_b), &
                                 info_b,iph,ith,iz,p_MStk(:,:,iz), &
                                 p_MProf(:,:,iz), &
                                 Norm(:,:,iph,ith,iz), &
                                 BStk(:,:,:,:,iph,ith,iz))

                  !
                  ! Calculate rest of integrals
                  !
                  call FInt_rest(Atom,MPID,Geom,Frec%omega, &
                                 Frec%W_freq,ip0l,ip1l, &
                                 Atmo%T(iz),info_b,iph,ith,WA, &
                                 p_MStk(:,:,iz),J00P(:,:,iz), &
                                 JKQC(:,:,:,iz),p_exu)

                  ! Nullify pointer
                  if (deal) deallocate(p_exu)
                  nullify(p_exu)

                end do ! heights
!$omp end do
              end do ! frequency domains
            enddo ! azimuthal directions
          enddo ! polar directions
!$omp single
          ! Nullify pointers
          nullify(p_MStk,p_MProf)
!$omp end single


          !
          ! Apply weights to JKQ, JKQS, and normalize
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
                  do itran=1,Atom(ia)%ntran

                    ! Apply atomic shift
                    jtran = itran + Atom(ia)%tshift

                    ! Get the weight
                    if (Norm(1,jtran,iph1,ith1,iz).gt.0d0) then

                      daux = WA/Norm(1,jtran,iph1,ith1,iz)

                      ! Integrate angle
                      JKQ(0:2,:,jtran,iz) = JKQ(0:2,:,jtran,iz) + &
                                   BStk(:,:,1,jtran,iph1,ith1,iz)*daux

                    end if

                    ! If there is stimulated emission
                    if (stm) then

                      ! Get the weight
                      if (Norm(2,jtran,iph1,ith1,iz).gt.0d0) then

                        daux = WA/Norm(2,jtran,iph1,ith1,iz)

                        ! Integrate angle
                        JKQS(0:2,:,jtran,iz) = JKQS(0:2,:,jtran,iz) + &
                                   BStk(:,:,2,jtran,iph1,ith1,iz)*daux

                      end if

                    end if ! Stimulated emission
                  end do ! transitions
                end do ! atoms
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

              end do ! b-f transitions
            end do ! atoms
          end do ! heights
!$omp end do
!$omp end parallel

          ! If not axial, need to complete JKQ
          if (.not.axial) then

            ! Continuum
            do K=1,Krad
              do iQ=1,K
                JKQC(-iQ,K,:,:) = Flgsg%sg(iQ)*conjg(JKQC(iQ,K,:,:))
              end do
            end do

            ! For each atom
            do ia=1,nA

              ! For each FS transition
              do itran=1,Atom(ia)%ntran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tshift

                ! K
                do K=1,Atom(ia)%Krad(itran)

                  ! Q
                  do iQ=1,K

                    ! Conjugate
                    JKQ(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                         conjg(JKQ(iQ,K,jtran,:))

                    ! Stimulated, conjugate
                    if (stm) &
                      JKQS(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                            conjg(JKQS(iQ,K,jtran,:))

                  end do ! Q
                end do ! K
              end do ! Transition
            end do ! Atom
          end if ! Not axial

          !
          ! Add the Ad-hoc asymmetries
          !
          if (Input%nasym.gt.0) &
            call addJKQasym(Bfield,Flgsg,JKQ_asym,JKQ,JKQS,JKQC)


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

            ! Determine the first and last height indexes to run over
            iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
            iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

            ! For each azimuthal direction
            do iph=1,Geom%nPh

              ! Error
              if (laborted) goto 3000

              !
              ! If angle-dependent, manage scattering angles
              if (.not.AV.and.PRD) &
                call scattering_manage(Geom,ith,iph)

              !
              ! First height
              !

              ! If going down, get top boundary
              if(diz.eq.1)then

                ! Call top boundary
                call top(MPID,data1M(:,:,5))

              ! If going up, get bottom boundary
              else

                ! Call bottom boundary
                call bottom(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                            Atmo%vy(iz0),Atmo%vz(iz0), &
                            Geom%V_mu(ith),Geom%V_mux(iph), &
                            Geom%V_muy(iph),MPID,data1M(:,:,5))

              endif ! propagation direction

              ! Identify current height
              o = iz0

              ! Index for Stokes
              if (PRD.and.ADD) op = o

              ! Calculate radiative coefficients
              call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                           Geom,o,ith,iph,if0,if1,JKQ_asym, &
                           JKQ(:,:,:,o),JKQC(:,:,:,o),Cont%ndir, &
                           Cont%c(:,:,:,o),Bfield, &
                           Stokes(:,:,:,:,op), data1M(:,:,0:4), &
                           data2O)
              if (laborted) goto 3000


              !
              ! Store in buffer
              !

              ! Stokes
              Stokes_s(:,:,o) = data1M(:,:,5)

              ! Profiles
              Prof_s(:,:,o) = data2O

              ! Identify next height
              p = iz0 + diz

              ! Index for Stokes
              if (PRD.and.ADD) op = p

              ! Calculate radiative coefficients
              call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                           Geom,p,ith,iph,if0,if1,JKQ_asym, &
                           JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                           Cont%c(:,:,:,p),Bfield, &
                           Stokes(:,:,:,:,op),data1O(:,:,0:4), &
                           data2O)
              if (laborted) goto 3000


              !
              ! Intermedium heights
              !

              ! For each height this CPU has assigned
              do iz=iz0,iz1,diz

                ! We treat the boundaries outside
                if(iz.eq.iz0.or.iz.eq.iz1)cycle

                ! Allocate P pointers
                allocate(data1P(0:3,if0:if1,0:5))
                allocate(data2P(Frec%ntfreq,2))

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
                if (PRD.and.ADD) op = p

                ! RT coefficients
                call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                             Geom,p,ith,iph,if0,if1,JKQ_asym, &
                             JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                             Cont%c(:,:,:,p),Bfield, &
                             Stokes(:,:,:,:,p),data1P(:,:,0:4), &
                             data2P)

                ! Point to the data
                p_K0M  => data1M(:,:,0)
                p_K1M  => data1M(:,:,1)
                p_K2M  => data1M(:,:,2)
                p_SM   => data1M(:,:,4)
                p_StkM => data1M(:,:,5)
                p_K0O  => data1O(:,:,0)
                p_K1O  => data1O(:,:,1)
                p_K2O  => data1O(:,:,2)
                p_SO   => data1O(:,:,4)
                p_StkO => data1O(:,:,5)
                p_K0P  => data1P(:,:,0)
                p_SP   => data1P(:,:,4)

                ! Apply short characteristics BESSER
                call RTStep(o,ith,iph,MPID%nf(pid), &
                            dsm,dsp,p_K0M,p_K1M,p_K2M, &
                            p_SM,p_K0O,p_K1O,p_K2O, &
                            p_SO,p_K0P,p_SP,p_StkM, &
                            p_StkO,.True.)

                !
                ! Store in buffer
                !

                ! Stokes
                Stokes_s(:,:,o) = data1O(:,:,5)

                ! Profiles
                Prof_s(:,:,o) = data2O

                ! Shift data (O->M, P->O)
                deallocate(data1M,data2O)
                data1M => data1O
                data1O => data1P
                data2O => data2P
                nullify(data1P,data2P)

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
              p_K0M  => data1M(:,:,0)
              p_K1M  => data1M(:,:,1)
              p_K2M  => data1M(:,:,2)
              p_SM   => data1M(:,:,4)
              p_StkM => data1M(:,:,5)
              p_K0O  => data1O(:,:,0)
              p_K1O  => data1O(:,:,1)
              p_K2O  => data1O(:,:,2)
              p_SO   => data1O(:,:,4)
              p_StkO => data1O(:,:,5)

              ! Apply short characteristics LINEAR
              call RTStep(o,ith,iph,MPID%nf(pid), &
                          dsm,dsp,p_K0M,p_K1M,p_K2M, &
                          p_SM,p_K0O,p_K1O,p_K2O, &
                          p_SO,p_K0P,p_SP,p_StkM, &
                          p_StkO,.False.)

              !
              ! Send to master if error
              !

              ! If had an error
3000          if (laborted) then

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

              ! Stokes
              Stokes_s(:,:,o) = data1O(:,:,5)

              ! Profiles
              Prof_s(:,:,o) = data2O

              ! Send indexes
              info_c = (/ pid, ith, iph /)
              do while (.True.)
                call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                              MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do

              ! Send Stokes
              do while (.True.)
                call MPI_SEND(Stokes_s(0,if0,Rz0), MPID%size4(pid), &
                              MPI_DOUBLE_PRECISION, 0, pid, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

              ! Send profiles
              do while (.True.)
                call MPI_SEND(Prof_s(1,1,Rz0), &
                             MPID%size5(pid), MPI_DOUBLE_PRECISION, &
                             0, 1+pid, MPI_COMM_RT, &
                             ierr)
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
                          MPID%recv, 1+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Receive JKQC
            call MPI_RECV(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                          MPI_DOUBLE_COMPLEX, &
                          MPID%recv, 2+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Receive Stokes if doing A-D PRD
            if (PRD.and.ADD) &
            call MPI_RECV(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                          MPI_DOUBLE_PRECISION, &
                          MPID%recv, 3+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)

            ! Bradcast J00 for b-f transitions
            call MPI_RECV(J00P(1,1,Rz0), MPID%sizei3(0), &
                          MPI_DOUBLE_PRECISION, &
                          MPID%recv, 4+pid, &
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
                           1+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,2), ierr)

            ! Send JKQC
            call MPI_ISEND(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                           MPI_DOUBLE_COMPLEX, &
                           MPID%lsend(istep), &
                           2+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,3), ierr)

            ! Send Stokes if doing A-D PRD
            if (PRD.and.ADD) &
            call MPI_ISEND(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                           MPI_DOUBLE_PRECISION, &
                           MPID%lsend(istep), &
                           3+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,4), ierr)

            ! Send J00 for b-f transitions
            call MPI_ISEND(J00P(1,1,Rz0), MPID%sizei3(0), &
                           MPI_DOUBLE_PRECISION, &
                           MPID%lsend(istep), &
                           4+MPID%lsend(istep), &
                           MPI_COMM_RT, &
                           MPID%requestA(istep,5), ierr)

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

          ! Share JKQC
          call MPI_BCAST(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

          ! Share Stokes if doing A-D PRD
          if (PRD.and.ADD) &
          call MPI_BCAST(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Share J00 for b-f transitions
          call MPI_BCAST(J00P(1,1,Rz0), MPID%sizei3(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Type of bcast


        !
        ! Solve SEE
        !
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,iter)
#endif

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,iphot,jphot,larmor,urou,umsg,tid) &
!$omp shared(nA,nz,Atom,JKQ,JKQS,J00P,Flgsg,Bfield,Rz0,Rz1) &
!$omp shared(laborted,vaborted)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tshift + 1
          jtran = itran + Atom(ia)%ntran - 1
          iphot = Atom(ia)%pshift + 1
          jphot = iphot + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Set magnetic data
            larmor = B2L*Bfield%Bstrength(iz)

            ! Solve the SEE
            call SEE(Atom(ia),JKQ(:,:,itran:jtran,iz), &
                     JKQS(:,:,itran:jtran,iz), &
                     J00P(iphot:jphot,:,iz),larmor,Flgsg,iz,tid)

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel
#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,iter)
#endif

        ! Control
        call control
        if (laborted) goto 2000


        !
        ! NG acceleration
        !

        ! Check if doing NG acceleration
        if(Input%NG.and.iter.gt.Input%NG_delay)then

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
            if (PRD) then

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
                                            Stokes(0,ifreq,iph,ith,iz)

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
                    NG_scratch(o,NG_entry) = dble(JKQC(0,0,ifreq,iz))

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD
            end if ! PRD

            ! Call NG and check if it should be processed
            call NG(NG_dim,Input%NG_ord,NG_scratch,NG_entry,doNG)

          ! Slave
          else

            ! If wrong order
            if (Input%NG_ord.lt.1.or.Input%NG_ord.gt.5) then
              doNG = .False.
            ! Valid order
            else
              ! Check if Master is in NG step
              if (NG_entry.gt.(Input%NG_ord+1)) then
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
                              MPID%recv, 7+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No master

              ! For each send
              do istep=1,MPID%nsend

                ! Send NG iteration
                call MPI_ISEND(NG_scratch(1,ing), NG_dim, &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               7+MPID%lsend(istep), &
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
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Advance index
                    o = o + 1

                    ! Accelerate rho
                    Atom(ia)%crho(p,iz) = &
                                        dcmplx(NG_scratch(o,ing), 0d0)

                  enddo ! Heights
                enddo ! Levels
              enddo ! Terms
            enddo ! Atoms

            ! If PRD
            if (PRD) then

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
                        Stokes(0,ifreq,iph,ith,iz) = NG_scratch(o,ing)

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
                    JKQC(0,0,ifreq,iz) = &
                                        dcmplx(NG_scratch(o,ing), 0d0)
                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD
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
          call MRC_sb(Atom,Rho_old,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5
          MRC%values(1,2) = Atmo%z(MRC%indexes(2,2))*1d-5

          ! Check exit criteria
          if (MRC%values(2,1).le.Input%mrc_i.and. &
              MRC%values(2,2).le.Input%mrc_p) &
            goout = .True.

        end if

        ! Only global master
        if (gpid.eq.0) then

          ! Write in stdout
          write(umsg, &
          '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
          '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,1x,i4,3x,i11,'// &
          '2x,f9.3,1x,i2,1x,i2)') &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(5,1),MRC%indexes(2,1), &
          MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
          MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
          MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
          MRC%indexes(6,2)
          call verbose

          ! If file
          if (Input%keep_MRC) then

            ! Write in MRC file
            open(800, file=trim(Input%folder)//'/MRC', &
                 iostat=ios,err=1000,position='append')
            write(800, &
            '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
            '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,2x,i4,2x,i11,'// &
            '2x,f9.3,1x,i2,1x,i2)', err=1100) &
            iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
            MRC%indexes(5,1),MRC%indexes(2,1), &
            MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
            MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
            MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
            MRC%indexes(6,2)
            close(800)

          end if
        end if

        ! We can swith now to AD if we had AV input
        if (tbAD) then
          AV = .False.
          tbAD = .False.
          PRAM = RPRAM
          ! Force at least two iterations
          goout = .False.
        end if

        !
        ! Save partial solution
        !

        ! Check if we want to store partial results
        if(Input%store.and.mod(iter,Input%store_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesol(Input,iterS,Frec%omega,Geom,Flgsg, &
                        Bfield,Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)
        ! Or have a check
        else

          call control
          if (laborted) goto 2000

        end if

        !
        ! Share if we are finished or not

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive goout
            call MPI_RECV(goout, 1, MPI_LOGICAL, &
                          MPID%recv, 5+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, ierr)

          end if ! No Master

          ! For each send
          do istep=1,MPID%nsend

            ! Send goout
            call MPI_ISEND(goout, 1, MPI_LOGICAL, &
                           MPID%lsend(istep), &
                           5+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,6), &
                           ierr)

          end do ! sends

        ! Normal bcast
        else

          ! Sjare goout
          call MPI_BCAST(goout, 1, MPI_LOGICAL, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Type of bcast

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

          end do ! Processors

        end if ! Domain decomposition

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
        if (ierr.ne.0) goto 1000

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = real(MRC%values(2,2))

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

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
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO)
      end if

      return

1000  umsg = 'Error opening MRC file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing MRC file'
      close(800)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine solver_alt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the NLTE problem of the second kind with lambda
      !! iteration with one CPU (serial).\n
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
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      subroutine solver_serial(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                               Red,Bfield,Geom,MPID,Input,Flgsg, &
                               JKQ_asym,Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Rhoc_class), dimension(:):: Rho_old
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      double precision, dimension(0:3,nfreq,Geom%nPh, &
                                  Geom%nTh,giz0:giz1), target:: Stokes
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

      ! Local

      type(MRC_class):: MRC

      character(LEN=20):: iterS

      logical:: goout, AD, ADD, doNG, laux, RPRAM, lp_exu

      integer:: ios,iz0,iz1,diz,m,o,p,op,itran,jtran,iphot,jphot
      integer:: iter,ith,iph,ia,iz,ifreq,if0,if1,K,iQ
      integer:: NG_dim,NG_entry,iterm,iJ,tid

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset

      double precision:: mu_inv,dsm,dsp,larmor
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P
      double precision, &
        dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
                                                     target:: Stokes_n

      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ_n
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQS_n
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC_n

      ! Pointers
      double precision, dimension(:,:,:), pointer:: data1M,data1O, &
                                                    data1P
      double precision, dimension(:,:), pointer:: data2O,data2P
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP
      double precision, dimension(:), pointer:: p_exu


      ! Routine name
      urou = 'solver_serial'

      ! Initialize converged flag
      goout = .False.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! If storing redistribution
      RPRAM = PRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AV = .True.
        PRAM = .False.
      end if


      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.

      ! If NG acceleration
      if (Input%NG) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + Rnz*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRD) then

          ! If we need Stokes
          if (ADD) then

            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            NG_dim = NG_dim + nfreq*RnZ

          end if !ADD
        end if ! PRD

        ! If it requires too much buffer
        if (dble(NG_dim)*8d-6.gt.maxbuffer_NG) then

          Input%NG = .False.

          umsg = ' # The buffer for NG acceleration '// &
                 'is too big. Not doing NG.'
          call verbose

        end if

        ! If finally doing it, allocate
        if (Input%NG) then

          allocate(NG_scratch(NG_dim, Input%NG_ord+2))

        end if ! NG
      end if ! NG

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = 1
      if1 = nfreq

      !
      ! Allocation
      !

      ! Common (Master and slave)
      ! Allocate O pointers
      allocate(data2O(Frec%ntfreq,2))

      ! Allocate M and O pointers for RT coeff
      allocate(data1M(0:3,nfreq,0:5))
      allocate(data1O(0:3,nfreq,0:5))

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

      ! Global Master
      if (gpid.eq.0) then

        ! Announce we are starting
        umsg = '   Iteration          MRC(rho^0_0) Atom_index '// &
               'Term_index    2J  Height_index Height(km) |   '// &
               '       MRC(rho^K_Q) Atom_index Term_index '// &
               "   2J   2J' Height_index Height(km)  K  Q"
        call verbose
      end if

      ! Open the file to store MRC
      if (gpid.eq.0.and.Input%keep_MRC) then
        ! If appending
        if (Input%appendMRC) then
          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRC', exist=laux)
          if (.not.laux) then
            open(800, file=trim(Input%folder)//'/MRC', &
                 action='write',iostat=ios,err=1000)
            write(800,'(A)') '! MRC file'
            write(800,'(A)') '!  Iteration          MRC(rho^0_0) '// &
                         'Atom_index Term_index    2J  '// &
                         'Height_index Height(km) |   '// &
                         '       MRC(rho^K_Q) Atom_index '// &
                         'Term_index    2J   2J'// &
                         "' Height_index Height(km)  K  Q"
            close(800)
          end if
        ! If not appending, create
        else
          open(800, file=trim(Input%folder)//'/MRC', &
               action='write',iostat=ios,err=1000)
          write(800,'(A)') '!  Iteration          MRC(rho^0_0) '// &
                       'Atom_index Term_index    2J  '// &
                       'Height_index Height(km) |   '// &
                       '       MRC(rho^K_Q) Atom_index '// &
                       'Term_index    2J   2J'// &
                       "' Height_index Height(km)  K  Q"
          close(800)
        end if
      end if

      ! Get the current rhoKQ into the 'old' structure
      do ia=1,nA
        Rho_old(ia)%crho = Atom(ia)%crho
      end do


      !
      ! Iterate
      !

      ! For each iteration between the limits specified
      do iter=Input%iter_min,Input%iter_max

        ! Flags for physics
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics
        if (iter.le.Input%allownphys_rho) then
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

        ! Reset radiation field variables
!$omp parallel workshare default(none) &
!$omp shared(Stokes_n,JKQ_n,JKQS_n,JKQC_n,J00P)
        Stokes_n = 0d0
        JKQ_n = cZero
        JKQS_n = cZero
        JKQC_n = cZero
        J00P = 0d0
!$omp end parallel workshare


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
            if (.not.AV.and.PRD) &
              call scattering_manage(Geom,ith,iph)

            !
            ! First height
            !

            ! If going down, get top boundary
            if(diz.eq.1)then

              ! Call top boundary
              call top(MPID,data1M(:,:,5))

            ! If going up, get bottom boundary
            else

              ! Call bottom boundary
              call bottom(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                          Atmo%vy(iz0),Atmo%vz(iz0), &
                          Geom%V_mu(ith),Geom%V_mux(iph), &
                          Geom%V_muy(iph),MPID,data1M(:,:,5))

            endif

            ! Identify current height
            o = iz0

            ! Index for Stokes
            if (PRD.and.ADD) op = o

            ! Calculate radiative coefficients
            call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                         Geom,o,ith,iph,if0,if1,JKQ_asym, &
                         JKQ(:,:,:,o),JKQC(:,:,:,o),Cont%ndir, &
                         Cont%c(:,:,:,o),Bfield, &
                         Stokes(:,:,:,:,op), &
                         data1M(:,:,0:4),data2O)
            if (laborted) goto 2000

            if (KSTK) Stokes_n(:,:,iph,ith,o) = data1M(:,:,5)

            ! Point to exu values
            if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

            !
            ! Calculate integrals
            !
            call Jcalc(Atom,Geom,Frec%omega,Frec%W_freq, &
                       Frec%lif0,Frec%lif1, &
                       Frec%pif0,Frec%pif1, &
                       Atmo%T(o),Atmo%ne(o),iph,ith,o, &
                       data1M(:,:,5), &
                       data2O,JKQ_n(:,:,:,o), &
                       JKQS_n(:,:,:,o),J00P(:,:,o), &
                       JKQC_n(:,:,:,o),p_exu)

            ! Identify next height
            p = iz0 + diz

            ! Index for Stokes
            if (PRD.and.ADD) op = p

            ! Calculate radiative coefficients
            call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                         Geom,p,ith,iph,if0,if1,JKQ_asym, &
                         JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                         Cont%c(:,:,:,p),Bfield, &
                         Stokes(:,:,:,:,op), &
                         data1O(:,:,0:4),data2O)
            if (laborted) goto 2000


            !
            ! Intermedium heights
            !

            ! For each height this CPU has assigned
            do iz=iz0,iz1,diz

              ! We treat the boundaries outside
              if(iz.eq.iz0.or.iz.eq.iz1)cycle

              ! Allocate P pointers
              allocate(data1P(0:3,nfreq,0:5))
              allocate(data2P(Frec%ntfreq,2))

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
              if (PRD.and.ADD) op = p

              ! RT coefficients
              call RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                           Geom,p,ith,iph,if0,if1,JKQ_asym, &
                           JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                           Cont%c(:,:,:,p),Bfield, &
                           Stokes(:,:,:,:,op), &
                           data1P(:,:,0:4),data2P)

              ! Point to the data
              p_K0M  => data1M(:,:,0)
              p_K1M  => data1M(:,:,1)
              p_K2M  => data1M(:,:,2)
              p_SM   => data1M(:,:,4)
              p_StkM => data1M(:,:,5)
              p_K0O  => data1O(:,:,0)
              p_K1O  => data1O(:,:,1)
              p_K2O  => data1O(:,:,2)
              p_SO   => data1O(:,:,4)
              p_StkO => data1O(:,:,5)
              p_K0P  => data1P(:,:,0)
              p_SP   => data1P(:,:,4)

              ! Apply short characteristics BESSER
              call RTStep(o,ith,iph,nfreq, &
                          dsm,dsp,p_K0M,p_K1M,p_K2M, &
                          p_SM,p_K0O,p_K1O,p_K2O, &
                          p_SO,p_K0P,p_SP,p_StkM, &
                          p_StkO,.True.)

              if (KSTK) Stokes_n(:,:,iph,ith,o) = data1O(:,:,5)

              ! Point to exu values
              if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

              ! Calculate contribution to JKQ
              call Jcalc(Atom,Geom,Frec%omega,Frec%W_freq, &
                         Frec%lif0,Frec%lif1, &
                         Frec%pif0,Frec%pif1, &
                         Atmo%T(o),Atmo%ne(o),iph,ith,o, &
                         data1O(:,:,5),data2O, &
                         JKQ_n(:,:,:,o),JKQS_n(:,:,:,o), &
                         J00P(:,:,o),JKQC_n(:,:,:,o),p_exu)

              ! Shift data (O->M, P->O)
              deallocate(data1M,data2O)
              data1M => data1O
              data1O => data1P
              data2O => data2P
              nullify(data1P,data2P)

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
            p_K0M  => data1M(:,:,0)
            p_K1M  => data1M(:,:,1)
            p_K2M  => data1M(:,:,2)
            p_SM   => data1M(:,:,4)
            p_StkM => data1M(:,:,5)
            p_K0O  => data1O(:,:,0)
            p_K1O  => data1O(:,:,1)
            p_K2O  => data1O(:,:,2)
            p_SO   => data1O(:,:,4)
            p_StkO => data1O(:,:,5)

            ! Apply short characteristics LINEAR
            call RTStep(o,ith,iph,nfreq, &
                        dsm,dsp,p_K0M,p_K1M,p_K2M, &
                        p_SM,p_K0O,p_K1O,p_K2O, &
                        p_SO,p_K0P,p_SP,p_StkM, &
                        p_StkO,.False.)
            if (laborted) goto 2000

            if (KSTK) Stokes_n(:,:,iph,ith,o) = data1O(:,:,5)

            ! Point to exu values
            if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

            !
            ! Calculate integrals
            !
            call Jcalc(Atom,Geom,Frec%omega,Frec%W_freq, &
                       Frec%lif0,Frec%lif1, &
                       Frec%pif0,Frec%pif1, &
                       Atmo%T(o),Atmo%ne(o),iph,ith,o, &
                       data1O(:,:,5),data2O,JKQ_n(:,:,:,o), &
                       JKQS_n(:,:,:,o),J00P(:,:,o), &
                       JKQC_n(:,:,:,o),p_exu)

          enddo ! azimuthal directions
        enddo ! polar directions


        ! If not axial, need to complete JKQ
        if (.not.axial) then

          ! Continuum
          do K=1,Krad
            do iQ=1,K
              JKQC_n(-iQ,K,:,:) = Flgsg%sg(iQ)*conjg(JKQC_n(iQ,K,:,:))
            end do
          end do

          ! For each atom
          do ia=1,nA

            ! For each FS transition
            do itran=1,Atom(ia)%ntran

              ! Apply atomic shift
              jtran = itran + Atom(ia)%tshift

              ! K
              do K=1,Atom(ia)%Krad(itran)

                ! Q
                do iQ=1,K

                  ! Conjugate
                  JKQ_n(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                         conjg(JKQ_n(iQ,K,jtran,:))

                  ! Stimulated, conjugate
                  if (stm) &
                    JKQS_n(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                           conjg(JKQS_n(iQ,K,jtran,:))

                end do ! Q
              end do ! K
            end do ! Transition
          end do ! Atom
        end if ! Not axial


        !
        ! Add the Ad-hoc asymmetries
        !
        if (Input%nasym.gt.0) &
          call addJKQasym(Bfield,Flgsg,JKQ_asym,JKQ_n,JKQS_n,JKQC)


        !
        ! Solve SEE
        !
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ_n,JKQS_n, &
                                  Input%folder,iter)
#endif

!$omp parallel default(none) &
!$omp private(ia,iz,itran,jtran,iphot,jphot,larmor,umsg,urou,tid) &
!$omp shared(nA,nz,Atom,JKQ_n,JKQS_n,J00P,Flgsg,Bfield,Rz0,Rz1) &
!$omp shared(laborted,vaborted)
#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
#else
        tid = -1
#endif
        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tshift + 1
          jtran = itran + Atom(ia)%ntran - 1
          iphot = Atom(ia)%pshift + 1
          jphot = iphot + Atom(ia)%nphot - 1

          ! For each height
!$omp do
          do iz=Rz0,Rz1

            ! Set magnetic data
            larmor = B2L*Bfield%Bstrength(iz)

            ! Solve the SEE
            call SEE(Atom(ia),JKQ_n(:,:,itran:jtran,iz), &
                     JKQS_n(:,:,itran:jtran,iz), &
                     J00P(iphot:jphot,:,iz),larmor,Flgsg,iz,tid)

          end do ! heights
!$omp end do nowait
        end do ! atoms
!$omp end parallel
#ifdef DEBUGRHOKQ
        call dump_rho(Atom,Input%folder,iter)
#endif

        ! Error
        if (laborted) goto 2000


        !
        ! NG acceleration
        !

        ! Check if doing NG acceleration
        if(Input%NG.and.iter.gt.Input%NG_delay)then

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
          if (PRD) then

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
                                          Stokes_n(0,ifreq,iph,ith,iz)

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
                  NG_scratch(o,NG_entry) = dble(JKQC_n(0,0,ifreq,iz))

                end do ! Frequencies
              end do ! Heights

            end if ! AV or AD
          end if ! PRD

          ! Call NG and check if it should be processed
          call NG(NG_dim,Input%NG_ord,NG_scratch,NG_entry,doNG)

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
                  p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Advance index
                    o = o + 1

                    ! Accelerate rho
                    Atom(ia)%crho(p,iz) = &
                                   dcmplx(NG_scratch(o,NG_entry), 0d0)

                  enddo ! Heights
                enddo ! Levels
              enddo ! Terms
            enddo ! Atoms

            ! PRD
            if (PRD) then

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
                        Stokes_n(0,ifreq,iph,ith,iz) = &
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
                    JKQC_n(0,0,ifreq,iz) = &
                                   dcmplx(NG_scratch(o,NG_entry), 0d0)

                  end do ! Frequencies
                end do ! Heights

              end if ! AV or AD
            end if ! PRD

            ! Send message if global Master
            if (gpid.eq.0) then
              umsg = 'NG acceleration'
              call verbose
            end if

            ! Reset entry
            NG_entry = 0

          endif ! Apply NG
        endif ! Doing NG acceleration

        !
        !  Calculate MRC
        !

        call MRC_sb(Atom,Rho_old,MRC)

        ! Convert cm into km
        MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5
        MRC%values(1,2) = Atmo%z(MRC%indexes(2,2))*1d-5

        ! Global Master
        if (gpid.eq.0) then

          ! Write in stdout
          write(umsg, &
          '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
          '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,1x,i4,3x,i11,'// &
          '2x,f9.3,1x,i2,1x,i2)') &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(5,1),MRC%indexes(2,1), &
          MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
          MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
          MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
          MRC%indexes(6,2)
          call verbose

          ! File
          if (Input%keep_MRC) then

            ! Write in MRC file
            open(800, file=trim(Input%folder)//'/MRC', &
                 iostat=ios,err=1000,position='append')
            write(800, &
            '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
            '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,2x,i4,2x,i11,'// &
            '2x,f9.3,1x,i2,1x,i2)', err=1100) &
            iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
            MRC%indexes(5,1),MRC%indexes(2,1), &
            MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
            MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
            MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
            MRC%indexes(6,2)
            close(800)

          end if
        end if


        ! Check exit criteria
        if (MRC%values(2,1).le.Input%mrc_i.and. &
            MRC%values(2,2).le.Input%mrc_p) &
          goout = .True.

        ! We can swith now to AD if we had AV input
        if (tbAD) then
          AV = .False.
          tbAD = .False.
          PRAM = RPRAM
          ! Force at least two iterations
          goout = .False.
        end if

        ! Save data of this steps if proceeds
        if(Input%store.and.mod(iter,Input%store_step).eq.0)then

          write(iterS,'(i0.4)') iter
          call writesol(Input,iterS,Frec%omega,Geom,Flgsg, &
                        Bfield,Atom,Atmo%z,Stokes_n,JKQ_n, &
                        JKQS_n,JKQC)

        end if

        ! Shift the new values into the proper variables
!$omp parallel default(none) &
!$omp shared(KSTK,Stokes,Stokes_n,JKQ,JKQS,JKQ_n,JKQS_n) &
!$omp shared(JKQC,JKQC_n)
        if (KSTK) then
!$omp workshare
          Stokes = Stokes_n
!$omp end workshare
        end if
!$omp workshare
        JKQ = JKQ_n
        JKQS = JKQS_n
        JKQC = JKQC_n
!$omp end workshare
!$omp end parallel

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
        if (ierr.ne.0) goto 1000

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = real(MRC%values(2,2))

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      ! Clean pointers
2000  deallocate(data1M,data1O,data2O)
      nullify(data1M,data1O,data2O)
      nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
      nullify(p_StkM,p_StkO)
      if (lp_exu) then
        nullify(p_exu)
      else
        deallocate(p_exu)
        nullify(p_exu)
      end if

      call control

      return

1000  umsg = 'Error opening MRC file'
      call aborted
      return
1010  umsg = 'Error seeking MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call aborted
      return
1100  umsg = 'Error writing MRC file'
      close(800)
      call aborted
      return
1300  umsg = 'Error writing MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call aborted
      return

      end subroutine solver_serial

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the last formal solution.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      subroutine emergent(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                          Geom,MPID,Input,Flgsg,JKQ_asym,Stokes, &
                          JKQ,JKQC,SolF)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      type(Solution_F_class):: SolF
      double precision, &
        dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                           giz0:giz1), target:: Stokes
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

      ! MPI
      if(MPID%mpi)then
        call emergence(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                       Geom,MPID,Input,Flgsg,JKQ_asym,Stokes, &
                       JKQ,JKQC,SolF)
      ! Serial
      else
        call emergence_serial(Atom,LTElines,Atmo,Cont,Frec,Red, &
                              Bfield,Geom,MPID,Input,Flgsg, &
                              JKQ_asym,Stokes,JKQ,JKQC,SolF)
      end if

      end subroutine emergent

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
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      subroutine emergence(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                           Geom,MPID,Input,Flgsg,JKQ_asym,Stokes, &
                           JKQ,JKQC,SolF)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      type(Solution_F_class):: SolF
      double precision, &
        dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                           giz0:giz1), target:: Stokes
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

      ! Local

      logical:: AD, ADD

      integer:: iaux,ncount,tau1size
      integer:: iz0,iz1,diz,m,o,p,op,sshift
      integer:: ith,iph,iz,ifreq,if0,if1,id,icount

      double precision:: mu_inv,dsm,dsp,dzm,dzp
      double precision, dimension(:,:), allocatable, target:: tau1
      double precision, dimension(:,:), allocatable, target:: tau
      double precision, dimension(:,:), allocatable:: Contr
      double precision, dimension(:,:,:), allocatable:: ContrG
      double precision, dimension(0:3,nfreq):: Stokes_out

      ! Buffers
      ! Receivers
      integer:: rpid
      double precision, dimension(:), allocatable:: Stokes_r
      double precision, dimension(:), allocatable:: Contr_r
      ! Senders
      double precision, dimension(:,:,:), allocatable:: Contr_s
      double precision, dimension(:,:), allocatable:: Stokes_s
      double precision, dimension(:,:), allocatable:: tau1_s
      ! Dual
      integer, dimension(2):: info_b

      ! Pointers
      double precision, dimension(:,:,:), pointer:: data1M,data1O, &
                                                    data1P
      double precision, dimension(:), pointer:: etaIM,etaIO
      double precision, dimension(:), pointer:: tauM
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP


      !
      ! Initializations
      !

      ! Routine name
      urou = 'emergence'

      ! Reset progress counter
      icount = 0

      ! Determine number of directions to do
      ncount = Geom%nThLOS*Geom%nPhLOS

      ! Initialize angle depended flag
      AD = .not.AV
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
        iaux = 4*MPID%nxfreq
        allocate(Stokes_r(iaux))

        ! If calculating height of tau=1, allocate
        if (Input%out_tau1) allocate(tau1(2,nfreq))

        ! If calculating contribution function, allocate
        if (Input%out_contr) then
          iaux = 4*MPID%nxfreq*Rnz
          allocate(Contr_r(iaux))
          allocate(ContrG(0:3,nfreq,Rz0:Rz1))
        end if

        ! If inverting, need to return the output
        if (run_mode.eq.-1) then
          allocate(SolF%e_Stk(0:3,nfreq,Geom%nPhLOS,Geom%nThLOS))
          if (Input%out_tau1) &
          allocate(SolF%e_tau1(Input%lim_tau%nn, &
                               Geom%nPhLOS,Geom%nThLOS))
          if (Input%out_contr) &
          allocate(SolF%e_Ctr(0:3,Input%lim_ctr%nn,Rz0:Rz1, &
                              Geom%nPhLOS,Geom%nThLOS))
        end if ! Inversion

      ! Slave
      else

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(0:3,if0:if1,0:5))
        allocate(data1O(0:3,if0:if1,0:5))

        ! To send Stokes chunks
        allocate(Stokes_s(0:3,if0:if1))

        ! If calculating tau 1 or contribution function, allocate
        if (input%out_tau1.or.input%out_contr) then

          allocate(tau(if0:if1,Rz0:Rz1))
          allocate(tau1(2,if0:if1))
          allocate(tau1_s(2,if0:if1))
          tau1size = MPID%nf(pid)*2
          allocate(etaIM(if0:if1))

          ! If calculating contribution function, allocate
          if (input%out_contr) then

            allocate(Contr(0:3,1:MPID%nf(pid)))
            allocate(Contr_s(0:3,1:MPID%nf(pid),Rz0:Rz1))

          end if ! contribution output
        else
          nullify(etaIM,etaIO,tauM)
        end if ! tau1 output

      end if ! Master or slave

      ! Control
      call control
      if (laborted) goto 2000


      !
      ! Formal solutions
      !


      ! For each LOS polar direction
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

        end if ! Slave

        ! For each LOS azimuthal direction
        do iph=1,Geom%nPhLOS

          !
          ! Master
          !
          if (pid.eq.0) then

            ! Advance the counter
            icount = icount + 1

            ! Write progress into stdout if global Master
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
                              MPI_ANY_SOURCE, 2, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)
                tau1size = MPID%nf(rpid)*2

                call MPI_recv(tau1(1,MPID%if0(rpid)), &
                              tau1size, &
                              MPI_DOUBLE_PRECISION, rpid, &
                              3+rpid, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)

              end do

              ! If inversion
              if (run_mode.eq.-1) then

                ! Keep tau1
                call settau(SolF%e_tau1(:,iph,ith),tau1(2,:), &
                            Input%lim_tau)

              ! If synthesis
              else

                ! Store the height where tau=1
                call writetau(Input%folder,iph,ith,Frec%omega,Geom, &
                              tau1(2,:),Input%lim_tau)
                call control
                if (laborted) goto 2000

              end if
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
                call MPI_recv(Contr_r(1), MPID%size3(info_b(1)), &
                              MPI_DOUBLE_PRECISION, info_b(1), &
                              1+info_b(1), MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)

                ! Reset shift in index
                sshift = 0

                ! For each height
                do iz=Rz0,Rz1

                  ! Rearrange the contribution function
                  do ifreq=0,MPID%nf(info_b(1))-1

                    ContrG(:,MPID%if0(info_b(1))+ifreq,iz) = &
                                         Contr_r(sshift+4*ifreq+1: &
                                                 sshift+4*(ifreq+1))

                  end do ! frequencies

                  ! Update the shift in the buffer
                  sshift = sshift + 4*MPID%nf(info_b(1))

                end do ! heights

              end if ! if contribution function

              ! If it is not a boundary, we are not receiving
              ! anything else
              if (info_b(2).lt.0) cycle

              ! Receive Stokes
              call MPI_recv(Stokes_r(1), MPID%size10(info_b(1)), &
                            MPI_DOUBLE_PRECISION, info_b(1), &
                            info_b(1), MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)

              ! Rearrange the Stokes
              do ifreq=0,MPID%nf(info_b(1))-1

                Stokes_out(:,MPID%if0(info_b(1))+ifreq) = &
                                     Stokes_r(4*ifreq+1:4*(ifreq+1))

              end do ! frequencies
            end do ! Frequenct blocks


          !
          ! Slave
          !
          else

            !
            ! Get geometry if PRD AD
            !
            if (PRD.and..not.AV) &
              call get_scattering_los(Geom,ith,iph)

            !
            ! Set geometrical tensors
            !
            call setTKQLOS(Geom,Flgsg,Bfield,ith,iph)

            ! Wait till last communication was received
            call MPI_WAIT(MPID%request3,MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%request4,MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%request5,MPI_STATUS_IGNORE,ierr)

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
              call RTAbs(Frec,Atom,LTElines,Atmo,Flgsg,Geom, &
                         o,ith,iph,if0,if1,Cont%ndir, &
                         Cont%c(:,:,:,o),Bfield,etaIM)

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
                call RTAbs(Frec,Atom,LTElines,Atmo,Flgsg,Geom, &
                           o,ith,iph,if0,if1,Cont%ndir, &
                           Cont%c(:,:,:,o),Bfield,etaIO)

                ! Accumulate tau
                call RTtau(MPID%nf(pid),dsm,Atmo%z(m),Atmo%z(o), &
                           etaIM(if0:if1),etaIO(if0:if1), &
                           tauM,tau(if0:if1,o),tau1(:,if0:if1))


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
              call RTAbs(Frec,Atom,LTElines,Atmo,Flgsg,Geom, &
                         o,ith,iph,if0,if1,Cont%ndir, &
                         Cont%c(:,:,:,o),Bfield,etaIO)

              ! Accumulate tau
              call RTtau(MPID%nf(pid),dsm,Atmo%z(m),Atmo%z(o), &
                         etaIM(if0:if1),etaIO(if0:if1),tauM, &
                         tau(if0:if1,o), tau1(:,if0:if1))

              ! Clean pointers
              deallocate(etaIO)
              nullify(etaIO)

              ! If outputting the height of tau=1
              if (input%out_tau1) then

                !
                ! Send to master the tau=1
                !

                ! Wait for last send to finish
                call MPI_WAIT(MPID%request7,MPI_STATUS_IGNORE, &
                              ierr)
                call MPI_WAIT(MPID%request8,MPI_STATUS_IGNORE, &
                              ierr)

                ! Send indexes
                call MPI_ISEND(pid,1,MPI_INTEGER, &
                               0,2,MPI_COMM_RT, &
                               MPID%request7,ierr)


                ! Send Intensity
                tau1_s = tau1
                call MPI_ISEND(tau1_s(1,if0), tau1size, &
                               MPI_DOUBLE_PRECISION, &
                               0,3+pid,MPI_COMM_RT, &
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
              call top(MPID,data1M(:,:,5))

            ! If going up, get bottom boundary
            else

              ! Call bottom boundary
              call bottom(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                          Atmo%vy(iz0),Atmo%vz(iz0), &
                          Geom%L_mu(ith),cos(Geom%L_phi(iph)), &
                          sin(Geom%L_phi(iph)),MPID,data1M(:,:,5))

            endif ! propagation direction

            ! Identify current height
            o = iz0

            ! Index for Stokes
            if (PRD.and.ADD) op = o

            ! Calculate radiative coefficients
            call RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom, &
                          o,ith,iph,if0,if1,JKQ_asym, &
                          JKQ(:,:,:,o),JKQC(:,:,:,o),Cont%ndir, &
                          Cont%c(:,:,:,o),Bfield, &
                          Stokes(:,:,:,:,op),data1M(:,:,0:4))

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Store in buffer. The first does not contribute
              contr_s(:,:,o) = 0d0

            end if ! computing contribution function

            ! Identify next height
            p = iz0 + diz

            ! Index for Stokes
            if (PRD.and.ADD) op = p

            ! Calculate radiative coefficients
            call RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom, &
                          p,ith,iph,if0,if1,JKQ_asym, &
                          JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                          Cont%c(:,:,:,p),Bfield, &
                          Stokes(:,:,:,:,op),data1O(:,:,0:4))


            !
            ! Intermedium heights
            !

            ! For each height this CPU has assigned
            do iz=iz0,iz1,diz

              ! We treat the boundaries outside
              if(iz.eq.iz0.or.iz.eq.iz1)cycle

              ! Allocate P pointers
              allocate(data1P(0:3,if0:if1,0:5))

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
              if (PRD.and.ADD) op = p

              ! RT coefficients
              call RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom, &
                            p,ith,iph,if0,if1,JKQ_asym, &
                            JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                            Cont%c(:,:,:,p),Bfield, &
                            Stokes(:,:,:,:,op),data1P(:,:,0:4))

              ! Point to the data
              p_K0M  => data1M(:,:,0)
              p_K1M  => data1M(:,:,1)
              p_K2M  => data1M(:,:,2)
              p_SM   => data1M(:,:,4)
              p_StkM => data1M(:,:,5)
              p_K0O  => data1O(:,:,0)
              p_K1O  => data1O(:,:,1)
              p_K2O  => data1O(:,:,2)
              p_SO   => data1O(:,:,4)
              p_StkO => data1O(:,:,5)
              p_K0P  => data1P(:,:,0)
              p_SP   => data1P(:,:,4)

              ! Apply short characteristics BESSER
              call RTStep(o,ith,iph,MPID%nf(pid), &
                          dsm,dsp,p_K0M,p_K1M,p_K2M, &
                          p_SM,p_K0O,p_K1O,p_K2O, &
                          p_SO,p_K0P,p_SP,p_StkM, &
                          p_StkO,.True.)

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
                call RTContr(MPID%nf(pid),dsm,dsp,dzm,dzp,p_K0M, &
                             p_K0O,p_K1O,p_K2O,p_SO,p_K0P, &
                             p_StkO,tau(if0:if1,o),Contr,.True.)

                ! Store in buffer
                contr_s(:,:,o) = contr

              end if ! calculating contribution function

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
            p_K0M  => data1M(:,:,0)
            p_K1M  => data1M(:,:,1)
            p_K2M  => data1M(:,:,2)
            p_SM   => data1M(:,:,4)
            p_StkM => data1M(:,:,5)
            p_K0O  => data1O(:,:,0)
            p_K1O  => data1O(:,:,1)
            p_K2O  => data1O(:,:,2)
            p_SO   => data1O(:,:,4)
            p_StkO => data1O(:,:,5)

            ! Apply short characteristics LINEAR
            call RTStep(o,ith,iph,MPID%nf(pid), &
                        dsm,dsp,p_K0M,p_K1M,p_K2M, &
                        p_SM,p_K0O,p_K1O,p_K2O, &
                        p_SO,p_K0P,p_SP,p_StkM, &
                        p_StkO,.False.)

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Calculate vertical distance between points
              if (ztau) then
                dzm = dsm/mu_inv
              else
                dzm = Atmo%z(o) - Atmo%z(m)
              end if

              ! Calculate contribution function
              call RTContr(MPID%nf(pid),dsm,dsp,dzm,dzp,p_K0M, &
                           p_K0O,p_K1O,p_K2O,p_SO,p_K0P, &
                           p_StkO,tau(if0:if1,o),Contr,.False.)

            end if ! calculating contribution function


            !
            ! Send to master
            !

            ! Send indexes
            info_b = (/ pid , 1 /)
            call MPI_ISEND(info_b(1),2,MPI_INTEGER, &
                           0,0, MPI_COMM_RT,MPID%request3, &
                           ierr)

            ! If calculating contribution function
            if (Input%out_contr) then

              ! Send contribution function
              contr_s(:,:,o) = contr
              call MPI_ISEND(contr_s(0,1,Rz0), &
                             MPID%size3(pid), MPI_DOUBLE_PRECISION, &
                             0, 1+pid, MPI_COMM_RT, &
                             MPID%request5, ierr)
            end if

            ! Send intensity
            Stokes_s = data1O(:,:,5)
            call MPI_ISEND(Stokes_s(0,MPID%if0(pid)), &
                           MPID%size10(pid), MPI_DOUBLE_PRECISION, &
                           0, pid, MPI_COMM_RT, &
                           MPID%request4, ierr)

          end if ! Master or slave


          ! Master
          if (pid.eq.0) then

            ! If inverting
            if (run_mode.eq.-1) then

              ! Keep Stokes
              call setstk(SolF%e_Stk(:,:,iph,ith),Stokes_out, &
                          Input%lim_stk,.False.)

              ! Keep contribution function
              if (Input%out_contr) &
                call setctr(SolF%e_Ctr(:,:,:,iph,ith),ContrG, &
                            Input%lim_ctr)

            ! Synthesis
            else

              ! Write stokes
              call writestk(Input%folder,iph,ith,Frec%omega,Geom, &
                            Stokes_out,Input%lim_stk)
              if (laborted) goto 2000

              ! Write contribution function
              if (Input%out_contr) then
                call writectr(Input%folder,iph,ith,Frec%omega,Geom, &
                              Atmo%z,ContrG,Input%lim_ctr)
                if (laborted) goto 2000
              end if

              ! Control
              call control
              if (laborted) goto 2000

            end if

            ! Say completed
            if (gpid.eq.0) then
              umsg = '   Completed'
              call verbose
            end if

          ! Slave in synthesis
          else if (run_mode.ne.-1) then

            call control
            if (laborted) goto 2000

          end if

        enddo ! azimuthal LOS directions
      enddo ! polar LOS directions

      !
      ! Clean slave pointers
      !
2000  if (pid.ne.0) then
        deallocate(data1M,data1O)
        nullify(data1M,data1O)
        nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
        nullify(p_StkM,p_StkO)
        if (associated(etaIM)) then
          deallocate(etaIM)
          nullify(etaIM,tauM)
        end if
      end if

      return

      end subroutine emergence

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
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      subroutine emergence_serial(Atom,LTElines,Atmo,Cont,Frec,Red, &
                                  Bfield,Geom,MPID,Input,Flgsg, &
                                  JKQ_asym,Stokes,JKQ,JKQC,SolF)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Geometry_class):: Geom
      type(Bfield_class):: Bfield
      type(Solution_F_class):: SolF
      double precision, &
        dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
                                                       target:: Stokes
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1):: JKQC

      ! Local

      logical:: AD,ADD

      integer:: ncount,iz0,iz1,diz,m,o,p,op
      integer:: ith,iph,iz,if0,if1,icount

      double precision:: mu_inv,dsm,dsp,dzm,dzp
      double precision, dimension(:,:), allocatable:: tau, tau1
      double precision, dimension(:,:,:), allocatable:: Contr

      ! Pointers
      double precision, dimension(:,:,:), pointer:: data1M,data1O, &
                                                    data1P
      double precision, dimension(:), pointer:: etaIM,etaIO
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP


      !
      ! Initializations
      !

      ! Reset progress counter
      icount = 0

      ! Determine number of directions to do
      ncount = Geom%nThLOS*Geom%nPhLOS

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! Initialize index of Stokes
      op = 1

      ! CPU limits
      if0 = 1
      if1 = nfreq

      ! Allocate M and O pointers for RT coeff
      allocate(data1M(0:3,nfreq,0:5))
      allocate(data1O(0:3,nfreq,0:5))

      ! If calculating height of tau=1, allocate
      if (Input%out_tau1.or.Input%out_contr) then
        allocate(tau(nfreq,Rz0:Rz1))
        allocate(etaIM(nfreq))
        allocate(tau1(2,nfreq))
        allocate(Contr(0:3,nfreq,Rz0:Rz1))
      else
        nullify(etaIM,etaIO)
      end if

      ! If inverting, need to return the output
      if (run_mode.eq.-1) then
        allocate(SolF%e_Stk(0:3,nfreq,Geom%nPhLOS,Geom%nThLOS))
        if (Input%out_tau1) &
        allocate(SolF%e_tau1(Input%lim_tau%nn,Geom%nPhLOS, &
                             Geom%nThLOS))
        if (Input%out_contr) &
        allocate(SolF%e_Ctr(0:3,Input%lim_ctr%nn,Rz0:Rz1, &
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

          !
          ! Get geometry if PRD AD
          !
          if (PRD.and..not.AV) &
            call get_scattering_los(Geom,ith,iph)

          !
          ! Set geometrical tensors
          !
          call setTKQLOS(Geom,Flgsg,Bfield,ith,iph)

          ! Advance the index
          icount = icount + 1

          ! If global Máster
          if (gpid.eq.0) then

            ! Communicate which direction we are doing
            write(umsg,'(A,i4,A,i4)') &
                        '   Doing direction ',icount,' of ',ncount
            call verbose

          end if

          ! If calculating contribution function
          if (Input%out_tau1.or.Input%out_contr) then

            ! Reset tau
            tau = 0
            tau1(1,:) = 0
            tau1(2,:) = Atmo%z(Rz0)

            ! Initial point
            o = Rz0

            ! Calculate absorptivity
            call RTAbs(Frec,Atom,LTElines,Atmo,Flgsg,Geom, &
                       o,ith,iph,if0,if1,Cont%ndir, &
                       Cont%c(:,:,:,o),Bfield,etaIM)

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
              call RTAbs(Frec,Atom,LTElines,Atmo,Flgsg,Geom, &
                         o,ith,iph,if0,if1,Cont%ndir, &
                         Cont%c(:,:,:,o),Bfield,etaIO)

              ! Calculate run of tau
              call RTtau(nfreq,dsm,Atmo%z(m),Atmo%z(o),etaIM, &
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
            call top(MPID,data1M(:,:,5))

          ! If going up, get bottom boundary
          else

            ! Call bottom boundary
            call bottom(Frec%omega,Atmo%T(iz0),Atmo%vx(iz0), &
                        Atmo%vy(iz0),Atmo%vz(iz0), &
                        Geom%L_mu(ith),cos(Geom%L_phi(iph)), &
                        sin(Geom%L_phi(iph)),MPID,data1M(:,:,5))

          endif ! propagation direction

          ! Identify current height
          o = iz0

          ! Index for Stokes
          if (PRD.and.ADD) op = o

          ! Calculate radiative coefficients
          call RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom, &
                        o,ith,iph,if0,if1,JKQ_asym, &
                        JKQ(:,:,:,o),JKQC(:,:,:,o),Cont%ndir, &
                        Cont%c(:,:,:,o),Bfield, &
                        Stokes(:,:,:,:,op),data1M(:,:,0:4))


          ! Identify next height
          p = iz0 + diz

          ! Index for Stokes
          if (PRD.and.ADD) op = p

          ! Calculate radiative coefficients
          call RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom, &
                        p,ith,iph,if0,if1,JKQ_asym, &
                        JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                        Cont%c(:,:,:,p),Bfield, &
                        Stokes(:,:,:,:,p),data1O(:,:,0:4))

          ! Contribution function at first point is 0
          if (Input%out_contr) Contr(:,:,o) = 0d0


          !
          ! Intermedium heights
          !

          ! Medium heights
          do iz=iz0,iz1,diz

            ! We treat the boundaries outside
            if(iz.eq.iz0.or.iz.eq.iz1)cycle

            ! Allocate P pointers
            allocate(data1P(0:3,if0:if1,0:5))

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
            if (PRD.and.ADD) op = p

            ! Calculate radiative coefficients
            call RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom, &
                          p,ith,iph,if0,if1,JKQ_asym, &
                          JKQ(:,:,:,p),JKQC(:,:,:,p),Cont%ndir, &
                          Cont%c(:,:,:,p),Bfield, &
                          Stokes(:,:,:,:,p),data1P(:,:,0:4))

            ! Point to the data
            p_K0M  => data1M(:,:,0)
            p_K1M  => data1M(:,:,1)
            p_K2M  => data1M(:,:,2)
            p_SM   => data1M(:,:,4)
            p_StkM => data1M(:,:,5)
            p_K0O  => data1O(:,:,0)
            p_K1O  => data1O(:,:,1)
            p_K2O  => data1O(:,:,2)
            p_SO   => data1O(:,:,4)
            p_StkO => data1O(:,:,5)
            p_K0P  => data1P(:,:,0)
            p_SP   => data1P(:,:,4)

            ! Apply short characteristics BESSER
            call RTStep(o,ith,iph,nfreq, &
                        dsm,dsp,p_K0M,p_K1M,p_K2M, &
                        p_SM,p_K0O,p_K1O,p_K2O, &
                        p_SO,p_K0P,p_SP,p_StkM, &
                        p_StkO,.True.)

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
              call RTContr(nfreq,dsm,dsp,dzm,dzp,p_K0M, &
                           p_K0O,p_K1O,p_K2O,p_SO,p_K0P, &
                           p_StkO,tau(:,o),Contr(:,:,o),.True.)

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
          p_K0M  => data1M(:,:,0)
          p_K1M  => data1M(:,:,1)
          p_K2M  => data1M(:,:,2)
          p_SM   => data1M(:,:,4)
          p_StkM => data1M(:,:,5)
          p_K0O  => data1O(:,:,0)
          p_K1O  => data1O(:,:,1)
          p_K2O  => data1O(:,:,2)
          p_SO   => data1O(:,:,4)
          p_StkO => data1O(:,:,5)

          ! Apply short characteristics LINEAR
          call RTStep(o,ith,iph,nfreq, &
                      dsm,dsp,p_K0M,p_K1M,p_K2M, &
                      p_SM,p_K0O,p_K1O,p_K2O, &
                      p_SO,p_K0P,p_SP,p_StkM, &
                      p_StkO,.False.)

          ! If calculating contribution function
          if (Input%out_contr) then

            ! Calculate vertical distance between points
            if (ztau) then
              dzm = dsm/mu_inv
            else
              dzm = Atmo%z(o) - Atmo%z(m)
            end if

            ! Calculate contribution function
            call RTContr(nfreq,dsm,dsp,dzm,dzp,p_K0M, &
                         p_K0O,p_K1O,p_K2O,p_SO, &
                         p_K0P,p_StkO,tau(:,o),Contr(:,:,o),.False.)

            ! Inversion
            if (run_mode.eq.-1) then

              ! Keep contribution function
              call setctr(SolF%e_Ctr(:,:,:,iph,ith),Contr, &
                          Input%lim_ctr)

            ! Synthesis
            else

              ! Write contribution function
              call writectr(Input%folder,iph,ith,Frec%omega,Geom, &
                            Atmo%z,Contr,Input%lim_ctr)

            end if ! Inversion
          end if ! calculating contribution function

          ! If inverting
          if (run_mode.eq.-1) then

            ! Keep Stokes
            call setstk(SolF%e_Stk(:,:,iph,ith),p_StkO, &
                        Input%lim_stk,.False.)

          ! Synthesis
          else

            ! Write stokes
            call writestk(Input%folder,iph,ith,Frec%omega, &
                          Geom,p_StkO,Input%lim_stk)

          end if ! Inverting

          ! Communicate we finished this direction
          if (gpid.eq.0) then
            umsg = '   Completed'
            call verbose
          end if

        enddo ! azimuthal LOS directions
      enddo ! polar LOS directions

      ! Clean pointers
      deallocate(data1M,data1O)
      nullify(data1M,data1O)
      nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
      nullify(p_StkM,p_StkO)
      if (associated(etaIM)) then
        deallocate(etaIM)
        nullify(etaIM)
      end if

      return

      end subroutine emergence_serial

!#####################################################################
!#####################################################################
!#####################################################################

      end module solver_mod
