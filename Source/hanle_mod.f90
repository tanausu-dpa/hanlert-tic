      !> Flow control of the code
      module hanle_mod
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
!     04/18/2017
!  Last version:
!     12/12/2023 V3.2.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     12/12/2023:    V3.2.3 - Call diabon_B0 in absence of magnetic
!                             fields (TdPA)
!
!     11/24/2023:    V3.2.2 - Do not normalize for intensity if there
!                             are no iterations to make (TdPA)
!
!     11/14/2023:    V3.2.1 - Bugfix: When failing to read Stokes
!                             from a solution file, checking "lio"
!                             can be insufficient to decide what
!                             to initialize (TdPA)
!
!     10/16/2023:    V3.2.0 - Split the hanle routine into calls
!                             of a number of routines: hanle_setup,
!                             hanle_reback, hanle_init,
!                             hanle_intensity, and
!                             hanle_polarization (TdPA)
!                           - Added the possibility of solving the
!                             polarized problem in two steps, one
!                             without magnetic field to initialize
!                             the magnetic problem (TdPA)
!                           - The line broadening is now computed
!                             when preparing the synthesis (TdPA)
!
!     10/04/2023:    V3.1.5 - Added new type of inversion (TdPA)
!
!     08/30/2023:    V3.1.4 - Added argument to initmemoJ call (TdPA)
!
!     08/28/2023:    V3.1.3 - Added argument to free_local call (TdPA)
!
!     08/17/2023:    V3.1.2 - Moved control call before memory
!                             freeing (TdPA)
!
!     08/07/2023:    V3.1.1 - Ensure the initialization of Stokes and
!                             JKQ continuum when cannot be read from
!                             a solution file (TdPA)
!                           - Added calls and arguments related to LTE
!                             lines (TdPA)
!                           - Moved most of the logic of the
!                             normalization of profiles to the
!                             normalization module (TdPA)
!                           - Moved the selection of the solvers to
!                             the respective modules (TdPA)
!                           - Reset counter for profile memory when
!                             removing normalization (TdPA)
!
!     07/03/2023:    V3.1.0 - Moved the beginning of the hanle
!                             routine which initialized populations,
!                             eq. of state, etc., outside of the
!                             routine, to facilitate the integration
!                             of the inversion module (TdPA)
!                           - Also moved the part about updating the
!                             electron density (TdPA)
!                           - The inversion sends a new argument
!                             where the solution is stored, instead
!                             of calling the reading/writing of a
!                             file (TdPA)
!                           - Added the call to the new chi_freq
!                             subroutine from background_mod (TdPA)
!                           - Added a subroutine to prepare the
!                             solution buffers when running an
!                             inversion (TdPA)
!                           - Changed calls to free memory to
!                             better integrate the inversion
!                             module (TdPA)
!
!     04/11/2023:    V3.0.8 - Cont%chi500 should be deallocated to
!                             redo the background (HL)
!
!     03/08/2023:    V3.0.7 - Do not compute the geometrical tensor
!                             components in the magnetic field
!                             reference frame if not doing
!                             polarization (TdPA)
!                           - Bugfix: omegabuildinI had Geom instead
!                             of GeomI in its arguments (TdPA)
!
!     02/14/2023:    V3.0.6 - Added GeomI exclusive for use in the
!                             intensity problem (TdPA)
!                           - Added option to force AA intensity in a
!                             AD polarization problem (TdPA)
!                           - It may be necessary to recompoute
!                             background if b-b transitions are
!                             present and the quadrature changes
!                             from intensity to polarization (TdPA)
!                           - Added dirty tricks to ignore
!                             velocities in the intensity solution
!                             even for dynamic models (TdPA)
!                           - Added dirty tricks to ignore
!                             the horizontal component of velocities
!                             in the intensity solution if forced to
!                             be axially symmetric (TdPA)
!                           - Do not deallocate Atmo%Pg, Atmo%Pe, and
!                             Atmo%rho in the inversion (HL)
!                           - If doing AD PRD, but the intensity is
!                             AA due to the user input, do not
!                             consider the existing Stokes as AA,
!                             but assume that was AD (TdPA)
!
!     11/10/2022:    V3.0.5 - Added JKQ_asym as argument for
!                             emergence and emergence_serial (TdPA)
!
!     10/25/2022:    V3.0.4 - Added JKQin as argument. This contains
!                             information about ad-hoc JKQ in the
!                             1.5D case (TdPA)
!                           - Added e_Stk as argument. This contains
!                             the emergent Stokes parameters in the
!                             inversion case (TdPA)
!                           - Changed the call to prepareatom to
!                             adjust for its changes (TdPA)
!                           - Changed the call to Initpopu to
!                             adjust for its changes (TdPA)
!                           - Changed the call to cleanFrecandRed to
!                             adjust for its changes (TdPA)
!                           - Added a call to setmpi_sizes to
!                             adjust for the new height limitation
!                             option (TdPA)
!                           - Added a call to new Initcrho (TdPA)
!                           - Changed the calls to  emergenceI,
!                             emergenceI_serial, emergence, and
!                             emergence_serial to add the e_Stk
!                             argument (TdPA)
!                           - Changed explicit cleaning of the
!                             Red structure with a proper
!                             cleaning call (TdPA)
!                           - Added call to restrict_zaxis (TdPA)
!                             adjust for its changes (TdPA)
!                           - Rearranged the order of some calls to
!                             adjust for the new height limitation
!                             option (TdPA)
!
!     07/27/2022:    V3.0.3 - Renamed MPI to MPID (TdPA)
!
!     07/13/2022:    V3.0.2 - Initializenlte now requires the Atmo
!                             argument (TdPA)
!                           - Input%resource is not required in the
!                             calls to eqstate, redo_ne, chemeq, and
!                             broad (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: Atomb cannot be deallocated
!                             unless we are running the single 1D
!                             synthesis case (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o The frequency structure is now an
!                                input and the output frequency
!                                axis is initialized elsewhere.
!                              o Now fudge factor and Kurucz lines
!                                are inputs to avoid reading the files
!                                repeatedly.
!                              o The flags that decide the path in
!                                hanle_mod are now inputs (lio,lie,lp,
!                                and lpe) because they are needed
!                                outside as well.
!                              o Some of the deallocations here are
!                                conditioned by the type of run via
!                                the free logical variable.
!                              o The deallocation of background atoms
!                                and molecules happens elsewhere.
!                              o The definition of the offset limit
!                                for MPI writing happens elsewhere.
!                              o The flagging of the atoms to fix
!                                their populations happens elsewhere.
!                              o Added call to prepareatom  and
!                                preparemol to allocate arrays in the
!                                Atom, Atomb, and Mol structures.
!                              o Added checks of need to abort that
!                                send the routine to the end.
!                              o The routine does not return without
!                                going via the new 1000 label anymore.
!                              o Added call to setTB to generate the
!                                geometrical tensors in the magnetic
!                                field reference frame.
!                              o The atoms are now diagonalized in
!                                this module.
!                              o The initialization of the
!                                photoionization quantities has been
!                                split into two, with the first part
!                                (cross section) initialized
!                                elsewhere and the second part called
!                                here by setphotoTEI.
!                              o The non-error verbosity has been
!                                limited to the global master and thus
!                                it should only happen in the pure 1D
!                                case.
!                              o Initcols is also called for passive
!                                atoms to initialize some arrays.
!                              o Added an argument to the calls to
!                                writecols, writedamp, writeatmo, and
!                                writeback.
!                              o The resizing of the frequency axis
!                                happens elsewhere.
!                              o Frec%exu is now nullified when not
!                                needed at all to avoid trying to
!                                deallocate it when freeing memory.
!                              o Background now also has fudge and
!                                kurucz as inputs.
!                              o Changed MPI communicator from
!                                MPI_COMM_WORLD to MPI_COMM_RT.
!                              o Change the first to arguments for
!                                writesol and writesolI for the Input
!                                structure.
!                              o wAtmo can only be called if in the
!                                pure 1D case.
!                              o Removed a superfluous if clause that
!                                was there to prepare for a change
!                                never implemented.
!                              o Added a call to clean the memory of
!                                height dependent quantities at the
!                                end of the 1D calculation.
!                             (TdPA)
!
!     05/24/2022:    V2.0.2 - Added a warning for when an intensity
!                             solution is loaded to compute emergent
!                             polarized profiles without iterating
!                             before (TdPA)
!
!     04/07/2022:    V2.0.1 - Added an argument to writesol/I (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added calls to new OpenMP related
!                             routines (TdPA)
!                           - Always call initialize_asym (TdPA)
!                           - Changed arguments in calls that do not
!                             need the MPI structure anymore (TdPA)
!
!     02/12/2021:   V1.2.23 - Added calls to report_time (TdPA)
!
!     02/04/2021:   V1.2.22 - Added argument to normalize due to
!                             changes in that routine (TdPA)
!
!     01/13/2021:   V1.2.21 - Call the initialization of asymmetry
!                             ad-hoc inputs and passes the new
!                             arguments into the solver (TdPA)
!
!     11/12/2020:   V1.2.20 - Adds protection agains chemical
!                             equilibrium if indicated (TdPA)
!
!     09/11/2020:   V1.2.19 - Added depar variable and passed to the
!                             relevant subroutines (TdPA)
!                           - Added control of allocated RAM (TdPA)
!                           - Initialize and clean PRD data in Frec
!                             and Red pointers (TdPA)
!                           - Always call normalization before
!                             input axis builders (TdPA)
!                           - Added the possibility to call a RAM
!                             use reporter (TdPA)
!
!     07/31/2020:   V1.2.18 - Always call frecresize (TdPA)
!                           - If not doing MPI and not storing
!                             photoionization quantities, allocate a
!                             dummy pointer Frec%exu to avoid calls
!                             with uninitialized pointers (TdPA)
!
!     06/26/2020:   V1.2.17 - Changed call to check_nchlt (TdPA)
!
!     06/02/2020:   V1.2.16 - Added lio and lp to the arguments of
!                             omegabuildin and omegabuildinI (TdPA)
!
!     06/01/2020:   V1.2.15 - Avoided unnecessary call to omegabuildin
!                             if not iterating but doing emergent
!                             solution (TdPA)
!                           - Added call to a routine to check the
!                             magnetic field regime of the lower
!                             term when using the non-coherent lower
!                             term approximation (TdPA)
!
!     05/11/2020:   V1.2.14 - Added Frec%omega as argument for the
!                             writeback subroutine (TdPA)
!
!     04/14/2020:   V1.2.13 - Bugfix: The call to correctpop before
!                             solverI must be avoided when reading
!                             (TdPA)
!
!     03/18/2020:   V1.2.12 - Removed initialization in case of no
!                             iterations. Radiation quantities get
!                             always initialized in their specific
!                             routine (TdPA)
!
!     03/10/2020:   V1.2.11 - If updating atmos and keeping the
!                             atmospheric data file, it will be
!                             updated also at the end (TdPA)
!
!     03/05/2020:   V1.2.10 - The input lH is not needed now (TdPA)
!                           - Define nlte and call routines to
!                             initialize and activate it (TdPA)
!                           - Changed order of the initial calls in
!                             this routine (TdPA)
!                           - Now populations are initialized in two
!                             steps. First, one reads the available
!                             files, and later the LTE is computed and
!                             the rho^0_0 properly initialized (TdPA)
!                           - Now the code tries to keep the hydrogen
!                             in the atmosphere and in the atom (if
!                             active or populations given) consistent
!                             between them (TdPA)
!                           - The normalization and de-normalization
!                             of populations is done in secondary
!                             routines called here (TdPA)
!                           - If you put iteration limits that result
!                             in no iterations, the solvers are not
!                             even called. In that case, if not
!                             loading, radiation tensors are
!                             initialized to zero (TdPA)
!                           - The code will not do polarization if you
!                             are updating the atmospheric model or
!                             changing the electrons (TdPA)
!
!     02/10/2020:    V1.2.9 - Added calls to routine to control
!                             inputs of hydrogen populations and
!                             for recalculating the electron number
!                             density (TdPA)
!
!     12/10/2019:    V1.2.8 - Now admits no LOS angles (TdPA)
!                           - Moved writeatmo after background,
!                             because it needs it to convert between
!                             heights and tau (TdPA)
!                           - Added routine to transform between tau
!                             and heights (TdPA)
!                           - Adjusted conditions under which JKQgen
!                             is called (TdPA)
!
!     11/19/2019:    V1.2.7 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - Split omegabuildin/I into two calls, one
!                             for quadrature and one for LOS (TdPA)
!                           - Initializes memory counters (TdPA)
!                           - Drops intensity redistribution when
!                             starting polarization calculations
!                             for normalizations (TdPA)
!
!     11/12/2019:    V1.2.6 - Changes to accomodate for the new
!                             normalization split between quadrature
!                             and LOS, as well as the Voigt data file
!                             option (TdPA)
!
!     10/18/2019:    V1.2.5 - Passing a variable between the Input and
!                             Atom structures (TdPA)
!
!     09/26/2019:    V1.2.4 - New calls due to additions to the code
!                             elsewhere (TdPA)
!
!     08/14/2019:    V1.2.3 - Added call to the new setNSCcoeff (TdPA)
!
!     06/03/2019:    V1.2.2 - Bugfix: Do not try to deallocate the
!                             input names if there were not any (TdPA)
!
!     05/31/2019:    V1.2.1 - Changed the call arguments of setmpi and
!                             added the call to setmpi_sizes, in
!                             agreement with changes in the former
!                             routine and the additin of the latter
!                             for the new MPI packages (TdPA)
!
!     05/08/2019:    V1.2.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!                           - Bugfix: Missing space in verbose (TdPA)
!
!     03/12/2019:    V1.1.1 - Added arguments to broad call (TdPA)
!                           - Bugfix: Message for frequency
!                             initialization (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!                           - Now controls allocations and
!                             deallocations (TdPA)
!
!     11/19/2018:   V1.0.26 - Updated writedamp call arguments (TdPA)
!
!     09/06/2018:   V1.0.25 - Added argument to writesolI (TdPA)
!                           - Bugfix: The polarization part must run
!                             with force='ALL' too (TdPA)
!                           - You can now keep the pure intensity
!                             solution even when solving for
!                             polarization (TdPA)
!
!     09/04/2018:   V1.0.24 - setmpi() needs some more inputs (TdPA)
!                           - Branched the calls to the different
!                             solvers due to the implementation of
!                             the alternative versions (TdPA)
!
!     08/08/2018:   V1.0.23 - Added some barriers for message
!                             synchronization (TdPA)
!
!     08/06/2018:   V1.0.22 - Call to new subroutine ramphoto and
!                             additional messages regarding RAM and
!                             profile allocations (TdPA)
!
!     07/10/2018:   V1.0.21 - Need to pass Flgsg to background (TdPA)
!
!     11/27/2017:   V1.0.20 - Changed JKQgen and JKQgen_serial calls
!                             to addapt for the changes within those
!                             subroutines (TdPA)
!
!     10/11/2017:   V1.0.19 - Moved the flags that control the flow
!                             of the code to the beginning of the
!                             subroutine (TdPA)
!                           - Passing Bfield%strength and lp.or.lpe
!                             into omegabuild (TdPA)
!
!     10/03/2017:   V1.0.18 - Changed arguments of JKQgen to be
!                             consistent with its changes (TdPA)
!
!     09/15/2017:   V1.0.17 - Added Input%resource to all the
!                             chemiqeq and broad calls (TdPA)
!
!     08/09/2017:   V1.0.16 - Added lload to solveri_serial and
!                             solveri arguments (TdPA)
!                           - Added a barrier after omegabuildin to
!                             not receive the message too soon (TdPA)
!                           - Force cannot be 'Q', 'P' is the correct
!                             value to force polarization (TdPA)
!
!     07/17/2017:   V1.0.15 - Bugfix: Checks if there are molecules
!                             to deallocate (TdPA)
!                           - Removed debug print (TdPA)
!
!     07/10/2017:   V1.0.14 - Bugfix: Now it is necessary to call
!                             frecresize even if the is no frequency
!                             splitting (TdPA)
!
!     07/06/2017:   V1.0.13 - Bugfix: JKQgen was being called in some
!                             situations when it should not (TdPA)
!
!     07/05/2017:   V1.0.12 - Not normalizing for terms if not
!                             going to calculate polarization (TdPA)
!
!     06/29/2017:   V1.0.11 - Put back resetWarr, but new simplified
!                             version (TdPA)
!
!     06/28/2017:   V1.0.10 - Declaring and passing Red (TdPA)
!                           - Removed the resetWarr and resetWarrI
!                             routines, not needed anymore (TdPA)
!
!     06/26/2017:    V1.0.9 - Bugfix: JKQgen was being called also
!                             for the serial case (TdPA)
!
!     06/23/2017:    V1.0.8 - Removed Atom from initialize and
!                             initializeI calls (TdPA)
!                           - Added argument to setmpi call (TdPA)
!
!     06/22/2017:    V1.0.7 - Now J00P has to be passed around (TdPA)
!                           - Changed how Jbar are initialized (TdPA)
!                           - Removed Atom from initialize and
!                             initializeI calls (TdPA)
!
!     06/19/2017:    V1.0.6 - Added Geom to omegabuildin call (TdPA)
!
!     06/15/2017:    V1.0.5 - Bugfix: omegabulidin was not being
!                             called for polarization (TdPA)
!
!     06/09/2017:    V1.0.4 - Added J iterations (TdPA)
!
!     05/12/2017:    V1.0.3 - Added omegabuildinI (TdPA)
!
!     05/05/2017:    V1.0.2 - Bugfix: Wrong condition for
!                             initialization (TdPA)
!
!     05/02/2017:    V1.0.1 - Bugfix: Only master should store
!                             a solution after ionization balance and
!                             that block should not be called if
!                             forcing intensity (TdPA)
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
!  hanle:
!    Manage the workflow to solve a RT problem
!
!  hanle_setup:
!    Set-up preliminar steps for any formal solution
!
!  hanle_reback:
!    Recalculate the background continuum quantities
!
!  hanle_init:
!    Initialize the RT problem
!
!  hanle_intensity:
!    Solve the NLTE problem for intensity
!
!  hanle_polarization:
!    Solve the polarized NLTE problem
!
!  prepare_buffers:
!    Allocate the arrays to store the solution of the forward problem
!  and the configure the synthesis mode for the inversion code,
!  depending on if the synthesis is a trial or a response function
!
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use background_mod
      use commons_mod
      use diagon_mod
      use free_mod
      use gauss_mod
      use getztau_mod
      use initialize_mod
      use initmemoj_mod
      use initphotoion_mod
      use initpopu_mod
      use iosolution_mod
      use normalizer_mod
      use omegabuild_mod
      use parameters_mod , only : TINYB
      use setmpi_mod
      use solver_mod
      use solveri_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Secondary main with that controls the execution flow.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!            Mol(Mol_class): Structure with the molecule data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!     GeomI(Geometry_class): Structure with geometry data for
      !!                            the intensity problem\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!          JKQin(double(:)): Data with JKQ asymmetries\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      !!            lload(logical): Reading solution file\n
      !!              lio(logical): Doing intensity formal solution\n
      !!              lie(logical): Doing intensity emergence\n
      !!               lp(logical): Doing polarized formal solution\n
      !!              lpe(logical): Doing polarized emergence\n
      !!             free(logical): If allowed to free some input
      !!                            data
      subroutine hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input, &
                       GeomI,Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                       JKQin,SolF,lload,lio,lie,lp,lpe,free)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI,Geom
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Solution_F_class):: SolF
      logical, intent(in)::lload,lio,lie,lp,lpe,free
      double precision, dimension(:), allocatable:: JKQin

      ! Local

      type(Bfield_class):: Bfield0
      type(Continuum_class):: Cont
      type(Rhoc_class), dimension(nA):: Rho_old

      logical:: rlimw = .True.
      logical:: csize
      logical:: rback,rdyn,raxial
      logical:: l1, l2, ofram

      integer:: iph,rnPh
      integer, dimension(:), allocatable:: nlte,depar

      double precision:: rVPhi,rVmux,rVmuy,rWmux

      double precision, dimension(:,:,:,:), allocatable:: StokesI
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P

      complex(kind=8), dimension(:,:,:), allocatable:: JKQ_asym
      complex(kind=8), dimension(:,:,:), allocatable:: JKQ_asym_fake

      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Initialize frequency and redistribution pointer
      nullify(Frec%dzao)

      ! If forcing the problem to be static for intensity, but it
      ! is dynamic
      if (dyn.and.Input%static_int) then

        ! Trick the problem
        rdyn = dyn
        dyn = .False.

        Atmo%vxa => Atmo%vx
        Atmo%vya => Atmo%vy
        Atmo%vza => Atmo%vz
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros
        Atmo%vz => Atmo%zeros

      ! If axial intensity and not polarization
      else if (axiali.and..not.axial) then

        ! Cheat the velocities
        Atmo%vxa => Atmo%vx
        Atmo%vya => Atmo%vy
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros

      ! Otherwise, null
      else

        ! Initialize
        nullify(Atmo%vxa)
        nullify(Atmo%vya)
        nullify(Atmo%vza)

      end if

      ! Call initialization
      call hanle_setup(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input, &
                       GeomI,Geom,Bfield,Frec,Flgsg,fudge, &
                       kurucz,SolF,Cont,Rho_old,free,rback)

      ! Control
      if (laborted) goto 1000


      !
      ! Check if we need to reevaluate sizes
      !
      if (Rnz.ne.nz) then
        l1 = .not.lload
        l2 = (lio.and..not.lie).or.(lload.and.(lp.or.lpe))
        if (MPID%mpi) &
          call setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,l1,l2,.True.)
      end if

      ! Control
      if (laborted) goto 1000

      !
      ! Initialize solution
      !
      call hanle_init(Atom,Atmo,MPID,Input,GeomI,Geom,Bfield, &
                            Frec,Flgsg,SolF,lload,lio,lie,lp,lpe, &
                            Stokes,JKQ,JKQS,JKQC, &
                            StokesI,J00,J00S,J00C,J00P)

      ! Control
      if (laborted) goto 1000


      !
      ! Intensity formal solution
      if (lio.or.lie) &
        call hanle_intensity(Atom,LTElines,Atmo,MPID,Input, &
                             GeomI,Bfield,Frec,Flgsg,SolF, &
                             Cont,Rho_old, &
                             StokesI,J00,J00S,J00C,J00P, &
                             lload,lio,lie,rlimw,ofram)

      ! Control
      if (laborted) goto 1000


      !
      ! Check if updating the model
      !
      if (Input%update_atmos.ge.0) then

        ! Check type of calculation
        if (pid.eq.0.and.(lp.or.lpe)) then
          umsg = 'The option update_atmos is not compatible '// &
                 'with computing the polarization. The code '// &
                 'is stopping after updating the atmosphere'
          call verbose
        end if

        ! Leave
        goto 1000

      else

        if (Input%redo_ne.eq.1.or.Input%redo_ne.eq.11) then
          if (pid.eq.0) then
            umsg = ' # Warning: You chose to redo electrons '// &
                   'at the end of iterations, but not '// &
                   'to update the atmosphere. Electrons are '// &
                   'not recomputed.'
            call verbose
          end if
        end if
      end if

      ! Deallocate atmospheric quantities if still in RAM
      if (run_mode.ne.-1) then
        if (allocated(Atmo%Pg)) deallocate(Atmo%Pg)
        if (allocated(Atmo%Pe)) deallocate(Atmo%Pe)
        if (allocated(Atmo%rho)) deallocate(Atmo%rho)
      end if
      if (allocated(nlte)) deallocate(nlte)
      if (allocated(depar)) deallocate(depar)
      if (free.and..not.rback) then
        if (allocated(Atomb)) then
          nAb = 0
          deallocate(Atomb)
        end if
      end if

      ! If forcing the problem to be static for intensity, but it
      ! is dynamic
      if (rdyn.and.Input%static_int) then

        ! Return the problem to its original state
        dyn = rdyn

        Atmo%vx => Atmo%vxa
        Atmo%vy => Atmo%vya
        Atmo%vz => Atmo%vza
        nullify(Atmo%vxa,Atmo%vya,Atmo%vza)

      ! If axial intensity and not polarization, un-cheat
      ! the velocities
      else if (axiali.and..not.axial) then

        Atmo%vx => Atmo%vxa
        Atmo%vy => Atmo%vya
        nullify(Atmo%vxa,Atmo%vya)

      end if


      !
      ! Polatization
      if (lp.or.lpe) then

        ! Initialize extra asymmetry
        call initialize_asym(Input,MPID,Flgsg,JKQin,JKQ_asym)

        !
        ! Update radiation RAM if coming from intensity
        !

        if ((lio.and..not.lie).or. &
            (lload.and..not.allocated(Stokes))) then

          ! Remove current size in radiation
          MPID%RAM = MPID%RAM - MPID%RRAM

          ! Pre-compute amount of RAM to fill with radiation
          if (AV.or..not.PRD) then
            MPID%RRAM = 8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*2)
          else
            MPID%RRAM = 8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*nz)
          end if
          MPID%RRAM = MPID%RRAM + 8d-6*dble(nz*(nxphot*2 + &
                                            15*2*(2*nxtran + nfreq)))

          ! Update RAM needed
          MPID%RAM = MPID%RAM + MPID%RRAM

        end if ! Refit RAM for radiation

        ! Redo-background?
        if (rback) &
          call hanle_reback(Atom,Atomb,Mol,Atmo,MPID,Input,Geom, &
                            Frec,Flgsg,fudge,kurucz,Cont,free)

        !
        ! If there is a magnetic field, and doing a two-step
        ! calculation
        if (maxval(Bfield%Bstrength).gt.TINYB.and. &
            Input%two_step_pol) then

          ! Set zero field structure
          allocate(Bfield0%Bstrength(nz))
          allocate(Bfield0%Btheta(nz))
          allocate(Bfield0%Bphi(nz))
          Bfield0%Bstrength = 0d0
          Bfield0%Btheta = 0d0
          Bfield0%Bphi = 0d0

          ! Get current status
          raxial = axial
          rnPh = Geom%nPh
          rVPhi = Geom%V_phi(1)
          rVmux = Geom%V_mux(1)
          rVmuy = Geom%V_muy(1)
          rWmux = Geom%W_mux(1)
          csize = .False.

          if(gpid.eq.0) then
            umsg = ' - Solving the polarized but '// &
                   'non-magnetic problem'
            call verbose
          end if

          ! Check geometry if intensity was static or axial
          if (Input%static_int.or.axiali) then

            ! Check axial velocity
            if (maxval(abs(Atmo%vx)).le.0d0.and. &
                maxval(abs(Atmo%vy)).le.0d0) then

              ! And fake it
              axial = .True.
              Geom%nPh = 1
              Geom%V_phi(1) = 0d0
              Geom%V_mux(1) = 1d0
              Geom%V_muy(1) = 1d0
              Geom%W_mux(1) = 1d0

              ! Fake sizes
              csize = .True.
              MPID%size4 = MPID%size4/rnPh
              MPID%size5 = MPID%size5/rnPh

            end if ! No horizontal velocity
          end if ! Intensity was axial or static

          if (axial) then

            ! Call polarization solution WITHOUT field or
            ! emergence with fake JKQ
            call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                    Geom,Bfield0,Frec,Flgsg,SolF, &
                                    Cont,Rho_old, &
                                    StokesI,J00,J00S,J00C,J00P, &
                                    Stokes,JKQ,JKQS,JKQC, &
                                    JKQ_asym_fake,rnPh,.False., &
                                    lload,lio,lie,lp,.False., &
                                    rlimw,ofram)

            ! If will not be axial
            if (.not.raxial) then

              ! Copy Stokes rest of azimuth
              do iph=2,rnPh
                Stokes(:,:,iph,:,:) = Stokes(:,:,1,:,:)
              end do

            end if ! Will not be axial

          else

            ! Call polarization solution WITHOUT field
            call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                    Geom,Bfield0,Frec,Flgsg,SolF, &
                                    Cont,Rho_old, &
                                    StokesI,J00,J00S,J00C,J00P, &
                                    Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                    rnPh,.False., &
                                    lload,lio,lie,lp,.False., &
                                    rlimw,ofram)
          end if

          ! Restore
          axial = raxial
          Geom%nPh = rnPh
          Geom%V_phi(1) = rVPhi
          Geom%V_mux(1) = rVmux
          Geom%V_muy(1) = rVmuy
          Geom%W_mux(1) = rWmux
          if (csize) then
            MPID%size4 = MPID%size4*rnPh
            MPID%size5 = MPID%size5*rnPh
          end if

          ! Solve the actual polarization problem
          call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                  Geom,Bfield,Frec,Flgsg,SolF, &
                                  Cont,Rho_old, &
                                  StokesI,J00,J00S,J00C,J00P, &
                                  Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                  Geom%nPh,.True., &
                                  lload,.False.,.False.,lp,lpe, &
                                  rlimw,ofram)

        ! Normal run
        else

          ! Solve the polarization problem
          call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                  Geom,Bfield,Frec,Flgsg,SolF, &
                                  Cont,Rho_old, &
                                  StokesI,J00,J00S,J00C,J00P, &
                                  Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                  Geom%nPh,.True., &
                                  lload,lio,lie,lp,lpe,rlimw,ofram)

        end if ! Zero field first step solution
      end if ! Polarization


      !
      ! Clean memory
      !

      ! Clean rest
1000  call free_local(Atom,Cont,Geom,Frec,LTElines)

      return

      end subroutine hanle

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up some preliminar steps for the formal solution.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!            Mol(Mol_class): Structure with the molecule data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!     GeomI(Geometry_class): Structure with geometry data for
      !!                            the intensity problem\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!          JKQin(double(:)): Data with JKQ asymmetries\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!             free(logical): If allowed to free some input
      !!                            data\n
      !!            rback(logical): Indicate if background needs to
      !!                            be recalculated here
      subroutine hanle_setup(Atom,Atomb,LTElines,Mol,Atmo,MPID, &
                             Input,GeomI,Geom,Bfield,Frec,Flgsg, &
                             fudge,kurucz,SolF,Cont,Rho_old,free, &
                             rback)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Continuum_class):: Cont
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI,Geom
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Solution_F_class):: SolF
      type(Rhoc_class), dimension(:):: Rho_old
      logical, intent(in)::free
      logical, intent(inout):: rback

      ! Local

      integer:: ia,ios


      ! Routine name
      urou = 'hanle_init'

      ! Control
      if (laborted) return

      ! Reset memory counters
      MPID%RAM = 0d0
      MPID%BRAM = 0d0
      MPID%PRAM = 0d0
      MPID%VRAM = 0d0
      MPID%WRAM = 0d0
      MPID%RRAM = 0d0

      ! Add size of SolF to radiation
      MPID%RRAM = sizeof(SolF)


      !
      ! Prepare OpenMP
      !
#ifdef _OPENMP
      call setomp
#endif

      ! Allocate space for chi scale
      if (.not.allocated(Atmo%chi500)) &
        allocate(Atmo%chi500(nz))

      !
      ! Calculate background continuum quantities
      !
      call background(Atom,Atomb,Mol,Atmo,fudge,kurucz, &
                      Input,Frec%omega,Cont,GeomI,MPID,Flgsg)

      ! Control
      if (laborted) return

      if (gpid.eq.0) then
        umsg = ' - Background quantities calculated'
        call verbose
      end if

      !
      ! Do we need to repeat the background calculation?
      !

      ! Initialize
      rback = .False.

      ! If bound-bound transitions in background
      if (Input%addbb) then

        ! If slave or not mpi
        if (.not.MPID%mpi.or.pid.gt.0) then

          ! If size of directions larger than 1
          if (size(Cont%c,3).gt.1) then

            ! If geometry changes between quadratures
            ! then we repeat
            rback = Geom%nPh.ne.GeomI%nPh.or. &
                    Geom%nTh.ne.GeomI%nTh
          end if
        end if

        ! Anyone found any true?
        if (MPID%mpi) &
          call MPI_ALLREDUCE(MPI_IN_PLACE,rback,1,MPI_LOGICAL, &
                             MPI_LOR,MPI_COMM_RT,ios)

      end if !B-B background

      !
      ! Calculate continuum opacity at reference frequency
      !
      call chi_freq(Atom,Atomb,Mol,Atmo,fudge,Input, &
                    Atmo%tfreq,Atmo%chi500,1,nz,MPID%mpi)

      ! Remove space for chi scale
      if (pid.eq.0.and.MPID%mpi) deallocate(Atmo%chi500)

      ! Control
      if (laborted) return


      !
      ! Compute missing height or tau
      !
      call getztau(Atmo,MPID,.True.)


      !
      ! Store atmosphere in ASCII
      !
      if (Input%keep_atmo) call writeatmo(Atmo,Bfield, &
                                          Input%folder, &
                                          Input%lim_atmo)

      ! Control
      if (laborted) return

      ! Check if we can drop chi500
      if (allocated(Atmo%chi500)) then
        if (MPID%mpi.and.pid.eq.0.and. &
            (Input%update_atmos.lt.0.and..not.Input%keep_atmo)) then
          deallocate(Atmo%chi500)
        else
          ! If not tau scale, drop it
          if (.not.ztau.and.Input%update_atmos.lt.0.and. &
              .not.Input%keep_atmo) &
            deallocate(Atmo%chi500)
        end if ! If Master in MPI
      end if ! Allocated chi500


      !
      ! You can remove pressures now, if not using them later
      !
      if (allocated(Atmo%Pg).and.Input%redo_ne.le.0.and. &
          (Input%update_atmos.lt.0.and.Atmo%typo.eq.0)) then
        ! The inversion needs this
        if (run_mode.ne.-1) then
          deallocate(Atmo%Pg)
          deallocate(Atmo%Pe)
          deallocate(Atmo%rho)
        end if
        if (free.and..not.rback) then
          nAb = 0
          deallocate(Atomb)
        end if
      end if

      ! Deallocate background atoms and molecules
      if (free.and..not.rback) then
        nM = 0
        if (allocated(Mol)) deallocate(Mol)
      end if

      ! Control
      call control

      ! Control
      if (laborted) return

      ! If keeping background, call writer
      if (Input%keep_back) call writeback(Cont,Frec%omega, &
                                          Input%folder,MPID, &
                                          Input%lim_back)

      ! Control
      if (laborted) return

      !
      ! Restrict heights
      !
      call restrict_zaxis(Atmo,Input,MPID)

      !
      ! If there are LTE lines
      !
      if (Input%nLTE.gt.0) &
        call restrict_LTE_lines(Atmo,LTElines,MPID)

      !
      ! Photoionization quantites, thermal part
      !
      do ia=1,nA
        call setphotoTEI(Atom(ia),Frec,Atmo%T,Atmo%ne,MPID,free)
      end do

      ! Control
      if (laborted) return

      ! Message
      if (gpid.eq.0) then
        umsg = ' - Initialized photoionization quantities (thermal)'
        call verbose
      end if

      ! Initialize density matrix
      do ia=1,nA
        call Initcrho_old(Atom(ia),Rho_old(ia))
      end do

      !
      ! Store quantities for epsIphoto
      !
      if (PIRAM) then
        call ramphoto(Atom,Frec,Atmo%T,MPID)
      else
        MPID%PRAM = 0d0
        ! Allocate a dummy index to avoid memory 'leak'
        if (.not.MPID%mpi) then
          allocate(Frec%exu(1,1))
        else
          nullify(Frec%exu)
        end if
      end if

      end subroutine hanle_setup

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up some preliminar steps for the formal solution.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!            Mol(Mol_class): Structure with the molecule data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!             free(logical): If allowed to free some input
      !!                            data
      subroutine hanle_reback(Atom,Atomb,Mol,Atmo,MPID,Input,Geom, &
                              Frec,Flgsg,fudge,kurucz,Cont,free)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(Continuum_class):: Cont
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(Frequency_class):: Frec
      type(Geometry_class):: Geom
      type(Input_class):: Input
      type(MPI_class):: MPID
      logical, intent(in)::free


      ! Clean
      if (allocated(Cont%c)) deallocate(Cont%c)

      ! Do it again
      call background(Atom,Atomb,Mol,Atmo,fudge,kurucz, &
                      Input,Frec%omega,Cont,Geom,MPID,Flgsg)
      ! Control
      if (laborted) return

      if (gpid.eq.0) then
        umsg = ' - Background quantities re-calculated'
        call verbose
      end if

      if (free) then
        nAb = 0
        deallocate(Atomb)
        nM = 0
        if (allocated(Mol)) deallocate(Mol)
      end if

      end subroutine hanle_reback

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the RT problem.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!     GeomI(Geometry_class): Structure with geometry data for
      !!                            the intensity problem\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      !!            lload(logical): Reading solution file\n
      !!              lio(logical): Doing intensity formal solution\n
      !!              lie(logical): Doing intensity emergence\n
      !!               lp(logical): Doing polarized formal solution\n
      !!              lpe(logical): Doing polarized emergence\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!  StokesI(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine hanle_init(Atom,Atmo,MPID,Input,GeomI,Geom,Bfield, &
                            Frec,Flgsg,SolF,lload,lio,lie,lp,lpe, &
                            Stokes,JKQ,JKQS,JKQC, &
                            StokesI,J00,J00S,J00C,J00P)
      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Fctsg_class):: Flgsg
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI,Geom
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Solution_F_class):: SolF
      logical, intent(in)::lload,lio,lie,lp,lpe
      double precision, dimension(:,:,:,:), allocatable:: StokesI
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      logical:: read_stokes

      integer:: ia


      !
      ! Load a previous solution
      !
      if (lload) then

        ! Inversion
        if (run_mode.eq.-1) then

          ! Get solution
          call getsol(SolF,GeomI,Geom,MPID,Flgsg,Bfield, &
                      Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                      J00,J00S,J00C,J00P,lio.or.lie)

          if (gpid.eq.0) then
            umsg = ' - Copied from solution'
            call verbose
          end if

        ! Not inversion
        else

          if (gpid.eq.0) then
            umsg = ' - Reading solution'
            call verbose
          end if

          ! Get solution
          call readsol(Input%solution,GeomI,Geom,MPID,Flgsg,Bfield, &
                       Atom,read_stokes, &
                       Stokes,JKQ,JKQS,JKQC, &
                       StokesI,J00,J00S,J00C,J00P)

          ! Control
          if (laborted) return

          if (lload.and.lpe.and..not.lp.and. &
              .not.lio.and..not.allocated(JKQ)) then
            umsg = ' # Warning: Loading an intensity solution '// &
                   'to compute polarized emergent profiles '// &
                   'without iterating'
            call verbose
          end if

          if (lio.and..not.allocated(J00)) then
            umsg = 'Cannot read multiterm solution as '// &
                   'initialization for multilevel'
            call aborted
            return
          end if

          ! If could not read Stokes
          if (.not.read_stokes) then

            ! Polarization
            if (allocated(JKQC)) then

              ! Initialize
              call initialize_failread(Frec,Atmo,Stokes,JKQC)

            ! Intensity
            else

              ! Initialize
              call initializeI_failread(Frec,Atmo,StokesI,J00C)

            end if ! Type of solution
          end if ! Could read Stokes
        end if ! Inversion/synthesis


      !
      ! Initialize radiation quantities
      !
      else

        if ((lio.or.lie).and..not.lload) then
          call initializeI(Frec,GeomI,Atmo,MPID, &
                           StokesI,J00,J00S,J00C,J00P,Input%mode)
        else if ((lp.or.lpe).or.lload) then
          call initialize(Frec,Geom,Atmo,MPID, &
                          Stokes,JKQ,JKQS,JKQC,J00P,Input%mode)
        end if

        ! Control
        if (laborted) return

        if (gpid.eq.0) then
          umsg = ' - Radiation Field Initialized'
          call verbose
        end if

        ! For every atom, normalize populations
        do ia=1,nA
          call correctpop(Atom(ia),0)
        end do

        ! Control
        if (laborted) return

      end if

      ! Initialize RAM
      MPID%RAM = MPID%PRAM + MPID%BRAM + MPID%RRAM

      end subroutine hanle_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the intensity RT NLTE problem.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!     GeomI(Geometry_class): Structure with geometry data for
      !!                            the intensity problem\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!  StokesI(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      !!            lload(logical): Reading solution file\n
      !!              lio(logical): Doing intensity formal solution\n
      !!              lie(logical): Doing intensity emergence\n
      !!            rlimw(logical): Write RAM limit message\n
      !!            ofram(logical): Indicates if out of RAM
      subroutine hanle_intensity(Atom,LTElines,Atmo,MPID,Input, &
                                 GeomI,Bfield,Frec,Flgsg,SolF, &
                                 Cont,Rho_old, &
                                 StokesI,J00,J00S,J00C,J00P, &
                                 lload,lio,lie,rlimw,ofram)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Fctsg_class):: Flgsg
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Solution_F_class):: SolF
      type(Continuum_class):: Cont
      type(Rhoc_class), dimension(:):: Rho_old
      logical, intent(in):: lload,lio,lie
      logical, intent(inout):: rlimw,ofram
      double precision, dimension(:,:,:,:), allocatable:: StokesI
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P

      ! Local

      type(Red_class):: Red

      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC


      ! Initialize redistribution pointer
      nullify(Red%dzao)

      !
      ! Formal solution and iterating
      !
      if (lio.and.Input%iteri_max.ge.Input%iteri_min) then

        ! Normalize first order
        call normalization(Atom,LTElines,Atmo,Atmo%zeros,GeomI, &
                           Frec,Input,Flgsg,MPID,rlimw,lio, &
                           .False.,.False.)

        ! Control
        if (laborted) goto 1000


        !
        ! PRD
        !

        ! If doing PRD, determine the input frequencies and weights
        if(Input%iter_ord.eq.2)then

          call omegabuildinI(Frec,Red,Atom,Atmo,Input,GeomI,MPID, &
                             lio,ofram,.False.)

          ! Control
          if (laborted) goto 1000

          ! If storing profiles
          if (IRAM) then

            ! If CPU went above ram
            if (ofram.and.rlimw) then

              write(umsg,'(A,i3,A)') ' # Processor',pid, &
                                     ' reached the limit of '// &
                                     'redistribution allocations.'
              call verbose
              rlimw = .False.

            end if

            call MPI_BARRIER(MPI_COMM_RT, ierr)

          end if ! Storing profiles

          if (gpid.eq.0) then
            umsg = ' - Input frequency axis initialized (intensity)'
            call verbose
          end if

#ifdef _OPENMP
          call setomp_2ord(Atom,Frec,GeomI,Bfield%Bstrength, &
                           .False.,.False.)
          ! Control
          if (laborted) goto 1000
#endif

        end if ! doing PRD

      end if ! Intensity formal solution

      ! If not loading
      if (.not.lload) then

        ! If J iterations specified
        if (Input%iter_j.gt.0) then

          ! If doing only intensity
          if (lio) then

            if (gpid.eq.0) then
              umsg = ' - Iterating radiation field'
              call verbose
            end if

            ! Solve
            call solveJ(Atmo,Cont,Frec,GeomI,MPID,Input, &
                        StokesI,J00C)

            ! Control
            if (laborted) goto 1000

          end if ! Only intensity
        end if ! J iterations
      end if ! Not load

      !
      ! Intensity formal solution
      !
      if (lio) then

        ! If iterating
        if (Input%iteri_max.ge.Input%iteri_min) then

          ! Solve the NLTE problem of the first kind
          if (gpid.eq.0) then
            umsg = ' - Starting intensity iterations'
            call verbose
          end if

          ! If we have to report RAM use
          if (Input%RAMreport) call RAMreport(MPID,Input%folder,0,1)

          ! If measuring performance
          if (Input%g_perf.and.pid.eq.0) &
            call report_time(Input%folder,Input%ID,.True.)

          ! Solution
          call solveI(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                      Red,GeomI,MPID,Input,lload,StokesI, &
                      J00,J00S,J00C,J00P)

          ! Control
          if (laborted) goto 1000

          ! If measuring performance
          if (Input%g_perf.and.pid.eq.0) &
            call report_time(Input%folder,Input%ID,.True.)

          ! If intensity was AA and pol is AD, treat as AD
          if (PRD.and..not.AV.and.AVI) tbAD = .False.

        end if ! Intensity iterations
      end if ! lio

      !
      ! Intensity emergence
      !

      if (lie) then

        ! Write the solution file
        if (gpid.eq.0) then
          umsg = ' - Saving intensity solution'
          call verbose
        end if

        ! If inversion
        if (run_mode.eq.-1) then

          ! Set solution
          if (SolF%keep_solution) &
            call setsol(SolF,Flgsg,Bfield, &
                        Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                        J00,J00S,J00C,J00P,.True.)

        ! If synthesis
        else

          ! Write to file
          call writesolI(Input,'NONE',Frec%omega,GeomI,Atom, &
                         Atmo%z,StokesI,J00,J00S,J00C,J00P,.False.)

        end if

        ! Control
        if (laborted) goto 1000

        ! If removing redistribution
        if (.not.(AVI.and..not.dyn.and.lio).or.(.not.AVI.or.dyn)) &
          MPID%RAM = MPID%RAM - MPID%WRAM

        ! If removing normalization
        if (.not.(.not.dyn.and.lio)) then
          MPID%RAM = MPID%RAM - MPID%VRAM
          MPID%VRAM = 0d0
        end if

        !
        ! Normalization for emergence
        !
        call normalization(Atom,LTElines,Atmo,Atmo%zeros,GeomI, &
                           Frec,Input,Flgsg,MPID,rlimw,lio, &
                           .False.,.True.)

        ! Control
        if (laborted) goto 1000


        !
        ! If doing PRD, determine the input frequencies and weights
        !

        ! Define the input frequency axis
        if(Input%iter_ord.eq.2)then

          call omegabuildinI(Frec,Red,Atom,Atmo, &
                             Input,GeomI,MPID,lio,ofram,.True.)

          ! Control
          if (laborted) goto 1000

          ! If storing profiles
          if (IRAM) then

            ! If CPU went above ram
            if (ofram.and.rlimw) then

              write(umsg,'(A,i3,A)') ' # Processor',pid, &
                                     ' reached the limit of '// &
                                     'redistribution allocations.'
              call verbose
              rlimw = .False.

            end if

            call MPI_BARRIER(MPI_COMM_RT, ierr)

          end if ! Storing profiles

          if (gpid.eq.0) then
            umsg = ' - Input frequency axis initialized (intensity)'
            call verbose
          end if

#ifdef _OPENMP
          call setomp_2ord(Atom,Frec,GeomI,Bfield%Bstrength, &
                           .False.,.True.)

          ! Control
          if (laborted) return
#endif

        end if ! PRD

        ! Calculate emergent solutions
        if (gpid.eq.0) then
          umsg = ' - Emergent intensity'
          call verbose
        end if

        ! If we have to report RAM use
        if (Input%RAMreport) call RAMreport(MPID,Input%folder,0,0)

        ! If measuring performance
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

        ! Call solver
        call emergentI(Atom,LTElines,Atmo,Cont,Frec,Red,GeomI, &
                       MPID,Input,StokesI,J00,J00C,SolF)

        ! Control
        if (laborted) goto 1000

        ! If measuring performance
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

      ! If we are doing polarization but we want to keep the
      ! intensity solution file, and we have a solution to store
      elseif (Input%keepIsol.and.lio) then

        ! If synthesis
        if (run_mode.ge.0) then

          ! Write the solution file
          if (gpid.eq.0) then
            umsg = ' - Saving intensity multi-level solution'
            call verbose
          end if

          ! Write to file
          call writesolI(Input,'NONE',Frec%omega,GeomI,Atom, &
                         Atmo%z,StokesI,J00,J00S,J00C,J00P,.True.)

        end if

        ! Control
        if (laborted) goto 1000

      end if

      ! Clean Frec and Red structures
1000  call cleanFrecandRed(Frec,Red,MPID)
      MPID%RAM = MPID%RAM - MPID%WRAM
      MPID%WRAM = 0d0

      return

      end subroutine hanle_intensity

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the polarized RT NLTE problem.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with settings data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!    SolF(Solution_F_class): Class to save the RT solution in
      !!                            RAM\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!       Rho_old(Rhoc_class): Structure to store rhoKQ
      !!                            quantities\n
      !!  StokesI(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !! JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                            tensors\n
      !!             rnPh(integer): Allocation size for Stokes\n
      !!           saving(logical): If creating a solution file\n
      !!            lload(logical): Reading solution file\n
      !!              lio(logical): Doing intensity formal solution\n
      !!              lie(logical): Doing intensity emergence\n
      !!               lp(logical): Doing polarized formal solution\n
      !!              lpe(logical): Doing polarized emergence\n
      !!            rlimw(logical): Write RAM limit message\n
      !!            ofram(logical): Indicates if out of RAM
      subroutine hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                    Geom,Bfield,Frec,Flgsg,SolF, &
                                    Cont,Rho_old, &
                                    StokesI,J00,J00S,J00C,J00P, &
                                    Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                    rnPh,saving, &
                                    lload,lio,lie,lp,lpe,rlimw,ofram)
          
      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Fctsg_class):: Flgsg
      type(Frequency_class):: Frec
      type(Geometry_class):: Geom
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(Solution_F_class):: SolF
      type(Continuum_class):: Cont
      type(Rhoc_class), dimension(:):: Rho_old
      logical, intent(in):: saving,lload,lio,lie,lp,lpe
      logical, intent(inout):: rlimw,ofram
      integer, intent(in):: rnPh
      double precision, dimension(:,:,:,:), allocatable:: StokesI
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC
      complex(kind=8), dimension(:,:,:), allocatable:: JKQ_asym

      ! Local

      type(Red_class):: Red

      logical:: l1

      integer:: ia


      ! Initialize redistribution pointer
      nullify(Red%dzao)

      !
      ! Prepare geometrical tensors
      !

      ! Get B geometrical tensors
      call setTB(Geom,Flgsg,Bfield)

      !
      ! Diagonalize Hamiltonian
      !

      ! For each atom, if there is magnetic field
      if (maxval(Bfield%Bstrength).gt.TINYB) then
        do ia=1,nA
          call diagon(Atom(ia),Bfield,Input%zeeman_mode,Flgsg)
        end do
        if(gpid.eq.0) then
          umsg = ' - Hamiltonian diagonalized'
          call verbose
        end if
      else
        do ia=1,nA
          call diagon_B0(Atom(ia))
        end do
      end if

      !
      ! Normalize
      !

      ! Remove Voigt from RAM
      MPID%RAM = MPID%RAM - MPID%VRAM

      ! Normalize first order
      call normalization(Atom,LTElines,Atmo,Bfield%Bstrength,Geom, &
                         Frec,Input,Flgsg,MPID,rlimw,lp, &
                         .True.,.False.)

      ! Control
      if (laborted) goto 1000

#ifdef _OPENMP
      ! Initialize memoization to use OpenMP later for each atom
      call initmemoJ(Atom,LTElines,Flgsg,Bfield%Bstrength,lp)

      ! Initialize OpenMP for magnetic routines
      call setomp_magn(Atom,Bfield%Bstrength)
#endif

      ! If we solved in intensity, convert multiterm Jbar into
      ! multilevel. Also, if we read a solution only in intensity
      if ((lio.and..not.lie.and.(lp.or.lpe)).or. &
          (lload.and.(lp.or.lpe).and..not.allocated(Stokes))) then

        ! Deal with the NCHLT when coming from intensity
        l1 = NCHLT
        if (NCHLT) NCHLT = .False.

        ! Call
        call JKQgenerate(Atom,Rho_old,Atmo,Frec,Geom, &
                         MPID,Flgsg,Input%Pcorr,Bfield,rnPh, &
                         StokesI,J00,J00S,J00C, &
                         Stokes,JKQ,JKQS,JKQC,J00P)

        ! Control
        if (laborted) goto 1000

        ! Only if synthesis
        if (run_mode.ne.-1) then

          if (gpid.eq.0) then
            umsg = ' - Saving intensity solution'
            call verbose
          end if

          ! Write the solution file
          call writesol(Input,'NONE',Frec%omega,Geom,Flgsg,Bfield, &
                        Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)

        end if ! Synthesis

        ! Deal with the NCHLT when coming from intensity
        NCHLT = l1

      end if


      ! If PRD
      if (Input%iter_ord.eq.2) then

        ! If non-coherent lower term, check the critical field
        ! condition
        if (NCHLT) call check_nchlt(Atom,JKQ,Bfield)

      end if ! PRD


      !
      ! Formal Solution
      !
      if (lp) then

        !
        ! Build input axis
        !

        ! Define the input frequency axis
        if (Input%iter_ord.eq.2) then

          call omegabuildin(Frec,Red,Atom,Atmo,Bfield%Bstrength, &
                            Input,Geom,MPID,lp,ofram,.False.)

          if (laborted) goto 1000

          ! If storing redistribution
          if (PRAM) then

            ! If CPU went above ram
            if (ofram.and.rlimw) then

              write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                  ' reached the limit of redistribution allocations.'
              call verbose

              rlimw = .False.

            end if

            call MPI_BARRIER(MPI_COMM_RT, ierr)

          end if ! Storing profiles

          if (gpid.eq.0) then
            umsg = ' - Input frequency axis initialized'
            call verbose
          end if ! Master

#ifdef _OPENMP
          call setomp_2ord(Atom,Frec,Geom,Bfield%Bstrength, &
                           .True.,.False.)
          if (laborted) goto 1000
          call setomp_magn_2ord(Atom,Flgsg,Bfield%Bstrength)
#endif

        end if ! PRD

        ! If iterating
        if (Input%iter_max.ge.Input%iter_min) then

          ! Solve the NLTE problem of the second kind
          if (gpid.eq.0) then
            umsg = ' - Starting iterations'
            call verbose
          end if

          ! If we have to report RAM use
          if (Input%RAMreport) call RAMreport(MPID,Input%folder,1,1)

          ! If measuring performance
          if (Input%g_perf.and.pid.eq.0) &
            call report_time(Input%folder,Input%ID,.True.)

          ! Solver
          call solve(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                     Red,Bfield,Geom,MPID,Input,Flgsg, &
                     JKQ_asym,Stokes,JKQ,JKQS,JKQC)

          ! Control
          if (laborted) goto 1000

          ! If measuring performance
          if (Input%g_perf.and.pid.eq.0) &
            call report_time(Input%folder,Input%ID,.True.)

          ! Store solution
          if (gpid.eq.0) then
            umsg = ' - Saving solution'
            call verbose
          end if

          ! If saving
          if (saving) then

            ! If inversion
            if (run_mode.eq.-1) then

              ! Set solution
              if (SolF%keep_solution) &
                call setsol(SolF,Flgsg,Bfield, &
                            Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                            J00,J00S,J00C,J00P,.False.)

            ! Synthesis
            else

              ! Write to file
              call writesol(Input,'NONE',Frec%omega,Geom,Flgsg, &
                            Bfield,Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)

            end if ! Inversion/synthesis
          end if ! Saving solution

        ! No iterations, but inversion, keeping solution, and allowed
        ! to do so
        else if (run_mode.eq.-1.and.SolF%keep_solution.and. &
                 saving) then

          ! Set solution
          call setsol(SolF,Flgsg,Bfield, &
                      Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                      J00,J00S,J00C,J00P,.False.)

        end if ! There are actual iterations
      endif ! Formal solution for polarization


      !
      ! Stokes emergence
      !

      if (lpe) then

        ! If removing redistribution
        if (.not.(AV.and..not.dyn.and.lp).or.(.not.AV.or.dyn)) &
          MPID%RAM = MPID%RAM - MPID%WRAM

        ! If removing normalization
        if (.not.(.not.dyn.and.lp)) then
          MPID%RAM = MPID%RAM - MPID%VRAM
          MPID%VRAM = 0d0
        end if

        ! Normalize first order
        call normalization(Atom,LTElines,Atmo,Bfield%Bstrength,Geom, &
                           Frec,Input,Flgsg,MPID,rlimw, &
                           lp,.True.,.True.)

        ! Control
        if (laborted) goto 1000

        !
        ! If doing PRD, determine the input frequencies and weights
        !
        if (Input%iter_ord.eq.2) then

          call omegabuildin(Frec,Red,Atom,Atmo,Bfield%Bstrength, &
                            Input,Geom,MPID,lp,ofram,.True.)

          if (laborted) goto 1000

          ! If storing profiles
          if (PRAM) then

            ! If CPU went above ram
            if (ofram.and.rlimw) then

              write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                  ' reached the limit of redistribution allocations.'
              call verbose

            end if

            call MPI_BARRIER(MPI_COMM_RT, ierr)

          end if ! Storing profiles

          if (gpid.eq.0) then
            umsg = ' - Input frequency axis initialized'
            call verbose
          end if ! Master
#ifdef _OPENMP
          call setomp_2ord(Atom,Frec,Geom,Bfield%Bstrength, &
                           .True.,.True.)
          if (laborted) goto 1000
          call setomp_magn_2ord(Atom,Flgsg,Bfield%Bstrength)
#endif
        end if ! PRD

        ! Calculate emergent solutions
        if (gpid.eq.0) then
          umsg = ' - Emergent Stokes'
          call verbose
        end if

        ! If we have to report RAM use
        if (Input%RAMreport) call RAMreport(MPID,Input%folder,1,0)

        ! If measuring performance
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

        ! Solve
        call emergent(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                      Geom,MPID,Input,Flgsg,JKQ_asym,Stokes, &
                      JKQ,JKQC,SolF)

        ! Control
        if (laborted) goto 1000

        ! If measuring performance
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

      end if ! Energence

      ! Control errors
      call control

      !
      ! Clean memory
      !

      ! Clean Frec and Red structures
1000  call cleanFrecandRed(Frec,Red,MPID)
      call free_local_geom(Geom)
      MPID%RAM = MPID%RAM - MPID%WRAM
      MPID%WRAM = 0d0

      return

      end subroutine hanle_polarization

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare solution structure for the inversion problem\n
      !!   SolF(Solution_F_class): Class with the full RT solution\n
      !!       Input(Input_class): Structure with settings data\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!    GeomI(Geometry_class): Structure with the geometry data
      !!                           for the intensity problem\n
      !!     Geom(Geometry_class): Structure with the geometry
      !!                           data\n
      !!              RF(logical): If calculating responses
      subroutine prepare_buffers(SolF,Input,Atom,GeomI,Geom,RF)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Geometry_class):: GeomI, Geom
      type(Input_class):: Input
      type(Solution_F_class):: SolF
      logical, intent(in):: RF

      ! Local

      integer:: ia


      !
      ! If response functions
      !
      if (RF) then

        ! Read
        if (Input%popuinit) then
          Input%mode = 'B'
        else
          Input%mode = 'W'
        end if
        SolF%keep_solution = .False.

      ! No response
      else

        ! Write
        Input%mode = 'W'
        SolF%keep_solution = .True.

      end if

      ! Only master continues
      if (pid.gt.0) return

      ! If not initialized the space for solutions
      if (SolF%no_initialized) then

        ! Initialize intensity?
        if (Input%Type_Inversion.eq.0.or. &
            (Input%Type_Inversion.eq.3.and.Input%force.eq.'I').or. &
            (Input%Type_Inversion.eq.4.and.Input%force.eq.'I')) then

          ! Allocate quantities
          if (KSTK) then
            allocate(SolF%i_StkI(nfreq,GeomI%nPh,GeomI%nTh,nz))
            allocate(SolF%i_StkI_b(nfreq,GeomI%nPh,GeomI%nTh,nz))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_StkI_t(nfreq,GeomI%nPh,GeomI%nTh,nz))
          end if
          allocate(SolF%i_J00(nxt,nz))
          allocate(SolF%i_J00C(nfreq,nz))
          allocate(SolF%i_J00P(nxphot,2,nz))
          allocate(SolF%i_J00_b(nxt,nz))
          allocate(SolF%i_J00C_b(nfreq,nz))
          allocate(SolF%i_J00P_b(nxphot,2,nz))
          if (Input%LM_Method.eq.1) then
            allocate(SolF%i_J00_t(nxt,nz))
            allocate(SolF%i_J00C_t(nfreq,nz))
            allocate(SolF%i_J00P_t(nxphot,2,nz))
          end if

          ! Allocate rhoes
          allocate(SolF%i_rhoes(na))
          allocate(SolF%i_rhoes_b(na))
          if (Input%LM_Method.eq.1) &
            allocate(SolF%i_rhoes_t(na))

          ! For each atom
          do ia=1,nA
            allocate(SolF%i_rhoes(ia)%rho(Atom(ia)%nlevel,nz))
            allocate(SolF%i_rhoes_b(ia)%rho(Atom(ia)%nlevel,nz))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_rhoes_t(ia)%rho(Atom(ia)%nlevel,nz))
          end do

        end if ! Intensity

        ! Initialize polarization
        if (Input%Type_inversion.eq.1.or. &
            Input%Type_inversion.eq.2.or. &
            (Input%Type_Inversion.eq.3.and.Input%force.eq.'N').or. &
            (Input%Type_Inversion.eq.4.and.Input%force.eq.'N')) then

          ! Allocate quantities
          if (KSTK) then
            allocate(SolF%i_Stk(0:3,nfreq,Geom%nPh,Geom%nTh,nz))
            allocate(SolF%i_Stk_b(0:3,nfreq,Geom%nPh,Geom%nTh,nz))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_Stk_t(0:3,nfreq,Geom%nPh,Geom%nTh,nz))
          end if
          allocate(SolF%i_JKQ(-2:2,0:2,nxtran,nz))
          allocate(SolF%i_JKQS(-2:2,0:2,nxtran,nz))
          allocate(SolF%i_JKQC(-2:2,0:2,nfreq,nz))
          allocate(SolF%i_JKQ_b(-2:2,0:2,nxtran,nz))
          allocate(SolF%i_JKQS_b(-2:2,0:2,nxtran,nz))
          allocate(SolF%i_JKQC_b(-2:2,0:2,nfreq,nz))
          if (Input%LM_Method.eq.1) then
            allocate(SolF%i_JKQ_t(-2:2,0:2,nxtran,nz))
            allocate(SolF%i_JKQS_t(-2:2,0:2,nxtran,nz))
            allocate(SolF%i_JKQC_t(-2:2,0:2,nfreq,nz))
          end if

          ! Allocate rhoes
          allocate(SolF%i_rhoes(na))
          allocate(SolF%i_rhoes_b(na))
          if (Input%LM_Method.eq.1) &
            allocate(SolF%i_rhoes_t(na))

          ! For each atom
          do ia=1,nA
            allocate(SolF%i_rhoes(ia)%crho(Atom(ia)%ndim,nz))
            allocate(SolF%i_rhoes_b(ia)%crho(Atom(ia)%ndim,nz))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_rhoes_t(ia)%crho(Atom(ia)%ndim,nz))
          end do

        end if ! Polarization

        ! Flag initialized
        SolF%no_initialized = .False.

      end if ! Initialize

      !
      ! Free space for outputs
      !
      if (allocated(SolF%e_Stk)) deallocate(SolF%e_Stk)
      if (allocated(SolF%e_Ctr)) deallocate(SolF%e_Ctr)
      if (allocated(SolF%e_tau1)) deallocate(SolF%e_tau1)

      end subroutine prepare_buffers

!#####################################################################
!#####################################################################
!#####################################################################

      end module hanle_mod
