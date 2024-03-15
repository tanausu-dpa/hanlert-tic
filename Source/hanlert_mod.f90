      !> HanleRT manager
      module hanlert_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Hao Li (IAC)
!     Roberto Casini (HAO)
!  Contributors:
!     Ricky Egeland (HAO)
!  Start:
!     06/22/2022
!  Last version:
!     03/15/2024 V3.0.25
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/15/2024:   V3.0.25 - Actually make the problem static if
!                             indicated for 1D synthesis (TdPA)
!
!     02/23/2024:   V3.0.24 - Added Input%fvmicro argument to rAtmo
!                             call (TdPA)
!
!     02/19/2024:   V3.0.23 - Added calls to new get_lims routine
!                             in 1.5D and CLE modes (TdPA)
!
!     02/14/2024:   V3.0.22 - In 1.5D mode, added option to skip
!                             specific pixels (TdPA)
!
!     01/09/2024:   V3.0.21 - In 1.5D mode, when the cache says that
!                             the pixel is solved, it is skipped
!                             regardless of mode (TdPA)
!
!     11/14/2023:   V3.0.20 - Bugfix: Initialize 0 maxB in CLE (TdPA)
!
!     10/17/2023:   V3.0.19 - Damping arrays were not being fred in
!                             TIC (TdPA)
!
!     10/16/2023:   V3.0.18 - Added argument to the call to the
!                             hanle subroutine (TdPA)
!                           - Ensure that maxB is initialized (TdPA)
!                           - Do not recycle the buffer_atmo variable
!                             in CLE as a dummy for omegabuild (TdPA)
!
!     09/08/2023:   V3.0.17 - Changed a return for a jump in case of
!                             error to ensure memory cleansing (TdPA)
!
!     08/30/2023:   V3.0.16 - Bugfix: The writing of the cache file
!                             was inconsistent when using a solution
!                             box (TdPA; pointed out by Hao in the
!                             tic_mod equivalent)
!
!     08/24/2023:   V3.0.15 - Added argument to omegabuild calls for
!                             compatibility reasons. All calls in this
!                             module have that argument as a dummy
!                             variable, just an allocatable double
!                             precision 1D array (TdPA)
!
!     08/17/2023:   V3.0.14 - Removed call to control (TdPA)
!
!     08/07/2023:   V3.0.13 - Added arguments for LTE lines (TdPA)
!
!     07/03/2023:   V3.0.12 - Added PRD as a condition to store
!                             the full array of Stokes par. (TdPA)
!                           - Added a call to prepare_syn in
!                             the synthesis in order to set-up
!                             the model atmosphere (TdPA)
!                           - Added a call to updateatmo at the
!                             exit of hanle (TdPA)
!                           - Added calls to free memory dependent
!                             on the pixel after hanle (TdPA)
!                           - Significantly rewriting of HanleRTTIC
!                             to better integrate the inversion
!                             module (TdPA)
!                           - Undid the removal of projection, as it
!                             can have its uses (TdPA)
!                           - Implemented a different version of the
!                             FWHM for multi-wavelength ranges (TdPA)
!                           - Bugfix: The serial 15DS was not reading
!                             the Barklem and kurucz files (TdPA)
!
!     06/13/2023:   V3.0.11 - Removed projection (HL)
!                           - Update the FWHM for multi-wavelength
!                             ranges (HL)
!
!     04/11/2023:   V3.0.10 - Update the scales for multi-wavelength
!                             ranges (HL)
!
!     03/21/2023:    V3.0.9 - Added argument to rBField for
!                             compatibility for a change due to the
!                             inversion mode (TdPA)
!                           - Removed some commented blocks remaining
!                             from the original TIC (TdPA)
!                           - The 1.5D synthesis can now work in
!                             serial mode, although I do not know why
!                             anybody would (TdPA)
!
!     03/08/2023:    V3.0.8 - Added HanleRTTIC to take manage the
!                             synthesis in the inversion module (TdPA)
!
!     02/14/2023:    V3.0.7 - The geometry is now defined in each
!                             code branch (TdPA)
!                           - Added GeomI exclusive for use in the
!                             intensity problem (TdPA)
!
!     11/24/2022:    V3.0.6 - Changes to make HanleCLE work, mainly
!                             adding calls to routines that were
!                             missing before and changing other to
!                             new ones specific for CLE (TdPA)
!
!     10/25/2022:    V3.0.5 - Added allocations in synthesis branches
!                             to adjust for new Hanle arguments (TdPA)
!                           - Changed call to setmpi_sizes to adjust
!                             for its changes (TdPA)
!                           - Changed call to hanle to adjust for its
!                             changes (TdPA)
!                           - Changed call to open_atm_and_cache and
!                             get_column to adjust for their changes
!                             (TdPA)
!                           - Changed call to set_io_buffers to
!                             adjust for its changes (TdPA)
!                           - Limited the call to rBarklem to not
!                             manager processes (TdPA)
!                           - Added reading and managing of ad-hoc
!                             JKQ to the 1.5D synthesis case (TdPA)
!                           - Implemented the option to limit the
!                             FoV to solve in 1.5 synthesis (TdPA)
!                           - Added HanleCLE routine for the CLE
!                             case (TdPA)
!
!     07/27/2022:    V3.0.4 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!
!     07/18/2022:    V3.0.3 - Bugfix: Added definition of maxV in
!                             HanleRT15DS (TdPA)
!                           - Added calls to adjust_IW (TdPA)
!
!     07/13/2022:    V3.0.2 - Added comment header to HanleRT1DS and
!                             HanleRT15DS (TdPA)
!                           - Added calls to rParfunAbund and
!                             rBarklem, which implies a dependence on
!                             the rbarklem_mod and rpfa_mod (TdPA)
!                           - The call to kurucz_get does not require
!                             the Input%resource argument, but Atmo is
!                             now required instead (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: There are quantities in Atom
!                             that require initialization and
!                             depend on the number of processes for
!                             RT. Added a call to prepareatomMPI
!                             for this (TdPA)
!                           - Bugfix: There are some quantities that
!                             the Master has to allocate that depend
!                             on the number of processes for RT. Added
!                             a omegainitmaster call for this (TdPA)
!                           - Bugfix: Fixed wrong tag when receiving
!                             and sending from/to RT groups (TdPA)
!                           - Bugfix: The cache index order was
!                             inverted (TdPA)
!                           - Bugfix: Fixed an issue where the global
!                             Master ignored the last working group
!                             because of a copied .ge. that now must
!                             be .gt. (TdPA)
!                           - Bugfix: Initialize aborting flag for
!                             the slaves (TdPA)
!                           - Bugfix: Fixed an error in the RT groups
!                             receiving coordinates and columns of the
!                             atmospheric model (TdPA)
!                           - Changed a condition on an output
!                             message (TdPA)
!                           - Reordered some commands (TdPA)
!                           - Added one informative message (TdPA)
!
!     06/29/2022:    V3.0.0 - First version. Taking from hanlert
!                             V2.0.1 and hanlert-15d V1.0.2 (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    Hanlert15D TODO: allow for inversion outputs as input model
!                     atmosphere
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!    Manages the model dependent initializations that used to be
!    in hanlert.
!
!    HanleRT1DS:
!      1D synthesis
!
!    HanleRT15DS:
!      1.5D synthesis
!
!    HanleRTTIC:
!      Handles the synthesis for the inversion module
!
!    HanleCLE:
!      Coronal Line Emission (NO OMP!!)
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use cle_mod
      use commons_mod
      use hanle_mod
      use initmodel_mod
      use inspect_mod
      use io_mod
#ifdef DEBUGATMO
      use iosolution_mod
#endif
      use kurucz_mod
      use normalizer_mod
      use omegabuild_mod
      use parameters_mod , only : c, pi
      use chianti_mod
      use psf_mod
      use ratmo_mod
      use ratom_mod
      use rbarklem_mod
      use rbfield_mod
      use rpfa_mod
      use setmpi_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Handles the 1D synthesis case\n
      !!    Input(Input_class): Structure with settings data\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!     Atomb(Atom_class): Structure with the atomic data for
      !!                        background opacities\n
      !!        Mol(Mol_class): Structure with the molecule data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine HanleRT1DS(Input,Atom,Atomb,Mol,Flgsg, &
                            fudge,MPID)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(MPI_class):: MPID

      ! Local
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI,Geom
      type(kurucz_class):: kurucz
      type(Solution_F_class):: dummy2
      type(LTEline_class), dimension(:), allocatable:: LTElines

      logical:: lDK,l1,l2
      logical:: lio,lie,lp,lpe,lporlpe,lload
      integer:: ia
      double precision:: maxB,DwTa
      double precision, dimension(:), allocatable:: dummy


      !
      ! Define communicator
      !
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, gpid, &
                          MPI_COMM_RT, ierr)
      call MPI_COMM_RANK(MPI_COMM_RT,pid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_RT,nproc,ierr)


      !
      ! Define the flags to route the code
      !
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'
      lio = (Input%mode.eq.'W'.and.Input%force.ne.'P').or. &
            (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'I'.or.Input%force.eq.'A'))
      lie = Input%force.eq.'I'.and.Input%nThLOS.gt.0
      lp = (Input%mode.eq.'W'.and.Input%force.ne.'I').or. &
           (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'N'.or.Input%force.eq.'P'.or. &
                      Input%force.eq.'A'))
      lpe = Input%force.ne.'I'.and.Input%nThLOS.gt.0
      lporlpe = lp.or.lpe

      !
      ! Allocate atomic MPI arrays
      !
      call prepareatomMPI(Atom)

      !
      ! Read atmospheric data
      !
      call rAtmo(Input%atmo,Input%source,Input%ID,Atmo,Input%fvmicro)
      nz = Atmo%nz

      ! Forcing static?
      if (dyn.and.Input%static) then

        ! Make the problem static
        dyn = .False.
        Atmo%vx = 0d0
        Atmo%vy = 0d0
        Atmo%vz = 0d0

      end if

      !
      ! Read magnetic field data
      !
      call rBField(Input%bfield,Input%source,Input%ID, &
                   Bfield,Atmo%nz,Input)

      !
      ! Read partition function data and abundances
      !
      call rParfunAbund(Input,Atmo)

      !
      ! Read barklem data
      !
      call rBarklem(Input,Atom,Atomb)

      !
      ! Get minimum and maximum temperatures, and maximum velocity
      !
      Input%minT = minval(Atmo%T)
      Input%maxT = maxval(Atmo%T)
      if (dyn) then
        Input%maxV = maxval(Atmo%vx*Atmo%vx + &
                            Atmo%vy*Atmo%vy + &
                            Atmo%vz*Atmo%vz)
        Input%maxV = sqrt(Input%maxV)
        Input%static = .False.
      else
        Input%maxV = 0d0
        Input%static = .True.
      end if

      !
      ! Check if unmagnetized
      !
      maxB = maxval(Bfield%Bstrength)
      if (maxB.le.0d0) then
        Input%unmagnetized = .True.
      else
        Input%unmagnetized = .False.
      end if

      !
      ! Decide if need to keep Stokes
      !
      if (pid.eq.0) then
        KSTK = KSTK.or.(PRD.and.(dyn.or..not.AV))
      else
        KSTK = PRD.and.(dyn.or..not.AV)
      end if

      !
      ! Set angular quadrature
      !
      call gauss(Input,GeomI,Geom,1,lp.or.lpe,lpe.or.lie,Flgsg)
      if(pid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbose
      end if

      !
      ! Check that the axial symmetry is consistent
      !
      if (axial) &
        call check_axial(Atmo%vx,Atmo%vy,Bfield%Btheta)


      !
      ! Define the output frequency axis
      !
      call omegabuild(Frec,Atom,Input,maxB,lporlpe,dummy)
      call omegainitmaster(Atom)
      if (pid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbose
      end if

      !
      ! Improve the weight determination?
      !
      if (nproc.gt.1) call adjust_IW(Input,Frec%IW_freq)

      !
      ! Organize the tasks splitting
      !
      call setmpi(MPID,Input,Frec%IW_freq)
      if (pid.eq.0.and.MPID%mpi) then
        umsg = ' - Tasks distributed'
        call verbose
      end if

      !
      ! Master remove fudge
      !
      if (pid.eq.0.and.MPID%mpi.and.allocated(fudge%fudge_v)) &
        deallocate(fudge%fudge_v)

      !
      ! Get Kurucz lines
      !
      if (pid.gt.0.or..not.MPID%mpi) then

        ! Read Kurucz data
        if (Input%NK.ge.1) then

          ! Compute Doppler precursor
          if(Input%dws.eq.'MAX')then
            lDK = .True.
            DwTa = sqrt(maxval(Atmo%T))
          else if(Input%dws.eq.'MIN')then
            lDK = .True.
            DwTa = sqrt(minval(Atmo%T))
          else if(Input%dws.eq.'NUM')then
            lDK = .False.
            DwTa = Input%dw*1d-9/c
          end if

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,lDK,kurucz)

        ! If there are not
        else

          kurucz%ntran = 0

        end if
      end if

      ! Control
      call gcontrol

      !
      ! Photoionization quantites
      !
      do ia=1,nA
        call setphoto(Atom(ia),Frec%omega,MPID)
      end do
      if (pid.eq.0) then
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbose
      end if

      !
      ! Resize some frequency quantities if doing MPI
      !
      call frecresize(Frec,Atom,Input,MPID)

      !
      ! Compute size for MPI messages in solvers
      !
      l1 = .not.lload
      l2 = (lio.and..not.lie).or.(lload.and.(lp.or.lpe))
      if (MPID%mpi) &
        call setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,l1,l2,.False.)

      ! Deallocate things not needed
      if (allocated(Input%atomback)) deallocate(Input%atomback)
      if (allocated(Input%mol)) deallocate(Input%mol)

#ifdef DEBUGATMO
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif
      ! Prepare for synthesis
      call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

      ! If error, skip
      if (laborted) return

      !
      ! Solution of the NLTE problem of the second kind
      !
      call hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input,GeomI, &
                 Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                 dummy,dummy2,lload,lio,lie,lp,lpe,.True.)

      ! If error, skip
      if (laborted) return

      ! Verbose
      if(pid.eq.0) then
        umsg = ' - Solver finished'
        call verbose
      end if

      ! Recalculate electrons
      call updateatmo(Atom,Atomb,Atmo,Bfield,Input)

      ! Clean memory
      call free_pix(Atom,Atomb,Mol,Bfield)
      call free_Atmo(Atmo,.True.)

      end subroutine HanleRT1DS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Handles the 1.5D synthesis case\n
      !!    Input(Input_class): Structure with settings data\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!     Atomb(Atom_class): Structure with the atomic data for
      !!                        background opacities\n
      !!        Mol(Mol_class): Structure with the molecule data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine HanleRT15DS(Input,Atom,Atomb,Mol, &
                             Flgsg,fudge,MPID)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(MPI_class):: MPID

      ! Local
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI, Geom
      type(kurucz_class):: kurucz
      type(Solution_F_class):: dummy
      type(LTEline_class), dimension(:), allocatable:: LTElines

      logical:: aborting,check,receiving,lcache
      logical:: lio,lie,lp,lpe,lporlpe,lload,l1,l2
      logical:: double,lexcl

      logical, dimension(:,:), allocatable:: cache

      integer:: unitA,unitC,unitJ
      integer:: ix1,iy1,sizeA,sizeJ
      integer:: ia,ix0,iy0,ix,iy,inod,NLOSr,NLOS,ip,iproc

      integer, dimension(3):: dims,out_dims,int_buff
      integer, dimension(:), allocatable:: cpu_free

      double precision:: maxB,DwTa

      double precision, dimension(:), allocatable,target:: buffer_atmo
      double precision, dimension(:), allocatable,target:: buffer_JKQ

      ! Deallocate things not needed
      if (allocated(Input%atomback)) deallocate(Input%atomback)
      if (allocated(Input%mol)) deallocate(Input%mol)

      ! Trick maxB if not unmagnetized by force
      if (Input%unmagnetized) then
        maxB = 0d0
      else
        maxB = 1d0
      end if

      ! Check known atmospheric limits
      if ((PRD.and.(Input%minT.lt.0d0.or.Input%maxT.lt.0d0)).or. &
          (Input%minT.lt.0d0.and.Input%dws.eq.'MIN').or. &
          (Input%maxT.lt.0d0.and.Input%dws.eq.'MAX').or. &
          (.not.Input%static.and.Input%maxV.lt.0d0)) then

        ! Read whole model to get limits
        call get_lims(Input,1,aborting)

        ! Check if could read
        laborted = aborting

        ! Check if aborting
        call gcontrol

      end if

      ! Correct maxV
      if (Input%static) then
        Input%maxV = 0d0
      else
        Input%maxV = Input%maxV*1d-6/c
      end if

      !
      ! Define the flags to route the code
      !
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'
      lio = (Input%mode.eq.'W'.and.Input%force.ne.'P').or. &
            (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'I'.or.Input%force.eq.'A'))
      lie = Input%force.eq.'I'.and.Input%nThLOS.gt.0
      lp = (Input%mode.eq.'W'.and.Input%force.ne.'I').or. &
           (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'N'.or.Input%force.eq.'P'.or. &
                      Input%force.eq.'A'))
      lpe = Input%force.ne.'I'.and.Input%nThLOS.gt.0
      lporlpe = lp.or.lpe

      !
      ! Define the output frequency axis
      !
      call omegabuild(Frec,Atom,Input,maxB,lporlpe,buffer_atmo)
      if (gpid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbose
      end if

      !
      ! Split in groups of tasks
      !
      call setmpi15D(MPID,Input)
      if (gpid.eq.0) then
        umsg = ' - Tasks distributed'
        call verbose
      end if

      !
      ! Decide if need to keep Stokes
      !
      if (.not.MPID%mpi15d.or.(pid.eq.0.and.gpid.ne.0)) then
        KSTK = KSTK.or.(PRD.and.(.not.Input%static.or..not.AV))
      else
        KSTK = PRD.and.(.not.Input%static.or..not.AV)
      end if

      !
      ! Set angular quadrature
      !
      call gauss(Input,GeomI,Geom,1,lp.or.lpe,lpe.or.lie,Flgsg)
      if(gpid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbose
      end if

      !
      ! Masters remove fudge
      !
      if (gpid.eq.0.and.allocated(fudge%fudge_v).and.MPID%mpi15d) &
        deallocate(fudge%fudge_v)


      !
      ! Pre-process
      !

      ! Initialize looping identifier
      ix1 = -1
      iy1 = -1

      ! Master
      if (gpid.eq.0) then

        ! Open files (ia and nz are a dummy variable here)
        call open_atm_and_cache(Input,1,unitA,unitC,aborting,dims, &
                                ia,double,nz,cache,lcache)

        ! The master does not need the background atoms or
        ! molecules
        if (MPID%mpi15d) deallocate(Atomb,Mol)

        ! Check if could read
        laborted = aborting

        ! Open JKQ asymmetry file
        if (Input%nasym.gt.0.and..not.laborted) then

          ! Open asymmetry file
          call open_asymm(Input,unitJ,aborting,dims)

          ! Check if could read
          laborted = aborting

        end if
      end if

      ! Check if aborting
      call gcontrol

      ! Share dims
      call MPI_BCAST(dims(1),3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ! Translate into common height variable
      nz = dims(3)

      !
      ! Now fix the real dimensions from solution box
      !

      ! Fix wildcards in input
      if (Input%sol_box(1).lt.1) Input%sol_box(1) = 1
      if (Input%sol_box(2).lt.1) Input%sol_box(2) = dims(1)
      if (Input%sol_box(3).lt.1) Input%sol_box(3) = 1
      if (Input%sol_box(4).lt.1) Input%sol_box(4) = dims(2)

      ! Get new dimensions
      if (gpid.eq.0) then

        ! Get output dimensions
        out_dims(1) = Input%sol_box(2) - Input%sol_box(1) + 1
        out_dims(2) = Input%sol_box(4) - Input%sol_box(3) + 1
        out_dims(3) = dims(3)

      else

        dims(1) = Input%sol_box(2) - Input%sol_box(1) + 1
        dims(2) = Input%sol_box(4) - Input%sol_box(3) + 1

      end if


      !
      ! Organize the tasks splitting (within groups)
      !
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! Improve the weight determination?
        if (nproc.gt.1) call adjust_IW(Input,Frec%IW_freq)

        ! Distribute
        call setmpi(MPID,Input,Frec%IW_freq)

        !
        ! Masters remove fudge
        !
        if ((pid.eq.0.and.MPID%mpi).and.allocated(fudge%fudge_v)) &
          deallocate(fudge%fudge_v)

        !
        ! Finish the atom
        !
        call prepareatomMPI(Atom)
        call omegainitmaster(Atom)

      end if

      ! Verbose
      if (gpid.eq.0) then

        umsg = ' - Tasks within groups distributed'
        call verbose

      end if

      ! Control
      call gcontrol


      !
      ! Photoionization quantites
      !
      if (.not.MPID%mpi15d.or.gpid.gt.0) then
        do ia=1,nA
          call setphoto(Atom(ia),Frec%omega,MPID)
        end do
      else
        call gcontrol
      end if

      ! Verbose
      if (gpid.eq.0) then
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbose
      end if

      ! Not global Master
      if (.not.MPID%mpi15d.or.gpid.gt.0)  then

        !
        ! Resize some frequency quantities if doing MPI
        !
        call frecresize(Frec,Atom,Input,MPID)

        !
        ! Compute size for MPI messages in solvers
        !
        l1 = .not.lload
        l2 = (lio.and..not.lie).or.(lload.and.(lp.or.lpe))
        if (MPID%mpi) &
          call setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,l1,l2,.False.)

      end if

      ! Sanity check the limits for the outputs
      call check_io_buffers_sanity_check(Input,Atom)

      ! Check if aborting
      call gcontrol

      ! Prepare buffer atmosphere
      sizeA = dims(3)*24
      allocate(buffer_atmo(sizeA))

      ! If asymmetry data
      if (Input%nasym.gt.0) then
        sizeJ = dims(3)*8
        allocate(buffer_JKQ(sizeJ))
      else
        sizeJ = 0
      end if

      ! Initialize buffer sizes
      call set_io_buffers(Input,0,Atom,Frec)

      ! Check if aborting
      call gcontrol

      !
      ! Check files exist
      !
      if (gpid.eq.0) then

        ! Check files exist
        if (lcache) call check_io_buffers_exist(Input,Atom, &
                                                Geom,aborting)
        laborted = aborting

      end if

      ! Check if aborting
      call gcontrol

      !
      ! Read partition function data and abundances
      !
      call rParfunAbund(Input,Atmo)

      ! Check if aborting
      call gcontrol

      !
      ! Read barklem data
      !
      if (.not.MPID%mpi15d.or.gpid.gt.0) &
        call rBarklem(Input,Atom,Atomb)

      ! Check if aborting
      call gcontrol

      !
      ! Get Kurucz lines
      !
      if ((.not.MPID%mpi15d.and.(pid.gt.0.or..not.MPID%mpi)).or. &
          (gpid.gt.0.and.(pid.gt.0.or..not.MPID%mpi))) then

        ! Read Kurucz data
        if (Input%NK.ge.1) then

          DwTa = Input%dw*1d-9/c

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,.False.,kurucz)

        ! If there are
        else

          kurucz%ntran = 0

        end if
      end if

      ! Check if aborting
      call gcontrol

      !
      ! Properly start now
      !
      if (gpid.eq.0) then
        umsg = ' - Starting to work on your model, this may '// &
               'take a while'
        call verbose
      end if

      !
      ! MPI version
      !
      if (MPID%mpi15d) then

        !
        ! Master
        !
        if (gpid.eq.0) then

          ! Open files so slaves can write later
          if (.not.lcache) call create_io_files(Input,Atom,out_dims, &
                                                Geom,GeomI,Frec)

          ! Allocate cpu_free with group status
          allocate(cpu_free(MPID%ngroup))
          cpu_free = 1

          ! Initialize indexes and sizes
          ix = 1
          iy = 0
          inod = 0
          NLOSr = 0
         !NLOS = dims(1)*dims(2)
          NLOS = out_dims(1)*out_dims(2)

          ! Work until exhausted
          do while (.True.)

            ! If aborting, send stop signal to everyone
            if (aborting) then

              ! For every CPU still sending
              do while (minval(cpu_free).lt..5d0)

                ! For every leader
                do ip=1,MPID%ngroup

                  ! If already free, skip
                  if (cpu_free(ip).gt..5d0) cycle

                  ! Test if slave in the group is sending something
                  call MPI_IPROBE(MPID%ltslave(ip), &
                                  4000000+MPID%ltslave(ip), &
                                  MPI_COMM_WORLD, receiving, &
                                  MPI_STATUS_IGNORE, ierr)

                  ! If slave is calling
                  if (receiving) then

                    ! Receive the ping
                    call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                  MPID%ltslave(ip), &
                                  4000000+MPID%ltslave(ip), &
                                  MPI_COMM_WORLD, &
                                  MPI_STATUS_IGNORE, ierr)

                    ! Free the group
                    cpu_free(ip) = 1

                  end if ! Receiving from a CPU

                end do ! Receive from everyone
              end do ! While there is someone working

            end if ! Aborting

            ! If there are LOS to do and at least one free CPU
            if (inod.lt.NLOS.and.maxval(cpu_free).gt..5d0) then

              ! Advance one node
              ix0 = ix
              iy0 = iy
              iy = iy + 1
              if (iy.gt.dims(2)) then
                ix = ix + 1
                iy = 1
                if (ix.gt.dims(1)) cycle
              end if

              ! Check if not looping
              if (ix.ne.ix1.or.iy.ne.iy1) then

                ! Get column from atmosphere
                call get_column(unitA, buffer_atmo, double, check)

                ! If JKQ data
                if (sizeJ.gt.0) &
                  call get_column(unitJ, buffer_JKQ, .True., check)

                ! Store last read
                ix1 = ix
                iy1 = iy

              end if ! Loopìng

              ! Check if out of box
              if (ix.lt.Input%sol_box(1).or. &
                  ix.gt.Input%sol_box(2).or. &
                  iy.lt.Input%sol_box(3).or. &
                  iy.gt.Input%sol_box(4)) then
                cycle
              end if

              ! Now advance LOS and set slave indexes
              inod = inod + 1
              icoords = (/ ix - Input%sol_box(1) + 1 , &
                           iy - Input%sol_box(3) + 1 , inod /)

              ! If done in cache, skip
              if (lcache) then
                if (cache(iy,ix)) then
                  NLOSr = NLOSr + 1
                  cycle
                end if
              end if

              ! If excluding
              if (Input%lexcl) then
                if (ix.ge.Input%excl(1,1).and. &
                    ix.le.Input%excl(1,Input%nexcl)) then
                  lexcl = .False.
                  do ia=1,Input%nexcl
                    if (Input%excl(1,ia).lt.ix) cycle
                    if (Input%excl(1,ia).gt.ix) exit
                    if (Input%excl(2,ia).eq.iy) then
                      lexcl = .True.
                      exit
                    end if
                  end do
                  if (lexcl) then
                    NLOSr = NLOSr + 1
                    cycle
                  end if
                end if
              end if

              ! Take a free cpu
              ip = maxloc(cpu_free, 1)

              ! Send signal to node
              call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                            MPID%ltslave(ip), &
                            2000000+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, ierr)

              ! If failed
              if (ierr.ne.0) then
                inod = inod - 1
                ix = ix0
                iy = iy0
                cycle
              end if

              ! Send data to node
              call MPI_SEND(buffer_atmo(1), sizeA, &
                            MPI_DOUBLE_PRECISION, &
                            MPID%ltslave(ip), &
                            3000000+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, ierr)

              ! If failed
              if (ierr.ne.0) then
                inod = inod - 1
                ix = ix0
                iy = iy0
                cycle
              end if

              ! If JKQ data
              if (sizeJ.gt.0) then

                ! Send JKQ data to node
                call MPI_SEND(buffer_JKQ(1), sizeJ, &
                              MPI_DOUBLE_PRECISION, &
                              MPID%ltslave(ip), &
                              5000000+MPID%ltslave(ip), &
                              MPI_COMM_WORLD, ierr)

                ! If failed
                if (ierr.ne.0) then
                  inod = inod - 1
                  ix = ix0
                  iy = iy0
                  cycle
                end if

              end if

              ! That CPU is now busy
              cpu_free(ip) = 0

            end if ! If there is work to do and free CPUs

            ! For every slave group
            ip = 1
            do while (.True.)

              ! Test if slave in the group is sending something
              call MPI_IPROBE(MPID%ltslave(ip), &
                              4000000+MPID%ltslave(ip), &
                              MPI_COMM_WORLD, receiving, &
                              MPI_STATUS_IGNORE, ierr)

              ! If failing
              if (ierr.ne.0) cycle

              ! If slave is calling
              if (receiving) then

                ! Try to receive
                do while (.True.)

                  ! Receive the ping
                  call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                MPID%ltslave(ip), &
                                4000000+MPID%ltslave(ip), &
                                MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! If failed, try again
                  if (ierr.ne.0) cycle

                  exit

                end do

                ! Convert ix coordinate into node coordinate
                int_buff(1) = &
                    (int_buff(1) - 2 + Input%sol_box(1))*dims(2) + &
                    int_buff(2) - 1 + Input%sol_box(3)

                ! Write in cache
                call write_cache(unitC,Input%cache,int_buff,check)
                aborting = .not.check

                ! Update NLOS received
                NLOSr = NLOSr + 1

                ! Free the group
                cpu_free(ip) = 1

              end if ! Receiving from a CPU group

              ip = ip + 1

              if (ip.gt.MPID%ngroup) exit

            end do ! Slaves

            ! If we went beyond the number of LOS, exit
            if (NLOSr.ge.NLOS) exit

          end do ! While there is work to do

          ! If we are done, notify to slaves
          icoords(1) = -1
          iproc = 1
          do while (.True.)

            ! send termination signal
            call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                          MPID%ltslave(iproc), &
                          2000000+MPID%ltslave(iproc), &
                          MPI_COMM_WORLD, ierr)

            ! If it fails
            if (ierr.ne.0) cycle

            iproc = iproc + 1

            if (iproc.gt.MPID%ngroup) exit

          end do ! slaves

        !
        ! Slaves
        !
        else

          ! Initialize
          aborting = .False.

          ! Work until further notice
          do while (.True.)

            ! If leader
            if (pid.eq.0) then

              ! Try receiving until success
              do while (.True.)

                ! Wait for signal
                call MPI_RECV(icoords(1), 3, MPI_INTEGER, 0, &
                              2000000+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! Check if it is the termination signal
                if (icoords(1).lt.1) then
                  aborting = .True.
                  exit
                end if

                ! Receive LOS
                call MPI_RECV(buffer_atmo(1), sizeA, &
                              MPI_DOUBLE_PRECISION, &
                              0, 3000000+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! If JKQ data
                if (sizeJ.gt.0) then

                  ! Receive LOS
                  call MPI_RECV(buffer_JKQ(1), sizeJ, &
                                MPI_DOUBLE_PRECISION, &
                                0, 5000000+gpid, MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! If it fails
                  if (ierr.ne.0) cycle

                end if

                ! Success
                exit

              end do

            end if

            ! If liutenant has friends
            if (nproc.gt.1) then

              ! Try until done
              do while (.True.)

                ! Broadcast
                call MPI_BCAST(icoords(1), 3, MPI_INTEGER, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! If aborting
                if (icoords(1).lt.1) then
                  aborting = .True.
                  exit
                end if

                ! Broadcast
                call MPI_BCAST(buffer_atmo(1), sizeA, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! If JKQ data
                if (sizeJ.gt.0) then

                  ! Broadcast
                  call MPI_BCAST(buffer_JKQ(1), sizeJ, &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                  ! If it fails
                  if (ierr.ne.0) cycle

                end if

                ! Success
                exit

              end do ! Try communicating until successfull

            end if ! MPI in RT

            ! Problem
            if (aborting) exit

            ! Reorder buffer
            call rAtmo_frombuffer(buffer_atmo, Input, &
                                  Atmo, Bfield, dims)

            ! Initialize aborting flag
            laborted = .False.

#ifdef DEBUGATMO
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif
            ! Prepare for synthesis
            call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

            ! If error, skip
            if (laborted) goto 1100

            ! Solve columns
            call hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input, &
                       GeomI,Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                       buffer_JKQ,dummy,lload,lio,lie,lp,lpe,.False.)

            ! If error, skip
            if (laborted) goto 1100

            ! Recalculate electrons
            call updateatmo(Atom,Atomb,Atmo,Bfield,Input)

            ! Clean memory
1100        call free_pix(Atom,Atomb,Mol,Bfield)
            call free_Atmo(Atmo,.False.)

            ! If liutenant, send message to grand master
            if (pid.eq.0) then

              ! Fail
              if (laborted) then
                icoords(3) = -1
              ! Success
              else
                icoords(3) = gpid
              end if

              ! Try until achieved
              do while (.True.)
                call MPI_SEND(icoords(1),3,MPI_INTEGER,0, &
                              4000000+gpid,MPI_COMM_WORLD,ierr)
                if (ierr.ne.0) cycle
                exit
              end do


            end if ! Send info back to grand master

          end do

        end if ! Master or slave

      ! Serial version
      else

        ! Open files so slaves can write later
        if (.not.lcache) call create_io_files(Input,Atom,out_dims, &
                                              Geom,GeomI,Frec)

        ! Initialize indexes and sizes
        ix = 1
        iy = 0
        inod = 0
        NLOSr = 0
       !NLOS = dims(1)*dims(2)
        NLOS = out_dims(1)*out_dims(2)

        ! Work until exhausted
        do while (.True.)

          ! Initialize
          aborting = .False.

          ! If aborting, send stop signal to everyone
          if (aborting) then

            ! Abort
            call aborted_silent

          end if ! Aborting

          ! Advance one node
          ix0 = ix
          iy0 = iy
          iy = iy + 1
          if (iy.gt.dims(2)) then
            ix = ix + 1
            iy = 1
            if (ix.gt.dims(1)) cycle
          end if

          ! Check if not looping
          if (ix.ne.ix1.or.iy.ne.iy1) then

            ! Get column from atmosphere
            call get_column(unitA, buffer_atmo, double, check)

            ! If JKQ data
            if (sizeJ.gt.0) &
              call get_column(unitJ, buffer_JKQ, .True., check)

            ! Store last read
            ix1 = ix
            iy1 = iy

          end if ! Loopìng

          ! Check if out of box
          if (ix.lt.Input%sol_box(1).or. &
              ix.gt.Input%sol_box(2).or. &
              iy.lt.Input%sol_box(3).or. &
              iy.gt.Input%sol_box(4)) then
            cycle
          end if

          ! Now advance LOS and set slave indexes
          inod = inod + 1
          icoords = (/ ix - Input%sol_box(1) + 1 , &
                       iy - Input%sol_box(3) + 1 , inod /)

          ! If done in cache, skip
          if (lcache) then
            if (cache(iy,ix)) then
              if (.not.lload) then
                NLOSr = NLOSr + 1
                cycle
              end if
            end if
          end if

          ! Reorder buffer
          call rAtmo_frombuffer(buffer_atmo, Input, &
                                Atmo, Bfield, dims)

          ! Initialize aborting flag
          laborted = .False.

#ifdef DEBUGATMO
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif
          ! Prepare for synthesis
          call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

          ! If error, skip
          if (laborted) goto 1200

          ! Solve columns
          call hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input,GeomI, &
                     Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                     buffer_JKQ,dummy,lload,lio,lie,lp,lpe,.False.)

          ! If error, skip
          if (laborted) goto 1200

          ! Recalculate electrons
          call updateatmo(Atom,Atomb,Atmo,Bfield,Input)

          ! Clean memory
1200      call free_pix(Atom,Atomb,Mol,Bfield)
          call free_Atmo(Atmo,.False.)

          ! Copy coords
          int_buff = icoords

          ! Convert ix coordinate into node coordinate
          int_buff(1) = &
              (int_buff(1) - 2 + Input%sol_box(1))*dims(2) + &
              int_buff(2) - 1 + Input%sol_box(3)

          ! Write in cache
          call write_cache(unitC,Input%cache,int_buff,check)
          aborting = .not.check

          ! Update NLOS received
          NLOSr = NLOSr + 1

          ! If we went beyond the number of LOS, exit
          if (NLOSr.ge.NLOS) exit

        end do ! While there is work to do

      end if ! MPI/serial

      !
      ! Close files
      !
      if (gpid.eq.0) call close_file(unitA)

      end subroutine HanleRT15DS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Handles the synthesis for the inversion module\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!        Atomb(Atom_class): Structure with the atomic data
      !!                           for background opacities\n
      !!           Mol(Mol_class): Structure with the molecule
      !!                           data\n
      !!    GeomI(Geometry_class): Structure with the geometry data
      !!                           for the intensity problem\n
      !!     Geom(Geometry_class): Structure with the geometry
      !!                           data\n
      !!       Flgsg(Fctsg_class): Structure with factorials and
      !!                           signs\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!       fudge(fudge_class): Structure with fudge data\n
      !!     kurucz(kurucz_class): Structure with Kurucz line data\n
      !!          MPID(MPI_class): Structure with MPI data
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Bfield(Bfield_class): Structure with the vertical
      !!                           magnetic field data\n
      !!       Input(Input_class): Structure with settings data\n
      !!      Sol(Solution_class): Class with the data of the RT
      !!                           solution\n
      !!   SolF(Solution_F_class): Class with the full RT solution\n
      !!              RF(logical): If calculating responses
      subroutine HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                            fudge,kurucz,MPID,Atmo,Bfield,Input, &
                            Sol,SolF,RF)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(Geometry_class):: GeomI, Geom
      type(Frequency_class):: Frec
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(MPI_class):: MPID
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Solution_class):: Sol
      type(Solution_F_class):: SolF
      logical, intent(in):: RF

      ! Local
      type(LTEline_class), dimension(:), allocatable:: LTElines

      logical:: PRAM_local, IRAM_local
      logical:: lio,lie,lp,lpe,lporlpe,lload
      logical:: tau1_local, contr_local

      integer:: ia,nfreq_local


      ! Store local copies
      PRAM_Local = PRAM
      IRAM_Local = IRAM
      nfreq_local = nfreq

      ! Be sure that recurrent atmospheric quantities are deallocated
      ! if the type of atmosphere is not zero
      if (allocated(Atmo%zalt)) deallocate(Atmo%zalt)
      if (associated(Atmo%Bx)) nullify(Atmo%Bz,Atmo%Bx,Atmo%By)


#ifdef DEBUGATMO
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif
      !
      ! Prepare the models for synthesis
      call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

      ! If failed, skip
      if (laborted) goto 1100

      !
      ! If projection is on
      if (Sol%Projection) then

        ! At every height
        do ia=1,Atmo%nZ

          ! If magnetic field strength is set negative
          if (Bfield%Bstrength(ia).lt.0) then

            ! Change to correct module but opposite inclination
            Bfield%Bstrength(ia) = -Bfield%Bstrength(ia)
            Bfield%Btheta(ia) = pi

          ! If positive, force vertical
          else

            Bfield%Btheta(ia) = 0d0

          end if ! Sign of Bstrength

        end do ! Heights

      end if ! Projection

      !
      ! Prepare buffers
      call prepare_buffers(SolF,Input,Atom,GeomI,Geom,RF)

      !
      ! Define the flags to route the code
      !
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'
      lio = (Input%mode.eq.'W'.and.Input%force.ne.'P').or. &
            (Input%mode.eq.'B'.and. &
             (Input%force.eq.'I'.or.Input%force.eq.'A'))
      lie = Input%force.eq.'I'.and.Geom%nThLOS.GT.0
      lp = (Input%mode.eq.'W'.and.Input%force.ne.'I').or. &
           (Input%mode.eq.'B'.and. &
            (Input%force.eq.'N'.or.Input%force.eq.'P'.or. &
             Input%force.eq.'A'))
      lpe = Input%force.ne.'I'.and.Geom%nThLOS.gt.0
      lporlpe = lp.or.lpe

      ! If RF, that means we do not need tau and contribution
      if (RF) then

        ! Copy to local
        tau1_local = Input%out_tau1
        contr_local = Input%out_contr

        ! Force to false
        Input%out_tau1 = .False.
        Input%out_contr = .False.

      end if

      ! Call synthesis
      call hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input,GeomI,Geom, &
                 Bfield,Frec,Flgsg,fudge,kurucz,Atmo%JKQin,SolF, &
                 lload,lio,lie,lp,lpe,.False.)

      ! If failed
      if (laborted) goto 1100


      ! Apply PSF
      call Profiles_Out(Frec,Sol,SolF%e_Stk,lpe, &
                        Input%lim_stk,Input%lim_fwhm)


      ! If diffuse light
      if (Sol%Diff_flag) then

        ! If only intensity
        if (Input%force.eq.'I') then

          ! Combine intensity
          Sol%Stokes_out(0,:) = Atmo%f_diff*Sol%Stokes_diff(0,:) + &
                               (1d0 - Atmo%f_diff)*Sol%Stokes_out(0,:)

        ! Polarization
        else

          ! Combine Stokes
          Sol%Stokes_out = Atmo%f_diff*Sol%Stokes_diff + &
                           (1d0 - Atmo%f_diff)*Sol%Stokes_out

        end if ! Polarized
      end if ! Diffuse light

      ! If active projection
      if (Sol%Projection) then

        ! For each height
        do ia=1,Atmo%nZ

          ! If inclination larger than 3rad
          if (Bfield%Btheta(ia).gt.3.0) then

            ! The field was just negative polarity
            Bfield%Bstrength(ia) = -Bfield%Bstrength(ia)
            Bfield%Btheta(ia) = 0d0

          end if

        end do

      end if

      ! If working with fractional Stokes parameters
      if (Sol%Fractional) then

        ! If only intensity
        if (Input%force.ne.'I') then

          ! Fractional polarization
          Sol%Stokes_out(1,:) = 1d2*Sol%Stokes_out(1,:)/ &
                                Sol%Stokes_out(0,:)
          Sol%Stokes_out(2,:) = 1d2*Sol%Stokes_out(2,:)/ &
                                Sol%Stokes_out(0,:)
          Sol%Stokes_out(3,:) = 1d2*Sol%Stokes_out(3,:)/ &
                                Sol%Stokes_out(0,:)

        end if ! Polarization

        ! For each range in wavelengths
        do ia=1,Sol%Num_Range

          ! Scale to the range just the intensity
          Sol%Stokes_out(0,Sol%Range(ia,1):Sol%Range(ia,2)) = &
                  Sol%Stokes_out(0,Sol%Range(ia,1):Sol%Range(ia,2))/ &
                  Sol%Scal_Stokes(ia)
        end do

      ! If non-fractional
      else

        ! If only intensity
        if (Input%force.eq.'I') then

          ! For each range in wavelengths
          do ia=1,Sol%Num_Range

            ! Scale to range just the intensity
            Sol%Stokes_out(0,Sol%Range(ia,1):Sol%Range(ia,2)) = &
              Sol%Stokes_out(0,Sol%Range(ia,1):Sol%Range(ia,2))/ &
              Sol%Scal_Stokes(ia)

          end do

        ! Full Stokes
        else

          ! For each range in wavelengths
          do ia=1,Sol%Num_Range

            ! Scale Stokes parameters in this range
            Sol%Stokes_out(0:3,Sol%Range(ia,1):Sol%Range(ia,2)) = &
              Sol%Stokes_out(0:3,Sol%Range(ia,1):Sol%Range(ia,2))/ &
              Sol%Scal_Stokes(ia)

          end do

        end if ! Intensity/Stokes
      end if ! Fractional or not

1100  nfreq = nfreq_local
      PRAM = PRAM_Local
      IRAM = IRAM_Local

      ! If RF
      if (RF) then

        ! Restore variables
        Input%out_tau1 = tau1_local
        Input%out_contr = contr_local

      end if

      ! Free memory
      call free_mol(Mol)
      call free_damp(Atom,nA)
      if (nAb.gt.0) call free_damp(Atomb,nAb)
      call free_gpop(Atom,Atomb,Mol)
      call free_lpop(Atom,Atomb)
      call free_cols(Atom)

      end subroutine HanleRTTIC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Handles the CLE synthesis case\n
      !!    Input(Input_class): Structure with settings data\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!     Atomb(Atom_class): Structure with the atomic data for
      !!                        background opacities\n
      !!        Mol(Mol_class): Structure with the molecule data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine HanleCLE(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(MPI_class):: MPID

      ! Local
      type(Atmo_class):: Atmo
      type(Frequency_class):: Frec
      type(Geometry_class):: Geom,Gdummy
      type(kurucz_class):: kurucz
      type(chianti_class):: chianti
      type(spect_class):: spect

      logical:: aborting,check,receiving,lcache
      logical:: lie,lpe,lload,double

      logical, dimension(:,:), allocatable:: cache

      integer:: unitA,unitC
      integer:: NIV,NIC
      integer:: iy1,iz1,sizeA,mode,norm,ioffset
      integer:: ia,iy0,iz0,iy,iz,inod,NLOSr,NLOS,ip,iproc

      integer, dimension(3):: dims,int_buff
      integer, dimension(:), allocatable:: cpu_free
      integer, dimension(:), allocatable:: unitI
      integer, dimension(nA):: ion_value_ind
      integer, dimension(nA):: ion_column_ind

      double precision:: maxB,DwTa,ycoor,zcoor
      double precision, dimension(nA):: ion_value

      double precision, dimension(:), allocatable,target:: buffer
      double precision, dimension(:), pointer:: buffer_atmo
      double precision, dimension(:), pointer:: buffer_ion
      double precision, dimension(:), allocatable:: x,y,z,dummy


      ! Deallocate things not needed
      if (allocated(Input%atomback)) deallocate(Input%atomback)
      if (allocated(Input%mol)) deallocate(Input%mol)

      ! Trick maxB if not unmagnetized by force
      if (Input%unmagnetized) then
        maxB = 0d0
      else
        maxB = 1d0
      end if

      ! Check known atmospheric limits
      if ((PRD.and.(Input%minT.lt.0d0.or.Input%maxT.lt.0d0)).or. &
          (Input%minT.lt.0d0.and.Input%dws.eq.'MIN').or. &
          (Input%maxT.lt.0d0.and.Input%dws.eq.'MAX').or. &
          (.not.Input%static.and.Input%maxV.lt.0d0)) then

        ! Read whole model to get limits
        call get_lims(Input,1,aborting)

        ! Check if could read
        laborted = aborting

        ! Check if aborting
        call gcontrol

      end if

      ! Correct maxV
      if (Input%static) then
        Input%maxV = 0d0
      else
        Input%maxV = Input%maxV*1d-6/c
      end if

      !
      ! Define the flags to route the code
      !
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'
      lie = Input%force.eq.'I'
      lpe = Input%force.ne.'I'

      !
      ! Set angular quadrature
      !
      call gauss(Input,Gdummy,Geom,2,.True.,.False.,Flgsg)
      if(pid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbose
      end if

      !
      ! Define the output frequency axis
      !
      call omegabuild(Frec,Atom,Input,maxB,lpe,dummy)
      if (gpid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbose
      end if

      !
      ! Split in groups of tasks
      !
      call setmpi15D(MPID,Input)
      if (gpid.eq.0) then
        umsg = ' - Tasks distributed'
        call verbose
      end if

      !
      ! Not keeping Stokes
      !
      if (pid.eq.0.and.gpid.ne.0) then
        KSTK = KSTK.or.(PRD.and.(.not.Input%static.or..not.AV))
      else
        KSTK = PRD.and.(.not.Input%static.or..not.AV)
      end if

      !
      ! Masters remove fudge
      !
      if (gpid.eq.0.and.allocated(fudge%fudge_v)) &
        deallocate(fudge%fudge_v)


      !
      ! Pre-process
      !

      ! Initialize looping identifier
      iy1 = -1
      iz1 = -1

      ! Master
      if (gpid.eq.0) then

        ! Open files
        call open_atm_and_cache(Input,2,unitA,unitC,aborting,dims, &
                                mode,double,norm,cache,lcache)

        ! The master does not need the background atoms or
        ! molecules
        deallocate(Atomb,Mol)

        ! Check if could read
        laborted = aborting

      end if

      ! Check if aborting
      call gcontrol

      ! Share mode
      call MPI_BCAST(mode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ! Share norm
      call MPI_BCAST(norm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ! Share dims
      call MPI_BCAST(dims(1),3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      !
      ! Deal with ionization fractions (loaded ones)
      !

      ! Initialize
      ion_value = -1d0
      ion_value_ind = -1
      ion_column_ind = -1
      NIC = 0
      NIV = 0

      ! For each atom
      do ia=1,NA

        ! Check if explicitly None
        if (Input%ionf(ia)%typ.lt.0) then

          cycle

        ! Number
        else if (Input%ionf(ia)%typ.eq.1) then

          NIV = NIV + 1
          ion_value_ind(ia) = NIV
          ion_value(NIV) = Input%ionf(ia)%val

        ! File
        else if (Input%ionf(ia)%typ.eq.0) then

          NIC = NIC + 1
          ion_column_ind(ia) = NIC

        ! What?!
        else

          umsg = 'Unexpected Atom_ion input'
          call verbose
          laborted = .True.

        end if

      end do

      ! Allocate units to read ion if it proceeds
      if (NIC.gt.0) then
        allocate(unitI(NIC))
        do ia=1,nA
          ip = ion_column_ind(ia)
          if (ip.lt.1) cycle
          if (ip.eq.1) then
            unitI(ip) = 10001
          else
            unitI(ip) = unitI(ip-1) + 1
          end if
          call open_file(unitI(ip), Input%ionf(ia)%str, &
                         0, .False., check)
        end do
      end if

      ! Check if aborting
      call gcontrol

      ! Initialize buffer sizes
      call set_io_CLE_buffers(Input,mode,dims,Frec)

      ! Check if aborting
      call gcontrol

      !
      ! Check files exist
      !
      if (gpid.eq.0) then

        ! Check files exist
        if (lcache) call check_io_CLE_buffers_exist(Input,aborting)
        laborted = aborting

      end if

      ! Check if aborting
      call gcontrol

      !
      ! SPECTRUM READ
      !
      call rinspect(Input,Atom,Geom,spect,Frec%omega)

      ! Check if aborting
      call gcontrol

      !
      ! Organize the tasks splitting (within groups)
      !
      if (gpid.gt.0) then

        ! Distribute tasks
        call setmpi_CLE(MPID,Input,Frec%IW_freq_in,Frec%IW_freq)

        ! Free IW_freq
        deallocate(Frec%IW_freq_in,Frec%IW_freq)

        !
        ! Finish the atom
        !
        call prepareatomMPI(Atom)
        call omegainitmaster(Atom)

      else if (Input%rt_group_n.gt.1) then

        umsg = ' - Tasks within groups distributed'
        call verbose

      end if

      ! Control
      call gcontrol

      !
      ! Photoionization quantites
      !
      if (gpid.gt.0) then
        do ia=1,nA
          call setphoto(Atom(ia),Frec%omega,MPID)
        end do
      else
        call gcontrol
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbose
      end if

      ! Only slaves
      if (gpid.gt.0) then

        ! Manage output frequency axis
        call refitfrec(Input,Frec,Atom,MPID)

        ! Allocate space for exu in Frec
        if (Frec%Mpif0(pid).ge.Frec%Mpif1(pid)) then
          allocate(Frec%exu(Frec%Mpif0(pid):Frec%Mpif1(pid),1))
        else
          nullify(Frec%exu)
        end if

      end if

      ! Check if aborting
      call gcontrol

      ! If cartesian or slab
      if (mode.eq.0.or.mode.eq.1) then

        ! Prepare buffers for data sending
        ioffset = dims(1)*22
        sizeA = ioffset + dims(1)*NIC
        if (mode.eq.1) then
          sizeA = sizeA+3
          ioffset = ioffset + 3
        end if
        allocate(buffer(sizeA))

        ! If cartesian
        if (mode.eq.0) then

          ! Get axes if master
          call get_axes(unitA,dims,x,y,z,double,gpid.eq.0)

          ! And share
          call MPI_BCAST(x(1),dims(1),MPI_DOUBLE_PRECISION,0, &
                         MPI_COMM_WORLD,ierr)
          call MPI_BCAST(y(1),dims(2),MPI_DOUBLE_PRECISION,0, &
                         MPI_COMM_WORLD,ierr)
          call MPI_BCAST(z(1),dims(3),MPI_DOUBLE_PRECISION,0, &
                         MPI_COMM_WORLD,ierr)

          ! Need to normalize?
          if (norm.le.0) then
            x = x/Input%R_star
            y = y/Input%R_star
            z = z/Input%R_star
          end if
        end if
      end if

      !
      ! Read partition function data and abundances
      !
      call rParfunAbund(Input,Atmo)

      ! Check if aborting
      call gcontrol

      !
      ! Read barklem data
      !
      if (gpid.ne.0) call rBarklem(Input,Atom,Atomb)

      ! Check if aborting
      call gcontrol

      !
      ! Get Kurucz lines
      !
      if (gpid.gt.0) then

        ! Read Kurucz data
        if (Input%NK.ge.1) then

          DwTa = Input%dw*1d-9/c

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,.False.,kurucz)

        ! If there are
        else

          kurucz%ntran = 0

        end if
      end if

      ! Check if aborting
      call gcontrol

      !
      ! CHIANTI READ
      !
      call rCHIANTI(Input,chianti)

      ! Check if aborting
      call gcontrol

      !
      ! Properly start now
      !
      if (gpid.eq.0) then
        umsg = ' - Starting to work on your model, this may '// &
               'take a while'
        call verbose
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Master
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (gpid.eq.0) then

        ! Open files so slaves can write later
        if (.not.lcache) call create_io_CLE_files(Input,mode,y,z, &
                                                  dims,Frec)

        ! Deallocate x,y,z
        if (allocated(x)) deallocate(x,y,z)

        ! Allocate cpu_free with group status
        allocate(cpu_free(MPID%ngroup))
        cpu_free = 1

        ! Initialize indexes and sizes
        iy = 1
        iz = 0
        inod = 0
        NLOSr = 0
        NLOS = dims(2)*dims(3)

        ! Work until exhausted
        do while (.True.)

          ! If aborting, send stop signal to everyone
          if (aborting) then

            ! For every CPU still sending
            do while (minval(cpu_free).lt..5d0)

              ! For every leader
              do ip=1,MPID%ngroup

                ! If already free, skip
                if (cpu_free(ip).gt..5d0) cycle

                ! Test if slave in the group is sending something
                call MPI_IPROBE(MPID%ltslave(ip), &
                                4000000+MPID%ltslave(ip), &
                                MPI_COMM_WORLD, receiving, &
                                MPI_STATUS_IGNORE, ierr)

                ! If slave is calling
                if (receiving) then

                  ! Receive the ping
                  call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                MPID%ltslave(ip), &
                                4000000+MPID%ltslave(ip), &
                                MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! Free the group
                  cpu_free(ip) = 1

                end if ! Receiving from a CPU

              end do ! Receive from everyone
            end do ! While there is someone working

          end if ! Aborting

          ! If there are LOS to do and at least one free CPU
          if (inod.lt.NLOS.and.maxval(cpu_free).gt..5d0) then

            ! Advance one node
            iy0 = iy
            iz0 = iz
            inod = inod + 1
            iz = iz + 1
            if (iz.gt.dims(3)) then
              iy = iy + 1
              iz = 1
              if (iy.gt.dims(2)) cycle
            end if
            icoords = (/ iy , iz , inod /)

            ! Check if not looping
            if (iy.ne.iy1.or.iz.ne.iz1) then

              ! If not cartesian, get coords
              if (mode.eq.2) then

                ! Get dimension and coordinates
                call get_point(unitA,dims(1),ycoor,zcoor,double)

                ! Correct icoords message
                icoords = (/ iy , iz , dims(1) /)

                ! Allocate buffer
                sizeA = dims(1)*(23 + NIC) + 2
                if (.not.allocated(buffer)) then
                  allocate(buffer(sizeA))
                else if (size(buffer).ne.sizeA) then
                  deallocate(buffer)
                  allocate(buffer(sizeA))
                end if

                ! Offset for ion
                ioffset = dims(1)*23 + 2

                ! Get column from atmosphere
                buffer(1) = ycoor
                buffer(2) = zcoor
                call get_column(unitA, buffer(3:ioffset), &
                                double, check)

              ! Cartesian or slab
              else

                ! Get column from atmosphere
                call get_column(unitA, buffer(1:ioffset), &
                                double, check)

              end if

              ! Get columns from ion
              if (NIC.gt.0) call get_column_ion(unitI,buffer, &
                                                dims(1),ioffset, &
                                                double,check)

              ! Store last read
              iy1 = iy
              iz1 = iz

            end if ! Loopìng

            ! If done in cache, skip
            if (lcache) then
              if (cache(iz,iy)) then
                if (.not.lload) then
                  NLOSr = NLOSr + 1
                  cycle
                end if
              end if
            end if

            ! Take a free cpu
            ip = maxloc(cpu_free, 1)

            ! Send signal to node
            call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                          MPID%ltslave(ip), &
                          2000000+MPID%ltslave(ip), &
                          MPI_COMM_WORLD, ierr)

            ! If failed
            if (ierr.ne.0) then
              inod = inod - 1
              iy = iy0
              iz = iz0
              cycle
            end if

            ! Send data to node
            call MPI_SEND(buffer(1), sizeA, &
                          MPI_DOUBLE_PRECISION, &
                          MPID%ltslave(ip), &
                          3000000+MPID%ltslave(ip), &
                          MPI_COMM_WORLD, ierr)

            ! If failed
            if (ierr.ne.0) then
              inod = inod - 1
              iy = iy0
              iz = iz0
              cycle
            end if

            ! That CPU is now busy
            cpu_free(ip) = 0

          end if ! If there is work to do and free CPUs

          ! For every slave group
          ip = 1
          do while (.True.)

            ! Test if slave in the group is sending something
            call MPI_IPROBE(MPID%ltslave(ip), &
                            4000000+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, receiving, &
                            MPI_STATUS_IGNORE, ierr)

            ! If failing
            if (ierr.ne.0) cycle

            ! If slave is calling
            if (receiving) then

              ! Try to receive
              do while (.True.)

                ! Receive the ping
                call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                              MPID%ltslave(ip), &
                              4000000+MPID%ltslave(ip), &
                              MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If failed, try again
                if (ierr.ne.0) cycle

                exit

              end do

              ! Convert iy coordinate into node coordinate
              int_buff(1) = (int_buff(1)-1)*dims(3) + int_buff(2)

              ! Write in cache
              call write_cache(unitC,Input%cache,int_buff,check)
              aborting = .not.check

              ! Update NLOS received
              NLOSr = NLOSr + 1

              ! Free the group
              cpu_free(ip) = 1

            end if ! Receiving from a CPU group

            ip = ip + 1

            if (ip.gt.MPID%ngroup) exit

          end do ! Slaves

          ! If we went beyond the number of LOS, exit
          if (NLOSr.ge.NLOS) exit

        end do ! While there is work to do

        ! If we are done, notify to slaves
        icoords(1) = -1
        iproc = 1
        do while (.True.)

          ! send termination signal
          call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                        MPID%ltslave(iproc), &
                        2000000+MPID%ltslave(iproc), &
                        MPI_COMM_WORLD, ierr)

          ! If it fails
          if (ierr.ne.0) cycle

          iproc = iproc + 1

          if (iproc.gt.MPID%ngroup) exit

        end do ! slaves

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Slaves
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else

        ! Initialize
        aborting = .False.

        ! Initialize atmosphere RAM
        call rAtmo_cle_prep(Input,Atmo,mode,norm)

        ! Normalize
        call normalize_cle(Atom)

        ! Work until further notice
        do while (.True.)

          ! If leader
          if (pid.eq.0) then

            ! Try receiving until success
            do while (.True.)

              ! Wait for signal
              call MPI_RECV(icoords(1), 3, MPI_INTEGER, 0, &
                            2000000+gpid, MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE, ierr)

              ! If it fails
              if (ierr.ne.0) cycle

              ! Check if it is the termination signal
              if (icoords(1).lt.1) then
                aborting = .True.
                exit
              end if

              ! If not-cartesian
              if (mode.eq.2) then

                ! Allocate buffer
                dims(1) = icoords(3)
                ioffset = dims(1)*23 + 2
                sizeA = ioffset + dims(1)*NIC
                if (allocated(buffer)) then
                  if (size(buffer).lt.sizeA) then
                    deallocate(buffer)
                    allocate(buffer(sizeA))
                  end if
                else
                  allocate(buffer(sizeA))
                end if

              end if ! Only non-cartesian grid

              ! Receive LOS
              call MPI_RECV(buffer(1), sizeA, &
                            MPI_DOUBLE_PRECISION, &
                            0, 3000000+gpid, MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE, ierr)

              ! If it fails
              if (ierr.ne.0) cycle

              ! Success
              exit

            end do

          end if

          ! If liutenant has friends
          if (nproc.gt.1) then

            ! Try until done
            do while (.True.)

              ! Broadcast
              call MPI_BCAST(icoords(1), 3, MPI_INTEGER, 0, &
                             MPI_COMM_RT, ierr)

              ! If it fails
              if (ierr.ne.0) cycle

              ! If aborting
              if (icoords(1).lt.1) then
                aborting = .True.
                exit
              end if

              ! If not-cartesian and a slave
              if (mode.eq.2.and.pid.ne.0) then

                ! Allocate buffer
                dims(1) = icoords(3)
                ioffset = dims(1)*23 + 2
                sizeA = ioffset + dims(1)*NIC
                if (allocated(buffer)) then
                  if (size(buffer).lt.sizeA) then
                    deallocate(buffer)
                    allocate(buffer(sizeA))
                  end if
                else
                  allocate(buffer(sizeA))
                end if

              end if ! Only non-cartesian grid

              ! Broadcast
              call MPI_BCAST(buffer(1), sizeA, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

              ! If it fails
              if (ierr.ne.0) cycle

              ! Success
              exit

            end do ! Try communicating until successfull

          end if ! MPI in RT

          ! Problem
          if (aborting) exit

          ! Buffer with ion data
          if (NIC.lt.1) then
            buffer_ion => buffer(1:1)
          else
            buffer_ion => buffer(ioffset+1:sizeA)
          end if

          !
          ! Call solver with correct arguments
          !

          ! If cartesian
          if (mode.eq.0) then

            buffer_atmo => buffer(1:ioffset)

            ! Solve columns
            call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Geom,Flgsg, &
                     fudge,kurucz,buffer_atmo,x, &
                     y(icoords(1)),z(icoords(2)), &
                     dims,buffer_ion,ion_column_ind, &
                     ion_value_ind,ion_value,spect,chianti)

          ! If slab
          else if (mode.eq.1) then

            buffer_atmo => buffer(1:ioffset)

            ! Solve columns
            call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Geom,Flgsg, &
                     fudge,kurucz,buffer_atmo,x,0d0,0d0, &
                     dims,buffer_ion,ion_column_ind, &
                     ion_value_ind,ion_value,spect,chianti)

          ! If non-cartesian
          else if (mode.eq.2) then

            ! Get true inode index
            icoords(3) = icoords(1) + (icoords(2)-1)*dims(2)

            buffer_atmo => buffer(3:ioffset)

            ! Solve columns
            call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Geom,Flgsg, &
                     fudge,kurucz,buffer_atmo,x,buffer(1), &
                     buffer(2),dims,buffer_ion,ion_column_ind, &
                     ion_value_ind,ion_value,spect,chianti)

          end if

          ! Nullify pointer
          nullify(buffer_ion)
          nullify(buffer_atmo)

          ! If liutenant, send message to grand master
          if (pid.eq.0) then

            ! Fail
            if (laborted) then
              icoords(3) = -1
            ! Success
            else
              icoords(3) = gpid
            end if

            ! Try until achieved
            do while (.True.)
              call MPI_SEND(icoords(1),3,MPI_INTEGER,0, &
                            4000000+gpid,MPI_COMM_WORLD,ierr)
              if (ierr.ne.0) cycle
              exit
            end do


          end if ! Send info back to grand master

        end do

      end if ! Master or slave

      !
      ! Close files
      !
      if (gpid.eq.0) then
        call close_file(unitA)
        do ip=1,NIC
          call close_file(unitI(ip))
        end do
      end if

      end subroutine HanleCLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module hanlert_mod
