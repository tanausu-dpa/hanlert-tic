      !> HanleRT manager
      module hanlert_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!  Start:
!     22/06/2022
!  Last version:
!     18/03/2025 V4.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/03/2025:    V4.0.2 - Added argument to the call to
!                             setmpi_sizes (TdPA)
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
!    - In Hanlert15D, allow for inversion outputs to be used as input
!      model atmosphere
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!    HanleRT1DS
!      Prepare the model atmosphere and run the solution in a single
!    1D model atmosphere
!
!    HanleRT15DS
!      Manage the solution pixel by pixel of a 3D model atmosphere
!
!    HanleRTTIC
!      Manage the solution in a model atmosphere in the inversion
!
!    HanleCLE
!      Manage the solution pixel by pixel of a 3D model atmosphere or
!    a list of slabs for the Coronal Line Emission problem
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
      use parameters_mod , only : c, pi, TINYB , TINYVEL
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

      !> Prepare the model atmosphere and run the solution in a single
      !! 1D model atmosphere\n
      !!    Input(Input_class): Structure with configuration data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine HanleRT1DS(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Input_class), intent(inout):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(inout):: fudge
      type(MPI_class), intent(inout):: MPID

      ! Local

      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI,Geom
      type(kurucz_class):: kurucz
      type(Solution_F_class):: dummy2
      type(LTEline_class), dimension(:), allocatable:: LTElines

      logical:: lDK,l1,l2,lio,lie,lp,lpe,lporlpe,lload

      integer:: ia

      double precision:: maxB,DwTa
      double precision, dimension(:), allocatable:: dummy


      !
      ! Define RT communicator. Here it is the same than WORLD
      ! because everyone works in the same (and single) pixel
      !
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, gpid, &
                          MPI_COMM_RT, ierr)
      call MPI_COMM_RANK(MPI_COMM_RT,pid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_RT,nproc,ierr)

      ! Count non-allocatable memory in local structures
      MRAMc = MRAMc + 1d-6*(sizeof(Atmo) + &
                            sizeof(Bfield) + &
                            sizeof(Frec) + &
                            sizeof(GeomI) + &
                            sizeof(Geom) + &
                            sizeof(kurucz))
      SRAMc = SRAMc + 1d-6*sizeof(dummy2)


      !
      ! Define the flags to decide the flow of the code
      !

      ! If loading file
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'

      ! If iterating intensity
      lio = (Input%mode.eq.'W'.and.Input%force.ne.'P').or. &
            (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'I'.or.Input%force.eq.'A'))

      ! If emerging intensity
      lie = Input%force.eq.'I'.and.Input%nThLOS.gt.0

      ! If iterating polarization
      lp = (Input%mode.eq.'W'.and.Input%force.ne.'I').or. &
           (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'N'.or.Input%force.eq.'P'.or. &
                      Input%force.eq.'A'))

      ! If emerging polarization
      lpe = Input%force.ne.'I'.and.Input%nThLOS.gt.0

      ! If doing anything polarization related
      lporlpe = lp.or.lpe


      !
      ! Allocate atomic MPI arrays
      call prepareatomMPI(Atom)

      !
      ! Read atmospheric data
      call rAtmo(Input%atmo,Input%source,Input%ID,Atmo,Input%fvmicro)

      ! Save dimension of read model in common variable
      nz = Atmo%nz

      ! If the model has velocities, but we are forcing static
      if (dyn.and.Input%static) then

        ! Make the problem static
        dyn = .False.
        Atmo%vx = 0d0
        Atmo%vy = 0d0
        Atmo%vz = 0d0

      end if ! Dynamic but forcing static

      !
      ! Read magnetic field data
      call rBField(Input%bfield,Input%source,Input%ID, &
                   Bfield,Atmo%nz,Input)

      !
      ! Index atom levels, sublevels, transitions, and
      ! magnetic components
      call set_atom_indexes(Atom,Input,lio.or.lie,lporlpe, &
                            minval(Bfield%Bstrength).le.TINYB, &
                            maxval(Bfield%Bstrength).gt.TINYB)

      !
      ! Read partition function data and abundances
      call rParfunAbund(Input,Atmo)

      !
      ! Read barklem data
      call rBarklem(Input,Atom,Atomb)

      !
      ! Get minimum and maximum temperatures, and maximum velocity
      !

      ! Temperature limits
      Input%minT = minval(Atmo%T)
      Input%maxT = maxval(Atmo%T)

      ! If dynamic in theory
      if (dyn) then

        ! Get maximum velocity
        Input%maxV = maxval(Atmo%vx*Atmo%vx + &
                            Atmo%vy*Atmo%vy + &
                            Atmo%vz*Atmo%vz)
        Input%maxV = sqrt(Input%maxV)

        ! Check if too small
        if (Input%maxV.le.TINYVEL) then

          ! Flag static
          Input%maxV = 0d0
          Input%static = .True.
          dyn = .False.

        ! Not too small
        else

          ! No static
          Input%static = .False.

        end if ! Maximum velocity is too small

      ! Model static by itself
      else

        ! Zero velocity
        Input%maxV = 0d0
        Input%static = .True.

      end if ! Determine maximum velocity

      ! Get maximum magnetic field
      maxB = maxval(Bfield%Bstrength)

      ! If smaller than small field
      if (maxB.le.TINYB) then

        ! Flag no magnetic
        Input%unmagnetized = .True.

      ! Larger than small field
      else

        ! Flag magnetic
        Input%unmagnetized = .False.

      end if ! Magnetic field

      !
      ! Decide if need to keep Stokes
      !

      ! Master
      if (pid.eq.0) then

        ! Keep Stokes if pre-imposed or doing PRD that
        ! is either dynamic or angle-dependent
        KSTK = KSTK.or.(PRD.and.(dyn.or..not.AV)).or.(lp.and.dyn)

      ! Slaves
      else

        ! Keep Stokes if doing PRD that is either dynamic
        ! or angle-dependent
        KSTK = (PRD.and.(dyn.or..not.AV)).or.(lp.and.dyn)

      end if ! Master/slave

      !
      ! Set angular quadrature
      call gauss(Input,GeomI,Geom,1,lp.or.lpe,lpe.or.lie,Flgsg)

      ! Verbose angular quadrature finished
      if (pid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbose
      end if

      !
      ! Check that the axial symmetry is consistent
      if (axial) &
        call check_axial(Atmo%vx,Atmo%vy,Bfield%Btheta)


      !
      ! Define the output frequency axis
      call omegabuild(Frec,Atom,Input,maxB,dummy)

      ! Initialize MPI arrays for Master
      call omegainitmaster(Atom)

      ! Verbose finished with output frequency axis
      if (pid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbose
      end if

      !
      ! Organize the tasks splitting
      call setmpi(MPID,Input,Frec%IW_freq)

      ! Verbose finished distributing tasks
      if (pid.eq.0.and.MPID%mpi) then
        umsg = ' - Tasks distributed'
        call verbose
      end if

      !
      ! Master can remove fudge data if doing MPI
      if (pid.eq.0.and.MPID%mpi.and.allocated(fudge%fudge_v)) then
        MRAMc = MRAMc - 1d-6*sizeof(fudge%fudge_v)
        deallocate(fudge%fudge_v)
      end if

      !
      ! Kurucz
      !

      !
      ! If slave or serial
      if (pid.gt.0.or..not.MPID%mpi) then

        ! There are Kurucz data to read
        if (Input%NK.ge.1) then

          ! If Doppler width for transformation is the maximum
          if(Input%dws.eq.'MAX')then

            ! Flag and compute maximum velocity (sqrt)
            lDK = .True.
            DwTa = sqrt(maxval(Atmo%T))

          ! If Doppler width for transformation is the minimum
          else if(Input%dws.eq.'MIN')then

            ! Flag and compute minimum velocity (sqrt)
            lDK = .True.
            DwTa = sqrt(minval(Atmo%T))

          ! If Doppler width for transformation is a number
          else if(Input%dws.eq.'NUM')then

            ! Flag and compute factor
            lDK = .False.
            DwTa = Input%dw*1d-9/c

          end if ! Type of Doppler width transformation

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,lDK,kurucz)

        ! If there are no Kurucz data to read
        else

          ! Set number of transitions in Kurucz structure
          ! to zero
          kurucz%ntran = 0

        end if ! There is Kurucz data to read
      end if ! Serial or slave

      ! Control
      call gcontrol

      !
      ! Photoionization quantites
      !

      ! For every active atom
      do ia=1,nA

        ! Set photoionization quantities
        call setphoto(Atom(ia),Frec%omega)

      end do ! Active atoms

      ! Verbose initialized photoionization cross sections
      if (pid.eq.0) then
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbose
      end if


      !
      ! Resize some frequency quantities
      call frecresize(Frec,Atom,Input,MPID)


      !
      ! Compute size for MPI messages in solvers
      !

      ! Flags to decide what sizes are necessary to define
      l1 = .not.lload
      l2 = (lio.and..not.lie).or.(lload.and.(lp.or.lpe))

      ! If doing MPI, setup the buffer sizes
      if (MPID%mpi) &
        call setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,l1,l2, &
                          Input%ALI_photo,.False.)

#ifdef DEBUGATMO
      ! Write in an ASCII file the model atmosphere
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif

      ! Prepare model atmosphere for synthesis
      call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

      ! If error, skip
      if (laborted) goto 2000

      !
      ! Solution of the NLTE problem
      !
      call hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input,GeomI, &
                 Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                 dummy,dummy2,lload,lio,lie,lp,lpe,.True.)

      ! If error, skip
      if (laborted) goto 2000

      ! Verbose that the synthesis finished
      if (pid.eq.0) then
        umsg = ' - Solver finished'
        call verbose
      end if

      ! Recalculate electrons
      call updateatmo(Atom,Atomb,Atmo,Bfield,Input)

      ! Clean memory
2000  call free_pix(Atom,Atomb,Mol,Bfield)
      call free_Atmo(Atmo,.True.)
      call free_LTElines_full(LTElines)

      end subroutine HanleRT1DS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage the solution pixel by pixel of a 3D model atmosphere\n
      !!    Input(Input_class): Structure with configuration data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine HanleRT15DS(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Input_class), intent(inout):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(inout):: fudge
      type(MPI_class), intent(inout):: MPID

      ! Local

      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI, Geom
      type(kurucz_class):: kurucz
      type(Solution_F_class):: dummy
      type(LTEline_class), dimension(:), allocatable:: LTElines

      logical:: aborting,check,receiving,lcache,double,lexcl,lDK
      logical:: lio,lie,lp,lpe,lporlpe,lload,l1,l2,warning
      logical, dimension(:,:), allocatable:: cache

      integer:: unitA,unitC,unitJ
      integer:: ix1,iy1,sizeA,sizeJ
      integer:: ia,ix0,iy0,ix,iy,inod,NLOSr,NLOS,ip,iproc
      integer, dimension(3):: dims,out_dims,int_buff
      integer, dimension(:), allocatable:: cpu_free

      double precision:: maxB,DwTa,lMRAMc
      double precision, dimension(:), allocatable,target:: buffer_atmo
      double precision, dimension(:), allocatable,target:: buffer_JKQ


      !
      ! Initialize
      !

      ! Can issue memory RAM count warning
      warning = .True.

      ! No info on last MRAMc
      lMRAMc = -1d0

      ! Count non-allocatable memory in local structures
      MRAMc = MRAMc + 1d-6*(sizeof(Atmo) + &
                            sizeof(Bfield) + &
                            sizeof(Frec) + &
                            sizeof(GeomI) + &
                            sizeof(Geom) + &
                            sizeof(kurucz))
      SRAMc = SRAMc + 1d-6*sizeof(dummy)

      ! If unmagnetized by configuration
      if (Input%unmagnetized) then

        ! Force maximum field to zero
        maxB = 0d0

      ! Not forcing zero magnetic field
      else

        ! Assume that there will be some magnetic field
        maxB = TINYB + 1d0

      end if ! Forcing magnetic field

      ! If (PRD and T limits unknown) or T limits are unknown and
      ! needed for the Doppler width to transform frequencies, or
      ! if not forcing static but the maximum velocity is unknown
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

      ! If forcing static
      if (Input%static) then

        ! Force the maximum velocity to zero
        Input%maxV = 0d0

      ! Not forcing static
      else

        ! Get maximum velocity in synthesis units
        Input%maxV = Input%maxV*1d-6/c

      end if ! Forcing static


      !
      ! Define the flags to route the code
      !

      ! If reading previous solutions
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'

      ! If iterating the intensity problem
      lio = (Input%mode.eq.'W'.and.Input%force.ne.'P').or. &
            (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'I'.or.Input%force.eq.'A'))

      ! If doing intensity emergence
      lie = Input%force.eq.'I'.and.Input%nThLOS.gt.0

      ! If iterating the polarization problem
      lp = (Input%mode.eq.'W'.and.Input%force.ne.'I').or. &
           (Input%mode.eq.'B'.and. &
                     (Input%force.eq.'N'.or.Input%force.eq.'P'.or. &
                      Input%force.eq.'A'))

      ! If doing polarization emergence
      lpe = Input%force.ne.'I'.and.Input%nThLOS.gt.0

      ! If doing anything in polarization
      lporlpe = lp.or.lpe


      !
      ! Index atom levels, sublevels, transitions, and
      ! magnetic components
      call set_atom_indexes(Atom,Input,lio.or.lie,lporlpe, &
                            .True.,maxB.gt.0d0)


      !
      ! Define the output frequency axis
      call omegabuild(Frec,Atom,Input,maxB,buffer_atmo)

      ! Verbose that the frequency axis is defined
      if (gpid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbose
      end if


      !
      ! Split in groups of tasks to distribute pixels
      call setmpi15D(MPID,Input)

      ! Verbose MPI splits
      if (gpid.eq.0) then
        umsg = ' - Tasks distributed'
        call verbose
      end if


      !
      ! Decide if need to keep Stokes
      !

      ! If not splitting in pixels or it is an RT master
      if (.not.MPID%mpi15d.or.(pid.eq.0.and.gpid.ne.0)) then

        ! Keep Stokes if already flagged or if doind PRD either
        ! dynamic or angle-dependent
        KSTK = KSTK.or.(PRD.and.(.not.Input%static.or..not.AV)).or. &
               (lp.and.dyn)

      ! If RT slave
      else

        ! Keep Stokes if doind PRD either dynamic or angle-dependent
        KSTK = (PRD.and.(.not.Input%static.or..not.AV)).or. &
               (lp.and.dyn)

      end if ! RT slave


      !
      ! Set angular quadrature
      call gauss(Input,GeomI,Geom,1,lp.or.lpe,lpe.or.lie,Flgsg)

      ! Verbose angular quadrature
      if(gpid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbose
      end if


      !
      ! Global master remove fudge if splitting pixels
      if (gpid.eq.0.and.allocated(fudge%fudge_v).and.MPID%mpi15d) then
        MRAMc = MRAMc - 1d-6*sizeof(fudge%fudge_v)
        deallocate(fudge%fudge_v)
      end if


      !
      ! Pre-process
      !

      ! Initialize looping identifier
      ix1 = -1
      iy1 = -1

      ! Global Master
      if (gpid.eq.0) then

        ! Open files (ia and nz are a dummy variable here)
        call open_atm_and_cache(Input,1,unitA,unitC,aborting,dims, &
                                ia,double,nz,cache,lcache)

        ! The master does not need the background atoms or
        ! molecules if just managing
        if (MPID%mpi15d) then
          call free_atom_full(Atomb)
          call free_mol_full(Mol)
        end if

        ! Check if could read
        laborted = aborting

        ! Open JKQ asymmetry file
        if (Input%nasym.gt.0.and..not.laborted) then

          ! Open asymmetry file
          call open_asymm(Input,unitJ,aborting,dims)

          ! Check if could read
          laborted = aborting

        end if ! JKQ asymmetry file
      end if ! Global master

      ! Check if aborting
      call gcontrol

      ! Share model dimensions
      call MPI_BCAST(dims(1),3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ! Save z axis dimension into common variable
      nz = dims(3)

      !
      ! Now fix the real dimensions from solution box
      !

      ! Fix wildcards in input
      if (Input%sol_box(1).lt.1) Input%sol_box(1) = 1
      if (Input%sol_box(2).lt.1) Input%sol_box(2) = dims(1)
      if (Input%sol_box(3).lt.1) Input%sol_box(3) = 1
      if (Input%sol_box(4).lt.1) Input%sol_box(4) = dims(2)

      ! Global Master
      if (gpid.eq.0) then

        ! Get output dimensions (actual sizes to save)
        out_dims(1) = Input%sol_box(2) - Input%sol_box(1) + 1
        out_dims(2) = Input%sol_box(4) - Input%sol_box(3) + 1
        out_dims(3) = dims(3)

      ! Rest
      else

        ! Correct dimensions
        dims(1) = Input%sol_box(2) - Input%sol_box(1) + 1
        dims(2) = Input%sol_box(4) - Input%sol_box(3) + 1

      end if ! Global Master


      !
      ! Organize the tasks splitting (within groups)
      !

      ! If not splitting pixels or not the distributer
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! Split taks in pixel RT
        call setmpi(MPID,Input,Frec%IW_freq)

        ! RT Masters remove fudge
        if ((pid.eq.0.and.MPID%mpi).and.allocated(fudge%fudge_v)) then
          MRAMc = MRAMc - 1d-6*sizeof(fudge%fudge_v)
          deallocate(fudge%fudge_v)
        end if

        !
        ! Finish preparing the atom
        call prepareatomMPI(Atom)
        call omegainitmaster(Atom)

      end if ! Not splitting pixels or not the distributer

      ! Verbose MPI RT distributed
      if (gpid.eq.0) then
        umsg = ' - Tasks within groups distributed'
        call verbose
      end if

      ! Control
      call gcontrol


      !
      ! Photoionization quantites
      !

      ! If not splitting pixels or not the distributer
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! For every active atom
        do ia=1,nA

          ! Setup photoionization cross sections
          call setphoto(Atom(ia),Frec%omega)

        end do ! Active atoms

      ! The distributer
      else

        ! Just call control to match the call in setphoto
        call gcontrol

      end if ! Not splitting pixels or not the distributer

      ! Verbose photoionization quantities initialized
      if (gpid.eq.0) then
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbose
      end if


      !
      ! MPI sizes
      !

      ! Not global Master or not distributer
      if (.not.MPID%mpi15d.or.gpid.gt.0)  then

        ! Resize some frequency quantities if doing MPI
        call frecresize(Frec,Atom,Input,MPID)

        !
        ! Compute size for MPI messages in solvers
        !

        ! Flags to decide what sizes are necessary to define
        l1 = .not.lload
        l2 = (lio.and..not.lie).or.(lload.and.(lp.or.lpe))

        ! If doing RT MPI, define the buffer sizes
        if (MPID%mpi) &
          call setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,l1,l2, &
                            Input%ALI_photo,.False.)

      end if ! Not global master or not distributer


      !
      ! Prepare buffers
      !

      ! Sanity check the limits for the outputs
      call check_io_buffers_sanity_check(Input,Atom)

      ! Check if aborting
      call gcontrol

      ! Prepare buffer atmosphere
      sizeA = dims(3)*24
      allocate(buffer_atmo(sizeA))
      MRAMc = MRAMc + 1d-6*sizeof(buffer_atmo)

      ! If asymmetry data
      if (Input%nasym.gt.0) then
        sizeJ = dims(3)*8
        allocate(buffer_JKQ(sizeJ))
        MRAMc = MRAMc + 1d-6*sizeof(buffer_JKQ)
      else
        sizeJ = 0
      end if

      ! Initialize buffer sizes for IO
      call set_io_buffers(Input,0,Atom,Frec)

      ! Check if aborting
      call gcontrol

      ! Global master
      if (gpid.eq.0) then

        ! Check files exist
        if (lcache) call check_io_buffers_exist(Input,Atom, &
                                                Geom,aborting)
        ! Check if success in checking
        laborted = aborting

      end if ! Global master

      ! Check if aborting
      call gcontrol


      !
      ! Read partition function data and abundances
      call rParfunAbund(Input,Atmo)

      ! Check if aborting
      call gcontrol


      !
      ! Read barklem data
      if (.not.MPID%mpi15d.or.gpid.gt.0) &
        call rBarklem(Input,Atom,Atomb)

      ! Check if aborting
      call gcontrol


      !
      ! Get Kurucz lines
      !

      ! If not distributing pixels and not the RT master doing MPI,
      ! or if not the distributer nor the RT master doing MPI
      if ((.not.MPID%mpi15d.and.(pid.gt.0.or..not.MPID%mpi)).or. &
          (gpid.gt.0.and.(pid.gt.0.or..not.MPID%mpi))) then

        ! If there is Kurucz data to read
        if (Input%NK.ge.1) then

          ! If Doppler width for transformation is the maximum
          if(Input%dws.eq.'MAX')then

            ! Flag and compute maximum velocity (sqrt)
            lDK = .True.
            DwTa = sqrt(Input%maxT)

          ! If Doppler width for transformation is the minimum
          else if(Input%dws.eq.'MIN')then

            ! Flag and compute minimum velocity (sqrt)
            lDK = .True.
            DwTa = sqrt(Input%minT)

          ! If Doppler width for transformation is a number
          else if(Input%dws.eq.'NUM')then

            ! Flag and compute factor
            lDK = .False.
            DwTa = Input%dw*1d-9/c

          end if ! Type of Doppler width transformation

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,lDK,kurucz)

        ! If there are not Kurucz data
        else

          ! Set number of transitions in Kurucz structure
          ! to zero
          kurucz%ntran = 0

        end if ! There is Kurucz data
      end if ! Could need Kurucz data

      ! Check if aborting
      call gcontrol

      !
      ! Properly start now
      !

      ! Global master anounce start
      if (gpid.eq.0) then
        umsg = ' - Starting to work on your model, this may '// &
               'take a while'
        call verbose
      end if

      !
      ! MPI version
      !

      ! Distributing tasks
      if (MPID%mpi15d) then

        !
        ! Global Master
        !
        if (gpid.eq.0) then

          ! Open files so slaves can write later
          if (.not.lcache) call create_io_files(Input,Atom,out_dims, &
                                                Geom,GeomI,Frec)

          ! Allocate cpu_free with group status
          allocate(cpu_free(MPID%ngroup))
          MRAMc = MRAMc + 1d-6*sizeof(cpu_free)
          cpu_free = 1

          ! Initialize indexes and sizes
          ix = 1
          iy = 0
          inod = 0
          NLOSr = 0
          NLOS = out_dims(1)*out_dims(2)

          ! Work until exhausted
          do while (.True.)

            ! If aborting
            if (aborting) then

              !
              ! Send stop signal to everyone
              !

              !
              ! First receive whatever is pending

              ! For every CPU still sending
              do while (minval(cpu_free).lt..5d0)

                ! For every leader
                do ip=1,MPID%ngroup

                  ! If already free, skip
                  if (cpu_free(ip).gt..5d0) cycle

                  ! Test if slave in the group is sending something
                  call MPI_IPROBE(MPID%ltslave(ip), &
                                  4+MPID%ltslave(ip), &
                                  MPI_COMM_WORLD, receiving, &
                                  MPI_STATUS_IGNORE, ierr)

                  ! If slave is calling
                  if (receiving) then

                    ! Receive the ping
                    call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                  MPID%ltslave(ip), &
                                  4+MPID%ltslave(ip), &
                                  MPI_COMM_WORLD, &
                                  MPI_STATUS_IGNORE, ierr)

                    ! Free the group
                    cpu_free(ip) = 1

                  end if ! Receiving from a CPU

                end do ! Receive from everyone
              end do ! While there is someone working

              ! And break the work loop
              exit

            end if ! Aborting

            ! If there are LOS to do and at least one free CPU
            if (inod.lt.NLOS.and.maxval(cpu_free).gt..5d0) then

              ! Save last values for pixel position
              ix0 = ix
              iy0 = iy

              ! Advance one pixel (X slow, Y fast)
              iy = iy + 1
              if (iy.gt.dims(2)) then
                ix = ix + 1
                iy = 1
                if (ix.gt.dims(1)) cycle
              end if

              ! Check if not repeating
              if (ix.ne.ix1.or.iy.ne.iy1) then

                ! Get column from atmosphere
                call get_column(unitA, buffer_atmo, double, check)

                ! If JKQ data, read
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

              ! Now advance LOS
              inod = inod + 1

              ! Set slave indexes
              icoords = (/ ix - Input%sol_box(1) + 1 , &
                           iy - Input%sol_box(3) + 1 , inod /)

              ! If done in cache, skip
              if (lcache) then
                if (cache(iy,ix)) then
                  NLOSr = NLOSr + 1
                  cycle
                end if
              end if

              ! If excluding pixels
              if (Input%lexcl) then

                ! Check if potentially excluded
                if (ix.ge.Input%excl(1,1).and. &
                    ix.le.Input%excl(1,Input%nexcl)) then

                  ! Initialize excluded flag
                  lexcl = .False.

                  ! For every exclusion
                  do ia=1,Input%nexcl

                    ! If below the current one, skip
                    if (Input%excl(1,ia).lt.ix) cycle

                    ! If larger than the current one, we went beyond
                    if (Input%excl(1,ia).gt.ix) exit

                    ! If y coincides
                    if (Input%excl(2,ia).eq.iy) then

                      ! Flag excluded
                      lexcl = .True.
                      exit

                    end if ! Y coincides

                  end do ! Exclusions

                  ! If excluded pixel, skip
                  if (lexcl) then
                    NLOSr = NLOSr + 1
                    cycle
                  end if

                end if ! Pixel potentially excluded
              end if ! Excluding pixels

              ! Take a free cpu
              ip = maxloc(cpu_free, 1)

              ! Send signal to node
              call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                            MPID%ltslave(ip), &
                            2+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, ierr)

              ! If failed
              if (ierr.ne.0) then

                ! Go back and try again
                inod = inod - 1
                ix = ix0
                iy = iy0
                cycle

              end if ! Message failed

              ! Send atmospheric data to node
              call MPI_SEND(buffer_atmo(1), sizeA, &
                            MPI_DOUBLE_PRECISION, &
                            MPID%ltslave(ip), &
                            3+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, ierr)

              ! If failed
              if (ierr.ne.0) then

                ! Go back and try again
                inod = inod - 1
                ix = ix0
                iy = iy0
                cycle

              end if ! Message failed

              ! If JKQ data
              if (sizeJ.gt.0) then

                ! Send JKQ data to node
                call MPI_SEND(buffer_JKQ(1), sizeJ, &
                              MPI_DOUBLE_PRECISION, &
                              MPID%ltslave(ip), &
                              5+MPID%ltslave(ip), &
                              MPI_COMM_WORLD, ierr)

                ! If failed
                if (ierr.ne.0) then

                  ! Go back and try again
                  inod = inod - 1
                  ix = ix0
                  iy = iy0
                  cycle

                end if ! Message failed
              end if ! There are JKQ data

              ! That CPU is now busy
              cpu_free(ip) = 0

            end if ! If there is work to do and free CPUs


            !
            ! Receive data
            !

            ! Initialize with first group
            ip = 1

            ! For every slave group
            do while (.True.)

              ! Test if slave in the group is sending something
              call MPI_IPROBE(MPID%ltslave(ip), &
                              4+MPID%ltslave(ip), &
                              MPI_COMM_WORLD, receiving, &
                              MPI_STATUS_IGNORE, ierr)

              ! If message failed, try again
              if (ierr.ne.0) cycle

              ! If slave is calling
              if (receiving) then

                ! Try to receive
                do while (.True.)

                  ! Receive the ping
                  call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                MPID%ltslave(ip), &
                                4+MPID%ltslave(ip), &
                                MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! If failed, try again
                  if (ierr.ne.0) cycle

                  ! Success, can leave
                  exit

                end do ! Try receiving

                ! Convert ix coordinate into node coordinate
                int_buff(1) = &
                    (int_buff(1) - 2 + Input%sol_box(1))*dims(2) + &
                    int_buff(2) - 1 + Input%sol_box(3)

                ! Write in cache
                call write_cache(unitC,Input%cache,int_buff,check)

                ! Check if could write
                aborting = .not.check

                ! Update NLOS received
                NLOSr = NLOSr + 1

                ! Free the group
                cpu_free(ip) = 1

              end if ! Receiving from a CPU group

              ! Advance group
              ip = ip + 1

              ! If checked every group, leave
              if (ip.gt.MPID%ngroup) exit

            end do ! Slaves

            ! If we went beyond the number of LOS, exit
            if (NLOSr.ge.NLOS) exit

          end do ! While there is work to do

          !
          ! Notify finished with work
          !

          ! Finished signal
          icoords(1) = -1

          ! First CPU
          iproc = 1

          ! Try until sucess
          do while (.True.)

            ! Send termination signal
            call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                          MPID%ltslave(iproc), &
                          2+MPID%ltslave(iproc), &
                          MPI_COMM_WORLD, ierr)

            ! If it fails try again
            if (ierr.ne.0) cycle

            ! Advance group
            iproc = iproc + 1

            ! If sent to all groups, leave
            if (iproc.gt.MPID%ngroup) exit

          end do ! Try sending messages until done

          ! Free cpu_free
          MRAMc = MRAMc - 1d-6*sizeof(cpu_free)
          deallocate(cpu_free)

        !
        ! Slaves
        !
        else

          ! Initialize
          aborting = .False.

          ! Work until further notice
          do while (.True.)

            ! If RT leader
            if (pid.eq.0) then

              ! Try receiving until success
              do while (.True.)

                ! Wait for signal
                call MPI_RECV(icoords(1), 3, MPI_INTEGER, 0, &
                              2+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! Check if it is the termination signal
                if (icoords(1).lt.1) then
                  aborting = .True.
                  exit
                end if

                ! Receive model atmosphere for pixel
                call MPI_RECV(buffer_atmo(1), sizeA, &
                              MPI_DOUBLE_PRECISION, &
                              0, 3+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! If JKQ data
                if (sizeJ.gt.0) then

                  ! Receive asymmetry JKQ for pixel
                  call MPI_RECV(buffer_JKQ(1), sizeJ, &
                                MPI_DOUBLE_PRECISION, &
                                0, 5+gpid, MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! If it fails, try again
                  if (ierr.ne.0) cycle

                end if ! There is JKQ data

                ! Success, leave loop
                exit

              end do ! Try receiving till success

            end if ! RT leader

            ! If MPI for the pixel
            if (nproc.gt.1) then

              ! Try until done
              do while (.True.)

                ! Broadcast coordinate
                call MPI_BCAST(icoords(1), 3, MPI_INTEGER, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! If aborting
                if (icoords(1).lt.1) then

                  ! Everyone leave
                  aborting = .True.
                  exit

                end if ! Aborting

                ! Broadcast model atmosphere for pixel
                call MPI_BCAST(buffer_atmo(1), sizeA, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! If JKQ data
                if (sizeJ.gt.0) then

                  ! Broadcast JKQ data
                  call MPI_BCAST(buffer_JKQ(1), sizeJ, &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                  ! If it fails, try again
                  if (ierr.ne.0) cycle

                end if

                ! Success, leave
                exit

              end do ! Try communicating until successfull

            end if ! MPI in RT

            ! Problem
            if (aborting) exit

            ! Build actual model atmosphere from buffer
            call rAtmo_frombuffer(buffer_atmo, Input, &
                                  Atmo, Bfield, dims)

            ! Initialize aborting flag
            laborted = .False.

#ifdef DEBUGATMO
            ! RT master write ASCII file with model atmosphere
            if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif
            ! Prepare model atmosphere for synthesis
            call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

            ! If error, skip
            if (laborted) goto 1100

            ! If no lMRAMc data
            if (lMRAMc.lt.0d0) then

              ! Set-up
              lMRAMc = MRAMc

            ! Check MRAMc
            else

              ! If different
              if (nint(1d6*abs(lMRAMc - MRAMc)).gt.1d0) then

                ! Warning
                if (warning) then

                  ! Deflag
                  warning = .False.

                  ! Write message
                  urou = 'HanleRT15DS'
                  write(umsg,'(2(A,es13.6),A)') &
                    'The miscellaneous RAM counter is different '// &
                    'between calls to the hanle function ',MRAMc, &
                    ' != ',lMRAMc,'. It is being corrected, but '// &
                    'this should not happen. Please, notify of '// &
                    'the issue providing your inputs'
                  call abortedS(umsg,urou,.False.,.True.)

                end if ! Can issue warning

                ! Correct
                MRAMc = lMRAMc

              end if ! Different
            end if ! lMRAMc data

            !
            ! Solve NLTE problem in pixel
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
            call free_LTElines_full(LTElines)

            ! If liutenant, send message to grand master
            if (pid.eq.0) then

              ! Failed
              if (laborted) then

                ! Send issue signal
                icoords(3) = -1

              ! Success
              else

                ! Send success signal
                icoords(3) = gpid

              end if ! Failure or success

              ! Try until achieved
              do while (.True.)

                ! Send status of pixel calculation
                call MPI_SEND(icoords(1),3,MPI_INTEGER,0, &
                              4+gpid,MPI_COMM_WORLD,ierr)

                ! If failed, try again
                if (ierr.ne.0) cycle

                ! Success, leave
                exit

              end do ! Try until success

            end if ! Send info back to grand master

          end do ! Work till finished

        end if ! Master or slave

      !
      ! Single pixel work version
      !
      else

        ! Open output files
        if (.not.lcache) call create_io_files(Input,Atom,out_dims, &
                                              Geom,GeomI,Frec)

        ! Initialize indexes and sizes
        ix = 1
        iy = 0
        inod = 0
        NLOSr = 0
        NLOS = out_dims(1)*out_dims(2)

        ! Initialize
        aborting = .False.

        ! Work until exhausted
        do while (.True.)

          ! If aborting, send stop signal to everyone
          if (aborting) then

            ! Abort
            call aborted_silent

          end if ! Aborting

          ! Save current position
          ix0 = ix
          iy0 = iy

          ! Advance one pixel
          iy = iy + 1
          if (iy.gt.dims(2)) then
            ix = ix + 1
            iy = 1
            if (ix.gt.dims(1)) cycle
          end if

          ! Check if not repeating
          if (ix.ne.ix1.or.iy.ne.iy1) then

            ! Get model atmosphere for this pixel
            call get_column(unitA, buffer_atmo, double, check)

            ! If JKQ data, get data for pixel
            if (sizeJ.gt.0) &
              call get_column(unitJ, buffer_JKQ, .True., check)

            ! Store last read
            ix1 = ix
            iy1 = iy

          end if ! Repeating pixel

          ! Skip if out of box
          if (ix.lt.Input%sol_box(1).or. &
              ix.gt.Input%sol_box(2).or. &
              iy.lt.Input%sol_box(3).or. &
              iy.gt.Input%sol_box(4)) then
            cycle
          end if

          ! Advance LOS
          inod = inod + 1

          ! Set slave indexes
          icoords = (/ ix - Input%sol_box(1) + 1 , &
                       iy - Input%sol_box(3) + 1 , inod /)

          ! If done in cache, skip
          if (lcache) then
            if (cache(iy,ix)) then
              NLOSr = NLOSr + 1
              cycle
            end if
          end if

          ! Build atmospheric model from buffer
          call rAtmo_frombuffer(buffer_atmo, Input, &
                                Atmo, Bfield, dims)

          ! Initialize aborting flag
          laborted = .False.

#ifdef DEBUGATMO
          ! RT master write atmospheric model in ASCII
          if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif
          ! Prepare atmospheric model for synthesis
          call prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

          ! If error, skip
          if (laborted) goto 1200

          ! If no lMRAMc data
          if (lMRAMc.lt.0d0) then

            ! Set-up
            lMRAMc = MRAMc

          ! Check MRAMc
          else

            ! If different
            if (nint(1d6*abs(lMRAMc - MRAMc)).gt.1d0) then

              ! Warning
              if (warning) then

                ! Deflag
                warning = .False.

                ! Write message
                urou = 'HanleRT15DS'
                write(umsg,'(2(A,es13.6),A)') &
                  'The miscellaneous RAM counter is different '// &
                  'between calls to the hanle function ',MRAMc, &
                  ' != ',lMRAMc,'. It is being corrected, but '// &
                  'this should not happen. Please, notify of '// &
                  'the issue providing your inputs'
                call abortedS(umsg,urou,.False.,.True.)

              end if ! Can issue warning

              ! Correct
              MRAMc = lMRAMc

            end if ! Different
          end if ! lMRAMc data

          !
          ! Solve NLTE problem
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
          call free_LTElines_full(LTElines)

          ! Copy coordinates
          int_buff = icoords

          ! Convert ix coordinate into node coordinate
          int_buff(1) = &
              (int_buff(1) - 2 + Input%sol_box(1))*dims(2) + &
              int_buff(2) - 1 + Input%sol_box(3)

          ! Write in cache
          call write_cache(unitC,Input%cache,int_buff,check)

          ! Check if could write in cache file
          aborting = .not.check

          ! Update NLOS received
          NLOSr = NLOSr + 1

          ! If we went beyond the number of LOS, exit
          if (NLOSr.ge.NLOS) exit

        end do ! While there is work to do

      end if ! Distribute pixels or single group

      !
      ! Close files
      !

      ! Global master
      if (gpid.eq.0) then

        ! Close model atmosphere
        call close_file(unitA)

        ! If JKQ data, close file
        if (sizeJ.gt.0) &
          call close_file(unitJ)

      end if ! Global master

      ! Clean memory
      if (allocated(buffer_atmo)) then
        MRAMc = MRAMc - 1d-6*sizeof(buffer_atmo)
        deallocate(buffer_atmo)
      end if
      if (allocated(buffer_JKQ)) then
        MRAMc = MRAMc - 1d-6*sizeof(buffer_JKQ)
        deallocate(buffer_JKQ)
      end if

      end subroutine HanleRT15DS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage the solution in a model atmosphere in the inversion\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!    Atomb(Atom_class(:)): Structures with atomic data for
      !!                          background atoms\n
      !!       Mol(Mol_class(:)): Structures with molecular data\n
      !!   GeomI(Geometry_class): Structure with geometric data for
      !!                          the intensity problem\n
      !!    Geom(Geometry_class): Structure with geometric data\n
      !!      Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                          and J-symbols\n
      !!   Frec(Frequency_class): Structure with frequency data\n
      !!      fudge(fudge_class): Structure with fudge data\n
      !!    kurucz(kurucz_class): Structure with Kurucz line data\n
      !!         MPID(MPI_class): Structure with MPI data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Bfield(Bfield_class): Structure with magnetic field data\n
      !!      Input(Input_class): Structure with configuration data\n
      !!     Sol(Solution_class): Structure with the frequency and
      !!                          synthetic Stokes parameters in the
      !!                          frequency range of the inverted
      !!                          data\n
      !!  SolF(Solution_F_class): Structure with the solution of
      !!                          the self-consistent problem and
      !!                          the corresponding emergent
      !!                          profiles, contribution function,
      !!                          and height for optical depth
      !!                          equal to one\n
      !!             RF(logical): If solving the NLTE problem to
      !!                          calculate the response function
      subroutine HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                            fudge,kurucz,MPID,Atmo,Bfield,Input, &
                            Sol,SolF,RF)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Input_class), intent(inout):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(Geometry_class), intent(inout):: GeomI, Geom
      type(Frequency_class), intent(inout):: Frec
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(MPI_class), intent(inout):: MPID
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: RF

      ! Local

      type(LTEline_class), dimension(:), allocatable:: LTElines

      logical:: PRAM_local, IRAM_local,tau1_local, contr_local
      logical:: lio,lie,lp,lpe,lporlpe,lload

      integer:: ia,nfreq_local


      ! Store local copies of global variables
      PRAM_Local = PRAM
      IRAM_Local = IRAM
      nfreq_local = nfreq

      ! Be sure that recurrent atmospheric quantities are deallocated
      if (allocated(Atmo%zalt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%zalt)
        deallocate(Atmo%zalt)
      end if
      if (associated(Atmo%Bx)) nullify(Atmo%Bz,Atmo%Bx,Atmo%By)

#ifdef DEBUGATMO
      ! RT master write model atmosphere in ASCII
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,0)
#endif

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
      ! Prepare Input and prepare the buffers to store solution
      call prepare_buffers(SolF,Input,Atom,GeomI,Geom,RF)

      !
      ! Define the flags to route the code
      !

      ! If loading from the solution
      lload = Input%mode.eq.'R'.or.Input%mode.eq.'B'

      ! If iterating in intensity
      lio = (Input%mode.eq.'W'.and.Input%force.ne.'P').or. &
            (Input%mode.eq.'B'.and. &
             (Input%force.eq.'I'.or.Input%force.eq.'A'))

      ! If emergent intensity
      lie = Input%force.eq.'I'.and.Geom%nThLOS.GT.0

      ! It iterating in polarization
      lp = (Input%mode.eq.'W'.and.Input%force.ne.'I').or. &
           (Input%mode.eq.'B'.and. &
            (Input%force.eq.'N'.or.Input%force.eq.'P'.or. &
             Input%force.eq.'A'))

      ! If emergent polarization
      lpe = Input%force.ne.'I'.and.Geom%nThLOS.gt.0

      ! If doing anything in polarization
      lporlpe = lp.or.lpe

      ! If calculating response function
      ! (that means we do not need tau and contribution)
      if (RF) then

        ! Copy to local
        tau1_local = Input%out_tau1
        contr_local = Input%out_contr

        ! Force to false
        Input%out_tau1 = .False.
        Input%out_contr = .False.

      end if


      !
      ! Solve the NLTE problem
      call hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input,GeomI,Geom, &
                 Bfield,Frec,Flgsg,fudge,kurucz,Atmo%JKQin,SolF, &
                 lload,lio,lie,lp,lpe,.False.)
      if (laborted) goto 1000


      !
      ! Apply PSF and select only the relevant frequencies
      call Profiles_Out(Frec,Sol,SolF%e_Stk,lpe, &
                        Input%lim_stk,Input%lim_fwhm)


      !
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

      !
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

      ! If active projection
1000  if (Sol%Projection) then

        ! For each height
        do ia=1,Atmo%nZ

          ! If inclination larger than 3rad
          if (Bfield%Btheta(ia).gt.3.0) then

            ! The field was just negative polarity
            Bfield%Bstrength(ia) = -Bfield%Bstrength(ia)
            Bfield%Btheta(ia) = 0d0

          end if ! Inclination

        end do ! Heights

      end if ! Active projection

      ! Restore global variables
1100  nfreq = nfreq_local
      PRAM = PRAM_Local
      IRAM = IRAM_Local

      ! If doing response function
      if (RF) then

        ! Restore variables
        Input%out_tau1 = tau1_local
        Input%out_contr = contr_local

      end if

      ! Free memory
      call free_mol(Mol)
      call free_damp(Atom)
      if (nAb.gt.0) call free_damp(Atomb)
      call free_gpop(Atom,Atomb,Mol)
      call free_lpop(Atom,Atomb)
      call free_cols(Atom)
      call free_LTElines_full(LTElines)

      end subroutine HanleRTTIC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage the solution pixel by pixel of a 3D model atmosphere
      !! or a list of slabs for the Coronal Line Emission problem\n
      !!    Input(Input_class): Structure with configuration data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine HanleCLE(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), allocatable, intent(inout):: Mol
      type(Input_class), intent(inout):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(inout):: fudge
      type(MPI_class), intent(inout):: MPID

      ! Local

      type(Atmo_class):: Atmo
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Geometry_class):: Geom,Gdummy
      type(kurucz_class):: kurucz
      type(chianti_class):: chianti
      type(spect_class):: spect

      logical:: aborting,check,receiving,lcache
      logical:: lie,lpe,lload,double,lDK
      logical, dimension(:,:), allocatable:: cache

      integer:: unitA,unitC,NIV,NIC
      integer:: iy1,iz1,sizeA,mode,norm,ioffset
      integer:: ia,iy0,iz0,iy,iz,inod,NLOSr,NLOS,ip,iproc
      integer, dimension(3):: dims,int_buff
      integer, dimension(nA):: ion_value_ind
      integer, dimension(nA):: ion_column_ind
      integer, dimension(:), allocatable:: cpu_free
      integer, dimension(:), allocatable:: unitI

      double precision:: maxB,DwTa,ycoor,zcoor
      double precision, dimension(nA):: ion_value
      double precision, dimension(:), allocatable, target:: buffer
      double precision, dimension(:), allocatable:: x,y,z,dummy

      ! Pointers
      double precision, dimension(:), pointer:: buffer_atmo
      double precision, dimension(:), pointer:: buffer_ion


      ! Count non-allocatable memory in local structures
      MRAMc = MRAMc + 1d-6*(sizeof(Atmo) + &
                            sizeof(Frec) + &
                            sizeof(Red) + &
                            sizeof(Geom) + &
                            sizeof(Gdummy) + &
                            sizeof(kurucz) + &
                            sizeof(chianti) + &
                            sizeof(spect))

      ! Nullify pointers
      nullify(buffer_atmo,buffer_ion)

      ! If forced to be unmagnetized
      if (Input%unmagnetized) then

        ! Say maximum B is zero
        maxB = 0d0

      ! Not forcing unmagnetized
      else

        ! Assume there will be magnetic field
        maxB = 1d0 + TINYB

      end if ! Force unmagnetized

      ! If (PRD and T limits unknown) or T limits are unknown and
      ! needed for the Doppler width to transform frequencies, or
      ! if not forcing static but the maximum velocity is unknown
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

      ! If forcing static
      if (Input%static) then

        ! Force the maximum velocity to zero
        Input%maxV = 0d0

      ! Not forcing static
      else

        ! Get maximum velocity in synthesis units
        Input%maxV = Input%maxV*1d-6/c

      end if ! Forcing static

      !
      ! Define the flags to route the code
      !

      ! Never reading previous solution
      lload = .False.

      ! If doing intensity emergence
      lie = Input%force.eq.'I'

      ! If doing polarization emergence
      lpe = Input%force.ne.'I'


      !
      ! Index atom levels, sublevels, transitions, and
      ! magnetic components
      call set_atom_indexes(Atom,Input,lie,lpe, &
                            .True.,maxB.gt.0d0)

      !
      ! Set angular quadrature
      call gauss(Input,Gdummy,Geom,2,.True.,.False.,Flgsg)

      ! Verbose angular quadrature defined
      if(pid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbose
      end if


      !
      ! Define the output frequency axis
      call omegabuild(Frec,Atom,Input,maxB,dummy)

      ! Verbose that the frequency axis is defined
      if (gpid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbose
      end if


      !
      ! Split in groups of tasks to distribute pixels
      call setmpi15D(MPID,Input)

      ! Verbose distribution of pixel tasks
      if (gpid.eq.0) then
        umsg = ' - Tasks distributed'
        call verbose
      end if


      !
      ! Decide if need to keep Stokes
      !

      ! If not splitting in pixels or it is an RT master
      if (.not.MPID%mpi15d.or.(pid.eq.0.and.gpid.ne.0)) then

        ! Keep Stokes if already flagged or if doind PRD either
        ! dynamic or angle-dependent
        KSTK = KSTK.or.(PRD.and.(.not.Input%static.or..not.AV))

      ! If RT slave
      else

        ! Keep Stokes if doind PRD either dynamic or angle-dependent
        KSTK = PRD.and.(.not.Input%static.or..not.AV)

      end if ! RT slave


      !
      ! Global master remove fudge if splitting pixels
      if (gpid.eq.0.and.allocated(fudge%fudge_v).and.MPID%mpi15d) then
        MRAMc = MRAMc - 1d-6*sizeof(fudge%fudge_v)
        deallocate(fudge%fudge_v)
      end if


      !
      ! Pre-process
      !

      ! Initialize looping identifier
      iy1 = -1
      iz1 = -1

      ! Global
      if (gpid.eq.0) then

        ! Open files
        call open_atm_and_cache(Input,2,unitA,unitC,aborting,dims, &
                                mode,double,norm,cache,lcache)

        ! The master does not need the background atoms or
        ! molecules
        if (MPID%mpi15d) then
          call free_atom_full(Atomb)
          call free_mol_full(Mol)
        end if

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

      ! Initialize ionization quantities
      ion_value = -1d0
      ion_value_ind = -1
      ion_column_ind = -1
      NIC = 0
      NIV = 0

      ! For each atom
      do ia=1,nA

        ! Check if is explicitly none
        if (Input%ionf(ia)%typ.lt.0) then

          ! Go to next atom
          cycle

        ! If a number was provided
        else if (Input%ionf(ia)%typ.eq.1) then

          ! Add count to number
          NIV = NIV + 1

          ! And store index and value
          ion_value_ind(ia) = NIV
          ion_value(NIV) = Input%ionf(ia)%val

        ! If a file was provided
        else if (Input%ionf(ia)%typ.eq.0) then

          ! Add count to files
          NIC = NIC + 1

          ! And store index
          ion_column_ind(ia) = NIC

        ! What?!
        else

          ! Error
          umsg = 'Unexpected Atom_ion input'
          call verbose
          laborted = .True.

        end if ! Type of ionization fraction input

      end do ! Active atoms

      !
      ! If ionization fraction files
      if (NIC.gt.0) then

        ! Allocate units
        allocate(unitI(NIC))

        ! Run over atoms
        do ia=1,nA

          ! Get rolling index for current atom
          ip = ion_column_ind(ia)

          ! If not file, skip
          if (ip.lt.1) cycle

          ! If first instance
          if (ip.eq.1) then

            ! First unit to be defined
            unitI(ip) = 10001

          ! Not first instance
          else

            ! Get new unit
            unitI(ip) = unitI(ip-1) + 1

          end if ! First instance

          ! Open ionization fraction file
          call open_file(unitI(ip), Input%ionf(ia)%str, &
                         0, .False., check)

        end do ! Atoms

      end if ! There are ionization fraction files

      ! Check if aborting
      call gcontrol


      !
      ! Initialize buffer sizes
      call set_io_CLE_buffers(Input,mode,dims,Frec)

      ! Check if aborting
      call gcontrol

      ! Global master
      if (gpid.eq.0) then

        ! Check files exist
        if (lcache) call check_io_CLE_buffers_exist(Input,aborting)
        laborted = aborting

      end if

      ! Check if aborting
      call gcontrol

      !
      ! Input spectrum read
      call rinspect(Input,Atom,Geom,spect,Frec%omega)

      ! Check if aborting
      call gcontrol

      !
      ! Organize the tasks splitting (within groups)
      !

      ! If not splitting pixels or not the distributer
      if (.not.MPID%mpi15d.or.gpid.gt.0)  then

        ! Distribute tasks
        call setmpi_CLE(MPID,Input,Frec%IW_freq_in,Frec%IW_freq)

        ! Free IW_freq
        deallocate(Frec%IW_freq_in,Frec%IW_freq)

        !
        ! Finish preparingthe atom
        call prepareatomMPI(Atom)
        call omegainitmaster(Atom)

      end if

      ! Verbose MPI RT distributed
      if (gpid.eq.0) then
        umsg = ' - Tasks within groups distributed'
        call verbose
      end if

      ! Control
      call gcontrol


      !
      ! Photoionization quantites
      !

      ! If not splitting pixels or not the distributer
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! For every active atom
        do ia=1,nA

          ! Setup photoionization cross sections
          call setphoto(Atom(ia),Frec%omega)

        end do ! Active atoms

      ! The distributer
      else

        ! Just call control to match the call in setphoto
        call gcontrol

      end if ! Not splitting pixels or not the distributer

      ! Verbose photoionization quantities initialized
      if (gpid.eq.0) then
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbose
      end if


      !
      ! MPI sizes
      !

      ! Not global Master or not distributer
      if (.not.MPID%mpi15d.or.gpid.gt.0)  then

        ! Manage output frequency axis
        call refitfrec(Input,Frec,Atom,MPID)

        ! If photoionization frequencies
        if (Frec%Mpif0(pid).le.Frec%Mpif1(pid)) then

          ! Allocate space for exu in Frec
          allocate(Frec%exu(Frec%Mpif0(pid):Frec%Mpif1(pid),1))
          PRAMc = PRAMc + 1d-6*sizeof(Frec%exu)

        ! No photoionization frequencies
        else

          ! Nullify pointer
          nullify(Frec%exu)

        end if ! Need to allocate photoionizations
      end if ! Not global master or not distributer

      ! Check if aborting
      call gcontrol


      !
      ! Grid sharing (for cartesian or slab)
      !

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
        MRAMc = MRAMc + 1d-6*sizeof(buffer)

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

        ! Slab
        else

          ! Initialize for valgrind
          allocate(x(1),y(1),z(1))
          MRAMc = MRAMc + 1d-6*sizeof(x)
          MRAMc = MRAMc + 1d-6*sizeof(y)
          MRAMc = MRAMc + 1d-6*sizeof(z)
          x = 0d0
          y = 0d0
          z = 0d0

        end if ! Cartesian/Slab
      end if ! Cartesian or slab)


      !
      ! Read partition function data and abundances
      call rParfunAbund(Input,Atmo)

      ! Check if aborting
      call gcontrol


      !
      ! Read barklem data
      if (.not.MPID%mpi15d.or.gpid.gt.0) &
        call rBarklem(Input,Atom,Atomb)

      ! Check if aborting
      call gcontrol

      !
      ! Get Kurucz lines
      !

      ! If not distributing pixels and not the RT master doing MPI,
      ! or if not the distributer nor the RT master doing MPI
      if ((.not.MPID%mpi15d.and.(pid.gt.0.or..not.MPID%mpi)).or. &
          (gpid.gt.0.and.(pid.gt.0.or..not.MPID%mpi))) then

        ! If there is Kurucz data to read
        if (Input%NK.ge.1) then

          ! If Doppler width for transformation is the maximum
          if(Input%dws.eq.'MAX')then

            ! Flag and compute maximum velocity (sqrt)
            lDK = .True.
            DwTa = sqrt(Input%maxT)

          ! If Doppler width for transformation is the minimum
          else if(Input%dws.eq.'MIN')then

            ! Flag and compute minimum velocity (sqrt)
            lDK = .True.
            DwTa = sqrt(Input%minT)

          ! If Doppler width for transformation is a number
          else if(Input%dws.eq.'NUM')then

            ! Flag and compute factor
            lDK = .False.
            DwTa = Input%dw*1d-9/c

          end if ! Type of Doppler width transformation

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,lDK,kurucz)

        ! If there are not Kurucz data
        else

          ! Set number of transitions in Kurucz structure
          ! to zero
          kurucz%ntran = 0

        end if ! There is Kurucz data
      end if ! Could need Kurucz data

      ! Check if aborting
      call gcontrol


      !
      ! CHIANTI read
      call rCHIANTI(Input,chianti)

      ! Check if aborting
      call gcontrol

      !
      ! Properly start now
      !

      ! Global master anounce start
      if (gpid.eq.0) then
        umsg = ' - Starting to work on your model, this may '// &
               'take a while'
        call verbose
      end if

      !
      ! MPI version
      !

      ! Distributing tasks
      if (MPID%mpi15d) then

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Master
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (gpid.eq.0) then

          ! Open files so slaves can write later
          if (.not.lcache) call create_io_CLE_files(Input,mode,y,z, &
                                                    dims,Frec)

          ! Deallocate x,y,z
          if (allocated(x)) then
            MRAMc = MRAMc - 1d-6*sizeof(x)
            MRAMc = MRAMc - 1d-6*sizeof(y)
            MRAMc = MRAMc - 1d-6*sizeof(z)
            deallocate(x,y,z)
          end if

          ! Allocate cpu_free with group status
          allocate(cpu_free(MPID%ngroup))
          MRAMc = MRAMc + 1d-6*sizeof(cpu_free)
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

              !
              ! Send stop signal to everyone
              !

              !
              ! First receive whatever is pending

              ! For every CPU still sending
              do while (minval(cpu_free).lt..5d0)

                ! For every leader
                do ip=1,MPID%ngroup

                  ! If already free, skip
                  if (cpu_free(ip).gt..5d0) cycle

                  ! Test if slave in the group is sending something
                  call MPI_IPROBE(MPID%ltslave(ip), &
                                  4+MPID%ltslave(ip), &
                                  MPI_COMM_WORLD, receiving, &
                                  MPI_STATUS_IGNORE, ierr)

                  ! If slave is calling
                  if (receiving) then

                    ! Receive the ping
                    call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                  MPID%ltslave(ip), &
                                  4+MPID%ltslave(ip), &
                                  MPI_COMM_WORLD, &
                                  MPI_STATUS_IGNORE, ierr)

                    ! Free the group
                    cpu_free(ip) = 1

                  end if ! Receiving from a CPU

                end do ! Receive from everyone
              end do ! While there is someone working

              ! And break the work loop
              exit

            end if ! Aborting

            ! If there are LOS to do and at least one free CPU
            if (inod.lt.NLOS.and.maxval(cpu_free).gt..5d0) then

              ! Save last values for pixel position
              iy0 = iy
              iz0 = iz

              ! Advance one pixel (Y slow, Z fast)
              iz = iz + 1
              if (iz.gt.dims(3)) then
                iy = iy + 1
                iz = 1
                if (iy.gt.dims(2)) cycle
              end if

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
                    MRAMc = MRAMc + 1d-6*sizeof(buffer)
                  else if (size(buffer).ne.sizeA) then
                    MRAMc = MRAMc - 1d-6*sizeof(buffer)
                    deallocate(buffer)
                    allocate(buffer(sizeA))
                    MRAMc = MRAMc + 1d-6*sizeof(buffer)
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

                  ! Setup-icoords
                  icoords = (/ iy , iz , inod+1 /)

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

              ! Advance LOS
              inod = inod + 1

              ! If done in cache, skip
              if (lcache) then
                if (cache(iz,iy)) then
                  NLOSr = NLOSr + 1
                  cycle
                end if
              end if

              ! Take a free cpu
              ip = maxloc(cpu_free, 1)

              ! Send signal to node
              call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                            MPID%ltslave(ip), &
                            2+MPID%ltslave(ip), &
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
                            3+MPID%ltslave(ip), &
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

            !
            ! Receive data
            !

            ! Initialize with first group
            ip = 1

            ! For every slave group
            do while (.True.)

              ! Test if slave in the group is sending something
              call MPI_IPROBE(MPID%ltslave(ip), &
                              4+MPID%ltslave(ip), &
                              MPI_COMM_WORLD, receiving, &
                              MPI_STATUS_IGNORE, ierr)

              ! If message failed, try again
              if (ierr.ne.0) cycle

              ! If slave is calling
              if (receiving) then

                ! Try to receive
                do while (.True.)

                  ! Receive the ping
                  call MPI_RECV(int_buff(1), 3, MPI_INTEGER, &
                                MPID%ltslave(ip), &
                                4+MPID%ltslave(ip), &
                                MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! If failed, try again
                  if (ierr.ne.0) cycle

                  ! Success, can leave
                  exit

                end do ! Try receiving

                ! Convert iy coordinate into node coordinate
                int_buff(1) = (int_buff(1)-1)*dims(3) + int_buff(2)

                ! Write in cache
                call write_cache(unitC,Input%cache,int_buff,check)

                ! Check if could write
                aborting = .not.check

                ! Update NLOS received
                NLOSr = NLOSr + 1

                ! Free the group
                cpu_free(ip) = 1

              end if ! Receiving from a CPU group

              ! Advance group
              ip = ip + 1

              ! If checked every group, leave
              if (ip.gt.MPID%ngroup) exit

            end do ! Slaves

            ! If we went beyond the number of LOS, exit
            if (NLOSr.ge.NLOS) exit

          end do ! While there is work to do

          !
          ! Notify finished with work
          !

          ! Finished signal
          icoords(1) = -1

          ! First CPU
          iproc = 1

          ! Try until sucess
          do while (.True.)

            ! send termination signal
            call MPI_SEND(icoords(1), 3, MPI_INTEGER, &
                          MPID%ltslave(iproc), &
                          2+MPID%ltslave(iproc), &
                          MPI_COMM_WORLD, ierr)

            ! If it fails try again
            if (ierr.ne.0) cycle

            ! Advance group
            iproc = iproc + 1

            ! If sent to all groups, leave
            if (iproc.gt.MPID%ngroup) exit

          end do ! slaves

          ! Free cpu_free
          MRAMc = MRAMc - 1d-6*sizeof(cpu_free)
          deallocate(cpu_free)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! Slaves
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else

          ! Initialize
          aborting = .False.

          ! Initialize atmosphere RAM
          call rAtmo_cle_prep(Input,Atmo,mode,norm)

          ! Normalize
          call normalize_cle(Atom,Red)

          ! Work until further notice
          do while (.True.)

            ! If RT leader
            if (pid.eq.0) then

              ! Try receiving until success
              do while (.True.)

                ! Wait for signal
                call MPI_RECV(icoords(1), 3, MPI_INTEGER, 0, &
                              2+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails, try again
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
                      MRAMc = MRAMc - 1d-6*sizeof(buffer)
                      deallocate(buffer)
                      allocate(buffer(sizeA))
                      MRAMc = MRAMc + 1d-6*sizeof(buffer)
                    end if
                  else
                    allocate(buffer(sizeA))
                    MRAMc = MRAMc + 1d-6*sizeof(buffer)
                  end if

                end if ! Only non-cartesian grid

                ! Receive LOS
                call MPI_RECV(buffer(1), sizeA, &
                              MPI_DOUBLE_PRECISION, &
                              0, 3+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! Success, leave loop
                exit

              end do ! Try receiving till success

            end if ! RT leader

            ! If liutenant has friends
            if (nproc.gt.1) then

              ! Try until done
              do while (.True.)

                ! Broadcast
                call MPI_BCAST(icoords(1), 3, MPI_INTEGER, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! If aborting
                if (icoords(1).lt.1) then

                  ! Everyone leave
                  aborting = .True.
                  exit

                end if ! Aborting

                ! If not-cartesian and a slave
                if (mode.eq.2.and.pid.ne.0) then

                  ! Allocate buffer
                  dims(1) = icoords(3)
                  ioffset = dims(1)*23 + 2
                  sizeA = ioffset + dims(1)*NIC
                  if (allocated(buffer)) then
                    if (size(buffer).lt.sizeA) then
                      MRAMc = MRAMc - 1d-6*sizeof(buffer)
                      deallocate(buffer)
                      allocate(buffer(sizeA))
                      MRAMc = MRAMc + 1d-6*sizeof(buffer)
                    end if
                  else
                    allocate(buffer(sizeA))
                    MRAMc = MRAMc + 1d-6*sizeof(buffer)
                  end if

                end if ! Only non-cartesian grid

                ! Broadcast model atmosphere for pixel
                call MPI_BCAST(buffer(1), sizeA, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails, try again
                if (ierr.ne.0) cycle

                ! Success, lave
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

            ! Clean state
            laborted = .False.

            !
            ! Call solver with correct arguments
            !

            ! If cartesian
            if (mode.eq.0) then

              ! Point to atmosphere part
              buffer_atmo => buffer(1:ioffset)

              ! Solve columns
              call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                       Flgsg,fudge,kurucz,buffer_atmo,x, &
                       y(icoords(1)),z(icoords(2)), &
                       dims,buffer_ion,ion_column_ind, &
                       ion_value_ind,ion_value,spect,chianti)

            ! If slab
            else if (mode.eq.1) then

              ! Point to atmosphere part
              buffer_atmo => buffer(1:ioffset)

              ! Solve columns
              call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                       Flgsg,fudge,kurucz,buffer_atmo,x,0d0,0d0, &
                       dims,buffer_ion,ion_column_ind, &
                       ion_value_ind,ion_value,spect,chianti)

            ! If non-cartesian
            else if (mode.eq.2) then

              ! Get true inode index
              icoords(3) = icoords(1) + (icoords(2)-1)*dims(2)

              ! Point to atmosphere part
              buffer_atmo => buffer(3:ioffset)

              ! Solve columns
              call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                       Flgsg,fudge,kurucz,buffer_atmo,x,buffer(1), &
                       buffer(2),dims,buffer_ion,ion_column_ind, &
                       ion_value_ind,ion_value,spect,chianti)

            end if

            ! Nullify pointer
            nullify(buffer_ion,buffer_atmo)

            ! If liutenant, send message to grand master
            if (pid.eq.0) then

              ! Failed
              if (laborted) then

                ! Send issue signal
                icoords(3) = -1
              ! Success
              else

                ! Send success signal
                icoords(3) = gpid

              end if ! Failure or success

              ! Try until achieved
              do while (.True.)

                ! Send status of pixel calculation
                call MPI_SEND(icoords(1),3,MPI_INTEGER,0, &
                              4+gpid,MPI_COMM_WORLD,ierr)

                ! If failed, try again
                if (ierr.ne.0) cycle

                ! Success, leave
                exit

              end do ! Try until success

            end if ! Send info back to grand master

          end do ! Work till finished

        end if ! Master or slave

      !
      ! Serial
      !
      else

        ! Open files so slaves can write later
        if (.not.lcache) call create_io_CLE_files(Input,mode,y,z, &
                                                  dims,Frec)

        ! Deallocate x,y,z
        if (allocated(x)) then
          MRAMc = MRAMc - 1d-6*sizeof(x)
          MRAMc = MRAMc - 1d-6*sizeof(y)
          MRAMc = MRAMc - 1d-6*sizeof(z)
          deallocate(x,y,z)
        end if

        ! Initialize indexes and sizes
        iy = 1
        iz = 0
        inod = 0
        NLOSr = 0
        NLOS = dims(2)*dims(3)

        ! Initialize
        aborting = .False.

        ! Initialize atmosphere RAM
        call rAtmo_cle_prep(Input,Atmo,mode,norm)

        ! Work until exhausted
        do while (.True.)

          ! If aborting, send stop signal to everyone
          if (aborting) then

            ! Abort
            call aborted_silent

          end if ! Aborting

          ! Save last values for pixel position
          iy0 = iy
          iz0 = iz

          ! Advance one pixel (Y slow, Z fast)
          iz = iz + 1
          if (iz.gt.dims(3)) then
            iy = iy + 1
            iz = 1
            if (iy.gt.dims(2)) cycle
          end if

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
                MRAMc = MRAMc + 1d-6*sizeof(buffer)
              else if (size(buffer).ne.sizeA) then
                MRAMc = MRAMc - 1d-6*sizeof(buffer)
                deallocate(buffer)
                allocate(buffer(sizeA))
                MRAMc = MRAMc + 1d-6*sizeof(buffer)
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

              ! Setup-icoords
              icoords = (/ iy , iz , inod+1 /)

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

          ! Advance LOS
          inod = inod + 1

          ! If done in cache, skip
          if (lcache) then
            if (cache(iz,iy)) then
              NLOSr = NLOSr + 1
              cycle
            end if
          end if

          ! Clean state
          laborted = .False.

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

            ! Point to atmosphere part
            buffer_atmo => buffer(1:ioffset)

            ! Solve columns
            call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                     Flgsg,fudge,kurucz,buffer_atmo,x, &
                     y(icoords(1)),z(icoords(2)), &
                     dims,buffer_ion,ion_column_ind, &
                     ion_value_ind,ion_value,spect,chianti)

          ! If slab
          else if (mode.eq.1) then

            ! Point to atmosphere part
            buffer_atmo => buffer(1:ioffset)

            ! Solve columns
            call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                     Flgsg,fudge,kurucz,buffer_atmo,x,0d0,0d0, &
                     dims,buffer_ion,ion_column_ind, &
                     ion_value_ind,ion_value,spect,chianti)

          ! If non-cartesian
          else if (mode.eq.2) then

            ! Get true inode index
            icoords(3) = icoords(1) + (icoords(2)-1)*dims(2)

            ! Point to atmosphere part
            buffer_atmo => buffer(3:ioffset)

            ! Solve columns
            call CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                     Flgsg,fudge,kurucz,buffer_atmo,x,buffer(1), &
                     buffer(2),dims,buffer_ion,ion_column_ind, &
                     ion_value_ind,ion_value,spect,chianti)

          end if

          ! Nullify pointer
          nullify(buffer_ion,buffer_atmo)

          ! Convert iy coordinate into node coordinate
          int_buff(1) = (int_buff(1)-1)*dims(3) + int_buff(2)

          ! Write in cache
          call write_cache(unitC,Input%cache,int_buff,check)

          ! Check if could write
          aborting = .not.check

          ! Update NLOS received
          NLOSr = NLOSr + 1

          ! If we went beyond the number of LOS, exit
          if (NLOSr.ge.NLOS) exit

        end do ! While there is work to do

      end if ! Distribute pixels or single group

      !
      ! Close files
      !
      if (gpid.eq.0) then
        call close_file(unitA)
        do ip=1,NIC
          call close_file(unitI(ip))
        end do
      end if

      ! Free memory
      if (allocated(buffer)) then
        MRAMc = MRAMc - 1d-6*sizeof(buffer)
        deallocate(buffer)
      end if
      if (allocated(x)) then
        MRAMc = MRAMc - 1d-6*sizeof(x)
        deallocate(x)
      end if
      if (allocated(y)) then
        MRAMc = MRAMc - 1d-6*sizeof(y)
        deallocate(y)
      end if
      if (allocated(z)) then
        MRAMc = MRAMc - 1d-6*sizeof(z)
        deallocate(z)
      end if

      end subroutine HanleCLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module hanlert_mod
