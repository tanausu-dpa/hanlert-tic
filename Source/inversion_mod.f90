      !> Inversion module
      module inversion_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC/NSSCC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     22/02/2023
!  Last version:
!     15/05/2025 V4.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/05/2025:    V4.0.2 - Generalized declarations of Atom to
!                             allow for empty arrays for any of
!                             them (TdPA)
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
!  Inversion
!    Manage the inversion of a given set of Stokes profiles
!
!  Fit
!    Manage the call to the minimization algorithm
!
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use bounds_mod
      use commons_mod
      use initinv_mod
      use lmfit_mod
      use model_mod
      use parameters_mod
      use regul_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage the inversion of a given set of Stokes profiles\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!       Atmo_in(Atmo_class): Structure with atmospheric data
      !!                            read from the input\n
      !!   Bfield_in(Bfield_class): Structure with magnetic field data
      !!                            read from the input\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!            imask(integer): Indicate if this pixel is masked
      !!                            when restarting the inversion
      subroutine Inversion(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec, &
                           fudge,kurucz,MPID,Atmo_in,Bfield_in, &
                           Input,Inf_Stokes,Inf_Nodes,Sol,imask)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Fctsg_class), intent(inout):: Flgsg
      type(Geometry_class), intent(inout):: GeomI, Geom
      type(Frequency_class), intent(inout):: Frec
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(MPI_class), intent(inout):: MPID
      type(Atmo_class), intent(in):: Atmo_in
      type(Bfield_class), intent(inout):: Bfield_in
      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      integer, intent(in):: imask

      ! Local

      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield,Bfield0

      integer:: i

      double precision, dimension(Input%nvar):: TMP_Weight


      ! Count memory in local structures
      MRAMc = MRAMc + 1d-6*sizeof(Atmo)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield0)

      ! Already broken, go back
      if (laborted) return

      ! If master
      if (pid.eq.0) then

        ! Print current pixel we are working on
        write(umsg, '(2(A,i5))') ''
        call verboseI(3)
        write(umsg, '(2(A,i5))') ''
        call verboseI(3)
        write(umsg, '(2(A,i5))') ' - Current pixel: ix = ', &
                                 icoords(1), &
                                 '; iy = ',  &
                                 icoords(2)
        call verboseI(0)
        if (vlevel.eq.0) call verboseI(3)

      end if ! Master

      ! Set up from common variable
      Atmo%nz = nz

      ! Store current regularization weights
      TMP_Weight = Inf_Nodes%Regul_weight

      !
      ! LOS geometry
      !

      ! Get LOS
      Inf_Nodes%mu = Inf_Stokes%mu
      Inf_Nodes%azimuth = Inf_Stokes%azimuth

      ! Set up geometry for LOS in synthesis structures
      GeomI%L_mu(1) = Inf_Stokes%mu
      GeomI%L_theta(1) = acos(GeomI%L_mu(1))
      GeomI%L_phi(1) = Inf_Stokes%azimuth
      if (Input%Type_inversion.ne.0) then
        Geom%L_mu(1) = Inf_Stokes%mu
        Geom%L_theta(1) = acos(Geom%L_mu(1))
        Geom%L_phi(1) = Inf_Stokes%azimuth
      end if

      !
      ! Initial model atmosphere
      !

      ! If asymmetry nodes, signal that they have dimension
      if (Inf_Nodes%Num_Asymmetry.gt.0) Input%nasym = 1

      ! If inverting from scratch and stratification from inputs
      if (trim(Input%Inv_init).eq.'INIT'.and. &
          Input%Atmo_Input.gt.0) then

        ! Generate a new stratification and interpolate the input
        ! model
        call Atmo_Stratify(Atmo_in,Atmo,Bfield_in,Bfield, &
                           Input%Atmo_strat_done)

      ! Restoring or copying stratification
      else

        ! Copy input atmosphere into actual model atmosphere
        call cAtmo(Atmo_in,Atmo)
        call cBfield(Bfield_in,Bfield)

      end if

      ! If Asymmetry not allocated
      if (.not.allocated(Atmo%JKQin)) then

        ! Allocate JKQin in model atmosphere
        allocate(Atmo%JKQin(8*nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%JKQin)

        ! And initialize
        Atmo%JKQin = 0d0

      end if ! Asymmetry not allocated

      ! If gas pressure from model atmosphere
      if (.not.Inf_Nodes%Pg_auto) then

        ! Boundary from model
        Inf_Nodes%Pg_bound = Atmo%Pg(1)

      end if ! Gas pressure from model atmosphere

      !
      ! Magnetic field and velocity vectors reference
      !

      ! If LOS field
      if (Inf_Nodes%Btype.eq.1) then

        ! Allocate magnetic field arrays for LOS
        allocate(Bfield%Bpos(Atmo%nZ),Bfield%Azimuth(Atmo%nZ))
        allocate(Bfield%Blos(Atmo%nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Bfield%Blos)
        MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bpos)
        MRAMc = MRAMc + 1d-6*sizeof(Bfield%Azimuth)

        ! If there is magnetic field already
        if (maxval(Bfield%Bstrength).gt.TINYB) then

          ! Convert to LOS
          call B2Blos(Atmo%nZ,GeomI%L_mu(1),GeomI%L_phi(1), &
                      Bfield%Bstrength, &
                      Bfield%Btheta, &
                      Bfield%Bphi, &
                      Bfield%Blos, &
                      Bfield%Bpos, &
                      Bfield%Azimuth)

        ! No field
        else

          ! Just set to zero
          Bfield%Blos = 0d0
          Bfield%Bpos = 0d0
          Bfield%Azimuth = 0d0

        end if ! Magnetic field exists
      end if ! LOS field

      ! If LOS velocity
      if (Inf_Nodes%vtype.eq.1) then

        ! Allocate velocity field arrays
        allocate(Atmo%vlos(Atmo%nZ),Atmo%vpos(Atmo%nZ))
        allocate(Atmo%vphi(Atmo%nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%vlos)
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%vpos)
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%vphi)

        ! If velocity
        if (dyn) then

          ! Convert to LOS
          call v2vlos(Atmo%nZ,GeomI%L_mu(1),GeomI%L_phi(1), &
                      Atmo%vx,Atmo%vy,Atmo%vz, &
                      Atmo%vlos,Atmo%vpos,Atmo%vphi)

        ! No velocity
        else

          ! Just zero
          Atmo%vlos = 0d0
          Atmo%vpos = 0d0
          Atmo%vphi = 0d0

        end if ! Dynamic
      end if ! LOS velocity

      !
      ! Set-up nodes
      !

      ! Initialize node positiones
      call Init_Nodes(Atmo,Input,Inf_Nodes)

      ! Initialize node values
      call Initialize_Nodes(Atmo,Bfield,Inf_Nodes)

      ! Set regularization constant for gas pressure
      if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg).and. &
          Inf_Nodes%hydroeq) &
        Inf_Nodes%Const(Inf_Nodes%index_Pg) = &
                             Inf_Nodes%Node(Inf_Nodes%index_Pg)%Var(1)

      ! Set regulatization constant for diffuse light
      Inf_Nodes%Const(Inf_Nodes%index_f) = Input%f_diff

      ! If retoring and not masking
      if (trim(Input%Inv_init).ne.'INIT'.and.imask.eq.0) then

        ! Magnetic field inclination index
        i = Inf_Nodes%index_Bt

        ! If inverting magnetic field inclination or BPOS
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If maximum value too small (1d-3)
          if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-3) then

            ! Master
            if (pid.eq.0) then

              ! If vertical coords.
              if (Inf_Nodes%Btype.eq.0) then

                ! Verbose
                write(umsg, '(A)') &
                  ' - Magnetic field inclination re-initialized: '
                call verboseI(3)

              ! If LOS coords.
              else

                ! Verbose
                write(umsg, '(A)') &
                  ' - Transversal magnetic field re-initialized: '
                call verboseI(3)

              end if ! Type of coordinates
            end if ! Master

            ! Get initial guess from Input structure
            call re_set_nodes(Inf_Nodes,i,Input%ini_Bpos)

          end if ! If too small inclination or Bpos
        end if ! If inversing Bpos

        ! Magnetic field azimuth index
        i = Inf_Nodes%index_Bp

        ! If inverting magnetic field azimuth
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If maximum value too small (1d-3)
          if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-3) then

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg, '(A)') &
                ' - Magnetic field azimuth re-initialized: '
              call verboseI(3)

            end if ! Master

            ! Get initial guess from Input structure
            call re_set_nodes(Inf_Nodes,i,Input%ini_Bazi)

          end if ! If too small inclination or Bpos
        end if ! If inversing Bpos

        ! Velocity x/pos
        i = Inf_Nodes%index_vx

        ! If inverting vx or vpos
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If maximum value too small (1d-3)
          if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-9/c) then

            ! Master
            if (pid.eq.0) then

              ! If vertical coords.
              if (Inf_Nodes%vtype.eq.0) then

                ! Verbose
                write(umsg, '(A)') &
                  ' - Velocity X component re-initialized: '
                call verboseI(3)

              ! If LOS coords.
              else

                ! Verbose
                write(umsg, '(A)') &
                  ' - Transversal velocity re-initialized: '
                call verboseI(3)

              end if ! Type of coordinates
            end if ! Master

            ! Get initial guess from Input structure
            call re_set_nodes(Inf_Nodes,i,Input%ini_vpos*1d-6/c)

          end if ! If too small inclination or vpos
        end if ! If inversing vpos

        ! Velocity y/azimuth
        i = Inf_Nodes%index_vy

        ! If inverting velocity field azimuth
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If azimuth
          if (Inf_Nodes%vtype.eq.1) then

            ! If maximum value too small (1d-3)
            if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-3) then

              ! Master
              if (pid.eq.0) then

                ! Verbose
                write(umsg, '(A)') &
                  ' - Velocity azimuth re-initialized: '
                call verboseI(3)

              end if ! Master

              ! Get initial guess from Input structure
              call re_set_nodes(Inf_Nodes,i,Input%ini_vazi)

            end if ! If too small vy

          ! If vy
          else

            ! If maximum value too small (1d-3)
            if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-9/c) then

              ! Master
              if (pid.eq.0) then

                ! Verbose
                write(umsg, '(A)') &
                  ' - Velocity Y component re-initialized: '
                call verboseI(3)

              end if ! Master

              ! Get initial guess from Input structure
              call re_set_nodes(Inf_Nodes,i,Input%ini_vazi*1d-6/c)

            end if ! If too small v azimuth
          end if ! vy or azimuth
        end if ! inverting vy or azimuth
      end if ! Type of restore file and mask

      ! If regularizing
      if (Inf_Nodes%Regul_Flag) then

        ! Initialize regulatization
        call Init_Regul(Inf_Nodes)
        if (laborted) goto 1000

      end if ! Regularizing

      ! Put node values within boundaries (magnetic field azimuth)
      i = Inf_Nodes%index_Bp
      call FoldBounds(Inf_Nodes%Node(i),Inf_Nodes%Num_Nodes(i))

      ! Put node values within boundaries (velocity azimuth)
      i = Inf_Nodes%index_vy
      if (Inf_Nodes%vtype.eq.1) &
        call FoldBounds(Inf_Nodes%Node(i),Inf_Nodes%Num_Nodes(i))

      ! Master
      if (pid.eq.0) then

        ! Verbose
        umsg = ' - Initialized model atmosphere/nodes'
        call verboseI(3)

      end if ! Master

      !
      ! Select depending on the type of inversion
      select case(Input%Type_Inversion)

        ! Only thermal
        case(0)

          ! If there are thermal nodes
          if (Inf_Nodes%Num_Thermal.gt.0) then

            ! Master
            if (pid.eq.0) then

              ! Verbose
              umsg = " - Fitting the thermal and dynamical parameters"
              call verboseI(3)

            end if ! Master

            ! Force inputs
            Inf_Nodes%Nodes_type = 0
            Input%force = 'I'

            ! If not masked, interpolate nodes into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                   Atomb,Mol,Input,fudge)

            ! Allocate mangetic arrays
            allocate(Bfield0%Bstrength(Atmo%nZ))
            allocate(Bfield0%Btheta(Atmo%nZ))
            allocate(Bfield0%Bphi(Atmo%nZ))
            allocate(Bfield0%Blos(Atmo%nZ))
            allocate(Bfield0%Bpos(Atmo%nZ))
            allocate(Bfield0%Azimuth(Atmo%nZ))
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bstrength)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Btheta)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bphi)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Blos)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bpos)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Azimuth)

            ! Initialize magnetic arrays
            Bfield0%Bstrength = 0d0
            Bfield0%Btheta = 0d0
            Bfield0%Bphi = 0d0
            Bfield0%Blos = 0d0
            Bfield0%Bpos = 0d0
            Bfield0%Azimuth = 0d0

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield0,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

            ! Deallocate magnetic arrays
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bstrength)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Btheta)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bphi)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Blos)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bpos)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Azimuth)
            deallocate(Bfield0%Bstrength,Bfield0%Btheta,Bfield0%Bphi)
            deallocate(Bfield0%Blos,Bfield0%Bpos,Bfield0%Azimuth)

          end if ! If there are thermal nodes

        ! Only magnetic
        case(1)

          ! If there are magnetic nodes
          if (Inf_Nodes%Num_Mag.gt.0) then

            ! Master
            if (pid.eq.0) then

              ! Verbose
              umsg = " - Fitting the of magnetic parameters"
              call verboseI(3)

            end if ! Master

            ! Force inputs
            Inf_Nodes%Nodes_type = 1
            Input%force = 'N'

            ! If not masked, interpolate nodes into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                   Atomb,Mol,Input,fudge)

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

          end if ! Magnetic nodes

        ! Fit all, together
        case(2)

          ! Master
          if (pid.eq.0) then

            ! Verbose
            umsg = " - Fitting the thermal, dynamical, "// &
                   "and magnetic parameters together"
            call verboseI(3)

          end if ! Master

          ! Force inputs
          Inf_Nodes%Nodes_type = 2
          Input%force = 'N'

          ! If not masked, interpolate nodes into the atmosphere
          if (imask.eq.0) &
            call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                 Atomb,Mol,Input,fudge)

          ! Fit the profiles
          call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                   kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                   Inf_Nodes,Sol,imask,.True.)

        ! Fit all, but not together
        case(3)

          ! If thermal nodes
          if (Inf_Nodes%Num_Thermal.gt.0) then

            ! Allocate magnetic field quantities
            allocate(Bfield0%Bstrength(Atmo%nZ))
            allocate(Bfield0%Btheta(Atmo%nZ))
            allocate(Bfield0%Bphi(Atmo%nZ))
            allocate(Bfield0%Blos(Atmo%nZ))
            allocate(Bfield0%Bpos(Atmo%nZ))
            allocate(Bfield0%Azimuth(Atmo%nZ))
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bstrength)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Btheta)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bphi)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Blos)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bpos)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Azimuth)

            ! Initialize magnetic field
            Bfield0%Bstrength = 0d0
            Bfield0%Btheta = 0d0
            Bfield0%Bphi = 0d0
            Bfield0%Blos = 0d0
            Bfield0%Bpos = 0d0
            Bfield0%Azimuth = 0d0

            ! Master
            if(pid.eq.0) then

              ! Verbose
              umsg = " - Fitting the thermal and dynamical parameters"
              call verboseI(3)

            end if ! Master

            ! Force inputs
            Inf_Nodes%Nodes_type = 0
            Input%force = 'I'

            ! If not masked, interpolate nodes into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                   Atomb,Mol,Input,fudge)

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield0,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.False.)

            ! Deallocate magnetic field quantities
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bstrength)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Btheta)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bphi)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Blos)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bpos)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Azimuth)
            deallocate(Bfield0%Bstrength,Bfield0%Btheta,Bfield0%Bphi)
            deallocate(Bfield0%Blos,Bfield0%Bpos,Bfield0%Azimuth)

          end if ! Thermal nodes

          ! Verbose merit function
          umsg = ' * '
          call verboseI(0)
          if (vlevel.eq.0) call verboseI(3)

          ! Verbose merit function
          write(umsg,'(A,i4,2(3x,A,es15.4))')  &
              ' * Thermal cycle ended, starting full cycle'

          ! Global master
          if (gpid.eq.0) then

            ! Verbose
            call verboseI(0)
            call verboseI(4)

          ! Slaves
          else

            ! Verbose
            call verboseI(0)
            if (vlevel.eq.0) call verboseI(3)

          end if ! Global master or other

          ! There are magnetic nodes
          if (Inf_Nodes%Num_Mag.gt.0.and..not.laborted) then

            ! Master
            if(pid.eq.0) then

              ! Verbose
              umsg = " - Fitting the thermal, dynamical, "// &
                     "and magnetic parameters together"
              call verboseI(1)

            end if ! Master

            ! Force inputs
            Inf_Nodes%Nodes_type = 2
            Input%force = 'N'

            ! If not masked, interpolate nodes into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                   Atomb,Mol,Input,fudge)

            ! Fit profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

          end if ! There are magnetic nodes

        ! Fit all, but not together, and second cycle is only magnetic
        case(4)

          ! If thermal nodes
          if (Inf_Nodes%Num_Thermal.gt.0) then

            ! Allocate magnetic field quantities
            allocate(Bfield0%Bstrength(Atmo%nZ))
            allocate(Bfield0%Btheta(Atmo%nZ))
            allocate(Bfield0%Bphi(Atmo%nZ))
            allocate(Bfield0%Blos(Atmo%nZ))
            allocate(Bfield0%Bpos(Atmo%nZ))
            allocate(Bfield0%Azimuth(Atmo%nZ))
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bstrength)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Btheta)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bphi)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Blos)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bpos)
            MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Azimuth)

            ! Initialize magnetic field
            Bfield0%Bstrength = 0d0
            Bfield0%Btheta = 0d0
            Bfield0%Bphi = 0d0
            Bfield0%Blos = 0d0
            Bfield0%Bpos = 0d0
            Bfield0%Azimuth = 0d0

            ! Master
            if(pid.eq.0) then

              ! Verbose
              umsg = " - Fitting the thermal and dynamical parameters"
              call verboseI(3)

            end if ! Master

            ! Force inputs
            Inf_Nodes%Nodes_type = 0
            Input%force = 'I'

            ! If not masked, interpolate nodes into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                   Atomb,Mol,Input,fudge)

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield0,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.False.)

            ! Deallocate magnetic field quantities
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bstrength)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Btheta)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bphi)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Blos)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Bpos)
            MRAMc = MRAMc - 1d-6*sizeof(Bfield0%Azimuth)
            deallocate(Bfield0%Bstrength,Bfield0%Btheta,Bfield0%Bphi)
            deallocate(Bfield0%Blos,Bfield0%Bpos,Bfield0%Azimuth)

          end if ! There are thermal nodes

          ! Verbose merit function
          umsg = ' * '
          call verboseI(0)
          if (vlevel.eq.0) call verboseI(3)

          ! Verbose merit function
          write(umsg,'(A,i4,2(3x,A,es15.4))')  &
              ' * Thermal cycle ended, starting only magnetic cycle'

          ! If global master
          if (gpid.eq.0) then

            ! Verbose
            call verboseI(0)
            call verboseI(4)

          ! Everyone else
          else

            ! Verbose
            call verboseI(0)
            if (vlevel.eq.0) call verboseI(3)

          end if ! Global master or other

          ! There are magnetic nodes
          if (Inf_Nodes%Num_Mag.gt.0.and..not.laborted) then

            ! Master
            if(pid.eq.0) then

              ! Verbose
              umsg = " - Fitting the thermal, dynamical, "// &
                     "and magnetic parameters together"
              call verboseI(1)

            end if ! Master

            ! Force inputs
            Inf_Nodes%Nodes_type = 1
            Input%force = 'N'

            ! If not masked, interpolate nodes into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom, &
                                   Atomb,Mol,Input,fudge)

            ! Fit profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

          end if ! There are magnetic nodes

      end select ! Type of inversion

      ! Restore regularization weights
1000  Inf_Nodes%Regul_weight = TMP_Weight

      ! Purge atmospheric model
      call free_Atmo(Atmo,.True.)
      call free_B(Bfield)

      ! Free inherent memory of local structures
      MRAMc = MRAMc - 1d-6*sizeof(Atmo)
      MRAMc = MRAMc - 1d-6*sizeof(Bfield)
      MRAMc = MRAMc - 1d-6*sizeof(Bfield0)

      ! If master
      if (pid.eq.0) then

        ! Message
        if (laborted) then
          umsg = " # Pixel inversion failed"
        else
          umsg = " - Pixel inversion finished"
        end if

        ! Global master
        if (gpid.eq.0) then

          ! Verbose
          call verboseI(0)
          call verboseI(4)

        ! Normal master
        else

          ! Verbose
          call verboseI(4)

        end if ! Global or local master
      end if ! Master

      ! Round
      MRAMc = 1d-6*nint(MRAMc*1d6)
      SRAMc = 1d-6*nint(SRAMc*1d6)

      ! Sanity
      if (abs(SRAMc).gt.0d0) then

        ! Warning
        if (Sol%warning) then

          ! Deflag
          Sol%warning = .False.

          ! Write message
          urou = 'Inversion'
          write(umsg,'(A,es13.6,A)') &
            'The Solution RAM counter is not equal to '// &
            'zero at the exit of the Inversion function ',SRAMc, &
            ' != 0. It is being corrected, but '// &
            'this should not happen. Please, notify of '// &
            'the issue providing your inputs'
          call abortedS(umsg,urou,.False.,.True.)

        end if ! Can issue warning

        ! Correct
        SRAMc = 0d0

      end if ! Sanify check

      return

      end subroutine Inversion

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage the call to the minimization algorithm\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!            imask(integer): Indicate if this pixel is masked
      !!                            when restarting the inversion\n
      !!           saving(logical): If the result is to be stored
      subroutine Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,saving)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Fctsg_class), intent(inout):: Flgsg
      type(Geometry_class), intent(inout):: GeomI, Geom
      type(Frequency_class), intent(inout):: Frec
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(MPI_class), intent(inout):: MPID
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      logical, intent(in):: saving
      integer, intent(in):: imask

      ! Local

      type(LMFIT_class):: LM_Stru


      ! Count inherent memory of local structure
      MRAMc = MRAMc + 1d-6*sizeof(LM_Stru)

      ! If thermal inversion
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Number of nodes to invert
        LM_Stru%Num = Inf_Nodes%Num_Thermal + Inf_Nodes%Num_glob

        ! If diffuse light
        if (Sol%Diff_flag) then

          ! First index is diffuse light
          Inf_Nodes%Indx_b = Inf_Nodes%index_f

        ! No diffuse light
        else

          ! First index is temperature
          Inf_Nodes%Indx_b = Inf_Nodes%index_T

        end if ! Diffuse light

        ! Last index is gas pressure
        Inf_Nodes%Indx_e = Inf_Nodes%index_Pg

      ! If magnetic inversion
      else if (Inf_Nodes%Nodes_Type.eq.1) then

        ! Number of nodes to invert
        LM_Stru%Num = Inf_Nodes%Num_Mag + Inf_Nodes%Num_Asymmetry + &
                      Inf_Nodes%Num_glob

        ! First index is B strength or LOS
        Inf_Nodes%Indx_b = Inf_Nodes%index_B

        ! If diffuse light
        if (Sol%Diff_flag) then

          ! Last index is diffuse light
          Inf_Nodes%Indx_e = Inf_Nodes%index_f

        ! No diffuse light
        else

          ! Last index is B azimuth
          Inf_Nodes%Indx_e = Inf_Nodes%index_Bp

        end if ! Diffuse light

      ! If full inversion
      else if (Inf_Nodes%Nodes_Type.eq.2) then

        ! Number of nodes to invert
        LM_Stru%Num = Inf_Nodes%Num_Fit

        ! All variables are included
        Inf_Nodes%Indx_b = 1
        Inf_Nodes%Indx_e = Input%nvar

      end if ! Inversion type

      ! Call LM fit routine
      call LMFIT(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge,kurucz, &
                 MPID,Atmo,Bfield,Input,Inf_Stokes,Inf_Nodes,Sol, &
                 LM_Stru,imask,saving)

      ! Memory in local structure
      MRAMc = MRAMc - 1d-6*sizeof(LM_Stru)

      return

      end subroutine Fit

!#####################################################################
!#####################################################################
!#####################################################################

      end module inversion_mod
