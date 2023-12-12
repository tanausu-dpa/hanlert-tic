      !> Initialization of the inversion
      module initinv_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/23/2023
!  Last version:
!     11/24/2023 V3.1.13
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/24/2023:   V3.1.13 - Added verbosity of mask file (TdPA)
!
!     10/04/2023:   V3.1.12 - Improved the verbosity of the PSF
!                             option (TdPA)
!                           - Updated verbosity of the type of
!                             inversion (TdPA)
!
!     10/03/2023:   V3.1.11 - Only master verbose the initialization
!                             of the ad-hoc asymmetry (HL)
!
!     09/28/2023:   V3.1.10 - Added re_set_nodes routine (TdPA)
!
!     08/30/2023:    V3.1.9 - Diffuse light factor node value not
!                             correctly initialized when initializing
!                             the nodes (TdPA)
!                           - Verbosity update (TdPA)
!
!     08/29/2023:    V3.1.8 - Diffuse light factor not correctly
!                             initialized when generating the
!                             stratification (TdPA)
!
!     08/24/2023:    V3.1.7 - Wrong formatting in new FWHM option
!                             verbosity (TdPA)
!
!     08/24/2023:    V3.1.6 - Updated the verbosity of the inversion
!                             set-up for the new possibilities of
!                             FWHM inputs (TdPA)
!
!     08/17/2023:    V3.1.5 - Initialize axially symmetric velocity
!                             pointers when creating the model
!                             atmosphere stratification (TdPA)
!
!     08/11/2023:    V3.1.4 - Inversion weights are no longer managed
!                             here. Changed the verbosity of that
!                             variable (TdPA)
!
!     07/31/2023:    V3.1.3 - Change the verbosity level in the
!                             inversion (HL)
!
!     07/06/2023:    V3.1.2 - Atmo_Stratify no longer processes the
!                             stratification inputs, that has been
!                             moved to rAtmo_mod in order to only do
!                             it once (TdPA)
!
!     07/05/2023:    V3.1.1 - Bugfix: Duplicated message when Pg_bound
!                             is not specified in hydrostatic
!                             equilibrium (TdPA)
!
!     07/03/2023:    V3.1.0 - The initial model atmosphere is always
!                             initialized when calling the routines
!                             in this module. Therefore, the routines
!                             Init_Guess_Thermal, Init_Guess_Mag,
!                             Restore_Thermal_Nodes, and
!                             Restore_Mag_Nodes have been
!                             removed (TdPA)
!                           - Added set_up_inversion and set_up_limits
!                             routines to initialize the inversion
!                             parameters. The code in them used to be
!                             in the TIC routine (TdPA)
!                           - Removed initial checks in Init_Nodes
!                             because they are elsewhere (TdPA)
!                           - Renamed Set_Nodes to Locate_Nodes (TdPA)
!                           - Moved Atmo_Stratify here, which also
!                             takes care of interpolating the input
!                             model into the new grid (TdPA)
!                           - Added set_nodes to do get the node
!                             values from the model atmosphere (TdPA)
!                           - Added Initialize_Nodes to manage the
!                             initialization of the node values (TdPA)
!
!     06/12/2023:    V3.0.3 - Rename the variable Inf_File (HL)
!
!     11/04/2023:    V3.0.2 - Bugfix: Vmi is a pointer. *1d-6/c will
!                             change the intialization for the next
!                             pixel. Now VmiC and VmiP are already
!                             in the code units (HL)
!
!     03/15/2023:    V3.0.1 - Commented hard-coded model variables not
!                             needed (TdPA)
!                           - Hard-coded model variables are no longer
!                             parameters, but pointers for additional
!                             flexibility (TdPA)
!                           - The initial model can be now specified
!                             via the input (TdPA)
!                           - Do not leave the initial guess routines
!                             just because there are no nodes, as we
!                             initialize the model there (TdPA)
!                           - Initializing atmosphere and magnetic
!                             field models in Init_Guess_* (TdPA)
!                           - Removed height branch in the node
!                             initialization (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/23/2023:    V0.0.0 - Started from 05/12/2020
!                             TIC@init_mod.f90 revision (TdPA)
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
!    set_up_inversion:
!      Copy from Input to the relevant inversion structures, set-up
!    and communicate the inverison parameters
!
!    set_up_limits:
!      Determine the limits of T, v, and B, as well as is the model
!    is static/dynamic or un/magnetized. A possible outcome is that
!    this needs to be checked pixelwise
!
!    Init_Guess_Thermal:
!      Guess the initial node values for thermodynamic variables
!
!    Init_Guess_Mag:
!      Guess the initial node values for magnetic variables
!
!    Init_Nodes:
!      Initialize the nodes
!
!    Locate_Nodes:
!      Equally spaces nodes positions
!
!    Atmo_Stratify:
!      Define the optical depth stratification from the input
!    parameters and interpolate the initial model into it
!
!    set_nodes:
!      Determine node values from the model atmosphere for a given
!    variable
!
!    re_set_nodes:
!      Changes the initial value of nodes in the transversal B or
!    v variables when initializing from other inversion and they
!    are too small
!
!    Initialize_Nodes:
!      Manage the initialization of the node values from the initial
!    model atmosphere
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use bounds_mod
      use commons_mod
      use inter_mod
      use model_mod
      use parameters_mod, only: c, RAD , TINYSP , PI , TINYB
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up the inversion structures\n
      !!       Input(Input_class): Structure with settings data\n
      !!   Inf_Nodes(Nodes_class): Structure with nodes data\n
      !! Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                           data\n
      !!      Sol(Solution_class): Class with the data of the RT
      !!                           solution
      subroutine set_up_inversion(Input,Inf_Nodes,Inf_Stokes,Sol)

      ! I/O
      type(Input_class):: Input
      type(Nodes_class):: Inf_Nodes
      type(Stokes_class):: Inf_Stokes
      type(Solution_class):: Sol

      ! Local
      character(len=10):: caux

      integer:: ii,jj,kk,ll,ia


      !
      ! Copy variables into node
      !
      Inf_Nodes%Node_type = Input%Node_type
      Inf_Nodes%Num_nodes = Input%Num_nodes
      do ii=1,Input%nvar
        if (allocated(Input%Node(ii)%H)) then
          allocate(Inf_Nodes%Node(ii)%H( &
                   Inf_Nodes%Num_nodes(ii)))
          Inf_Nodes%Node(ii)%H = Input%Node(ii)%H
          deallocate(Input%Node(ii)%H)
        end if
        Inf_Nodes%Node(ii)%Bounds = Input%Node(ii)%Bounds
        Inf_Nodes%Node(ii)%nebound = Input%Node(ii)%nebound
        if (allocated(Input%Node(ii)%ebound)) then
          allocate(Inf_Nodes%Node(ii)%ebound(4, &
                       Input%Node(ii)%nebound))
          Inf_Nodes%Node(ii)%ebound = Input%Node(ii)%ebound
          deallocate(Input%Node(ii)%ebound)
        end if
      end do
      Inf_Nodes%Nodes_flags = Input%Nodes_flags
      Inf_Nodes%Nodes_Regul = Input%Nodes_Regul
      Inf_Nodes%Indx_regul = Input%Indx_regul
      Inf_Nodes%Regul_weight = Input%Regul_weight
      Inf_Nodes%Scal = Input%Scal
      Inf_Nodes%Perturb = Input%Perturb
      Inf_Nodes%min_rel_Pert = Input%min_rel_Pert
      deallocate(Input%Node_type)
      deallocate(Input%Node)
      deallocate(Input%Num_nodes)
      deallocate(Input%Nodes_flags)
      deallocate(Input%Nodes_Regul)
      deallocate(Input%Indx_regul)
      deallocate(Input%Regul_weight)
      deallocate(Input%Scal)
      deallocate(Input%Perturb)
      deallocate(Input%min_rel_Pert)
      Inf_Nodes%Interpolation = Input%Interpolation
      Inf_Nodes%Btype = Input%Btype
      Inf_Nodes%vtype = Input%vtype
      Inf_Nodes%Pos_Correction = Input%Pos_Correction
      Inf_Nodes%Threshold_svd = Input%Threshold_svd
      Inf_Nodes%hydroeq = Input%hydroeq
      Inf_Nodes%Pg_bound = Input%Pg_bound
      Inf_Nodes%Max_step = Input%Max_step


      !
      ! Hard-coded type for pressure (if hydrostatic eq) and diffuse
      ! light
      !
      if (Inf_Nodes%hydroeq) &
        Inf_Nodes%Node_Type(Inf_Nodes%index_Pg) = 0
      Inf_Nodes%Node_Type(Inf_Nodes%index_f) = 0


      !
      ! Check if need to repeat hydrostatic equilibrium
      !
      Inf_Nodes%hydros = (Inf_Nodes% &
                               Nodes_flags(Inf_Nodes%index_T).or. &
                          Inf_Nodes% &
                               Nodes_flags(Inf_Nodes%index_Pg)).and. &
                         Inf_Nodes%hydroeq


      !
      ! Copy variable into Stokes
      !
      Inf_Stokes%auto_weight = Input%auto_weight

      !
      ! Copy variable into solution
      !
      Sol%fractional = Input%fractional
      Sol%Projection = Input%Projection


      !
      ! Sanity checks
      !

      ! There is scattering polarization
      if (Krad.gt.0.or.Kcut.gt.0) then

        ! If inverting the magnetic field vector or if it is
        ! in the LOS
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Bt).or. &
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Bp).or. &
            (Input%Btype.eq.1.and. &
             Inf_Nodes%Nodes_Flags(Inf_Nodes%index_B))) then

          ! Cannot be axial
          if (Input%nPh.le.0) then

            ! Abort
            umsg = 'If scattering polarization is included and '// &
                   'you are inverting the magnetic field '// &
                   'vector or its longitudinal component, then '// &
                   'the polarization cannot use an axial '// &
                   'quadrature'
            urou = 'set_up_inv'
            call gabortedv

          end if ! Axial polarization
        end if ! Inverting B vector or along LOS
      end if ! Scattering polarization included

      ! If inverting the velocity vector or if it is in the LOS
      if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vx).or. &
          Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vy).or. &
          (Input%vtype.eq.1.and. &
           Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vz))) then

        ! Cannot be axial
        if (Input%nPhI.le.0.or. &
            (Input%force.ne.'I'.and.Input%nPh.le.0)) then

          ! Abort
          umsg = 'You are inverting the horizontal velocity '// &
                 'vector or its longitudinal component, then '// &
                 'the you cannot use an axial quadrature'
          urou = 'set_up_inv'
          call gabortedv

        end if ! Axial polarization
      end if ! Inverting B vector or along LOS


      !
      ! Manage regularization inputs
      !

      ! If regularization limit larger than zero
      if (Input%Regul_Limit.ge.0d0) then

        ! For each variable
        do ii=1,Input%nvar

          ! If inverting the variable, with regularization and
          ! some weight
          Inf_Nodes%Nodes_Regul(ii) = Inf_Nodes%Nodes_Flags(ii).and. &
                                  Inf_Nodes%Indx_regul(ii).gt.0.and. &
                                  Inf_Nodes%Regul_Weight(ii).gt.0d0
        end do ! Variables

        ! Check if any true
        Inf_Nodes%Regul_flag = any(Inf_Nodes%Nodes_Regul)

      ! The limit is negative
      else

        ! No regularization
        Inf_Nodes%Nodes_Regul = .False.
        Inf_Nodes%Regul_flag = .False.

      end if ! Limit in regularization


      !
      ! Prepare for hydrostatic equilibrium
      !

      ! Initialize constant value for everyone
      Inf_Nodes%Const = 0d0

      ! If hydrostatic equilibrium
      if (Inf_Nodes%hydroeq) then

        ! If inverting Pgas
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg)) then

          ! Hard-code number of nodes, scale, and
          ! perturbation
          Inf_Nodes%Num_Nodes(Inf_Nodes%index_Pg) = 1
          Inf_Nodes%Scal(Inf_Nodes%index_Pg) = 2d0
          Inf_Nodes%Perturb(Inf_Nodes%index_Pg) = 0.05d0
          Inf_Nodes%Node_Type(Inf_Nodes%index_Pg) = 0

          ! Allocate the single node and set to 0
          if (allocated(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H)) &
            deallocate(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H)
          allocate(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H(1))
          Inf_Nodes%Node(Inf_Nodes%index_Pg)%H(1) = 0d0

          ! Automatic Pg_bound
          Inf_Nodes%Pg_Auto = Inf_Nodes%Pg_Bound.gt.0d0

          ! If pressure has a regularization weight
          if (Inf_Nodes%Regul_weight(Inf_Nodes%index_Pg).gt.0d0) then

            ! Set regularization to yes
            Inf_Nodes%Nodes_Regul(Inf_Nodes%index_Pg) = .True.

            ! If the regularization is other than constant,
            ! there is aproblem
            if (Inf_Nodes%Indx_regul(Inf_Nodes%index_Pg).ne.0.and. &
                Inf_Nodes%Indx_regul(Inf_Nodes%index_Pg).ne.2.and. &
                Inf_Nodes%Indx_regul(Inf_Nodes%index_Pg).ne.5) then
              umsg = 'Wrong regularization for pressure in '// &
                     'hydrostatic equilibrium'
              urou = 'set_up_inv'
              call gabortedv
            end if

          ! No regularization weight
          else
            ! Flag to no regulatization
            Inf_Nodes%Nodes_Regul(Inf_Nodes%index_Pg) = .False.
          end if

        ! Not inverting Pgas
        else

            ! Flag everything to no
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg) = .False.
            Inf_Nodes%Pg_Auto = .False.
            Inf_Nodes%Num_Nodes(Inf_Nodes%index_Pg) = 0
            Inf_Nodes%Nodes_Regul(Inf_Nodes%index_Pg) = .False.

        end if ! Inverting Pgas
      end if ! Hydrostatic equilibrium


      !
      ! Prepare for diffuse light
      !

      ! If inverting diffuse light
      if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_f)) then

        ! Hard-code number of nodes, scale, and
        ! perturbation
        Inf_Nodes%Num_Nodes(Inf_Nodes%index_f) = 1
        Inf_Nodes%Node_Type(Inf_Nodes%index_f) = 0

        ! Allocate the single node and set to 0
        if (allocated(Inf_Nodes%Node(Inf_Nodes%index_f)%H)) &
          deallocate(Inf_Nodes%Node(Inf_Nodes%index_f)%H)
        allocate(Inf_Nodes%Node(Inf_Nodes%index_f)%H(1))
        Inf_Nodes%Node(Inf_Nodes%index_f)%H(1) = 0d0

        ! If it has a regularization weight
        if (Inf_Nodes%Regul_weight(Inf_Nodes%index_f).gt.0d0) then

          ! Set regularization to yes
          Inf_Nodes%Nodes_Regul(Inf_Nodes%index_f) = .True.

          ! If the regularization is other than constant,
          ! there is aproblem
          if (Inf_Nodes%Indx_regul(Inf_Nodes%index_f).ne.0.and. &
              Inf_Nodes%Indx_regul(Inf_Nodes%index_f).ne.2.and. &
              Inf_Nodes%Indx_regul(Inf_Nodes%index_f).ne.5) then
            umsg = 'Wrong regularization for diffuse light'
            urou = 'set_up_inv'
            call gabortedv
          end if

        ! No regularization weight
        else
          ! Flag to no regulatization
          Inf_Nodes%Nodes_Regul(Inf_Nodes%index_f) = .False.
        end if

      ! Not inverting Pgas
      else

          ! Flag everything to no
          Inf_Nodes%Num_Nodes(Inf_Nodes%index_f) = 0
          Inf_Nodes%Nodes_Regul(Inf_Nodes%index_f) = .False.

      end if ! Inverting Pgas


      !
      ! Number of changing nodes
      !

      ! For each variable
      do ii=1,Input%nvar

        ! If inverting
        if (Inf_Nodes%Nodes_Flags(ii)) then

          ! Check type of node
          select case(Inf_Nodes%Node_Type(ii))

            ! Value all nodes
            case(0)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)
              Inf_Nodes%Node_Vary(1,ii) = 1

            ! Value all nodes except last
            case(1)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)-1
              Inf_Nodes%Node_Vary(1,ii) = 1

            ! Value all nodes except first
            case(2)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)-1
              Inf_Nodes%Node_Vary(1,ii) = 2

            ! Value all nodes except extremes
            case(3)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)-2
              Inf_Nodes%Node_Vary(1,ii) = 2

            ! Correction
            case(4)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)
              Inf_Nodes%Node_Vary(1,ii) = 1

            ! Correction all nodes execept last
            case(5)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)-1
              Inf_Nodes%Node_Vary(1,ii) = 1

            ! Correction all nodes except first
            case(6)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)-1
              Inf_Nodes%Node_Vary(1,ii) = 2

            ! Correction all nodes except extremes
            case(7)

              ! Get number of changing nodes and first to change
              Inf_Nodes%Num_Vary(ii) = Inf_Nodes%Num_Nodes(ii)-2
              Inf_Nodes%Node_Vary(1,ii) = 2

          end select ! Type of node

          ! Compute last node that changes
          Inf_Nodes%Node_Vary(2,ii) = Inf_Nodes%Node_Vary(1,ii) + &
                                      Inf_Nodes%Num_Vary(ii) - 1

        ! Not inverting
        else

            ! No nodes and no change
            Inf_Nodes%Num_Nodes(ii) = 0
            Inf_Nodes%Num_Vary(ii) = 0
            Inf_Nodes%Node_Vary(:,ii) = 0

        end if ! Inverting or not

      end do ! Variables


      !
      ! Count nodes
      !

      ! Initialize
      Inf_Nodes%Num_Fit = 0
      Inf_Nodes%Num_Mag = 0
      Inf_Nodes%Num_Thermal = 0
      Inf_Nodes%Num_Asymmetry = 0
      Inf_Nodes%Num_glob = 0

      ! Magnetic variables
      ii = Inf_Nodes%index_B
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Mag = Inf_Nodes%Num_Mag + &
                            Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_Bt
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Mag = Inf_Nodes%Num_Mag + &
                            Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_Bp
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Mag = Inf_Nodes%Num_Mag + &
                            Inf_Nodes%Num_Vary(ii)
      end if

      ! Thermal variables
      ii = Inf_Nodes%index_T
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_vx
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_vy
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_vz
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_vm
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_Pg
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      ! Global variables
      ii = Inf_Nodes%index_f
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_glob = Inf_Nodes%Num_glob + &
                             Inf_Nodes%Num_Vary(ii)
      end if

      ! Asymmetry variables
      ii = Inf_Nodes%index_J21R
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_J21I
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_J22R
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if
      ii = Inf_Nodes%index_J22I
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if

      ! Total nodes
      Inf_Nodes%Num_Fit = Inf_Nodes%Num_Mag + &
                          Inf_Nodes%Num_Thermal + &
                          Inf_Nodes%Num_glob + &
                          Inf_Nodes%Num_Asymmetry

      ! No Nodes?
      if (Inf_Nodes%Num_Fit.le.0) then
        umsg = 'There are no nodes to fit'
        urou = 'set_up_inv'
        call gabortedv
      end if



      !
      ! Show the input to the user
      !
      if (gpid.eq.0) then

        ! Information message
        umsg = ' - Called the inversion module with the '// &
               'following parameters:'
        call verboseI(1)


        write(umsg,'(A,i1)') '   o Verbosity level:',vlevel
        call verboseI(1)


        write(umsg,'(A,A)') '   o Input data file: ', &
                             trim(Input%Filename_ob)
        call verboseI(1)

        ! If mask
        if (trim(Input%Inv_mask).ne.'NONE'.and. &
            trim(Input%Inv_init).ne.'INIT') then
          write(umsg,'(A,A)') '   o Mask file: ', &
                              trim(Input%Inv_mask)
          call verboseI(1)
        end if

        if (Input%Type_Inversion.eq.0) then
          write(umsg,'(A)') '   o Thermodynamic inversion'
        else if (Input%Type_Inversion.eq.1) then
          write(umsg,'(A)') '   o Magnetic inversion'
        else if (Input%Type_Inversion.eq.2) then
          write(umsg,'(A)') '   o Full inversion together'
        else if (Input%Type_Inversion.eq.3) then
          write(umsg,'(A)') '   o Full inversion sequential'
        else
          write(umsg,'(A)') '   o Full inversion sequential, '// &
                            'first non-magnetic, then only '// &
                            'magnetic'
        end if
        call verboseI(1)


        if (Input%auto_weight) then
          write(umsg,'(A)') '   o Automatic Stokes weights'
          call verboseI(1)
        else
          ! File
          if (Input%linv_weight) then
            write(umsg,'(A,A)') &
              '   o Stokes weights from file ',trim(Input%inv_weight)
            call verboseI(1)
          ! Numbers
          else
            do ia=1,Input%Num_weight
              write(umsg,'(A,2(1x,f9.3),A,4(1x,es10.3))') &
                '   o Stokes weights between range', &
                Input%weight(4,ia),Input%weight(5,ia), &
                ' nm :',Input%weight(0:3,ia)
              call verboseI(1)
            end do
          end if
        end if


        if (trim(Input%Inv_init).ne.'NONE') then
          write(umsg,'(A,A)') '   o Restore from file ', &
                               Input%Inv_init
          call verboseI(1)
        end if


        if (Input%centered) then
          write(umsg,'(A)') '   o Centered derivative for RF'
        else
          write(umsg,'(A)') '   o Non-centered derivative for RF'
        end if
        call verboseI(1)


        write(umsg,'(A,i5)') '   o Maximum iterations ', &
                             Input%Num_Iter
        call verboseI(1)


        if (Input%Interpolation.eq.0) then
          write(umsg,'(A)') '   o Linear interpolation'
        else if (Input%Interpolation.eq.1) then
          write(umsg,'(A)') '   o Quadratic bezier interpolation'
        else if (Input%Interpolation.eq.2) then
          write(umsg,'(A)') '   o Cubic bezier interpolation'
        end if
        call verboseI(1)


        if (Input%btype.eq.0) then
          write(umsg,'(A)') '   o Magnetic field in vertical ref.'
        else if (Input%btype.eq.1) then
          write(umsg,'(A)') '   o Magnetic field in LOS ref.'
        end if
        call verboseI(1)


        if (Input%vtype.eq.0) then
          write(umsg,'(A)') '   o Cartesian velocity'
        else if (Input%vtype.eq.1) then
          write(umsg,'(A)') '   o Velocity in LOS ref.'
        end if
        call verboseI(1)


        if (Input%Pos_Correction) then
          write(umsg,'(A)') '   o Correct node positions'
        else
          write(umsg,'(A)') '   o Do not correct node positions'
        end if
        call verboseI(1)


        ! Variables
        do ii=1,Input%nvar

          if (ii.eq.Inf_Nodes%index_B) then

            if (Input%btype.eq.0) then
              write(umsg,'(A,i2,A)') '   o |B| (',ii,'):'
            else if (Input%btype.eq.1) then
              write(umsg,'(A,i2,A)') '   o Bpar (',ii,'): '
            end if

          else if (ii.eq.Inf_Nodes%index_Bt) then

            if (Input%btype.eq.0) then
              write(umsg,'(A,i2,A)') '   o B_incl (',ii,'):'
            else if (Input%btype.eq.1) then
              write(umsg,'(A,i2,A)') '   o Btrv (',ii,'): '
            end if

          else if (ii.eq.Inf_Nodes%index_Bp) then

            if (Input%btype.eq.0) then
              write(umsg,'(A,i2,A)') '   o B_azim (',ii,'):'
            else if (Input%btype.eq.1) then
              write(umsg,'(A,i2,A)') '   o B_pos-azim (',ii,'):'
            end if

          else if (ii.eq.Inf_Nodes%index_T) then

            write(umsg,'(A,i2,A)') '   o Temperature: (',ii,')'

          else if (ii.eq.Inf_Nodes%index_vx) then

            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle
            if (Input%vtype.eq.0) then
              write(umsg,'(A,i2,A)') '   o X axis velocity (',ii,'):'
            else if (Input%vtype.eq.1) then
              write(umsg,'(A,i2,A)') '   o POS velocity (',ii,'):'
            end if

          else if (ii.eq.Inf_Nodes%index_vy) then

            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle
            if (Input%vtype.eq.0) then
              write(umsg,'(A,i2,A)') '   o Y axis velocity (',ii,'):'
            else if (Input%vtype.eq.1) then
             write(umsg,'(A,i2,A)') '   o Velocity POS azim (',ii,'):'
            end if

          else if (ii.eq.Inf_Nodes%index_vz) then

            if (Input%vtype.eq.0) then
             write(umsg,'(A,i2,A)') '   o Vertical velocity (',ii,'):'
            else if (Input%vtype.eq.1) then
              write(umsg,'(A,i2,A)') '   o LOS velocity (',ii,'):'
            end if

          else if (ii.eq.Inf_Nodes%index_vm) then

            write(umsg,'(A,i2,A)') '   o Turbulent velocity (',ii,'):'

          else if (ii.eq.Inf_Nodes%index_Pg) then

            write(umsg,'(A,i2,A)') '   o Gas pressure (',ii,'):'

          else if (ii.eq.Inf_Nodes%index_f) then

            write(umsg,'(A,i2,A)') '   o Diffuse light (',ii,'):'

          else if (ii.eq.Inf_Nodes%index_J21R) then

            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle
            write(umsg,'(A,i2,A)') '   o J21 real (',ii,'):'

          else if (ii.eq.Inf_Nodes%index_J21I) then

            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle
            write(umsg,'(A,i2,A)') '   o J21 imag (',ii,'):'

          else if (ii.eq.Inf_Nodes%index_J22R) then

            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle
            write(umsg,'(A,i2,A)') '   o J22 real (',ii,'):'

          else if (ii.eq.Inf_Nodes%index_J22I) then

            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle
            write(umsg,'(A,i2,A)') '   o J22 imag (',ii,'):'

          else

            write(umsg,'(A)') '   o ??'
            cycle

          end if
          call verboseI(1)

          if (Inf_Nodes%Nodes_Flags(ii)) then
            write(umsg,'(A)') '     + Inverting'
          else
            write(umsg,'(A)') '     + Fixed'
          end if
          call verboseI(1)

          ! Pressure
          if (ii.eq.Inf_Nodes%index_Pg.and.Input%hydroeq) then

            write(umsg,'(A)') '     + Hydrostatic equilibrium'
            call verboseI(1)
            if (Input%Pg_bound.gt.0d0) then
              write(umsg,'(A,es10.3)') '     + Boundary value ', &
                                       Input%Pg_bound
              call verboseI(1)
            end if

          ! Non-pressure
          else

            if (.not.Inf_Nodes%Nodes_Flags(ii)) cycle

            if (Inf_Nodes%Node_type(ii).eq.0) then

              write(umsg,'(A)') '     + Nodes have values'

            else if (Inf_Nodes%Node_type(ii).eq.1) then

              write(umsg,'(A)') '     + Nodes have values '// &
                                'with last node fixed'

            else if (Inf_Nodes%Node_type(ii).eq.2) then

              write(umsg,'(A)') '     + Nodes have values '// &
                                'with first node fixed'

            else if (Inf_Nodes%Node_type(ii).eq.3) then

              write(umsg,'(A)') '     + Nodes have values '// &
                                'with extreme nodes fixed'

            else if (Inf_Nodes%Node_type(ii).eq.4) then

              write(umsg,'(A)') '     + Nodes have corrections'

            else if (Inf_Nodes%Node_type(ii).eq.5) then

              write(umsg,'(A)') '     + Nodes have corrections '// &
                                'with last node fixed'

            else if (Inf_Nodes%Node_type(ii).eq.6) then

              write(umsg,'(A)') '     + Nodes have corrections '// &
                                'with first node fixed'

            else if (Inf_Nodes%Node_type(ii).eq.7) then

              write(umsg,'(A)') '     + Nodes have corrections '// &
                                'with extreme nodes fixed'

            end if
            call verboseI(1)

            write(umsg,'(A,i3)') '     + Number of nodes ', &
                                 Inf_Nodes%Num_nodes(ii)
            call verboseI(1)

            if (allocated(Inf_Nodes%Node(ii)%H)) then
              write(caux,'(i10)') Inf_Nodes%Num_nodes(ii)
              kk = 0
              do jj=1,10
                if (caux(jj:jj).ne.' ') then
                  kk = jj
                  exit
                end if
              end do
              ll = 10 - kk + 1
              do jj=1,ll
                caux(jj:jj) = caux(kk+jj-1:kk+jj-1)
              end do
              do jj=ll+1,10
                caux(jj:jj) = ' '
              end do
              write(umsg,'(A,'//trim(caux)//'(1x,es10.3))') &
                 '     + Fixed positions ',Inf_Nodes%Node(ii)%H
              call verboseI(1)
            end if

            write(umsg,'(A,2(1x,es10.3))') &
              '     + Scale ',Inf_Nodes%Scal(ii)
            call verboseI(1)

            write(umsg,'(A,2(1x,es10.3))') &
              '     + Perturbation ',Inf_Nodes%Perturb(ii)
            call verboseI(1)

            write(umsg,'(A,2(1x,es10.3))') &
              '     + Bound values ',Inf_Nodes%Node(ii)%Bounds
            call verboseI(1)

            ! Special limits
            do jj=1,Inf_Nodes%Node(ii)%nebound
              write(umsg,'(A,2(1x,es10.3),A,2(1x,es10.3))') &
                '     + Bound values between', &
                Inf_Nodes%Node(ii)%ebound(1:2,jj),': ', &
                Inf_Nodes%Node(ii)%ebound(3:4,jj)
              call verboseI(1)
            end do ! special limits

          end if ! Pressure


          if (Inf_Nodes%Nodes_regul(ii)) then

            if (Inf_Nodes%Indx_regul(ii).eq.1) then

              write(umsg,'(A,es10.3)') &
                '     + Mean regularization with weight ', &
                Inf_Nodes%Regul_weight(ii)

            else if (Inf_Nodes%Indx_regul(ii).eq.2) then

              write(umsg,'(A,es10.3)') &
                '     + Constant regularization with weight ', &
                Inf_Nodes%Regul_weight(ii)

            else if (Inf_Nodes%Indx_regul(ii).eq.3) then

              write(umsg,'(A,es10.3)') &
                '     + First derivative regularization with '// &
                'weight ',Inf_Nodes%Regul_weight(ii)

            else if (Inf_Nodes%Indx_regul(ii).eq.4) then

              write(umsg,'(A,es10.3)') &
                '     + Second derivative regularization with '// &
                'weight ',Inf_Nodes%Regul_weight(ii)

            end if
            call verboseI(1)

          end if ! Regularization

        end do ! Variables


        write(umsg,'(A,i4)') '   o Total number of nodes ', &
                             Inf_Nodes%Num_Fit
        call verboseI(1)


        if (Inf_Nodes%Num_Mag.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of magnetic '// &
                               'nodes ',Inf_Nodes%Num_Mag
          call verboseI(1)
        end if


        if (Inf_Nodes%Num_Thermal.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of thermal '// &
                               'nodes ',Inf_Nodes%Num_Thermal
          call verboseI(1)
        end if


        if (Inf_Nodes%Num_glob.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of global '// &
                               'nodes ',Inf_Nodes%Num_glob
          call verboseI(1)
        end if


        if (Inf_Nodes%Num_Asymmetry.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of symmetry '// &
                               'nodes ',Inf_Nodes%Num_Asymmetry
          call verboseI(1)
        end if


        if (maxval(Inf_Nodes%Indx_regul).gt.0) then
          write(umsg,'(A,es10.3)') &
            '   o Regularization limit',Input%Regul_Limit
          call verboseI(1)
        end if

        write(umsg,'(A,es10.3)') '   o Chi2 threshold ', &
                                 Input%Threshold_chisq
        call verboseI(1)

        write(umsg,'(A,es10.3)') '   o Chi2 MRC ', &
                                 Input%Chisq_fraction
        call verboseI(1)

        if (Input%SVD_type.eq.0) then
          write(umsg,'(A,es10.3)') '   o SVD method: tradiational'
        else if (Input%SVD_type.eq.2) then
          write(umsg,'(A,es10.3)') '   o SVD method: SIR-like'
        end if
        call verboseI(1)


        write(umsg,'(A,es10.3)') '   o SVD threshold ', &
                                 Input%Threshold_svd
        call verboseI(1)


        write(umsg,'(A,A)') '   o Initial atmospheric model file: ', &
                             trim(Input%atmo)
        call verboseI(1)


        write(umsg,'(A,i4)') '   o Atmospheric nodes in '// &
                             'synthesis: ',Input%Atmo_Input
        call verboseI(1)


        write(umsg,'(A,es10.3)') '   o SVD maximum step ', &
                                 Input%Max_Step
        call verboseI(1)


        if (Input%Err_type.eq.0) then

          write(umsg,'(A)') '   o Error from Hessian'
          call verboseI(1)

        else if (Input%Err_type.eq.1) then

          write(umsg,'(A)') '   o Error not from Hessian'
          call verboseI(1)

        else if (Input%Err_type.eq.2) then

          write(umsg,'(A)') '   o Error from RF'
          call verboseI(1)

        else if (Input%Err_type.eq.3) then

          write(umsg,'(A)') '   o Error worst from Hessian or RF'
          call verboseI(1)

        end if


        if (allocated(Input%lim_fwhm)) then
          if (Input%lim_fwhm(1)%nn.gt.1) then
            write(umsg,'(A,i2)') '   o PSF number of ranges ', &
                                 Input%lim_fwhm(1)%nn
          end if
          do ii=1,Input%lim_fwhm(1)%nn
            ! Gaussian
            if (Input%lim_fwhm(ii)%gaussian) then
              if (Input%lim_fwhm(ii)%doub(2).lt.1d99) then
                write(umsg,'(A,es10.3," (",es10.3,",",es10.3,")")') &
                       '   o PSF FWHM ',Input%lim_fwhm(ii)%doub(3), &
                        Input%lim_fwhm(ii)%doub(1:2)
              else
                write(umsg,'(A,es10.3)') &
                       '   o PSF FWHM ',Input%lim_fwhm(ii)%doub(3)
              end if
              call verboseI(1)
            ! No gaussian
            else
              if (Input%lim_fwhm(ii)%doub(2).lt.1d99) then
                write(umsg,'(A,A," (",es10.3,",",es10.3,")")') &
                    '   o PSF from file',trim(Input%fwhm_fil(ii)%str), &
                     Input%lim_fwhm(ii)%doub(1:2)
              else
                write(umsg,'(A,A)') &
                    '   o PSF from file',trim(Input%fwhm_fil(ii)%str)
              end if
              call verboseI(1)
            end if
          end do
        end if


        if (Input%PopuInit) then
          write(umsg,'(A)') '   o Initialize RF with solution'
        else
          write(umsg,'(A)') '   o Do not initialize RF with solution'
        end if
        call verboseI(1)


        if (Input%Fractional) then
          write(umsg,'(A)') '   o Use fractional Stokes parameters'
        else
          write(umsg,'(A)') '   o Use absolute Stokes parameters'
        end if
        call verboseI(1)


        write(umsg,'(A,2(1x,es10.3))') &
          '   o Range of taus ',Input%Tau_Range
        call verboseI(1)


        if (Input%Broyden) then

          write(umsg,'(A)') '   o LM Broyden'
          call verboseI(1)

        end if


        if (Input%LM_Method.eq.0) then

          write(umsg,'(A,es10.3)') '   o Traditional LM'

        else if (Input%LM_Method.eq.1) then

          write(umsg,'(A,es10.3)') '   o LM with backtracking'

        end if
        call verboseI(1)

        if (Input%Sigma_neglect) then

          write(umsg,'(A)') '   o Neglect sigma'
          call verboseI(1)

        end if
      end if ! Spit out the inputs


      !
      ! Transform units in boundaries, scale, and perturbation
      ! for velocities
      !

      ! For vx
      ii = Inf_Nodes%index_vx

      ! Transform units
      Inf_Nodes%Perturb(ii) = Inf_Nodes%Perturb(ii)*1d-6/c
      Inf_Nodes%Scal(ii) = Inf_Nodes%Scal(ii)*1d-6/c
      Inf_Nodes%Node(ii)%bounds = Inf_Nodes%Node(ii)%bounds*1d-6/c

      ! For each special limit
      do jj=1,Inf_Nodes%Node(ii)%nebound

        ! Transform units
        Inf_Nodes%Node(ii)%ebound(3:4,jj) = &
                           Inf_Nodes%Node(ii)%ebound(3:4,jj)*1d-6/c

      end do ! special limit

      ! For vy if cartesian
      if (Input%vtype.eq.0) then

        ii = Inf_Nodes%index_vy

        ! Transform units
        Inf_Nodes%Perturb(ii) = Inf_Nodes%Perturb(ii)*1d-6/c
        Inf_Nodes%Scal(ii) = Inf_Nodes%Scal(ii)*1d-6/c
        Inf_Nodes%Node(ii)%bounds = Inf_Nodes%Node(ii)%bounds*1d-6/c

        ! For each special limit
        do jj=1,Inf_Nodes%Node(ii)%nebound

          ! Transform units
          Inf_Nodes%Node(ii)%ebound(3:4,jj) = &
                             Inf_Nodes%Node(ii)%ebound(3:4,jj)*1d-6/c

        end do ! special limit

      end if

      ! For vz
      ii = Inf_Nodes%index_vz

      ! Transform units
      Inf_Nodes%Perturb(ii) = Inf_Nodes%Perturb(ii)*1d-6/c
      Inf_Nodes%Scal(ii) = Inf_Nodes%Scal(ii)*1d-6/c
      Inf_Nodes%Node(ii)%bounds = Inf_Nodes%Node(ii)%bounds*1d-6/c

      ! For each special limit
      do jj=1,Inf_Nodes%Node(ii)%nebound

        ! Transform units
        Inf_Nodes%Node(ii)%ebound(3:4,jj) = &
                           Inf_Nodes%Node(ii)%ebound(3:4,jj)*1d-6/c

      end do ! special limit

      ! For micro
      ii = Inf_Nodes%index_vm

      ! Transform units
      Inf_Nodes%Perturb(ii) = Inf_Nodes%Perturb(ii)*1d-6/c
      Inf_Nodes%Scal(ii) = Inf_Nodes%Scal(ii)*1d-6/c
      Inf_Nodes%Node(ii)%bounds = Inf_Nodes%Node(ii)%bounds*1d-6/c

      ! For each special limit
      do jj=1,Inf_Nodes%Node(ii)%nebound

        ! Transform units
        Inf_Nodes%Node(ii)%ebound(3:4,jj) = &
                           Inf_Nodes%Node(ii)%ebound(3:4,jj)*1d-6/c

      end do ! special limit


      !
      ! Trick in azimuth if limits are complete
      !

      ! If inverting magnetic field azimuth
      ii = Inf_Nodes%index_Bp
      if (Inf_Nodes%Nodes_flags(ii)) then

        ! If full range
        if (abs(Inf_Nodes%Node(ii)%bounds(2) - &
                Inf_Nodes%Node(ii)%bounds(1)).ge.(2d0*PI-TINYA)) then

          ! Move limits by pi
          Inf_Nodes%Node(ii)%bounds(2) = &
                                     Inf_Nodes%Node(ii)%bounds(2) + PI
          Inf_Nodes%Node(ii)%bounds(1) = &
                                     Inf_Nodes%Node(ii)%bounds(1) - PI

        end if ! Full azimuth range

        ! Special ranges
        do jj=1,Inf_Nodes%Node(ii)%nebound

          ! If full range
          if (abs(Inf_Nodes%Node(ii)%ebound(4,jj) - &
                  Inf_Nodes%Node(ii)%ebound(3,jj)).ge. &
              (2d0*PI-TINYA)) then

            ! Move limits by pi
            Inf_Nodes%Node(ii)%ebound(4,jj) = &
                                  Inf_Nodes%Node(ii)%ebound(4,jj) + PI
            Inf_Nodes%Node(ii)%ebound(3,jj) = &
                                  Inf_Nodes%Node(ii)%ebound(3,jj) - PI

          end if ! Full azimuth range

        end do ! special ranges

      end if ! Inverting azimuth

      ! If inverting velocity azimuth
      ii = Inf_Nodes%index_vy
      if (Inf_Nodes%Nodes_flags(ii).and.Input%vtype.eq.1) then

        ! If full range
        if (abs(Inf_Nodes%Node(ii)%bounds(2) - &
                Inf_Nodes%Node(ii)%bounds(1)).ge.(2d0*PI-TINYA)) then

          ! Move limits by pi
          Inf_Nodes%Node(ii)%bounds(2) = &
                                     Inf_Nodes%Node(ii)%bounds(2) + PI
          Inf_Nodes%Node(ii)%bounds(1) = &
                                     Inf_Nodes%Node(ii)%bounds(1) - PI

        end if ! Full azimuth range

        ! Special ranges
        do jj=1,Inf_Nodes%Node(ii)%nebound

          ! If full range
          if (abs(Inf_Nodes%Node(ii)%ebound(4,jj) - &
                  Inf_Nodes%Node(ii)%ebound(3,jj)).ge. &
              (2d0*PI-TINYA)) then

            ! Move limits by pi
            Inf_Nodes%Node(ii)%ebound(4,jj) = &
                                  Inf_Nodes%Node(ii)%ebound(4,jj) + PI
            Inf_Nodes%Node(ii)%ebound(3,jj) = &
                                  Inf_Nodes%Node(ii)%ebound(3,jj) - PI

          end if ! Full azimuth range

        end do ! special ranges

      end if ! Inverting azimuth

      end subroutine set_up_inversion

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up the limits for some variables\n
      !!       Input(Input_class): Structure with settings data\n
      !!   Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!      Atmo_in(Atmo_class): Structure with atmospheric data
      !!                           read from model atmosphere\n
      !!  Bfield_in(Bfield_class): Structure with the magnetic field
      !!                           data read from the input\n
      !!             maxB(double): Maximum magnetic field\n
      !!     update_tlim(logical): If T limit has to be updated
      !!                           at every pixel\n
      !!     update_vlim(logical): If v limits have to be updated
      !!                           at every pixel\n
      !!     update_blim(logical): If B limit has to be updated
      !!                           at every pixel
      subroutine set_up_limits(Input,Inf_Nodes,Atmo_in, &
                               Bfield_in,maxB,update_tlim, &
                               update_vlim,update_blim)

      ! I/O
      type(Input_class):: Input
      type(Nodes_class):: Inf_Nodes
      type(Atmo_class), intent(in):: Atmo_in
      type(Bfield_class), intent(inout):: Bfield_in
      logical, intent(out):: update_tlim,update_vlim,update_blim
      double precision, intent(out):: maxB

      ! Local
      integer:: ii,jj

      double precision:: maxvx,maxvy,maxvz

      !
      ! Get minimum and maximum temperatures, and maximum velocity
      !

      ! If inverting temperature
      ii = Inf_Nodes%index_T
      if (Inf_Nodes%Nodes_Flags(ii)) then

        ! Temperature limits
        Input%minT = Inf_Nodes%Node(ii)%Bounds(1)
        Input%maxT = Inf_Nodes%Node(ii)%Bounds(2)
        ! Temperature special limits
        do jj=1,Inf_Nodes%Node(ii)%nebound
          if (Inf_Nodes%Node(ii)%ebound(1,jj).lt.Input%minT) &
            Input%minT = Inf_Nodes%Node(ii)%ebound(1,jj)
          if (Inf_Nodes%Node(ii)%ebound(2,jj).gt.Input%maxT) &
            Input%maxT = Inf_Nodes%Node(ii)%ebound(2,jj)
        end do

        ! No need to update
        update_tlim = .False.

      ! Not inverting temperature
      else

        ! If atmosphere already set
        if (Input%atmoin_type.eq.0) then

          ! Temperature limits
          Input%minT = minval(Atmo_in%T)
          Input%maxT = maxval(Atmo_in%T)

          ! No need to update
          update_tlim = .False.

        else

          ! Temperature limits (guess)
          Input%minT = 3d3
          Input%maxT = 5d5

          ! Need to update
          update_tlim = .True.

        end if ! Type of model input
      end if ! Inverting T or not

      ! Initialize
      maxvx = 0d0
      maxvy = 0d0
      maxvz = 0d0

      ! Inverting velocities
      if ((Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vx).or. &
           Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vy).or. &
           Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vz)).and. &
           (Input%vtype.eq.0.or. &
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vx).or. &
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vz))) then

        ! vx
        ii=Inf_Nodes%index_vx
        if (Inf_Nodes%Nodes_Flags(ii)) then

          ! Maximum from boundaries
          maxvx = maxval(abs(Inf_Nodes%Node(ii)%Bounds))

          ! Special limits
          do jj=1,Inf_Nodes%Node(ii)%nebound

            ! New maximum
            maxvx = max(maxvx, &
                       maxval(abs(Inf_Nodes%Node(ii)%ebound(3:4,jj))))

          end do ! special limits

        end if ! Inverting

        ! vy
        ii=Inf_Nodes%index_vy
        if (Inf_Nodes%Nodes_Flags(ii).and.Input%vtype.eq.0) then

          ! Maximum from boundaries
          maxvy = maxval(abs(Inf_Nodes%Node(ii)%Bounds))

          ! Special limits
          do jj=1,Inf_Nodes%Node(ii)%nebound

            ! New maximum
            maxvy = max(maxvy, &
                       maxval(abs(Inf_Nodes%Node(ii)%ebound(3:4,jj))))

          end do ! special limits

        end if ! Inverting

        ! vz
        ii=Inf_Nodes%index_vz
        if (Inf_Nodes%Nodes_Flags(ii)) then

          ! Maximum from boundaries
          maxvz = maxval(abs(Inf_Nodes%Node(ii)%Bounds))

          ! Special limits
          do jj=1,Inf_Nodes%Node(ii)%nebound

            ! New maximum
            maxvz = max(maxvz, &
                       maxval(abs(Inf_Nodes%Node(ii)%ebound(3:4,jj))))

          end do ! special limits

        end if ! Inverting

        ! No need to update
        update_vlim = .False.

      ! No inverting velocities
      else

        ! If atmosphere already set
        if (Input%atmoin_type.eq.0) then

          ! Velocity limit
          maxvz = maxval(sqrt(Atmo_in%vx*Atmo_in%vx + &
                              Atmo_in%vy*Atmo_in%vy + &
                              Atmo_in%vz*Atmo_in%vz))

          ! No need to update
          update_vlim = .False.

        else

          ! Need to update
          update_vlim = .True.

        end if ! Type of model input
      end if ! Inverting velocity or not

      ! If to update later
      if (update_vlim) then

        ! Asumme dynamic 5 km/s to start with
        dyn = .True.
        Input%static = .False.
        Input%maxV = 5d-6/c

      ! No need to update
      else

        ! Combine
        maxvz = sqrt(maxvx*maxvx + maxvy*maxvy + maxvz*maxvz)

        ! If maximum velocity
        if (maxvz.gt.0d0) then
          dyn = .True.
          Input%maxV = maxvz
          Input%static = .False.
        else
          dyn = .False.
          Input%maxV = 0d0
          Input%static = .True.
        end if

      end if ! To update later


      !
      ! Check if unmagnetized
      !

      ! Only thermal inversion
      if (Input%Type_inversion.eq.0) then

        ! No field
        maxB = 0d0
        update_blim = .False.

      ! Not only thermal
      else

        ! Inverting B
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_B).or. &
            (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Bt).and. &
             Input%Btype.eq.1)) then

          ! Bstrength or BLOS
          ii = Inf_Nodes%index_B
          maxB = maxval(abs(Inf_Nodes%Node(ii)%Bounds))

          ! Special limits
          do jj=1,Inf_Nodes%Node(ii)%nebound

            ! Maximum
            maxB = max(maxB, &
                       maxval(abs(Inf_Nodes%Node(ii)%ebound(3:4,jj))))

          end do

          ! If magnetic field in LOS frame
          if (Input%Btype.eq.1) then

            ! Bpos
            ii = Inf_Nodes%index_Bt
            maxB = max(maxB,maxval(abs(Inf_Nodes%Node(ii)%Bounds)))

            ! Special limits
            do jj=1,Inf_Nodes%Node(ii)%nebound

              ! Maximum
              maxB = max(maxB, &
                       maxval(abs(Inf_Nodes%Node(ii)%ebound(3:4,jj))))

            end do

          end if ! LOS reference frame

          ! No need to update
          update_blim = .False.

        ! Not inverting B
        else

          ! If atmosphere already set
          if (Input%atmoin_type.eq.0) then

            ! Bield limit
            maxB = maxval(Bfield_in%Bstrength)

            ! No need to update
            update_blim = .False.

          ! Atmo not set
          else

            ! Need to update
            update_blim = .True.

          end if ! Type of model input
        end if ! Inverting B
      end if ! Type of inversion


      ! If to update later
      if (update_blim) then

        ! Asumme magnetic to start with
        Input%unmagnetized = .False.
        maxB = 1d0

      ! No need to update
      else

        ! Magnetic field?
        if (maxB.le.TINYB) then
          Input%unmagnetized = .True.
        else
          Input%unmagnetized = .False.
        end if

      end if ! To update later

      end subroutine set_up_limits

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the nodes\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Input(Input_class): Structure with settings data\n
      !!   Inf_Nodes(Nodes_class): Structure with nodes data
      subroutine Init_Nodes(Atmo,Input,Inf_Nodes)


      ! IO
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input
      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local
      integer:: i, j, k

      double precision:: Delta_Tau, Mini, H_min, H_max
      double precision, dimension(:), allocatable:: TAU_LG

      !
      ! Get height limits
      H_min = log10(Atmo%z(nZ))
      H_max = log10(Atmo%z(1))


      ! For each variable
      do i=1,Input%nvar

        ! If there are nodes and inverting
        if (Inf_Nodes%Num_Nodes(i).gt.0.and. &
            Inf_Nodes%Nodes_Flags(i)) then

          ! Reset values and errors
          if (allocated(Inf_Nodes%Node(i)%Var)) &
            deallocate(Inf_Nodes%Node(i)%Var)
          if (allocated(Inf_Nodes%Node(i)%Errors)) &
            deallocate(Inf_Nodes%Node(i)%Errors)
          allocate(Inf_Nodes%Node(i)%Var(Inf_Nodes%Num_Nodes(i)))
          allocate(Inf_Nodes%Node(i)%Errors(Inf_Nodes%Num_Nodes(i)))

          ! If pressure and hydrostatic
          if ((i.eq.Inf_Nodes%index_Pg).and.Inf_Nodes%hydroeq) then

            ! Forget about the input nodes
            if (allocated(Inf_Nodes%Node(i)%H)) &
              deallocate(Inf_Nodes%Node(i)%H)
            allocate(Inf_Nodes%Node(i)%H(Inf_Nodes%Num_Nodes(i)))

            ! Set in upper boundary
            Inf_Nodes%Node(i)%H(1) = Input%Tau_range(1)

          ! If diffuse light
          else if (i.eq.Inf_Nodes%index_f) then

            ! Forget about the input nodes
            if (allocated(Inf_Nodes%Node(i)%H)) &
              deallocate(Inf_Nodes%Node(i)%H)
            allocate(Inf_Nodes%Node(i)%H(Inf_Nodes%Num_Nodes(i)))

            ! Set in upper boundary
            Inf_Nodes%Node(i)%H(1) = Input%Tau_range(1)

          ! No stratification specified for nodes
          else if (.not.allocated(Inf_Nodes%Node(i)%H)) then

            ! Set equally spaced nodes
            call Locate_Nodes(H_Min,H_Max, &
                              Inf_Nodes%Num_Nodes(i), &
                              Inf_Nodes%Node(i)%H)

          end if ! Input stratification

        ! No consistent for invertion
        else

          ! Flag and set to 0
          Inf_Nodes%Nodes_flags(i) = .False.
          Inf_Nodes%Num_Nodes(i) = 0

        end if ! Inverting variable

      end do ! Variables

      ! If correct position
      if (Inf_Nodes%Pos_Correction) then

        ! Allocate auxiliar
        allocate(TAU_LG(Atmo%nZ))
        TAU_LG = log10(Atmo%z)

        ! Variables
        do i=1,Input%nvar

          ! Skip if not inverting
          if (Inf_Nodes%Num_Nodes(i).le.0.or. &
              .not.Inf_Nodes%Nodes_Flags(i)) cycle

          ! Skip diffuse light or hydrostatic Pg
          if (Inf_Nodes%index_f.eq.i.or. &
              (Inf_Nodes%index_Pg.eq.i.and.Inf_Nodes%hydroeq)) cycle

          ! Reset tau index
          if (allocated(Inf_Nodes%Node(i)%Tau_Indx)) &
            deallocate(Inf_Nodes%Node(i)%Tau_Indx)
          allocate(Inf_Nodes%Node(i)%Tau_Indx( &
                                            Inf_Nodes%Num_Nodes(i)))

          ! For each node
          do j=1,Inf_Nodes%Num_Nodes(i)

            ! Initialize minimum and index
            Mini = abs(Inf_Nodes%Node(i)%H(j)-TAU_LG(1))
            Inf_Nodes%Node(i)%Tau_Indx(j) = 1

            ! For each height in the atmosphere
            do k=2,Atmo%nZ

              ! Get difference
              Delta_Tau = abs(Inf_Nodes%Node(i)%H(j) - TAU_LG(k))

              ! If new minimum
              if (Delta_Tau.lt.Mini) then

                ! Store minimum and index
                Mini = Delta_Tau
                Inf_Nodes%Node(i)%Tau_Indx(j) = k

              end if ! New minimum

            end do ! Atmosphere nodes

            ! Set node to closest atmosphere node
            Inf_Nodes%Node(i)%H(j) = &
                               TAU_LG(Inf_Nodes%Node(i)%Tau_Indx(j))

          end do ! For each node
        end do ! Variables

        ! Deallocate auxiliar
        deallocate(TAU_LG)

      end if ! Correct position

      ! Reset indexing for nodes
      if (allocated(Inf_Nodes%Inf_Inv)) deallocate(Inf_Nodes%Inf_Inv)
      allocate(Inf_Nodes%Inf_Inv(2,Inf_Nodes%Num_Fit))

      ! Initialize counter
      j = 1

      ! For each variable
      do k=1,Input%nvar

        ! If inverting
        if (Inf_Nodes%Nodes_Flags(k)) then

          ! Count nodes with variation
          do i=Inf_Nodes%Node_Vary(1,k),Inf_Nodes%Node_Vary(2,k)

            ! Index nodes linearly
            Inf_Nodes%Inf_Inv(1,j) = k
            Inf_Nodes%Inf_Inv(2,j) = i

            ! Advance counter
            j = j+1

          end do ! Nodes with variation

        end if ! Inverting

      end do ! Variables

      return

      end subroutine Init_Nodes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Equally spaces nodes positions\n
      !!      H_min(double): Minimum position of the nodes\n
      !!      H_max(double): Maximum position of the nodes\n
      !!       Num(integer): Number of nodes\n
      !!  Height(double(:)): Node positions
      subroutine Locate_Nodes(H_min,H_max,Num,Height)

      ! IO
      integer, intent(in):: Num
      double precision, intent(in):: H_min, H_max
      double precision, dimension(:), allocatable, intent(out):: &
                                                                Height
      ! Local
      integer:: i

      double precision:: dx


      ! Allocate positions
      allocate(Height(Num))

      ! If more than two positions
      if (Num.ge.2) then

        ! Top boundary
        Height(1) = H_max

        ! Get equal step
        dx = (H_max-H_min)/(Num-1)

        ! Setup rest of nodes except the last
        do i=2,Num-1
          Height(i) = H_max-(i-1)*dx
        end do

      end if ! More than two positions

      ! The last is the other boundary
      Height(Num) = H_min

      return

      end subroutine Locate_Nodes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Create an optical depth stratification for the model
      !! atmosphere in the inversion\n
      !!     Atmo_in(Atmo_class): Structure with input atmospheric
      !!                          data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !! Bfield_in(Bfield_class): Structure with the magnetic field
      !!                          data read from the input\n
      !!    Bfield(Bfield_class): Structure with the magnetic field
      !!                          data\n
      !!            z(double(:)): Stratification
      subroutine Atmo_Stratify(Atmo_in,Atmo,Bfield_in,Bfield,z)

      ! IO
      type(Atmo_class), intent(in):: Atmo_in
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(in):: Bfield_in
      type(Bfield_class), intent(inout):: Bfield
      double precision, dimension(:), intent(in):: z

      ! Local
      integer:: i

      double precision, dimension(Atmo_in%nz):: LTAUI


      ! Get log optical depth of input
      LTAUI = log10(Atmo_in%z)

      ! Set optical depth scale, 500 nm ref, solar gravity and
      ! electron density type
      Atmo%scal = 'T'
      Atmo%typo = Atmo_in%typo
      Atmo%tfreq = Atmo_in%tfreq
      Atmo%logg = Atmo_in%logg

      ! Set global variables
      ztau = .True.

      ! Allocate variables
      allocate(Atmo%z(nZ),Atmo%T(nZ),Atmo%ne(nZ),Atmo%Pg(nZ))
      allocate(Atmo%vmi(nZ),Atmo%vx(nZ),Atmo%vy(nZ),Atmo%vz(nZ))
      allocate(Atmo%nh(nZ,6),Atmo%nht(nZ),Atmo%nha(nZ),Atmo%nhm(nZ))
      allocate(Atmo%nhe(1,1),Atmo%zeros(nZ))
      allocate(Bfield%Bstrength(nZ))
      allocate(Bfield%Btheta(nZ))
      allocate(Bfield%Bphi(nZ))
      nullify(Atmo%Bx,Atmo%By,Atmo%Bz)
      nullify(Atmo%vxa,Atmo%vya,Atmo%vza)

      ! Allocs
      Atmo%alloc_a = .True.
      Atmo%alloc_b = .False.

      ! Copy rest of variables
      Atmo%NT = Atmo_in%NT
      Atmo%nele = Atmo_in%nele
      allocate(Atmo%ele(Atmo%nele),Atmo%pT(Atmo%NT))
      allocate(Atmo%abund(size(Atmo_in%abund)))
      Atmo%ele = Atmo_in%ele
      Atmo%abund = Atmo_in%abund
      Atmo%pT = Atmo_in%pT

      ! Initialize everything
      Atmo%nhe(1,1) = -1
      Atmo%zeros = 0d0
      Atmo%T = 0d0
      Atmo%ne = 0d0
      Atmo%vmi = 0d0
      Atmo%vx = 0d0
      Atmo%vy = 0d0
      Atmo%vz = 0d0
      Atmo%nh = 0d0
      Atmo%nht = 0d0
      Atmo%nha = 0d0
      Atmo%nhm = 0d0
      Atmo%Pg = 0d0
      Atmo%f_diff = Atmo_in%f_diff

      ! Copy stratification
      Atmo%z = z


      !
      ! Interpolate variables
      !

      ! Temperature
      call Intpol(LTAUI, Atmo_in%T, Atmo_in%nZ, &
                  Atmo%z, Atmo%T, Atmo%nZ, 2, 1)

      ! Pgas
      call Intpol(LTAUI, Atmo_in%Pg, Atmo_in%nZ, &
                  Atmo%z, Atmo%Pg, Atmo%nZ, 2, 1)

      ! Micro
      call Intpol(LTAUI, Atmo_in%vmi, Atmo_in%nZ, &
                  Atmo%z, Atmo%vmi, Atmo%nZ, 2, 1)

      ! vx
      call Intpol(LTAUI, Atmo_in%vx, Atmo_in%nZ, &
                  Atmo%z, Atmo%vx, Atmo%nZ, 2, 1)

      ! vy
      call Intpol(LTAUI, Atmo_in%vy, Atmo_in%nZ, &
                  Atmo%z, Atmo%vy, Atmo%nZ, 2, 1)

      ! vz
      call Intpol(LTAUI, Atmo_in%vz, Atmo_in%nZ, &
                  Atmo%z, Atmo%vz, Atmo%nZ, 2, 1)

      ! If magnetic field
      if (maxval(Bfield_in%Bstrength).gt.TINYB) then

        ! B
        call Intpol(LTAUI, Bfield_in%Bstrength, Atmo_in%nZ, &
                    Atmo%z, Bfield%Bstrength, Atmo%nZ, 2, 1)
        ! Btheta
        call Intpol(LTAUI, Bfield_in%Btheta, Atmo_in%nZ, &
                    Atmo%z, Bfield%Btheta, Atmo%nZ, 2, 1)
        ! Bphi
        call Intpol(LTAUI, Bfield_in%Bphi, Atmo_in%nZ, &
                    Atmo%z, Bfield%Bphi, Atmo%nZ, 2, 1)

        !
        ! Sanity
        do i=1,Atmo%nz

          ! No field
          if (Bfield%Bstrength(i).le.TINYB) then

            ! Just zero
            Bfield%Bstrength(i) = 0d0
            Bfield%Btheta(i) = 0d0
            Bfield%Bphi(i) = 0d0

          ! Check angles
          else

            if (Bfield%Btheta(i).lt.0d0) Bfield%Btheta(i) = 0d0
            if (Bfield%Btheta(i).gt.PI) Bfield%Btheta(i) = PI
            if (Bfield%Bphi(i).lt.0d0) Bfield%Bphi(i) = &
                                                    Bfield%Bphi(i) + &
                                                    2d0*PI
            if (Bfield%Bphi(i).gt.2d0*PI) Bfield%Bphi(i) = &
                                                    Bfield%Bphi(i) - &
                                                    2d0*PI

          end if ! Bfield value

        end do


      ! No field
      else

        ! Just zero
        Bfield%Bstrength = 0d0
        Bfield%Btheta = 0d0
        Bfield%Bphi = 0d0

      end if

      ! If JKQ in input
      if (allocated(Atmo_in%JKQin)) then

        ! J21R
        if (maxval(abs(Atmo_in%JKQin(4*Atmo_in%nz+1: &
                                     5*Atmo_in%nz))).gt.0d0) &
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(4*Atmo_in%nz+1:5*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(4*Atmo%nz+1:5*Atmo%nz), &
                      Atmo%nZ, 2, 1)
        ! J21I
        if (maxval(abs(Atmo_in%JKQin(5*Atmo_in%nz+1: &
                                     6*Atmo_in%nz))).gt.0d0) &
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(5*Atmo_in%nz+1:6*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(5*Atmo%nz+1:6*Atmo%nz), &
                      Atmo%nZ, 2, 1)
        ! J22R
        if (maxval(abs(Atmo_in%JKQin(6*Atmo_in%nz+1: &
                                     7*Atmo_in%nz))).gt.0d0) &
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(6*Atmo_in%nz+1:7*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(6*Atmo%nz+1:7*Atmo%nz), &
                      Atmo%nZ, 2, 1)
        ! J22I
        if (maxval(abs(Atmo_in%JKQin(7*Atmo_in%nz+1: &
                                     8*Atmo_in%nz))).gt.0d0) &
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(7*Atmo_in%nz+1:8*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(7*Atmo%nz+1:8*Atmo%nz), &
                      Atmo%nZ, 2, 1)
      end if

      ! Make the scale linear
      Atmo%z = 10d0**Atmo%z

      return

      end subroutine Atmo_Stratify

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate incoming variable into nodes
      !!         tau(double(:)): Log10 optical depth scale\n
      !!         var(double(:)): Stratificaiton of any variable\n
      !!             nn(iteger): Dimension of tau and var\n
      !! Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!           indx(iteger): Index of current variable
      subroutine set_nodes(tau,var,nn,Inf_Nodes,indx)

      ! I/O
      type(Nodes_class), intent(inout):: Inf_Nodes
      integer, intent(in):: nn, indx
      double precision, dimension(:), intent(in):: tau,var

      ! Local
      character(3):: length
      character(30):: fmt

      integer::j


      ! Inverting variable
      if (Inf_Nodes%Nodes_Flags(indx)) then

        ! If by value
        if (Inf_Nodes%Node_Type(indx).le.3) then

          ! Interpolate into nodes
          call Intpol(tau, var, nn, &
                      Inf_Nodes%Node(indx)%H, &
                      Inf_Nodes%Node(indx)%Var, &
                      Inf_Nodes%Num_Nodes(indx), 2, 1)

          ! Check boundaries
          call CheckBounds(Inf_Nodes%Node(indx), &
                           Inf_Nodes%Num_Nodes(indx))

        ! If by correction
        else

          ! Initialize
          Inf_Nodes%Node(indx)%Var = 0d0

        end if ! Type of node

        ! Master
        if (pid.eq.0) then

          ! Verbose
          write(umsg, '(A,i2)') "   Parameter index = ",indx
          call verboseI(3)

          ! Get format
          write(length, "(i3)") Inf_Nodes%Num_Nodes(indx)
          fmt = '(A,'//trim(adjustl(length))//'es15.4)'
          fmt = trim(adjustl(fmt))

          ! Verbose positions
          write(umsg, FMT=fmt) "   Position: ", &
            (Inf_Nodes%Node(indx)%H(j), j=1, &
             Inf_Nodes%Num_Nodes(indx))
          call verboseI(3)

          !
          ! Verbose values
          !
          ! If velocity
          if (indx.eq.Inf_Nodes%index_vz.or. &
              indx.eq.Inf_Nodes%index_vx.or. &
              (indx.eq.Inf_Nodes%index_vy.and. &
               Inf_Nodes%vtype.eq.0).or. &
              indx.eq.Inf_Nodes%index_vm) then

            write(umsg, FMT=fmt) "     Values: ", &
             (Inf_Nodes%Node(indx)%Var(j)*1d6*c, j=1, &
              Inf_Nodes%Num_Nodes(indx))
            call verboseI(3)

          ! No velocity
          else

            write(umsg, FMT=fmt) "     Values: ", &
             (Inf_Nodes%Node(indx)%Var(j), j=1, &
              Inf_Nodes%Num_Nodes(indx))
            call verboseI(3)

          end if ! Velocity
        end if ! Master
      end if ! Inverting variable

      end subroutine set_nodes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Changed the node values into a given value\n
      !! Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!           indx(iteger): Index of current variable\n
      !!            val(double): Value to set
      subroutine re_set_nodes(Inf_Nodes,indx,val)

      ! I/O
      type(Nodes_class), intent(inout):: Inf_Nodes
      integer, intent(in):: indx
      double precision, intent(in):: val

      ! Local
      character(3):: length
      character(30):: fmt

      integer::j

      ! Not inverting variable, skip
      if (.not.Inf_Nodes%Nodes_Flags(indx)) return

      ! If not by value, skip
      if (Inf_Nodes%Node_Type(indx).gt.3) return

      ! Set value
      Inf_Nodes%Node(indx)%var = val

      ! Master
      if (pid.eq.0) then

        ! Verbose
        write(umsg, '(A,i2)') "   Parameter index = ",indx
        call verboseI(3)

        ! Get format
        write(length, "(i3)") Inf_Nodes%Num_Nodes(indx)
        fmt = '(A,'//trim(adjustl(length))//'es15.4)'
        fmt = trim(adjustl(fmt))

        ! Verbose positions
        write(umsg, FMT=fmt) "   Position: ", &
          (Inf_Nodes%Node(indx)%H(j), j=1, &
           Inf_Nodes%Num_Nodes(indx))
        call verboseI(3)

        !
        ! Verbose values
        !
        ! If velocity
        if (indx.eq.Inf_Nodes%index_vz.or. &
            indx.eq.Inf_Nodes%index_vx.or. &
            (indx.eq.Inf_Nodes%index_vy.and. &
             Inf_Nodes%vtype.eq.0).or. &
            indx.eq.Inf_Nodes%index_vm) then

          write(umsg, FMT=fmt) "     Values: ", &
           (Inf_Nodes%Node(indx)%Var(j)*1d6*c, j=1, &
            Inf_Nodes%Num_Nodes(indx))
          call verboseI(3)

        ! No velocity
        else

          write(umsg, FMT=fmt) "     Values: ", &
           (Inf_Nodes%Node(indx)%Var(j), j=1, &
            Inf_Nodes%Num_Nodes(indx))
          call verboseI(3)

        end if ! Velocity
      end if ! Master

      end subroutine re_set_nodes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up the node values from the initial model atmosphere\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Bfield(Bfield_class): Structure with the magnetic field\n
      !! Inf_Nodes(Nodes_class): Structure with nodes data
      subroutine Initialize_Nodes(Atmo,Bfield,Inf_Nodes)

      ! I/O
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local

      integer:: ii

      double precision, dimension(Atmo%nz):: TAU_LG


      ! Get log tau
      TAU_LG = log10(Atmo%z)


      !
      ! Magnetic field
      !
      if (Inf_Nodes%Num_Mag.gt.0) then

        ! Master
        if (pid.eq.0) then

          ! Verbose
          write(umsg, '(A)') &
            ' - Magnetic parameters initialized: '
          call verboseI(3)
          write(umsg, '(A,i1)') &
            '   Btype is ', Inf_Nodes%Btype
          call verboseI(3)

        end if

        ! Type of magnetic field
        select case (Inf_Nodes%Btype)

          ! In vertical
          case(0)

            ! Bstrength
            ii = Inf_Nodes%index_B

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Bstrength, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Bstrength
            ii = Inf_Nodes%index_Bt

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Btheta, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Bstrength
            ii = Inf_Nodes%index_Bp

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Bphi, &
                           Atmo%nz,Inf_Nodes,ii)

          ! In LOS
          case(1)

            ! Bstrength
            ii = Inf_Nodes%index_B

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Blos, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Bstrength
            ii = Inf_Nodes%index_Bt

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Bpos, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Bstrength
            ii = Inf_Nodes%index_Bp

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Azimuth, &
                           Atmo%nz,Inf_Nodes,ii)
        end select

      end if ! Magnetic nodes


      !
      ! Thermal variables
      !
      if (Inf_Nodes%Num_Thermal.gt.0) then

        ! Verbose
        if (pid.eq.0) then
          write(umsg, '(A)') &
            ' - Thermodynamic parameters initialized: '
          call verboseI(3)
          write(umsg, '(A,i1)') &
            '   vtype is ', Inf_Nodes%vtype
          call verboseI(3)
        end if

        ! Temperature
        ii = Inf_Nodes%index_T

        ! Get nodes
        call set_nodes(TAU_LG,Atmo%T,Atmo%nz, &
                       Inf_Nodes,ii)

        ! Type of velocity
        select case (Inf_Nodes%vtype)

          ! In vertical
          case(0)

            ! X velocity
            ii = Inf_Nodes%index_vx

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vx,Atmo%nz, &
                           Inf_Nodes,ii)

            ! Y velocity
            ii = Inf_Nodes%index_vy

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vy,Atmo%nz, &
                           Inf_Nodes,ii)

            ! Vertical velocity
            ii = Inf_Nodes%index_vz

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vz,Atmo%nz, &
                           Inf_Nodes,ii)

          ! In LOS
          case(1)

            ! POS velocity
            ii = Inf_Nodes%index_vx

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vpos,Atmo%nz, &
                           Inf_Nodes,ii)

            ! Azimuth velocity
            ii = Inf_Nodes%index_vy

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vphi,Atmo%nz, &
                           Inf_Nodes,ii)

            ! LOS velocity
            ii = Inf_Nodes%index_vz

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vlos,Atmo%nz, &
                           Inf_Nodes,ii)

        end select

        ! Micro. velocity
        ii = Inf_Nodes%index_vm

        ! Get nodes
        call set_nodes(TAU_LG,Atmo%vmi,Atmo%nz, &
                       Inf_Nodes,ii)

        ! Gas pressure
        ii = Inf_Nodes%index_Pg

        ! Hydrostatic
        if (Inf_Nodes%hydroeq.and.Inf_Nodes%Nodes_Flags(ii)) then

          ! If value
          if (Inf_Nodes%Node_Type(ii).le.0) then

            ! Get boundary
            Inf_Nodes%Node(ii)%Var(1) = Inf_Nodes%Pg_Bound

          ! If correction
          else

            ! Set to zero
            Inf_Nodes%Node(ii)%Var(1) = 0d0

          end if ! Type of node

          ! Master
          if (pid.eq.0) then

            ! Verbose
            write(umsg, '(A,i2)') "   Parameter index = ",ii
            call verboseI(3)

            ! Verbose values
            write(umsg, '(A,es15.4)') "     Values: ", &
                                 Inf_Nodes%Node(ii)%Var(1)
            call verboseI(3)

          end if

        ! No hydrostatic
        else

          ! Get nodes
          call set_nodes(TAU_LG,Atmo%Pg,Atmo%nz, &
                         Inf_Nodes,ii)

        end if ! Hydrostatic eq.

      end if ! Thermal


      !
      ! Global
      !
      if (Inf_Nodes%Num_glob.gt.0) then

        ! Diffuse light
        ii = Inf_Nodes%index_f

        ! Get nodes
        if (Inf_Nodes%Nodes_flags(ii).and.pid.eq.0) then

          ! Set from atmosphere
          Inf_Nodes%Node(ii)%var(1) = Atmo%f_diff

          ! Verbose
          write(umsg, '(A)') &
            ' - Global parameters initialized: '
          call verboseI(1)

          ! Verbose
          write(umsg, '(A,i2)') "   Parameter index = ",ii
          call verboseI(1)

          ! Verbose values
          write(umsg, '(A,es15.4)') "     Values: ", &
                               Inf_Nodes%Node(ii)%Var(1)
          call verboseI(1)

        end if

      end if


      !
      ! Asymmetry variables
      !
      if (Inf_Nodes%Num_asymmetry.gt.0) then

        ! Verbose
        if (pid.eq.0) then
          write(umsg, '(A)') &
            ' - Asymmetry parameters initialized: '
          call verboseI(1)
        end if

        ! J21R
        ii = Inf_Nodes%index_J21R

        ! Get nodes
        call set_nodes(TAU_LG, &
                       Atmo%JKQin(4*Atmo%nz+1:5*Atmo%nz), &
                       Atmo%nz, Inf_Nodes,ii)

        ! J21I
        ii = Inf_Nodes%index_J21I

        ! Get nodes
        call set_nodes(TAU_LG, &
                       Atmo%JKQin(5*Atmo%nz+1:6*Atmo%nz), &
                       Atmo%nz, Inf_Nodes,ii)

        ! J22R
        ii = Inf_Nodes%index_J22R

        ! Get nodes
        call set_nodes(TAU_LG, &
                       Atmo%JKQin(6*Atmo%nz+1:7*Atmo%nz), &
                       Atmo%nz, Inf_Nodes,ii)

        ! J22I
        ii = Inf_Nodes%index_J22I

        ! Get nodes
        call set_nodes(TAU_LG, &
                       Atmo%JKQin(7*Atmo%nz+1:8*Atmo%nz), &
                       Atmo%nz, Inf_Nodes,ii)

      end if

      end subroutine Initialize_Nodes

!#####################################################################
!#####################################################################
!#####################################################################

      end module initinv_mod
