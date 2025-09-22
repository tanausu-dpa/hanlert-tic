      !> Initialization of the inversion
      module initinv_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC/NSSCC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     23/02/2023
!  Last version:
!     29/08/2025 V4.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     29/08/2025:    V4.0.5 - Made changes to accomodate the new
!                             input to determine the extrapolation
!                             mode (TdPA)
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
!    set_up_inversion
!      Set-up the inversion structures by copying the information
!    in the structure with the configuration data
!
!    set_up_limits
!      Determine the limits of T, v, and B, as well as if the model
!    is static/dynamic or un/magnetized. A possible outcome is that
!    this needs to be checked pixelwise
!
!    guess_polarity
!      Make a first guess of the magnetic field based on WFA
!
!    Init_Nodes
!      Initialize the location, values, and errors of the nodes
!
!    Locate_Nodes
!      Generate an array with equally spaced nodes within the given
!    limits
!
!    Atmo_Stratify
!      Create a model atmosphere by interpolating a given model into
!    a given optical depth stratification
!
!    set_nodes
!      Determine node values from the model atmosphere for a given
!    variable
!
!    re_set_nodes
!      Change the initial value of nodes in the transversal B or v
!    variables when initializing from other inversion and they are too
!    small
!
!    Initialize_Nodes
!      Initialize the node values from a given model atmosphere
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use bounds_mod
      use commons_mod
      use inter_mod
      use model_mod
      use parameters_mod, only: c, RAD , TINYSP , PI , TINYB , &
                                TINYVEL , wfac
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up the inversion structures by copying the information
      !! in the structure with the configuration data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data
      subroutine set_up_inversion(Input,Inf_Nodes,Inf_Stokes,Sol)

      ! I/O

      type(Input_class), intent(inout):: Input
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Solution_class), intent(inout):: Sol

      ! Local

      character(len=10):: caux

      logical:: check

      integer:: ii,jj,kk,ll,ia


      !
      ! Copy variables into node
      ! There are no changes in memory count because we are copying
      ! and removing the original
      !

      ! Type of node data and number of nodes
      Inf_Nodes%Node_type = Input%Node_type
      Inf_Nodes%Num_nodes = Input%Num_nodes

      ! For every variable
      do ii=1,Input%nvar

        ! If there is location data
        if (allocated(Input%Node(ii)%H)) then

          ! Copy to node structure and free input
          allocate(Inf_Nodes%Node(ii)%H( &
                   Inf_Nodes%Num_nodes(ii)))
          Inf_Nodes%Node(ii)%H = Input%Node(ii)%H
          deallocate(Input%Node(ii)%H)

        end if ! There is location data

        ! Copy bounds info
        Inf_Nodes%Node(ii)%Bounds = Input%Node(ii)%Bounds
        Inf_Nodes%Node(ii)%nebound = Input%Node(ii)%nebound

        ! If there is special bound info
        if (allocated(Input%Node(ii)%ebound)) then

          ! Copy to node structure and free input
          allocate(Inf_Nodes%Node(ii)%ebound(4, &
                       Input%Node(ii)%nebound))
          Inf_Nodes%Node(ii)%ebound = Input%Node(ii)%ebound
          deallocate(Input%Node(ii)%ebound)

        end if ! There is special bound info

      end do ! Variables

      ! Copy flags, regularization data, scaling, perturbation size,
      ! minimum relative perturbation, extrapolation
      Inf_Nodes%Nodes_flags = Input%Nodes_flags
      Inf_Nodes%Nodes_Regul = Input%Nodes_Regul
      Inf_Nodes%Indx_regul = Input%Indx_regul
      Inf_Nodes%Regul_weight = Input%Regul_weight
      Inf_Nodes%Scal = Input%Scal
      Inf_Nodes%Perturb = Input%Perturb
      Inf_Nodes%min_rel_Pert = Input%min_rel_Pert
      Inf_Nodes%extrapolation = Input%extrapolation

      ! Free the already copied information from the input
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
      deallocate(Input%extrapolation)

      ! Copy type of interpolation, type of magnetic field vector,
      ! type of velocity vector, if correcting node position, the
      ! SVD threshold, if doing hydrostatic equilibrium, the value
      ! of the gas pressure at the top boundary, and the maximum
      ! step in SVD
      Inf_Nodes%Interpolation = Input%Interpolation
      Inf_Nodes%Btype = Input%Btype
      Inf_Nodes%vtype = Input%vtype
      Inf_Nodes%Pos_Correction = Input%Pos_Correction
      Inf_Nodes%Threshold_svd = Input%Threshold_svd
      Inf_Nodes%hydroeq = Input%hydroeq
      Inf_Nodes%Pg_bound = Input%Pg_bound
      Inf_Nodes%Max_step = Input%Max_step


      ! The type of node is hard-coded type for pressure (if in
      ! hydrostatic eq) and for diffuse light
      if (Inf_Nodes%hydroeq) &
        Inf_Nodes%Node_Type(Inf_Nodes%index_Pg) = 0
      Inf_Nodes%Node_Type(Inf_Nodes%index_f) = 0


      ! Check if need to repeat hydrostatic equilibrium when
      ! updating the model, i.e., if temperature or gas pressure
      ! are free variables and the model is in hydrostatic
      ! equilibrium
      Inf_Nodes%hydros = (Inf_Nodes% &
                               Nodes_flags(Inf_Nodes%index_T).or. &
                          Inf_Nodes% &
                               Nodes_flags(Inf_Nodes%index_Pg)).and. &
                         Inf_Nodes%hydroeq


      ! Copy weights in Stokes structure
      Inf_Stokes%auto_weight = Input%auto_weight

      ! Copy if profiles are fractional and the kind of magnetic
      ! field projection in the solution structure
      Sol%fractional = Input%fractional
      Sol%Projection = Input%Projection


      !
      ! Sanity checks
      !

      ! There is scattering polarization
      if (Krad.gt.0.or.Kcut.gt.0) then

        ! If inverting the non-vertical components of the magnetic
        ! field vector or if inverting the magnetic field along the
        ! LOS
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Bt).or. &
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Bp).or. &
            (Input%Btype.eq.1.and. &
             Inf_Nodes%Nodes_Flags(Inf_Nodes%index_B))) then

          ! If the problem is axial
          if (Input%nPh.le.0) then

            ! Abort because it is not possible
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

      ! If inverting the horizontal components of the velocity vector
      ! or inverting in the LOS reference frame
      if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vx).or. &
          Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vy).or. &
          (Input%vtype.eq.1.and. &
           Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vz))) then

        ! If the problem is set to exial
        if (Input%nPhI.le.0.or. &
            (Input%force.ne.'I'.and.Input%nPh.le.0)) then

          ! Abort because it is not possible
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

      ! If regularization limit is larger than zero
      if (Input%Regul_Limit.ge.0d0) then

        ! For each variable
        do ii=1,Input%nvar

          ! Flag to regularize if inverting the variable, with
          ! regularization and it has some weight
          Inf_Nodes%Nodes_Regul(ii) = Inf_Nodes%Nodes_Flags(ii).and. &
                                  Inf_Nodes%Indx_regul(ii).gt.0.and. &
                                  Inf_Nodes%Regul_Weight(ii).gt.0d0
        end do ! Variables

        ! Check if there is any regulatization to do
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

          ! Hard-code number of nodes, scale, perturbation, and node
          ! type
          Inf_Nodes%Num_Nodes(Inf_Nodes%index_Pg) = 1
          Inf_Nodes%Scal(Inf_Nodes%index_Pg) = 2d0
          Inf_Nodes%Perturb(Inf_Nodes%index_Pg) = 0.05d0
          Inf_Nodes%Node_Type(Inf_Nodes%index_Pg) = 0

          ! If the node location is already set
          if (allocated(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H)) then

            ! Free location data
            MRAMc = MRAMc - &
                    1d-6*sizeof(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H)
            deallocate(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H)

          end if ! Already allocate

          ! Allocate single height
          allocate(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H(1))
          MRAMc = MRAMc + &
                  1d-6*sizeof(Inf_Nodes%Node(Inf_Nodes%index_Pg)%H)

          ! Set location at top
          Inf_Nodes%Node(Inf_Nodes%index_Pg)%H(1) = 0d0

          ! Flag if automatic Pg_bound
          Inf_Nodes%Pg_Auto = Inf_Nodes%Pg_Bound.gt.0d0

          ! If pressure has a regularization weight
          if (Inf_Nodes%Regul_weight(Inf_Nodes%index_Pg).gt.0d0) then

            ! Set regularization to yes
            Inf_Nodes%Nodes_Regul(Inf_Nodes%index_Pg) = .True.

            ! If the regularization is other than constant
            if (Inf_Nodes%Indx_regul(Inf_Nodes%index_Pg).ne.0.and. &
                Inf_Nodes%Indx_regul(Inf_Nodes%index_Pg).ne.2.and. &
                Inf_Nodes%Indx_regul(Inf_Nodes%index_Pg).ne.5) then

              ! Abort because it makes no sense
              umsg = 'Wrong regularization for pressure in '// &
                     'hydrostatic equilibrium'
              urou = 'set_up_inv'
              call gabortedv

            end if ! Regularization is anything but constant

          ! No regularization weight
          else

            ! Flag no regulatization
            Inf_Nodes%Nodes_Regul(Inf_Nodes%index_Pg) = .False.

          end if ! Check regulatization weight

        ! Not inverting Pgas
        else

            ! Flag everything to no and set no nodes
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

        ! Hard-code number of nodes and node type
        Inf_Nodes%Num_Nodes(Inf_Nodes%index_f) = 1
        Inf_Nodes%Node_Type(Inf_Nodes%index_f) = 0

        ! If the node location is already set
        if (allocated(Inf_Nodes%Node(Inf_Nodes%index_f)%H)) then

          ! Free location data
          MRAMc = MRAMc - &
                  1d-6*sizeof(Inf_Nodes%Node(Inf_Nodes%index_f)%H)
          deallocate(Inf_Nodes%Node(Inf_Nodes%index_f)%H)

        end if ! Already allocate

        ! Allocate single height
        allocate(Inf_Nodes%Node(Inf_Nodes%index_f)%H(1))
        MRAMc = MRAMc + &
                1d-6*sizeof(Inf_Nodes%Node(Inf_Nodes%index_f)%H)

        ! Set location at top
        Inf_Nodes%Node(Inf_Nodes%index_f)%H(1) = 0d0

        ! If it has a regularization weight
        if (Inf_Nodes%Regul_weight(Inf_Nodes%index_f).gt.0d0) then

          ! Set regularization to yes
          Inf_Nodes%Nodes_Regul(Inf_Nodes%index_f) = .True.

          ! If the regularization is other than constant
          if (Inf_Nodes%Indx_regul(Inf_Nodes%index_f).ne.0.and. &
              Inf_Nodes%Indx_regul(Inf_Nodes%index_f).ne.2.and. &
              Inf_Nodes%Indx_regul(Inf_Nodes%index_f).ne.5) then

            ! Abort because it makes no sense
            umsg = 'Wrong regularization for diffuse light'
            urou = 'set_up_inv'
            call gabortedv

          end if ! Regularization is anything but constant

        ! No regularization weight
        else

          ! Flag to regulatization
          Inf_Nodes%Nodes_Regul(Inf_Nodes%index_f) = .False.

        end if ! Check regulatization weight

      ! Not inverting diffuse light
      else

          ! Flag no regulatization and no nodes
          Inf_Nodes%Num_Nodes(Inf_Nodes%index_f) = 0
          Inf_Nodes%Nodes_Regul(Inf_Nodes%index_f) = .False.

      end if ! Inverting diffuse light


      !
      ! Count number of changing nodes per variable
      !

      ! For each variable
      do ii=1,Input%nvar

        ! If inverting the variable
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

      ! Initialize count
      Inf_Nodes%Num_Fit = 0
      Inf_Nodes%Num_Mag = 0
      Inf_Nodes%Num_Thermal = 0
      Inf_Nodes%Num_Asymmetry = 0
      Inf_Nodes%Num_glob = 0

      !
      ! Magnetic variables

      ! B strength or Blos
      ii = Inf_Nodes%index_B
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Mag = Inf_Nodes%Num_Mag + &
                            Inf_Nodes%Num_Vary(ii)
      end if

      ! B theta or Btrans
      ii = Inf_Nodes%index_Bt
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Mag = Inf_Nodes%Num_Mag + &
                            Inf_Nodes%Num_Vary(ii)
      end if

      ! B chi or B POS azimuth
      ii = Inf_Nodes%index_Bp
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Mag = Inf_Nodes%Num_Mag + &
                            Inf_Nodes%Num_Vary(ii)
      end if

      !
      ! Thermal variables

      ! Temperature
      ii = Inf_Nodes%index_T
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      ! X component of the velocity or vtrans
      ii = Inf_Nodes%index_vx
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      ! Y component of the velocity or v POS phi
      ii = Inf_Nodes%index_vy
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      ! Z component of the velocity or vlos
      ii = Inf_Nodes%index_vz
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      ! Microturbulent velocity
      ii = Inf_Nodes%index_vm
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      ! Gas pressure
      ii = Inf_Nodes%index_Pg
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Thermal = Inf_Nodes%Num_Thermal + &
                                Inf_Nodes%Num_Vary(ii)
      end if

      !
      ! Global variables

      ! Diffuse light
      ii = Inf_Nodes%index_f
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_glob = Inf_Nodes%Num_glob + &
                             Inf_Nodes%Num_Vary(ii)
      end if

      !
      ! Asymmetry variables

      ! Real J21
      ii = Inf_Nodes%index_J21R
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if

      ! Imaginary J21
      ii = Inf_Nodes%index_J21I
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if

      ! Real J22
      ii = Inf_Nodes%index_J22R
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if

      ! Imaginary J22
      ii = Inf_Nodes%index_J22I
      if (Inf_Nodes%Nodes_Flags(ii)) then
        Inf_Nodes%Num_Asymmetry = Inf_Nodes%Num_Asymmetry + &
                                  Inf_Nodes%Num_Vary(ii)
      end if

      ! Total number of nodes
      Inf_Nodes%Num_Fit = Inf_Nodes%Num_Mag + &
                          Inf_Nodes%Num_Thermal + &
                          Inf_Nodes%Num_glob + &
                          Inf_Nodes%Num_Asymmetry

      ! If no nodes to invert
      if (Inf_Nodes%Num_Fit.le.0) then

        ! Leave this absurd run
        umsg = 'There are no nodes to fit'
        urou = 'set_up_inv'
        call gabortedv

      end if ! No nodes

      !
      ! Check backtracking
      !

      ! If backtracking mehotd
      if (Input%LM_Method.eq.1) then

        ! If desparate backtracking mode
        if (Input%LM_Back_Mode.eq.1) then

          ! Check that all parameters are within the limits
          check = Input%LM_lam_big_test.le.Input%Lam_Range(1).or. &
                  Input%LM_lam_big_prove.le.Input%Lam_Range(1).or. &
                  Input%LM_lam_small_test.le.Input%Lam_Range(1).or. &
                  Input%LM_lam_small_prove.le.Input%Lam_Range(1).or. &
                  Input%LM_lam_big_test.ge.Input%Lam_Range(2).or. &
                  Input%LM_lam_big_prove.ge.Input%Lam_Range(2).or. &
                  Input%LM_lam_small_test.ge.Input%Lam_Range(2).or. &
                  Input%LM_lam_small_prove.ge.Input%Lam_Range(2).or. &
                  Input%LM_lam_small_prove.ge.Input%LM_lam_big_prove

          ! Wrong ranges
          if (check) then

            ! Issue error
            umsg = 'Found a consistency issue between the ranges '// &
                   'of the lambda parameter and the test and '// &
                   'prove parameters of the backtracking '// &
                   'desperate shake-up. Check that the "small" '// &
                   'variables are smaller than the corresponding '// &
                   '"big" variables and that all of them are '// &
                   'strictly between the limits for lambda'
            urou = 'set_up_inv'
            call gabortedv

          end if ! Correct ranges
        end if ! Desperate backtracking
      end if ! Backtracking


      !
      ! Show the whole configuration to the user
      !

      ! Only the global master
      if (gpid.eq.0) then

        ! Information message
        umsg = ' - Called the inversion module with the '// &
               'following parameters:'
        call verboseI(1)

        ! Verbosity level
        write(umsg,'(A,i1)') '   o Verbosity level:',vlevel
        call verboseI(1)

        ! Data file
        write(umsg,'(A,A)') '   o Input data file: ', &
                             trim(Input%Filename_ob)
        call verboseI(1)

        ! Mask file
        if (trim(Input%Inv_mask).ne.'NONE'.and. &
            trim(Input%Inv_init).ne.'INIT') then
          write(umsg,'(A,A)') '   o Mask file: ', &
                              trim(Input%Inv_mask)
          call verboseI(1)
        end if

        ! Type of inversion
        if (Input%Type_Inversion.eq.0) then
          write(umsg,'(A)') '   o Thermodynamic inversion'
        else if (Input%Type_Inversion.eq.1) then
          write(umsg,'(A)') '   o Magnetic inversion'
        else if (Input%Type_Inversion.eq.2) then
          write(umsg,'(A)') '   o Full inversion together'
        else if (Input%Type_Inversion.eq.3) then
          write(umsg,'(A)') '   o Full inversion sequential'
        else if (Input%Type_Inversion.eq.4) then
          write(umsg,'(A)') '   o Thermodynamic followed by '// &
                            'magnetic inversion'
        else if (Input%Type_Inversion.eq.5) then
          write(umsg,'(A)') '   o Thermodynamic, then magnetic, '// &
                            'and full inversion'
        else
          write(umsg,'(A)') '   o Inversion type not recognized'
        end if
        call verboseI(1)

        ! Stokes parameters weights
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
              if (Input%weight(5,ia).lt.1d99) then
                write(umsg,'(A,2(1x,f9.3),A,4(1x,es10.3))') &
                  '   o Stokes weights between range', &
                  Input%weight(4,ia),Input%weight(5,ia), &
                  ' nm :',Input%weight(0:3,ia)
              else
                write(umsg,'(A,4(1x,es10.3))') &
                  '   o Stokes weights :', Input%weight(0:3,ia)
              end if
              call verboseI(1)
            end do
          end if
        end if

        ! If initializing from previous solution
        if (trim(Input%Inv_init).ne.'NONE') then
          write(umsg,'(A,A)') '   o Restore from file ', &
                               Input%Inv_init
          call verboseI(1)
        end if

        ! If centered derivatives
        if (Input%centered) then
          write(umsg,'(A)') '   o Centered derivative for RF'
        else
          write(umsg,'(A)') '   o Non-centered derivative for RF'
        end if
        call verboseI(1)

        ! Maximum number of inversion iterations
        write(umsg,'(A,i5)') '   o Maximum iterations ', &
                             Input%Num_Iter
        call verboseI(1)

        ! Type of interpolation to build the model atmosphere
        if (Input%Interpolation.eq.0) then
          write(umsg,'(A)') '   o Linear interpolation'
        else if (Input%Interpolation.eq.1) then
          write(umsg,'(A)') '   o Quadratic bezier interpolation'
        else if (Input%Interpolation.eq.2) then
          write(umsg,'(A)') '   o Cubic bezier interpolation'
        end if
        call verboseI(1)

        ! Type of magnetic field vector
        if (Input%btype.eq.0) then
          write(umsg,'(A)') '   o Magnetic field in vertical ref.'
        else if (Input%btype.eq.1) then
          write(umsg,'(A)') '   o Magnetic field in LOS ref.'
        end if
        call verboseI(1)

        ! Type of velocity vector
        if (Input%vtype.eq.0) then
          write(umsg,'(A)') '   o Cartesian velocity'
        else if (Input%vtype.eq.1) then
          write(umsg,'(A)') '   o Velocity in LOS ref.'
        end if
        call verboseI(1)

        ! If correct the node positions to coincide with positions
        ! in the model atmosphere
        if (Input%Pos_Correction) then
          write(umsg,'(A)') '   o Correct node positions'
        else
          write(umsg,'(A)') '   o Do not correct node positions'
        end if
        call verboseI(1)


        ! For every variable
        do ii=1,Input%nvar

          ! If magnetic field strength or Blos
          if (ii.eq.Inf_Nodes%index_B) then

            ! Verbose
            if (Input%btype.eq.0) then
              write(umsg,'(A,i2,A)') '   o |B| (',ii,'):'
            else if (Input%btype.eq.1) then
              write(umsg,'(A,i2,A)') '   o Bpar (',ii,'): '
            end if

          ! If magnetic field inclination or Btrans
          else if (ii.eq.Inf_Nodes%index_Bt) then

            ! Verbose
            if (Input%btype.eq.0) then
              write(umsg,'(A,i2,A)') '   o B_incl (',ii,'):'
            else if (Input%btype.eq.1) then
              write(umsg,'(A,i2,A)') '   o Btrv (',ii,'): '
            end if

          ! If magnetic field azimuth or POS azimuth
          else if (ii.eq.Inf_Nodes%index_Bp) then

            ! Verbose
            if (Input%btype.eq.0) then
              write(umsg,'(A,i2,A)') '   o B_azim (',ii,'):'
            else if (Input%btype.eq.1) then
              write(umsg,'(A,i2,A)') '   o B_pos-azim (',ii,'):'
            end if

          ! If temperature
          else if (ii.eq.Inf_Nodes%index_T) then

            ! Verbose
            write(umsg,'(A,i2,A)') '   o Temperature: (',ii,')'

          ! If X component of the velocity or vtrans
          else if (ii.eq.Inf_Nodes%index_vx) then

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle

            ! Verbose
            if (Input%vtype.eq.0) then
              write(umsg,'(A,i2,A)') '   o X axis velocity (',ii,'):'
            else if (Input%vtype.eq.1) then
              write(umsg,'(A,i2,A)') '   o POS velocity (',ii,'):'
            end if

          ! If Y component of the velocity or POS azimuth
          else if (ii.eq.Inf_Nodes%index_vy) then

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle

            ! Verbose
            if (Input%vtype.eq.0) then
              write(umsg,'(A,i2,A)') '   o Y axis velocity (',ii,'):'
            else if (Input%vtype.eq.1) then
             write(umsg,'(A,i2,A)') '   o Velocity POS azim (',ii,'):'
            end if

          ! If vertical velocity or vlos
          else if (ii.eq.Inf_Nodes%index_vz) then

            ! Verbose
            if (Input%vtype.eq.0) then
             write(umsg,'(A,i2,A)') '   o Vertical velocity (',ii,'):'
            else if (Input%vtype.eq.1) then
              write(umsg,'(A,i2,A)') '   o LOS velocity (',ii,'):'
            end if

          ! If microturbulent velocity
          else if (ii.eq.Inf_Nodes%index_vm) then

            ! Verbose
            write(umsg,'(A,i2,A)') '   o Turbulent velocity (',ii,'):'

          ! If gas pressure
          else if (ii.eq.Inf_Nodes%index_Pg) then

            ! Verbose
            write(umsg,'(A,i2,A)') '   o Gas pressure (',ii,'):'

          ! If diffuse light
          else if (ii.eq.Inf_Nodes%index_f) then

            ! Verbose
            write(umsg,'(A,i2,A)') '   o Diffuse light (',ii,'):'

          ! If real J21
          else if (ii.eq.Inf_Nodes%index_J21R) then

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle

            ! Verbose
            write(umsg,'(A,i2,A)') '   o J21 real (',ii,'):'

          ! If imaginary J21
          else if (ii.eq.Inf_Nodes%index_J21I) then

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle

            ! Verbose
            write(umsg,'(A,i2,A)') '   o J21 imag (',ii,'):'

          ! If real J22
          else if (ii.eq.Inf_Nodes%index_J22R) then

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle

            ! Verbose
            write(umsg,'(A,i2,A)') '   o J22 real (',ii,'):'

          ! If imaginary J22
          else if (ii.eq.Inf_Nodes%index_J22I) then

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_flags(ii)) cycle

            ! Verbose
            write(umsg,'(A,i2,A)') '   o J22 imag (',ii,'):'

          ! Variable not recognized
          else

            ! Verbose no idea and skip rest
            write(umsg,'(A)') '   o ??'
            cycle

          end if
          call verboseI(1)

          ! Verbose if inverting
          if (Inf_Nodes%Nodes_Flags(ii)) then
            write(umsg,'(A)') '     + Inverting'
          else
            write(umsg,'(A)') '     + Fixed'
          end if
          call verboseI(1)

          ! If variable is pressure and doing hydrostatic Eq.
          if (ii.eq.Inf_Nodes%index_Pg.and.Input%hydroeq) then

            ! Verbose
            write(umsg,'(A)') '     + Hydrostatic equilibrium'
            call verboseI(1)

            ! If specified value
            if (Input%Pg_bound.gt.0d0) then

              ! Verbose
              write(umsg,'(A,es10.3)') '     + Boundary value ', &
                                       Input%Pg_bound
              call verboseI(1)

            end if ! Specified gas pressure

          ! Any other variable
          else

            ! Skip if not inverting
            if (.not.Inf_Nodes%Nodes_Flags(ii)) cycle

            ! If node value
            if (Inf_Nodes%Node_type(ii).eq.0) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have values'

            ! If node values with fixed last node
            else if (Inf_Nodes%Node_type(ii).eq.1) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have values '// &
                                'with last node fixed'

            ! If node value with fixed first node
            else if (Inf_Nodes%Node_type(ii).eq.2) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have values '// &
                                'with first node fixed'

            ! If node value with fixed extremes
            else if (Inf_Nodes%Node_type(ii).eq.3) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have values '// &
                                'with extreme nodes fixed'

            ! If correction
            else if (Inf_Nodes%Node_type(ii).eq.4) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have corrections'

            ! If correcting with fixed last node
            else if (Inf_Nodes%Node_type(ii).eq.5) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have corrections '// &
                                'with last node fixed'

            ! If correcting with fixed first node
            else if (Inf_Nodes%Node_type(ii).eq.6) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have corrections '// &
                                'with first node fixed'

            ! If correcting with fixed extreme nodes
            else if (Inf_Nodes%Node_type(ii).eq.7) then

              ! Verbose
              write(umsg,'(A)') '     + Nodes have corrections '// &
                                'with extreme nodes fixed'

            end if ! Type of node
            call verboseI(1)

            ! Number of nodes
            write(umsg,'(A,i3)') '     + Number of nodes ', &
                                 Inf_Nodes%Num_nodes(ii)
            call verboseI(1)

            ! If locations are already set
            if (allocated(Inf_Nodes%Node(ii)%H)) then

              ! Write number of nodes to message
              write(caux,'(i10)') Inf_Nodes%Num_nodes(ii)

              ! Find space
              kk = 0
              do jj=1,10
                if (caux(jj:jj).ne.' ') then
                  kk = jj
                  exit
                end if
              end do

              ! Move number to beginning to trim it in both directions
              ll = 10 - kk + 1
              do jj=1,ll
                caux(jj:jj) = caux(kk+jj-1:kk+jj-1)
              end do
              do jj=ll+1,10
                caux(jj:jj) = ' '
              end do

              ! Write locations in message and verbose it
              write(umsg,'(A,'//trim(caux)//'(1x,es10.3))') &
                 '     + Fixed positions ',Inf_Nodes%Node(ii)%H
              call verboseI(1)

            end if ! Locations already set

            ! Value scale
            write(umsg,'(A,2(1x,es10.3))') &
              '     + Scale ',Inf_Nodes%Scal(ii)
            call verboseI(1)

            ! Perturbation for response function calculation
            write(umsg,'(A,2(1x,es10.3))') &
              '     + Perturbation ',Inf_Nodes%Perturb(ii)
            call verboseI(1)

            ! Boundary values
            write(umsg,'(A,2(1x,es10.3))') &
              '     + Bound values ',Inf_Nodes%Node(ii)%Bounds
            call verboseI(1)

            ! For each special limit
            do jj=1,Inf_Nodes%Node(ii)%nebound

              ! Write the limit data
              write(umsg,'(A,2(1x,es10.3),A,2(1x,es10.3))') &
                '     + Bound values between', &
                Inf_Nodes%Node(ii)%ebound(1:2,jj),': ', &
                Inf_Nodes%Node(ii)%ebound(3:4,jj)
              call verboseI(1)

            end do ! special limits

          end if ! Pressure

          ! If regularizing this variable
          if (Inf_Nodes%Nodes_regul(ii)) then

            ! If average regularization
            if (Inf_Nodes%Indx_regul(ii).eq.1) then

              ! Verbose
              write(umsg,'(A,es10.3)') &
                '     + Mean regularization with weight ', &
                Inf_Nodes%Regul_weight(ii)

            ! If constant regularization
            else if (Inf_Nodes%Indx_regul(ii).eq.2) then

              ! Verbose
              write(umsg,'(A,es10.3)') &
                '     + Constant regularization with weight ', &
                Inf_Nodes%Regul_weight(ii)

            ! If first derivative regularization
            else if (Inf_Nodes%Indx_regul(ii).eq.3) then

              ! Verbose
              write(umsg,'(A,es10.3)') &
                '     + First derivative regularization with '// &
                'weight ',Inf_Nodes%Regul_weight(ii)

            ! If second derivative regularization
            else if (Inf_Nodes%Indx_regul(ii).eq.4) then

              ! Verbose
              write(umsg,'(A,es10.3)') &
                '     + Second derivative regularization with '// &
                'weight ',Inf_Nodes%Regul_weight(ii)

            end if ! Type of regularization

            ! Write message
            call verboseI(1)

          end if ! Regularization

          ! If extrapolating this variable
          if (Inf_Nodes%extrapolation(ii).gt.0) then

            if (Inf_Nodes%extrapolation(ii)/4.gt.0) then

              select case (Inf_Nodes%extrapolation(ii)/4)

                ! Zero
                case(1)

                  ! Verbose
                  write(umsg,'(A)') &
                      '     + Setting to zero toward the top'

                ! Constant
                case(2)

                  ! Verbose
                  write(umsg,'(A)') &
                      '     + Extending value toward the top'

                ! Linear
                case(3)

                  ! Verbose
                  write(umsg,'(A)') &
                      '     + Linear extrapolation toward the top'

              end select ! Top case

              ! Write message
              call verboseI(1)

            end if ! Extrapolating top

            if (mod(Inf_Nodes%extrapolation(ii),4).gt.0) then

              select case (mod(Inf_Nodes%extrapolation(ii),4))

                ! Zero
                case(1)

                  ! Verbose
                  write(umsg,'(A)') &
                      '     + Setting to zero toward the bottom'

                ! Constant
                case(2)

                  ! Verbose
                  write(umsg,'(A)') &
                      '     + Extending value toward the bottom'

                ! Linear
                case(3)

                  ! Verbose
                  write(umsg,'(A)') &
                      '     + Linear extrapolation toward the bottom'

              end select ! Top case

              ! Write message
              call verboseI(1)

            end if ! Extrapolating top
          end if ! Extrapolating variable

        end do ! Variables

        ! Total number of variable nodes
        write(umsg,'(A,i4)') '   o Total number of nodes ', &
                             Inf_Nodes%Num_Fit
        call verboseI(1)

        ! Total number of variable nodes for magnetic quantities if
        ! larger than 0
        if (Inf_Nodes%Num_Mag.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of magnetic '// &
                               'nodes ',Inf_Nodes%Num_Mag
          call verboseI(1)
        end if

        ! Total number of variable nodes for thermal quantities if
        ! larger than 0
        if (Inf_Nodes%Num_Thermal.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of thermal '// &
                               'nodes ',Inf_Nodes%Num_Thermal
          call verboseI(1)
        end if

        ! Total number of variable nodes for global quantities if
        ! larger than 0
        if (Inf_Nodes%Num_glob.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of global '// &
                               'nodes ',Inf_Nodes%Num_glob
          call verboseI(1)
        end if

        ! Total number of variable nodes for asymmetry quantities if
        ! larger than 0
        if (Inf_Nodes%Num_Asymmetry.gt.0) then
          write(umsg,'(A,i4)') '   o Total number of symmetry '// &
                               'nodes ',Inf_Nodes%Num_Asymmetry
          call verboseI(1)
        end if

        ! Limit for regularization contribution to merit function
        if (maxval(Inf_Nodes%Indx_regul).gt.0) then
          write(umsg,'(A,es10.3)') &
            '   o Regularization limit',Input%Regul_Limit
          call verboseI(1)
        end if

        ! Merit function threshold for solution
        write(umsg,'(A,es10.3)') '   o Chi2 threshold ', &
                                 Input%Threshold_chisq
        call verboseI(1)

        ! Merit function MRC to consider convergence
        write(umsg,'(A,es10.3)') '   o Chi2 MRC ', &
                                 Input%Chisq_fraction
        call verboseI(1)

        ! Type of SVD solution
        if (Input%SVD_type.eq.0) then
          write(umsg,'(A,es10.3)') '   o SVD method: tradiational'
        else if (Input%SVD_type.eq.2) then
          write(umsg,'(A,es10.3)') '   o SVD method: SIR-like'
        end if
        call verboseI(1)

        ! Threshold for SVD
        write(umsg,'(A,es10.3)') '   o SVD threshold ', &
                                 Input%Threshold_svd
        call verboseI(1)

        ! Initial model atmosphere
        write(umsg,'(A,A)') '   o Initial atmospheric model file: ', &
                             trim(Input%atmo)
        call verboseI(1)

        ! Number of nodes in model atmosphere during synthesis
        write(umsg,'(A,i4)') '   o Atmospheric nodes in '// &
                             'synthesis: ',Input%Atmo_Input
        call verboseI(1)

        ! Maximum step for SVD
        write(umsg,'(A,es10.3)') '   o SVD maximum step ', &
                                 Input%Max_Step
        call verboseI(1)

        ! If Hessian error
        if (Input%Err_type.eq.0) then

          ! Verbose
          write(umsg,'(A)') '   o Error from Hessian'
          call verboseI(1)

        ! If non-Hessian error
        else if (Input%Err_type.eq.1) then

          ! Verbose
          write(umsg,'(A)') '   o Error not from Hessian'
          call verboseI(1)

        ! If error from response function
        else if (Input%Err_type.eq.2) then

          ! Verbose
          write(umsg,'(A)') '   o Error from RF'
          call verboseI(1)

        ! If worst error
        else if (Input%Err_type.eq.3) then

          ! Verbose
          write(umsg,'(A)') '   o Error worst from Hessian or RF'
          call verboseI(1)

        end if ! Type of error

        ! If there is a PSF
        if (allocated(Input%lim_fwhm)) then

          ! If there is more than one range
          if (Input%lim_fwhm(1)%nn.gt.1) &
            write(umsg,'(A,i2)') '   o PSF number of ranges ', &
                                 Input%lim_fwhm(1)%nn

          ! For each range
          do ii=1,Input%lim_fwhm(1)%nn

            ! If Gaussian PSF
            if (Input%lim_fwhm(ii)%gaussian) then

              ! If actual wavelength limits
              if (Input%lim_fwhm(ii)%doub(2).lt.1d99) then

                ! Verbose limits and value
                write(umsg,'(A,es10.3," (",es10.3,",",es10.3,")")') &
                       '   o PSF FWHM ',Input%lim_fwhm(ii)%doub(3), &
                        Input%lim_fwhm(ii)%doub(1:2)

              ! Fake limits
              else

                ! Verbose value
                write(umsg,'(A,es10.3)') &
                       '   o PSF FWHM ',Input%lim_fwhm(ii)%doub(3)

              end if ! Type of limits

              ! Verbose
              call verboseI(1)

            ! No gaussian, so it is file
            else

              ! If actual wavelength limits
              if (Input%lim_fwhm(ii)%doub(2).lt.1d99) then

                ! Verbose limits and filename
                write(umsg,'(A,A," (",es10.3,",",es10.3,")")') &
                    '   o PSF from file', &
                     trim(Input%fwhm_fil(ii)%str), &
                     Input%lim_fwhm(ii)%doub(1:2)

              ! Fake limits
              else

                ! Verbose filename
                write(umsg,'(A,A)') &
                    '   o PSF from file',trim(Input%fwhm_fil(ii)%str)
              end if ! Type of limits

              ! Verbose
              call verboseI(1)

            end if ! Type of PSF

          end do ! List of PSF

        end if ! If there are PSF

        ! If initialize from last solution when calculating response
        ! functions
        if (Input%PopuInit) then
          write(umsg,'(A)') '   o Initialize RF with solution'
        else
          write(umsg,'(A)') '   o Do not initialize RF with solution'
        end if
        call verboseI(1)

        ! If initialize from last solution when calculating trials
        if (Input%trialinit) then
          if (Input%trialtpinit) then
            write(umsg,'(A)') '   o Initialize Trials with '// &
                                'solution if T and Pg are fixed'
          else
            write(umsg,'(A)') '   o Initialize Trials with '// &
                                'solution if thermal '// &
                                'parameters are fixed'
          end if
          call verboseI(1)
        end if

        ! If fractional polarization profiles
        if (Input%Fractional) then
          write(umsg,'(A)') '   o Use fractional Stokes parameters'
        else
          write(umsg,'(A)') '   o Use absolute Stokes parameters'
        end if
        call verboseI(1)

        ! Verbose the optical depth range
        write(umsg,'(A,2(1x,es10.3))') &
          '   o Range of taus ',Input%Tau_Range
        call verboseI(1)

        ! If using Broyden method
        if (Input%Broyden) then

          ! Verbose
          write(umsg,'(A)') '   o LM Broyden'
          call verboseI(1)

        end if ! Broyden method


        ! If tradiational LM
        if (Input%LM_Method.eq.0) then

          ! Verbose
          write(umsg,'(A,es10.3)') '   o Traditional LM'

        ! If backtracking LM
        else if (Input%LM_Method.eq.1) then

          ! Verbose
          write(umsg,'(A,es10.3)') '   o LM with backtracking'

        end if ! Type of LM
        call verboseI(1)

        ! If neglecting Stokes errors
        if (Input%Sigma_neglect) then

          ! Verbose
          write(umsg,'(A)') '   o Neglect sigma'
          call verboseI(1)

        end if ! Neglecting Stokes errors
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

      ! If cartesian
      if (Input%vtype.eq.0) then

        ! For vy
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

      end if ! Cartesian

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

      ! Magnetic field azimuth
      ii = Inf_Nodes%index_Bp

      ! If inverting magnetic field azimuth
      if (Inf_Nodes%Nodes_flags(ii)) then

        ! If full range
        if (abs(Inf_Nodes%Node(ii)%bounds(2) - &
                Inf_Nodes%Node(ii)%bounds(1)).ge.(2d0*PI-TINYA)) then

          ! Expand limits by pi
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

            ! Expand limits by pi
            Inf_Nodes%Node(ii)%ebound(4,jj) = &
                                  Inf_Nodes%Node(ii)%ebound(4,jj) + PI
            Inf_Nodes%Node(ii)%ebound(3,jj) = &
                                  Inf_Nodes%Node(ii)%ebound(3,jj) - PI

          end if ! Full azimuth range

        end do ! special ranges

      end if ! Inverting azimuth

      ! Velocity azimuth
      ii = Inf_Nodes%index_vy

      ! If inverting velocity azimuth
      if (Inf_Nodes%Nodes_flags(ii).and.Input%vtype.eq.1) then

        ! If full range
        if (abs(Inf_Nodes%Node(ii)%bounds(2) - &
                Inf_Nodes%Node(ii)%bounds(1)).ge.(2d0*PI-TINYA)) then

          ! Expand limits by pi
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

            ! Expand limits by pi
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

      !> Determine the limits of T, v, and B, as well as if the model
      !! is static/dynamic or un/magnetized. A possible outcome is
      !! that this needs to be checked pixelwise\n
      !!       Input(Input_class): Structure with configuration data\n
      !!   Inf_Nodes(Nodes_class): Structure with inversion node
      !!                           data\n
      !!      Atmo_in(Atmo_class): Structure with atmospheric data
      !!                           read from model atmosphere\n
      !!  Bfield_in(Bfield_class): Structure with magnetic field data
      !!                           read from the input\n
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
      type(Input_class), intent(inout):: Input
      type(Nodes_class), intent(in):: Inf_Nodes
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

      ! Temperature
      ii = Inf_Nodes%index_T

      ! If inverting temperature
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

        ! Atmosphere not already set
        else

          ! Need to update
          update_tlim = .True.

        end if ! Type of model input
      end if ! Inverting T or not

      ! Initialize
      maxvx = 0d0
      maxvy = 0d0
      maxvz = 0d0

      ! Inverting any velocity with values
      if ((Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vx).or. &
           Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vy).or. &
           Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vz)).and. &
           (Input%vtype.eq.0.or. &
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vx).or. &
            Inf_Nodes%Nodes_Flags(Inf_Nodes%index_vz))) then

        ! vx
        ii=Inf_Nodes%index_vx

        ! Inverting vx or vtrans
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

        ! Inverting vy
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

        ! Inverting vz or vlos
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

        ! Atmosphere not set yet
        else

          ! Need to update
          update_vlim = .True.

        end if ! Type of model input
      end if ! Inverting velocity or not

      ! If to update later
      if (update_vlim) then

        ! Asumme dynamic to start with
        dyn = .True.
        Input%static = .False.

      ! No need to update
      else

        ! Combine
        maxvz = sqrt(maxvx*maxvx + maxvy*maxvy + maxvz*maxvz)

        ! If maximum velocity if large enough
        if (maxvz.gt.TINYVEL) then

          ! Set dynamic
          dyn = .True.
          Input%maxV = maxvz
          Input%static = .False.

        ! Small maximum velocity
        else

          ! Set static
          dyn = .False.
          Input%maxV = 0d0
          Input%static = .True.

        end if ! Velocity magnitude
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

        ! If inverting B
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_B).or. &
            (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Bt).and. &
             Input%Btype.eq.1)) then

          ! Bstrength or BLOS
          ii = Inf_Nodes%index_B

          ! Get maximum value
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

            ! Check maximum value
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

          ! Atmosphere not set
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

        ! If small field
        if (maxB.le.TINYB) then

          ! Set unmagnetized
          Input%unmagnetized = .True.

        ! If significant field
        else

          ! Set not unmagnetized
          Input%unmagnetized = .False.

        end if ! Magnetic field value
      end if ! To update later

      end subroutine set_up_limits

!#####################################################################
!#####################################################################
!#####################################################################

      !> Make a first guess of the magnetic field based on WFA\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data
      subroutine guess_polarity(Input,GeomI,Inf_Stokes,Sol,Bfield)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(inout):: GeomI
      type(Stokes_class), intent(in):: Inf_Stokes
      type(Solution_class), intent(in):: Sol
      type(Bfield_class), intent(inout):: Bfield

      ! Local

      logical:: homogeneous

      integer:: iz,il0,il1,il
      integer, dimension(1):: ill

      double precision:: sig,Blos
      double precision, dimension(:), allocatable:: Ider
      double precision, dimension(:), allocatable:: aBl,aBtr,aBa
      double precision, dimension(:), allocatable:: aBm,aBth,aBp


      !
      ! Check if homogeneous
      !

      ! Initialize
      homogeneous = .True.

      ! Check heights
      do iz=2,nz

        ! Check difference strength
        if (abs(real(Bfield%Bstrength(iz) - Bfield%Bstrength(1))).gt.&
            TINYSP) then

          ! No homogeneous
          homogeneous = .False.
          exit

        end if ! Check 

        ! Check difference theta
        if (abs(real(Bfield%Btheta(iz) - Bfield%Btheta(1))).gt.&
            TINYSP) then

          ! No homogeneous
          homogeneous = .False.
          exit

        end if ! Check 

        ! Check difference phi
        if (abs(real(Bfield%Bphi(iz) - Bfield%Bphi(1))).gt.&
            TINYSP) then

          ! No homogeneous
          homogeneous = .False.
          exit

        end if ! Check 

      end do ! Rest of heights

      ! If no homogeneous, then there was intent in the input
      if (.not.homogeneous) then

        ! Verbose
        write(umsg, '(A)') ' - The input says to guess the '// &
          'polarity, but the input magnetic field is not '// &
          'homogeneous'
        call verboseI(2)
        return

      end if

      !
      ! If here, that means that the magnetic field is homogeneous
      ! and we can proceed with our guess
      !

      ! Find limits in omega
      ill = minloc(abs(Sol%omega_input - Input%gp_l))
      il0 = ill(1)
      ill = minloc(abs(Sol%omega_input - Input%gp_r))
      il1 = ill(1)

      ! Move indexes
      if (Sol%omega_input(il0).lt.Input%gp_l) il0 = il0 + 1
      if (Sol%omega_input(il1).gt.Input%gp_r) il1 = il1 + 1

      ! Sanity
      if (il1 - il0.lt.3) then
        ! Verbose
        write(umsg, '(A)') ' - Not enough wavelengths to guess '// &
          'the polarity'
        call verboseI(3)
        return
      end if
 
      ! Allocate derivative
      allocate(Ider(il0:il1))

      ! First point derivative
      Ider(il0) = 0.5*(Inf_Stokes%Stokes_Ob(0,il0+1) - &
                       Inf_Stokes%Stokes_Ob(0,il0))/ &
                  (Sol%omega_input(il0+1) - Sol%omega_input(il0))

      ! Intermediate points
      do il=il0+1,il1-1

        ! Derivative
        Ider(il) = 0.5*(Inf_Stokes%Stokes_Ob(0,il+1) - &
                        Inf_Stokes%Stokes_Ob(0,il-1))/ &
                   (Sol%omega_input(il+1) - Sol%omega_input(il-1))

      end do

      ! Last point derivative
      Ider(il1) = 0.5*(Inf_Stokes%Stokes_Ob(0,il1-1) - &
                       Inf_Stokes%Stokes_Ob(0,il1))/ &
                  (Sol%omega_input(il1-1) - Sol%omega_input(il1))

      ! If guessing strength as well
      if (Input%gp_g.gt.-1d85) then

        ! Scale derivative
        ! 1d1 = 1d-1 derivate in A * (1e1)^2 wavelength in A
        Ider = Ider*1d1*Input%gp_w*Input%gp_w*Input%gp_g

        ! If fractional, scale by intensity
        if (Sol%Fractional) &
          Ider = Ider/Inf_Stokes%Stokes_Ob(0,il0:il1)

        ! Guess Blos
        Blos = -1d0*(1d0/wfac)*sum(Ider* &
                                   Inf_Stokes%Stokes_ob(3,il0:il1))/ &
                               sum(Ider*Ider)

        ! Allocate auxiliars
        allocate(aBl(1),aBtr(1),aBa(1),aBm(1),aBth(1),aBp(1))

        ! Get values from input
        aBm(1) = Bfield%Bstrength(1)
        aBth(1) = Bfield%Btheta(1)
        aBp(1) = Bfield%Bphi(1)

        ! Transform into LOS
        call B2Blos(1,GeomI%L_mu(1),GeomI%L_phi(1), &
                    aBm,aBth,aBp,aBl,aBtr,aBa)

        ! Change Blos
        aBl(1) = Blos

        ! Transform back
        call Bconversion(1,GeomI%L_mu(1),GeomI%L_phi(1), &
                         aBl,aBtr,aBa,aBm,aBth,aBp)

        ! Copy
        Bfield%Bstrength = aBm(1)
        Bfield%Btheta = aBth(1)
        Bfield%Bphi = aBp(1)

        ! Deallocate auxialiars
        deallocate(aBl,aBtr,aBa,aBm,aBth,aBp)

      ! If only guessing polarity
      else

        ! Sign
        sig = -1d0*sign(1d0,sum(Ider/Inf_Stokes%Stokes_Ob(3,il0:il1)))

        ! If positive polarity and negative sign
        if (Bfield%Btheta(1).lt.0.5*PI.and.sig.lt.0d0) then

          ! Switch
          Bfield%Btheta = PI-Bfield%Btheta

        ! If positive polarity and negative sign
        else if (Bfield%Btheta(1).gt.0.5*PI.and.sig.gt.0d0) then

          ! Switch
          Bfield%Btheta = PI-Bfield%Btheta

        end if ! Non-coincident polarities
      end if ! What are we guessing

      end subroutine guess_polarity

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the location, values, and errors of the nodes\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Input(Input_class): Structure with configuration data\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data
      subroutine Init_Nodes(Atmo,Input,Inf_Nodes)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input
      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local

      integer:: i, j, k

      double precision:: Delta_Tau, Mini, H_min, H_max
      double precision, dimension(:), allocatable:: TAU_LG

      !
      ! Get height limits from model atmosphere
      H_min = log10(Atmo%z(nZ))
      H_max = log10(Atmo%z(1))

      !
      ! For each variable
      do i=1,Input%nvar

        ! If there are nodes and inverting them
        if (Inf_Nodes%Num_Nodes(i).gt.0.and. &
            Inf_Nodes%Nodes_Flags(i)) then

          ! If values are allocated
          if (allocated(Inf_Nodes%Node(i)%Var)) then

            ! Free them
            MRAMc = MRAMc - 1d-6*sizeof(Inf_Nodes%Node(i)%Var)
            deallocate(Inf_Nodes%Node(i)%Var)

          end if ! Allocated values

          ! If errors are allocated
          if (allocated(Inf_Nodes%Node(i)%Errors)) then

            ! Free them
            MRAMc = MRAMc - 1d-6*sizeof(Inf_Nodes%Node(i)%Errors)
            deallocate(Inf_Nodes%Node(i)%Errors)

          end if ! Allcoated errors

          ! Allocate values and errors
          allocate(Inf_Nodes%Node(i)%Var(Inf_Nodes%Num_Nodes(i)))
          allocate(Inf_Nodes%Node(i)%Errors(Inf_Nodes%Num_Nodes(i)))
          MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Node(i)%Var)
          MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Node(i)%Errors)

          ! If variable is gas pressure and using hydrostatic Eq.
          if ((i.eq.Inf_Nodes%index_Pg).and.Inf_Nodes%hydroeq) then

            ! If allocated input locations
            if (allocated(Inf_Nodes%Node(i)%H)) then

              ! Free locations
              MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Node(i)%H)
              deallocate(Inf_Nodes%Node(i)%H)

            end if ! Allocated locations

            ! Allocate new array of locations
            allocate(Inf_Nodes%Node(i)%H(Inf_Nodes%Num_Nodes(i)))
            MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Node(i)%H)

            ! Set in upper boundary
            Inf_Nodes%Node(i)%H(1) = Input%Tau_range(1)

          ! If variable is diffuse light
          else if (i.eq.Inf_Nodes%index_f) then

            ! If allocated input locations
            if (allocated(Inf_Nodes%Node(i)%H)) then

              ! Free locations
              MRAMc = MRAMc - 1d-6*sizeof(Inf_Nodes%Node(i)%H)
              deallocate(Inf_Nodes%Node(i)%H)

            end if ! Allocated locations

            ! Allocate new array of locations
            allocate(Inf_Nodes%Node(i)%H(Inf_Nodes%Num_Nodes(i)))
            MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Node(i)%H)

            ! Set in upper boundary
            Inf_Nodes%Node(i)%H(1) = Input%Tau_range(1)

          ! No stratification specified for this variable
          else if (.not.allocated(Inf_Nodes%Node(i)%H)) then

            ! Set equally spaced nodes between the stablished limits
            call Locate_Nodes(H_Min,H_Max, &
                              Inf_Nodes%Num_Nodes(i), &
                              Inf_Nodes%Node(i)%H)

          end if ! Input stratification

        ! Not inverting this variable
        else

          ! Flag and set to 0
          Inf_Nodes%Nodes_flags(i) = .False.
          Inf_Nodes%Num_Nodes(i) = 0

        end if ! Inverting variable

      end do ! Variables

      ! If to correct position
      if (Inf_Nodes%Pos_Correction) then

        ! Allocate auxiliar optical depth scale
        allocate(TAU_LG(Atmo%nZ))
        TAU_LG = log10(Atmo%z)

        ! For each variable
        do i=1,Input%nvar

          ! Skip if not inverting this variable
          if (Inf_Nodes%Num_Nodes(i).le.0.or. &
              .not.Inf_Nodes%Nodes_Flags(i)) cycle

          ! Skip diffuse light or gas pressure if hydrostatic Eq.
          if (Inf_Nodes%index_f.eq.i.or. &
              (Inf_Nodes%index_Pg.eq.i.and.Inf_Nodes%hydroeq)) cycle

          ! If allocated Tau_Indx array
          if (allocated(Inf_Nodes%Node(i)%Tau_Indx)) then

            ! Free array
            MRAMc = MRAMc - 1d-6*sizeof(Inf_Nodes%Node(i)%Tau_Indx)
            deallocate(Inf_Nodes%Node(i)%Tau_Indx)

          end if ! Allocated Tau_Indx

          ! Allocate Tau_Indx
          allocate(Inf_Nodes%Node(i)%Tau_Indx( &
                                            Inf_Nodes%Num_Nodes(i)))
          MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Node(i)%Tau_Indx)

          ! For each node for this variable
          do j=1,Inf_Nodes%Num_Nodes(i)

            ! Initialize minimum distance and its index
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

            ! Set node position at closest atmosphere node
            Inf_Nodes%Node(i)%H(j) = &
                               TAU_LG(Inf_Nodes%Node(i)%Tau_Indx(j))

          end do ! For each node
        end do ! Variables

        ! Deallocate auxiliar
        deallocate(TAU_LG)

      end if ! Correct position

      ! If indexing is allocated
      if (allocated(Inf_Nodes%Inf_Inv)) then

        ! Free indexing
        MRAMc = MRAMc - 1d-6*sizeof(Inf_Nodes%Inf_Inv)
        deallocate(Inf_Nodes%Inf_Inv)

      end if ! Allocated indexing

      ! Allocate indexing
      allocate(Inf_Nodes%Inf_Inv(2,Inf_Nodes%Num_Fit))
      MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes%Inf_Inv)

      ! Initialize counter
      j = 1

      ! For each variable
      do k=1,Input%nvar

        ! If inverting this variable
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

      !> Generate an array with equally spaced nodes within the given
      !! limits\n
      !!      H_min(double): Minimum position of the nodes\n
      !!      H_max(double): Maximum position of the nodes\n
      !!       Num(integer): Number of nodes\n
      !!  Height(double(:)): Node positions
      subroutine Locate_Nodes(H_min,H_max,Num,Height)

      ! I/O

      integer, intent(in):: Num
      double precision, intent(in):: H_min, H_max
      double precision, dimension(:), &
                        allocatable, intent(out):: Height

      ! Local

      integer:: i

      double precision:: dx


      ! Allocate positions (this variable is Inf_Nodes%Node(i)%H)
      allocate(Height(Num))
      MRAMc = MRAMc + 1d-6*sizeof(Height)

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

      !> Create a model atmosphere by interpolating a given model into
      !! a given optical depth stratification\n
      !!      Atmo_in(Atmo_class): Structure with atmospheric data
      !!                           read from model atmosphere\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Bfield_in(Bfield_class): Structure with magnetic field data
      !!                           read from the input\n
      !!     Bfield(Bfield_class): Structure with magnetic field
      !!                           data\n
      !!             z(double(:)): Final optical depth stratification
      subroutine Atmo_Stratify(Atmo_in,Atmo,Bfield_in,Bfield,z)

      ! I/O

      type(Atmo_class), intent(in):: Atmo_in
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(in):: Bfield_in
      type(Bfield_class), intent(inout):: Bfield
      double precision, dimension(:), intent(in):: z

      ! Local

      integer:: i

      double precision, dimension(Atmo_in%nz):: LTAUI


      ! Get logarithm of input optical depth
      LTAUI = log10(Atmo_in%z)

      ! Set optical depth scale, 500 nm ref, solar gravity, and
      ! electron density type of new model
      Atmo%scal = 'T'
      Atmo%typo = Atmo_in%typo
      Atmo%tfreq = Atmo_in%tfreq
      Atmo%logg = Atmo_in%logg

      ! Set global variables
      ztau = .True.

      ! Allocate variables and nullify pointers
      allocate(Atmo%z(nZ),Atmo%T(nZ),Atmo%ne(nZ),Atmo%Pg(nZ))
      allocate(Atmo%vmi(nZ),Atmo%vx(nZ),Atmo%vy(nZ),Atmo%vz(nZ))
      allocate(Atmo%nh(nZ,6),Atmo%nht(nZ),Atmo%nha(nZ),Atmo%nhm(nZ))
      allocate(Atmo%nhe(1,1),Atmo%zeros(nZ))
      allocate(Bfield%Bstrength(nZ))
      allocate(Bfield%Btheta(nZ))
      allocate(Bfield%Bphi(nZ))
      nullify(Atmo%Bx,Atmo%By,Atmo%Bz)
      nullify(Atmo%vxa,Atmo%vya,Atmo%vza)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%z)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%T)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%ne)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vmi)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vx)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vy)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vz)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nh)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nht)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nha)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhm)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%zeros)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bstrength)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Btheta)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bphi)

      ! Type of model allocation
      Atmo%alloc_a = .True.
      Atmo%alloc_b = .False.

      ! Copy rest of variables in model, i.e., partition function
      ! and abundance tabulations
      Atmo%NT = Atmo_in%NT
      Atmo%nele = Atmo_in%nele
      allocate(Atmo%ele(Atmo%nele),Atmo%pT(Atmo%NT))
      allocate(Atmo%abund(size(Atmo_in%abund)))
      Atmo%ele = Atmo_in%ele
      Atmo%abund = Atmo_in%abund
      Atmo%pT = Atmo_in%pT

      ! Count memory in these arrays
      do i=1,size(Atmo%ele)
        if (allocated(Atmo%ele(i)%Ei)) &
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%ele(i)%Ei)
        if (allocated(Atmo%ele(i)%pf)) &
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%ele(i)%pf)
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%ele(i))
      end do
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%pT)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%abund)

      ! Initialize everything to zero
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
                  Atmo%z, Atmo%T, Atmo%nZ, 2, 10)

      ! Pgas
      call Intpol(LTAUI, Atmo_in%Pg, Atmo_in%nZ, &
                  Atmo%z, Atmo%Pg, Atmo%nZ, 2, 10)

      ! Micro
      call Intpol(LTAUI, Atmo_in%vmi, Atmo_in%nZ, &
                  Atmo%z, Atmo%vmi, Atmo%nZ, 2, 10)

      ! vx
      call Intpol(LTAUI, Atmo_in%vx, Atmo_in%nZ, &
                  Atmo%z, Atmo%vx, Atmo%nZ, 2, 10)

      ! vy
      call Intpol(LTAUI, Atmo_in%vy, Atmo_in%nZ, &
                  Atmo%z, Atmo%vy, Atmo%nZ, 2, 10)

      ! vz
      call Intpol(LTAUI, Atmo_in%vz, Atmo_in%nZ, &
                  Atmo%z, Atmo%vz, Atmo%nZ, 2, 10)

      ! If magnetic field
      if (maxval(Bfield_in%Bstrength).gt.TINYB) then

        ! B
        call Intpol(LTAUI, Bfield_in%Bstrength, Atmo_in%nZ, &
                    Atmo%z, Bfield%Bstrength, Atmo%nZ, 2, 10)
        ! Btheta
        call Intpol(LTAUI, Bfield_in%Btheta, Atmo_in%nZ, &
                    Atmo%z, Bfield%Btheta, Atmo%nZ, 2, 10)
        ! Bphi
        call Intpol(LTAUI, Bfield_in%Bphi, Atmo_in%nZ, &
                    Atmo%z, Bfield%Bphi, Atmo%nZ, 2, 10)

        !
        ! Sanity
        !

        ! For each height
        do i=1,Atmo%nz

          ! No field
          if (Bfield%Bstrength(i).le.TINYB) then

            ! Just zero
            Bfield%Bstrength(i) = 0d0
            Bfield%Btheta(i) = 0d0
            Bfield%Bphi(i) = 0d0

          ! Check angles
          else

            ! Theta between 0 and pi
            if (Bfield%Btheta(i).lt.0d0) Bfield%Btheta(i) = 0d0
            if (Bfield%Btheta(i).gt.PI) Bfield%Btheta(i) = PI

            ! Azimuth between 0 and 2pi
            if (Bfield%Bphi(i).lt.0d0) Bfield%Bphi(i) = &
                                                    Bfield%Bphi(i) + &
                                                    2d0*PI
            if (Bfield%Bphi(i).gt.2d0*PI) Bfield%Bphi(i) = &
                                                    Bfield%Bphi(i) - &
                                                    2d0*PI

          end if ! Bfield value

        end do ! Heights

      ! No field
      else

        ! Just zero
        Bfield%Bstrength = 0d0
        Bfield%Btheta = 0d0
        Bfield%Bphi = 0d0

      end if ! Field or no field

      ! If JKQ in input
      if (allocated(Atmo_in%JKQin)) then

        ! If there is J21R
        if (maxval(abs(Atmo_in%JKQin(4*Atmo_in%nz+1: &
                                     5*Atmo_in%nz))).gt.0d0) &

          ! Interpolate
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(4*Atmo_in%nz+1:5*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(4*Atmo%nz+1:5*Atmo%nz), &
                      Atmo%nZ, 2, 10)

        ! If there is J21I
        if (maxval(abs(Atmo_in%JKQin(5*Atmo_in%nz+1: &
                                     6*Atmo_in%nz))).gt.0d0) &

          ! Interpolate
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(5*Atmo_in%nz+1:6*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(5*Atmo%nz+1:6*Atmo%nz), &
                      Atmo%nZ, 2, 10)

        ! If there is J22R
        if (maxval(abs(Atmo_in%JKQin(6*Atmo_in%nz+1: &
                                     7*Atmo_in%nz))).gt.0d0) &

          ! Interpolate
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(6*Atmo_in%nz+1:7*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(6*Atmo%nz+1:7*Atmo%nz), &
                      Atmo%nZ, 2, 10)

        ! If there is J22I
        if (maxval(abs(Atmo_in%JKQin(7*Atmo_in%nz+1: &
                                     8*Atmo_in%nz))).gt.0d0) &

          ! Interpolate
          call Intpol(LTAUI, &
                      Atmo_in%JKQin(7*Atmo_in%nz+1:8*Atmo_in%nz), &
                      Atmo_in%nZ, Atmo%z, &
                      Atmo%JKQin(7*Atmo%nz+1:8*Atmo%nz), &
                      Atmo%nZ, 2, 10)

      end if ! Ad-hoc JKQ

      ! Make the optical depth scale linear
      Atmo%z = 10d0**Atmo%z

      return

      end subroutine Atmo_Stratify

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine node values from the model atmosphere for a given
      !! variable\n
      !!          tau(double(:)): Log10 optical depth scale\n
      !!          var(double(:)): Stratificaiton of a variable\n
      !!              nn(iteger): Dimension of tau and var\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!            indx(iteger): Index of current variable
      subroutine set_nodes(tau,var,nn,Inf_Nodes,indx)

      ! I/O

      type(Nodes_class), intent(inout):: Inf_Nodes
      integer, intent(in):: nn, indx
      double precision, dimension(:), intent(in):: tau,var

      ! Local

      character(3):: length
      character(30):: fmt

      integer::j,i0,nnode
      double precision:: ff


      ! If inverting variable
      if (Inf_Nodes%Nodes_Flags(indx)) then

        ! If by value
        if (Inf_Nodes%Node_Type(indx).le.3) then

          ! Interpolate into nodes
          call Intpol(tau, var, nn, &
                      Inf_Nodes%Node(indx)%H, &
                      Inf_Nodes%Node(indx)%Var, &
                      Inf_Nodes%Num_Nodes(indx), 2, 10)

          ! Check boundaries
          call CheckBounds(Inf_Nodes%Node(indx), &
                           Inf_Nodes%Num_Nodes(indx))

        ! If by correction
        else

          ! Initialize
          Inf_Nodes%Node(indx)%Var = 0d0

        end if ! Type of node

        ! If Master
        if (pid.eq.0) then

          ! Verbose
          write(umsg, '(A,i2)') "   Parameter index = ",indx
          call verboseI(3)

          ! Small enough nodes for a single line
          if (Inf_Nodes%Num_Nodes(indx).lt.33) then

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

          ! Too many nodes
          else

            ! If velocity
            if (indx.eq.Inf_Nodes%index_vz.or. &
                indx.eq.Inf_Nodes%index_vx.or. &
                (indx.eq.Inf_Nodes%index_vy.and. &
                 Inf_Nodes%vtype.eq.0).or. &
                indx.eq.Inf_Nodes%index_vm) then

              ! Multiplicative factor
              ff = 1d6*c

            ! No velocity
            else

              ! Multiplicative factor
              ff = 1d0

            end if

            ! Save
            nnode = Inf_Nodes%Num_Nodes(indx)
            i0 = 1

            ! Get format
            write(length, "(i3)") 32
            fmt = '(A,'//trim(adjustl(length))//'es15.4)'
            fmt = trim(adjustl(fmt))

            ! Verbose positions
            write(umsg, FMT=fmt) "   Position: ", &
              (Inf_Nodes%Node(indx)%H(j), j=1, 32)
            call verboseI(3)

            ! Verbose values
            write(umsg, FMT=fmt) "     Values: ", &
             (Inf_Nodes%Node(indx)%Var(j)*ff, j=1, 32)
            call verboseI(3)

            ! Update
            i0 = i0 + 32
            nnode = nnode - 32

            ! If still more
            do while (nnode.gt.32)

              ! Get format
              write(length, "(i3)") 32
              fmt = '(A,'//trim(adjustl(length))//'es15.4)'
              fmt = trim(adjustl(fmt))

              ! Verbose positions
              write(umsg, FMT=fmt) "   Position: ", &
                (Inf_Nodes%Node(indx)%H(j), j=i0, i0+32-1)
              call verboseI(3)

              ! Verbose values
              write(umsg, FMT=fmt) "     Values: ", &
               (Inf_Nodes%Node(indx)%Var(j)*ff, j=i0, i0+32-1)
              call verboseI(3)

              ! Update
              i0 = i0 + 32
              nnode = nnode - 32

            end do

            ! If remaining
            if (nnode.gt.0) then

              ! Get format
              write(length, "(i3)") nnode
              fmt = '(A,'//trim(adjustl(length))//'es15.4)'
              fmt = trim(adjustl(fmt))

              ! Verbose positions
              write(umsg, FMT=fmt) "   Position: ", &
                (Inf_Nodes%Node(indx)%H(j), j=i0, i0+nnode-1)
              call verboseI(3)

              ! Verbose values
              write(umsg, FMT=fmt) "     Values: ", &
               (Inf_Nodes%Node(indx)%Var(j)*ff, j=i0, i0+nnode-1)
              call verboseI(3)

            end if ! There are nodes to write
          end if ! Number of nodes
        end if ! Master
      end if ! Inverting variable

      end subroutine set_nodes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Change the initial value of nodes in the transversal B or v
      !! variables when initializing from other inversion and they are
      !! too small\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!            indx(iteger): Index of current variable\n
      !!             val(double): Value to set the node to
      subroutine re_set_nodes(Inf_Nodes,indx,val)

      ! I/O

      type(Nodes_class), intent(inout):: Inf_Nodes
      integer, intent(in):: indx
      double precision, intent(in):: val

      ! Local

      character(3):: length
      character(30):: fmt

      integer::j


      ! If not inverting this variable, skip
      if (.not.Inf_Nodes%Nodes_Flags(indx)) return

      ! If not set by value, skip
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

      !> Initialize the node values from a given model atmosphere\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Bfield(Bfield_class): Structure the magnetic field data\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      subroutine Initialize_Nodes(Atmo,Bfield,Inf_Nodes)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local

      integer:: ii

      double precision, dimension(Atmo%nz):: TAU_LG


      ! Get logarithmic scale for optical depth
      TAU_LG = log10(Atmo%z)


      !
      ! If inverting any magnetic field node
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

        end if ! Master

        ! Type of magnetic field
        select case (Inf_Nodes%Btype)

          ! In vertical
          case(0)

            ! Bstrength
            ii = Inf_Nodes%index_B

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Bstrength, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Btheta
            ii = Inf_Nodes%index_Bt

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Btheta, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Bphi
            ii = Inf_Nodes%index_Bp

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Bphi, &
                           Atmo%nz,Inf_Nodes,ii)

          ! In LOS
          case(1)

            ! Blos
            ii = Inf_Nodes%index_B

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Blos, &
                           Atmo%nz,Inf_Nodes,ii)

            ! Btrans
            ii = Inf_Nodes%index_Bt

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Bpos, &
                           Atmo%nz,Inf_Nodes,ii)

            ! B POS phi
            ii = Inf_Nodes%index_Bp

            ! Get nodes
            call set_nodes(TAU_LG,Bfield%Azimuth, &
                           Atmo%nz,Inf_Nodes,ii)

        end select ! Type of magnetic field vector

      end if ! Magnetic nodes


      !
      ! If inverting any thermal node
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

            ! POS azimuth velocity
            ii = Inf_Nodes%index_vy

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vphi,Atmo%nz, &
                           Inf_Nodes,ii)

            ! LOS velocity
            ii = Inf_Nodes%index_vz

            ! Get nodes
            call set_nodes(TAU_LG,Atmo%vlos,Atmo%nz, &
                           Inf_Nodes,ii)

        end select ! Type of velocity vector

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

        ! If inverting and master
        if (Inf_Nodes%Nodes_flags(ii)) then

          ! Set from atmosphere
          Inf_Nodes%Node(ii)%var(1) = Atmo%f_diff

          ! If Master
          if (pid.eq.0) then

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

          end if ! Master
        end if ! Inverting
      end if ! Inverting global variable


      !
      ! Asymmetry variables
      !
      if (Inf_Nodes%Num_asymmetry.gt.0) then

        ! Master verbose
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

      end if ! Inverting ad-hoc asymmetry variables

      end subroutine Initialize_Nodes

!#####################################################################
!#####################################################################
!#####################################################################

      end module initinv_mod
