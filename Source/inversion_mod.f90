      !> Inversion module
      module inversion_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/22/2023
!  Last version:
!     11/29/2023 V3.1.10
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/29/2023:   V3.1.10 - Masked pixels do not change their
!                             model atmosphere now (TdPA)
!
!     11/24/2023:    V3.1.9 - Added inversion argument to pass to
!                             LMFIT (TdPA)
!
!     10/04/2023:    V3.1.8 - Changed priority of the current pixel
!                             verbosity (TdPA)
!                           - Added saving argument to Fit and
!                             LMfit call (TdPA)
!                           - Added message between the two cycles
!                             in sequential mode (TdPA)
!                           - Added the 'sequential magnetic' type of
!                             inversion (TdPA)
!
!     09/28/2023:    V3.1.7 - When changing the value of the Bt, Bp,
!                             vx, or vy nodes, call re_set_nodes and
!                             added verbosity to indicate it (TdPA)
!
!     09/08/2023:    V3.1.6 - Added extra lines in verbosity to better
!                             distinguish between pixels in the
!                             verbosity file (TdPA)
!                           - Verbosity update (TdPA)
!
!     08/18/2023:    V3.1.5 - Error in "if" nesting did not allow for
!                             restoring previous cycles (TdPA)
!
!     08/11/2023:    V3.1.4 - Added verbosity about current pixel
!                             to be inverted (TdPA)
!
!     08/11/2023:    V3.1.3 - Instead of counting number of Stokes
!                             parameters with weights, count degree
!                             of freedom of the data accounting for
!                             zero weights (TdPA)
!
!     07/31/2023:    V3.1.2 - Change the verbosity level in the
!                             inversion (HL)
!                           - Revise the initial values of vx and vy
!                             variables (HL)
!
!     07/06/2023:    V3.1.1 - Changed arguments in Atmo_Stratify call
!                             to adjust to its changes (TdPA)
!
!     07/03/2023:    V3.1.0 - The initial model atmosphere is fully
!                             determined at entry. This implies
!                             a significant rewrite of the logic
!                             in the inversion routine (TdPA)
!                           - No need to rename files anymore (TdPA)
!                           - Results are writen in lmfit and not
!                             in fit (TdPA)
!
!     06/12/2023:    V3.0.4 - Rename the variable Inf_File (HL)
!
!     05/16/2023:    V3.0.3 - Added running index in Inf_File for the
!                             binary files (TdPA)
!
!     04/11/2023:    V3.0.2 - Update the weights for multi-wavelength
!                             ranges (HL)
!                           - Remove the solutions after the inversion
!                             is finished (HL)
!
!     03/15/2023:    V3.0.1 - The restore file can only be INIT or a
!                             path to a file (TdPA)
!                           - Inversion now accepts an input model
!                             atmosphere for initializing (TdPA)
!                           - New argument in Init_nodes (TdPA)
!                           - Intpol_Atmo_all does not need the Flgsg
!                             argument (TdPA)
!                           - When restoring node values, the
!                             values of Bpos/theta and Bphi/azimuth
!                             are forced to be not exactly zero (TdPA)
!                           - Removed some commented lines (TdPA)
!                           - The Blos variables are in the same
!                             structure than the polar ones (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/22/2023:    V0.0.0 - Started from 05/12/2020
!                             TIC@inversion_mod.f90 revision from
!                             Hao (TdPA)
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
!    Inversion:
!      Carry out the inversion of Stokes profiles (manager)
!
!    Fit:
!      Actually fit Stokes profiles
!
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use bounds_mod
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

      !> Carry out the inversion of Stokes profiles\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!        Atomb(Atom_class): Structure with the atomic data for
      !!                           background opacities\n
      !!           Mol(Mol_class): Structure with the molecule data\n
      !!     Geom(Geometry_class): Structure with the geometry data\n
      !!    GeomI(Geometry_class): Structure with the geometry data
      !!                           for the intensity problem\n
      !!       Flgsg(Fctsg_class): Structure with factorials and
      !!                           signs\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!       fudge(fudge_class): Structure with fudge data\n
      !!     kurucz(kurucz_class): Structure with Kurucz line data\n
      !!          MPID(MPI_class): Structure with MPI data
      !!      Atmo_in(Atmo_class): Structure with atmospheric data
      !!                           read from model atmosphere\n
      !!  Bfield_in(Bfield_class): Structure with the magnetic field
      !!                           data read from the input\n
      !!       Input(Input_class): Structure with settings data\n
      !! Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                           data\n
      !!   Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!      Sol(Solution_class): Class with the data of the RT
      !!                           solution\n
      !!           imask(integer): Indicate if this pixel is masked
      !!                           in the restart
      subroutine Inversion(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec, &
                           fudge,kurucz,MPID,Atmo_in, &
                           Bfield_in,Input,Inf_Stokes, &
                           Inf_Nodes,Sol,imask)


      ! IO
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Fctsg_class):: Flgsg
      type(Geometry_class):: GeomI, Geom
      type(Frequency_class):: Frec
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(MPI_class):: MPID
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

      integer:: i, j

      double precision, dimension(Input%nvar):: TMP_Weight


      ! Print current pixel
      if (pid.eq.0) then
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
      end if

      ! Set up from common variable
      Atmo%nz = nz

      ! Store current regularization weights
      TMP_Weight = Inf_Nodes%Regul_weight

      ! Get LOS
      Inf_Nodes%mu = Inf_Stokes%mu
      Inf_Nodes%azimuth = Inf_Stokes%azimuth
      GeomI%L_mu(1) = Inf_Stokes%mu
      GeomI%L_theta(1) = acos(GeomI%L_mu(1))
      GeomI%L_phi(1) = Inf_Stokes%azimuth
      if (Input%Type_inversion.ne.0) then
        Geom%L_mu(1) = Inf_Stokes%mu
        Geom%L_theta(1) = acos(Geom%L_mu(1))
        Geom%L_phi(1) = Inf_Stokes%azimuth
      end if

      ! Initialize number of data points
      Inf_Stokes%Num_freedomI = 0
      Inf_Stokes%Num_freedom = 0

      ! Wavelengths
      do j=1,Inf_Stokes%Num_wavelength

        ! If non-zero intensity weight, add to count
        if (Inf_Stokes%Weight(0,j).gt.0d0) then

          ! Total
          Inf_Stokes%Num_freedom = 1 + Inf_Stokes%Num_freedom
          Inf_Stokes%Num_freedomI = 1 + Inf_Stokes%Num_freedomI

        end if ! Non-zero weight

        ! For each Stokes parameter
        do i=1,3

          ! If non-zero weight, add to count
          if (Inf_Stokes%Weight(i,j).gt.0d0) &
            Inf_Stokes%Num_freedom = 1 + Inf_Stokes%Num_freedom

        end do ! Stokes parameters
      end do ! Wavelengths

      ! If asymmetry nodes
      if (Inf_Nodes%Num_Asymmetry.gt.0) Input%nasym = 1

      ! If not restoring and stratification from inputs
      if (trim(Input%Inv_init).eq.'INIT'.and. &
          Input%Atmo_Input.gt.0) then

        ! Generate stratification and interpolate
        call Atmo_Stratify(Atmo_in,Atmo, &
                           Bfield_in,Bfield, &
                           Input%Atmo_strat_done)

      ! Restoring or copying stratification
      else

        ! Copy input into actual model atmosphere
        call cAtmo(Atmo_in,Atmo)
        Bfield = Bfield_in

      end if

      ! If Asymmetry not allocated, do it
      if (.not.allocated(Atmo%JKQin)) then
        allocate(Atmo%JKQin(8*nz))
        Atmo%JKQin = 0d0
      end if

      ! If gas pressure from model atmosphere
      if (.not.Inf_Nodes%Pg_auto) then

        ! Boundary from model
        Inf_Nodes%Pg_bound = Atmo%Pg(1)

      end if

      ! If LOS field
      if (Inf_Nodes%Btype.eq.1) then

        ! Allocate magnetic field arrays
        allocate(Bfield%Bpos(Atmo%nZ),Bfield%Azimuth(Atmo%nZ))
        allocate(Bfield%Blos(Atmo%nZ))

        ! If magnetic field
        if (maxval(Bfield%Bstrength).gt.TINYB) then

          ! Convert to LOS
          call B2Blos(Atmo%nZ, &
                      GeomI%L_mu(1),  &
                      GeomI%L_phi(1), &
                      Bfield%Bstrength, &
                      Bfield%Btheta, &
                      Bfield%Bphi, &
                      Bfield%Blos, &
                      Bfield%Bpos, &
                      Bfield%Azimuth)

        ! No field
        else

          ! Just zero
          Bfield%Blos = 0d0
          Bfield%Bpos = 0d0
          Bfield%Azimuth = 0d0

        end if
      end if

      ! If LOS velocity
      if (Inf_Nodes%vtype.eq.1) then

        ! Allocate velocity field arrays
        allocate(Atmo%vlos(Atmo%nZ),Atmo%vpos(Atmo%nZ))
        allocate(Atmo%vphi(Atmo%nZ))

        ! If velocity
        if (dyn) then

          ! Convert to LOS
          call v2vlos(Atmo%nZ, GeomI%L_mu(1), GeomI%L_phi(1), &
                      Atmo%vx, Atmo%vy, Atmo%vz, &
                      Atmo%vlos, Atmo%vpos, Atmo%vphi)

        ! No field
        else

          ! Just zero
          Atmo%vlos = 0d0
          Atmo%vpos = 0d0
          Atmo%vphi = 0d0

        end if
      end if


      ! Initialize node positiones
      call Init_Nodes(Atmo,Input,Inf_Nodes)

      ! Initialize node values
      call Initialize_Nodes(Atmo,Bfield,Inf_Nodes)

      ! Set regularization constants
      if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg).and. &
          Inf_Nodes%hydroeq) &
        Inf_Nodes%Const(Inf_Nodes%index_Pg) = &
                             Inf_Nodes%Node(Inf_Nodes%index_Pg)%Var(1)
      Inf_Nodes%Const(Inf_Nodes%index_f) = Input%f_diff

      ! If retoring
      if (trim(Input%Inv_init).ne.'INIT'.and.imask.eq.0) then

        ! Magnetic field inclination index
        i = Inf_Nodes%index_Bt

        ! If inverting magnetic field inclination or BPOS
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If maximum value too small (1d-3)
          if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-3) then

            ! Verbose
            if (Inf_Nodes%Btype.eq.0) then
              write(umsg, '(A)') &
                ' - Magnetic field inclination re-initialized: '
              call verboseI(3)
            else
              write(umsg, '(A)') &
                ' - Transversal magnetic field re-initialized: '
              call verboseI(3)
            end if

            ! Get from Input structure
            call re_set_nodes(Inf_Nodes,i,Input%ini_Bpos)

          end if ! If too small inclination or Bpos
        end if ! If inversing Bpos

        ! Magnetic field azimuth index
        i = Inf_Nodes%index_Bp

        ! If inverting magnetic field azimuth
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If maximum value too small (1d-3)
          if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-3) then

            ! Verbose
            write(umsg, '(A)') &
              ' - Magnetic field azimuth re-initialized: '
            call verboseI(3)

            ! Get from Input structure
            call re_set_nodes(Inf_Nodes,i,Input%ini_Bazi)

          end if ! If too small inclination or Bpos
        end if ! If inversing Bpos

        ! Velocity x/pos
        i = Inf_Nodes%index_vx

        ! If inverting vx or vpos
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If maximum value too small (1d-3)
          if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-9/c) then

            ! Verbose
            if (Inf_Nodes%vtype.eq.0) then
              write(umsg, '(A)') &
                ' - Velocity X component re-initialized: '
              call verboseI(3)
            else
              write(umsg, '(A)') &
                ' - Transversal velocity re-initialized: '
              call verboseI(3)
            end if

            ! Get from Input structure
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

              ! Verbose
              write(umsg, '(A)') &
                  ' - Velocity azimuth re-initialized: '
                call verboseI(3)

              ! Get from Input structure
              call re_set_nodes(Inf_Nodes,i,Input%ini_vazi)

            end if ! If too small vy

          ! If vy
          else

            ! If maximum value too small (1d-3)
            if (maxval(abs(Inf_Nodes%Node(i)%Var)).lt.1d-9/c) then

              write(umsg, '(A)') &
                  ' - Velocity Y component re-initialized: '
                  call verboseI(3)

              ! Get from Input structure
              call re_set_nodes(Inf_Nodes,i,Input%ini_vazi*1d-6/c)

            end if ! If too small v azimuth
          end if ! vy or azimuth
        end if ! inverting vy or azimuth
      end if ! Type of restore file and mask

      ! If regularizing, initialize regularization
      if (Inf_Nodes%Regul_Flag) then
        call Init_Regul(Inf_Nodes)
        if (laborted) goto 1000
      end if

      ! Put node values within boundaries (magnetic field azimuth)
      i = Inf_Nodes%index_Bp
      call FoldBounds(Inf_Nodes%Node(i),Inf_Nodes%Num_Nodes(i))

      ! Put node values within boundaries (velocity azimuth)
      i = Inf_Nodes%index_vy
      if (Inf_Nodes%vtype.eq.1) &
        call FoldBounds(Inf_Nodes%Node(i),Inf_Nodes%Num_Nodes(i))


      ! Master verbose
      if (pid.eq.0) then
        umsg = ' - Initialized model atmosphere/nodes'
        call verboseI(3)
      end if

      !
      ! Select depending on the type of inversion
      select case(Input%Type_Inversion)

        ! Only thermal
        case(0)

          ! If there are thermal nodes
          if (Inf_Nodes%Num_Thermal.gt.0) then

            ! Master verbose
            if (pid.eq.0) then
              umsg = " - Fitting the thermal and dynamical parameters"
              call verboseI(3)
            end if

            ! Force inputs
            Inf_Nodes%Nodes_type = 0
            Input%force = 'I'

            ! Interpolate into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                   Atomb, Mol, Input, fudge)

            ! Allocate mangetic arrays
            allocate(Bfield0%Bstrength(Atmo%nZ))
            allocate(Bfield0%Btheta(Atmo%nZ))
            allocate(Bfield0%Bphi(Atmo%nZ))
            allocate(Bfield0%Blos(Atmo%nZ))
            allocate(Bfield0%Bpos(Atmo%nZ))
            allocate(Bfield0%Azimuth(Atmo%nZ))

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
            deallocate(Bfield0%Bstrength,Bfield0%Btheta,Bfield0%Bphi)
            deallocate(Bfield0%Blos,Bfield0%Bpos,Bfield0%Azimuth)

          end if ! If there are thermal nodes

        ! Only magnetic
        case(1)

          ! If there are magnetic nodes
          if (Inf_Nodes%Num_Mag.gt.0) then

            ! Master verbose
            if (pid.eq.0) then
              umsg = " - Fitting the of magnetic parameters"
              call verboseI(3)
            end if

            ! Force inputs
            Inf_Nodes%Nodes_type = 1
            Input%force = 'N'

            ! Interpolate into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                   Atomb, Mol, Input, fudge)

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

          end if ! Magnetic nodes

        ! Fit all, together
        case(2)

          ! Master verbose
          if (pid.eq.0) then
            umsg = " - Fitting the thermal, dynamical, "// &
                   "and magnetic parameters together"
            call verboseI(3)
          end if

          ! Force inputs
          Inf_Nodes%Nodes_type = 2
          Input%force = 'N'

          ! Interpolate into the atmosphere
          if (imask.eq.0) &
            call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                 Atomb, Mol, Input, fudge)

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

            ! Initialize magnetic field
            Bfield0%Bstrength = 0d0
            Bfield0%Btheta = 0d0
            Bfield0%Bphi = 0d0
            Bfield0%Blos = 0d0
            Bfield0%Bpos = 0d0
            Bfield0%Azimuth = 0d0

            ! Master verbose
            if(pid.eq.0) then
              umsg = " - Fitting the thermal and dynamical parameters"
              call verboseI(3)
            end if

            ! Force inputs
            Inf_Nodes%Nodes_type = 0
            Input%force = 'I'

            ! Interpolate into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                   Atomb, Mol, Input, fudge)

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield0,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.False.)

            ! Deallocate magnetic field quantities
            deallocate(Bfield0%Bstrength,Bfield0%Btheta,Bfield0%Bphi)
            deallocate(Bfield0%Blos,Bfield0%Bpos,Bfield0%Azimuth)

          end if ! Thermal nodes

          ! Verbose merit function if in output
          umsg = ' * '
          call verboseI(0)
          if (vlevel.eq.0) call verboseI(3)

          ! Verbose merit function
          write(umsg,'(A,i4,2(3x,A,es15.4))')  &
              ' * Thermal cycle ended, starting full cycle'

          if (gpid.eq.0) then
            call verboseI(0)
            call verboseI(4)
          else
            call verboseI(0)
            if (vlevel.eq.0) call verboseI(3)
          end if

          ! There are magnetic nodes
          if (Inf_Nodes%Num_Mag.gt.0.and..not.laborted) then

            ! Master verbose
            if(pid.eq.0) then
              umsg = " - Fitting the thermal, dynamical, "// &
                     "and magnetic parameters together"
              call verboseI(1)
            end if

            ! Force inputs
            Inf_Nodes%Nodes_type = 2
            Input%force = 'N'

            ! Interpolate into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                   Atomb, Mol, Input, fudge)

            ! Fit profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

          end if

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

            ! Initialize magnetic field
            Bfield0%Bstrength = 0d0
            Bfield0%Btheta = 0d0
            Bfield0%Bphi = 0d0
            Bfield0%Blos = 0d0
            Bfield0%Bpos = 0d0
            Bfield0%Azimuth = 0d0

            ! Master verbose
            if(pid.eq.0) then
              umsg = " - Fitting the thermal and dynamical parameters"
              call verboseI(3)
            end if

            ! Force inputs
            Inf_Nodes%Nodes_type = 0
            Input%force = 'I'

            ! Interpolate into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                   Atomb, Mol, Input, fudge)

            ! Fit the profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield0,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.False.)

            ! Deallocate magnetic field quantities
            deallocate(Bfield0%Bstrength,Bfield0%Btheta,Bfield0%Bphi)
            deallocate(Bfield0%Blos,Bfield0%Bpos,Bfield0%Azimuth)

          end if ! Thermal nodes

          ! Verbose merit function if in output
          umsg = ' * '
          call verboseI(0)
          if (vlevel.eq.0) call verboseI(3)

          ! Verbose merit function
          write(umsg,'(A,i4,2(3x,A,es15.4))')  &
              ' * Thermal cycle ended, starting only magnetic cycle'

          if (gpid.eq.0) then
            call verboseI(0)
            call verboseI(4)
          else
            call verboseI(0)
            if (vlevel.eq.0) call verboseI(3)
          end if

          ! There are magnetic nodes
          if (Inf_Nodes%Num_Mag.gt.0.and..not.laborted) then

            ! Master verbose
            if(pid.eq.0) then
              umsg = " - Fitting the thermal, dynamical, "// &
                     "and magnetic parameters together"
              call verboseI(1)
            end if

            ! Force inputs
            Inf_Nodes%Nodes_type = 1
            Input%force = 'N'

            ! Interpolate into the atmosphere
            if (imask.eq.0) &
              call Intpol_Atmo_all(Inf_Nodes, Atmo, Bfield, Atom, &
                                   Atomb, Mol, Input, fudge)

            ! Fit profiles
            call Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,.True.)

          end if

      end select ! Type of inversion

      ! Restore regularization weights
1000  Inf_Nodes%Regul_weight = TMP_Weight

      ! Purge atmospheric model
      call free_Atmo(Atmo,.True.)
      call free_B(Bfield)

      ! If master
      if (pid.eq.0) then

        ! Message
        if (laborted) then
          umsg = " # Pixel inversion failed"
        else
          umsg = " - Pixel inversion finished"
        end if

        ! Master verbose
        if (gpid.eq.0) then
          call verboseI(0)
          call verboseI(4)
        else
          call verboseI(4)
        end if
      end if

      return

      end subroutine Inversion

!#####################################################################
!#####################################################################
!#####################################################################

      !> Fit the Stokes profiles\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!         Input(Input_class): Structure with settings data\n
      !!   Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                             data\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!             imask(integer): Indicate if this pixel is masked
      !!                             in the restart\n
      !!            saving(logical): If the result is to be stored
      subroutine Fit(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                     kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                     Inf_Nodes,Sol,imask,saving)


      ! IO
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Fctsg_class):: Flgsg
      type(Geometry_class):: GeomI, Geom
      type(Frequency_class):: Frec
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(MPI_class):: MPID
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


      ! If thermal inversion
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Set thermal indexes
        LM_Stru%Num = Inf_Nodes%Num_Thermal + Inf_Nodes%Num_glob

        ! If diffuse light
        if (Sol%Diff_flag) then
          Inf_Nodes%Indx_b = Inf_Nodes%index_f
        else
          Inf_Nodes%Indx_b = Inf_Nodes%index_T
        end if

        Inf_Nodes%Indx_e = Inf_Nodes%index_Pg

      ! If magnetic inversion
      else if (Inf_Nodes%Nodes_Type.eq.1) then

        ! Set Magnetic indexes
        LM_Stru%Num = Inf_Nodes%Num_Mag + Inf_Nodes%Num_Asymmetry + &
                      Inf_Nodes%Num_glob

        Inf_Nodes%Indx_b = Inf_Nodes%index_B

        ! If diffuse light
        if (Sol%Diff_flag) then
          Inf_Nodes%Indx_e = Inf_Nodes%index_f
        else
          Inf_Nodes%Indx_e = Inf_Nodes%index_Bp
        end if

      ! If full inversion
      else if (Inf_Nodes%Nodes_Type.eq.2) then

        ! Set full indexes
        LM_Stru%Num = Inf_Nodes%Num_Fit
        Inf_Nodes%Indx_b = 1
        Inf_Nodes%Indx_e = Input%nvar

      end if ! Inversion type

      ! Call LM fit routine
      call LMFIT(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge,kurucz, &
                 MPID,Atmo,Bfield,Input,Inf_Stokes,Inf_Nodes,Sol, &
                 LM_Stru,imask,saving)

      return

      end subroutine Fit

!#####################################################################
!#####################################################################
!#####################################################################

      end module inversion_mod
