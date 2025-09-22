      !> Inversion main
      module tic_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!  Start:
!     16/02/2023
!  Last version:
!     28/08/2025 V4.0.7
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/08/2025:    V4.0.7 - Initialize the JKQ tensors if the
!                             result being restored had them even if
!                             not inverting them (TdPA)
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
!  TIC
!    Perform the inversion of the given spectropolarimetric data
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use background_mod
      use commons_mod
      use free_mod
      use gauss_mod
      use initinv_mod
      use inversion_mod
      use io_mod
      use iotic_mod
      use omegabuild_mod
      use parameters_mod , only: c , TINYA , TINYSP , TINYB , PI
      use psf_mod
      use ratmo_mod
      use ratom_mod
      use rbarklem_mod
      use rpfa_mod
      use types_mod
      use w2freq_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Perform the inversion of the given spectropolarimetric data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       MPID(MPI_class): Structure with MPI data
      subroutine TIC(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), allocatable, intent(inout):: Mol
      type(Input_class), intent(inout):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(inout):: fudge
      type(MPI_class), intent(inout):: MPID

      ! Local

      type(Atmo_class):: Atmo_in
      type(Bfield_class):: Bfield_in
      type(Frequency_class):: Frec
      type(Geometry_class):: GeomI,Geom
      type(kurucz_class):: kurucz
      type(Stokes_class):: Inf_Stokes
      type(Nodes_class):: Inf_Nodes
      type(Solution_class):: Sol

      logical:: aborting,check,lcache,double,restoring
      logical:: update_tlim,update_vlim,update_blim
      logical:: double_jkq,receiving,atmojkq,lexcl,warning
      logical, dimension(:,:), allocatable:: cache

      integer:: unitD,unitC,unitA,unitJ,unitM
      integer:: ix0,iy0,ix1,iy1,ix2,iy2,ix,iy,inod,ip,iproc,imask
      integer:: ia,ii,jj,NLOS,NLOSr,aindex,type_atmo_size
      integer:: s_data_buffer,s_atmo_buffer
      integer:: s_jkq_buffer,s_transfer_buffer
      integer, dimension(3):: dims,out_dims,dims_atmo
      integer, dimension(4):: int_buff,finfo
      integer, dimension(:), allocatable:: cpu_free

      double precision:: maxB,DwTa,lMRAMc

      ! Pointers

      double precision, dimension(:), pointer:: p_transfer_buffer


      ! Memory count
      MRAMc = MRAMc + 1d-6*(sizeof(Atmo_in) + &
                            sizeof(Bfield_in) + &
                            sizeof(Frec) + &
                            sizeof(GeomI) + &
                            sizeof(Geom) + &
                            sizeof(kurucz) + &
                            sizeof(Inf_Stokes) + &
                            sizeof(Inf_Nodes) + &
                            sizeof(Sol))

      ! Initialize miscellaneous memory warning
      Sol%warning = .True.
      warning = .True.
      lMRAMc = -2.5d0


      !
      ! From now on we will distinguish between the inversion
      ! verbosity and the synthesis verbosity
      !

      ! Copy the verbosity file name to the inversion
      verbosefv = verbosef

      ! Set up the inversion structures
      call set_up_inversion(Input,Inf_Nodes,Inf_Stokes,Sol)


      ! Get index atoms
      call set_atom_indexes(Atom,Input,.True., &
                            Input%Type_Inversion.ne.0, &
                            .True., &
                            Input%Type_Inversion.ne.0)


#ifdef FITSSUP
      ! Determine if fits file
      call name_check(Input%Filename_ob, Input%Fits_Index)

      ! If its a fits
      Input%FITSFILE = Input%Fits_Index.gt.0
#else
      ! No fits allowed
      Input%FITSFILE = .False.
#endif


      ! Split in groups of tasks
      call setmpi15D(MPID,Input)
      if (gpid.eq.0) then
        umsg = ' - Tasks distributed'
        call verbosev
      end if

      ! Initialize units
      unitA = -1
      unitD = -1
      unitC = -1
      unitJ = -1
      unitM = -1

      ! Slaves, add identifiers for verbosity if 1.5D
      if (gpid.gt.0.and.MPID%mpi15d) &
        write(verbosefv,'(A,"_",i5.5)') trim(verbosefv),gpid


      !
      ! Get inversion data from file
      !

      ! Master
      if (gpid.eq.0) then

        ! Open files (ia and nz are a dummy variable here)
        call open_data_and_cache(Input,unitD,unitC, &
                                 aborting,dims,finfo, &
                                 cache,lcache)

        ! The master does not need the background atoms or
        ! molecules
        if (MPID%mpi15d) then
          call free_atom_full(Atomb)
          call free_mol_full(Mol)
        end if

        ! Check if could read
        laborted = aborting

        ! Open mask
        call open_mask(Input,unitM,aborting,dims)

      end if ! G. master

      ! Check if aborting
      call gcontrol

      ! Share dims
      ! dims(1) = nx (slow axis)
      ! dims(2) = ny (fast axis)
      ! dims(3) = number_wavelengths
      call MPI_BCAST(dims(1),3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ! Share file info
      ! finfo(1) = 0 If only intensity, 1 if full Stokes
      ! finfo(2) = 0 If only one LOS, 1 if a LOS per position
      ! finfo(3) = 0 no sigma
      !            1 constant sigma common for every pix
      !            2 wavelength dependent sigma common for every pix
      !            3 constant sigma at each pix
      !            4 wavelength dependent sigma at each pix
      ! finfo(4) = 0 no diffuse light
      !            1 single diffuse light only intensity
      !            2 single diffuse light polarized
      !            3 intensity diffuse light at each pixel
      !            4 polarized diffuse light at each pixel
      call MPI_BCAST(finfo(1),4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


      !
      ! Get lambda
      !

      ! Save size in the right structures
      Sol%Num_Wavelength = dims(3)
      Inf_Stokes%Num_Wavelength = dims(3)

      ! Allocate Stokes and frequency
      allocate(Sol%omega_input(Inf_Stokes%Num_Wavelength))
      allocate(Sol%Stokes_out(0:3,Inf_Stokes%Num_Wavelength))
      allocate(Inf_Stokes%Stokes_Ob(0:3,Inf_Stokes%Num_Wavelength))
      MRAMc = MRAMc + 1d-6*sizeof(Sol%omega_input)
      MRAMc = MRAMc + 1d-6*sizeof(Sol%Stokes_out)
      MRAMc = MRAMc + 1d-6*sizeof(Inf_Stokes%Stokes_ob)

      ! Master
      if (gpid.eq.0) then

        ! Get wavelength axis
        call get_data_wavelength(unitD,aborting,Sol,Input%FITSFILE)

        ! Check if could read
        laborted = aborting

      end if ! Master

      ! Check if aborting
      call gcontrol

      ! Share lambda
      call MPI_BCAST(Sol%omega_input(1),Sol%Num_Wavelength, &
                     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      ! Convert wavelengths
      call Range_Check(Sol%omega_input,Inf_Stokes)

      ! Set-up weights if not automatic
      if (.not.Inf_Stokes%auto_weight) &
        call inversion_weights(Input,Inf_Stokes,Sol%omega_input)

      ! Check if aborting
      call gcontrol

      !
      ! Limit the outputs for Stokes parameters
      call prepare_lambda_limits(Input,Inf_Stokes,Sol%omega_input)


      !
      ! Calculate expected buffer size
      !

      ! Basic size
      s_data_buffer = dims(3)

      ! If LOS included
      if (finfo(2).eq.1) s_data_buffer = s_data_buffer + 2

      ! Only intensity
      if (finfo(1).eq.0) then

        ! If sigma constant
        if (finfo(3).eq.3) then

          s_data_buffer = s_data_buffer + 1

        ! If sigma wavelength dependent
        else if (finfo(3).eq.4) then

          s_data_buffer = s_data_buffer + dims(3)

        end if

        ! If polarized diffuse light
        if (finfo(4).eq.2.or.finfo(4).eq.4) then

          ! Error
          umsg = 'The data file is only intensity, but has '// &
                 'polarized diffuse light, the file is not '// &
                 'valid'
          urou = 'TIC'
          call gabortedv

        end if

        ! If there is pixel dependent diffuse light
        if (finfo(4).eq.1) s_data_buffer = s_data_buffer + dims(3)

      ! Polarization
      else

        ! Add Stokes
        s_data_buffer = s_data_buffer + 3*dims(3)

        ! If sigma constant
        if (finfo(3).eq.3) then

          ! Add four numbers
          s_data_buffer = s_data_buffer + 4

        ! If sigma wavelength dependent
        else if (finfo(3).eq.4) then

          ! Add four profiles
          s_data_buffer = s_data_buffer + dims(3)*4

        end if ! Constant or profile sigma

        ! If intensity diffuse light
        if (finfo(4).eq.3) then

          ! Add one profile
          s_data_buffer = s_data_buffer + dims(3)

        ! If polarized diffuse light
        else if (finfo(4).eq.4) then

          ! Add four profiles
          s_data_buffer = s_data_buffer + dims(3)*4

        end if ! Only intensity or polarized diffuse light
      end if ! Intensity/Polarization


      ! Get constant LOS now
      if (finfo(2).eq.0) then

        ! Master
        if (gpid.eq.0) then

          call get_data_los(unitD,aborting,Inf_Stokes,Input%FITSFILE)

        endif ! Master

        ! Check if aborting
        call gcontrol

        ! Share LOS
        call MPI_BCAST(Inf_Stokes%mu,1,MPI_DOUBLE_PRECISION,0, &
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Inf_Stokes%azimuth,1,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)

      end if ! If constant LOS


      ! Get "constant" sigma now
      if (finfo(3).eq.1.or.finfo(3).eq.2) then

        ! Allocate input sigma
        allocate(Inf_Stokes%Sigma_in(0:3,Inf_Stokes%Num_Wavelength))
        MRAMc = MRAMc + 1d-6*sizeof(Inf_Stokes%Sigma_in)

        ! Flag
        Inf_Stokes%Sigma_ct = .True.

        ! Master
        if (gpid.eq.0) then

          ! Get sigma
          call get_data_sigma(unitD,aborting,finfo, &
                              Inf_Stokes,Input%FITSFILE)

          ! If to enhance Sigma
          if (allocated(Input%Sigma_factor)) then

            ! Enhance the input directly
            call enhance_sigma(Input%Sigma_factor, &
                               Inf_Stokes%Sigma_in, &
                               Sol%omega_input)

            ! Free space
            MRAMc = MRAMc - 1d-6*sizeof(Input%Sigma_factor)
            deallocate(Input%Sigma_factor)

          end if ! Enhancing Sigma

        endif

        ! Check if aborting
        call gcontrol

        ! Share sigma
        call MPI_BCAST(Inf_Stokes%Sigma_in(0,1), &
                       Inf_Stokes%Num_Wavelength*4, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        ! Master manager
        if (gpid.eq.0.and.MPID%mpi15d) then

          ! Forget
          MRAMc = MRAMc - 1d-6*sizeof(Inf_Stokes%Sigma_in)
          deallocate(Inf_Stokes%Sigma_in)

        end if ! Master manager

      ! Not constant
      else

        ! Flag
        Inf_Stokes%Sigma_ct = .False.

      end if ! If constant sigma


      ! Get "constant" diffuse light now (assuming FITS are going
      ! to read it pixel by pixel)
      if (finfo(4).eq.1.or.finfo(4).eq.2) then

        ! Allocate input sigma
        allocate(Inf_Stokes%Diff_in(0:3,Inf_Stokes%Num_Wavelength))
        MRAMc = MRAMc + 1d-6*sizeof(Inf_Stokes%Diff_in)

        ! Flag
        Inf_Stokes%Diff_ct = .True.

        ! Master
        if (gpid.eq.0) then

          ! Read diffuse light
          call get_data_diff(unitD,aborting,finfo,Inf_Stokes, &
                             Input%FITSFILE)

        endif ! Master

        ! Check if aborting
        call gcontrol

        ! Share LOS
        call MPI_BCAST(Inf_Stokes%Diff_in(0,1), &
                       Inf_Stokes%Num_Wavelength*4, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        ! Master manager
        if (gpid.eq.0.and.MPID%mpi15d) then

          ! Forget
          MRAMc = MRAMc - 1d-6*sizeof(Inf_Stokes%Diff_in)
          deallocate(Inf_Stokes%Diff_in)

        end if ! Master manager

      ! Not constant
      else

        ! Flag no constant diffuse light
        Inf_Stokes%Diff_ct = .False.

      end if ! Constant diffuse light


      ! Slaves allocate sigma if needed
      if ((gpid.gt.0.or..not.MPID%mpi15d).and.finfo(3).gt.0) then
        allocate(Inf_Stokes%Sigma_W(0:3,Inf_Stokes%Num_wavelength))
        MRAMc = MRAMc + 1d-6*sizeof(Inf_Stokes%Sigma_W)
      end if


      ! Slaves allocate diffuse light if needed
      if (gpid.gt.0.or..not.MPID%mpi15d) then

        ! There is diffuse light in the file
        if (finfo(4).gt.0) then

          ! Allocate
          allocate(Sol%Stokes_diff(0:3,Inf_Stokes%Num_wavelength))
          MRAMc = MRAMc + 1d-6*sizeof(Sol%Stokes_diff)

          ! If there is diffuse light
          if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_f).or. &
              Input%f_diff.gt.0d0) then

            ! There is diffuse light
            Sol%Diff_flag = .True.

          ! No diffuse light
          else

            ! There is no diffuse light afterall
            Sol%Diff_flag = .False.

          end if ! Diffuse light?

        ! De-flag diffuse light
        else

          ! Do not apply
          Sol%Diff_flag = .False.

          ! If inverting diffuse light
          if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_f)) then

            ! Issue error
            umsg = 'Cannot invert diffuse light if not '// &
                   'provided in the data file'
            urou = 'TIC'
            call aborted

          end if ! Inverting diffuse light
        end if ! If there is diffuse light
      end if ! Slaves or no-manager master

      ! Control
      call gcontrol


      ! Read partition function data and abundances
      if (gpid.gt.0.or..not.MPID%mpi15d) &
        call rParfunAbund(Input,Atmo_in)


      ! Initialize J from atmosphere
      atmojkq = .False.


      !
      ! Deal with input model
      !

      ! Get input Atmosphere from restoring file
      if (trim(Input%Inv_init).ne.'INIT') then

        ! Flag restoring
        restoring = .True.

        ! Check label
        if (gpid.eq.0) &
          call get_atmo_type(Input%Inv_init, &
                             Input%atmoin_type)

        ! Share type
        call MPI_BCAST(Input%atmoin_type,1,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)

        ! Overwrite atmo
        Input%atmo = Input%Inv_init

        ! If not output, abort
        if (Input%atmoin_type.ne.2) then

          ! Error
          umsg = 'The restore file '//trim(Input%Inv_init)// &
                 ' is not a valid output inversion file'
          urou = 'TIC'
          call gabortedv

        end if ! Error

        ! No automatic gas pressure
        Inf_Nodes%Pg_auto = .False.

      ! Initialize
      else

        ! Flag no restoring
        restoring = .False.

        ! No input Atmosphere
        if (trim(Input%atmo).eq.'NONE') then

          ! Call hard-coded model
          call gAtmo(Atmo_in,Input%Init_Thermal)

          ! If no input field, generate one
          if (trim(Input%bfield).eq.'NONE'.and. &
              .not.Input%bfieldn) then

            ! Hard-code input
            Input%bfieldn = .True.

            ! If thermal
            if (Input%Type_Inversion.eq.0) then

              ! Zero field
              Input%bfieldv = (/ 0d0, 0d0, 0d0 /)

            ! Magnetic
            else

              ! Depending on the inverted quantities
              ii = Inf_Nodes%index_Bt
              jj = Inf_Nodes%index_Bp

              ! If inverting Bt or Bp
              if (Inf_Nodes%Nodes_Flags(ii).and. &
                  Inf_Nodes%Nodes_Flags(jj)) then

                ! Add an arbitrary direction field
                Input%bfieldv = (/ 1d0, PI*0.25d0, 5.7d0 /)

              ! If inverting Bt and not Bp
              else if (Inf_Nodes%Nodes_Flags(ii)) then

                ! Add an inclined field at azimuth 0
                Input%bfieldv = (/ 1d0, PI*0.25d0, 0.0d0 /)

              ! If inverting Bp and not Bt
              else if (Inf_Nodes%Nodes_Flags(jj)) then

                ! Add a vertical magnetic field
                Input%bfieldv = (/ 1d0,       0d0, 5.7d0 /)

              ! Not inverting either
              else

                ! Add a vertical magnetic field
                Input%bfieldv = (/ 1d0,       0d0,   0d0 /)

              end if ! Inverting Bt or Bp
            end if ! Thermal or magnetic inversion
          end if ! No magnetic field model and no input

          ! Set-up initial field (global master does not care)
          call rBField(Input%bfield,Input%source, &
                       Input%ID,Bfield_in,Atmo_in%nz,Input)

          ! Flag input
          Input%atmoin_type = 0

        ! Get input Atmosphere, if not restoring
        else

          ! Check label
          if (gpid.eq.0) &
            call get_atmo_type(Input%atmo,Input%atmoin_type)

          ! Share
          call MPI_BCAST(Input%atmoin_type,1,MPI_INTEGER,0, &
                         MPI_COMM_WORLD,ierr)

          ! If no type, try with 1D
          if (Input%atmoin_type.lt.0) then

            ! Read input atmosphere
            call rAtmo(Input%atmo,Input%source, &
                       Input%ID, Atmo_in, -1d0)

            ! If non-numeric input
            if (.not.Input%bfieldn) then

              ! If no input field, generate one
              if (trim(Input%bfield).eq.'NONE') then

                ! Hard-code input
                Input%bfieldn = .True.

                ! If thermal
                if (Input%Type_Inversion.eq.0) then

                  ! No magnetic field
                  Input%bfieldv = (/ 0d0, 0d0, 0d0 /)

                ! Magnetic
                else

                  ! Depending on the inverted quantities
                  ii = Inf_Nodes%index_Bt
                  jj = Inf_Nodes%index_Bp

                  ! If inverting Bt or Bp
                  if (Inf_Nodes%Nodes_Flags(ii).and. &
                      Inf_Nodes%Nodes_Flags(jj)) then

                    ! Add an arbitrary direction field
                    Input%bfieldv = (/ 1d0, PI*0.25d0, 5.7d0 /)

                  ! If inverting Bt and not Bp
                  else if (Inf_Nodes%Nodes_Flags(ii)) then

                    ! Add an inclined field at azimuth 0
                    Input%bfieldv = (/ 1d0, PI*0.25d0, 0.0d0 /)

                  ! If inverting Bp and not Bt
                  else if (Inf_Nodes%Nodes_Flags(jj)) then

                    ! Add a vertical magnetic field
                    Input%bfieldv = (/ 1d0,       0d0, 5.7d0 /)

                  ! Not inverting either
                  else

                    ! Add a vertical magnetic field
                    Input%bfieldv = (/ 1d0,       0d0,   0d0 /)

                  end if ! Inverting Bt or Bp
                end if ! Thermal or magnetic inversion
              end if ! No magnetic field model
            end if ! Non-numeric input

            ! Set-up initial field
            call rBField(Input%bfield,Input%source, &
                         Input%ID,Bfield_in,Atmo_in%nz,Input)

            ! Flag initialized
            Input%atmoin_type = 0

          end if ! No 3D

          ! send to the global master
          call MPI_BCAST(Atmo_in%nz,1,MPI_INTEGER,0, &
                         MPI_COMM_WORLD,ierr)

        end if ! Input atmosphere
      end if ! Restore


      !
      ! Now fix the real dimensions from solution box
      !

      ! Fix wildcards in input
      if (Input%sol_box(1).lt.1) Input%sol_box(1) = 1
      if (Input%sol_box(2).lt.1) Input%sol_box(2) = dims(1)
      if (Input%sol_box(2).gt.dims(1)) Input%sol_box(2) = dims(1)
      if (Input%sol_box(3).lt.1) Input%sol_box(3) = 1
      if (Input%sol_box(4).lt.1) Input%sol_box(4) = dims(2)
      if (Input%sol_box(4).gt.dims(2)) Input%sol_box(4) = dims(2)

      ! Global master
      if (gpid.eq.0) then

        ! Get output dimensions
        out_dims(1) = Input%sol_box(2) - Input%sol_box(1) + 1
        out_dims(2) = Input%sol_box(4) - Input%sol_box(3) + 1
        out_dims(3) = dims(3)

        ! Check output dimensions
        if (out_dims(1).lt.1.or.out_dims(2).lt.1) then

          ! Issue error
          write(umsg,'(A,i4,",",i4)') &
            ' # Error: X-Y zero size after applying '// &
            'the SOLUTION_BOX limits: ',out_dims(1:2)
          call verbosev
          laborted = .True.

        end if ! Final output dimensions

      ! Others
      else

        ! Overwrite dimensions
        dims(1) = Input%sol_box(2) - Input%sol_box(1) + 1
        dims(2) = Input%sol_box(4) - Input%sol_box(3) + 1

      end if ! Global master or not

      ! Control
      call gcontrol

      ! If reading from a 3D model, open the file
      if (Input%atmoin_type.gt.0) then

        ! Check known atmospheric limits
        if ((PRD.and.(Input%minT.lt.0d0.or.Input%maxT.lt.0d0)).or. &
            (Input%minT.lt.0d0.and.Input%dws.eq.'MIN').or. &
            (Input%maxT.lt.0d0.and.Input%dws.eq.'MAX').or. &
            (.not.Input%static.and.Input%maxV.lt.0d0)) then

          ! Read whole model to get limits
          call get_lims(Input,-1,aborting)

          ! Check if could read
          laborted = aborting

        end if ! Unknown atmospheric limits

        ! Master
        if (gpid.eq.0.and..not.laborted) then

          ! Open files (ii is a dummy variable here)
          call open_atm(Input,1,unitA,aborting,dims_atmo, &
                        aindex,double,ii)

          ! Check if could read
          laborted = aborting

          ! Check dimensions
          if (.not.laborted) then

            ! Check the data sizes
            if (dims_atmo(1).eq.dims(1).and. &
                dims_atmo(2).eq.dims(2)) then

              ! Model atmosphere same dimensions original data
              type_atmo_size = 1

            ! Check the output sizes
            else if (dims_atmo(1).eq.out_dims(1).and. &
                     dims_atmo(2).eq.out_dims(2)) then

              ! Model atmosphere same dimensions than output
              type_atmo_size = 2

            ! Problem with sizes
            else

              ! Issue error
              write(umsg,'(3(A,i4,",",i4))') &
                ' # Error: X-Y size in data ', &
                dims(1:2), &
                ' and output ',out_dims(1:2), &
                ' different in model atmosphere ', &
                dims_atmo(1:2)
              call verbosev
              laborted = .True.

            end if ! Right sizes
          end if ! Success getting dimensions
        end if ! Master

        ! Check if aborting
        call gcontrol

        ! Share data
        call MPI_BCAST(aindex,1,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dims_atmo(1),3,MPI_INTEGER,0, &
                       MPI_COMM_WORLD,ierr)

        ! If inversion result, get if JKQ data
        if (Input%atmoin_type.eq.2) atmojkq = aindex.gt.7


        ! Initialize pointers
        call iAtmo_p(Atmo_in)

        !
        ! Prepare buffer atmosphere
        !

        ! If 1.5DS model
        if (Input%atmoin_type.eq.1) then

          ! Full synthesis size
          s_atmo_buffer = dims_atmo(3)*24

        ! If inversion model
        else

          ! Yes JKQ
          if (atmojkq) then

            ! Full inversion size
            s_atmo_buffer = dims_atmo(3)*27 + 1

          ! No JKQ
          else

            ! Reduced inversion size
            s_atmo_buffer = dims_atmo(3)*19 + 1

          end if ! JKQ data
        end if ! Synthesis or inversion model

      ! Just 1D model
      else

        ! Flag type of atmosphere
        s_atmo_buffer = 0

      end if ! Reading from a 3D model


      !
      ! JKQ asymmetry file?
      !

      ! If restoring file with JKQ
      if (restoring.and.atmojkq) then

        ! Buffer size
        s_jkq_buffer = 0

        ! Flag restore
        double_jkq = .False.
        unitJ = -1

      ! JKQ asymmetry file?
      else if (Input%nasym.gt.0) then

        ! Master
        if (gpid.eq.0) then

          ! Open asymmetry file
          call open_asymm(Input,unitJ,aborting,dims_atmo)

          ! Check if could read
          laborted = aborting

        end if ! Master

        ! Check if aborting
        call gcontrol

        ! Buffer size
        s_jkq_buffer = dims_atmo(3)*8

        ! Flag restore
        double_jkq = .True.

      ! No asymmetry file
      else

        ! No buffer size
        s_jkq_buffer = 0
        double_jkq = .False.
        unitJ = -1

      end if ! Asymmetry file


      ! Allocate transfer buffer
      s_transfer_buffer = s_data_buffer + s_atmo_buffer + s_jkq_buffer
      allocate(p_transfer_buffer(s_transfer_buffer))
      MRAMc = MRAMc + 1d-6*sizeof(p_transfer_buffer)


      !
      ! Set up initial atmospheric model
      !

      ! If restoring
      if (restoring) then

        ! Get size from input model
        Atmo_in%nz = dims_atmo(3)
        nz = Atmo_in%nz

        ! Initialize
        Atmo_in%tfreq = Input%omega_ref
        Atmo_in%scal = 'T'
        Atmo_in%logg = 4.44d0

      ! If the input was to initialize or no file
      else if (trim(Input%Inv_init).eq.'INIT') then

        ! If input model atmosphere not 1D
        if (Input%atmoin_type.gt.0) then

          ! Reference
          Atmo_in%tfreq = Input%omega_ref
          Atmo_in%scal = Input%atm_scale
          Atmo_in%logg = 4.44d0

        end if ! Not a 1D input model

        ! If more than 0 nodes in the atmosphere
        if (Input%Atmo_Input.gt.0) then

          ! Set up height dimension
          nz = Input%Atmo_Input

          ! Build stratification now
          call gAtmo_strat(Input%Tau_Range(1), &
                           Input%Tau_Range(2), &
                           Input%Atmo_strat, &
                           Input%Atmo_input, &
                           Input%Atmo_strat_done)

        ! If 0 nodes, then read the input atmosphere
        else

          ! If reading a 1D model
          if (Input%atmoin_type.eq.0) then

            ! Set up height dimension
            nz = Atmo_in%nz

          ! Copy the stratification
          else

            ! Set up height dimension
            nz = dims_atmo(3)

          end if ! 1D or 3D initialization
        end if ! Number of nodes
      end if ! File not specified

      ! If nodes in JKQin
      if ((Inf_Nodes%Num_Asymmetry.gt.0.or.s_jkq_buffer.gt.0.or. &
           atmojkq).and.(gpid.gt.0.or..not.MPID%mpi15d)) then

        ! Allocate and initialize
        allocate(Atmo_in%JKQin(8*nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo_in%JKQin)
        Atmo_in%JKQin = 0d0

      end if ! JKQ nodes


      ! Read barklem data
      if (.not.MPID%mpi15d.or.gpid.gt.0) &
        call rBarklem(Input,Atom,Atomb)

      ! Control
      call gcontrol


      ! Get limits for T, v, and B (global master does not care)
      if (.not.MPID%mpi15d.or.gpid.gt.0) &
        call set_up_limits(Input,Inf_Nodes,Atmo_in,Bfield_in, &
                           maxB,update_Tlim,update_vlim,update_Blim)


      ! Decide if need to keep Stokes
      if (.not.MPID%mpi15d.or.(pid.eq.0.and.gpid.ne.0)) then
        KSTK = KSTK.or.(PRD.and.(dyn.or..not.AV)).or. &
               (dyn.and.Input%Type_inversion.gt.0)
      else
        KSTK = (PRD.and.(dyn.or..not.AV)).or. &
               (dyn.and.Input%Type_inversion.gt.0)
      end if


      !
      ! Set angular quadrature
      !

      ! If only thermal inversion
      if (Input%Type_inversion.eq.0) then

        call gauss(Input,GeomI,Geom,1,.False.,.True.,Flgsg)

      ! Not only thermal inversion
      else

        call gauss(Input,GeomI,Geom,1,.True.,.True.,Flgsg)

      end if ! Thermal or not-only-thermal inversion

      ! Master verbose
      if (gpid.eq.0) then
        umsg = ' - Angular quadrature initialized'
        call verbosev
      end if


      ! Define the output frequency axis
      call omegabuild(Frec,Atom,Input,maxB,Sol%omega_input)

      ! Control
      call gcontrol

      ! Master verbose
      if (gpid.eq.0) then
        write(umsg,'(" - Frequency axis initialized with",'// &
                   '1x,i6," frequencies")') nfreq
        call verbosev
      end if


      !
      ! Prepare limits for the PSF
      !
      if (gpid.gt.0.or..not.MPID%mpi15d) &
        call set_psf_ranges(Inf_Stokes,Input%lim_fwhm, &
                            Input%fwhm_fil,Sol%omega_input,Frec%omega)


      !
      ! Organize the tasks splitting
      !

      ! Not manager
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! Distribute
        call setmpi(MPID,Input,Frec%IW_freq)

        ! Allocate atomic MPI arrays
        call prepareatomMPI(Atom)
        call omegainitmaster(Atom)

      ! Global master
      else

        ! Remove fudge
        if (allocated(fudge%fudge_v)) then
          MRAMc = MRAMc - 1d-6*sizeof(fudge%fudge_v)
          deallocate(fudge%fudge_v)
        end if

        ! Verbose
        umsg = ' - Tasks within groups distributed'
        call verbosev

      end if ! Manager or not


      !
      ! Get Kurucz lines
      !
      if ((.not.MPID%mpi15d.and.(pid.gt.0.or..not.MPID%mpi)).or. &
          (gpid.gt.0.and.(pid.gt.0.or..not.MPID%mpi))) then

        ! Read Kurucz data
        if (Input%NK.ge.1) then

          ! Compute Doppler precursor
          DwTa = Input%dw*1d-9/c

          ! Get Kurucz data
          call kurucz_get(Atom,Atomb,Atmo_in,Input%LTEline, &
                          Input%kurucz,Input%NK, &
                          Frec%omega,MPID,DwTa,.False.,kurucz)

        ! If there are not
        else

          ! No Kurucz input
          kurucz%ntran = 0

        end if ! Kurucz data
      end if ! May need Kurucz data

      ! Control
      call gcontrol

      !
      ! Photoionization quantites
      !

      ! If not manager
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! For every atom
        do ia=1,nA

          ! Set constant part of photoionization information
          call setphoto(Atom(ia),Frec%omega)

        end do ! Atoms

      end if ! Not manager

      ! Master verbose
      if (gpid.eq.0.and.nA.gt.0) then
        umsg = ' - Initialized photoionization quantities '//&
               '(cross section)'
        call verbosev
      end if


      ! Not global master
      if (.not.MPID%mpi15d.or.gpid.gt.0) then

        ! Resize some frequency quantities if doing MPI
        call frecresize(Frec,Atom,Input,MPID)

        ! If RT MPI
        if (MPID%mpi) then

          ! If only thermal
          if (Input%Type_inversion.eq.0) then

            ! Compute size for MPI messages in solvers
            call setmpi_sizes(MPID,GeomI,Geom,Frec,.True.,.False., &
                              .True.,.False.,Input%ALI_photo,.False.)

          ! If full Stokes
          else

            ! Compute size for MPI messages in solvers
            call setmpi_sizes(MPID,GeomI,Geom,Frec,.True.,.True., &
                              .True.,.True.,Input%ALI_photo,.False.)

          end if ! Intensity/Stokes
        end if ! Freq. MPI
      end if ! Global master

      ! Determine if JKQ in the output
      Input%out_jkqa = s_jkq_buffer.gt.0.or. &
                       atmojkq.or. &
                       any(Inf_Nodes% &
                             Nodes_Flags(Inf_Nodes%index_J21R: &
                                         Inf_Nodes%index_J22I))

      ! Check if aborting
      call gcontrol

      !
      ! If there is input model here, prepare it
      !
      if ((.not.MPID%mpi15d.or.gpid.gt.0).and. &
          Input%atmoin_type.eq.0) then

        ! Set up input model
        call setup_Atmo_ininv(Atom,Atomb,Mol,Atmo_in,Input, &
                              fudge,Atmo_in%zalt,.True.)

        ! Initialize diffuse light
        if (Input%f_diff.gt.0d0) then
          Atmo_in%f_diff = Input%f_diff
        else
          Atmo_in%f_diff = 0d0
        end if

      end if ! Input model 1D

      ! Control
      call gcontrol

      ! Initialize buffer sizes (synthesis)
      call set_io_buffers(Input,0,Atom,Frec)

      ! Initialize buffer sizes (inversion)
      if (gpid.eq.0) then
        call set_inv_io_buffers(Input,Inf_Nodes,out_dims,nz)
      else
        call set_inv_io_buffers(Input,Inf_Nodes,dims,nz)
      end if

      ! If global master
      if (gpid.eq.0) then

        ! If existing cache
        if (lcache) then

          ! Check files exist
          call check_io_inv_buffer_exists(Input,aborting)

        end if ! Existing cache

        ! Get abortion flag
        laborted = aborting

      end if ! Master

      ! Check if aborting
      call gcontrol

      ! Tell verbosity that synthesis is in inversion
      ninv_mode = .False.

      ! Announce readiness
      if (gpid.eq.0) then
        umsg = ' - Ready to begin the inversion'
        call verbosev
      end if

      ! From now on, the synthesis output points to a different file
      verbosef = trim(verbosef)//'_syn'

      !
      ! Carry out inversion
      !

      !
      ! MPI version
      !
      if (MPID%mpi15d) then

        !
        ! Master
        !
        if (gpid.eq.0) then

          ! Open files so slaves can write later
          if (.not.lcache) then

            ! If constant mu
            if (finfo(2).eq.0) then

              ! Create files
              call create_io_inv_files(Input,Sol,Inf_Nodes, &
                                       acos(Inf_Stokes%mu)*180d0/PI, &
                                       Inf_Stokes%azimuth*180d0/PI, &
                                       out_dims,nz,Frec)

            ! Variable mu
            else

              ! Create files
              call create_io_inv_files(Input,Sol,Inf_Nodes, &
                                       -1d0,-1d0, &
                                       out_dims,nz,Frec)

            end if ! Type of mu
          end if ! No cache

          ! Allocate cpu_free with group status
          allocate(cpu_free(MPID%ngroup))
          MRAMc = MRAMc + 1d-6*sizeof(cpu_free)
          cpu_free = 1

          ! Initialize indexes and sizes
          ix = 1
          iy = 0
          ix1 = -1
          iy1 = -1
          ix2 = -1
          iy2 = -1
          inod = 0
          NLOSr = 0
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

                ! Get data
                call get_data_column(Input,unitD,dims, &
                              p_transfer_buffer(1:s_data_buffer), &
                                     finfo,ix,iy,check)

                ! Get column from atmosphere (input size)
                if (Input%atmoin_type.ge.1.and. &
                    s_atmo_buffer.gt.0.and.type_atmo_size.eq.1) then

                  ! Get input atmosphere data
                  ii = s_data_buffer + 1
                  jj = s_data_buffer + s_atmo_buffer
                  call get_column(unitA,p_transfer_buffer(ii:jj), &
                                  double,check)

                end if ! If input atmosphere data

                ! If JKQ data
                if (s_jkq_buffer.gt.0) then

                  ! Get input ad-hoc JKQ data
                  ii = s_data_buffer + s_atmo_buffer + 1
                  jj = s_transfer_buffer
                  call get_column(unitJ,p_transfer_buffer(ii:jj), &
                                  double_jkq,check)

                end if ! If JKQ data

                ! If mask unit
                if (unitM.gt.0) then

                  ! Get mask
                  read(unitM) imask

                ! No masks
                else

                  ! Flag unmasked
                  imask = 0

                end if ! Mask unit

                ! Store last read
                ix1 = ix
                iy1 = iy

              end if ! Loopìng

              ! Skip if out of box
              if (ix.lt.Input%sol_box(1).or. &
                  ix.gt.Input%sol_box(2).or. &
                  iy.lt.Input%sol_box(3).or. &
                  iy.gt.Input%sol_box(4)) cycle

              ! Check if not looping
              if (ix.ne.ix2.or.iy.ne.iy2) then

                ! Get column from atmosphere (output size)
                if (Input%atmoin_type.ge.1.and. &
                    s_atmo_buffer.gt.0.and.type_atmo_size.eq.2) then

                  ! Get input atmosphere data
                  ii = s_data_buffer + 1
                  jj = s_data_buffer + s_atmo_buffer
                  call get_column(unitA,p_transfer_buffer(ii:jj), &
                                  double,check)

                end if ! If input atmosphere data

                ! Store last read
                ix2 = ix
                iy2 = iy

              end if ! Not looping

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

              ! If excluding data
              if (Input%lexcl) then

                ! Check x coordinate could be excluded
                if (ix.ge.Input%excl(1,1).and. &
                    ix.le.Input%excl(1,Input%nexcl)) then

                  ! Initialize as no
                  lexcl = .False.

                  ! For all exclusion entries
                  do ia=1,Input%nexcl

                    ! If x coordinate below minimum x excluded, skip
                    if (Input%excl(1,ia).lt.ix) cycle

                    ! If x coordinate above maximum x excluded, leave
                    if (Input%excl(1,ia).gt.ix) exit

                    ! If y coordinate is excluded
                    if (Input%excl(2,ia).eq.iy) then

                      ! Flag as excluded
                      lexcl = .True.
                      exit

                    end if ! Excluded pixel

                  end do ! Exclusion entries

                  ! Skip if excluded pixel
                  if (lexcl) then
                    NLOSr = NLOSr + 1
                    cycle
                  end if ! Excluded pixel
                end if ! x coordinate could be excluded
              end if ! Exclusion data

              ! Take a free cpu
              ip = maxloc(cpu_free, 1)

              ! Message
              int_buff(1:3) = icoords
              int_buff(4) = imask

              ! Send signal to node
              call MPI_SEND(int_buff(1), 4, MPI_INTEGER, &
                            MPID%ltslave(ip), &
                            2+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, ierr)

              ! If failed, try again
              if (ierr.ne.0) then
                inod = inod - 1
                ix = ix0
                iy = iy0
                cycle
              end if

              ! Send data to node
              call MPI_SEND(p_transfer_buffer(1), &
                            s_transfer_buffer, &
                            MPI_DOUBLE_PRECISION, &
                            MPID%ltslave(ip), &
                            3+MPID%ltslave(ip), &
                            MPI_COMM_WORLD, ierr)

              ! If failed, try again
              if (ierr.ne.0) then
                inod = inod - 1
                ix = ix0
                iy = iy0
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
                              4+MPID%ltslave(ip), &
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
                                4+MPID%ltslave(ip), &
                                MPI_COMM_WORLD, &
                                MPI_STATUS_IGNORE, ierr)

                  ! If failed, try again
                  if (ierr.ne.0) cycle

                  exit

                end do ! Until received

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

              ! Advance group
              ip = ip + 1

              ! Leave if checked all groups
              if (ip.gt.MPID%ngroup) exit

            end do ! Slaves

            ! If we went beyond the number of LOS, exit
            if (NLOSr.ge.NLOS) exit

          end do ! While there is work to do

          ! If we are done, notify to slaves
          int_buff(1) = -1
          iproc = 1
          do while (.True.)

            ! send termination signal
            call MPI_SEND(int_buff(1), 4, MPI_INTEGER, &
                          MPID%ltslave(iproc), &
                          2+MPID%ltslave(iproc), &
                          MPI_COMM_WORLD, ierr)

            ! If it fails
            if (ierr.ne.0) cycle

            ! Advance group
            iproc = iproc + 1

            ! Leave if sent to all groups
            if (iproc.gt.MPID%ngroup) exit

          end do ! slaves

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

            ! If leader
            if (pid.eq.0) then

              ! Try receiving until success
              do while (.True.)

                ! Wait for signal
                call MPI_RECV(int_buff(1), 4, MPI_INTEGER, 0, &
                              2+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! Check if it is the termination signal
                if (int_buff(1).lt.1) then
                  aborting = .True.
                  exit
                end if

                ! Receive LOS
                call MPI_RECV(p_transfer_buffer(1), &
                              s_transfer_buffer, &
                              MPI_DOUBLE_PRECISION, &
                              0, 3+gpid, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! Success
                exit

              end do ! Until received successfully

            end if ! Leader

            ! If liutenant has friends
            if (nproc.gt.1) then

              ! Try until done
              do while (.True.)

                ! Broadcast
                call MPI_BCAST(int_buff(1), 4, MPI_INTEGER, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! If aborting
                if (int_buff(1).lt.1) then
                  aborting = .True.
                  exit
                end if

                ! Broadcast
                call MPI_BCAST(p_transfer_buffer(1), &
                               s_transfer_buffer, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

                ! If it fails
                if (ierr.ne.0) cycle

                ! Success
                exit

              end do ! Try communicating until successfull

            end if ! MPI in RT

            ! Get coords
            icoords = int_buff(1:3)
            imask = int_buff(4)

            ! Problem
            if (aborting) exit

            ! Initialize aborting flag
            laborted = .False.

            ! Get data from buffer
            call set_up_data_frombuffer(finfo,Inf_Stokes,Sol, &
                                        Input%Sigma_factor, &
                                   p_transfer_buffer(1:s_data_buffer))

            ! If restoring
            if (restoring) then

              ! Put model in atmo
              ii = s_data_buffer + 1
              jj = s_data_buffer + s_atmo_buffer
              call set_up_atmo_frombuffer(.True.,atmojkq,.True., &
                                          p_transfer_buffer(ii:jj), &
                                          Atmo_in,Bfield_in)

            ! Not restoring, but with input model
            else if (s_atmo_buffer.gt.0.or. &
                     s_jkq_buffer.gt.0) then

              ! Buffer limits
              ii = s_data_buffer + 1
              jj = s_data_buffer + s_atmo_buffer

              ! Standard 1.5DS model input
              if (Input%atmoin_type.eq.1) then

                ! Update atmospheric type
                Atmo_in%typo = Input%atmo_char

                ! Put model in atmo
                call set_up_atmo_frombuffer(.False.,.False., &
                                            Atmo_in%scal.eq.'T', &
                                            p_transfer_buffer(ii:jj),&
                                            Atmo_in,Bfield_in)
                !
                ! Set-up input model
                !
                call setup_Atmo_ininv(Atom,Atomb,Mol,Atmo_in, &
                                      Input,fudge, &
                                      Atmo_in%zalt,.False.)

              ! Inversion result
              else if (s_atmo_buffer.gt.0) then

                ! Put model in atmo
                call set_up_atmo_frombuffer(.True.,atmojkq,.True., &
                                            p_transfer_buffer(ii:jj),&
                                            Atmo_in,Bfield_in)

              end if

              ! Initialize diffuse light
              if (Input%f_diff.gt.0d0) then
                Atmo_in%f_diff = Input%f_diff
              else
                Atmo_in%f_diff = 0d0
              end if

              ! If input JKQ from file
              if (s_jkq_buffer.gt.0) then

                ! Buffer limits
                ii = s_data_buffer + s_atmo_buffer + 1
                jj = s_transfer_buffer

                ! Put JKQ in model Atmosphere
                call set_up_JKQ_frombuffer(p_transfer_buffer(ii:jj), &
                                           Atmo_in)

              end if
            end if ! Type of input atmosphere

            ! Need to update T
            if (update_tlim) then

              ! Temperature limits
              Input%minT = minval(Atmo_in%T)
              Input%maxT = maxval(Atmo_in%T)

            end if

            ! Need to update v
            if (update_vlim) then

              ! Velocity limit
              DwTa = maxval(sqrt(Atmo_in%vx*Atmo_in%vx + &
                                 Atmo_in%vy*Atmo_in%vy + &
                                 Atmo_in%vz*Atmo_in%vz))

              ! If maximum velocity
              if (DwTa.gt.0d0) then
                dyn = .True.
                Input%maxV = DwTa
                Input%static = .False.
              else
                dyn = .False.
                Input%maxV = 0d0
                Input%static = .True.
              end if ! Determine dynamic
            end if ! If we need to update vlim

            ! Need to update B
            if (update_blim) then

              ! Bield limit
              maxB = maxval(Bfield_in%Bstrength)

              ! Magnetic field?
              if (maxB.le.TINYB) then
                Input%unmagnetized = .True.
              else
                Input%unmagnetized = .False.
              end if

            end if

            ! Skip first
            if (lMRAMc.lt.-2d0) then

              ! Upgrade
              lMRAMc = -1.5d0

            ! If no lMRAMc data
            else if (lMRAMc.lt.0d0) then

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
                  urou = 'TIC'
                  write(umsg,'(2(A,es13.6),A)') &
                    'The miscellaneous RAM counter is different '// &
                    'between calls to the Inversion function ', &
                    MRAMc,' != ',lMRAMc,'. It is being '// &
                    'corrected, but this should not happen. '// &
                    'Please, notify of the issue providing '// &
                    'your inputs'
                  call abortedS(umsg,urou,.False.,.True.)

                end if ! Can issue warning

                ! Correct
                MRAMc = lMRAMc

              end if ! Different
            end if ! lMRAMc data

            ! Carry out the inversion
            call Inversion(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec, &
                           fudge,kurucz,MPID,Atmo_in, &
                           Bfield_in,Input,Inf_Stokes, &
                           Inf_Nodes,Sol,imask)

            ! Clean atmosphere?
            if (s_atmo_buffer.gt.0.or.restoring) &
              call free_Atmo(Atmo_in,.False.)

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
                              4+gpid,MPI_COMM_WORLD,ierr)
                if (ierr.ne.0) cycle
                exit
              end do

            end if ! Send info back to grand master

          end do ! Working loop

        end if ! Master or slave

      ! Serial version
      else

        ! Open files so slaves can write later
        if (.not.lcache) then

          ! If constant mu
          if (finfo(2).eq.0) then

            ! Create files
            call create_io_inv_files(Input,Sol,Inf_Nodes, &
                                     acos(Inf_Stokes%mu)*180d0/PI, &
                                     Inf_Stokes%azimuth*180d0/PI, &
                                     out_dims,nz,Frec)

          ! Variable mu
          else

            ! Create files
            call create_io_inv_files(Input,Sol,Inf_Nodes, &
                                     -1d0,-1d0, &
                                     out_dims,nz,Frec)

          end if ! Type of mu
        end if ! Files must exist

        ! Initialize indexes and sizes
        ix = 1
        iy = 0
        inod = 0
        NLOSr = 0
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

          ! If there are no more LOS to do, leave
          if (inod.ge.NLOS) exit

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

            ! Get data
            call get_data_column(Input,unitD,dims, &
                              p_transfer_buffer(1:s_data_buffer), &
                              finfo,ix,iy,check)

            ! Get column from atmosphere
            if (Input%atmoin_type.ge.1.and. &
                s_atmo_buffer.gt.0) then

              ii = s_data_buffer+1
              jj = s_data_buffer + s_atmo_buffer
              call get_column(unitA,p_transfer_buffer(ii:jj), &
                              double,check)

            end if

            ! If JKQ data
            if (s_jkq_buffer.gt.0) then

              ii = s_data_buffer + s_atmo_buffer + 1
              jj = s_transfer_buffer
              call get_column(unitJ,p_transfer_buffer(ii:jj), &
                              double_jkq,check)
            end if

            ! Get mask
            if (unitM.gt.0) then
              read(unitM) imask
            else
              imask = 0
            end if

            ! Store last read
            ix1 = ix
            iy1 = iy

          end if ! Loopìng

          ! Check if out of box
          if (ix.lt.Input%sol_box(1).or. &
              ix.gt.Input%sol_box(2).or. &
              iy.lt.Input%sol_box(3).or. &
              iy.gt.Input%sol_box(4)) cycle

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

          ! Initialize aborting flag
          laborted = .False.

          ! Get data from buffer
          call set_up_data_frombuffer(finfo,Inf_Stokes,Sol, &
                                      Input%Sigma_factor, &
                                   p_transfer_buffer(1:s_data_buffer))

          ! If restoring
          if (restoring) then

            ! Put model in atmo
            ii = s_data_buffer + 1
            jj = s_data_buffer + s_atmo_buffer
            call set_up_atmo_frombuffer(.True.,atmojkq,.True., &
                                        p_transfer_buffer(ii:jj), &
                                        Atmo_in,Bfield_in)

          ! Not restoring, but with input model
          else if (s_atmo_buffer.gt.0.or. &
                   s_jkq_buffer.gt.0) then

            ! Buffer limits
            ii = s_data_buffer + 1
            jj = s_data_buffer + s_atmo_buffer

            ! Standard 1.5DS model input
            if (Input%atmoin_type.eq.1) then

              ! Update atmospheric type
              Atmo_in%typo = Input%atmo_char

              ! Put model in atmo
              call set_up_atmo_frombuffer(.False.,.False., &
                                          Atmo_in%scal.eq.'T', &
                                          p_transfer_buffer(ii:jj), &
                                          Atmo_in,Bfield_in)
              !
              ! Set-up input model
              !
              call setup_Atmo_ininv(Atom,Atomb,Mol,Atmo_in, &
                                    Input,fudge, &
                                    Atmo_in%zalt,.False.)

            ! Inversion result
            else if (s_atmo_buffer.gt.0) then

              ! Put model in atmo
              call set_up_atmo_frombuffer(.True.,atmojkq,.True., &
                                          p_transfer_buffer(ii:jj), &
                                          Atmo_in,Bfield_in)

            end if

            ! Initialize diffuse light
            if (Input%f_diff.gt.0d0) then
              Atmo_in%f_diff = Input%f_diff
            else
              Atmo_in%f_diff = 0d0
            end if

            ! If input JKQ from file
            if (s_jkq_buffer.gt.0) then

              ! Buffer limits
              ii = s_data_buffer + s_atmo_buffer + 1
              jj = s_transfer_buffer

              ! Put JKQ in model Atmosphere
              call set_up_JKQ_frombuffer(p_transfer_buffer(ii:jj), &
                                         Atmo_in)

            end if

          end if ! Type of input atmosphere

          ! Need to update T
          if (update_tlim) then

            ! Temperature limits
            Input%minT = minval(Atmo_in%T)
            Input%maxT = maxval(Atmo_in%T)

          end if

          ! Need to update v
          if (update_vlim) then

            ! Velocity limit
            DwTa = maxval(sqrt(Atmo_in%vx*Atmo_in%vx + &
                               Atmo_in%vy*Atmo_in%vy + &
                               Atmo_in%vz*Atmo_in%vz))

            ! If maximum velocity
            if (DwTa.gt.0d0) then
              dyn = .True.
              Input%maxV = DwTa
              Input%static = .False.
            else
              dyn = .False.
              Input%maxV = 0d0
              Input%static = .True.
            end if
          end if

          ! Need to update B
          if (update_blim) then

            ! Bield limit
            maxB = maxval(Bfield_in%Bstrength)

            ! Magnetic field?
            if (maxB.le.TINYB) then
              Input%unmagnetized = .True.
            else
              Input%unmagnetized = .False.
            end if

          end if

          ! Skip first
          if (lMRAMc.lt.-2d0) then

            ! Upgrade
            lMRAMc = -1.5d0

          ! If no lMRAMc data
          else if (lMRAMc.lt.0d0) then

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
                urou = 'TIC'
                write(umsg,'(2(A,es13.6),A)') &
                  'The miscellaneous RAM counter is different '// &
                  'between calls to the Inversion function ', &
                  MRAMc,' != ',lMRAMc,'. It is being '// &
                  'corrected, but this should not happen. '// &
                  'Please, notify of the issue providing '// &
                  'your inputs'
                call abortedS(umsg,urou,.False.,.True.)

              end if ! Can issue warning

              ! Correct
              MRAMc = lMRAMc

            end if ! Different
          end if ! lMRAMc data

          ! Carry out the inversion
          call Inversion(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec, &
                         fudge,kurucz,MPID,Atmo_in, &
                         Bfield_in,Input,Inf_Stokes, &
                         Inf_Nodes,Sol,imask)

          ! Clean Atmo?
          if (s_atmo_buffer.gt.0.or.restoring) &
            call free_Atmo(Atmo_in,.False.)

          ! Copy coords
          int_buff(1:3) = icoords

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
      if (gpid.eq.0.and.unitA.gt.0) call close_file(unitA)
      if (gpid.eq.0.and.unitD.gt.0) call close_file(unitD)
      if (gpid.eq.0.and.unitC.gt.0) call close_file(unitC)
      if (gpid.eq.0.and.unitJ.gt.0) call close_file(unitJ)
      if (gpid.eq.0.and.unitM.gt.0) call close_file(unitM)


      !
      ! Restore verbosity name
      !
      verbosef = verbosefv

      ! Switch off inversion mode
      ninv_mode = .True.

      end subroutine TIC

!#####################################################################
!#####################################################################
!#####################################################################

      end module tic_mod
