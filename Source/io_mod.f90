      !> Input/Ouput for 1.5D
      module io_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     29/06/2022
!  Last version:
!     11/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  get_atmo_type:
!    Determine to what format of model atmosphere corresponds to the
!  specified file
!
!  open_atm:
!    Open the atmospheric model file and determine the success
!
!  open_atm_and_cache:
!    Open the atmospheric model file and the cache file if exists,
!  and determine the success
!
!  open_asymm
!    Open the ad-hoc radiation field tensor file and determine the
!  success
!
!  open_file
!    Open any binary file and determine success
!
!  close_file
!    Close the specified unit
!
!  get_dims
!    Read three integers from a file that are expected to be the
!  model dimensions and check that they are reasonable
!
!  get_dims_info
!    Read the dimension and the data properties in a data file for
!  the inversion
!
!  check_reasonable
!    Check if a given integer is reasonable for being the size of a
!  discretized axis
!
!  get_axes
!    Read the cartesian axes of a model atmosphere for the CLE mode
!
!  get_point
!    Read pixel data for a non-cartesian CLE model atmosphere
!
!  get_column
!    Read data from an opened file into the given buffer
!
!  get_lims
!    Check the temperature and velocity limits for an atmospheric
!  model
!
!  get_column_ion
!    Read data about ionization in the given buffer after the
!  specified offset
!
!  get_cache
!    Read the cache from an existing cache file
!
!  start_cache
!    Initialize the cache file
!
!  write_cache
!    Write entry into the cache file
!
!  prepare_lambda_limits
!    Allocate and setup the index limits in the wavelength dependent
!  variables to cut the synthesis input based on the wavelengths in
!  the inversion data
!
!  set_lambda_limit
!    Specify the index limit to choose what to write in the output
!  of a wavelength dependent variable
!
!  set_cols_limit
!    Specify the buffer size to write in the output for collisions
!
!  set_damp_limit
!    Specify the buffer size to write in the output for damping
!  parameters
!
!  set_pop_limit
!    Specify the buffer size to write in the output for populations
!  and departure coefficients
!
!  set_io_buffers
!    Setup buffer sizes and limits to write the output in 1.5D
!  synthesis mode
!
!  check_cols_limit
!    Checks that the specified limits for the collisional output
!  comply with the atomic models
!
!  check_damp_limit
!    Checks that the specified limits for the damping parameter output
!  comply with the atomic models
!
!  check_qel_limit
!    Checks that the specified limits for the elastic rates output
!  comply with the atomic models
!
!  check_pop_limit
!    Checks that the specified limits for the population output
!  comply with the atomic models
!
!  check_io_buffers_sanity_check
!    Call the check_*_limit routines in sequence to check that
!  the specified limits make sense for the actual atomic models
!
!  check_io_buffers_exist
!    Check if the output files for the 1.5D synthesis already exist
!
!  create_io_files
!    Create the output files for the 1.5D synthesis and write their
!  headers
!
!  set_io_CLE_buffers
!    Setup buffer sizes and limits to write the output in CLE
!  synthesis mode
!
!  check_io_CLE_buffers_exist
!    Check if the output files for the CLE synthesis already exist
!
!  create_io_CLE_files:
!    Create the output files for the CLE synthesis and write their
!  headers
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : PI
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine to what format of model atmosphere corresponds to
      !! the specified file\n
      !!  filename(character(500)): Filename\n
      !!      atmoin_type(integer): Type of model\n
      !!                             1: 1.5D synthesis model\n
      !!                             2: 1.5D inversion model\n
      !!                            -1: Not a valid model
      subroutine get_atmo_type(filename,atmoin_type)

      ! I/O

      character(len=500), intent(in):: filename
      integer, intent(out):: atmoin_type

      ! Local

      character(len=4):: label

      logical:: check

      integer:: unitA


      ! Unit (local scope)
      unitA = 16

      ! Try opening file
      call open_file(unitA, filename, 0, .False., check)

      ! Read label
      read(unitA, err=1100) label

      ! If normal 3D
      if (label.eq.'2Dat') then

        ! Flag
        atmoin_type = 1

      ! If inversion 3D
      else if (label.eq.'invo') then

        ! Flag
        atmoin_type = 2

      ! None of the above
      else

        ! Flag
        atmoin_type = -1

      end if

      ! Close
      call close_file(unitA)

      return

      ! Error
1100  atmoin_type = -1
      return

      end subroutine get_atmo_type

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open the atmospheric model file and determine the success\n
      !!  Input(Input_class): Structure with configuration data\n
      !!   run_mode(integer): Mode of the run\n
      !!      unitA(integer): Unit where to open atmospheric model\n
      !!   aborting(logical): Indicate failure at output\n
      !!    dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!       mode(integer): Type of atmospheric model\n
      !!     double(logical): If data in double precision\n
      !!       norm(integer): If the geometrical axes are normalized
      subroutine open_atm(Input,run_mode,unitA,aborting,dims,mode, &
                          double,norm)

      ! I/O

      type(Input_class), intent(in):: Input
      logical, intent(out):: aborting,double
      integer, intent(in):: run_mode
      integer, intent(out):: unitA,mode,norm
      integer, dimension(:), intent(out):: dims

      ! Local

      logical:: check


      !
      ! Give unit a number
      !
      unitA = 16

      !
      ! Open atmospheric file to read
      !
      call open_file(unitA, Input%atmo, 0, .False., check)

      ! Check could open
      if (.not.check) then

        ! Flag error and return
        aborting = .True.
        return

      end if

      ! Get dimensions from atmosphere file
      call get_dims(unitA,run_mode,mode,double,norm,dims,check)

      ! Check could read
      if (.not.check) then

        ! Flag error and return
        aborting = .True.
        return

      end if

      ! Verbose
      write(umsg,'(A)') ' - Reading atmospheric file '// &
                        trim(Input%atmo)
      call verbose

      return

      end subroutine open_atm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open the atmospheric model file and the cache file if
      !! exists, and determine the success\n
      !!   Input(Input_class): Structure with configuration data\n
      !!    run_mode(integer): Mode of the run\n
      !!       unitA(integer): Unit where to open atmospheric model\n
      !!       unitC(integer): Unit where to open cache\n
      !!    aborting(logical): Indicate failure at output\n
      !!     dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!        mode(integer): Type of atmospheric model\n
      !!      double(logical): If data in double precision\n
      !!        norm(integer): If the geometrical axes are
      !!                       normalized\n
      !!  cache(logical(:,:)): Cache of already done columns\n
      !!      lcache(logical): If there is a cache at output
      subroutine open_atm_and_cache(Input,run_mode,unitA,unitC, &
                                    aborting,dims,mode,double,norm, &
                                    cache,lcache)

      ! I/O

      type(Input_class), intent(in):: Input
      logical, intent(out):: aborting,lcache,double
      logical, dimension(:,:), allocatable, intent(out):: cache
      integer, intent(in):: run_mode
      integer, intent(out):: unitA,unitC,mode,norm
      integer, dimension(:), intent(out):: dims

      ! Local

      logical:: check


      !
      ! Give units a number and order to open
      !

      ! Unit of atmosphere file
      unitA = 16

      ! Unit of cache file
      unitC = 17

      !
      ! Open atmospheric file to read
      !
      call open_file(unitA, Input%atmo, 0, .False., check)

      ! Check could open
      if (.not.check) then

        ! Flag error and return
        aborting = .True.
        return

      end if

      ! Get dimensions from atmosphere file
      call get_dims(unitA,run_mode,mode,double,norm,dims,check)

      ! Check could open
      if (.not.check) then

        ! Flag error and return
        aborting = .True.
        return

      end if

      !
      ! Check if there is a cache
      !

      ! Check if there is a cache file
      inquire(file=trim(Input%cache), exist=lcache)

      ! If there is a cache
      if (lcache) then

        ! 1.5D
        if (run_mode.eq.1) then

          ! Allocate cache data
          allocate(cache(dims(2),dims(1)))

          ! Read cache
          call get_cache(unitC,Input%cache,dims(1:2),cache,check)

        ! CLE
        else if (run_mode.eq.2) then

          ! Allocate cache data
          allocate(cache(dims(3),dims(2)))

          ! Read cache
          call get_cache(unitC,Input%cache,dims(2:3),cache,check)

        end if

        ! Check if valid cache
        lcache = check

        ! If there is a cache file
        if (lcache) then

          ! If nothing actually done, ignore cache
          if (.not.any(cache)) lcache = .False.

        end if

        ! If there is cache afterall
        if (lcache) then

          ! Count memory
          MRAMc = MRAMc + 1d-6*sizeof(cache)

        ! There is no cache
        else

          ! Free its space
          deallocate(cache)

        end if ! There is cache, actually
      endif ! If there is cache

      ! 1.5D
      if (run_mode.eq.1) then

        ! Open or initialize cache
        call start_cache(unitC,Input%cache,dims(1:2),lcache,check)

        ! Check success
        aborting = .not.check

      ! CLE
      else if (run_mode.eq.2) then

        ! Open or initialize cache
        call start_cache(unitC,Input%cache,dims(2:3),lcache,check)

        ! Check success
        aborting = .not.check

      end if ! 1.5D or CLE

      ! Verbose
      write(umsg,'(A)') ' - Reading atmospheric file '// &
                        trim(Input%atmo)
      call verbose

      return

      end subroutine open_atm_and_cache

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open the ad-hoc radiation field tensor file and determine the
      !! success\n
      !!   Input(Input_class): Structure with configuration data\n
      !!       unitH(integer): Unit where to open asymmetry\n
      !!    aborting(logical): Indicate failure at output\n
      !!     dims(integer(:)): Grid dimensions (X,Y,Z)
      subroutine open_asymm(Input,unitJ,aborting,dims)

      ! I/O

      type(Input_class), intent(in):: Input
      logical, intent(out):: aborting
      integer, intent(out):: unitJ
      integer, dimension(:), intent(in):: dims

      ! Local

      character(len=3):: label

      logical:: check

      integer:: ii,prec
      integer, dimension(3):: ldims

      !
      ! Give units a number and order to open
      !

      ! Unit of JKQ asymmetry file
      unitJ = 18

      !
      ! Open file to read
      !
      call open_file(unitJ, Input%asym_fil(1)%str, 0, .False., check)

      ! Check could open
      if (.not.check) then

        ! Flag error and return
        aborting = .True.
        return

      end if

      ! Read label
      read(unitJ, err=1100) label

      ! Check label
      if (label.ne.'JKQ') then

        ! Issue error
        urou = 'open_asymm'
        write(umsg,'(A)') ' # Wrong file identifier in '// &
                          'the asymmetry file, '// &
                          'expected "JKQ" and got "'// &
                          label//'"'
        call verbose
        aborting = .True.
        return

      end if

      ! Read precision
      read(unitJ, err=1100) prec

      ! Check if not double
      if (prec.ne.8) then

        ! Issue error
        urou = 'open_asymm'
        write(umsg,'(A)') ' # The asymmetry file only '// &
                          'admits double precision right '// &
                          'now'
        call verbose
        aborting = .True.
        return

      end if

      ! Read integers
      read(unitJ, err=1100) ldims(1:3)

      ! For each dimension
      do ii=1,3

        ! Check if different size
        if (dims(ii).ne.ldims(ii)) then

          ! Issue error
          urou = 'open_asymm'
          write(umsg,'(A,3(1x,i4),A,3(1x,i4),A)') &
                            ' # The asymmetry file has '// &
                            'different dimensions (',ldims, &
                            ') than the model atmosphere (', &
                            dims,')'
          call verbose
          aborting = .True.
          return

        end if ! Different size

      end do ! Dimensions

      ! Verbose
      write(umsg,'(A)') ' - Reading asymmetry file '// &
                        trim(Input%asym_fil(1)%str)
      call verbose

      return

1100  aborting = .True.
      write(umsg,'(A)') ' # Error reading from asymmetry file'
      call verbose

      return

      end subroutine open_asymm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open any binary file and determine success\n
      !!     unit_index(integer): Unit where to open the file\n
      !!  filename(character(:)): Name of the file to open\n
      !!     write_flag(integer): Specifies if it is opening for
      !!                          writing\n
      !!         append(logical): Specifies if it is appending\n
      !!          check(logical): Success flag
      subroutine open_file(unit_index,filename,write_flag,append, &
                           check)

      ! I/O

      character(len=500), intent(in):: filename
      logical, intent(in):: append
      logical, intent(out):: check
      integer, intent(in):: unit_index, write_flag

      ! Local

      character(len=5), allocatable:: caction
      character(len=6), allocatable:: cpos

      integer:: ios


      !
      ! Define caction and cpos
      !

      ! If not writing
      if (write_flag.eq.0) then

        ! Action is read
        caction = 'read '

        ! If appending
        if (append) then

          ! Set append
          cpos = 'APPEND'

        ! Otherwise
        else

          ! Not appending
          cpos = 'ASIS  '

        end if ! Appending

      ! Writing
      else

        ! Set action
        caction = 'write'
        cpos = 'ASIS'

      end if ! Read or writing

      ! Open file
      open (unit_index, file=trim(filename), status='unknown', &
            iostat=ios, err=1000, access='stream', &
            action=trim(caction), form='unformatted', &
            position=trim(cpos))

      ! If reached this, it was a success
      check = .True.

      return

1000  check = .False.
      write(umsg,'(A)') ' # Error opening file '//trim(filename)
      call verbose

      return

      end subroutine open_file

!#####################################################################
!#####################################################################
!#####################################################################

      !> Close the specified unit\n
      !!   unit_index(integer): Unit to close
      subroutine close_file(unit_index)

      ! I/O

      integer, intent(in):: unit_index


      ! Close the unit
      close(unit_index)

      end subroutine close_file

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read three integers from a file that are expected to be the
      !! model dimensions and check that they are reasonable\n
      !!  unit_index(integer): Unit to read from\n
      !!    run_mode(integer): Mode of the run\n
      !!        mode(integer): Type of atmospheric model\n
      !!      double(logical): If data in double precision\n
      !!        norm(integer): If the geometrical axes are
      !!                       normalized\n
      !!     dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!       check(logical): Success flag
      subroutine get_dims(unit_index,run_mode,mode,double,norm,dims, &
                          check)

      ! I/O

      logical, intent(out):: check,double
      integer, intent(in):: unit_index,run_mode
      integer, intent(out):: mode,norm
      integer, dimension(:), intent(out):: dims

      ! Local

      character(len=3):: lavel
      character(len=4):: label

      integer:: prec


      ! 1.5D synthesis or inversion
      if (run_mode.eq.1.or.run_mode.eq.-1) then

        ! Read label
        read(unit_index, err=1100) label

        ! Check label
        if (label.ne.'2Dat'.and.label.ne.'invo') then

          ! Issue error
          check = .False.
          urou = 'get_dims'
          write(umsg,'(A)') ' # Wrong file identifier in '// &
                            'the atmospheric model, '// &
                            'expected "2Dat" or "invo" and got "'// &
                            label//'"'
          call verbose
          return

        end if

        ! If 1.5DS model
        if (label.eq.'2Dat') then

          ! Read precision
          read(unit_index, err=1100) prec

          ! Unused variables
          mode = -1
          norm = -1

        ! Inversion output
        else

          ! Read type data
          read(unit_index, err=1100) mode

          ! Force single precision
          prec = 4

          ! Unused variable
          norm = -1

        end if ! Type of input

        ! Check precision
        if (prec.ne.4.and.prec.ne.8) then

          ! Issue error
          check = .False.
          urou = 'get_dims'
          write(umsg,'(A,i6)') ' # Precision must be 4 or 8 '// &
                               'and got ',prec
          call verbose
          return

        end if

        ! Read integers with sizes
        read(unit_index, err=1100) dims(1:3)

        !
        ! Check dimensions
        !

        ! X
        call check_reasonable(dims(1),'NX',check)
        if (.not.check) return

        ! Y
        call check_reasonable(dims(2),'NY',check)
        if (.not.check) return

        ! Z
        call check_reasonable(dims(3),'NZ',check)
        if (.not.check) return

      ! CLE synthesis
      else if (run_mode.eq.2) then

        ! Read label
        read(unit_index, err=1100) lavel

        ! Check label
        if (lavel.ne.'CLE') then

          ! Issue error
          check = .False.
          urou = 'get_dims'
          write(umsg,'(A)') ' # Wrong file identifier in '// &
                            'the atmospheric model, '// &
                            'expected "CLE" and got "'// &
                            lavel//'"'
          call verbose
          return

        end if

        ! Read precision
        read(unit_index, err=1100) prec

        ! Check precision
        if (prec.ne.4.and.prec.ne.8) then

          ! Issue error
          check = .False.
          urou = 'get_dims'
          write(umsg,'(A,i6)') ' # Precision must be 4 or 8 and '// &
                               'got ',prec
          call verbose
          return

        end if

        ! Read type of atmospheric file
        read(unit_index, err=1100) mode
        read(unit_index, err=1100) norm

        ! LOS cartesian model
        if (mode.eq.0) then

          ! Announce
          umsg = ' - LOS cartesian mode of calculation'
          call verbose

          ! Read dimensions
          read(unit_index, err=1100) dims

          !
          ! Check dimensions
          !

          ! X
          call check_reasonable(dims(1),'NX',check)
          if (.not.check) return

          ! Y
          call check_reasonable(dims(2),'NY',check)
          if (.not.check) return

          ! Z
          call check_reasonable(dims(3),'NZ',check)
          if (.not.check) return

        ! Slab mode
        else if (mode.eq.1) then

          ! Announce
          umsg = ' - Slab mode of calculation'
          call verbose

          ! Get dimension
          read(unit_index,err=1100) dims(2)

          ! Unused dimensions
          dims(1) = 1
          dims(3) = 1

          ! Check dimension
          call check_reasonable(dims(2),'NS',check)
          if (.not.check) return

        ! LOS non-cartesian
        else if (mode.eq.2) then

          ! Announce
          umsg = ' - LOS non-cartesian mode of calculation'
          call verbose

          ! Get dimensions
          read(unit_index, err=1100) dims(2)

          ! Unused variables
          dims(1) = 1
          dims(3) = 1

          ! Check dimension
          call check_reasonable(dims(2),'NL',check)
          if (.not.check) return

        ! No recognized mode
        else

          ! Issue error
          check = .False.
          umsg = ' # Atmospheric model mode not recognized'
          call verbose
          return

        end if ! CLE model type

      ! Unexpected running mode
      else

        ! Return failure
        check = .False.
        urou = 'get_dims'
        write(umsg,'(A,i1,A)') ' # Wrong running mode in '// &
                               'get_dims()',run_mode,' != 1 or 2'
        call verbose
        return

      end if ! Running mode

      ! Set double precision flag
      double = prec.eq.8

      ! Success if reached this
      check = .True.

      ! Return
      return

1100  check = .False.
      write(umsg,'(A)') ' # Error getting dimensions from'// &
                        ' the atmospheric file'
      call verbose

      return

      end subroutine get_dims

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the dimension and the data properties in a data file for
      !! the inversion\n
      !!  unit_index(integer): Unit to read from\n
      !!     dims(integer(:)): Dimensions in data file\n
      !!    finfo(integer(:)): Input data information\n
      !!       check(logical): Success flag
      subroutine get_dims_info(unit_index,dims,finfo,check)

      ! I/O

      logical, intent(out):: check
      integer, intent(in):: unit_index
      integer, dimension(:), intent(out):: dims,finfo

      ! Local
      character(len=4):: label


      ! Read label
      read(unit_index, err=1100) label

      ! Check label
      if (label.ne.'invi') then

        ! Issue error
        check = .False.
        urou = 'get_dims_info'
        write(umsg,'(A)') ' # Wrong file identifier in '// &
                          'the input inversion file, '// &
                          'expected "invi" and got "'// &
                          label//'"'
        call verbose
        return

      end if

      ! Read integers (nx,ny,nl)
      read(unit_index, err=1100) dims(1:3)

      !
      ! Check dimensions
      !

      ! X
      call check_reasonable(dims(1),'NX',check)
      if (.not.check) return

      ! Y
      call check_reasonable(dims(2),'NY',check)
      if (.not.check) return

      ! Wavelength
      call check_reasonable(dims(3),'NL',check)
      if (.not.check) return

      ! Read Stokes type, LOS type, sigma type, and diffuse light type
      read(unit_index, err=1100) finfo(1:4)

      !
      ! Check values
      !

      ! Stokes
      if (finfo(1).lt.0.or.finfo(1).gt.1) then

        ! Problem
        urou = 'get_dims_info'
        umsg = ' # Stokes label must be 0 or 1'
        call verbose
        check = .False.
        return

      end if

      ! LOS type
      if (finfo(2).lt.0.or.finfo(2).gt.1) then

        ! Problem
        urou = 'get_dims_info'
        umsg = ' # LOS label must be 0 or 1'
        call verbose
        check = .False.
        return

      end if

      ! Sigma type
      if (finfo(3).lt.0.or.finfo(3).gt.4) then

        ! Problem
        urou = 'get_dims_info'
        umsg = ' # Sigma label must be 0, 1, 2, 3, or 4'
        call verbose
        check = .False.
        return

      end if

      ! Diffuse light type
      if (finfo(4).lt.0.or.finfo(4).gt.4) then

        ! Problem
        urou = 'get_dims_info'
        umsg = ' # Diffuse light label must be 0, 1, 2, 3, or 4'
        call verbose
        check = .False.
        return

      end if

      ! Success if reached this
      check = .True.

      ! Return
      return

1100  check = .False.
      write(umsg,'(A)') ' # Error getting dimensions and file '// &
                        'information from the input inversion file'
      call verbose

      return

      end subroutine get_dims_info

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if a given integer is reasonable for being the size of
      !! a discretized axis\n
      !!         val(integer): Value of the integer to evaluate\n
      !!  label(character(2)): Axis name for verbosity\n
      !!       check(logical): Success flag
      subroutine check_reasonable(val,label,check)

      ! I/O

      character(len=2), intent(in):: label
      integer, intent(in):: val
      logical, intent(out):: check

      ! Local

      integer, parameter:: reason=10000


      ! Wrong because too small or too big
      if (val.lt.1.or.val.gt.reason) then

        ! Issue error
        check = .False.
        urou = 'get_dims'

        ! Inversion
        if (run_mode.eq.-1) then

          ! Message
          umsg = ' # Dimension in the input inversion '// &
                 'file out of reasonable limits, probably '// &
                 'wrong file.'

        ! Synthesis
        else

          ! Message
          umsg = ' # Dimension in '// &
                 'the atmospheric file out of '// &
                 'reasonable limits, probably '// &
                 'wrong model.'
        end if

        ! Verbose message
        call verbose

        ! Complete message with wrong dimension
        write(umsg,'(A,A,A,1x,i6)') ' ',label,' =',val
        call verbose

        ! Leave
        return

      end if

      ! Correct if reached this
      check = .True.

      return

      end subroutine check_reasonable

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the cartesian axes of a model atmosphere for the CLE
      !! mode\n
      !!  unit_index(integer): Unit to read from\n
      !!     dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!         x(double(:)): X axis\n
      !!         y(double(:)): Y axis\n
      !!         z(double(:)): Z axis\n
      !!      double(logical): If data in double precision\n
      !!     reading(logical): If actually reading the axes
      subroutine get_axes(unit_index, dims, x, y, z, double, reading)

      ! I/O

      logical, intent(in):: reading,double
      integer, intent(in):: unit_index
      integer, dimension(3), intent(in):: dims
      double precision, dimension(:), allocatable, intent(out):: x
      double precision, dimension(:), allocatable, intent(out):: y
      double precision, dimension(:), allocatable, intent(out):: z

      ! Local

      real, dimension(:), allocatable:: xf
      real, dimension(:), allocatable:: yf
      real, dimension(:), allocatable:: zf

      !
      ! Allocate axes
      !

      ! x
      allocate(x(dims(1)))
      MRAMc = MRAMc + 1d-6*sizeof(x)

      ! y
      allocate(y(dims(2)))
      MRAMc = MRAMc + 1d-6*sizeof(y)

      ! z
      allocate(z(dims(3)))
      MRAMc = MRAMc + 1d-6*sizeof(z)

      ! If actually reading, only the master should have this
      if (reading) then

        ! If double precision
        if (double) then

          ! Read directly into output variables
          read(unit_index, err=1100) x
          read(unit_index, err=1100) y
          read(unit_index, err=1100) z

        ! If single precision
        else

          ! Allocate auxiliar reals
          allocate(xf(dims(1)))
          allocate(yf(dims(2)))
          allocate(zf(dims(3)))

          ! Read reals
          read(unit_index, err=1100) xf
          read(unit_index, err=1100) yf
          read(unit_index, err=1100) zf

          ! Convert ot double
          x = dble(xf)
          y = dble(yf)
          z = dble(zf)

          ! Free
          deallocate(xf,yf,zf)

        end if ! Double precision
      end if ! Actually reading

      return

1100  umsg = ' # Error reading axes from the atmospheric file'
      call verbose

      return

      end subroutine get_axes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read pixel data for a non-cartesian CLE model atmosphere\n
      !!  unit_index(integer): Unit to read from\n
      !!          nx(integer): X axis size\n
      !!            y(double): Y position\n
      !!            z(double): Z position\n
      !!      double(logical): If data in double precision
      subroutine get_point(unit_index,nx,y,z,double)

      ! I/O

      logical, intent(in):: double
      integer, intent(in):: unit_index
      integer, intent(out):: nx
      double precision, intent(out):: y,z

      ! Local

      real:: yf,zf


      !
      ! Read y and z coordinates
      !

      ! If double precision
      if (double) then

        ! Read directly into double
        read(unit_index, err=1100) y
        read(unit_index, err=1100) z

      ! If single
      else

        ! Read into float
        read(unit_index, err=1100) yf
        read(unit_index, err=1100) zf

        ! Transform into double
        y = dble(yf)
        z = dble(zf)

      end if ! Precision

      ! Read size of LOS axis
      read(unit_index, err=1100) nx

      return

1100  umsg = ' # Error reading position from the atmospheric file'
      call verbose

      return

      end subroutine get_point

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read data from an opened file into the given buffer\n
      !!  unit_index(integer): Unit to read from\n
      !!    buffer(double(:)): Buffer to read the data into\n
      !!      double(logical): If data in double precision\n
      !!       check(logical): Success flag
      subroutine get_column(unit_index,buffer,double,check)

      ! I/O

      logical, intent(in):: double
      logical, intent(out):: check
      integer, intent(in):: unit_index
      double precision, dimension(:), intent(out):: buffer

      ! Local

      real, dimension(:), allocatable:: bufferr


      ! If double precision
      if (double) then

        ! Read directly into buffer
        read(unit_index, err=1100) buffer

      ! If single precision
      else

        ! Allocate float buffer
        allocate(bufferr(size(buffer)))

        ! Read into float buffer
        read(unit_index, err=1100) bufferr

        ! Convert to double precision
        buffer = dble(bufferr)

        ! Free float buffer
        deallocate(bufferr)

      end if ! Precision

      ! Success
      check = .True.

      return

1100  check = .False.
      umsg = ' # Error reading column from the atmospheric file'
      call verbose

      return

      end subroutine get_column

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check the temperature and velocity limits for an atmospheric
      !! model\n
      !!  Input(Input_class): Structure with configuration data\n
      !!   run_mode(integer): Mode of the run\n
      !!   aborting(logical): Indicate failure at output
      subroutine get_lims(Input,run_mode,aborting)

      ! I/O

      type(Input_class), intent(inout):: Input
      logical, intent(out):: aborting
      integer, intent(in):: run_mode

      ! Local

      logical:: double,check,readx,readxy

      integer:: unitA,mode,norm,sizeA,jump,ix,iy,iz,i0,i1
      integer:: d0,di0,di1,it,iv
      integer, dimension(3):: dims
      integer, dimension(4):: sol_box

      double precision:: minT,maxT,maxV,lv,ycoor,zcoor
      double precision, dimension(:), allocatable:: x,y,z
      double precision, dimension(:), allocatable, target:: buffer
      double precision, dimension(:), pointer:: p_T,p_vx,p_vy,p_vz

      !
      ! Only Master
      !
      if (gpid.eq.0) then

        ! Initialize good
        aborting = .False.

        ! Dummy loop
        do while (.True.)

          !
          ! Give unit a number (local scope)
          !
          unitA = 16

          !
          ! Open atmospheric file to read
          !
          call open_file(unitA, Input%atmo, 0, .False., check)

          ! Check could open
          if (.not.check) then

            ! Issue error
            aborting = .True.
            exit

          end if

          ! Get dimensions from atmosphere file
          call get_dims(unitA,run_mode,mode,double,norm,dims,check)

          ! Check could read
          if (.not.check) then

            ! Issue error
            aborting = .True.
            exit

          end if

          ! Initialize limits
          minT = 1d299
          maxT = -1d299
          if (Input%static) then
            maxV = 0d0
          else
            maxV = -1d299
          end if

          ! Nullify pointers
          nullify(p_T,p_vx,p_vy,p_vz)

          ! If 1.5D or inversion with 1.5D model
          if (run_mode.eq.1.or.(run_mode.eq.-1.and.mode.lt.0)) then

            ! If 1.5D synthesis
            if (run_mode.eq.1) then

              ! Get copy of solution box and fix wildcard indexes
              sol_box = Input%sol_box
              if (sol_box(1).lt.1) sol_box(1) = 1
              if (sol_box(2).lt.1) sol_box(2) = dims(1)
              if (sol_box(3).lt.1) sol_box(3) = 1
              if (sol_box(4).lt.1) sol_box(4) = dims(2)

            ! Otherwise
            else

              ! Set box from dimension
              sol_box = (/ 1, dims(1), 1, dims(2) /)

            end if ! Type of run

            ! Allocate column buffer
            sizeA = dims(3)*24
            allocate(buffer(sizeA))

            ! Specify amount to jump to skip columns
            jump = sizeA*8

            ! For each X
            do ix=1,dims(1)

              ! Check error
              if (aborting) exit

              ! Read?
              readx = ix.ge.sol_box(1).and.ix.le.sol_box(2)

              ! If below lower limit
              if (ix.lt.sol_box(1)) then

                ! Skip full line
                call fseek(unitA,jump*dims(2),1)
                cycle

              end if ! Below X lower limit

              ! If beyond limit, stop already
              if (ix.gt.sol_box(2)) exit

              ! For each Y
              do iy=1,dims(2)

                ! Read?
                readxy = readx.and. &
                         iy.ge.sol_box(3).and.iy.le.sol_box(4)

                ! If reading
                if (readxy) then

                  ! Get column
                  call get_column(unitA,buffer,double,check)

                  ! Check could read
                  if (.not.check) then

                    ! Issue error
                    aborting = .True.
                    exit

                  end if

                  ! Point to temperature and check limit
                  p_T => buffer(3*dims(3)+1:4*dims(3))
                  minT = min(minT,minval(p_T))
                  maxT = max(maxT,maxval(p_T))

                  ! Non-static
                  if (.not.Input%static) then

                    ! Point to velocity and update limit
                    p_vx => buffer(9*dims(3)+1:10*dims(3))
                    p_vy => buffer(10*dims(3)+1:11*dims(3))
                    p_vz => buffer(11*dims(3)+1:12*dims(3))
                    maxV = max(maxV,maxval(sqrt(p_vx*p_vx + &
                                                p_vy*p_vy + &
                                                p_vz*p_vz)))
                  end if

                ! Not reading
                else

                  ! Skip column
                  call fseek(unitA,jump,1)

                end if ! Reasing column

              end do ! Y
            end do ! X

          ! If inversion and previous solution
          else if (run_mode.eq.-1) then

            ! If jkq
            if (mode.gt.7) then

              ! Size atmosphere
              sizeA = dims(3)*27 + 1

            ! No jkq
            else

              ! Size atmosphere
              sizeA = dims(3)*19 + 1

            end if

            ! Allocate column buffer
            jump = sizeA*8
            allocate(buffer(sizeA))

            ! For each X
            do ix=1,dims(1)

              ! For each Y
              do iy=1,dims(2)

                ! Get column
                call get_column(unitA,buffer,double,check)

                ! Check could read
                if (.not.check) then

                  ! Issue error
                  aborting = .True.
                  exit

                end if

                ! Point to temperature and update limits
                p_T => buffer(  dims(3)+1:2*dims(3))
                minT = min(minT,minval(p_T))
                maxT = max(maxT,maxval(p_T))

                ! If non-static
                if (.not.Input%static) then

                  ! Point to velocity and update limit
                  p_vx => buffer(6*dims(3)+1:7*dims(3))
                  p_vy => buffer(7*dims(3)+1:8*dims(3))
                  p_vz => buffer(8*dims(3)+1:9*dims(3))
                  maxV = max(maxV,maxval(sqrt(p_vx*p_vx + &
                                              p_vy*p_vy + &
                                              p_vz*p_vz)))

                end if ! Non-static

              end do ! Y
            end do ! X

          ! If CLE
          else if (run_mode.eq.2) then

            ! If cartesian or slab
            if (mode.eq.0.or.mode.eq.1) then

              ! Get axes if cartesian
              if (mode.eq.0) &
                call get_axes(unitA,dims,x,y,z,double,gpid.eq.0)

              ! Size
              sizeA = dims(1)*22
              if (mode.eq.1) sizeA = sizeA + 3

              ! Cartesian
              if (mode.eq.0) then

                ! Indexes for buffer call
                d0 = dims(1)
                di0 = 1
                di1 = dims(1)
                it = 1
                iv = 7

              ! Slab
              else

                ! Indexes for buffer call
                d0 = 1
                di0 = 0
                di1 = 0
                it = 5
                iv = 11

              end if

            ! If not cartesian
            else if (mode.eq.2) then

              ! Size
              sizeA = dims(1)*23

              ! Indexes for buffer call
              d0 = dims(1)
              di0 = 1
              di1 = dims(1)
              it = 2
              iv = 8

            end if

            ! Allocate buffer
            allocate(buffer(sizeA))

            ! Static
            if (Input%static) then

              ! For each column
              do iy=1,dims(2)

                ! Aborting
                if (aborting) exit

                do iz=1,dims(3)

                  ! Not cartesian
                  if (mode.eq.2) &
                    call get_point(unitA,dims(1),ycoor,zcoor,double)

                  ! Get data
                  call get_column(unitA,buffer,double,check)

                  ! Check could read
                  if (.not.check) then
                    aborting = .True.
                    exit
                  end if

                  ! For each point
                  do ix=1,dims(1)

                    ! Coordinate
                    i1 = di0*ix

                    ! Update temperature limits
                    i0 = it*d0
                    minT = min(minT,buffer(i0+iz))
                    maxT = max(maxT,buffer(i0+iz))

                  end do
                end do
              end do

            ! Not static
            else

              ! For each column
              do iy=1,dims(2)

                ! Aborting
                if (aborting) exit

                do iz=1,dims(3)

                  ! Not cartesian
                  if (mode.eq.2) &
                    call get_point(unitA,dims(1),ycoor,zcoor,double)

                  ! Get data
                  call get_column(unitA,buffer,double,check)

                  ! Check could read
                  if (.not.check) then
                    aborting = .True.
                    exit
                  end if

                  ! For each point
                  do ix=1,dims(1)

                    ! Coordinate
                    i1 = di0*ix

                    ! Update temperature limits
                    i0 = it*d0
                    minT = min(minT,buffer(i0+iz))
                    maxT = max(maxT,buffer(i0+iz))

                    ! Update velocity limit
                    lv = 0d0
                    i0 = iv*d0
                    lv = lv + buffer(i0+iz)*buffer(i0+iz)
                    i0 = i0 + d0
                    lv = lv + buffer(i0+iz)*buffer(i0+iz)
                    i0 = i0 + d0
                    lv = lv + buffer(i0+iz)*buffer(i0+iz)
                    maxV = max(maxV,sqrt(lv))

                  end do
                end do
              end do

            end if ! Static
          end if ! 15D/CLE

        ! Dummy loop
        exit
        end do

        ! Close
        call close_file(unitA)

        ! Clean
        nullify(p_T,p_vx,p_vy,p_vz)

      end if ! Master

      ! Master share status
      call MPI_BCAST(aborting,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      ! Aborting?
      if (aborting) return

      ! Minimun temperature
      if (Input%minT.lt.0d0) then

        ! Master update input
        if (gpid.eq.0) Input%minT = minT

        ! Share
        call MPI_BCAST(Input%minT,1,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)

      end if

      ! Maximum temperature
      if (Input%maxT.lt.0d0) then

        ! Master update input
        if (gpid.eq.0) Input%maxT = maxT

        ! Share
        call MPI_BCAST(Input%maxT,1,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)

      end if

      ! If static
      if (Input%static) then

        ! Just 0
        Input%maxV = 0d0

      ! Not explicitly static
      else if (Input%maxV.lt.0d0) then

        ! Master update
        if (gpid.eq.0) Input%maxV = maxV

        ! Share
        call MPI_BCAST(Input%maxV,1,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)

      end if ! V

      return

      end subroutine get_lims

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read data about ionization in the given buffer after the
      !! specified offset\n
      !!  unit_index(integer(:)): Units to read from\n
      !!       buffer(double(:)): Buffer to store the data\n
      !!             nx(integer): Size of LOS axis\n
      !!         offset(integer): Positions already filled with the
      !!                          model atmosphere in buffer\n
      !!         double(logical): If data in double precision\n
      !!          check(logical): Success flag
      subroutine get_column_ion(unit_index,buffer,nx,offset, &
                                double,check)

      ! I/O

      logical, intent(in):: double
      logical, intent(out):: check
      integer, intent(in):: nx,offset
      integer, dimension(:), intent(in):: unit_index
      double precision, dimension(:), intent(out):: buffer

      ! Local

      integer:: ion,i0,i1

      real, dimension(:), allocatable:: bufferr


      ! Indexes to read into
      i0 = offset + 1
      i1 = offset + nx

      ! If single precision, allocate float buffer
      if (.not.double) allocate(bufferr(nx))

      ! For each ion
      do ion=1,size(unit_index)

        ! Double
        if (double) then

          ! Read directly into double
          read(unit_index(ion), err=1100) buffer(i0:i1)

        ! Single
        else

          ! Read into float buffer
          read(unit_index(ion), err=1100) bufferr

          ! Transform into double precision
          buffer(i0:i1) = dble(bufferr)

        end if ! Double precision

        ! Shift position in buffer
        i0 = i0 + nx
        i1 = i1 + nx

      end do ! Ions we need to read

      ! If single, deallocate float buffer
      if (.not.double) deallocate(bufferr)

      ! Success
      check = .True.

      return

1100  check = .False.
      umsg = ' # Error reading column from the ion file'
      call verbose

      return

      end subroutine get_column_ion

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the cache from an existing cache file\n
      !!     unit_index(integer): Unit to read from\n
      !!  filename(character(:)): Name of the file to open\n
      !!        dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!     cache(logical(:,:)): Cache of already done columns\n
      !!          check(logical): Success flag
      subroutine get_cache(unit_index, filename, dims, cache, check)

      ! I/O

      character(len=500), intent(in):: filename
      logical, intent(out):: check
      logical, dimension(:,:), intent(inout):: cache
      integer, intent(in):: unit_index
      integer, dimension(:), intent(in):: dims

      ! Local

      character(len=5):: label

      integer:: ios,buffer,ix,iy
      integer, dimension(2):: ldims


      ! Initialize success flag
      check = .False.

      ! Open file
      open (unit_index, file=trim(filename), status='unknown', &
            iostat=ios, err=1000, access='stream', action='read', &
            form='unformatted')

      ! Read label
      read(unit_index,err=1100,end=1100) label

      ! Wrong label
      if (label.ne.'cache') then

        ! Issue error
        umsg = ' # Cache file has wrong label'
        call verbose
        goto 1100

      ! Correct label
      else

        ! Verbose
        umsg = ' # Read cache file: '//trim(filename)
        call verbose

      end if

      ! Read dimensions in cache
      read(unit_index,err=1100) ldims

      ! Check dimensions
      if (ldims(1).ne.dims(1).or.ldims(2).ne.dims(2)) then

        ! Issue error
        umsg = ' # Cache file has wrong dimensions'
        call verbose
        goto 1100

      end if

      ! Success if got here
      check = .True.

      ! Read cache
      do ix=1,dims(1)
        do iy=1,dims(2)

          ! Read data and stop if failure to do so
          read(unit_index,err=1100,end=1100) buffer

          ! Transform to bool
          cache(iy,ix) = buffer.gt.0

        end do ! Y
      end do ! X

      ! Close unit and leave
1100  close(unit_index)
1000  return

      end subroutine get_cache

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the cache file\n
      !!     unit_index(integer): Unit where to open the file\n
      !!  filename(character(:)): Name of the file to open\n
      !!        dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!         lcache(logical): If cache exists\n
      !!          check(logical): Success flag
      subroutine start_cache(unit_index,filename,dims,lcache,check)

      ! I/O

      character(len=500), intent(in):: filename
      logical, intent(in):: lcache
      logical, intent(out):: check
      integer, intent(in):: unit_index
      integer, dimension(:), intent(in):: dims

      ! Local

      character(len=6), allocatable:: cpos

      integer:: ios,ic,NC


      ! Initialize success flag
      check = .False.

      ! Define cpos
      if (lcache) then

        ! If the cache exists, we are appending
        cpos = 'APPEND'

      ! If cache does not exist
      else

        ! We are initializing
        cpos = 'ASIS  '

      end if

      ! Open file
      open (unit_index, file=trim(filename), status='unknown', &
            iostat=ios, err=1000, access='stream', action='write', &
            form='unformatted',position=trim(cpos))

      ! If file does not exist
      if (.not.lcache) then

        ! Write header
        write(unit_index,err=1000) 'cache'
        write(unit_index,err=1000) dims(1:2)

        ! Initialize full file
        NC = dims(1)*dims(2)

        ! For all columns
        do ic=1,NC

          ! Write 0
          write(unit_index,err=1000) 0

        end do ! Columns

      end if ! File does not exist

      ! Success
      check = .True.

      ! Close unit
      close(unit_index)

      return

1000  umsg = ' # Error initializing cache file'
      call verbose
      return

      end subroutine start_cache

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write entry into the cache file\n
      !!     unit_index(integer): Unit to write to\n
      !!  filename(character(:)): Name of the file to open\n
      !!    register(integer(:)): Coordinates of column and result\n
      !!          cache(integer): Status of solution\n
      !!          check(logical): Success flag
      subroutine write_cache(unit_index,filename,register,check)

      ! I/O

      character(len=500), intent(in):: filename
      logical, intent(out):: check
      integer, intent(in):: unit_index
      integer, dimension(:), intent(in):: register

      ! Local

      integer:: ios,offset,whence


      ! Initialize success flag
      check = .False.

      ! Open file
      open (unit_index, file=trim(filename), status='unknown', &
            iostat=ios, err=1000, access='stream', action='write', &
            form='unformatted',position='APPEND')

      ! Seek position in file
      offset = 13 + (register(1)-1)*4
      whence = 0
      call fseek(unit_index, offset, whence, ios)

      ! Check error
      if (ios.ne.0) goto 1000

      ! Write result
      write(unit_index,err=1000) register(3)

      ! Success if here
      check = .True.

      ! Close
      close(unit_index)

      return

1000  umsg = ' # Error writing in cache file'
      return

      end subroutine write_cache

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate and setup the index limits in the wavelength
      !! dependent variables to cut the synthesis input based on the
      !! wavelengths in the inversion data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!          omega(double(:)): Wavelength axis in the observation
      subroutine prepare_lambda_limits(Input,Inf_Stokes,omega)

      ! I/O

      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      double precision, dimension(:), intent(in):: omega

      ! Local

      integer:: iran,jran


      ! Copy number of ranges in data
      Input%lim_stk%nran = Inf_Stokes%Num_Range

      ! Allocate Stokes ranges
      allocate(Input%lim_stk%doub(2,Input%lim_stk%nran))
      MRAMc = MRAMc + 1d-6*sizeof(Input%lim_stk%doub)

      !
      ! Set-up Stokes ranges
      !

      ! For each range
      do iran=1,Input%lim_stk%nran

        ! Ranges are usually reversed, so get correct index
        jran = Input%lim_stk%nran - iran + 1

        ! Get from input axis
        Input%lim_stk%doub(1,jran) = &
                                   1d2/omega(Inf_Stokes%Range(iran,2))
        Input%lim_stk%doub(2,jran) = &
                                   1d2/omega(Inf_Stokes%Range(iran,1))

      end do ! Ranges

      ! Output tau
      if (Input%out_tau1) then

        ! Copy number of ranges in data
        Input%lim_tau%nran = Inf_Stokes%Num_Range

        ! Allocate Tau ranges
        allocate(Input%lim_tau%doub(2,Input%lim_tau%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_tau%doub)

        !
        ! Set-up Stokes ranges
        !

        ! For each range
        do iran=1,Input%lim_tau%nran

          ! Ranges are usually reversed, so get correct index
          jran = Input%lim_stk%nran - iran + 1

          ! Get from input axis
          Input%lim_tau%doub(1,jran) = &
                                   1d2/omega(Inf_Stokes%Range(iran,2))
          Input%lim_tau%doub(2,jran) = &
                                   1d2/omega(Inf_Stokes%Range(iran,1))

        end do ! Ranges

      end if ! Output tau

      ! Output contribution function
      if (Input%out_contr) then

        ! Copy number of ranges in data
        Input%lim_ctr%nran = Inf_Stokes%Num_Range

        ! Allocate contribution function ranges
        allocate(Input%lim_ctr%doub(2,Input%lim_ctr%nran))
        MRAMc = MRAMc + 1d-6*sizeof(Input%lim_ctr%doub)

        !
        ! Set-up Stokes ranges
        !

        ! For each range
        do iran=1,Input%lim_ctr%nran

          ! Ranges are usually reversed, so get correct index
          jran = Input%lim_stk%nran - iran + 1

          ! Get from input axis
          Input%lim_ctr%doub(1,jran) = &
                                   1d2/omega(Inf_Stokes%Range(iran,2))
          Input%lim_ctr%doub(2,jran) = &
                                   1d2/omega(Inf_Stokes%Range(iran,1))

        end do ! Ranges

      end if ! Output contribution function

      end subroutine prepare_lambda_limits

!#####################################################################
!#####################################################################
!#####################################################################

      !> Specify the index limit to choose what to write in the output
      !! of a wavelength dependent variable\n
      !!  buff(IO_helper_class): Structure with IO limit data\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine set_lambda_limit(buff,Frec)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      type(IO_helper_class), intent(inout):: buff

      ! Local

      integer:: iran


      ! If limits specified
      if (buff%nran.gt.0) then

        ! Allocate indexes limits and range size
        allocate(buff%indx(2,buff%nran))
        MRAMc = MRAMc + 1d-6*sizeof(buff%indx)
        allocate(buff%nbuff(buff%nran))
        MRAMc = MRAMc + 1d-6*sizeof(buff%nbuff)

        ! Initialize size
        buff%nn = 0

        ! Look for indexes in each range
        do iran=1,buff%nran

          ! Check location of closest frequencies
          buff%indx(1,iran) = minloc(abs(Frec%omega - &
                                         buff%doub(1,iran)),1)
          buff%indx(2,iran) = minloc(abs(Frec%omega - &
                                         buff%doub(2,iran)),1)

          ! Get size of this range
          buff%nbuff(iran) = buff%indx(2,iran) - buff%indx(1,iran) + 1

          ! Get cumulative size
          buff%nn = buff%nn + buff%nbuff(iran)

        end do ! Ranges

        ! Deallocate doubles
        MRAMc = MRAMc - 1d-6*sizeof(buff%doub)
        deallocate(buff%doub)

      ! Not specified
      else

        ! Just write the full range
        buff%nn = nfreq

      end if ! Specified limits

      end subroutine set_lambda_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Specify the buffer size to write in the output for
      !! collisions\n
      !!  buff1(IO_helper_class): Structure with IO data\n
      !!  buff2(IO_helper_class): Structure with IO data\n
      !!        Atom(Atom_class): Structure with the atomic data
      subroutine set_cols_limit(buff1,buff2,Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff1
      type(IO_helper_class), intent(inout):: buff2

      ! Local

      integer:: ia


      !
      ! First set (term to term)
      !

      ! If specified
      if (buff1%nran.gt.0) then

        ! Get size
        buff1%nn = buff1%nran*2

      ! Not specified
      else

        ! Initialize sizes
        buff1%nn = 0

        ! For each atom
        do ia=1,na

          ! Add everything to size
          buff1%nn = buff1%nn + Atom(ia)%nMulti*Atom(ia)%nMulti

        end do ! Atoms

      end if ! Specified

      !
      ! Second set (level to level)
      !

      ! If specified
      if (buff2%nran.gt.0) then

        ! Get size
        buff2%nn = buff2%nran*2

      ! Not specified
      else

        ! Initialize sizes
        buff2%nn = 0

        ! For each atom
        do ia=1,na

          ! Add everything to size
          buff2%nn = buff2%nn + Atom(ia)%nlevel*Atom(ia)%nlevel

        end do ! Atoms

      end if ! Specified

      end subroutine set_cols_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Specify the buffer size to write in the output for damping
      !! parameters\n
      !!  buff(IO_helper_class): Structure with IO data\n
      !!       Atom(Atom_class): Structure with the atomic data
      subroutine set_damp_limit(buff,Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff

      ! Local

      integer:: ia


      ! If specified
      if (buff%nran.gt.0) then

        ! Get size
        buff%nn = buff%nran

      ! Not specified
      else

        ! Initialize sizes
        buff%nn = 0

        ! For each atom
        do ia=1,na

          ! Add everythingto size
          buff%nn = buff%nn + Atom(ia)%ntran

        end do ! Atoms

      end if ! Specified

      end subroutine set_damp_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Specify the buffer size to write in the output for
      !! populations and departure coefficients\n
      !!  buff(IO_helper_class): Structure with IO data\n
      !!       Atom(Atom_class): Structure with the atomic data
      subroutine set_pop_limit(buff,Atom)

      ! I/O
      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff

      ! Local
      integer:: ia,iran

      ! Allocate buffer per atom
      allocate(buff%nbuff(nA))
      MRAMc = MRAMc + 1d-6*sizeof(buff%nbuff)

      ! Initialize size
      buff%nbuff = 0

      ! If specified
      if (buff%nran.gt.0) then

        ! For each atom
        do ia=1,nA

          ! For each range
          do iran=1,buff%nran

            ! If atom coincide coincide, add
            if (ia.eq.buff%indx(1,iran)) &
              buff%nbuff(ia) = buff%nbuff(ia) + 1

          end do ! Range
        end do ! Atoms

      ! Not specified
      else

        ! For each atom
        do ia=1,na

          ! Add everything to size
          buff%nbuff(ia) = Atom(ia)%nlevel

        end do ! Atoms

      end if ! Specified

      end subroutine set_pop_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Setup buffer sizes and limits to write the output in 1.5D
      !! synthesis mode\n
      !!     Input(Input_class): Structure with configuration data\n
      !!          mode(integer): Type of atmospheric model\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine set_io_buffers(Input,mode,Atom,Frec)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Input_class), intent(inout):: Input
      integer, intent(in):: mode


      !
      ! Stokes
      !

      ! Look for wavelength indexes
      call set_lambda_limit(Input%lim_stk,Frec)

      ! Header and buffer sizes
      Input%lim_stk%head_size = 32 + Input%lim_stk%nn*8
      Input%lim_stk%buffer_size = 16*Input%lim_stk%nn


      !
      ! Tau
      !
      if (Input%out_tau1) then

        ! Look for wavelength indexes
        call set_lambda_limit(Input%lim_tau,Frec)

        ! Header and buffer sizes
        Input%lim_tau%head_size = 32 + Input%lim_tau%nn*8
        Input%lim_tau%buffer_size = 4*Input%lim_tau%nn

      end if

      ! If mode is not cartesian
      if (mode.ne.0) then

        ! Skip rest
        Input%out_contr = .False.
        Input%keep_cols = .False.
        Input%keep_damp = .False.
        Input%keep_qel = .False.
        Input%keep_back = .False.
        Input%keep_pop = .False.
        Input%keep_dep = .False.
        Input%keep_atmo = .False.
        return

      end if


      !
      ! Contribution function
      !
      if (Input%out_contr) then

        ! Look for wavelength indexes
        call set_lambda_limit(Input%lim_ctr,Frec)

        ! Header and buffer sizes
        Input%lim_ctr%head_size = 36 + Input%lim_ctr%nn*8
        Input%lim_ctr%buffer_size = 16*Input%lim_ctr%nn*nz

      end if


      !
      ! Collisions
      !
      if (Input%keep_cols) then

        ! Set collision limits
        call set_cols_limit(Input%lim_cols_tt,Input%lim_cols_ll,Atom)

        ! Header and buffer sizes
        Input%lim_cols_tt%head_size = 24
        Input%lim_cols_tt%buffer_size = Input%lim_cols_tt%nn*4*nz
        Input%lim_cols_ll%head_size = 24
        Input%lim_cols_ll%buffer_size = Input%lim_cols_ll%nn*4*nz

      end if


      !
      ! Damping parameter
      !
      if (Input%keep_damp) then

        ! Set damping parameter limits
        call set_damp_limit(Input%lim_damp,Atom)

        ! Header and buffer sizes
        Input%lim_damp%head_size = 24
        Input%lim_damp%buffer_size = Input%lim_damp%nn*4*nz

      end if


      !
      ! Elastic rates
      !
      if (Input%keep_qel) then

        ! Set elastic rates limits
        call set_damp_limit(Input%lim_qel,Atom)

        ! Header and buffer sizes
        Input%lim_qel%head_size = 24
        Input%lim_qel%buffer_size = Input%lim_qel%nn*4*nz

      end if


      !
      ! Background
      !
      if (Input%keep_back) then

        ! Look for wavelength indexes
        call set_lambda_limit(Input%lim_back,Frec)

        ! Header and buffer sizes
        Input%lim_back%head_size = 24 + Input%lim_back%nn*8
        Input%lim_back%buffer_size = 12*Input%lim_back%nn*nz

      end if


      !
      ! Populations
      !
      if (Input%keep_pop.or.Input%keep_dep) then

        ! Set population limits
        call set_pop_limit(Input%lim_pop,Atom)

        ! Header and buffer sizes
        Input%lim_pop%head_size = 20
        Input%lim_pop%nbuff = Input%lim_pop%nbuff*nz*4

      end if

      !
      ! Atmosphere
      !
      if (Input%keep_atmo) then
        Input%lim_atmo%head_size = 20
        Input%lim_atmo%buffer_size = 24*8*nz
      end if

      end subroutine set_io_buffers

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks that the specified limits for the collisional output
      !! comply with the atomic models\n
      !!  buff1(IO_helper_class): Structure with IO data\n
      !!  buff2(IO_helper_class): Structure with IO data\n
      !!     Atom(Atom_class(:)): Structures with atomic data
      subroutine check_cols_limit(buff1,buff2,Atom)

      ! I/O
      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff1
      type(IO_helper_class), intent(inout):: buff2

      ! Local

      integer:: ia,i1,i2,iran


      !
      ! Term-term collisions
      !

      ! If specified
      if (buff1%nran.gt.0) then

        ! For each range
        do iran=1,buff1%nran

          ! Get atom and terms
          ia = buff1%indx(1,iran)
          i1 = buff1%indx(2,iran)
          i2 = buff1%indx(3,iran)

          ! Check atom
          if (ia.gt.nA) then

            ! Error
            umsg = 'You have specified an atomic index to '//&
                   'output term to term collisions larger '//&
                   'than the number of atoms'
            urou = 'check_cols_limit'
            call aborted
            return

          end if

          ! Check terms
          if (i1.gt.Atom(ia)%nmulti.or.i2.gt.Atom(ia)%nmulti) then

            ! Error
            umsg = 'You have specified a term index to '//&
                   'output term to term collisions larger '//&
                   'than the number of terms in the atom'
            urou = 'check_cols_limit'
            call aborted
            return

          end if

        end do ! Atoms

      end if

      !
      ! Level-level collisions
      !

      ! If specified
      if (buff2%nran.gt.0) then

        ! For each range
        do iran=1,buff2%nran

          ! Get atom and levels
          ia = buff2%indx(1,iran)
          i1 = buff2%indx(2,iran)
          i2 = buff2%indx(3,iran)

          ! Check atom
          if (ia.gt.nA) then

            ! Error
            umsg = 'You have specified an atomic index to '//&
                   'output level to level collisions larger '//&
                   'than the number of atoms'
            urou = 'check_cols_limit'
            call aborted
            return

          end if

          ! Check level
          if (i1.gt.Atom(ia)%nlevel.or.i2.gt.Atom(ia)%nlevel) then

            ! Error
            umsg = 'You have specified a level index to '//&
                   'output level to level collisions larger '//&
                   'than the number of levels in the atom'
            urou = 'check_cols_limit'
            call aborted
            return

          end if

        end do ! Atoms

      end if ! Specified limits

      end subroutine check_cols_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks that the specified limits for the damping parameter
      !! output comply with the atomic models\n
      !!  buff(IO_helper_class): Structure with IO data\n
      !!    Atom(Atom_class(:)): Structures with atomic data
      subroutine check_damp_limit(buff,Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff

      ! Local

      integer:: ia,iran,it


      ! If specified
      if (buff%nran.gt.0) then

        ! For each range
        do iran=1,buff%nran

          ! Get atom and levels
          ia = buff%indx(1,iran)
          it = buff%indx(2,iran)

          ! Check atom
          if (ia.gt.nA) then

            ! Issue error
            umsg = 'You have specified an atomic index to '//&
                   'output damping parameter larger '//&
                   'than the number of atoms'
            urou = 'check_damp_limit'
            call aborted
            return

          end if

          ! Check terms
          if (it.gt.Atom(ia)%ntran) then

            ! Issue error
            umsg = 'You have specified a transition index to '//&
                   'output damping parameter larger '//&
                   'than the number of transitions in the atom'
            urou = 'check_damp_limit'
            call aborted
            return

          end if

        end do ! Atoms

      end if ! Specified limits

      end subroutine check_damp_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks that the specified limits for the elastic rates output
      !! comply with the atomic models\n
      !!  buff(IO_helper_class): Structure with IO data\n
      !!    Atom(Atom_class(:)): Structures with atomic data
      subroutine check_qel_limit(buff,Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff

      ! Local

      integer:: ia,iran,it


      ! If specified
      if (buff%nran.gt.0) then

        ! For each range
        do iran=1,buff%nran

          ! Get atom and levels
          ia = buff%indx(1,iran)
          it = buff%indx(2,iran)

          ! Check atom
          if (ia.gt.nA) then

            ! Issue error
            umsg = 'You have specified an atomic index to '//&
                   'output elastic rates larger '//&
                   'than the number of atoms'
            urou = 'check_qel_limit'
            call aborted
            return

          end if

          ! Check terms
          if (it.gt.Atom(ia)%ntran) then

            ! Issue error
            umsg = 'You have specified a transition index to '//&
                   'output elastic rates larger '//&
                   'than the number of transitions in the atom'
            urou = 'check_qel_limit'
            call aborted
            return

          end if

        end do ! Atoms

      end if ! Specified limits

      end subroutine check_qel_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks that the specified limits for the population output
      !! comply with the atomic models\n
      !!  buff(IO_helper_class): Structure with IO data\n
      !!    Atom(Atom_class(:)): Structures with atomic data
      subroutine check_pop_limit(buff,Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(inout):: buff

      ! Local

      integer:: iran,ia,il


      ! If specified
      if (buff%nran.gt.0) then

        ! For each range
        do iran=1,buff%nran

          ! Get atom and levels
          ia = buff%indx(1,iran)
          il = buff%indx(2,iran)

          ! Check atom
          if (ia.gt.nA) then

            ! Issue error
            umsg = 'You have specified an atomic index to '//&
                   'output populations larger '//&
                   'than the number of atoms'
            urou = 'check_pop_limit'
            call aborted
            return

          end if

          ! Check levels
          if (il.gt.Atom(ia)%nlevel) then

            ! Issue error
            umsg = 'You have specified a level index to '//&
                   'output populations larger '//&
                   'than the number of levels in the atom'
            urou = 'check_pop_limit'
            call aborted
            return

          end if

        end do ! Atoms

      end if ! Specified limits

      end subroutine check_pop_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the check_*_limit routines in sequence to check that
      !! the specified limits make sense for the actual atomic
      !! models\n
      !!     Input(Input_class): Structure with configuration data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine check_io_buffers_sanity_check(Input,Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Input_class), intent(inout):: Input


      ! Collisions
      if (Input%keep_cols) &
        call check_cols_limit(Input%lim_cols_tt, &
                              Input%lim_cols_ll,Atom)

      ! Damping parameter
      if (Input%keep_damp) call check_damp_limit(Input%lim_damp,Atom)

      ! Elastic rates
      if (Input%keep_qel) call check_qel_limit(Input%lim_qel,Atom)

      ! Populations
      if (Input%keep_pop.or.Input%keep_dep) &
        call check_pop_limit(Input%lim_pop,Atom)

      end subroutine check_io_buffers_sanity_check

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if the output files for the 1.5D synthesis already
      !! exist\n
      !!    Input(Input_class): Structure with configuration data\n
      !!   Atom(Atom_class(:)): Structures with the atomic data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!     aborting(logical): If something goes wrong
      subroutine check_io_buffers_exist(Input,Atom,Geom,aborting)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Input_class), intent(inout):: Input
      logical, intent(out):: aborting

      ! Local

      character(len=4):: cth,cph
      integer:: ith,iph,ios,ia


      ! Initialize
      aborting = .False.

      !
      ! Stokes, Contribution, and Tau
      !

      ! For each LOS polar direction
      do ith=1,Geom%nThLOS

        ! Convert the integers into appropriate length strings
        if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
        if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
        if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith

        ! For each LOS azimuth
        do iph=1,Geom%nPhLOS

          ! Convert the integers into appropriate length strings
          if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
          if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
          if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

          ! If intensity calc.
          if (Input%force.eq.'I') then

            ! Try opening Stokes
            open (200,file=trim(Input%folder)//'/StokesI_'//&
                  trim(cth)//'_'//trim(cph), status='old', &
                  iostat=ios, access='stream', action='read', &
                  form='unformatted')

          ! If polarization calc.
          else

            ! Try opening Stokes
            open (200,file=trim(Input%folder)//'/Stokes_'//&
                  trim(cth)//'_'//trim(cph), status='old', &
                  iostat=ios, access='stream', action='read', &
                  form='unformatted')

          end if ! Intensity or polarization

          ! If could not open
          if (ios.ne.0) then

            ! Issue error
            umsg = 'There is no existing Stokes files'
            call verbose
            aborting = .True.
            return

          end if

          ! Close
          close(200)

          ! If contribution function output
          if (Input%out_contr) then

            ! Try opening Contribution function
            open (200,file=trim(Input%folder)//'/Contribution_'//&
                  trim(cth)//'_'//trim(cph), status='old', &
                  iostat=ios, access='stream', action='read', &
                  form='unformatted')

            ! If could not open
            if (ios.ne.0) then

              ! Issue error
              umsg = 'There is no existing contribution '// &
                     'function files'
              call verbose
              aborting = .True.
              return

            end if

            ! Close
            close(200)

          end if ! Contribution function output

          ! If tau output
          if (Input%out_tau1) then

            ! Try opening tau
            open (200,file=trim(Input%folder)//'/Tau_'//&
                  trim(cth)//'_'//trim(cph), status='old', &
                  iostat=ios, access='stream', action='read', &
                  form='unformatted')

            ! If could not open
            if (ios.ne.0) then

              ! Issue error
              umsg = 'There is no existing tau_1 files'
              call verbose
              aborting = .True.
              return

            end if

            ! Close
            close(200)

          end if ! Tau output

        end do ! Azimuth (LOS)
      end do ! Polar (LOS)

      ! If storing collisions
      if (Input%keep_cols) then

        ! Try opening collisions file term-term
        open (200,file=trim(Input%folder)//'/cols-TT', status='old', &
              iostat=ios, access='stream', action='read', &
              form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing cols-TT file'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

        ! Try opening collisions file level-level
        open (200,file=trim(Input%folder)//'/cols-LL', status='old', &
              iostat=ios, access='stream', action='read', &
              form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing cols-LL file'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

      end if ! Storing collisions

      ! If storing damping parameter
      if (Input%keep_damp) then

        ! Try opening damping parameter file
        open (200,file=trim(Input%folder)//'/damping', status='old', &
              iostat=ios, access='stream', action='read', &
              form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing damping file'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

      end if ! Storing damping parameter

      ! If storing elastic rates
      if (Input%keep_qel) then

        ! Try opening elastic rates file
        open (200,file=trim(Input%folder)//'/qel', status='old', &
              iostat=ios, access='stream', action='read', &
              form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing elastic rates file'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

      end if ! Storing elastic rates

      ! If storing background
      if (Input%keep_back) then

        ! Try opening background file
        open (200,file=trim(Input%folder)//'/background', &
              status='old', iostat=ios, access='stream', &
              action='read', form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing background file'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

      end if ! Storing background quantities

      ! If storing populations or departure coeff.
      if (Input%keep_pop.or.Input%keep_dep) then

        ! For each atom
        do ia=1,nA

          ! If buffer does not have size, skip
          if (Input%lim_pop%nbuff(ia).lt.1) cycle

          ! If storing pop
          if (Input%keep_pop) then

            ! Try opening population file
            open (200,file=trim(Input%folder)//'/'//&
                           trim(Atom(ia)%file_label)//'.pop', &
                  status='old', iostat=ios, access='stream', &
                  action='read', form='unformatted')

            ! If could not open
            if (ios.ne.0) then

              ! Issue error
              umsg = 'There is no existing population file'
              call verbose
              aborting = .True.
              return

            end if

            ! Close
            close(200)

          end if ! Storing population

          ! If storing dep
          if (Input%keep_dep) then

            ! Try opening departure file
            open (200,file=trim(Input%folder)//'/'//&
                             trim(Atom(ia)%file_label)//'.dep', &
                  status='old', iostat=ios, access='stream', &
                  action='read', form='unformatted')

            ! If could not open
            if (ios.ne.0) then

              ! Issue error
              umsg = 'There is no existing departura coefficient '//&
                     'file'
              call verbose
              aborting = .True.
              return

            end if

            ! Close
            close(200)

          end if ! Storing departure coefficient

        end do ! Atoms

      end if ! Storing population or departure coefficient

      ! If storing atmosphere
      if (Input%keep_atmo) then

        ! Try opening atmosphetic file
        open (200,file=trim(Input%folder)//'/atmo.hrt', &
              status='old', iostat=ios, access='stream', &
              action='read', form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing atmosphere file'
          call verbose
          aborting = .True.
          return
        end if

        ! Close
        close(200)

      end if ! Storing atmosphere

      ! If storing MRC
      if (Input%keep_MRC) then

        ! Try opening MRC file
        open (200,file=trim(Input%folder)//'/MRC', &
              status='old', iostat=ios, access='stream', &
              action='read', form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing MRC file'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

      end if ! Storing atmosphere

      end subroutine check_io_buffers_exist

!#####################################################################
!#####################################################################
!#####################################################################

      !> Create the output files for the 1.5D synthesis and write
      !! their headers\n
      !!     Input(Input_class): Structure with configuration data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!  GeomI(Geometry_class): Structure with geometric data for the
      !!                         intensity problem\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine create_io_files(Input,Atom,dims,Geom,GeomI,Frec)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom, GeomI
      type(Input_class), intent(inout):: Input
      integer, dimension(:), intent(in):: dims

      ! Local

      character(len=4):: cth,cph

      integer:: ith,iph,ios,ia,iran,nThLOS,nPhLOS

      double precision, dimension(:), allocatable:: Tlos, PLos


      ! If intensity only
      if (Input%force.eq.'I') then

        ! Get number of LOS directions from GeomI
        nThLOS = GeomI%nThLOS
        nPhLOS = GeomI%nPhLOS

        ! If more than zero
        if (nThLOS.gt.0) then

          ! Allocate space and copy
          allocate(Tlos(nThLOS),Plos(nPhLOS))
          Tlos = GeomI%L_theta
          Plos = GeomI%L_phi

        end if

      ! If polarization
      else

        ! Get from Geom
        nThLOS = Geom%nThLOS
        nPhLOS = Geom%nPhLOS

        ! If more than zero
        if (nThLOS.gt.0) then

          ! Allocate space and copy
          allocate(Tlos(nThLOS),Plos(nPhLOS))
          Tlos = Geom%L_theta
          Plos = Geom%L_phi

        end if

      end if ! intensity or polarization


      !
      ! Stokes, Contribution, Tau
      !

      ! For each polar LOS direction
      do ith=1,nThLOS

        ! Convert the integers into appropriate length strings
        if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
        if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
        if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith

        ! For each azimuths LOS direction
        do iph=1,nPhLOS

          ! Convert the integers into appropriate length strings
          if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
          if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
          if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

          ! Intensity calc.
          if (Input%force.eq.'I') then

            ! Open intensity file
            open (200,file=trim(Input%folder)//'/StokesI_'//&
                  trim(cth)//'_'//trim(cph), status='unknown', &
                  iostat=ios, access='stream', action='write', &
                  form='unformatted')

          ! Polarization calc.
          else

            ! Open Stokes
            open (200,file=trim(Input%folder)//'/Stokes_'//&
                  trim(cth)//'_'//trim(cph), status='unknown', &
                  iostat=ios, access='stream', action='write', &
                  form='unformatted')

          end if ! Intensity/polarization

          ! Write header
          write(200) '2Dbe'
          write(200) Input%lim_stk%nn
          write(200) dims(1:2)
          write(200) Tlos(ith)*180d0/PI
          write(200) Plos(iph)*180d0/PI

          ! If specified ranges
          if (Input%lim_stk%nran.gt.0) then

            ! For each range
            do iran=1,Input%lim_stk%nran

              ! Write relevant frequencies
              write(200) Frec%omega(Input%lim_stk%indx(1,iran): &
                                    Input%lim_stk%indx(2,iran))

            end do ! Ranges

          ! Full range
          else

            ! Write frequency
            write(200) Frec%omega

          end if ! Specified ranges

          ! Close
          close(200)

          ! If contribution function output
          if (Input%out_contr) then

            ! Open Contribution function
            open (200,file=trim(Input%folder)//'/Contribution_'//&
                  trim(cth)//'_'//trim(cph), status='unknown', &
                  iostat=ios, access='stream', action='write', &
                  form='unformatted')

            ! Write header
            write(200) '2Dbc'
            write(200) Input%lim_ctr%nn
            write(200) dims
            write(200) Tlos(ith)*180d0/PI
            write(200) Plos(iph)*180d0/PI

            ! If specified ranges
            if (Input%lim_ctr%nran.gt.0) then

              ! For each range
              do iran=1,Input%lim_ctr%nran

                ! Write relevant frequencies
                write(200) Frec%omega(Input%lim_ctr%indx(1,iran): &
                                      Input%lim_ctr%indx(2,iran))
              end do ! Ranges

            ! Full range
            else

              ! Write frequency
              write(200) Frec%omega

            end if ! Specified ranges

            ! Close
            close(200)

          end if

          ! If tau output
          if (Input%out_tau1) then

            ! Open tau
            open (200,file=trim(Input%folder)//'/Tau_'//&
                  trim(cth)//'_'//trim(cph), status='unknown', &
                  iostat=ios, access='stream', action='write', &
                  form='unformatted')

            ! Write header
            write(200) '2Dbt'
            write(200) Input%lim_tau%nn
            write(200) dims(1:2)
            write(200) Tlos(ith)*180d0/PI
            write(200) Plos(iph)*180d0/PI

            ! If specified ranges
            if (Input%lim_tau%nran.gt.0) then

              ! For each range
              do iran=1,Input%lim_tau%nran

                ! Write relevant frequencies
                write(200) Frec%omega(Input%lim_tau%indx(1,iran): &
                                      Input%lim_tau%indx(2,iran))
              end do ! Ranges

            ! Full range
            else

              ! Write frequency
              write(200) Frec%omega

            end if ! Specified ranges

            ! Close
            close(200)

          end if ! Output tau

        end do ! Azimuth (LOS)
      end do ! Polar (LOS)


      !
      ! Collisions
      !
      if (Input%keep_cols) then

        ! Open collisions file term-term
        open (200,file=trim(Input%folder)//'/cols-TT', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dct'
        write(200) na
        write(200) dims
        write(200) Input%lim_cols_tt%nn

        ! Close
        close(200)

        ! Open collisions file level-level
        open (200,file=trim(Input%folder)//'/cols-LL', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dcl'
        write(200) na
        write(200) dims
        write(200) Input%lim_cols_ll%nn

        ! Close
        close(200)

      end if ! Output collisions


      !
      ! Damping parameter
      !
      if (Input%keep_damp) then

        ! Open damping parameter file
        open (200,file=trim(Input%folder)//'/damping', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dda'
        write(200) na
        write(200) dims
        write(200) Input%lim_damp%nn

        ! Close
        close(200)

      end if ! Output damping parameter


      !
      ! Elastic rates
      !
      if (Input%keep_qel) then

        ! Open elastic rates file
        open (200,file=trim(Input%folder)//'/qel', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dqe'
        write(200) na
        write(200) dims
        write(200) Input%lim_qel%nn

        ! Close
        close(200)

      end if ! Output elastic rates


      !
      ! Background
      !
      if (Input%keep_back) then

        ! Try opening background file
        open (200,file=trim(Input%folder)//'/background', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dba'
        write(200) dims
        write(200) Input%lim_back%nn
        write(200) 1

        ! If specified ranges
        if (Input%lim_back%nran.gt.0) then

          ! For each range
          do iran=1,Input%lim_back%nran

            ! Write relevant frequencies
            write(200) Frec%omega(Input%lim_back%indx(1,iran): &
                                  Input%lim_back%indx(2,iran))
          end do ! Ranges

        ! Full range
        else

          ! Write frequency
          write(200) Frec%omega

        end if ! Specified ranges

        ! Close
        close(200)

      end if ! Output background quantities


      !
      ! Populations
      !
      if (Input%keep_pop) then

        ! Atoms
        do ia=1,na

          ! If no output, skip
          if (Input%lim_pop%nbuff(ia).lt.1) cycle

          ! Open file
          open (200,file=trim(Input%folder)//'/'//&
                         trim(Atom(ia)%file_label)//'.pop', &
                status='unknown', iostat=ios, access='stream', &
                action='write', form='unformatted')

          ! Write header
          write(200) '2Dbp'
          write(200) dims
          write(200) Input%lim_pop%nbuff(ia)/4/dims(3)

          ! Close
          close(200)

        end do ! Atoms

      end if ! Output populations


      !
      ! Departure
      !
      if (Input%keep_dep) then

        ! Atoms
        do ia=1,na

          ! If no output, skip
          if (Input%lim_pop%nbuff(ia).lt.1) cycle

          ! Open file
          open (200,file=trim(Input%folder)//'/'//&
                         trim(Atom(ia)%file_label)//'.dep', &
                status='unknown', iostat=ios, access='stream', &
                action='write', form='unformatted')

          ! Write header
          write(200) '2Dbb'
          write(200) dims
          write(200) Input%lim_pop%nbuff(ia)/4/dims(3)

          ! Close
          close(200)

        end do ! Atoms

      end if ! Output departure coefficients


      !
      ! Atmosphere
      !
      if (Input%keep_atmo) then

        ! Open file
        open (200,file=trim(Input%folder)//'/atmo.hrt', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dat'
        write(200) int(8)
        write(200) dims

        ! Close
        close(200)

      end if ! Output atmosphere


      !
      ! MRC
      !
      if (Input%keep_MRC) then

        ! Open file
        open (200,file=trim(Input%folder)//'/MRC', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) 'MRC'
        write(200) dims(1:2)

        ! Close
        close(200)

      end if ! Output MRC

      ! Free
      if (allocated(Tlos)) deallocate(Tlos,Plos)

      end subroutine create_io_files

!#####################################################################
!#####################################################################
!#####################################################################

      !> Setup buffer sizes and limits to write the output in CLE
      !! synthesis mode\n
      !!     Input(Input_class): Structure with configuration data\n
      !!          mode(integer): Type of atmospheric model\n
      !!       dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine set_io_CLE_buffers(Input,mode,dims,Frec)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      type(Input_class), intent(inout):: Input
      integer, intent(in):: mode
      integer, dimension(3), intent(in):: dims


      !
      ! Stokes
      !

      ! Look for wavelength indexes
      call set_lambda_limit(Input%lim_stk,Frec)

      ! Header and buffer sizes
      Input%lim_stk%buffer_size = 16*Input%lim_stk%nn

      !
      ! Set geometry part size
      !

      ! Cartesian
      if (mode.eq.0) then

        ! Header and column size
        Input%lim_stk%head_size = 20 + Input%lim_stk%nn*8
        Input%lim_stk%geom_size = (dims(2) + dims(3))*8

      ! Slab or non-cartesian
      else if (mode.eq.1.or.mode.eq.2) then

        ! Header and column size
        Input%lim_stk%head_size = 16 + Input%lim_stk%nn*8
        Input%lim_stk%geom_size = dims(2)*16

      end if

      end subroutine set_io_CLE_buffers

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if the output files for the CLE synthesis already
      !! exist\n
      !!    Input(Input_class): Structure with configuration data\n
      !!     aborting(logical): Signals if something goes wrong
      subroutine check_io_CLE_buffers_exist(Input,aborting)

      ! I/O

      type(Input_class), intent(inout):: Input
      logical, intent(out):: aborting

      ! Local

      integer:: ios


      ! Initialize
      aborting = .False.

      !
      ! Stokes
      !

      ! Try opening Stokes
      open (200,file=trim(Input%folder)//'/Stokes', &
            status='old', iostat=ios, access='stream', &
            action='read', form='unformatted')

      ! If could not open
      if (ios.ne.0) then

        ! Issue error
        umsg = 'There is no existing Stokes files'
        call verbose
        aborting = .True.
        return

      end if

      ! Close
      close(200)

      !
      ! Tau output
      !
      if (Input%out_tau1) then

        ! Try opening tau
        open (200,file=trim(Input%folder)//'/Tau', &
              status='old', iostat=ios, access='stream', &
              action='read', form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing tau_1 files'
          call verbose
          aborting = .True.
          return

        end if

        ! Close
        close(200)

      end if

      end subroutine check_io_CLE_buffers_exist

!#####################################################################
!#####################################################################
!#####################################################################

      !> Create the output files for the CLE synthesis and write their
      !! headers\n
      !!     Input(Input_class): Structure with configuration data\n
      !!          mode(integer): Type of atmospheric model\n
      !!           y(double(:)): Y axis\n
      !!           z(double(:)): Z axis\n
      !!       dims(integer(:)): Grid dimensions (X,Y,Z)\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine create_io_CLE_files(Input,mode,y,z,dims,Frec)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      type(Input_class), intent(inout):: Input
      integer, intent(in):: mode
      integer, dimension(:), intent(in):: dims
      double precision, dimension(:), intent(in):: y,z

      ! Local

      integer:: ios,iran


      !
      ! Stokes
      !

      ! Open Stokes
      open (200,file=trim(Input%folder)//'/Stokes', &
            status='unknown', iostat=ios, access='stream', &
            action='write', form='unformatted')

      ! Write header
      write(200) 'CLEe'
      write(200) mode
      write(200) Input%lim_stk%nn

      ! If specified ranges
      if (Input%lim_stk%nran.gt.0) then

        ! For each range
        do iran=1,Input%lim_stk%nran

          ! Write relevant frequency
          write(200) Frec%omega(Input%lim_stk%indx(1,iran): &
                                Input%lim_stk%indx(2,iran))

        end do ! Ranges

      ! No specified ranges
      else

        ! Write whole axis
        write(200) Frec%omega

      end if ! Specified ranges

      ! Cartesian
      if (mode.eq.0) then

        ! Write Y
        write(200) dims(2)
        write(200) y

        ! Write Z
        write(200) dims(3)
        write(200) z

      ! Slab or non-cartesian
      else if (mode.eq.1.or.mode.eq.2) then

        ! Write dim
        write(200) dims(2)

      end if ! Type of model

      ! Close
      close(200)

      !
      ! Tau output
      !
      if (Input%out_tau1) then

        ! Open tau
        open (200,file=trim(Input%folder)//'/Tau', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) 'CLEc'
        write(200) mode
        write(200) Input%lim_stk%nn

        ! If specified range
        if (Input%lim_stk%nran.gt.0) then

          ! For each range
          do iran=1,Input%lim_stk%nran

            ! Write relevant frequencies
            write(200) Frec%omega(Input%lim_stk%indx(1,iran): &
                                  Input%lim_stk%indx(2,iran))

          end do ! Ranges

        ! No range specified
        else

          ! Write full axis
          write(200) Frec%omega

        end if ! Specified ranges

        ! Cartesian
        if (mode.eq.0) then

          ! Write Y
          write(200) dims(2)
          write(200) y

          ! Write Z
          write(200) dims(3)
          write(200) z

        ! Slab or non-cartesian
        else if (mode.eq.1.or.mode.eq.2) then

          ! Write dim
          write(200) dims(2)

        end if ! Type of model

        ! Close
        close(200)

      end if ! Output tau

      end subroutine create_io_CLE_files

!#####################################################################
!#####################################################################
!#####################################################################

      end module io_mod
