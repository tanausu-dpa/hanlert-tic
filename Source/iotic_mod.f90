      !> IO for TIC (no fits)
      module iotic_mod
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
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  name_check
!    Determine if a file is in fits format
!
!  Verbose_Model
!    Verbose the current values in the nodes
!
!  open_data_and_cache
!    Open the data file and the cache file, if present, and return the
!  dimensions, the metadata, and the success of the task
!
!  open_mask
!    Open the mask file and check dimensions
!
!  get_data_wavelength
!    Read the wavelength axis from the data file
!
!  get_data_los
!    Read the LOS information, constant for the whole FoV, from the
!  data file
!
!  get_data_sigma
!    Read the sigma (error) information, constant for the whole FoV,
!  from the data file
!
!  get_data_diff
!    Read the diffuse light profile, constant for the whole FoV, from
!  the data file
!
!  get_data_column
!    Read one pixel in the data file
!
!  set_up_data_frombuffer
!    Transform data in buffer into manageable variables
!
!  set_up_atmo_frombuffer
!    Transform atmosphere in buffer into the structure with
!  atmospheric data
!
!  set_up_JKQ_frombuffer
!    Transform JKQ in buffer into the structure with atmospheric data
!
!  set_inv_io_buffers
!    Calculate and store the sizes of the data to navigate the output
!  files
!
!  check_io_inv_buffer_exists
!    Check if the output files already exist
!
!  create_io_inv_files
!    Open and write the header for the output files
!
!  Write_Result
!    Write the results of the inversion of one pixel into a file
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use commons_mod
      use io_mod
      use model_mod
      use parameters_mod , only : c , TINYB
      use types_mod
      use w2freq_mod

      contains

#ifdef FITSSUP

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine if a file is in fits format\n
      !!  filename(character(:)): Path to the file\n
      !!     fits_index(integer): Indicate if it is a fits file
      subroutine name_check(filename,fits_index)

      ! I/O

      character(len=500), intent(in):: filename
      integer, intent(out):: fits_index

      ! Local

      logical:: anynull

      integer:: unitt,stat,rwmode,blocksize,group,bitpix
      integer:: fpixel,hdutype,hdunum,num,naxis
      integer, dimension(4):: naxes

      double precision:: nullval
      double precision, dimension(:), allocatable:: dummy

      !
      ! Look for the string with the format in the name
      ! and return if any of them is found, trusting the
      ! user naming convenction
      !
      fits_index = index(filename, '.fits', .True.)
      if (fits_index.gt.0) return
      fits_index = index(filename, '.FITS', .True.)
      if (fits_index.gt.0) return
      fits_index = index(filename, '.fts', .True.)
      if (fits_index.gt.0) return
      fits_index = index(filename, '.FTS', .True.)
      if (fits_index.gt.0) return

      ! Initialize assuming it is not a fits file
      fits_index = 0

      ! Only the master tries this
      if (pid.eq.0) then

        ! Initialize variables to read fits
        rwmode = 0
        group = 1
        fpixel = 1
        blocksize = 1
        nullval = -1d0
        naxes = 0
        stat = 0

        ! Do until done, the code aborts at the first error that
        ! is found when trying to read the fits
        do while (.True.)

          !
          ! Give input unit a number
          !
          unitt = 19

          ! Open file as fits
          call FTOPEN(unitt,trim(filename),rwmode,blocksize,stat)
          if (stat.ne.0) exit

          ! Read the total number of HDUS in the fits file
          call FTTHDU(unitt, hdunum, stat)
          if (stat.ne.0) exit
          if (hdunum.lt.3) exit

          ! Get the data type
          ! 8, 16, 32, 64, -32, or -64 corresponding to unsigned byte,
          ! signed 2-byte integer, signed 4-byte integer, signed
          ! 8-byte integer, real, and double.
          call FTGIDT(unitt,  bitpix, stat)
          if (stat.ne.0) exit
          if (bitpix.ne.-64) exit

          ! Get the dimension of the image
          call FTGIDM(unitt, naxis, stat)
          if (stat.ne.0) exit
          if (naxis.gt.4) exit

          !Move to hdu1, the wavelength dhu
          call FTMRHD(unitt, 1, hdutype, stat)
          if (stat.ne.0) exit

          ! Get the data type
          call FTGIDT(unitt,  bitpix, stat)
          if (stat.ne.0) exit
          if (bitpix.ne.-64) exit

          ! Get the dimension of the image
          call FTGIDM(unitt, naxis, stat)
          if (stat.ne.0) exit
          if (naxis.ne.1) exit

          ! Get the size of all dimensions of the image
          call FTGISZ(unitt, naxis, naxes, stat)
          if (stat.ne.0) exit
          Num = naxes(1)
          allocate(dummy(naxes(1)))

          ! Read lambda
          call FTGPVD(unitt, group, fpixel, naxes(1), nullval, &
                      dummy(1), anynull, stat)
          deallocate(dummy)

          ! If we reached this, we can kind of assume it is a fits
          ! file, but the name does not have the extension
          fits_index = len(trim(filename)) + 1

          ! Stop checking
          exit

        end do ! Do till done

        ! Close the file
        call FTCLOS(unitt, stat)

      end if ! Master

      ! Share index with everyone
      CALL MPI_BCAST(fits_index, 1, MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, ierr)

      return

      end subroutine name_check
#endif

!#####################################################################
!#####################################################################
!#####################################################################

      !> Verbose the current values in the nodes\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !           trial(logical): If called from a trial solution
      subroutine Verbose_Model(Inf_Nodes,trial)

      ! I/O

      type(Nodes_class), intent(in):: Inf_Nodes
      logical, intent(in):: trial

      ! Local

      character(3):: length
      character(30):: fmt

      integer:: i, j


      !
      ! Head verbosity of the model
      !

      ! If called from trial
      if (trial) then

        ! Verbose
        umsg = ' - Node parameters for trial'
        call verboseI(3)

      ! If not called from trial
      else

        ! Verbose
        umsg = ' - Current node parameters'
        call verboseI(3)

      end if ! From trial or not

      ! For each variable
      do i=1,Inf_Nodes%nvar

        ! If inverting variable
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! Write into message
          write(umsg, '(A,i2)') "   Parameter index = ",i
          call verboseI(3)

          ! Write number of nodes in string and prepare format
          write(length, "(i3)") Inf_Nodes%Num_Nodes(i)
          fmt = '   ('//trim(adjustl(length))//'es15.5)'
          fmt = trim(adjustl(fmt))

          ! If vertical or microturbulent velocity
          if (i.eq.Inf_Nodes%index_vx.or. &
              (i.eq.Inf_Nodes%index_vy.and. &
               Inf_Nodes%vtype.eq.0).or. &
              i.eq.Inf_Nodes%index_vz.or. &
              i.eq.Inf_Nodes%index_vm) then

            ! Write transformed velocity in verbosity file
            write(umsg, FMT=fmt) &
              (Inf_Nodes%Node(i)%Var(j)*c*1d6, j = 1, &
                                               Inf_Nodes%Num_Nodes(i))
            call verboseI(3)

          ! Other variables
          else

            ! Write to verbosity file
            write(umsg, FMT=fmt) &
              (Inf_Nodes%Node(i)%Var(j), j = 1, &
                                         Inf_Nodes%Num_Nodes(i))
            call verboseI(3)

          end if ! If velocity or not
        end if ! If inverting

      end do ! Variables

      return

      end subroutine Verbose_Model

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open the data file and the cache file, if present, and return
      !! the dimensions, the metadata, and the success of the task\n
      !!   Input(Input_class): Structure with configuration data\n
      !!       unitD(integer): Unit where to open open data\n
      !!       unitC(integer): Unit where to open cache\n
      !!    aborting(logical): Indicate failure at output\n
      !!     dims(integer(:)): Grid dimensions (X,Y,L)\n
      !!    finfo(integer(:)): Data file metadata\n
      !!  cache(logical(:,:)): Cache of already done columns\n
      !!      lcache(logical): If there is a cache at output
      subroutine open_data_and_cache(Input,unitD,unitC, &
                                     aborting,dims,finfo, &
                                     cache,lcache)

      ! I/O

      type(Input_class), intent(in):: Input
      logical, intent(out):: aborting,lcache
      logical, dimension(:,:), allocatable, intent(out):: cache
      integer, intent(out):: unitD,unitC
      integer, dimension(:), intent(out):: dims,finfo

      ! Local

      logical:: check

#ifdef FITSSUP
      integer:: stat,rwmode,blocksize,naxis,hdunum,hdutype,bitpix
      integer, dimension(4):: naxes
#endif


      ! Give cache unit a number
      unitC = 17

#ifdef FITSSUP

      !
      ! If FITS file
      !
      if (Input%FITSFILE) then

        !
        ! Give input unit a number
        !
        unitD = 19
        naxes = 0
        blocksize = 1
        rwmode = 0
        stat = 0

        ! Open file
        call FTOPEN(unitD,Input%Filename_Ob,rwmode,blocksize,stat)

        ! Read the total number of hdus in the fits file
        call FTTHDU(unitD,hdunum,stat)

        ! Get the data type
        call FTGIDT(unitD,bitpix,stat)

        ! Get the dimension of the image
        call FTGIDM(unitD,naxis,stat)

        ! Get the size of all dimensions of the image
        call FTGISZ(unitD,naxis,naxes,stat)

        dims(3) = naxes(1)

        ! Check the dimension
        if (naxis.eq.3) then

          ! Only intensity
          finfo(1) = 0
          dims(1) = naxes(3)
          dims(2) = naxes(2)
          dims(3) = naxes(1)

        else if (naxis.eq.4) then

          ! 4 Stokes
          finfo(1) = 1
          dims(1) = naxes(4)
          dims(2) = naxes(3)
          dims(3) = naxes(1)
          if (naxes(2).ne.4) goto 1000

        else

          ! If the dimension of the data is smaller than 3,
          ! there is a problem
          goto 1000

        end if

        ! Move to hdu1, the wavelength dhu
        call FTMRHD(unitD,1,hdutype,stat)

        ! Get the dimension of the image
        call FTGIDM(unitD,naxis,stat)

        ! Get the size of all dimensions of the image
        call FTGISZ(unitD,naxis,naxes,stat)

        ! If multi-dimension, leave
        if (naxis.gt.1) goto 1000

        ! problem with the size of Wavelength
        if (dims(3).ne.naxes(1)) goto 1000

        ! Move to hdu2 the mu dhu
        call FTMRHD(unitD,1,hdutype,stat)

        ! Get the dimension of the mu
        call FTGIDM(unitD,naxis,stat)

        ! Get the size of all dimensions
        call FTGISZ(unitD,naxis,naxes,stat)

        ! Check the mu size
        select case(naxis)
          case(1)
            if (naxes(1).eq.2) then
              finfo(2) = 0
            else
              goto 1000
            end if

          case(3)
            if (naxes(1).eq.2.and.naxes(2).eq.1 &
                .and.naxes(3).eq.1) then
              finfo(2) = 0
            else if (naxes(1).eq.2.and.naxes(2).eq.dims(2) &
                .and.naxes(3).eq.dims(1)) then
              finfo(2) = 1
            else
              goto 1000
            end if

          case default
            goto 1000

        end select ! Type of node

        ! if more than 3 hdu
        if (hdunum.gt.3) then

          ! Move to hdu3 the sigma dhu
          call FTMRHD(unitD,1,hdutype,stat)

          ! Get the dimension of the sigma
          call FTGIDM(unitD,naxis,stat)

          if (naxis.gt.0) then

            ! Get the size of all dimensions
            call FTGISZ(unitD,naxis,naxes,stat)

            ! Get the date information
            ! if only intensity
            if (finfo(1).eq.0) then

              ! Check the size
              select case(naxis)
                case(1)
                  !data(1)
                  if (naxes(1).eq.1) then
                    finfo(3) = 1
                  !data(nl)
                  else if (naxes(1).eq.dims(3)) then
                    finfo(3) = 2
                  else
                    finfo(3) = 0
                  end if

                case(2)
                  !data(1,1)
                  if (naxes(1).eq.1.and.naxes(2).eq.1) then
                    finfo(3) = 1
                  !data(nl,1)
                  else if (naxes(1).eq.dims(3) &
                      .and.naxes(2).eq.1) then
                    finfo(3) = 2
                  !data(ny,nx)
                  else if (naxes(1).eq.dims(2) &
                      .and.naxes(2).eq.dims(1)) then
                    finfo(3) = 3
                  else
                    finfo(3) = 0
                  end if

                case(3)
                  !data(1,ny,nx)
                  if (naxes(1).eq.1.and.naxes(2).eq.dims(2) &
                      .and.naxes(3).eq.dims(1)) then
                    finfo(3) = 3
                  !data(nl,ny,nx)
                  else if (naxes(1).eq.dims(3).and.naxes(2) &
                      .eq.dims(2).and.naxes(3).eq.dims(1)) then
                    finfo(3) = 4
                  else
                    finfo(3) = 0
                  end if

                case(4)
                  !data(1,1,ny,nx)
                  if (naxes(1).eq.1.and.naxes(2).eq.1.and.naxes(3) &
                      .eq.dims(2).and.naxes(4).eq.dims(1)) then
                    finfo(3) = 3
                  !data(nl,1,ny,nx)
                  else if (naxes(1).eq.dims(3).and.naxes(2).eq.1 &
                      .and.naxes(3).eq.dims(2) &
                      .and.naxes(4).eq.dims(1)) then
                    finfo(3) = 4
                  else
                    finfo(3) = 0
                  end if

                case default
                  finfo(3) = 0

              end select ! Type of node

            ! 4 stokes
            else

              ! Check the size
              select case(naxis)
                case(1)
                  !data(4)
                  if (naxes(1).eq.4) then
                    finfo(3) = 1
                  else
                    finfo(3) = 0
                  end if

                case(2)
                  !data(1,4)
                  if (naxes(1).eq.1.and.naxes(2).eq.4) then
                    finfo(3) = 1
                  !data(nl,4)
                  else if (naxes(1).eq.dims(3) &
                      .and.naxes(2).eq.4) then
                    finfo(3) = 2
                  else
                    finfo(3) = 0
                  end if

                case(3)
                  !data(4,ny,nx)
                  if (naxes(1).eq.4.and.naxes(2).eq.dims(2) &
                      .and.naxes(3).eq.dims(1)) then
                    finfo(3) = 3
                  else
                    finfo(3) = 0
                  end if

                case(4)
                  !data(1,4,ny,nx)
                  if (naxes(1).eq.1.and.naxes(2).eq.4.and.naxes(3) &
                      .eq.dims(2).and.naxes(4).eq.dims(1)) then
                    finfo(3) = 3
                  !data(nl,4,ny,nx)
                  else if (naxes(1).eq.dims(3).and.naxes(2).eq.4 &
                      .and.naxes(3).eq.dims(2) &
                      .and.naxes(4).eq.dims(1)) then
                    finfo(3) = 4
                  else
                    finfo(3) = 0
                  end if

              case default
                finfo(3) = 0

              end select
            end if

          else

            finfo(3) = 0

          end if

          ! if more than 4 hdu
          if (hdunum.gt.4) then

          ! Move to hdu4 the diffuse dhu
            call FTMRHD(unitD,1,hdutype,stat)

            ! Get the dimension of the diffuse
            call FTGIDM(unitD,naxis,stat)

            if (naxis.gt.0) then

              ! Get the size of all dimensions
              call FTGISZ(unitD,naxis,naxes,stat)

              ! Get the date information
              ! Check the size
              select case(naxis)
                case(1)
                  !data(n1)
                  if (naxes(1).eq.dims(3)) then
                    finfo(4) = 1
                  else
                    finfo(4) = 0
                  end if

                case(2)
                  !data(nl,1)
                  if (naxes(1).eq.dims(3).and.naxes(2).eq.1) then
                    finfo(4) = 1
                  !data(nl,4)
                  else if (naxes(1).eq.dims(3) &
                      .and.naxes(2).eq.4) then
                    finfo(4) = 2
                  else
                    finfo(4) = 0
                  end if

                case(3)
                  !data(nl,ny,nx)
                  if (naxes(1).eq.dims(3).and.naxes(2).eq.dims(2) &
                      .and.naxes(3).eq.dims(1)) then
                    finfo(4) = 3
                  else
                    finfo(4) = 0
                  end if

                case(4)
                  !data(nl,1,ny,nx)
                  if (naxes(1).eq.dims(3).and. &
                      naxes(2).eq.1.and. &
                      naxes(3).eq.dims(2).and. &
                      naxes(4).eq.dims(1)) then
                    finfo(4) = 3
                  !data(nl,4,ny,nx)
                  else if (naxes(1).eq.dims(3).and. &
                      naxes(2).eq.4.and. &
                      naxes(3).eq.dims(2).and. &
                      naxes(4).eq.dims(1)) then
                    finfo(4) = 4
                  else
                    finfo(4) = 0
                  end if
                case default
                  finfo(4) = 0

              end select ! Type of node
            else
              finfo(4) = 0
            end if

          else
            finfo(4) = 0
          end if

        else
          finfo(3) = 0
          finfo(4) = 0

        end if

      ! Binary
      else
#endif

        ! Give data unit a number
        unitD = 19

        ! Open data file to read
        call open_file(unitD,Input%Filename_Ob,0,.False.,check)

        ! Check could open
        if (.not.check) then

          ! Issue error
          aborting = .True.
          return

        end if ! Could not open

        ! Get dimensions and info from input file
        call get_dims_info(unitD,dims,finfo,check)

        ! Check could read
        if (.not.check) then

          ! Issue error
          aborting = .True.
          return

        end if ! Could not read

#ifdef FITSSUP
      end if ! Type of file
#endif

      !
      ! Check if there is a cache
      !

      ! Check if there is a cache file
      inquire(file=trim(Input%cache), exist=lcache)

      ! If there is a cache
      if (lcache) then

        ! Allocate cache data
        allocate(cache(dims(2),dims(1)))

        ! Read cache
        call get_cache(unitC,Input%cache,dims(1:2),cache,check)

        ! The existence of the cache depends on if it could be read
        lcache = check

        ! If there is a cache
        if (lcache) then

          ! Check if there is something actually done
          if (.not.any(cache)) lcache = .False.

        end if

        ! If there is cache
        if (lcache) then

          ! Count memory
          MRAMc = MRAMc + 1d-6*sizeof(cache)

        ! No cache
        else

          ! No need to allocate it
          deallocate(cache)

        end if ! There is cache actually
      endif ! If there is cache

      ! Open or initialize cache
      call start_cache(unitC,Input%cache,dims(1:2),lcache,check)

      ! Check success
      aborting = .not.check

      ! Verbose
      write(umsg,'(A)') ' - Reading input inversion file '// &
                        trim(Input%Filename_Ob)
      call verbosev

      return

#ifdef FITSSUP
1000  call FTCLOS(unitD, stat)
      aborting = .True.
      return
#endif

      end subroutine open_data_and_cache

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open the mask file and check dimensions\n
      !!  Input(Input_class): Structure with configuration data\n
      !!      unitM(integer): Unit where to open mask\n
      !!   aborting(logical): Indicate failure at output\n
      !!    dims(integer(:)): Grid dimensions (X,Y,L)
      subroutine open_mask(Input,unitM,aborting,dims)

      ! I/O

      type(Input_class), intent(in):: Input
      logical, intent(out):: aborting
      integer, intent(out):: unitM
      integer, dimension(:), intent(in):: dims

      ! Local

      logical:: check

      integer, dimension(2):: ldims


      ! If Mask is none, return
      if (trim(Input%Inv_mask).eq.'NONE') return

      ! If inverting from scratch, return
      if (trim(Input%Inv_init).eq.'INIT') return

      ! Give mask unit a number
      unitM = 20

      ! Open data file to read
      call open_file(unitM,Input%Inv_mask,0,.False.,check)

      ! Check could open
      if (.not.check) then

        ! Issue error
        aborting = .True.
        return

      end if ! Coul open

      ! Read dimensions
      read(unitM,err=1000) ldims

      ! Check dimensions
      if (ldims(1).ne.dims(1).or.ldims(2).ne.dims(2)) then

        ! Issue error
        aborting = .True.
        umsg = ' # Mask file has wrong dimensions'
        call verbose

      end if

      return

      ! Issue error
1000  aborting = .True.
      return

      end subroutine open_mask

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the wavelength axis from the data file\n
      !!       unitD(integer): Unit where to read from\n
      !!    aborting(logical): Indicate failure at output\n
      !!  Sol(Solution_class): Structure with the frequency and
      !!                       synthetic Stokes parameters in the
      !!                       frequency range of the inverted data\n
      !!    fitsfile(logical): If the input is in a fits file
      subroutine get_data_wavelength(unitD,aborting,Sol,fitsfile)

      ! I/O

      type(Solution_class), intent(inout):: Sol
      logical, intent(in):: fitsfile
      logical, intent(out):: aborting
      integer, intent(in):: unitD

      ! Local

      logical:: anynull

#ifdef FITSSUP
      integer:: stat,hdutype,group,fpixel

      double precision:: nullval

      ! If fits
      if (fitsfile) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITS !!

        nullval = -1d0
        group = 1
        fpixel = 1
        stat = 0

        ! move to hdu1 (wavelength hdu)
        call FTMAHD(unitD,2,hdutype,stat)
        if(stat.ne.0) goto 1100

        ! read wavelength
        call FTGPVD(unitD,group,fpixel,Sol%Num_Wavelength,nullval, &
                    Sol%omega_input(1), anynull, stat)
        if(stat.ne.0) goto 1100
        if(anynull) goto 1100

      ! If binary
      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BINARY !!!!
#endif

        ! Read lambda
        read(unitD,err=1100) Sol%omega_input

#ifdef FITSSUP
      end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

      ! Success
      aborting = .False.

      return

      ! Issue error
1100  aborting = .True.
      write(umsg,'(A)') ' # Error getting wavelength '// &
                        'from the input inversion file'
      call verbose

      return

      ! Deceive compiler
      anynull = fitsfile

      end subroutine get_data_wavelength

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the LOS information, constant for the whole FoV, from
      !! the data file\n
      !!            unitD(integer): Unit where to read from\n
      !!         aborting(logical): Indicate failure at output\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!         fitsfile(logical): If the input is in a fits file
      subroutine get_data_los(unitD,aborting,Inf_Stokes,fitsfile)

      ! I/O

      type(Stokes_class), intent(inout):: Inf_Stokes
      logical, intent(in):: fitsfile
      logical, intent(out):: aborting
      integer, intent(in):: unitD

      ! Local

      logical:: anynull

#ifdef FITSSUP
      integer:: stat, hdutype,group,fpixel

      double precision:: nullval


      ! If fits
      if (fitsfile) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITS !!

        nullval = -1d0
        group = 1
        stat = 0

        ! move to hdu2 (mu hdu)
        call FTMAHD(unitD,3,hdutype,stat)
        if(stat.ne.0) goto 1100

        ! read mu
        fpixel = 1
        call FTGPVD(unitD,group,fpixel,1,nullval, &
                    Inf_Stokes%mu, anynull, stat)
        Inf_Stokes%mu = cos(Inf_Stokes%mu)
        if(stat.ne.0) goto 1100

        ! read azimuth
        fpixel = 2
        call FTGPVD(unitD,group,fpixel,1,nullval, &
                    Inf_Stokes%azimuth, anynull, stat)
        if(stat.ne.0) goto 1100

      ! If binary
      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BINARY !!!!
#endif

        ! Read LOS theta
        read(unitD,err=1100) Inf_Stokes%mu

        ! Save cosine in structure
        Inf_Stokes%mu = cos(Inf_Stokes%mu)

        ! Read LOS azimuth
        read(unitD,err=1100) Inf_Stokes%azimuth

#ifdef FITSSUP
      end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

      ! Success
      aborting = .False.

      return

      ! Issue error
1100  aborting = .True.
      write(umsg,'(A)') ' # Error getting constant LOS '// &
                        'from the input inversion file'
      call verbose

      ! Deceive compiler
      anynull = fitsfile

      end subroutine get_data_los

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the sigma (error) information, constant for the whole
      !! FoV, from the data file\n
      !!            unitD(integer): Unit where to read from\n
      !!         aborting(logical): Indicate failure at output\n
      !!         finfo(integer(:)): Data file metadata\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!         fitsfile(logical): If the input is in a fits file
      subroutine get_data_sigma(unitD,aborting,finfo, &
                                Inf_Stokes,fitsfile)

      ! I/O

      type(Stokes_class), intent(inout):: Inf_Stokes
      logical, intent(in):: fitsfile
      logical, intent(out):: aborting
      integer, intent(in):: unitD
      integer, dimension(3), intent(in):: finfo

      ! Local

      double precision:: daux

      logical:: anynull

#ifdef FITSSUP
      integer:: stat, hdutype,group,fpixel

      double precision:: nullval

      ! If fits
      if (fitsfile) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITS !!

        ! If only intensity
        if (finfo(1).eq.0) Inf_Stokes%Sigma_in(1:3,:) = 0d0

        nullval = -1d0
        group = 1
        fpixel = 1
        stat = 0

        ! move to hdu3 (mu hdu)
        call FTMAHD(unitD,4,hdutype,stat)
        if(stat.ne.0) goto 1100

        ! If constant
        if (finfo(3).eq.1) then

          ! Read sigmaI
          call FTGPVD(unitD,group,fpixel,1,nullval, &
                      daux, anynull, stat)
          if(stat.ne.0) goto 1100
          Inf_Stokes%Sigma_in(0,:) = daux

          ! If Polarization
          if (finfo(1).eq.1) then

            ! Read sigmaQ
            fpixel = 2
            call FTGPVD(unitD,group,fpixel,1,nullval, &
                        daux, anynull, stat)
            if(stat.ne.0) goto 1100
            Inf_Stokes%Sigma_in(1,:) = daux

            ! Read sigmaU
            fpixel = 3
            call FTGPVD(unitD,group,fpixel,1,nullval, &
                        daux, anynull, stat)
            if(stat.ne.0) goto 1100
            Inf_Stokes%Sigma_in(2,:) = daux

            ! Read sigmaV
            fpixel = 4
            call FTGPVD(unitD,group,fpixel,1,nullval, &
                        daux, anynull, stat)
            if(stat.ne.0) goto 1100
            Inf_Stokes%Sigma_in(3,:) = daux

          end if

        ! Variable
        else

          ! Read sigmaI
          call FTGPVD(unitD,group,fpixel,Inf_Stokes%Num_Wavelength, &
                      nullval,Inf_Stokes%Sigma_in(0,1), anynull, stat)
          if(stat.ne.0) goto 1100

          ! If Polarization
          if (finfo(1).eq.1) then

            ! Read sigmaQ
            fpixel = fpixel+Inf_Stokes%Num_Wavelength
            call FTGPVD(unitD,group,fpixel, &
                        Inf_Stokes%Num_Wavelength,nullval, &
                        Inf_Stokes%Sigma_in(1,1), anynull, stat)
            if(stat.ne.0) goto 1100

            ! Read sigmaU
            fpixel = fpixel+Inf_Stokes%Num_Wavelength
            call FTGPVD(unitD,group,fpixel, &
                        Inf_Stokes%Num_Wavelength,nullval, &
                        Inf_Stokes%Sigma_in(2,1), anynull, stat)
            if(stat.ne.0) goto 1100

            ! Read sigmaV
            fpixel = fpixel+Inf_Stokes%Num_Wavelength
            call FTGPVD(unitD,group,fpixel, &
                        Inf_Stokes%Num_Wavelength,nullval, &
                        Inf_Stokes%Sigma_in(3,1), anynull, stat)
            if(stat.ne.0) goto 1100

          end if

        end if ! Constant/wavelength dependent

      ! If binary
      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BINARY !!!!
#endif

        ! If only intensity
        if (finfo(1).eq.0) Inf_Stokes%Sigma_in(1:3,:) = 0d0

        ! If constant
        if (finfo(3).eq.1) then

          ! Read sigmaI and save in structure
          read(unitD,err=1100) daux
          Inf_Stokes%Sigma_in(0,:) = daux

          ! If Polarization
          if (finfo(1).eq.1) then

            ! Read sigmaQ and save in structure
            read(unitD,err=1100) daux
            Inf_Stokes%Sigma_in(1,:) = daux

            ! Read sigmaU and save in structure
            read(unitD,err=1100) daux
            Inf_Stokes%Sigma_in(2,:) = daux

            ! Read sigmaV and save in structure
            read(unitD,err=1100) daux
            Inf_Stokes%Sigma_in(3,:) = daux

          end if

        ! Constant but wavelength dependent
        else

          ! Read sigmaI
          read(unitD,err=1100) Inf_Stokes%Sigma_in(0,:)

          ! If Polarization
          if (finfo(1).eq.1) then

            ! Read sigmaQUV
            read(unitD,err=1100) Inf_Stokes%Sigma_in(1,:)
            read(unitD,err=1100) Inf_Stokes%Sigma_in(2,:)
            read(unitD,err=1100) Inf_Stokes%Sigma_in(3,:)

          end if ! Polarization
        end if ! Constant/wavelength dependent

#ifdef FITSSUP
      end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

      ! Success
      aborting = .False.

      return

      ! Issue error
1100  aborting = .True.
      write(umsg,'(A)') ' # Error getting constant sigma '// &
                        'from the input inversion file'
      call verbose

      return

      ! Deceive compiler
      anynull = fitsfile

      end subroutine get_data_sigma

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the diffuse light profile, constant for the whole FoV,
      !! from the data file\n
      !!            unitD(integer): Unit where to read from\n
      !!         aborting(logical): Indicate failure at output\n
      !!         finfo(integer(:)): Data file metadata\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!         fitsfile(logical): If the input is in a fits file
      subroutine get_data_diff(unitD,aborting,finfo,Inf_Stokes, &
                               fitsfile)

      ! I/O

      type(Stokes_class), intent(inout):: Inf_Stokes
      logical, intent(in):: fitsfile
      logical, intent(out):: aborting
      integer, intent(in):: unitD
      integer, dimension(4), intent(in):: finfo

      ! Local

      logical:: anynull

#ifdef FITSSUP
      integer:: stat, hdutype,group,fpixel

      double precision:: nullval


      ! If fits
      if (fitsfile) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITS !!

        ! If only intensity
        if (finfo(4).eq.1) Inf_Stokes%Diff_in(1:3,:) = 0d0

        nullval = -1d0
        group = 1
        fpixel = 1
        stat = 0

        ! move to hdu4 (diffuse light hdu)
        call FTMAHD(unitD,5,hdutype,stat)
        if(stat.ne.0) goto 1100

        ! read intensity diffuse light
        call FTGPVD(unitD,group,fpixel,Inf_Stokes%Num_Wavelength, &
                    nullval,Inf_Stokes%Diff_in(0,1), anynull, stat)
        if(stat.ne.0) goto 1100

        ! If Polarization
        if (finfo(4).eq.2) then

          ! Read diffuse Q
          fpixel = fpixel+Inf_Stokes%Num_Wavelength
          call FTGPVD(unitD,group,fpixel, &
                      Inf_Stokes%Num_Wavelength,nullval, &
                      Inf_Stokes%Diff_in(1,1), anynull, stat)
          if(stat.ne.0) goto 1100

          ! Read diffuse U
          fpixel = fpixel+Inf_Stokes%Num_Wavelength
          call FTGPVD(unitD,group,fpixel, &
                      Inf_Stokes%Num_Wavelength,nullval, &
                      Inf_Stokes%Diff_in(2,1), anynull, stat)
          if(stat.ne.0) goto 1100

          ! Read diffuse V
          fpixel = fpixel+Inf_Stokes%Num_Wavelength
          call FTGPVD(unitD,group,fpixel, &
                      Inf_Stokes%Num_Wavelength,nullval, &
                      Inf_Stokes%Diff_in(3,1), anynull, stat)
          if(stat.ne.0) goto 1100

        end if

      ! If binary
      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BINARY !!!!
#endif

        ! If only intensity, set polarization to zero
        if (finfo(4).eq.1) Inf_Stokes%Diff_in(1:3,:) = 0d0

        ! Read diffuse light intensity profile
        read(unitD,err=1100) Inf_Stokes%Diff_in(0,:)

        ! If Polarization
        if (finfo(4).eq.2) then

          ! Read diffuse light QUV profiles
          read(unitD,err=1100) Inf_Stokes%Diff_in(1,:)
          read(unitD,err=1100) Inf_Stokes%Diff_in(2,:)
          read(unitD,err=1100) Inf_Stokes%Diff_in(3,:)

        end if ! Polarization

#ifdef FITSSUP
      end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

      ! Success
      aborting = .False.

      return

      ! Issue error
1100  aborting = .True.
      write(umsg,'(A)') ' # Error getting constant diffuse light '// &
                        'from the input inversion file'
      call verbose

      return

      ! Deceive compiler
      anynull = fitsfile

      end subroutine get_data_diff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read one pixel in the data file\n
      !!  Input(Input_class): Structure with configuration data\n
      !!      unitD(integer): Unit where to read from\n
      !!   buffer(double(:)): Buffer to store the data\n
      !!    dims(integer(:)): Grid dimensions (X,Y,L)\n
      !!   finfo(integer(:)): Data file metadata\n
      !!         ix(integer): Pixel index in x dimension\n
      !!         iy(integer): Pixel index in y dimension\n
      !!      check(logical): If read is a success
      subroutine get_data_column(Input,unitD,dims,buffer,finfo, &
                                 ix,iy,check)

      ! I/O

      type(Input_class), intent(in):: Input
      logical, intent(out):: check
      integer, intent(in):: unitD,ix,iy
      integer, dimension(:), intent(in):: finfo,dims
      double precision, dimension(:), intent(out):: buffer

      ! Local

      integer:: i0,i1,nwp

      double precision:: daux

#ifdef FITSSUP
      logical:: anynull

      integer:: stat,hdutype,group,fpixel

      double precision:: nullval
#endif

      ! Size of wavelength package
      if (finfo(1).eq.0) then

        ! Only intensity
        nwp = dims(3)

      ! Polarization
      else

        ! Polarization
        nwp = dims(3)*4

      end if ! Size of wavelength package

#ifdef FITSSUP

      !
      ! If FITS file
      !
      if (Input%FITSFILE) then !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FITS !!

        nullval = -1d0
        group = 1
        stat = 0

        ! If LOS
        if (finfo(2).eq.1) then

          ! move to hdu2 (mu hdu)
          call FTMAHD(unitD,3,hdutype,stat)
          if(stat.ne.0) goto 1100

          ! read mu
          fpixel = 1+2*((ix-1)*dims(2)+(iy-1))
          call FTGPVD(unitD,group,fpixel,1,nullval, &
                      daux, anynull, stat)
          if(stat.ne.0) goto 1100
          buffer(1) = cos(daux)

          ! read azimuth
          fpixel = fpixel+1
          call FTGPVD(unitD,group,fpixel,1,nullval, &
                      daux, anynull, stat)
          if(stat.ne.0) goto 1100
          buffer(2) = daux

          ! Receiving indexes
          i0 = 3
          i1 = nwp + i0 - 1

        ! No LOS
        else

          ! Receiving indexes
          i0 = 1
          i1 = nwp

        end if

        ! move to hdu0 (stokes hdu)
        call FTMAHD(unitD,1,hdutype,stat)
        if(stat.ne.0) goto 1100

        ! read stokes
        fpixel = 1+nwp*((ix-1)*dims(2)+(iy-1))

        call FTGPVD(unitD,group,fpixel,nwp,nullval,buffer(i0), &
                    anynull,stat)
        if(stat.ne.0) goto 1100

        ! There is sigma here
        if (finfo(3).ge.3) then

          ! move to hdu3 (stokes hdu)
          call FTMAHD(unitD,4,hdutype,stat)
          if(stat.ne.0) goto 1100

          ! Shift initial index
          i0 = i1 + 1

          ! If constant
          if (finfo(3).eq.3) then

            ! Intensity
            if (finfo(1).eq.0) then

              ! Final index
              i1 = i0
              fpixel = 1+dims(3)*((ix-1)*dims(2)+(iy-1))
              call FTGPVD(unitD,group,fpixel,1,nullval,buffer(i0), &
                          anynull,stat)
              if(stat.ne.0) goto 1100

            ! Polarization
            else
              ! Final index
              i1 = i0+3
              fpixel = 1+4*((ix-1)*dims(2)+(iy-1))
              call FTGPVD(unitD,group,fpixel,4,nullval,buffer(i0), &
                          anynull,stat)
              if(stat.ne.0) goto 1100

            end if

          ! Wavelength dependent
          else

            ! Final index
            i1 = nwp + i0 - 1
            fpixel = 1+nwp*((ix-1)*dims(2)+(iy-1))

            ! Read
            call FTGPVD(unitD,group,fpixel,nwp,nullval,buffer(i0), &
                        anynull,stat)
            if(stat.ne.0) goto 1100

          end if ! Type of pixel sigma

        end if

        ! There is diffuse light here
        if (finfo(4).ge.3) then

          ! Shift initial index
          i0 = i1 + 1

          ! If only intensity
          if (finfo(4).eq.3) then

            ! Read
            fpixel = 1+dims(3)*((ix-1)*dims(2)+(iy-1))
            call FTGPVD(unitD,group,fpixel,dims(3),nullval, &
                        buffer(i0),anynull,stat)
            if(stat.ne.0) goto 1100

          ! Polarized
          else

            ! Read
            fpixel = 1+nwp*((ix-1)*dims(2)+(iy-1))
            call FTGPVD(unitD,group,fpixel,nwp,nullval,buffer(i0), &
                        anynull,stat)
            if(stat.ne.0) goto 1100

          end if ! Polarization?
        end if ! If pixel diffuse light

      ! Binary
      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BINARY !!!!
#endif

        ! If LOS is pixel dependent
        if (finfo(2).eq.1) then

          ! Read mu
          read(unitD,err=1100) daux
          buffer(1) = cos(daux)

          ! Read azimuth
          read(unitD,err=1100) daux
          buffer(2) = daux

          ! Receiving indexes for Stokes
          i0 = 3
          i1 = nwp + i0 - 1

        ! No pixel-wise LOS
        else

          ! Receiving indexes for Stokes
          i0 = 1
          i1 = nwp

        end if ! LOS pixel dependent

        ! Read Stokes
        read(unitD,err=1100) buffer(i0:i1)

        ! Sigma is pixel-wise
        if (finfo(3).ge.3) then

          ! Shift initial index
          i0 = i1 + 1

          ! If not wavelength dependent
          if (finfo(3).eq.3) then

            ! Only intensity
            if (finfo(1).eq.0) then

              ! Final index
              i1 = i0

            ! Polarization
            else

              ! Final index
              i1 = i0 + 3

            end if ! Intensity/Polarization

            ! Read sigma
            read(unitD,err=1100) buffer(i0:i1)

          ! Wavelength dependent
          else

            ! Final index
            i1 = nwp + i0 - 1

            ! Read
            read(unitD,err=1100) buffer(i0:i1)

          end if ! Type of pixel sigma
        end if ! If pixel sigma

        ! Diffuse light is pixel-wise
        if (finfo(4).ge.3) then

          ! Shift initial index
          i0 = i1 + 1

          ! If only intensity
          if (finfo(4).eq.3) then

            ! Final index
            i1 = dims(3) + i0 - 1

            ! Read intensity
            read(unitD,err=1100) buffer(i0:i1)

          ! Polarized
          else

            ! Final index
            i1 = nwp + i0 - 1

            ! Read
            read(unitD,err=1100) buffer(i0:i1)

          end if ! Polarization?
        end if ! If pixel diffuse light
#ifdef FITSSUP
      end if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

      ! Success
      check = .True.

      return

      ! Issue error
1100  check = .False.
      write(umsg,'(A)') ' # Error getting data from pixel '// &
                        'from the input inversion file'
      call verbose

      return

      ! Deceive compiler
      i1 = Input%iter_min
      i1 = ix
      i1 = iy

      end subroutine get_data_column

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform data in buffer into manageable arrays\n
      !!         finfo(integer(:)): Data file metadata\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!         buffer(double(:)): Data buffer
      subroutine set_up_data_frombuffer(finfo,Inf_Stokes,Sol,buffer)

      ! I/O

      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Solution_class), intent(inout):: Sol
      integer, dimension(:), intent(in):: finfo
      double precision, dimension(:), intent(in):: buffer

      ! Local

      integer:: i0,i1,istk,nstk

      ! Pixel-wise LOS
      if (finfo(2).eq.1) then

        ! Get LOS from buffer
        Inf_Stokes%mu = buffer(1)
        Inf_Stokes%azimuth = buffer(2)

        ! Set final index of this section
        i1 = 2

      ! Known constant LOS
      else

        ! Nothing read from buffer
        i1 = 0

      end if ! Type of LOS

      ! If only intensity
      if (finfo(1).eq.0) then

        ! Last Stokes index
        nstk = 0

        ! Zero Pol
        Inf_Stokes%Stokes_Ob(1:3,:) = 0d0

      ! Polarization
      else

        ! Last Stokes index
        nstk = 3

      end if

      ! For each stokes parameter
      do istk=0,nstk

        ! Set indexes
        i0 = i1 + 1
        i1 = i0 + Inf_Stokes%Num_Wavelength - 1

        ! Get Stokes parameter
        Inf_Stokes%Stokes_Ob(istk,:) = buffer(i0:i1)

      end do ! Stokes parameters

      !
      ! Sigma values
      !

      ! No sigma
      if (finfo(3).eq.0) then

        ! Flag as false
        Inf_Stokes%Sigma_Flag = .False.

      ! Yes sigma
      else

        ! Flag as True
        Inf_Stokes%Sigma_Flag = .True.

        ! Constant sigma
        if (finfo(3).eq.1.or.finfo(3).eq.2) then

          ! Copy from the input sigma
          Inf_Stokes%Sigma_W = Inf_Stokes%Sigma_in

        ! Pixel sigma not wavelength dependent
        else if (finfo(3).eq.3) then

          ! For each Stokes parameter
          do istk=0,nstk

            ! Set indexes
            i0 = i1 + 1
            i1 = i0

            ! Copy replicating
            Inf_Stokes%Sigma_W(istk,:) = buffer(i0)

          end do ! Stokes parameters

        ! Pixel sigma wavelength dependent
        else if (finfo(3).eq.4) then

          ! For each Stokes parameter
          do istk=0,nstk

            ! Set indexes
            i0 = i1 + 1
            i1 = i0 + Inf_Stokes%Num_Wavelength - 1

            ! Copy
            Inf_Stokes%Sigma_W(istk,:) = buffer(i0:i1)

          end do ! Stokes paramters

        end if ! Type of sigma

        ! For each Stokes parameter
        do istk=0,nstk

          ! Check if there are negatives
          if (minval(Inf_Stokes%Sigma_W(istk,:)).le.0d0) then

            ! I
            if (istk.eq.0) then

                ! Message
                umsg = 'At least one "sigma" value for '// &
                       'intensity is not positive'

            ! Q
            else if (istk.eq.1) then

                ! Message
                umsg = 'At least one "sigma" value for '// &
                       'Q is not positive'

            ! U
            else if (istk.eq.2) then

                ! Message
                umsg = 'At least one "sigma" value for '// &
                       'U is not positive'

            ! V
            else if (istk.eq.3) then

                ! Message
                umsg = 'At least one "sigma" value for '// &
                       'V is not positive'

            end if

            ! Issue error
            urou = 'set_up_data_frombuff'
            call aborted
            exit

          end if ! Negative sigma

        end do ! For each Stokes parameter

        ! Share
        call control

      end if ! If there is sigma at all


      !
      ! Diffuse light values
      !

      ! None
      if (finfo(4).eq.0) then

        ! Flag as false
        Inf_Stokes%Diff_Flag = .False.
        Sol%Diff_Flag = .False.

      ! There is diffuse light
      else

        ! Flag as True
        Inf_Stokes%Diff_Flag = .True.
        Sol%Diff_Flag = .True.

        ! Constant
        if (finfo(4).eq.1.or.finfo(4).eq.2) then

          ! Copy from the input sigma
          Sol%Stokes_diff = Inf_Stokes%Diff_in

        ! Pixel-wise only intensity
        else if (finfo(4).eq.3) then

          ! Set indexes
          i0 = i1 + 1
          i1 = i0 + Inf_Stokes%Num_Wavelength - 1

          ! Copy
          Sol%Stokes_diff(0,:) = buffer(i0:i1)
          Sol%Stokes_diff(1:3,:) = 0d0

        ! Pixel-wise and polarized
        else if (finfo(4).eq.4) then

          ! For each Stokes parameter
          do istk=0,nstk

            ! Set indexes
            i0 = i1 + 1
            i1 = i0 + Inf_Stokes%Num_Wavelength - 1

            ! Copy
            Sol%Stokes_diff(istk,:) = buffer(i0:i1)

          end do ! Stokes paramters

        end if ! Type of diffuse light
      end if ! If there is diffuse light at all

      ! Convert profiles to correct format
      call Profile_Conversion(Inf_Stokes, Sol)

      end subroutine set_up_data_frombuffer

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform atmosphere in buffer into the structure with
      !! atmospheric data\n
      !!         finv(logical): If the model was read from an
      !!                        inversion file\n
      !!        jkqin(logical): If there are JKQ\n
      !!      tauscal(logical): If input is in tau scale\n
      !!     buffer(double(:)): Atmosphere buffer\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine set_up_atmo_frombuffer(finv,jkqin,tauscal, &
                                        buffer,Atmo,Bfield)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(out):: Bfield
      logical, intent(in):: finv,tauscal,jkqin
      double precision, dimension(:), intent(in), target:: buffer

      ! Local

      integer:: lnz,iz

      double precision:: ikbcgs


      ! nz local
      lnz = Atmo%nz

      ! Constant
      ikbcgs = 1d-7/kb

      ! Ensure necessary free space in atmospheric structure
      if (allocated(Atmo%Pg)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pg)
        deallocate(Atmo%Pg)
      end if
      if (allocated(Atmo%nH)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nH)
        deallocate(Atmo%nH)
      end if
      if (allocated(Atmo%nht)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nht)
        deallocate(Atmo%nht)
      end if
      if (allocated(Atmo%nha)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nha)
        deallocate(Atmo%nha)
      end if
      if (allocated(Atmo%nhm)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nhm)
        deallocate(Atmo%nhm)
      end if
      if (allocated(Atmo%ne)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%ne)
        deallocate(Atmo%ne)
      end if

      ! Initialize zero array if not allocated
      if (.not.associated(Atmo%zeros)) then
        allocate(Atmo%zeros(lnz))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%zeros)
        Atmo%zeros = 0d0
      end if

      ! Set type of allocations for sets of pointers
      Atmo%alloc_a = .False.
      Atmo%alloc_b = .False.

      ! Allocation
      allocate(Atmo%nH(lnz,6),Atmo%nht(lnz),Atmo%nha(lnz))
      allocate(Atmo%nhm(lnz),Atmo%ne(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nH)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nht)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nha)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhm)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%ne)

      ! If from inversion file
      if (finv) then

        ! Type of scale must be tau
        ztau = .True.

        ! Tau
        Atmo%z   => buffer(1:lnz)

        ! Temperature
        Atmo%T   => buffer(   lnz+1:2*lnz)

        ! Microturbulence
        Atmo%vmi => buffer( 9*lnz+1:10*lnz)

        ! Velocity
        Atmo%vx  => buffer( 6*lnz+1: 7*lnz)
        Atmo%vy  => buffer( 7*lnz+1: 8*lnz)
        Atmo%vz  => buffer( 8*lnz+1: 9*lnz)

        ! Reallocate gas pressure
        if (allocated(Atmo%Pg)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pg)
          deallocate(Atmo%Pg)
        end if
        allocate(Atmo%Pg(lnz))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)

        ! Gas pressure
        Atmo%Pg = buffer( 2*lnz+1: 3*lnz)

        ! Type of atmosphere
        Atmo%typo = 4

        ! Initialize
        Atmo%ne = 0d0
        Atmo%nht = 0d0
        Atmo%nha = 0d0
        Atmo%nH = 0d0
        Atmo%nhm = 0d0

      ! Not inversion file
      else

        ! Tau scale
        if (tauscal) then

          ! Point to tau
          Atmo%z   => buffer( 2*lnz+1:3*lnz)
          ztau = .True.

        ! Height scale
        else

          ! Point to z
          Atmo%z   => buffer(   lnz+1:2*lnz)
          ztau = .False.

        end if ! Type of scale

        ! Temperature
        Atmo%T   => buffer( 3*lnz+1:4*lnz)

        ! Microturbulence
        Atmo%vmi => buffer(12*lnz+1:13*lnz)

        ! Velocity
        Atmo%vx  => buffer( 9*lnz+1:10*lnz)
        Atmo%vy  => buffer(10*lnz+1:11*lnz)
        Atmo%vz  => buffer(11*lnz+1:12*lnz)

        ! Type of atmosphere
        if (Atmo%typo.eq.0.or.Atmo%typo.eq.1) then

          ! Copy electron number density
          Atmo%ne = buffer(14*lnz+1:15*lnz)

          ! Full density
          if (Atmo%typo.eq.0) &
            Atmo%nH = reshape(buffer(18*nz+1:24*nz), (/ nz, 6 /))

        ! Electron pressure
        else if (Atmo%typo.eq.2) then

          ! Copy electron pressure
          Atmo%ne = buffer(13*lnz+1:14*lnz)

          ! Transform to electron number density
          Atmo%ne = Atmo%ne*ikbcgs/Atmo%T
          Atmo%typo = 1

        ! Gas pressure
        else if (Atmo%typo.eq.4) then

          ! Reallocate gas pressure
          if (allocated(Atmo%Pg)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pg)
            deallocate(Atmo%Pg)
          end if
          allocate(Atmo%Pg(lnz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)

          ! Initialize electron number density
          Atmo%ne = 0d0

          ! Copy from buffer
          Atmo%Pg = buffer( 4*lnz+1:5*lnz)

        ! Mass density
        else if (Atmo%typo.eq.5) then

          ! Reallocate mass density
          if (allocated(Atmo%rho)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%rho)
            deallocate(Atmo%rho)
          end if
          allocate(Atmo%rho(lnz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%rho)

          ! Initialize electron number density
          Atmo%ne = 0d0

          ! Copy from buffer
          Atmo%rho = buffer(5*lnz+1:6*lnz)

        end if ! Type of model atmosphere
      end if ! Input model format

      ! Just make a flag to know that there is no helium input
      allocate(Atmo%nhe(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
      Atmo%nhe(1,1) = -1

      ! Divide velocities by c (1d5*1d-11/cbar)
      Atmo%vx = Atmo%vx*1d-6/c
      Atmo%vy = Atmo%vy*1d-6/c
      Atmo%vz = Atmo%vz*1d-6/c
      Atmo%vmi = Atmo%vmi*1d-6/c


      !
      ! Magnetic field
      !

      ! Allocate
      allocate(Bfield%Bstrength(lnz))
      allocate(Bfield%Btheta(lnz))
      allocate(Bfield%Bphi(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bstrength)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Btheta)
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bphi)

      ! Point to buffer
      Atmo%Bx => buffer( 3*lnz+1:7*lnz)
      Atmo%By => buffer( 4*lnz+1:8*lnz)
      Atmo%Bz => buffer( 5*lnz+1:9*lnz)

      ! Compute module and angles for B at each height
      do iz=1,lnz

        ! Module
        Bfield%Bstrength(iz) = sqrt(Atmo%Bx(iz)*Atmo%Bx(iz) + &
                                    Atmo%By(iz)*Atmo%By(iz) + &
                                    Atmo%Bz(iz)*Atmo%Bz(iz))

        ! There is a magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          ! Get angles
          Bfield%Btheta(iz) = acos(Atmo%Bz(iz)/Bfield%Bstrength(iz))
          Bfield%Bphi(iz) = atan2(Atmo%By(iz),Atmo%Bx(iz))

        ! There is no field
        else

          ! Trivial angles
          Bfield%Bstrength(iz) = 0d0
          Bfield%Btheta(iz) = 0d0
          Bfield%Bphi(iz) = 0d0

        end if ! If there is magnetic field

      end do ! Heights

      ! If there are ad-hoc JKQ
      if (jkqin) then

        ! Get from buffer
        Atmo%JKQin(       1:  lnz) = buffer(19*lnz+1:20*lnz)
        Atmo%JKQin(   lnz+1:2*lnz) = buffer(20*lnz+1:21*lnz)
        Atmo%JKQin( 2*lnz+1:3*lnz) = buffer(21*lnz+1:22*lnz)
        Atmo%JKQin( 3*lnz+1:4*lnz) = buffer(22*lnz+1:23*lnz)
        Atmo%JKQin( 4*lnz+1:5*lnz) = buffer(23*lnz+1:24*lnz)
        Atmo%JKQin( 5*lnz+1:6*lnz) = buffer(24*lnz+1:25*lnz)
        Atmo%JKQin( 6*lnz+1:7*lnz) = buffer(25*lnz+1:26*lnz)
        Atmo%JKQin( 7*lnz+1:8*lnz) = buffer(26*lnz+1:27*lnz)

      end if ! Ad-hoc JKQ

      ! Inversion file
      if (finv) then

        ! If there are JKQ
        if (jkqin) then

          ! Diffuse light
          Atmo%f_diff = buffer(27*lnz+1)

        ! No JKQ
        else

          ! Diffuse light
          Atmo%f_diff = buffer(19*lnz+1)

        end if ! If there are JKQ
      end if ! If it is a model from an inversion file

      ! Nullify Bx,y, and z
      nullify(Atmo%Bx,Atmo%By,Atmo%Bz)

      end subroutine set_up_atmo_frombuffer

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform JKQ in buffer into the structure with atmospheric
      !! data\n
      !!  buffer(double(:)): JKQin buffer\n
      !!   Atmo(Atmo_class): Structure with atmospheric data
      subroutine set_up_JKQ_frombuffer(buffer,Atmo)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      double precision, dimension(:), intent(in), target:: buffer

      ! Copy from buffer
      Atmo%JKQin = buffer

      end subroutine set_up_JKQ_frombuffer

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate and store the sizes of the data to navigate the
      !! output files\n
      !!      Input(Input_class): Structure with configuration data\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!        dims(integer(:)): Grid dimensions (X,Y,L)\n
      !!             nz(integer): Vertical dimension of the model
      !!                          atmosphere
      subroutine set_inv_io_buffers(Input,Inf_Nodes,dims,nz)

      ! I/O

      type(Input_class), intent(inout):: Input
      type(Nodes_class), intent(in):: Inf_Nodes
      integer, intent(in):: nz
      integer, dimension(:), intent(in):: dims

      ! Local

      integer:: ivar


      ! Size of common header
      Input%s_inv_h= 20

      !
      ! Size of atmosphere block
      !

      ! Column per variable
      Input%s_inv_atmo_c = 4*nz

      ! If output JKQ
      if (Input%out_jkqa) then

        ! Size of each column of atmosphere
        Input%s_inv_atmo_c = Input%s_inv_atmo_c*27

      ! No output JKQ
      else

        ! Size of each column of atmosphere
        Input%s_inv_atmo_c = Input%s_inv_atmo_c*19

      end if ! Output JKQ

      ! Diffuse light
      Input%s_inv_atmo_c = Input%s_inv_atmo_c + 4*1

      ! Total block
      Input%s_inv_atmo = dims(1)*dims(2)*Input%s_inv_atmo_c


      !
      ! Sizes for results
      !

      ! Head (11 number of variables)         v 4 + 4
      Input%s_inv_res_h = 4 + 4 + 8*dims(3) + 8*Input%nvar

      ! If thermal
      if (Input%Type_Inversion.eq.0) then

        ! Each pixel
        Input%s_inv_res_c = 8 + 8*dims(3)

      ! Magnetic
      else

        ! Each pixel
        Input%s_inv_res_c = 8 + 32*dims(3)

      end if ! Thermal/magnetic

      ! Run over variables
      do ivar=1,Input%nvar

        ! If there are nodes
        if (Inf_Nodes%Num_nodes(ivar).gt.0) then

          ! If inverting
          if (Inf_Nodes%Nodes_Flags(ivar)) then

            ! Add node size
            Input%s_inv_res_c = Input%s_inv_res_c + &
                                   12*Inf_Nodes%Num_nodes(ivar)

          ! No inverting
          else

            ! Add node size without error
            Input%s_inv_res_c = Input%s_inv_res_c + &
                                   8*Inf_Nodes%Num_nodes(ivar)

          end if ! Inverting?
        end if ! Nodes?

      end do ! Variables

      ! Total
      Input%s_inv_res = Input%s_inv_res_c*dims(1)*dims(2)

      ! If output response functions
      if (Input%Keep_RF) then

        ! Head (number of variables with RF)
        Input%s_inv_RF_h = 4

        ! Column
        Input%s_inv_RF_c = 0

        ! Thermal
        if (Input%Type_Inversion.eq.0) then

          ! Run over variables
          do ivar=1,Input%nvar

            ! Skip not inverting
            if (.not.Inf_Nodes%Nodes_Flags(ivar)) cycle

            ! Add to header index and size
            Input%s_inv_RF_h = Input%s_inv_RF_h + 4 + 4

            ! If not varying
            if (Inf_Nodes%Num_vary(ivar).le.0) cycle

            ! Add node size
            Input%s_inv_RF_c = Input%s_inv_RF_c + &
                                  (4 + 4*dims(3))* &
                                  Inf_Nodes%Num_vary(ivar)

          end do ! Variables

        ! Polarization
        else

          ! Run over variables
          do ivar=1,Input%nvar

            ! Skip not inverting
            if (.not.Inf_Nodes%Nodes_Flags(ivar)) cycle

            ! Add to header index and size
            Input%s_inv_RF_h = Input%s_inv_RF_h + 4 + 4

            ! If not varying
            if (Inf_Nodes%Num_vary(ivar).le.0) cycle

            ! Add node size
            Input%s_inv_RF_c = Input%s_inv_RF_c + &
                                  (4 + 4*4*dims(3))* &
                                  Inf_Nodes%Num_vary(ivar)

          end do ! Variables

        end if ! If polarization
      end if ! Output RF?

      end subroutine set_inv_io_buffers

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if the output files already exist\n
      !!  Input(Input_class): Structure with configuration data\n
      !!   aborting(logical): Signals if something goes wrong
      subroutine check_io_inv_buffer_exists(Input,aborting)

      ! I/O

      type(Input_class), intent(inout):: Input
      logical, intent(out):: aborting

      ! Local

      integer:: ios


      ! Initialize
      aborting = .False.

      !
      ! Stokes, Contribution, Tau
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

      end if ! Could not open

      ! Close
      close(200)

      ! If contribution function output
      if (Input%out_contr) then

        ! Try opening Contribution function
        open (200,file=trim(Input%folder)//'/Contribution', &
              status='old', iostat=ios, access='stream', &
              action='read', form='unformatted')

        ! If could not open
        if (ios.ne.0) then

          ! Issue error
          umsg = 'There is no existing contribution '// &
                 'function files'
          call verbose
          aborting = .True.
          return

        end if ! Could not open

        ! Close
        close(200)

      end if ! Output contribution function

      ! If tau output
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

        end if ! Could not open

        ! Close
        close(200)

      end if ! Output tau

      ! Try opening result
      open (200,file=trim(Input%folder)//'/Result', &
            status='old', iostat=ios, access='stream', &
            action='read', form='unformatted')

      ! If could not open
      if (ios.ne.0) then

        ! Issue error
        umsg = 'There is no existing Result file'
        call verbose
        aborting = .True.
        return

      end if ! Could not open

      ! Close
      close(200)

      end subroutine check_io_inv_buffer_exists

!#####################################################################
!#####################################################################
!#####################################################################

      !> Open and write the header for the output files\n
      !!      Input(Input_class): Structure with configuration data\n
      !!     Sol(Solution_class): Structure with the frequency and
      !!                          synthetic Stokes parameters in the
      !!                          frequency range of the inverted
      !!                          data\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!              th(double): LOS polar angle\n
      !!             azi(double): LOS azimuth angle\n
      !!        dims(integer(:)): Grid dimensions (X,Y,L)\n
      !!             nz(integer): Vertical dimension of the model
      !!                          atmosphere\n
      !!   Frec(Frequency_class): Structure with frequency data
      subroutine create_io_inv_files(Input,Sol,Inf_Nodes,th,azi, &
                                     dims,nz,Frec)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Nodes_class), intent(in):: Inf_Nodes
      type(Frequency_class), intent(in):: Frec
      type(Solution_class), intent(in):: Sol
      integer, intent(in):: nz
      integer, dimension(:), intent(in):: dims
      double precision, intent(in):: th,azi

      ! Local

      integer:: i,j,ios,iran


      !
      ! Stokes, Contribution, Tau
      !

      ! Open Stokes
      open (200,file=trim(Input%folder)//'/Stokes', &
            status='unknown', iostat=ios, access='stream', &
            action='write', form='unformatted')

      ! Write header
      write(200) '2Dbe'
      write(200) Input%lim_stk%nn
      write(200) dims(1:2)
      write(200) th
      write(200) azi

      ! If specific ranges
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

      end if ! Specific ranges

      ! Close
      close(200)

      ! If contribution function output
      if (Input%out_contr) then

        ! Open Contribution function
        open (200,file=trim(Input%folder)//'/Contribution', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dbc'
        write(200) Input%lim_ctr%nn
        write(200) dims(1:2)
        write(200) nz
        write(200) th
        write(200) azi

        ! If specific ranges
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

        end if ! Specific ranges

        ! Close
        close(200)

      end if ! Output contribution function

      ! If tau output
      if (Input%out_tau1) then

        ! Open tau
        open (200,file=trim(Input%folder)//'/Tau', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write header
        write(200) '2Dbt'
        write(200) Input%lim_tau%nn
        write(200) dims(1:2)
        write(200) th
        write(200) azi

        ! If specific ranges
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

        end if ! Specific ranges

        ! Close
        close(200)

      end if ! Output tau

      !
      ! Result
      !

      ! Open result file
      open (200,file=trim(Input%folder)//'/Result', &
            status='unknown', iostat=ios, access='stream', &
            action='write', form='unformatted')

      !
      ! Write header
      !

      ! Label
      write(200) 'invo'

      ! Initialize info index
      i = 0

      ! If magnetic
      if (Input%type_inversion.gt.0) i = i + 1

      ! Type of field
      i = i + Input%btype*2

      ! Type of velocity
      i = i + Input%vtype*4

      ! JKQ
      if (Input%out_jkqa) i = i + 8

      ! Write info index
      write(200) i

      ! Dimensions
      write(200) dims(1:2)
      write(200) nz

      ! Skip atmo
      call fseek(200,Input%s_inv_atmo,1)

      ! Second header
      write(200) Input%nvar
      write(200) dims(3)
      write(200) Sol%omega_input(1:Sol%Num_Wavelength)

      ! For each variable
      do i=1,Input%nvar

        ! If inverting
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! Write positive flag
          write(200) int(1)

        ! Not inverting
        else

          ! Write negative flag
          write(200) int(0)

        end if ! Inverting variable

        ! Number of nodes
        write(200) Inf_Nodes%Num_Nodes(i)

      end do ! Variables

      ! Output RF
      if (Input%Keep_RF) then

        ! Jump result block
        call fseek(200,Input%s_inv_res,1)

        ! Count number of variables in the inversion
        j = 0

        ! For each variable
        do i=1,Input%nvar

          ! In inverting, add one
          if (Inf_Nodes%Nodes_Flags(i)) j = j + 1

        end do ! Variables

        ! Write number of variables in inversion
        write(200) j

        ! Run over variables
        do i=1,Input%nvar

          ! Skip not inverting
          if (.not.Inf_Nodes%Nodes_Flags(i)) cycle

          ! Write index
          write(200) i

          ! Write nodes
          write(200) Inf_Nodes%Num_vary(i)

        end do ! Variables

        ! Skip to last four bytes
        call fseek(200,dims(1)*dims(2)*Input%s_inv_RF_c-4,1)

        ! Write a zero to set the filesize
        write(200) real(0.0)

      ! Not output of RF
      else

        ! Jump result block but 1 position
        call fseek(200,Input%s_inv_res-4,1)

        ! Write a zero to set the filesize
        write(200) real(0.0)

      end if ! Output RF?

      ! Close file
      close(200)

      end subroutine create_io_inv_files

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write the results of the inversion of one pixel into a file\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data
      subroutine Write_Result(Inf_Stokes,Sol,LM_Stru, &
                              Inf_Nodes,Atmo,Bfield,Input)

      ! I/O

      type(Stokes_class), intent(in):: Inf_Stokes
      type(Solution_class), intent(in):: Sol
      type(LMFIT_class), intent(in):: LM_Stru
      type(Nodes_class), intent(in):: Inf_Nodes
      type(Atmo_class), Intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ii,jj,ir,is
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset
      double precision, dimension(:,:), allocatable:: RF
      double precision, dimension(:,:), allocatable:: Stokes_help


      ! Routine name
      urou = 'Write_Result'

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(Input%folder)//'/Result', &
                         MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                         ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Atmosphere
      !

      ! Allocate buffer
      allocate(buffer(Input%s_inv_atmo_c/4))

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(Input%s_inv_atmo_c) + &
                dble(Input%s_inv_h)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      !
      ! Copy into buffer

      ! Initialize
      ii = 0

      ! Tau
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%z)
      ii = ii + Atmo%nz

      ! T
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%T)
      ii = ii + Atmo%nz

      ! Pg
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%Pg)
      ii = ii + Atmo%nz

      ! Bx
      buffer(ii+1:ii+Atmo%nz) = real(Bfield%Bstrength* &
                                     sin(Bfield%Btheta)* &
                                     cos(Bfield%Bphi))
      ii = ii + Atmo%nz

      ! By
      buffer(ii+1:ii+Atmo%nz) = real(Bfield%Bstrength* &
                                     sin(Bfield%Btheta)* &
                                     sin(Bfield%Bphi))
      ii = ii + Atmo%nz

      ! Bz
      buffer(ii+1:ii+Atmo%nz) = real(Bfield%Bstrength* &
                                     cos(Bfield%Btheta))
      ii = ii + Atmo%nz

      ! vx
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%vx*c*1d6)
      ii = ii + Atmo%nz

      ! vy
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%vy*c*1d6)
      ii = ii + Atmo%nz

      ! vz
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%vz*c*1d6)
      ii = ii + Atmo%nz

      ! vmi
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%vmi*c*1d6)
      ii = ii + Atmo%nz

      ! ne
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%ne)
      ii = ii + Atmo%nz

      ! nHt
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%nht)
      ii = ii + Atmo%nz

      ! nHa
      buffer(ii+1:ii+Atmo%nz) = real(Atmo%nha)
      ii = ii + Atmo%nz

      ! nH
      buffer(ii+1:ii+6*Atmo%nz) = real(reshape(Atmo%nh, (/ nz*6 /)))
      ii = ii + 6*Atmo%nz

      ! JKQ?
      if (Input%out_jkqa) then

        ! Init
        jj = 0

        ! J10
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J11R
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J11I
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J20
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J21R
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J21I
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J22R
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz
        jj = jj + Atmo%nz

        ! J22I
        buffer(ii+1:ii+Atmo%nz) = real(Atmo%JKQin(jj+1:jj+Atmo%nz))
        ii = ii + Atmo%nz

      end if ! JKQ

      ! Diffuse light
      buffer(ii+1:ii+1) = real(Atmo%f_diff)
      ii = ii + 1

      ! Write buffer
      call MPI_FILE_WRITE(funit,buffer(1),Input%s_inv_atmo_c/4, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100


      !
      ! Result
      !

      ! Allocate buffer
      deallocate(buffer)
      allocate(buffer(Input%s_inv_res_c/4))

      !
      ! Column offset
      !

      ! Reset offset to zero
      offset = 0
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_SET, ierr)

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(Input%s_inv_res_c) + &
                dble(Input%s_inv_h) + &
                dble(Input%s_inv_atmo) + &
                dble(Input%s_inv_res_h)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      !
      ! Copy into buffer

      ! chi2 original
      buffer(1) = real(LM_Stru%Chisq_0)

      ! chi2 final
      buffer(2) = real(LM_Stru%Chisq)

      ! Current position of buffer
      ii = 2

      ! If polarized
      if (Inf_Nodes%Nodes_Type.ne.0) then

        ! Allocate help
        allocate(Stokes_help(0:3,Sol%Num_Wavelength))

      ! If only intensity
      else

        ! Allocate help
        allocate(Stokes_help(0:0,Sol%Num_Wavelength))

      end if ! If polarization


      !
      ! Observation
      !

      ! Convert the units of I for each range
      do ir=1,Sol%Num_Range

        ! Observation
        Stokes_help(0,Sol%Range(ir,1):Sol%Range(ir,2)) = &
          Inf_Stokes%Stokes_Ob(0,Sol%Range(ir,1):Sol%Range(ir,2))* &
          Sol%Scal_Stokes(ir)*1d-14/c

      end do ! Ranges

      ! Intensity (obs)
      buffer(ii+1:ii+Sol%Num_Wavelength) = &
                             real(Stokes_help(0,1:Sol%Num_Wavelength))
      ii = ii + Sol%Num_Wavelength

      ! Polarization?
      if (Inf_Nodes%Nodes_Type.ne.0) then

        ! If fractional polarization
        if (Sol%Fractional) then

          ! Polarization parameters
          do is=1,3

            ! Convert to absolute the observation
            Stokes_help(is,:) = Inf_Stokes%Stokes_Ob(is,:)* &
                                Stokes_help(0,:)*1d-2

          end do ! Polarization parameters

        ! Non-fractional
        else

          ! Convert the units of I for each range
          do ir=1,Sol%Num_Range

            ! The observation
            Stokes_help(1:3,Sol%Range(ir,1): &
                            Sol%Range(ir,2)) = &
                           Inf_Stokes%Stokes_Ob(1:3,Sol%Range(ir,1): &
                                                   Sol%Range(ir,2))* &
                                           Sol%Scal_Stokes(ir)*1d-14/c

          end do ! Ranges

        end if ! Fractional

        ! -Q
        buffer(ii+1:ii+Sol%Num_Wavelength) = &
                        real(-1d0*Stokes_help(1,1:Sol%Num_Wavelength))
        ii = ii + Sol%Num_Wavelength

        ! -U
        buffer(ii+1:ii+Sol%Num_Wavelength) = &
                        real(-1d0*Stokes_help(2,1:Sol%Num_Wavelength))
        ii = ii + Sol%Num_Wavelength

        ! V
        buffer(ii+1:ii+Sol%Num_Wavelength) = &
                             real(Stokes_help(3,1:Sol%Num_Wavelength))
        ii = ii + Sol%Num_Wavelength

      end if ! Polarization


      !
      ! Fit
      !

      ! Convert the units of I for each range
      do ir=1,Sol%Num_Range

        ! Emergent
        Stokes_help(0,Sol%Range(ir,1):Sol%Range(ir,2)) = &
          Sol%Stokes_out(0,Sol%Range(ir,1):Sol%Range(ir,2))* &
          Sol%Scal_Stokes(ir)*1d-14/c

      end do

      ! Intensity (obs)
      buffer(ii+1:ii+Sol%Num_Wavelength) = &
                             real(Stokes_help(0,1:Sol%Num_Wavelength))
      ii = ii + Sol%Num_Wavelength

      ! Polarization?
      if (Inf_Nodes%Nodes_Type.ne.0) then

        ! If fractional polarization
        if (Sol%Fractional) then

          ! Polarization parameters
          do is=1,3

            ! Convert to absolute the synthesis
            Stokes_help(is,:) = Sol%Stokes_out(is,:)* &
                                Stokes_help(0,:)*1d-2

          end do ! Polarization parameters

        ! Non-fractional
        else

          ! Convert the units of I for each range
          do ir=1,Sol%Num_Range

            ! The synthesis
            Stokes_help(1:3,Sol%Range(ir,1): &
                            Sol%Range(ir,2)) = &
                                 Sol%Stokes_out(1:3,Sol%Range(ir,1): &
                                                   Sol%Range(ir,2))* &
                                           Sol%Scal_Stokes(ir)*1d-14/c

          end do ! Ranges

        end if ! Fractional

        ! -Q
        buffer(ii+1:ii+Sol%Num_Wavelength) = &
                        real(-1d0*Stokes_help(1,1:Sol%Num_Wavelength))
        ii = ii + Sol%Num_Wavelength

        ! -U
        buffer(ii+1:ii+Sol%Num_Wavelength) = &
                        real(-1d0*Stokes_help(2,1:Sol%Num_Wavelength))
        ii = ii + Sol%Num_Wavelength

        ! V
        buffer(ii+1:ii+Sol%Num_Wavelength) = &
                             real(Stokes_help(3,1:Sol%Num_Wavelength))
        ii = ii + Sol%Num_Wavelength

      end if ! Polarization

      ! Deallocate helper
      deallocate(Stokes_help)


      !
      ! Nodes
      !

      ! Run over variables
      do jj=1,Inf_Nodes%nvar

        ! If nodes
        if (Inf_Nodes%Num_nodes(jj).gt.0) then

          ! Positions
          buffer(ii+1:ii+Inf_Nodes%Num_nodes(jj)) = &
                                      real(10**Inf_Nodes%Node(jj)%H)
          ii = ii + Inf_Nodes%Num_nodes(jj)

          !
          ! Values

          ! If velocity
          if (jj.eq.Inf_Nodes%index_vz.or. &
              jj.eq.Inf_Nodes%index_vx.or. &
              (jj.eq.Inf_Nodes%index_vy.and. &
               Inf_Nodes%vtype.eq.1).or. &
              jj.eq.Inf_Nodes%index_vm) then

            ! Transform to km/s
            buffer(ii+1:ii+Inf_Nodes%Num_nodes(jj)) = &
                                  real(Inf_Nodes%Node(jj)%Var*c*1d6)
            ii = ii + Inf_Nodes%Num_nodes(jj)

            ! If inverting
            if (Inf_Nodes%Nodes_Flags(jj)) then

              ! Transform to km/s
              buffer(ii+1:ii+Inf_Nodes%Num_nodes(jj)) = &
                               real(Inf_Nodes%Node(jj)%Errors*c*1d6)
              ii = ii + Inf_Nodes%Num_nodes(jj)

            end if ! Inverting

          ! Others
          else

            ! Copy
            buffer(ii+1:ii+Inf_Nodes%Num_nodes(jj)) = &
                                        real(Inf_Nodes%Node(jj)%Var)
            ii = ii + Inf_Nodes%Num_nodes(jj)

            ! If inverting
            if (Inf_Nodes%Nodes_Flags(jj)) then

              ! Copy
              buffer(ii+1:ii+Inf_Nodes%Num_nodes(jj)) = &
                                     real(Inf_Nodes%Node(jj)%Errors)
              ii = ii + Inf_Nodes%Num_nodes(jj)

            end if ! Inverting
          end if ! Particular variable
        end if ! Nodes?

      end do ! Variables

      ! Write buffer
      call MPI_FILE_WRITE(funit,buffer(1),Input%s_inv_res_c/4, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      !
      ! If keeping RF
      !
      if (Input%Keep_RF) then

        ! Allocate RF
        allocate(RF(Sol%Num_Wavelength,0:3))

        ! Allocate buffer
        deallocate(buffer)
        allocate(buffer(Input%s_inv_RF_c/4))

        !
        ! Column offset
        !

        ! Reset offset to zero
        offset = 0
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_SET, ierr)

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(Input%s_inv_RF_c) + &
                  dble(Input%s_inv_h) + &
                  dble(Input%s_inv_atmo) + &
                  dble(Input%s_inv_res_h) + &
                  dble(Input%s_inv_res) + &
                  dble(Input%s_inv_RF_h)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer offset
        ii = 0

        ! For each variable
        do jj=1,11

          ! Skip not inverting
          if (.not.Inf_Nodes%Nodes_Flags(jj)) cycle

          ! If not varying
          if (Inf_Nodes%Num_vary(jj).le.0) cycle

          ! Add to buffer
          buffer(ii+1:ii+Inf_Nodes%Num_vary(jj)) = &
            real(10d0** &
                 Inf_Nodes%Node(jj)%H(Inf_Nodes%Node_Vary(1,jj): &
                                      Inf_Nodes%Node_Vary(2,jj)))
          ii = ii + Inf_Nodes%Num_vary(jj)

        end do ! Variables

        ! Thermal
        if (Inf_Nodes%Nodes_Type.eq.0) then

          ! For each node
          do jj=1,LM_Stru%Num

            ! Get RF
            RF(:,0) = LM_Stru%JacobianI(:,jj)

            ! For each range in wavelength
            do ir=1,Sol%Num_Range

              ! De-scale RF Stokes parameters
              RF(Sol%Range(ir,1):Sol%Range(ir,2),0) = &
                            RF(Sol%Range(ir,1):Sol%Range(ir,2),0)* &
                            Sol%Scal_Stokes(ir)

            end do ! Wavelength ranges

            ! De-scale parameter scale
            RF(:,0) = RF(:,0)*Inf_Nodes%Scal(Inf_Nodes%Inf_Inv(1,jj))

            ! Add to buffer
            buffer(ii+1:ii+Sol%Num_Wavelength) = real(RF(:,0))
            ii = ii + Sol%Num_Wavelength

          end do ! Nodes

        ! Mag. or all
        else

          ! For each node
          do jj=1,LM_Stru%Num

            ! Get RF
            do is=0,3
              RF(:,is) = LM_Stru%Jacobian(is,:,jj)
            end do

            ! If fractional
            if (Sol%Fractional) then

              ! For each range in wavelength
              do ir=1,Sol%Num_Range

                ! De-scale RF Stokes parameters
                RF(Sol%Range(ir,1):Sol%Range(ir,2),0) = &
                            RF(Sol%Range(ir,1):Sol%Range(ir,2),0)* &
                            Sol%Scal_Stokes(ir)

              end do ! Wavelength ranges

              ! Factor 1e2 for fractional
              RF(:,1:3) = RF(:,1:3)*1d-2

            ! If not fractional
            else

              ! For each Stokes parameter
              do is=0,3

                ! For each range in wavelength
                do ir=1,Sol%Num_Range

                    ! De-scale RF Stokes parameters
                    RF(Sol%Range(ir,1):Sol%Range(ir,2),is) = &
                             RF(Sol%Range(ir,1):Sol%Range(ir,2),is)* &
                             Sol%Scal_Stokes(ir)

                end do ! Wavelength ranges
              end do ! Stokes paramters

            end if ! Fracional?

            ! Scale the response function
            RF = RF*Inf_Nodes%Scal(Inf_Nodes%Inf_Inv(1,jj))

            ! For each Stokes parameter
            do is=0,3

              ! Add to buffer
              buffer(ii+1:ii+Sol%Num_Wavelength) = real(RF(:,is))
              ii = ii + Sol%Num_Wavelength

            end do ! Stokes parameters
          end do ! Nodes

        end if ! Type of inversion

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),Input%s_inv_RF_c/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1100

        ! Deallocate RF
        deallocate(RF)

      end if ! Keep RF

      ! Deallocate buffer
      deallocate(buffer)

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      return

1000  umsg = 'Error opening Result file to write'
      call aborted
      call control
      return
1010  umsg = 'Error seeking Result file'
      call MPI_FILE_CLOSE(funit, ierr)
      call aborted
      call control
      return
1100  umsg = 'Error writing Result file'
      call MPI_FILE_CLOSE(funit, ierr)
      call aborted
      call control
      return

      end subroutine Write_Result

!#####################################################################
!#####################################################################
!#####################################################################

      end module iotic_mod
