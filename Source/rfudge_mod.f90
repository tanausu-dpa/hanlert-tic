      !> Reading of fudge data
      module rfudge_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     27/06/2022
!  Last version:
!     17/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     17/12/2024:    V4.0.0 - Updated headers (TdPA)
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
!  rfudge
!    Read a fudge tabulation from the specified file
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read a fudge tabulation from the specified file\n
      !!   filename(character(:)): Name of the file to read\n
      !!     source(character(:)): Path to the source code\n
      !!         ID(character(:)): ID of this run\n
      !!       fudge(fudge_class): Structure with fudge data
      subroutine rFudge(filename,source,ID,fudge)

      ! I/O

      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      type(fudge_class), intent(inout):: fudge

      ! Local

      integer:: ios


      ! Routine name
      urou = 'rFudge'

      ! If there is a fudge factor file
      if (trim(filename).ne.'NONE') then

        ! Master translate the fudge file with python
        if (pid.eq.0) call system('python '//trim(source)// &
                                  'rfudge.py '//trim(filename)// &
                                  ' '//ID//' '//verbosef)

        ! Wait for the master to finish with python
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ! Open temporal file with fudge data
        open (100,file='tmp_fud_'//ID,status='old',iostat=ios, &
              err=1000)

        ! Success
        read (100,*,err=1100) ios

        ! If no correct file
        if (ios.lt.0) then

          ! Issue error
          umsg = 'Problem translating the fudge file'
          goto 1200

        end if ! Wrong file

        ! How many frequencies with data
        read (100,*,err=1100) fudge%nfreq_f

        ! If there is any data
        if (fudge%nfreq_f.gt.0) then

          ! Allocate to read
          allocate(fudge%fudge_v(fudge%nfreq_f,4))
          MRAMc = MRAMc + 1d-6*sizeof(fudge%fudge_v)

          ! For each input frequency
          do ios=1,fudge%nfreq_f

            ! Read the data
            read (100,*) fudge%fudge_v(ios,:)

          end do ! Input frequencies

          ! Master
          if (pid.eq.0) then

            ! Verbose
            umsg = ' - Fudge '//trim(filename)//' read'
            call verbose

          end if ! Master
        end if ! fudge data

        ! Cloase fudge file
        close (100)

        ! Control
        call control

        ! Master remove temporal file
        if (pid.eq.0) CALL SYSTEM('rm tmp_fud_'//ID)

      ! There is no specified file to read
      else

        ! No fudge frequencies
        fudge%nfreq_f = 0

      end if ! There is a fudge file

      return

1000  umsg = 'Error opening fudge file'
      call aborted
      return
1100  umsg = 'Error reading fudge file'
1200  close(100)
      call aborted
      return

      end subroutine rFudge

!#####################################################################
!#####################################################################
!#####################################################################

      end module rfudge_mod
