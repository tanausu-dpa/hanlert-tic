      !> Reading of fudge data
      module rfudge_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     06/27/2022
!  Last version:
!     06/29/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/29/2022:    V3.0.0 - Initial version (TdPA)
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
!  This subroutine reads the fudge data
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

      !> Reads a file with the fudge data.\n
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

        ! Read the fudge data
        if (pid.eq.0) call system('python '//trim(source)// &
                                  'rfudge.py '//trim(filename)// &
                                  ' '//ID//' '//verbosef)

        ! Wait for the master to finish with python
        call MPI_BARRIER(MPI_COMM_RT, ierr)

        ! Open temporal file with fudge data
        open (100,file='tmp_fud_'//ID,status='old',iostat=ios, &
              err=1000)

        ! Success
        read (100,*,err=1100) ios

        ! If no correct file, abort
        if (ios.lt.0) then

          umsg = 'Problem translating the fudge file'
          goto 1200

        end if

        ! How many frequencies with data
        read (100,*,err=1100) fudge%nfreq_f

        ! If there is any data
        if (fudge%nfreq_f.gt.0) then

          ! Allocate to read
          allocate(fudge%fudge_v(fudge%nfreq_f,4))

          ! For each input frequency, read the data
          do ios=1,fudge%nfreq_f
            read (100,*) fudge%fudge_v(ios,:)
          end do

          if (pid.eq.0) then
            umsg = ' - Fudge '//trim(filename)//' read'
            call verbose
          end if

        end if ! fudge data

        ! Cloase fudge file
        close (100)

        ! Control
        call control

        ! Remove temporal file
        if (pid.eq.0) CALL SYSTEM('rm tmp_fud_'//ID)

      ! There is no specified file to read
      else

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
