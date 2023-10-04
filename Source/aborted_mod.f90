      !> Message printing routine and error handling
      module aborted_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/17/2017
!  Last version:
!     09/08/2023 V3.0.8
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/08/2023:    v3.0.8 - Added the possibility to shut-up the
!                             verbosity output depending on the level
!                             in the INPUT (TdPA)
!                           - main_verbosity now receives an integer
!                             with the message level, instead of the
!                             result of the comparison (TdPA)
!
!     08/17/2023:    v3.0.7 - Control allreduce is performed in place
!                             now (TdPA)
!
!     08/07/2023:    v3.0.6 - Added pixel information to the output
!                             of the error in abortedS (TdPA)
!
!     07/03/2023:    v3.0.5 - Changed call to aborted in the MPI
!                             case to accomodate for parallelization
!                             in pixels (TdPA)
!                           - Cugfix: the parallel approach in
!                             aborteds led to issues in the output
!                             message due to reusing a common
!                             variable (TdPA)
!
!     05/25/2023:    V3.0.4 - Fixed verbosity for abortedS serial
!                             when not crashing (TdPA)
!
!     05/16/2023:    V3.0.3 - Added explicit verbosity in abortedS
!                             for the serial case, because further
!                             control was needed for the inversion
!                             runs (TdPA)
!
!     03/08/2023:    V3.0.2 - Added gabortedv subroutine, which is
!                             the error control called from the
!                             inversion part of the TIC module (TdPA)
!                           - Added verbosev to control the verbosity
!                             in the inversion part of the TIC
!                             module (TdPA)
!                           - Added verboseI to control the verbosity
!                             in the inversion part of the TIC when
!                             it shoudl be decided where to write it
!                             depending on an input (TdPA)
!                           - Added a branch in aborted, gaborted, and
!                             abortedS for the inversion mode (TdPA)
!                           - Added a branch in aborted for the
!                             inversion mode (TdPA)
!
!     07/27/2022:    V3.0.1 - Changed size of input in abortedS for
!                             consistency with global variable (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case I had to do
!                             some dirty edits in this module.
!                              o Now aborted branches to other
!                                routines depending on the running
!                                mode.
!                              o Created gaborted (g == global) which
!                                does what aborted used to do.
!                              o Created gcontrol (g == global) which
!                                does what control used to do.
!                              o Now control only finalizes the run
!                                if in the pure 1D case.
!                              o AbortedS now needs to check on the
!                                new gnproc (g == global) variable
!                                to determine if the run is parallel.
!                                It also needs to call gaborted in the
!                                serial case, instead of aborted.
!                             (TdPA)
!
!     03/23/2021:    V2.0.1 - Added arguments to abortedS to be use it
!                             with OpenMP. This leaft abortedS_w
!                             obsolete and was thus removed (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added correct call to time when using
!                             OpenMP (TdPA)
!
!     02/12/2021:    V1.2.4 - Added report_time, report_mpi_timeI, and
!                             report_mpi_time routines (TdPA)
!
!     09/20/2020:    V1.2.3 - The name of the error file is taken
!                             from a common variable (TdPA)
!
!     05/08/2019:    V1.2.2 - Each CPU writes in standard output that
!                             it had an error just once (TdPA)
!
!     03/13/2019:    V1.2.1 - A CPU only informs once about error
!                             file generation (TdPA)
!                           - Fixed aborted message (TdPA)
!
!     02/20/2019:    V1.2.0 - Completely changed aborted subroutine
!                             and added aborted_silent, abortedS,
!                             verbose, and control routines (TdPA)
!
!     02/14/2019:    V1.1.0 - Made changes to avoid the waiting for
!                             buffer filling when no verbosity (TdPA)
!
!     04/17/2017:    V1.0.0 - First version (TdPA)
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
!    aborted
!      Handles the abortion in case of error depending on the type
!    of run
!
!    gaborted
!      Aborts the code in a controlled way with all CPU
!
!    gabortedv
!      Aborts the code in a controlled way with all CPU. Inversion
!      version
!
!    aborted_silent
!      Aborts the code in a controlled way with all CPU, but in
!      silence
!
!    abortedS
!      Generates error file and flags for abortion
!
!    verbose
!      For message output
!
!    verbosev
!      For message output for the inversion
!
!    verboseI
!      Special verbosity for additional inversion information
!
!    report_time
!      Report CPU time into time file
!
!    report_mpi_timeI
!      Report CPU time called from solveri
!
!    report_mpi_time
!      Report CPU time called from solver.\n
!
!    control
!      Controls if any CPU has crashed and stops if needed depending
!      on the type of run
!
!    gcontrol
!      Controls if any CPU has crashed and stops if needed
!
!#####################################################################
!#####################################################################
!#####################################################################

      use commons_mod
      use omp_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Handles the abortion depending on the case running
      subroutine aborted

      !
      ! If 1D synthesis
      !
      if (run_mode.eq.0) then
        call gaborted
      !
      ! Inversion
      !
      else if (run_mode.eq.-1) then
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)
        call abortedS(umsg,urou,-1,.True.,.True.)
      !
      ! If 1.5D synthesis
      !
      else if (run_mode.eq.1) then
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)
        call abortedS(umsg,urou,-1,.True.,.True.)
      !
      ! If CLE
      !
      else if (run_mode.eq.2) then
        call gaborted
      end if

      end subroutine aborted

!#####################################################################
!#####################################################################
!#####################################################################

      !> Sends a message and waits for every other process before
      !! stopping.
      subroutine gaborted

      !
      ! Output message
      !
      if (pid.eq.0) then

        ! Build message
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)

        ! Verbose the error depending on the running mode
        if (run_mode.eq.-1) then
          call verbosev
        else
          call verbose
        end if
      end if

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      end subroutine gaborted

!#####################################################################
!#####################################################################
!#####################################################################

      !> Sends a message and waits for every other process before
      !! stopping. Inversion version
      subroutine gabortedv

      !
      ! Output message
      !
      if (pid.eq.0) then
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)
        call verbosev
      end if

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      end subroutine gabortedv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Waits for every other process before stopping.
      subroutine aborted_silent

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      end subroutine aborted_silent

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generates error file and returns.\n
      !!      lumsg(character(650)): Error message\n
      !!       lurou(character(20)): Name of calling routine\n
      subroutine abortedS(lumsg,lurou,tid,flag,inform)

      ! I/O
      character(len=500):: lumsg
      character(len=20):: lurou
      logical:: flag,inform
      integer:: tid

      ! Local
      logical:: exists
      character(len=2):: stid
      character(len=500):: cumsg

      !
      ! If serial, just crash
      if (gnproc.eq.1) then


        ! If crashing
        if (flag) then

          ! Dump into globals
          umsg = lumsg
          urou = lurou

          ! Abort
          call gaborted

        end if

        ! Copy message
        cumsg = lumsg

        ! Dump in standard
        if (vaborted.and.inform) then
          umsg = ' ## In routine'//trim(urou)//': '
          ! Call verbose depending on running mode
          if (run_mode.eq.-1) then
            call verbosev
          else
            call verbose
          end if
          ! If 1.5D
          if (run_mode.ne.0) then
            if (icoords(3).gt.0) then
              write(umsg,'(A,1x,i7,1x,"(",i4,",",i4,")")') &
                                  ' For pixel',icoords(3),icoords(1:2)
              ! Call verbose depending on running mode
              if (run_mode.eq.-1) then
                call verbosev
              else
                call verbose
              end if
            end if
          end if
          ! Reset
          umsg = ' Error: '//trim(cumsg)
          ! Call verbose depending on running mode
          if (run_mode.eq.-1) then
            call verbosev
          else
            call verbose
          end if
        end if
      end if ! Serial

      !
      ! If threaded
      if (tid.gt.0) then

        ! Get character for thread
        write(stid,'(I0.2)') tid

        !
        ! ERROR
        inquire(file=trim(errorf)//'-'//trim(stid), exist=exists)
        if(.not.exists)then
          open(800,file=trim(errorf)//'-'//trim(stid))
        else
          open(800,file=trim(errorf)//'-'//trim(stid), &
                   position='append')
        endif

      ! Not threaded
      else

        !
        ! ERROR
        inquire(file=trim(errorf), exist=exists)
        if(.not.exists)then
          open(800,file=trim(errorf))
        else
          open(800,file=trim(errorf),position='append')
        endif

      end if

      !
      ! Inform about error
      if (inform) then

        ! 1D
        if (run_mode.eq.0) then
          write(800,'(A)') ' ## In routine '//trim(lurou)//': '// &
                           trim(lumsg)
        else
          write(800,'(A,1x,i7,1x,"(",i4,",",i4,")",A)') &
                      ' ## In routine '//trim(lurou)// &
                      ' for pixel ',icoords(3),icoords(1:2),': '// &
                      trim(lumsg)
        end if
      end if

      !
      ! Close
      close(800)

      ! And in standard too
      if (vaborted.and.inform) then
        vaborted = .False.
        ! Threaded
        if (tid.gt.0) then
          write(umsg,'(A)') ' ## One of the processes found an '// &
                            'error, check '//trim(errorf)//'-'// &
                            trim(stid)//' file'
        ! Not threaded
        else
          write(umsg,'(A)') ' ## One of the processes found an '// &
                            'error, check '//trim(errorf)//' file'
        end if

        ! Call verbose depending on running mode
        if (run_mode.eq.-1) then
          call verbosev
        else
          call verbose
        end if
      end if

      ! Flag it
      if (flag) laborted = .True.

      ! Return
      return

      end subroutine abortedS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Output messages.\n
      subroutine verbose

      ! Local

      logical:: exists


      !
      ! If output to terminal
      !
      if (verbosity.and.ninv_mode) then

        write(*,'(A)') trim(umsg)

      !
      ! If output to file
      !
      else

        inquire(file=trim(verbosef), exist=exists)
        if(.not.exists)then
          open(900,file=trim(verbosef))
        else
          open(900,file=trim(verbosef),position='append')
        endif
        write(900,'(A)') trim(umsg)
        close(900)

      end if

      end subroutine verbose

!#####################################################################
!#####################################################################
!#####################################################################

      !> Output messages for the inversion.\n
      subroutine verbosev

      ! Local

      logical:: exists


      !
      ! If output to terminal
      !
      if (verbosity) then

        write(*,'(A)') trim(umsg)

      !
      ! If output to file
      !
      else

        inquire(file=trim(verbosefv), exist=exists)
        if(.not.exists)then
          open(900,file=trim(verbosefv))
        else
          open(900,file=trim(verbosefv),position='append')
        endif
        write(900,'(A)') trim(umsg)
        close(900)

      end if

      end subroutine verbosev

!#####################################################################
!#####################################################################
!#####################################################################

      !> Special verbosity for additional inversion information\n
      !!   main_verbosity(integer): Level of the message
      subroutine verboseI(level)

      ! IO
      integer, intent(in):: level

      ! Local
      logical:: exists


      ! Shut-up?
      if (level.gt.slevel) return

      !
      ! If to main file
      !
      if (level.le.vlevel) then

        !
        ! If output to terminal
        !
        if (verbosity) then

          write(*,'(A)') trim(umsg)

        !
        ! If output to file
        !
        else

          inquire(file=trim(verbosefv), exist=exists)
          if(.not.exists)then
            open(900,file=trim(verbosefv))
          else
            open(900,file=trim(verbosefv),position='append')
          endif
          write(900,'(A)') trim(umsg)
          close(900)

        end if

      !
      ! Into extra file
      !
      else

        inquire(file=trim(verbosefv)//'_extra', exist=exists)
        if(.not.exists)then
          open(900,file=trim(verbosefv)//'_extra')
        else
          open(900,file=trim(verbosefv)//'_extra',position='append')
        endif
        write(900,'(A)') trim(umsg)
        close(900)

      end if

      end subroutine verboseI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time into time file.\n
      !!   folder(character(500)): Path to the folder where to write\n
      !!         ID(characteR(9)): ID of the current run\n
      !!          append(logical): Determine if creating the file or
      !!                           writing into an existing one
      subroutine report_time(folder,ID,append)

      ! I/O
      character(len=500), intent(in):: folder
      character(len=9), intent(in):: ID
      logical, intent(in):: append

      ! Local
      real:: tt

      ! Get initial time
      call cpu_time(tt)

      !
      ! If appending
      !
      if (append) then

        open(800,file=trim(folder)//'/times'//ID, &
             position='append')

      !
      ! If not appending
      !
      else

        open(800,file=trim(folder)//'/times'//ID)

      end if

      ! Write time
      write(800,*) tt

      ! And close
      close(800)

      end subroutine report_time

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time called from solveri.\n
      !!   folder(character(500)): Path to the folder where to write\n
      !!         ID(characteR(9)): ID of the current run\n
      !!             cpu(integer): CPU that sent to master\n
      !!            iter(integer): Current SEE iteration\n
      !!           iterr(integer): Current radiation sub-iteration\n
      !!          append(logical): Determine if creating the file or
      !!                           writing into an existing one
      subroutine report_mpi_timeI(folder,ID,cpu,iter,iterr,append)

      ! I/O
      character(len=500), intent(in):: folder
      character(len=9), intent(in):: ID
      logical, intent(in):: append
      integer, intent(in):: cpu, iter,iterr

      ! Local
      double precision:: tt

      ! Get initial time
#ifdef _OPENMP
      tt = omp_get_wtime()
#else
      call cpu_time(tt)
#endif

      !
      ! If appending
      !
      if (append) then

        open(800,file=trim(folder)//'/times_mpiI'//ID, &
             position='append')

      !
      ! If not appending
      !
      else

        open(800,file=trim(folder)//'/times_mpiI'//ID)

      end if

      ! Write time
      write(800,'(i5,1x,i5,1x,i5,1x,es15.8)') iter,iterr,cpu,tt

      ! And close
      close(800)

      end subroutine report_mpi_timeI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time called from solver.\n
      !!   folder(character(500)): Path to the folder where to write\n
      !!         ID(characteR(9)): ID of the current run\n
      !!             cpu(integer): CPU that sent to master\n
      !!            iter(integer): Current SEE iteration\n
      !!          append(logical): Determine if creating the file or
      !!                           writing into an existing one
      subroutine report_mpi_time(folder,ID,cpu,iter,append)

      ! I/O
      character(len=500), intent(in):: folder
      character(len=9), intent(in):: ID
      logical, intent(in):: append
      integer, intent(in):: cpu, iter

      ! Local
      real:: tt

      ! Get initial time
      call cpu_time(tt)

      !
      ! If appending
      !
      if (append) then

        open(800,file=trim(folder)//'/times_mpi'//ID, &
             position='append')

      !
      ! If not appending
      !
      else

        open(800,file=trim(folder)//'/times_mpi'//ID)

      end if

      ! Write time
      write(800,'(i5,1x,i5,1x,es15.8)') iter,cpu,tt

      ! And close
      close(800)

      end subroutine report_mpi_time

!#####################################################################
!#####################################################################
!#####################################################################

      !> Controls if any CPU has crashed.\n
      subroutine control


      call MPI_ALLREDUCE(MPI_IN_PLACE,laborted,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_CTRL,ierr)

      ! If no failure, continue
      if (.not.laborted.or.run_mode.ne.0) return

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      return

      end subroutine control

!#####################################################################
!#####################################################################
!#####################################################################

      !> Controls if any CPU has crashed.\n
      subroutine gcontrol

      ! Local

      logical:: slaborted

      slaborted = laborted

      call MPI_ALLREDUCE(slaborted, laborted, 1, MPI_LOGICAL, &
                         MPI_LOR, MPI_COMM_WORLD, ierr)

      ! If no failure, continue
      if (.not.laborted) return

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      return

      end subroutine gcontrol

!#####################################################################
!#####################################################################
!#####################################################################

      end module aborted_mod
