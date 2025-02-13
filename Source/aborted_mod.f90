      !> Message writing and error handling
      module aborted_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     17/04/2017
!  Last version:
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    v4.0.0 - Removed references to threads in all
!                             routines (TdPA)
!                           - Cleaned code (TdPA)
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
!    aborted
!      Handles the abortion in case of error checked in the RT group
!    depending on the type of run
!
!    gaborted
!      Write error message and terminate because all CPU are in
!    failure status. This is only called by synthesis modules
!
!    gabortedv
!      Write error message and terminate because all CPU are in
!    failure status. This is only called by inversion modules
!
!    aborted_silent
!      Terminate MPI and stop without writing messages
!
!    abortedS
!      Flag failure or inform of an error from the caller CPU
!
!    verbose
!      Message output in synthesis mode
!
!    verbosev
!      Message output in inversion mode
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
!      Report CPU time called from solver
!
!    control
!      Check if any of the CPU in the control (CTRL) group  is in
!    failure status
!
!    gcontrol
!      Check if any of the CPU (WORLD) is in failure status
!
!#####################################################################
!#####################################################################
!#####################################################################

      use commons_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Handles the abortion in case of error checked in the RT group
      !! depending on the type of run
      subroutine aborted

      !
      ! If 1D synthesis
      !
      if (run_mode.eq.0) then

        ! Call aborted with all CPU
        call gaborted

      !
      ! Inversion
      !
      else if (run_mode.eq.-1) then

        ! Called with RT group, write in log and return failure
        ! signal
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)
        call abortedS(umsg,urou,.True.,.True.)

      !
      ! If 1.5D synthesis
      !
      else if (run_mode.eq.1) then

        ! Called with RT group, write in log and return failure
        ! signal
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)
        call abortedS(umsg,urou,.True.,.True.)

      !
      ! If CLE
      !
      else if (run_mode.eq.2) then

        ! Call aborted with all CPU
        call gaborted

      end if ! Type of run

      end subroutine aborted

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write error message and terminate because all CPU are in
      !! failure status. This is only called by synthesis modules
      subroutine gaborted

      !
      ! Output message
      !

      ! Only master
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

      end if ! Master


      !
      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      end subroutine gaborted

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write error message and terminate because all CPU are in
      !! failure status. This is only called by inversion modules
      subroutine gabortedv

      !
      ! Output message
      !

      ! Only Master
      if (pid.eq.0) then

        ! Write message
        umsg = ' ## Controlled abortion called with message: '// &
               trim(umsg(1:416))//' in routine '//urou(1:15)
        call verbosev

      end if ! Master


      !
      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      end subroutine gabortedv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Terminate MPI and stop without writing messages
      subroutine aborted_silent

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      end subroutine aborted_silent

!#####################################################################
!#####################################################################
!#####################################################################

      !! Flag failure or inform of an error from the caller CPU\n
      !!  lumsg(character(650)): Error message\n
      !!   lurou(character(20)): Name of failing routine\n
      !!          flag(logical): If this call is to be considered
      !!                         a termination error\n
      !!        inform(logical): If the error message has to be
      !!                         written in the appropriate buffer
      subroutine abortedS(lumsg,lurou,flag,inform)

      ! I/O

      character(len=500), intent(in):: lumsg
      character(len=20), intent(in):: lurou
      logical, intent(in):: flag,inform

      ! Local

      character(len=500):: cumsg

      logical:: exists


      !
      ! Serial
      !
      if (gnproc.eq.1) then

        ! If crashing
        if (flag) then

          ! Dump message and routine name into the global
          ! variables
          umsg = lumsg
          urou = lurou

          ! Abort right here
          call gaborted

        end if ! Crashing

        ! Copy message
        cumsg = lumsg

        ! If I have not already send a message and I was requested
        ! to inform
        if (vaborted.and.inform) then

          ! Write routine name in message beginning
          umsg = ' ## In routine'//trim(urou)//': '

          !
          ! Call verbose depending on running mode

          ! Inversion
          if (run_mode.eq.-1) then
            call verbosev
          ! Synthesis
          else
            call verbose
          end if

          !
          ! If not 1D synthesis
          if (run_mode.ne.0) then

            ! If non-trivial coordinates
            if (icoords(3).gt.0) then

              ! Write information about the failing pixel
              write(umsg,'(A,1x,i7,1x,"(",i4,",",i4,")")') &
                                  ' For pixel',icoords(3),icoords(1:2)

              !
              ! Call verbose depending on running mode

              ! Inversion
              if (run_mode.eq.-1) then
                call verbosev
              ! Synthesis
              else
                call verbose
              end if ! Inversion/synthesis
            end if ! Non-trivial coordinates
          end if ! Not 1D synthesis

          ! Save the error message in global message buffer
          umsg = ' Error: '//trim(cumsg)

          !
          ! Call verbose depending on running mode

          ! Inversion
          if (run_mode.eq.-1) then
            call verbosev
          ! Synthesis
          else
            call verbose
          end if ! Inversion/synthesis
        end if ! Requested to inform and first error
      end if ! Serial

      !
      ! Inquire if the error file exists
      inquire(file=trim(errorf), exist=exists)

      ! If does not exist
      if(.not.exists)then

        ! Create new file
        open(800,file=trim(errorf))

      ! If it exists
      else

        ! Append to error file
        open(800,file=trim(errorf),position='append')

      endif ! Error file existence


      !
      ! Inform about error
      !
      if (inform) then

        ! 1D synthesis mode
        if (run_mode.eq.0) then

          ! Write the routine name and the error message
          write(800,'(A)') ' ## In routine '//trim(lurou)//': '// &
                           trim(lumsg)

        ! Any other mode
        else

          ! Write the routine name, the failing pixel, and the
          ! error message
          write(800,'(A,1x,i7,1x,"(",i4,",",i4,")",A)') &
                      ' ## In routine '//trim(lurou)// &
                      ' for pixel ',icoords(3),icoords(1:2),': '// &
                      trim(lumsg)

        end if ! 1D or any other mode
      end if ! Informing about the error

      ! Close error file
      close(800)

      !
      ! Write to stdout
      !

      ! If I have not already send a message and I was requested
      ! to inform
      if (vaborted.and.inform) then

        ! This CPU is not allowed to send more notifications
        ! about error to stdout
        vaborted = .False.

        ! Write message in global variable
        write(umsg,'(A)') ' ## One of the processes found an '// &
                          'error, check '//trim(errorf)//' file'

        !
        ! Call verbose depending on running mode

        ! Inversion
        if (run_mode.eq.-1) then
          call verbosev
        ! Synthesis
        else
          call verbose
        end if ! Inversion/synthesis
      end if ! Requested to inform and first error

      ! Flag failure if requested
      if (flag) laborted = .True.

      ! Return
      return

      end subroutine abortedS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Output messages in synthesis mode
      subroutine verbose

      ! Local

      logical:: exists


      ! If I have the permission to write in terminal and this is not
      ! an inversion
      if (verbosity.and.ninv_mode) then

        ! Output to terminal
        write(*,'(A)') trim(umsg)

      ! Cannot write in terminal
      else

        ! Inquire verbosity file existence
        inquire(file=trim(verbosef), exist=exists)

        ! If there is no verbosity file
        if(.not.exists)then

          ! Create a new one
          open(900,file=trim(verbosef))

        ! If there is a verbosity file
        else

          ! Open to append
          open(900,file=trim(verbosef),position='append')

        endif ! Verbosity file existence

        ! Write message in file
        write(900,'(A)') trim(umsg)

        ! Close file
        close(900)

      end if ! Output to terminal or to file

      end subroutine verbose

!#####################################################################
!#####################################################################
!#####################################################################

      !> Message output in inversion mode
      subroutine verbosev

      ! Local

      logical:: exists


      ! If I have the permission to write in terminal
      if (verbosity) then

        ! Output to terminal
        write(*,'(A)') trim(umsg)

      ! Cannot write in terminal
      else

        ! Inquire verbosity file existence
        inquire(file=trim(verbosefv), exist=exists)

        ! If there is no verbosity file
        if(.not.exists)then

          ! Create a new one
          open(900,file=trim(verbosefv))

        ! If there is a verbosity file
        else

          ! Open to append
          open(900,file=trim(verbosefv),position='append')

        endif ! Verbosity file existence

        ! Write message in file
        write(900,'(A)') trim(umsg)

        ! Close file
        close(900)

      end if ! Output to terminal or to file

      end subroutine verbosev

!#####################################################################
!#####################################################################
!#####################################################################

      !> Special verbosity for additional inversion information\n
      !!   level(integer): Level of the issued message
      subroutine verboseI(level)

      ! I/O

      integer, intent(in):: level

      ! Local

      logical:: exists


      ! If the level of the issued message is above the current
      ! shut-up limit, ignore it
      if (level.gt.slevel) return

      !
      ! If the level of the issued message is important enough to go
      ! into the main file
      if (level.le.vlevel) then

        ! If I have the permission to write in terminal
        if (verbosity) then

          ! Output to terminal
          write(*,'(A)') trim(umsg)

        ! Cannot write in terminal
        else

          ! Inquire verbosity file existence
          inquire(file=trim(verbosefv), exist=exists)

          ! If there is no verbosity file
          if(.not.exists)then

            ! Create a new one
            open(900,file=trim(verbosefv))

          ! If there is a verbosity file
          else

            ! Open to append
            open(900,file=trim(verbosefv),position='append')

          endif ! Verbosity file existence

          ! Write message in file
          write(900,'(A)') trim(umsg)

          ! Close file
          close(900)

        end if ! Output to terminal or to file

      !
      ! If the level of the issued message is NOT important enough to
      ! go into the main file
      else

        ! Inquire verbosity file existence
        inquire(file=trim(verbosefv)//'_extra', exist=exists)

        ! If there is no verbosity file
        if(.not.exists)then

          ! Create a new one
          open(900,file=trim(verbosefv)//'_extra')

        ! If there is a verbosity file
        else

          ! Open to append
          open(900,file=trim(verbosefv)//'_extra',position='append')

        endif ! Verbosity file existence

        ! Write message in file
        write(900,'(A)') trim(umsg)

        ! Close file
        close(900)

      end if ! Level of issued message

      end subroutine verboseI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time called from solveri.\n
      !! folder(character(500)): Path to the output folder\n
      !!        append(logical): If appending or creating the file
      subroutine report_ram(folder,append)

      ! I/O

      character(len=500), intent(in):: folder
      logical, intent(in):: append

      ! Local

      character(len=5):: CPUC


      ! Get string CPU ID
      write(CPUC,'(I0.5)') gpid

      !
      ! If appending
      !
      if (append) then

        ! Open existing file to append message
        open(800,file=trim(folder)//'/ram_'//CPUC, &
             position='append')

      !
      ! If not appending
      !
      else

        ! Open new file to write message
        open(800,file=trim(folder)//'/ram_'//CPUC)

        ! Header
        write(800,'(A)') ' | Background '// &
                         ' | Photoioni. '// &
                         ' |  Radiation '// &
                         ' |   Solution '// &
                         ' | Memoizati. '// &
                         ' | Voigt prf. '// &
                         ' | Redistrib. '// &
                         ' | 1st o. PRD '// &
                         ' | Rad. Tran. '// &
                         ' |      Misc. '// &
                         ' |'

      end if ! Appending/writing

      ! Write RAM count
      write(800,'(10(" | ",2x,f9.3)," |",2x,A)') &
            BRAMc,PRAMc,RRAMc, &
            SRAMc,ERAMc,VRAMc, &
            WRAMc,ORAMc,TRAMc, &
            MRAMc,trim(umsg)

      ! And close
      close(800)

      end subroutine report_ram

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time into time file.\n
      !!  folder(character(500)): Path to the output folder\n
      !!        ID(character(9)): ID of the current run\n
      !!         append(logical): If appending or creating the file
      subroutine report_time(folder,ID,append)

      ! I/O

      character(len=500), intent(in):: folder
      character(len=9), intent(in):: ID
      logical, intent(in):: append

      ! Local

      real:: tt


      ! Get time
      call cpu_time(tt)

      !
      ! If appending
      !
      if (append) then

        ! Open existing file to append message
        open(800,file=trim(folder)//'/times'//ID, &
             position='append')

      !
      ! If not appending
      !
      else

        ! Open new file to write message
        open(800,file=trim(folder)//'/times'//ID)

      end if ! Appending/writing

      ! Write time
      write(800,*) tt

      ! And close
      close(800)

      end subroutine report_time

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time called from solveri.\n
      !! folder(character(500)): Path to the output folder\n
      !!       ID(character(9)): ID of the current run\n
      !!           cpu(integer): CPU that sent to master\n
      !!          iter(integer): Current solver iteration\n
      !!         iterr(integer): Current radiation sub-iteration\n
      !!        append(logical): If appending or creating the file
      subroutine report_mpi_timeI(folder,ID,cpu,iter,iterr,append)

      ! I/O

      character(len=500), intent(in):: folder
      character(len=9), intent(in):: ID
      logical, intent(in):: append
      integer, intent(in):: cpu, iter,iterr

      ! Local

      double precision:: tt


      ! Get time
      call cpu_time(tt)

      !
      ! If appending
      !
      if (append) then

        ! Open existing file to append message
        open(800,file=trim(folder)//'/times_mpiI'//ID, &
             position='append')

      !
      ! If not appending
      !
      else

        ! Open new file to write message
        open(800,file=trim(folder)//'/times_mpiI'//ID)

      end if ! Appending/writing

      ! Write iteration and CPU information and time
      write(800,'(i5,1x,i5,1x,i5,1x,es15.8)') iter,iterr,cpu,tt

      ! And close
      close(800)

      end subroutine report_mpi_timeI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report CPU time called from solver\n
      !! folder(character(500)): Path to the output folder\n
      !!       ID(character(9)): ID of the current run\n
      !!           cpu(integer): CPU that sent to master\n
      !!          iter(integer): Current solver iteration\n
      !!        append(logical): If appending or creating the file
      subroutine report_mpi_time(folder,ID,cpu,iter,append)

      ! I/O

      character(len=500), intent(in):: folder
      character(len=9), intent(in):: ID
      logical, intent(in):: append
      integer, intent(in):: cpu, iter

      ! Local

      real:: tt


      ! Get time
      call cpu_time(tt)

      !
      ! If appending
      !
      if (append) then

        ! Open existing file to append message
        open(800,file=trim(folder)//'/times_mpi'//ID, &
             position='append')

      !
      ! If not appending
      !
      else

        ! Open new file to write message
        open(800,file=trim(folder)//'/times_mpi'//ID)

      end if ! Appending/writing

      ! Write iteration and CPU information and time
      write(800,'(i5,1x,i5,1x,es15.8)') iter,cpu,tt

      ! And close
      close(800)

      end subroutine report_mpi_time

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if any of the CPU in the control (CTRL) group  is in
      !! failure status
      subroutine control

      ! Check if any aborted status is activated
      call MPI_ALLREDUCE(MPI_IN_PLACE,laborted,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_CTRL,ierr)

      ! If no failure or if doing anything but 1D synthesis
      if (.not.laborted.or.run_mode.ne.0) return

      ! Only 1D synthesis goes through here

      ! Try to exit MPI
      call MPI_FINALIZE(ierr)

      ! And stop
      stop

      return

      end subroutine control

!#####################################################################
!#####################################################################
!#####################################################################

      !> Controls if any CPU (WORLD) is in failure status
      subroutine gcontrol

      ! Local

      logical:: slaborted


      ! Copy local value
      slaborted = laborted

      ! Check if any aborted status is activated
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
