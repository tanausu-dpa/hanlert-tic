      !> Main
      program hanlert
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!     Hao Li (IAC/NSSCC)
!  Start:
!     17/04/2017
!  Last version:
!     06/06/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/06/2025:    V4.0.3 - Added the size of a new logical in
!                             commons to the misc. counter (TdPA)
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
!    - Improve equation of state (mid priority)
!    - Introduce the possibility of forbidden lines, following always
!      multilevel formalism (low priority).
!    - Generalize frequency axis build to take into account magnetic
!      splitting (maybe irrelevant, low priority)
!    - Implement the modification to the Lambda operator for blended
!      transitions (haven't found problems yet, low priority).
!
!#####################################################################
!#####################################################################
!
!  Main
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use chemicaux_mod
      use fctsg_mod
      use free_mod
      use rfudge_mod
      use gauss_mod
      use hanlert_mod
      use ratom_mod
      use rinput_mod
      use rmol_mod
      use strength_mod
      use tic_mod
      use types_mod

      ! Variables

      type(Atom_class), dimension(:), allocatable:: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(MPI_class):: MPID
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge

      character(len=MPI_MAX_PROCESSOR_NAME) :: pname

      logical:: isPRD, lH

      integer:: ia,iab,i,slen
      integer(kind=MPI_OFFSET_KIND):: offset

      double precision:: t0,t1,day,hour,minute,second


      ! Initialize MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,gpid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,gnproc,ierr)
      MPID%gpid = gpid
      MPID%gnproc = gnproc
      pid = gpid
      nproc = gnproc
      MPID%pid = pid
      MPID%nproc = nproc

      ! Initialize control communicator as world
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, gpid, &
                          MPI_COMM_CTRL, ierr)

      !
      ! Initialize RAM counters
      !
      PRAMc = 0d0
      VRAMc = 0d0
      WRAMc = 0d0
      RRAMc = 0d0
      BRAMc = 0d0
      MRAMc = 0d0
      TRAMc = 0d0
      ORAMc = 0d0
      FRAMc = 0d0
      ERAMc = 0d0
      SRAMc = 0d0
      DRAMc = 0d0
      DRAM2c = 0d0

      !
      ! Add RAM in commons
      !
      MRAMc = 1d-6*dble(500*3 + &
                        650 + &
                        20 + &
                        4*29 + &
                        4 + &
                        4*2 + &
                        4*2 + &
                        4*3 + &
                        4*2 + &
                        4*24 + &
                        4*3 + &
                        4*2 + &
                        4*3 + &
                        8*12)

      !
      ! Add already used RAM
      !

      ! Miscellaneous
      MRAMc = MRAMc + 1d-6*sizeof(MPID) + &
                      1d-6*sizeof(Input) + &
                      1d-6*sizeof(Flgsg) + &
                      1d-6*sizeof(fudge)

      ! Start timer
      if (pid.eq.0) call cpu_time(t0)

      ! Initial verbosity
      verbosity = .True.
      vaborted = .True.
      ninv_mode = .True.


      !
      ! Check if the input file is specified explicitly
      !
      ia = command_argument_count()

      ! There is an argument
      if (ia.ge.1) then

        ! Get command line argument
        call get_command_argument(1, Input%input)

        ! Check argument length
        if (len_trim(Input%input).eq.0) Input%input = 'INPUT'

      ! No argument
      else

        ! Default name
        Input%input = 'INPUT'

      end if ! Command line argument

      ! Check if input file exists
      open(200, file=trim(Input%input), status='old', iostat=ierr)
      close(200)

      ! If failed to open input file
      if (ierr.ne.0) then

        ! Issue error
        umsg = 'Could not find input file with name: '// &
               trim(Input%input)//''
        urou = 'hanlert'
        call gaborted

      end if ! Error


      !
      ! Identify source folder
      !

      ! Try opening cheat file
      open(200, file='spath', status='old', iostat=ierr)

      ! If there is a path file
      if (ierr.eq.0) then

        ! Read path
        read(200,'(A)') Input%source

        ! Close
        close(200)

        ! Try opening source file
        open(200, file=trim(Input%source)//'rinput.py', &
             status='old', iostat=ierr)
        close(200)

        ! Could not find source file
        if (ierr.ne.0) then

          ! Issue error
          umsg = 'Wrong source path specified in spath file'
          urou = 'hanlert'
          call gaborted

        end if ! Error

      ! If there is not a path file try some typical locations
      else

        ! Try finding the source file up to four folders above
        open(200, file='Source/rinput.py', &
             status='old', iostat=ierr)
        close(200)
        if (ierr.eq.0) then
          Input%source = 'Source/'
        else
          open(200, file='../Source/rinput.py', &
               status='old', iostat=ierr)
          close(200)
          if (ierr.eq.0) then
            Input%source = '../Source/'
          else
            open(200, file='../../Source/rinput.py', &
                 status='old', iostat=ierr)
            close(200)
            if (ierr.eq.0) then
              Input%source = '../../Source/'
            else
              open(200, file='../../../Source/rinput.py', &
                   status='old', iostat=ierr)
              close(200)
              if (ierr.eq.0) then
                Input%source = '../../../Source/'
              else
                open(200, file='../../../../Source/rinput.py', &
                     status='old', iostat=ierr)
                close(200)
                if (ierr.eq.0) then
                  Input%source = '../../../../Source/'
                else
                  umsg = 'Source path not specified and '// &
                         'could not find it'
                  urou = 'hanlert'
                  call gaborted
                end if
              end if
            end if
          end if
        end if
      end if


      !
      ! Specify resource folder
      !

      ! Try opening cheat file
      open(200, file='sreso', status='old', iostat=ierr)

      ! If there is a path file
      if (ierr.eq.0) then

        ! Read path
        read(200,'(A)') Input%resource

        ! Close
        close(200)

        ! Try opening source file
        open(200, file=trim(Input%resource)//'partfunc', &
             status='old', iostat=ierr)
        close(200)

        ! Could not find source file
        if (ierr.ne.0) then

          ! Issue error
          umsg = 'Wrong resource path specified in '// &
                 'sreso file'
          urou = 'hanlert'
          call gaborted

        end if ! Error

      ! If there is not a path file try some typical locations
      else

        ! Try finding the resource files up to four folders above
        open(200, file='Resources/partfunc', &
             status='old', iostat=ierr)
        close(200)
        if (ierr.eq.0) then
          Input%resource = 'Resources/'
        else
          open(200, file='../Resources/partfunc', &
               status='old', iostat=ierr)
          close(200)
          if (ierr.eq.0) then
            Input%resource = '../Resources/'
          else
            open(200, file='../../Resources/partfunc', &
                 status='old', iostat=ierr)
            close(200)
            if (ierr.eq.0) then
              Input%resource = '../../Resources/'
            else
              open(200, file='../../../Resources/partfunc', &
                   status='old', iostat=ierr)
              close(200)
              if (ierr.eq.0) then
                Input%resource = '../../../Resources/'
              else
                open(200, file='../../../../Resources/partfunc', &
                     status='old', iostat=ierr)
                close(200)
                if (ierr.eq.0) then
                  Input%resource = '../../../../Resources/'
                else
                  umsg = 'Resource path not specified and '// &
                         'could not find it'
                  urou = 'hanlert'
                  call gaborted
                end if
              end if
            end if
          end if
        end if
      end if


      !
      ! Create run ID
      !

      ! Master
      if (pid.eq.0) then

        ! Get unique ID
        call fgetpid(ia)
        call MPI_GET_PROCESSOR_NAME(pname, slen, ierr)
        do i=1,slen
          ia = ia + ichar(pname(i:i))
        end do

      end if ! Master

      ! Share ID
      call MPI_BCAST(ia, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ! Write string ID
      write(Input%ID,'(i9.9)') ia


      ! Read input data
      call rInput(Input)


      ! If measuring performance, write now
      if (Input%g_perf.and.pid.eq.0) &
        call report_time(Input%folder,Input%ID,.False.)


      !
      ! Read atomic data
      !

      ! Allocate the array of atoms and initialize common variables
      call allocateatom(Atom, NA)
      nxtran = 0 ; nxt = 0 ; nxb = 0 ; nxphot = 0
      nkx = 0 ; nxJ = 0 ; nxS = 0 ; nxL = 0 ; nJs = 0
      isPRD = .False.
      tbAD = .False.

      ! For each active atom
      do ia=1,nA

        ! Read atom with index ia
        call rAtom(Input%atom(ia)%str,Input%source, &
                   Input%ID,Input%skip_wave(ia), &
                   Input%Kcut_input,Atom(ia),ia, &
                   isPRD,.True.)

        ! If not the first atom, add shifts for rolling indexes
        if (ia.gt.1) then
          Atom(ia)%tshift = Atom(ia-1)%tshift + Atom(ia-1)%ntran
          Atom(ia)%tfshift = Atom(ia-1)%tfshift + Atom(ia-1)%nftran
          Atom(ia)%pshift = Atom(ia-1)%pshift + Atom(ia-1)%nphot
        end if

      end do ! Active atoms

      ! Control dimensions, need minimum of 1 to avoid allocations
      ! with zero size
      if (nxb.lt.1) nxb = 1
      if (nxphot.lt.1) nxphot = 1

      ! Update size for signs and factorials
      nxdim = max(4,nxJ*4,nint(dble(nxJ+nxL+nxS)*0.5d0))

      ! Account for LTE lines
      do ia=1,nLTEl

        ! Update size
        nxdim = max(nxdim,nint(4*Input%LTEline(ia)%Ju), &
                          nint(4*Input%LTEline(ia)%Jl))

      end do

      ! Check that there is some PRD line if we chose PRD calculation
      if(.not.isPRD.and.PRD) then
        if (verbosity) then
          umsg = 'PRD is on, but there are not PRD lines'
          urou = 'hanlert'
          call aborted
        else
          umsg = ' ## PRD is on, but there are not PRD lines'
          urou = 'hanlert'
          if (gpid.eq.0) call verbose
          laborted = .True.
          call gcontrol
        end if
      end if

      ! Check if dynamic and angle-averaged
      if(PRD.and.dyn.and.AV.and.pid.eq.0) then
        umsg = ' - Warning: You have selected Angle averaged PRD '// &
               'with velocity fields. An hybrid approach will '// &
               'be used to deal with this'
        call verbose
      end if

      ! Renorm the abundance modifiers
      call abund(Atom)

      ! Unless unmagnetized is already specified, initialize the
      ! Mblock sizes
      if (.not.Input%unmagnetized) call initMblock(Atom)

      ! If there are background atoms to load
      if (Input%nAb.gt.0) then

        ! Allocate Atomb array
        call allocateatom(Atomb, nAb)

        ! For each background atom to read
        do iab=1,Input%nAb

          ! Read atom with index iab
          call rAtom(Input%atomback(iab)%str,Input%source, &
                     Input%ID,.True.,Input%Kcut_input, &
                     Atomb(iab),iab,isPRD,.False.)

        end do ! background atoms to read

      end if ! Background atoms to read

      ! If there were custom K cut inputs
      if (allocated(Input%Kcut_input)) then

        ! Free used Input data if allocated
        MRAMc = MRAMc - 1d-6*sizeof(Input%Kcut_input)
        deallocate(Input%Kcut_input)

        ! Update Kradl
        do ia=1,nA

          ! Get maximum K rad
          Kradl = max(Kradl,maxval(Atom(ia)%Krad))

        end do ! Atoms

      end if ! Custom K cut inputs


      !
      ! Check that H existed in the background and, if not, create
      ! a hard-coded one if its not active neither
      !

      ! Initialize
      lH = .False.

      ! Check between active atoms if there is Hydrogen
      do ia=1,nA
        if (Atom(ia)%Element.eq.' H') then
          lH = .True.
          exit
        end if
      end do

      ! Check between non active atoms if there is Hydrogen
      ! and have not found it yet
      if (.not.lH) then
        if (nAb.gt.0) then
          do ia=1,nAb
            if (Atomb(ia)%Element.eq.' H') then
              lH = .True.
              exit
            end if ! Atom is hydrogen
          end do ! Background atoms
        end if ! There are background atoms
      end if ! Not found Hydrogen yet

      ! Get file labels for atomic stuff
      call set_atom_label(Atom,nA)

      ! If there was no H
      if (.not.lH) then

        ! Shift all background atoms one position to the
        ! right
        call shiftatoms(Atomb,Input%popuback)

        ! Update Input variable (nAb added one in shiftatoms)
        Input%nAb = nAb

        ! Get hard-coded Hydrogen in position 1
        call AtomH(Atomb(1))

      end if ! No Hydrogen model to read

      !
      ! Check there is no repeated atom both in active and
      ! in background
      do ia=1,nA
        do iab=1,nAb
          if (Atom(ia)%Element.eq.Atomb(iab)%Element) then
            umsg = 'You have the same atom in the list '// &
                   'of input and background atoms'
            urou = 'hanlert'
            call aborted
          end if
        end do
      end do

      !
      ! If there are LTE lines
      !
      if (nLTEl.gt.0) then

        ! Run over active atoms
        do ia=1,nA

          ! Run over LTE transitions
          do iab=1,nLTEl

            ! If LTE line from active atom
            if (atom_char2index(Atom(ia)%Element).eq. &
                Input%LTEline(iab)%ele) then

              ! Check if ion for LTE line included in model
              if (minval(Atom(ia)%stage).le. &
                  Input%LTEline(iab)%stage.and. &
                  maxval(Atom(ia)%stage).ge. &
                  Input%LTEline(iab)%stage) then

                ! Abort because the ion is in the model
                write(umsg,'(A,i2,A)') &
                  'You have the same atom in the list of input '// &
                  'and in LTE lines, '//Atom(ia)%Element// &
                  ' ,sharing stages, ',Input%LTEline(iab)%stage, &
                  ', this is not allowed'
                urou = 'hanlert'
                call aborted

              end if ! Ion in the model
            end if ! LTE line from an atom with model

          end do ! LTE lines
        end do ! Active atoms

        ! Run over background atoms
        do ia=1,NAb

          ! Run over LTE transitions
          do iab=1,nLTEl

            ! If LTE line from background atom
            if (atom_char2index(Atomb(ia)%Element).eq. &
                Input%LTEline(iab)%ele) then

              ! Save the model for LTE line
              Input%LTEline(iab)%is_passive = .True.
              Input%LTEline(iab)%ia = ia

              ! If the LTE line stage is in the model
              if (minval(Atomb(ia)%stage).le. &
                  Input%LTEline(iab)%stage.and. &
                  maxval(Atomb(ia)%stage).ge. &
                  Input%LTEline(iab)%stage) &
                ! Remove the LTE transition from the
                ! atomic model
                call remove_LTE_transition(Atomb(ia), &
                                           Input%LTEline(iab))

            end if ! LTE line in background model atom

          end do ! LTE lines
        end do ! Passive atoms

        ! For each LTE line
        do ia=1,nLTEl

          ! Setup atomic quantities
          call setup_LTE_transition(Atomb,Input%LTEline(ia))

        end do ! LTE lines

      end if ! If there are LTE lines


      !
      ! Molecules
      !

      ! Allocate molecules array
      call allocatemol(Mol,Input%nM)

      ! If there are molecules to read
      if (Input%nM.gt.0) then

        ! For each molecule to read
        do ia=1,Input%nM

          ! Read molecule model
          call rMol(Input%mol(ia)%str,Input%source,Input%ID,Mol(ia))

        end do ! Molecules

      end if ! Molecules to read


      !
      ! Initialize Racah algebra
      !

      ! Tell structure if doing memoization
      Flgsg%memo = Input%memo

      ! Initialize factorials and signs
      call fctsg(Flgsg)

      ! Verbose
      if (pid.eq.0) then
        umsg = ' - Factorials and signs initialized'
        call verbose
      end if


      !
      ! Initialize FS Einstein coefficients
      !

      ! For every active atom
      do ia=1,nA

        ! Setup Einstein coefficients for FS transitions
        call setFScoeff(Atom(ia),Flgsg)

      end do ! Active atoms

      ! For every background atom
      do ia=1,NAb

        ! Skip hard-coded Hydrogen
        if (.not.lH.and.Atomb(ia)%Element.eq.' H') cycle

        ! Setup Einstein coefficients for FS transitions
        call setFScoeff(Atomb(ia),Flgsg)

      end do ! Background atoms


      !
      ! Calculate dipole strengths
      !

      ! For every active atom
      do ia=1,nA

        ! Call strength routine
        call strength(Atom(ia),Flgsg)

      end do ! Active atoms

      ! Verbose
      if (pid.eq.0.and.nA.gt.0) then
        umsg = ' - Dipole strengths calculated'
        call verbose
      end if


      !
      ! Read fudge data
      call rFudge(Input%fudge,Input%source,Input%ID,fudge)


      !
      ! Tell the atoms if their populations are fixed
      !

      ! For each active atom
      do ia=1,nA

        ! Copy from input
        Atom(ia)%fixp = Input%fixp(ia)
        Atom(ia)%fixplt = Input%fixplt(ia)
        Atom(ia)%zero_ion = Input%zero_ion(ia)

      end do ! Active atoms

      ! Deallocate fixp, fixplt, zero_ion, and filenames
      if (allocated(Input%fixp)) &
        deallocate(Input%fixp,Input%zero_ion,Input%atom,Input%fixplt)
      if (allocated(Input%atomback)) deallocate(Input%atomback)
      if (allocated(Input%mol)) deallocate(Input%mol)

      ! Maximum value of offset MPI kind
      offlimit = 2d0**(MIN(kind(ia),kind(offset))*8d0 - 1d0) - 1d0


      !
      ! Branch the code here for different running modes
      !

      ! Inversion mode
      if (run_mode.eq.-1) then

        ! Verbose inversion mode
        umsg = ' - Inversion mode'
        if (gpid.eq.0) call verbose

         ! Call the Tenerife Inversion Code
        call TIC(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! 1D Synthesis mode
      else if (run_mode.eq.0) then

        ! Verbose 1D synthesis mode
        umsg = ' - 1D synthesis mode'
        if (gpid.eq.0) call verbose

        ! Call the manager for a single 1D synthesis
        call HanleRT1DS(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! 1.5D Sythesis mode
      else if (run_mode.eq.1) then

        ! Verbose 1.5D synthesis mode
        umsg = ' - 1.5D synthesis mode'
        if (gpid.eq.0) call verbose

        ! Call the manager for 1.5D synthesis in a 3D model
        call HanleRT15DS(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! CLE Sythesis mode
      else if (run_mode.eq.2) then

        ! Verbose CLE synthesis mode
        umsg = ' - CLE synthesis mode'
        if (gpid.eq.0) call verbose

        ! Call the manager for the coronal line emission
        call HanleCLE(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      end if

      !
      ! Cleaning up
      !
      call free_atom_full(Atom)
      call free_atom_full(Atomb)
      call free_mol_full(Mol)


      ! Global master
      if (gpid.eq.0) then

        ! Final timer
        call cpu_time(t1)

        ! Initialize time components
        second = t1 - t0
        minute = 0d0
        hour = 0d0
        day = 0d0

        ! Transform overflowing second to minutes
        if (second.ge.60.0) then
          minute = real(floor(second/60.))
          second = second - minute*60.0
        end if

        ! Transform overflowing minutes to hours
        if (minute.ge.60.0) then
          hour = real(floor(minute/60.))
          minute = minute - hour*60.0
        end if

        ! Transform overflowing hours to days
        if (hour.ge.24.) then
          day = real(floor(hour/24.))
          hour = hour - day*24.
        end if

        !
        ! Program finished
        !

        ! Verbose finished
        umsg = ' - Hanle+RT finished without '// &
               'technical problems'
        call verbose

        !
        ! Verbose time elapsed
        !

        ! If took days
        if (nint(day).gt.0) then

          ! Write d-h-m-s
          write(umsg,'(" - Time: ",i2," days, ",i2," hours, "'// &
                    ',i2," minutes, and ",i2," seconds")') &
                nint(day),nint(hour),nint(minute),nint(second)

        ! If took hours
        else if (nint(hour).gt.0) then

          ! Write h-m-s
          write(umsg,'(" - Time: ",i2," hours, ",i2,'// &
                    '" minutes, and ",i2," seconds")') &
                nint(hour),nint(minute),nint(second)

        ! If took minutes
        else if (nint(minute).gt.0) then

          ! Write m-s
          write(umsg,'(" - Time: ",i2," minutes and "'// &
                    ',i2," seconds")') nint(minute),nint(second)

        ! If took seconds
        else

          ! Write s
          write(umsg,'(" - Time: ",i2," seconds")') nint(second)

        end if

        ! Output time elapsed
        call verbose

        !
        ! If measuring performance
        !
        if (Input%g_perf) &
          call report_time(Input%folder,Input%ID,.True.)

      endif ! Master

      !
      ! Finalize MPI
      !
      call MPI_FINALIZE(ierr)

      end program hanlert
