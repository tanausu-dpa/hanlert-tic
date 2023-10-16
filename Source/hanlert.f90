      !> Main
      program hanlert
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!     Hao Li (IAC)
!  Contributors:
!     Ricky Egeland (HAO)
!  Start:
!     04/17/2017
!  Last version:
!     09/29/2023 V3.0.12
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:   V3.0.12 - Updated calls to rAtom (TdPA)
!                           - Added cleaning of the array
!                             Input%Kcut_input (TdPA)
!
!     09/25/2023:   V3.0.11 - Added call to set_atom_label (TdPA)
!
!     08/17/2023:   V3.0.10 - Improved error message when a LTE line
!                             is in an ion with model atom (TdPA)
!
!     08/07/2023:    V3.0.9 - Added logic to deal with LTE lines that
!                             exist in model atoms (TdPA)
!
!     07/03/2023:    V3.0.8 - Updated the call to rAtom with the
!                             Input%skip_wave flag (TdPA)
!                           - Added the parsing of Input%fixplt to
!                             Atom%fixplt (TdPA)
!
!     03/23/2023:    V3.0.7 - Fixed call to MPI_INIT with OpenMP with
!                             newest MPI standards (TdPA)
!
!     03/08/2023:    V3.0.6 - Added TIC module (TdPA)
!
!     02/14/2023:    V3.0.5 - Setting up the geometry is now full
!                             responsibility for each branch (TdPA)
!
!     11/10/2022:    V3.0.4 - Added logic to pass zero_ion from
!                             the Input structure to the Atom
!                             structure (TdPA)
!
!     10/25/2022:    V3.0.3 - Added branch for CLE (TdPA)
!                           - Added some memory cleaning (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: It is necessary to initialize
!                             the M blocks in the atoms if there is
!                             potential of having a non-zero
!                             magnetic field (TdPA)
!                           - Bugfix: The timer can only be stopped by
!                             the global Master (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o The old nproc and pid have been
!                                changed to gnproc and gpid (g ==
!                                global).
!                              o A new MPI_COMM_CTRL communicator
!                                has been defined, but not sure
!                                how necessary it is.
!                              o Changed the calls to aborted into
!                                calls to gaborted.
!                              o Atmosphere and magnetic field are
!                                read elsewhere now.
!                              o The flag KSTK is defined elsewhere.
!                              o rMol does not requiere the argument
!                                Atmo anymore.
!                              o Gauss does not requiere the argument
!                                Bfield anymore.
!                              o check_axial is called elsewhere now.
!                              o The Hamiltonians are diagonalized
!                                elsewhere now.
!                              o The fudge factors and the Kurucz
!                                lines are read here just once.
!                              o Atoms are flagged to fix their
!                                populations here now.
!                              o New routines are called to perform
!                                the calculations depending on the
!                                type of run.
!                             (TdPA)
!                           - Fixed wrong indentation for the section
!                             that creates the run ID (TdPA)
!
!     09/30/2021:    V2.0.1 - nthread asked from a single CPU (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Change MPI initialization when compiling
!                             with OpenMP (TdPA)
!                           - Fix number of OpenMP threads and make
!                             it static (TdPA)
!                           - Changed time calls depending on OpenMP
!                             compilation (TdPA)
!                           - Moved check of axial compliance into
!                             a module (TdPA)
!
!     02/12/2021:   V1.2.15 - Added calls to report_time (TdPA)
!
!     01/12/2021:   V1.2.14 - Added definition of KSTK (TdPA)
!
!     09/28/2020:   V1.2.13 - Added a new method to get an unique
!                             job ID that does not collide for
!                             fast sends of jobs (RE)
!
!     09/11/2020:   V1.2.12 - The search for the Source and
!                             Resources folders reaches one more step
!                             up in the folder hierarchy (TdPA)
!
!     07/22/2020:   V1.2.11 - Changed initialization of nxb and
!                             nxphot. There was wasted memory in the
!                             J00P arrays. Added a control on the
!                             dimensions instead (TdPA)
!
!     03/05/2020:   V1.2.10 - Passing another parameter to shiftatoms
!                             routine (TdPA)
!                           - Not passing lH to hanle routine (TdPA)
!
!     12/17/2019:    V1.2.9 - Passing a new parameter to diagon (TdPA)
!
!     12/11/2019:    V1.2.8 - Communication between Input and Flgsg
!                             regarding memoization (TdPA)
!
!     09/26/2019:    V1.2.7 - Changes in calls due to changes inside
!                             those modules (TdPA)
!
!     09/13/2019:    V1.2.6 - Changed aborting message when asking for
!                             PRD angle-averaged with dynamics,
!                             because now it is allowed (TdPA)
!
!     08/08/2019:    V1.2.5 - Added argument to setFScoeff (TdPA)
!
!     06/12/2019:    V1.2.4 - Fixed two abortion messages with wrong
!                             initial 'tick' (TdPA)
!
!     06/03/2019:    V1.2.3 - Molecule allocation moved outside of
!                             conditional block (TdPA)
!
!     05/08/2019:    V1.2.2 - Added new needed initializations (TdPA)
!
!     03/18/2019:    V1.2.1 - Bugfix: Wrong logic when checking file
!                             existence (TdPA)
!
!     02/20/2019:    V1.2.0 - New verbosity (TdPA)
!                           - Using unit 200 now (TdPA)
!
!     09/27/2018:    V1.1.1 - Added a level higher when checking the
!                             location of the source and resources due
!                             to my own convenience (TdPA)
!
!     08/09/2018:    V1.1.0 - Allows a command line argument to
!                             specify explicitly the name of the input
!                             file (TdPA)
!
!     08/08/2018:    V1.0.8 - Bugfix: time variables were not
!                             initialized if the time was shorter than
!                             the unit of the variable (TdPA)
!
!     08/06/2018:    V1.0.7 - Outputs invested time (TdPA)
!
!     08/03/2018:    V1.0.6 - Removed warning for Kcut and PRD (TdPA)
!
!     07/27/2018:    V1.0.5 - Added warning for Kcut and PRD (TdPA)
!
!     05/17/2018:    V1.0.4 - Added a check for the vectorial
!                             quantities when axial symmetry is
!                             explicit (TdPA)
!
!     09/15/2017:    V1.0.3 - Added run ID (TdPA)
!                           - Added search of Source and Resources
!                             folder (TdPA)
!
!     07/14/2017:    V1.0.2 - Bugfix: one loop that was suppose to
!                             go through the background atoms was
!                             going through the active atoms (TdPA)
!                           - Bugfix: Initialize background atom
!                             FS Einstein coefficients (TdPA)
!
!     07/05/2017:    V1.0.1 - Checks if AV PRD and dynamics to cancel
!                             the calculation (TdPA)
!
!     04/17/2017:    V1.0.0 - First version (TdPA)
!
!#####################################################################
!#####################################################################
!
!  TODO:
!
!    - Generalize frequency axis build to take into account magnetic
!      splitting.
!    - Implement the modification to the Lambda operator for blended
!      transitions.
!    - (Maybe) introduce the possibility of forbidden lines, following
!      always multilevel formalism.
!    - Implement a better alternative solveri and solver that splits
!      frequencies within each CPU.
!    - Optimize the input frequency axis in such a way that one avoids
!      computing zero redistribution profiles
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
!    Main routine
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use chemicaux_mod
      use fctsg_mod
      use rfudge_mod
      use gauss_mod
      use hanlert_mod
      use omp_mod
      use ratom_mod
      use rinput_mod
      use rmol_mod
      use strength_mod
      use tic_mod
      use types_mod

      ! Variables

      logical:: isPRD, lH

      integer:: ia, iab

      type(Atom_class), dimension(:), allocatable:: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(MPI_class):: MPID
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge

      character(len=MPI_MAX_PROCESSOR_NAME) :: pname

      integer :: i, slen

      integer(kind=MPI_OFFSET_KIND):: offset

      double precision:: t0, t1
      double precision:: day, hour, minute, second


      !
      ! Initialize MPI
      !
#ifdef _OPENMP
#ifdef oldmpi
      call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED)
#else
      call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,ierr)
#endif
#else
      call MPI_INIT(ierr)
#endif
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
      ! Initialize OpenMP num. threads
      !
      nthread = 1
#ifdef _OPENMP
!$omp parallel shared(nthread)
      call omp_set_dynamic(.False.)
!$omp single
      nthread = omp_get_num_threads()
!$omp end single
!$omp end parallel
#endif

      ! Start timer
#ifdef _OPENMP
      t0 = omp_get_wtime()
#else
      if (pid.eq.0) call cpu_time(t0)
#endif

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

        call get_command_argument(1, Input%input)

        ! Check argument length
        if (len_trim(Input%input).eq.0) Input%input = 'INPUT'

      ! Default name
      else

        Input%input = 'INPUT'

      end if

      ! Check if input file exists
      open(200, file=trim(Input%input), status='old', iostat=ierr)
      close(200)
      if (ierr.ne.0) then
        umsg = 'Could not find input file with name: '// &
               trim(Input%input)//''
        urou = 'hanlert'
        call gaborted
      end if


      !
      ! Identify source folder
      !
      open(200, file='spath', status='old', iostat=ierr)
      ! If there is a path file
      if (ierr.eq.0) then
        read(200,'(A)') Input%source
        close(200)
        open(200, file=trim(Input%source)//'rinput.py', &
             status='old', iostat=ierr)
        close(200)
        if (ierr.ne.0) then
          umsg = 'Wrong source path specified in spath file'
          urou = 'hanlert'
          call gaborted
        end if
      ! If there is not a path file try some typical locations
      else
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
      open(200, file='sreso', status='old', iostat=ierr)
      ! If there is a path file
      if (ierr.eq.0) then
        read(200,'(A)') Input%resource
        close(200)
        open(200, file=trim(Input%resource)//'partfunc', &
             status='old', iostat=ierr)
        close(200)
        if (ierr.ne.0) then
          umsg = 'Wrong resource path specified in '// &
                 'sreso file'
          urou = 'hanlert'
          call gaborted
        end if
      ! If there is not a path file try some typical locations
      else
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
      if (pid.eq.0) then
        call fgetpid(ia)
        call MPI_GET_PROCESSOR_NAME(pname, slen, ierr)

        do i = 1, slen
          ia = ia + ichar(pname(i:i))
        end do
      end if
      call MPI_BCAST(ia, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      write(Input%ID,'(i9.9)') ia


      !
      ! Read input data
      !
      call rInput(Input)


      !
      ! If measuring performance
      !
      if (Input%g_perf.and.pid.eq.0) &
        call report_time(Input%folder,Input%ID,.False.)


      !
      ! Read atomic data
      !

      ! Allocate the array of atoms and initialize common variables
      call allocateatom(Atom, NA)
      nxtran = 0 ; nxt = 0 ; nxb = 0 ; nxphot = 0
      nkx = 0 ; nxJ = 0 ; nxS = 0 ; nxL = 0
      isPRD = .False.
      tbAD = .False.

      ! For each active atom
      do ia=1,nA
        call rAtom(Input%atom(ia)%str,Input%source, &
                   Input%ID,Input%skip_wave(ia), &
                   Input%Kcut_input,Atom(ia),ia, &
                   isPRD,.True.)

        ! If not the first atom, add shifts
        if (ia.gt.1) then
          Atom(ia)%tshift = Atom(ia-1)%tshift + Atom(ia-1)%ntran
          Atom(ia)%tfshift = Atom(ia-1)%tfshift + Atom(ia-1)%nftran
          Atom(ia)%pshift = Atom(ia-1)%pshift + Atom(ia-1)%nphot
        end if

      end do

      ! Control dimensions
      if (nxb.lt.1) nxb = 1
      if (nxphot.lt.1) nxphot = 1

      ! Update size for signs and factorials
      nxdim = max(nxJ*4,nint(dble(nxJ+nxL+nxS)*0.5d0))

      ! Check that there is some PRD line if we chose PRD calculation
      if(.not.isPRD.and.PRD) then
        umsg = 'PRD is on, but there are not PRD lines'
        urou = 'hanlert'
        call aborted
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

      ! Read atomic background data
      if (Input%nAb.gt.0) then
        ! Call to allocate
        call allocateatom(Atomb, nAb)
        do iab=1,Input%nAb
          call rAtom(Input%atomback(iab)%str,Input%source, &
                     Input%ID,.True.,Input%Kcut_input, &
                     Atomb(iab),iab,isPRD,.False.)
        end do
      end if

      ! If there were custom K cut inputs
      if (allocated(Input%Kcut_input)) then

        ! Free used Input data if allocated
        deallocate(Input%Kcut_input)

        ! Update Kradl
        do ia=1,nA

          ! Get maximum K rad
          Kradl = max(Kradl,maxval(Atom(ia)%Krad))

        end do

      end if ! Custom K cut inputs

      !
      ! Check that H existed in the background and, if not, create
      ! a hardwired one if its not active neither
      !

      ! Initialize
      lH = .False.

      ! Check between active atoms
      do ia=1,nA
        if (Atom(ia)%Element.eq.' H') then
          lH = .True.
          exit
        end if
      end do

      ! Check between non active atoms
      if (nAb.gt.0) then
        do ia=1,nAb
          if (Atomb(ia)%Element.eq.' H') then
            lH = .True.
            exit
          end if
        end do
      end if

      ! Get file labels for atomic stuff
      call set_atom_label(Atom,nA)

      ! If there was no H, add H to the beginning
      if (.not.lH) then
        call shiftatoms(Atomb,Input%popuback)
        Input%nAb = nAb
        call AtomH(Atomb(1))
      end if

      !
      ! Check there is no atom both in calculation and input
      do ia=1,NA
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

        ! Check if active as well
        do ia=1,NA
          do iab=1,nLTEl

            ! Check element
            if (atom_char2index(Atom(ia)%Element).eq. &
                Input%LTEline(iab)%ele) then

              ! Check stage
              if (minval(Atom(ia)%stage).le. &
                  Input%LTEline(iab)%stage.and. &
                  maxval(Atom(ia)%stage).ge. &
                  Input%LTEline(iab)%stage) then

                write(umsg,'(A,i2,A)') &
                  'You have the same atom in the list of input '// &
                  'and in LTE lines, '//Atom(ia)%Element// &
                  ' ,sharing stages, ',Input%LTEline(iab)%stage, &
                  ', this is not allowed'
                urou = 'hanlert'
                call aborted

              end if
            end if

          end do ! LTE lines
        end do ! Active atoms

        ! Check if passive as well
        do ia=1,NAb
          do iab=1,nLTEl

            ! If the passive atom has LTE line
            if (atom_char2index(Atomb(ia)%Element).eq. &
                Input%LTEline(iab)%ele) then

              ! Save the model for LTE line
              Input%LTEline(iab)%is_passive = .True.
              Input%LTEline(iab)%ia = ia

              ! Remove the transition in the atom, if present and in
              ! same stage
              if (minval(Atomb(ia)%stage).le. &
                  Input%LTEline(iab)%stage.and. &
                  maxval(Atomb(ia)%stage).ge. &
                  Input%LTEline(iab)%stage) &
                call remove_LTE_transition(Atomb(ia), &
                                           Input%LTEline(iab))

            end if

          end do ! LTE lines
        end do ! Passive atoms

        ! For each LTE line, get atomic quantities
        do ia=1,nLTEl
          call setup_LTE_transition(Atomb,Input%LTEline(ia))
        end do

      end if ! If there are LTE lines


      !
      ! Read Molecules
      !
      call allocatemol(Mol,Input%nM)
      if (Input%nM.gt.0) then
        do ia=1,Input%nM
          call rMol(Input%mol(ia)%str,Input%source,Input%ID,Mol(ia))
        end do
      end if


      !
      ! Initialize Racah algebra
      !
      ! Tell structure if doing memoization
      Flgsg%memo = Input%memo
      call fctsg(Flgsg)
      if(pid.eq.0) then
        umsg = ' - Factorials and signs initialized'
        call verbose
      end if


      !
      ! Initialize FS Einstein coefficients
      !
      do ia=1,NA
        call setFScoeff(Atom(ia),Flgsg)
      end do
      do ia=1,NAb
        if (.not.lH.and.Atomb(ia)%Element.eq.' H') cycle
        call setFScoeff(Atomb(ia),Flgsg)
      end do


      !
      ! Calculate dipole strengths
      !
      do ia=1,nA
        call strength(Atom(ia),Flgsg)
      end do
      if(pid.eq.0) then
        umsg = ' - Dipole strengths calculated'
        call verbose
      end if


      !
      ! Read fudge data
      !
      call rFudge(Input%fudge,Input%source,Input%ID,fudge)

      !
      ! Tell the atoms if their populations are fixed
      !

      ! For each atom
      do ia=1,na

        ! Copy from input
        Atom(ia)%fixp = Input%fixp(ia)
        Atom(ia)%fixplt = Input%fixplt(ia)
        Atom(ia)%zero_ion = Input%zero_ion(ia)

      end do ! Atoms

      ! Deallocate fixp, fixplt, zero_ion, and filenames
      deallocate(Input%fixp,Input%zero_ion,Input%atom,Input%fixplt)

      ! Maximum value of offset MPI kind
      offlimit = 2d0**(MIN(kind(ia),kind(offset))*8d0 - 1d0) - 1d0

      !
      ! Branch the code here for different running modes
      !

      ! Inversion mode
      if (run_mode.eq.-1) then

        umsg = ' - Inversion mode'
        if (gpid.eq.0) call verbose

        call TIC(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! 1D Synthesis mode
      else if (run_mode.eq.0) then

        umsg = ' - 1D synthesis mode'
        if (gpid.eq.0) call verbose

        call HanleRT1DS(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! 1.5D Sythesis mode
      else if (run_mode.eq.1) then

        umsg = ' - 1.5D synthesis mode'
        if (gpid.eq.0) call verbose

        call HanleRT15DS(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      ! CLE Sythesis mode
      else if (run_mode.eq.2) then

        umsg = ' - CLE synthesis mode'
        if (gpid.eq.0) call verbose

        call HanleCLE(Input,Atom,Atomb,Mol,Flgsg,fudge,MPID)

      end if

      !
      ! Cleaning up
      !
      if (allocated(Atom)) deallocate(Atom)
      if (allocated(Atomb)) deallocate(Atomb)
      if (allocated(Mol)) deallocate(Mol)


      ! Stop timer, handle time, and finish
      if (gpid.eq.0) then

#ifdef _OPENMP
        t1 = omp_get_wtime()
#else
        call cpu_time(t1)
#endif

        second = t1 - t0
        minute = 0d0
        hour = 0d0
        day = 0d0

        if (second.ge.60.0) then
          minute = real(floor(second/60.))
          second = second - minute*60.0
        end if

        if (minute.ge.60.0) then
          hour = real(floor(minute/60.))
          minute = minute - hour*60.0
        end if

        if (hour.ge.24.) then
          day = real(floor(hour/24.))
          hour = hour - day*24.
        end if


        !
        ! Program finished
        !

        umsg = ' - Hanle+RT finished without '// &
               'technical problems'
        call verbose

        ! Output time it took
        if (nint(day).gt.0) then

          write(umsg,'(" - Time: ",i2," days, ",i2," hours, "'// &
                    ',i2," minutes, and ",i2," seconds")') &
                nint(day),nint(hour),nint(minute),nint(second)

        else if (nint(hour).gt.0) then

          write(umsg,'(" - Time: ",i2," hours, ",i2,'// &
                    '" minutes, and ",i2," seconds")') &
                nint(hour),nint(minute),nint(second)

        else if (nint(minute).gt.0) then

          write(umsg,'(" - Time: ",i2," minutes and "'// &
                    ',i2," seconds")') nint(minute),nint(second)

        else

          write(umsg,'(" - Time: ",i2," seconds")') nint(second)

        end if

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
