      !> Distribution of MPI tasks
      module setmpi_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/19/2017
!  Last version:
!     02/19/2024 V3.0.14
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     02/19/2024:   V3.0.14 - Bugfix: The wrong geometrical structure
!                             was being used to re-evaluate the size
!                             sizei8 (TdPA)
!
!     01/09/2024:   V3.0.13 - Bugfix: There was a termination
!                             condition when splitting 1.5D tasks that
!                             was too strict (TdPA)
!
!     09/29/2023:   V3.0.12 - Only output MPI_info files by user
!                             request (TdPA)
!
!     08/22/2023:   V3.0.11 - Specify the running mode in a warning
!                             message (TdPA)
!
!     07/03/2023:   V3.0.10 - Bugfix: Typo in warning message (TdPA)
!
!     03/21/2023:    V3.0.9 - Bugfix: When using a single CPU for
!                             1.5D, or selecting a size of RT groups
!                             larger or equal than the number of
!                             available slaves, the implementation
!                             was wrong (TdPA)
!
!     02/14/2023:    V3.0.8 - Split Geom and GeomI in setmpi_sizes
!                             routine (TdPA)
!                           - Split AV and AVI (TdPA)
!
!     11/24/2022:    V3.0.7 - Added setmpi_CLE (TdPA)
!
!     10/26/2022:    V3.0.6 - Changed the storage structure of the
!                             rdip variable (TdPA)
!
!     10/25/2022:    V3.0.5 - Added argument to setmpi_sizes, which
!                             now can re-adjust the sizes which
!                             depend on the dimension of the vertical
!                             axis (TdPA)
!                           - Implemented restriction of the height
!                             axis (TdPA)
!                           - Implemented in setomp_magn the
!                             management of repeated calls (TdPA)
!
!     07/27/2022:    V3.0.4 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!
!     07/18/2022:    V3.0.3 - Added tge adjust_IW subroutine (TdPA)
!
!     07/13/2022:    V3.0.2 - Changed intent of Input argument from in
!                             to inout in setmpi15D (TdPA)
!                           - Added some control in case the user ask
!                             for more RT processes than available
!                             in setmpi15D (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: Always initialize MPI%pid,
!                             MPI%nproc, and MPI%mpi15d (TdPA)
!                           - Bugfix: The intent of MPI must be
!                             inout, and not out, in setmpi (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Added the setmpi15D routine.
!                              o Added returns for when we need to
!                                abort the run.
!                              o Changed the MPI communicator from
!                                MPI_COMM_WORLD to MPI_COMM_RT.
!                              o Only the global master can write
!                                outputs here.
!                             (TdPA)
!
!     06/21/2022:    V2.0.5 - Added detailed output for the split
!                             of tasks (TdPA)
!
!     12/16/2021:    V2.0.4 - Bugfix: Removed wrong nxtran factor when
!                             predicting the size of the JKQgen MPI
!                             package, resulting in the impossibility
!                             of solving for big atomic models without
!                             a real reson (TdPA)
!
!     03/24/2021:    V2.0.3 - Bugfix: Wrong initialization for the
!                             Atom%omp_2c%mxnU1 array (TdPA)
!
!     03/23/2021:    V2.0.2 - Removed calculation of minto and maxto
!                             for output transitions that was not
!                             needed (TdPA)
!
!     03/19/2021:    V2.0.1 - Bugfix: In setomp_2ord, the limits in
!                             the transition loop must be limited
!                             to the existing Frec%dzao dimensions.
!                             It was including all transitions (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Added routines to prepare the OpenmP
!                             split: setomp, setomp_magn,
!                             setomp_magn_2ord, and setomp_2ord, and
!                             set_index (TdPA)
!
!     09/20/2020:    V1.5.4 - Missing "/" in RAM file names (TdPA)
!
!     09/11/2020:    V1.5.3 - Added RAMreport routine (TdPA)
!
!     11/19/2019:    V1.5.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     08/19/2019:    V1.5.1 - Bugfix: Moved altern variable definition
!                             outside of the conditionals, so it is
!                             initialized for every case (TdPA)
!
!     05/31/2019:    V1.5.0 - Some MPI sizes require running the
!                             resize routine, which requires running
!                             setmpi. Therefore, the size
!                             determination has been removed from
!                             setmpi and now there is a second routine
!                             setmpi_sizes that computes them (TdPA)
!
!     04/15/2019:    V1.4.1 - The alternative bcast is optional (TdPA)
!                           - Modified sizes due to the change on
!                             radiation field tensor dimensions (TdPA)
!
!     02/20/2019:    V1.4.0 - New verbosity (TdPA)
!                           - Now uses unit 200 (TdPA)
!
!     09/25/2018:    V1.3.1 - Added a column to MPI%requestA (TdPA)
!
!     09/04/2018:    V1.3.0 - Added a correction to the sizes of the
!                             packages for MPI in the solvers if the
!                             size of the buffers is excesive (TdPA)
!                           - Added a stop if the buffers for
!                             communication are too big even with the
!                             alternative method (TdPA)
!
!     08/04/2018:    V1.2.5 - Initialize RAM parameter (TdPA)
!
!     06/29/2018:    V1.2.4 - Bugfix: I apparently forgot to define
!                             senders and receivers for the two CPU
!                             case. Therefore, it crashed as soon
!                             as there was a bcast like section in the
!                             source (TdPA)
!
!     05/16/2018:    V1.2.3 - Changed the size of Stokes and StokesI
!                             to take into account that are dimension
!                             nPh, and not nPh2 (TdPA)
!
!     10/31/2017:    V1.2.2 - Bugfix: The send_tree algorithm failed
!                             for exactly 4 CPU (TdPA)
!                           - Added message for cases when the spread
!                             algorithm reaches maximum number of
!                             tries (TdPA)
!
!     10/18/2017:    V1.2.1 - Introduced a maximum number of tries to
!                             optimize the frequency spread. There
!                             could be some combiations were it
!                             just oscillates (TdPA)
!
!     09/27/2017:    V1.2.0 - Created variables to apply send_tree
!                             algorithm (TdPA)
!
!     08/31/2017:    V1.1.5 - Refined frequency distribution
!                             algorithm (TdPA)
!
!     08/30/2017:    V1.1.4 - The format for the MPI_info output was
!                             missing numerals next to 'x' (TdPA)
!
!     08/24/2017:    V1.1.3 - Bugfix in the frequency distributer
!                             algorithm (TdPA)
!
!     08/22/2017:    V1.1.2 - Added more sizes (TdPA)
!
!     08/21/2017:    V1.1.1 - Changed sizes to comply with the
!                             changes in the solver routines (TdPA)
!                           - Bugfix: When only domain decomposition
!                             and no frequency split, the CPUs were
!                             assigned only one frequency (TdPA)
!                           - Allocated new sizes
!
!     08/01/2017:    V1.1.0 - Revamp of frequency distribution (TdPA)
!
!     07/31/2017:    V1.0.4 - Bugfix: Four TW in conditionals should
!                             be TWb. The routine could enter in
!                             an infinite loop (TdPA)
!                           - Bugfix: The CPU could only check
!                             backwards when optimizing the split of
!                             frequencies (TdPA)
!
!     07/10/2017:    V1.0.3 - One of the dimensions for MPI was wrong
!                             for the two CPU case, what leaded to
!                             the slave trying to send more than was
!                             allocated for one variable (TdPA)
!
!     06/23/2017:    V1.0.2 - Changed the way of splitting in
!                             frequencies to a weighted distribution
!                             of the nodes (TdPA)
!
!     06/22/2017:    V1.0.1 - Changed two sizes of MPI (TdPA)
!                           - Added request11 (TdPA)
!
!     04/19/2017:    V1.0.0 - First version (TdPA)
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
!  setmpi15D
!    Split tasks for the 1.5D case, creating groups of CPU that will
!  work independently in the RT problem. It creates a new communicator
!  to minimize changes in the 1D pure case.
!
!  setmpi
!    Split domains for the different processors
!
!  setmpi_CLE
!    Split domains for the different processors for CLE, that needs
!  to split two different axes
!
!  adjust_IW
!    Read weigths and times from a previous run to try readjust the
!  expected weights
!
!  setmpi_sizes
!    Compute the sizes of MPI messages
!
!  reset_mpirequest
!    Set all request variables to MPI_REQUEST_NULL
!
!  RAMreport
!    Writes in ASCII format the amomunt of RAM used which is
!  controlled by the MPI structure
!
!  setomp
!    Just checks the number of threads and initializes the global omp
!
!  setomp_magn
!    Decides how the first order RT coefficients are going to be
!  split in threads in those nodes where there is a magnetic field
!
!  setomp_magn_2ord
!    Prepares the thread split in magnetic components for emiss2ord
!
!  setomp_2ord
!    Prepares the thread split in emissI2ord and emiss2ordNB
!
!  set_index
!    Routine that returns the initial and final indexes to split a
!  domain in approximately equal amount of jobs, given the weight of
!  each node in the domain
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use funnj_mod
      use omp_mod
      use parameters_mod , only : TINYSP , TINYB , cZero , TINYCO , &
                                  TINYEV , TINYJS
      use types_mod

      ! Parameters

      ! Maximum size of solver buffer in MB
      double precision, parameter:: maxbuffer_personal = 5000
      double precision, parameter:: maxbuffer_supercom = 900

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Split the CPU in groups if doing 1.5D\n
      !!       MPID(MPI_class): Structure with MPI data\n
      !!    Input(Input_class): Structure with settings data\n
      subroutine setmpi15D(MPID,Input)

      ! I/O

      type(MPI_class), intent(out):: MPID
      type(Input_class), intent(inout):: Input

      ! Local
      integer:: nslave,ngroup,indx,jndx,lid,icolor
      integer, dimension(:), allocatable:: nelem
      integer, dimension(0:gnproc-1):: colors
#ifdef DEBUGSYN
      ! Only one group doing synthesis
      Input%rt_group_n = gnproc - 1
#endif


      !
      ! If only 1 CPU or whole MPI in RT, no groups
      !
      if (gnproc.eq.1) then

        ! Warn this may not be a good idea
        if (run_mode.eq.-1) then
          umsg = ' # It is recommended to use MPI for the 1.5D '// &
                 'inversion'
        else
          umsg = ' # It is recommended to use MPI for the 1.5D '// &
                 'synthesis'
        end if
        call verbose

        ! Trivial
        pid = gpid
        nproc = gnproc
        MPID%pid = pid
        MPID%nproc = nproc
        MPID%mpi15d = .False.
        MPID%ngroup = 1

        !
        ! Define communicator
        !
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, gpid, &
                            MPI_COMM_RT, ierr)
        call MPI_COMM_RANK(MPI_COMM_RT,pid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_RT,nproc,ierr)

      ! We are going to split columns
      else

        ! Number of slaves
        nslave = gnproc - 1

        ! Master by itself
        colors(0) = 0

        ! Control expectations
        if (Input%rt_group_n.gt.nslave) then
          Input%rt_group_n = nslave
          if (gpid.eq.0) then
            write(umsg,'(A,i5)') &
              ' # You choose too many RT slaves for the '// &
              'number of processes, set to ',nslave
            call verbose
          end if
        end if

        ! Number of groups
        ngroup = nslave/Input%rt_group_n
        MPID%ngroup = ngroup

        ! Allocate number of elements per group
        allocate(nelem(ngroup))
        nelem = Input%rt_group_n

        ! Add elements until all split
        indx = ngroup
        do while (sum(nelem).lt.nslave)
          if (indx.lt.1) indx = ngroup
          nelem(indx) = nelem(indx) + 1
          indx = indx-1
        end do

        ! Master allocate leaders
        if (gpid.eq.0) allocate(MPID%ltslave(ngroup))

        !
        ! Mark leaders and give colors
        !

        ! Initialize
        lid = 0
        icolor = 0

        ! For each group
        do indx=1,ngroup

          ! Advance leader and color
          lid = lid + 1
          icolor = icolor + 1
          colors(lid) = icolor
          if (gpid.eq.0) MPID%ltslave(indx) = lid

          ! Rest of members of the group
          do jndx=2,nelem(indx)

            ! Advance who and give color
            lid = lid + 1
            colors(lid) = icolor

          end do ! Members
        end do ! Groups

        !
        ! Define new communicators
        !
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, colors(gpid), gpid, &
                            MPI_COMM_RT, ierr)
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, colors(gpid), gpid, &
                            MPI_COMM_CTRL, ierr)
        call MPI_COMM_RANK(MPI_COMM_RT,pid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_RT,nproc,ierr)

        ! And set variables in MPI structure
        MPID%pid = pid
        MPID%nproc = nproc
        MPID%mpi15d = .True.

      end if ! Setting MPI environments

      end subroutine setmpi15D

!#####################################################################
!#####################################################################
!#####################################################################

      !> Assigns to each CPU a block of frequencies/heights, and
      !! initializes the isendtree algorithm for broadcasting.\n
      !!       MPID(MPI_class): Structure with MPI data\n
      !!    Input(Input_class): Structure with settings data\n
      !!   IW_freq(integer(:)): Vector with work weights for each
      !!                        frequency node
      subroutine setmpi(MPID,Input,IW_freq)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(in):: Input
      integer, dimension(:), allocatable:: IW_freq

      ! Local

      logical:: change, found

      integer:: ios,nnd,iaux,iproc,cproc,coun,maxcoun
      integer:: ifbl,ifblb,ii,is0,istep
      integer:: Mwg,Mwg1,Mwg2,mMwg,MTW,TTW,TW,TWM,TTWb,TWb,LW
      integer, dimension(:,:), allocatable:: itmp

      double precision:: daux


      ! Routine name
      urou = 'setmpi'


      !
      ! If we have more than one slave
      !
      if (nproc.gt.2) then

        ! We are doing MPI
        MPID%mpi = .True.

        ! Maximum number of tries to optimize frequency spread
        maxcoun = 1000000

        ! Pass bcast to MPI type
        MPID%altbcast = Input%altbcast


        !
        ! Allocate variables
        !

        ! Number of frequencies assigned
        allocate(MPID%nf(0:nproc-1))
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:nproc-1))
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:nproc-1))


        !
        ! Split frequencies between CPUs
        !

        ! Check if the axis is splittable itself
        nnd = nproc - 1
        MPID%nnd = nnd
        iaux = nfreq/nnd

        ! We need at least one frequency per processor
        if(iaux.lt.1)then
          umsg = 'Too many processors for this frequency grid'
          call aborted
          return
        end if

        ! Flag it if we are splitting
        MPID%lf = .True.

        ! Compute Maximum weight
        Mwg = maxval(IW_freq)

        ! Split frequencies between CPUs
        iaux = nfreq/nnd

        ! Determine number of frequencies per CPU
        MPID%nf(0) = 0
        do iproc=1,nnd
          MPID%nf(iproc) = iaux
        enddo

        ! Complex splitting if they do not coincide in size
        if(iaux*nnd.ne.nfreq)then
          iaux = nfreq - iaux*nnd
          do ifbl=1,iaux
            MPID%nf(ifbl) = MPID%nf(ifbl) + 1
          end do
        end if

        ! Determine the frequency limits for each processor in
        ! one group

        ! First CPU
        MPID%if0(1) = 1
        MPID%if1(1) = MPID%nf(1)

        ! For each group
        do ifbl=2,nnd
          MPID%if0(ifbl) = MPID%if1(ifbl-1) + 1
          MPID%if1(ifbl) = MPID%if0(ifbl) + MPID%nf(ifbl) - 1
        end do

        !
        ! Now optimize how they are distributed
        !

        ! Set logical variable
        change = .True.
        coun = 0

        ! Master look for limits
        if (pid.eq.0) then

          ! While there is change, repeat
          do while (change)

            ! Reset flag
            change = .False.

            MTW = sum(IW_freq(MPID%if0(1):MPID%if1(1)))
            do ifbl=2,nnd
              LW = sum(IW_freq(MPID%if0(ifbl):MPID%if1(ifbl)))
              if (LW.gt.MTW) MTW = LW
            end do

            ! First, check if there are very big differences and try
            ! big shifts
            do ifbl=1,nnd-1


              TW = sum(IW_freq(MPID%if0(ifbl):MPID%if1(ifbl)))
              Mwg1 = maxval(IW_freq(MPID%if0(ifbl):MPID%if1(ifbl)))


              do ifblb=ifbl+1,nnd

                TWb = sum(IW_freq(MPID%if0(ifblb):MPID%if1(ifblb)))
                Mwg2 = maxval(IW_freq(MPID%if0(ifblb): &
                                      MPID%if1(ifblb)))


                ! Minimum if maximum weights
                mMwg = min(Mwg1,Mwg2)

                if (abs(TW - TWb).gt.mMwg) then

                  if (TW.gt.TWb) then

                    if (MPID%nf(ifbl).gt.1) then

                      TTW = TW
                      TTWb = TWb

                      TWM = TW

                      do ii=ifbl,ifblb
                        LW = sum(IW_freq(MPID%if0(ii):MPID%if1(ii)))
                        if (LW.gt.TTW) TTW = LW
                      end do

                      LW = sum(IW_freq(MPID%if0(ifbl): &
                                       MPID%if1(ifbl)-1))
                      if (LW.gt.TWb) TWb = LW
                      if (LW.gt.TTWb) TTWb = LW
                      do ii=ifbl+1,ifblb-1
                        LW = sum(IW_freq(MPID%if0(ii)-1: &
                                         MPID%if1(ii)-1))
                        if (LW.gt.TTWb) TTWb = LW
                      end do
                      LW = sum(IW_freq(MPID%if0(ifblb)-1: &
                                       MPID%if1(ifblb)))
                      if (LW.gt.TWb) TWb = LW
                      if (LW.gt.TTWb) TTWb = LW

                      if (TWb.lt.TWM.and. &
                          TTWb.le.TTW+abs(TWb-TWM)) then
                        MPID%if1(ifbl) = MPID%if1(ifbl)-1
                        MPID%nf(ifbl) = MPID%nf(ifbl)-1
                        do ii=ifbl+1,ifblb-1
                          MPID%if0(ii) = MPID%if0(ii)-1
                          MPID%if1(ii) = MPID%if1(ii)-1
                        end do
                        MPID%if0(ifblb) = MPID%if0(ifblb)-1
                        MPID%nf(ifblb) = MPID%nf(ifblb)+1
                        TW = sum(IW_freq(MPID%if0(ifbl): &
                                         MPID%if1(ifbl)))
                        Mwg1 = maxval(IW_freq(MPID%if0(ifbl): &
                                              MPID%if1(ifbl)))
                        change = .True.
                      end if
                    end if

                  else

                    if (MPID%nf(ifblb).gt.1) then

                      TTW = TW
                      TTWb = TWb

                      TWM = TWb

                      do ii=ifbl,ifblb
                        LW = sum(IW_freq(MPID%if0(ii):MPID%if1(ii)))
                        if (LW.gt.TTW) TTW = LW
                      end do

                      TWb = sum(IW_freq(MPID%if0(ifbl): &
                                        MPID%if1(ifbl)+1))
                      do ii=ifbl+1,ifblb-1
                        LW = sum(IW_freq(MPID%if0(ii)+1: &
                                         MPID%if1(ii)+1))
                        if (LW.gt.TWb) TWb = LW
                      end do
                      LW = sum(IW_freq(MPID%if0(ifblb)+1: &
                                       MPID%if1(ifblb)))
                      if (LW.gt.TWb) TWb = LW

                      if (TWb.lt.TWM.and. &
                          TTWb.le.TTW+abs(TWb-TWM)) then
                        MPID%if1(ifbl) = MPID%if1(ifbl)+1
                        MPID%nf(ifbl) = MPID%nf(ifbl)+1
                        do ii=ifbl+1,ifblb-1
                          MPID%if0(ii) = MPID%if0(ii)+1
                          MPID%if1(ii) = MPID%if1(ii)+1
                        end do
                        MPID%if0(ifblb) = MPID%if0(ifblb)+1
                        MPID%nf(ifblb) = MPID%nf(ifblb)-1
                        TW = sum(IW_freq(MPID%if0(ifbl): &
                                         MPID%if1(ifbl)))
                        Mwg1 = maxval(IW_freq(MPID%if0(ifbl): &
                                              MPID%if1(ifbl)))
                        change = .True.
                      end if

                    end if
                  end if

                end if

              end do
            end do

            ! Cycle of changes
            coun = coun + 1

            ! If we have changed more than the limit
            if (coun.gt.maxcoun) then
              change = .False.
              if (pid.eq.0) then
                umsg = ' # Limit of iterations in task distributer'
                call verbose
              end if
            end if

          end do ! while changing

          ! Master have everything
          MPID%if0(0) = 1
          MPID%if1(0) = nfreq
          MPID%nf(0) = nfreq

        end if ! Master

        ! Send to everyone
        call MPI_BCAST(MPID%nf(0), nproc, MPI_INTEGER, 0, &
                       MPI_COMM_RT, ierr)
        call MPI_BCAST(MPID%if0(0), nproc, MPI_INTEGER, 0, &
                       MPI_COMM_RT, ierr)
        call MPI_BCAST(MPID%if1(0), nproc, MPI_INTEGER, 0, &
                       MPI_COMM_RT, ierr)

        ! Define maximums of sizes

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = maxval(MPID%nf(1:nproc-1))

        !
        ! Define variables for custom BCAST
        !

        if (Input%altbcast) then

          ! Calculate steps
          daux = log(dble(nproc))/log(2d0)
          if (abs(daux - dble(int(daux))).gt.1d-16) then
            MPID%steps = int(daux) + 1
          else
            MPID%steps = int(daux)
          end if

          ! Allocations
          allocate(itmp(0:MPID%steps,0:nproc-1))
          itmp = -1
          itmp(0,0) = 0

          !
          ! Mount vector
          !

          ! Initialize number of sends
          MPID%nsend = 0

          ! Initialize recv for master
          if (pid.eq.0) MPID%recv = 0

          ! Initialize processor to acknowledge
          ii = 0

          ! Master sending

          ! For every step
          do istep=1,MPID%steps

            ! Advance the processor index
            ii = ii + 1

            ! If we have split everything, finish
            if (ii.eq.nproc) exit

            ! In istep, master will send to ii
            itmp(istep,0) = ii

            ! If the master is looking at this, store
            ! the data
            if (pid.eq.0) MPID%nsend = MPID%nsend + 1

            ! If ii is looking at this, store that will receive
            if (pid.eq.ii) MPID%recv = 0

          end do

          ! Slaves sending and reeciving

          ! Do until we run out of processors
          do while (.True.)

            ! If we ran out of processors, exit
            if (ii.eq.nproc) exit

            ! Search for the next processor to be sender

            ! Check in every step
            do istep=1,MPID%steps
              ! Check over all processors
              do iproc=0,nproc-1

                ! Initialize flag
                found = .False.

                ! If the processor has received before
                if (itmp(istep,iproc).ge.0) then

                  ! This would be the new sender
                  cproc = itmp(istep,iproc)

                  ! And it would send in the next step
                  is0 = istep + 1

                  ! Check if it has not been the sender before
                  if (itmp(0,cproc).lt.0) then

                    ! If not, we found it
                    found = .True.
                    exit

                  end if ! was not sender before

                end if ! has received already
              end do ! processors
              if (found) exit
            end do ! steps

            ! If we did not find sender, finish
            if (.not.found) exit

            ! For all remaining steps
            do istep=is0,MPID%steps

              ! Advance the processor to receive
              ii = ii + 1

              ! If ran out of processors, we are done
              if (ii.eq.nproc) exit

              ! Store who is sending to
              itmp(istep,cproc) = ii

              ! Flag as been sender already
              itmp(0,cproc) = 0

              ! Update number of sends
              if (pid.eq.cproc) MPID%nsend = MPID%nsend + 1

              ! Update the receiver
              if (pid.eq.ii) MPID%recv = cproc

            end do ! steps

            ! If ran out of processors, finish
            if (ii.eq.nproc-1) exit

          end do ! While there is processes to assign

          ! If the processor is sending
          if (MPID%nsend.gt.0) then

            ! Allocate list of sends
            allocate(MPID%lsend(MPID%nsend))

            ! Initialize index
            ii = 0

            ! For every step
            do istep=1,MPID%steps

              ! If there is a send
              if (itmp(istep,pid).gt.0) then

                ! Shift index
                ii = ii+1

                ! Store the send
                MPID%lsend(ii) = itmp(istep,pid)

              end if ! There is a send

            end do ! steps

            ! Deallocate the temporal array
            deallocate(itmp)

          end if ! There is sending from this processor

        ! Normal bcast
        else

          MPID%nsend = 0

        end if ! Alternative bcast


        !
        ! Allocate request
        !

        ! Alternative bcast
        if (MPID%altbcast) then

          allocate(MPID%requestA(MPID%nsend,8))

        end if ! Domain decomposition


      !
      ! If there is only one slave
      !
      else if (nproc.eq.2) then

        ! We are doing MPI
        MPID%mpi = .True.

        ! Pass bcast to MPI type
        MPID%altbcast = Input%altbcast

        ! The number of processors
        nnd = 1
        MPID%nnd = 1

        !
        ! Allocate variables
        !

        ! Number of frequencies assigned
        allocate(MPID%nf(0:nproc-1))
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:nproc-1))
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:nproc-1))

        ! There is no split in frequencies
        MPID%lf = .False.

        ! Determine number of frequencies per CPU
        MPID%nf(0) = nfreq
        MPID%nf(1) = nfreq

        ! Define maximums of sizes

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = nfreq

        ! Set the limits of each processor
        MPID%if0(0) = 1
        MPID%if0(1) = 1
        MPID%if1(0) = nfreq
        MPID%if1(1) = nfreq


        !
        ! Define variables for custom BCAST
        !

        if (Input%altbcast) then

          ! Number of steps
          MPID%steps = 1

          ! Receive from
          MPID%recv = 0

          ! Master
          if (pid.eq.0) then

            MPID%nsend = 1
            allocate(MPID%lsend(1))
            allocate(MPID%requestA(1,8))

            MPID%lsend(1) = 1

          ! slave
          else

            MPID%nsend = 0
            MPID%recv = 0

          end if ! master or slave

        ! Normal bcast
        else

          MPID%nsend = 0

        end if ! alternative bcast


      !
      ! If we are in serial
      !
      else

        ! There is no MPI in any form
        MPID%mpi = .False.
        MPID%lf = .False.

        ! Number of frequencies assigned
        allocate(MPID%nf(0:0))
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:0))
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:0))

        ! The only processor sees everything
        MPID%nf(0) = nfreq
        MPID%if0(0) = 1
        MPID%if1(0) = nfreq

        !
        ! Define maximums of sizes
        !

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = nfreq

      end if ! Many, 2, or 1 CPU

      ! Reset the requests to null
      call reset_mpirequest(MPID)

      ! Write information about the distribution of tasks
      if (gpid.eq.0.and.Input%keep_mpil) then

        ! Open file
        open (200,file=trim(Input%folder)//'/MPI_info', &
              status='unknown', iostat=ios, action='write')

        write(200,'("Running with",1x,i4,1x,"threads and",'// &
                   '1x,i4,1x,"processes")') &
              nthread,nproc

        ! Write MPI info
        do iproc=1,nproc-1

          write(200,'("Process",1x,i4,1x,"received",1x,i6,1x,'// &
                    '"frequencies, with a total weight of",'// &
                    '1x,i7,1x," with maximum of",1x,i7)')&
                iproc,MPID%nf(iproc), &
                sum(IW_freq(MPID%if0(iproc):MPID%if1(iproc))), &
                maxval(IW_freq(MPID%if0(iproc):MPID%if1(iproc)))

        end do

        ! Close file
        close(200)

      end if

      ! Write detailed information about the distribution of tasks
      if (gpid.eq.0.and.Input%keep_mpidl) then

        ! Open file
        open (200,file=trim(Input%folder)//'/MPI_info_detail', &
              status='unknown', iostat=ios, access='stream', &
              action='write', form='unformatted')

        ! Write number of frequencies and weights
        write(200) nfreq
        write(200) IW_freq

        ! Now write the data in MPI_info but in binary
        write(200) nproc
        write(200) MPID%if0
        write(200) MPID%if1
        write(200) MPID%nf

        ! Close file
        close(200)

      end if ! Master

      ! Initialize flags for alternative MPI
      MPID%alternI = .False.
      MPID%alternP = .False.
      MPID%alternJ = .False.
      MPID%alternJgen = .False.

      ! Deallocate IW_freq
      deallocate(IW_freq)

      ! Initialize RAM parameter
      MPID%RAM = 0d0
      MPID%PRAM = 0d0

      ! Check if everything is fine
      call control

      return

      end subroutine setmpi

!#####################################################################
!#####################################################################
!#####################################################################

      !> Assigns to each CPU a block of frequencies/heights, and
      !! initializes the isendtree algorithm for broadcasting. This
      !! one is special for CLE because there is not a true master
      !! process.\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!     Input(Input_class): Structure with settings data\n
      !! IW_freq_ou(integer(:)): Vector with work weights for each
      !!                         frequency node for the output\n
      !! IW_freq_in(integer(:)): Vector with work weights for each
      !!                         frequency node for the input
      subroutine setmpi_CLE(MPID,Input,IW_freq_ou,IW_freq_in)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(in):: Input
      integer, dimension(:), target:: IW_freq_in,IW_freq_ou

      ! Local

      logical:: change, found

      integer:: ios,iaux,iproc,cproc,coun,maxcoun
      integer:: ifbl,ifblb,ii,is0,istep,iaxis,mfreq,iran,jfreq,ifreq
      integer:: Mwg,Mwg1,Mwg2,mMwg,MTW,TTW,TW,TWM,TTWb,TWb,LW
      integer, dimension(0:nproc-1):: nf,if0,if1
      integer, dimension(:,:), allocatable:: itmp

      double precision:: daux

      integer, dimension(:), pointer:: IW_freq


      ! Routine name
      urou = 'setmpi'


      !
      ! If we have more than one CPU
      !
      if (nproc.gt.1) then

        ! We are doing MPI
        MPID%mpi = .True.

        ! Maximum number of tries to optimize frequency spread
        maxcoun = 1000000

        ! Pass bcast to MPI type
        MPID%altbcast = Input%altbcast

        ! Flag it if we are splitting
        MPID%lf = .True.


        !
        ! Allocate variables
        !

        ! Number of frequencies assigned
        allocate(MPID%nf(0:nproc-1))
        allocate(MPID%inf(0:nproc-1))
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:nproc-1))
        allocate(MPID%iif0(0:nproc-1))
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:nproc-1))
        allocate(MPID%iif1(0:nproc-1))


        !
        ! Run twice, for input and output axes
        !
        do iaxis=1,2

          ! iaxis 1 ->  INPUT PART
          ! iaxis 2 -> OUTPUT PART

          ! If axis 1
          if (iaxis.eq.1) then

            ! Get input frequencies (size)
            mfreq = nfreq

            ! Point to IW_freq
            IW_freq => IW_freq_in

          ! Axis is 2
          else

            ! Check range is limited
            if (Input%lim_stk%nran.le.0) then

              ! Copy in and finish
              MPID%nf = MPID%inf
              MPID%if0 = MPID%iif0
              MPID%if1 = MPID%iif1
              exit

            end if ! Limited range

            ! Get input frequencies (size)
            mfreq = Input%lim_stk%nn

            ! Allocate IW_freq
            allocate(IW_freq(mfreq))

            !
            ! Let's copy only relevant IW_freq
            !

            ! Initialize rolling frequency index
            ifreq = 0

            ! For each range
            do iran=1,Input%lim_stk%nran

              ! For each frequency within the range
              do jfreq=Input%lim_stk%indx(1,iran), &
                       Input%lim_stk%indx(2,iran)

                ! Advance rolling index
                ifreq = ifreq + 1

                ! Save IW_freq for input in continuous order
                IW_freq(ifreq) = IW_freq_ou(jfreq)

              end do ! For each frequency within the range
            end do ! For each range

          end if

          !
          ! Split frequencies between CPUs
          !

          ! Check if the axis is splittable itself
          MPID%nnd = nproc
          iaux = mfreq/nproc

          ! We need at least one frequency per processor
          if(iaux.lt.1.and.iaxis.eq.1)then
            umsg = 'Too many processors for this frequency grid'
            call aborted
            return
          end if

          ! Compute Maximum weight
          Mwg = maxval(IW_freq)

          ! Split frequencies between CPUs
          iaux = mfreq/nproc

          ! Determine number of frequencies per CPU
          nf = iaux

          ! Complex splitting if they do not coincide in size
          if(iaux*nproc.ne.mfreq)then
            iaux = mfreq - iaux*nproc
            do ifbl=0,iaux-1
              nf(ifbl) = nf(ifbl) + 1
            end do
          end if

          ! Determine the frequency limits for each processor in
          ! one group

          ! First CPU
          if0(0) = 1
          if1(0) = nf(0)

          ! For each group
          do ifbl=1,nproc-1
            if0(ifbl) = if1(ifbl-1) + 1
            if1(ifbl) = if0(ifbl) + nf(ifbl) - 1
          end do

          !
          ! Now optimize how they are distributed
          !

          ! Set logical variable
          change = .True.
          coun = 0

          ! First CPU looks for limits
          if (pid.eq.0) then

            ! While there is change, repeat
            do while (change)

              ! Reset flag
              change = .False.

              MTW = sum(IW_freq(if0(0):if1(0)))
              do ifbl=1,nproc-1
                LW = sum(IW_freq(if0(ifbl):if1(ifbl)))
                if (LW.gt.MTW) MTW = LW
              end do

              ! First, check if there are very big differences and try
              ! big shifts
              do ifbl=0,nproc-2


                TW = sum(IW_freq(if0(ifbl):if1(ifbl)))
                Mwg1 = maxval(IW_freq(if0(ifbl):if1(ifbl)))


                do ifblb=ifbl+1,nproc-1

                  TWb = sum(IW_freq(if0(ifblb):if1(ifblb)))
                  Mwg2 = maxval(IW_freq(if0(ifblb): &
                                        if1(ifblb)))


                  ! Minimum if maximum weights
                  mMwg = min(Mwg1,Mwg2)

                  if (abs(TW - TWb).gt.mMwg) then

                    if (TW.gt.TWb) then

                      if (nf(ifbl).gt.1) then

                        TTW = TW
                        TTWb = TWb

                        TWM = TW

                        do ii=ifbl,ifblb
                          LW = sum(IW_freq(if0(ii):if1(ii)))
                          if (LW.gt.TTW) TTW = LW
                        end do

                        LW = sum(IW_freq(if0(ifbl): &
                                         if1(ifbl)-1))
                        if (LW.gt.TWb) TWb = LW
                        if (LW.gt.TTWb) TTWb = LW
                        do ii=ifbl+1,ifblb-1
                          LW = sum(IW_freq(if0(ii)-1: &
                                           if1(ii)-1))
                          if (LW.gt.TTWb) TTWb = LW
                        end do
                        LW = sum(IW_freq(if0(ifblb)-1: &
                                         if1(ifblb)))
                        if (LW.gt.TWb) TWb = LW
                        if (LW.gt.TTWb) TTWb = LW

                        if (TWb.lt.TWM.and. &
                            TTWb.le.TTW+abs(TWb-TWM)) then
                          if1(ifbl) = if1(ifbl)-1
                          nf(ifbl) = nf(ifbl)-1
                          do ii=ifbl+1,ifblb-1
                            if0(ii) = if0(ii)-1
                            if1(ii) = if1(ii)-1
                          end do
                          if0(ifblb) = if0(ifblb)-1
                          nf(ifblb) = nf(ifblb)+1
                          TW = sum(IW_freq(if0(ifbl): &
                                           if1(ifbl)))
                          Mwg1 = maxval(IW_freq(if0(ifbl): &
                                                if1(ifbl)))
                          change = .True.
                        end if
                      end if

                    else

                      if (nf(ifblb).gt.1) then

                        TTW = TW
                        TTWb = TWb

                        TWM = TWb

                        do ii=ifbl,ifblb
                          LW = sum(IW_freq(if0(ii):if1(ii)))
                          if (LW.gt.TTW) TTW = LW
                        end do

                        TWb = sum(IW_freq(if0(ifbl): &
                                          if1(ifbl)+1))
                        do ii=ifbl+1,ifblb-1
                          LW = sum(IW_freq(if0(ii)+1: &
                                           if1(ii)+1))
                          if (LW.gt.TWb) TWb = LW
                        end do
                        LW = sum(IW_freq(if0(ifblb)+1: &
                                         if1(ifblb)))
                        if (LW.gt.TWb) TWb = LW

                        if (TWb.lt.TWM.and. &
                            TTWb.le.TTW+abs(TWb-TWM)) then
                          if1(ifbl) = if1(ifbl)+1
                          nf(ifbl) = nf(ifbl)+1
                          do ii=ifbl+1,ifblb-1
                            if0(ii) = if0(ii)+1
                            if1(ii) = if1(ii)+1
                          end do
                          if0(ifblb) = if0(ifblb)+1
                          nf(ifblb) = nf(ifblb)-1
                          TW = sum(IW_freq(if0(ifbl):if1(ifbl)))
                          Mwg1 = maxval(IW_freq(if0(ifbl): &
                                                if1(ifbl)))
                          change = .True.
                        end if

                      end if
                    end if

                  end if

                end do
              end do

              ! Cycle of changes
              coun = coun + 1

              ! If we have changed more than the limit
              if (coun.gt.maxcoun) then
                change = .False.
                if (pid.eq.0) then
                  umsg = ' # Limit of iterations in task distributer'
                  call verbose
                end if
              end if

            end do ! while changing

          end if ! Task 0

          ! Send to everyone
          call MPI_BCAST(nf(0), nproc, MPI_INTEGER, 0, &
                         MPI_COMM_RT, ierr)
          call MPI_BCAST(if0(0), nproc, MPI_INTEGER, 0, &
                         MPI_COMM_RT, ierr)
          call MPI_BCAST(if1(0), nproc, MPI_INTEGER, 0, &
                         MPI_COMM_RT, ierr)

          ! Input axis
          if (iaxis.eq.1) then

            ! Copy to relevant variable in structure and free
            MPID%inf = nf
            MPID%iif0 = if0
            MPID%iif1 = if1
            nullify(IW_freq)

          ! Output axis
          else

            ! Copy to relevant variable in structure and free
            MPID%nf = nf
            MPID%if0 = if0
            MPID%if1 = if1
            deallocate(IW_freq)
            nullify(IW_freq)

          end if ! Input or output axis

        end do ! Axes

        ! Define maximums of sizes

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = maxval(MPID%inf(0:nproc-1))

        !
        ! Define variables for custom BCAST
        !

        if (Input%altbcast) then

          ! Calculate steps
          daux = log(dble(nproc))/log(2d0)
          if (abs(daux - dble(int(daux))).gt.1d-16) then
            MPID%steps = int(daux) + 1
          else
            MPID%steps = int(daux)
          end if

          ! Allocations
          allocate(itmp(0:MPID%steps,0:nproc-1))
          itmp = -1
          itmp(0,0) = 0

          !
          ! Mount vector
          !

          ! Initialize number of sends
          MPID%nsend = 0

          ! Initialize recv for master
          if (pid.eq.0) MPID%recv = 0

          ! Initialize processor to acknowledge
          ii = 0

          ! Master sending

          ! For every step
          do istep=1,MPID%steps

            ! Advance the processor index
            ii = ii + 1

            ! If we have split everything, finish
            if (ii.eq.nproc) exit

            ! In istep, master will send to ii
            itmp(istep,0) = ii

            ! If the master is looking at this, store
            ! the data
            if (pid.eq.0) MPID%nsend = MPID%nsend + 1

            ! If ii is looking at this, store that will receive
            if (pid.eq.ii) MPID%recv = 0

          end do

          ! Slaves sending and reeciving

          ! Do until we run out of processors
          do while (.True.)

            ! If we ran out of processors, exit
            if (ii.eq.nproc) exit

            ! Search for the next processor to be sender

            ! Check in every step
            do istep=1,MPID%steps
              ! Check over all processors
              do iproc=0,nproc-1

                ! Initialize flag
                found = .False.

                ! If the processor has received before
                if (itmp(istep,iproc).ge.0) then

                  ! This would be the new sender
                  cproc = itmp(istep,iproc)

                  ! And it would send in the next step
                  is0 = istep + 1

                  ! Check if it has not been the sender before
                  if (itmp(0,cproc).lt.0) then

                    ! If not, we found it
                    found = .True.
                    exit

                  end if ! was not sender before

                end if ! has received already
              end do ! processors
              if (found) exit
            end do ! steps

            ! If we did not find sender, finish
            if (.not.found) exit

            ! For all remaining steps
            do istep=is0,MPID%steps

              ! Advance the processor to receive
              ii = ii + 1

              ! If ran out of processors, we are done
              if (ii.eq.nproc) exit

              ! Store who is sending to
              itmp(istep,cproc) = ii

              ! Flag as been sender already
              itmp(0,cproc) = 0

              ! Update number of sends
              if (pid.eq.cproc) MPID%nsend = MPID%nsend + 1

              ! Update the receiver
              if (pid.eq.ii) MPID%recv = cproc

            end do ! steps

            ! If ran out of processors, finish
            if (ii.eq.nproc-1) exit

          end do ! While there is processes to assign

          ! If the processor is sending
          if (MPID%nsend.gt.0) then

            ! Allocate list of sends
            allocate(MPID%lsend(MPID%nsend))

            ! Initialize index
            ii = 0

            ! For every step
            do istep=1,MPID%steps

              ! If there is a send
              if (itmp(istep,pid).gt.0) then

                ! Shift index
                ii = ii+1

                ! Store the send
                MPID%lsend(ii) = itmp(istep,pid)

              end if ! There is a send

            end do ! steps

            ! Deallocate the temporal array
            deallocate(itmp)

          end if ! There is sending from this processor

        ! Normal bcast
        else

          MPID%nsend = 0

        end if ! Alternative bcast


        !
        ! Allocate request
        !

        ! Alternative bcast
        if (MPID%altbcast) then

          allocate(MPID%requestA(MPID%nsend,8))

        end if ! Domain decomposition

      !
      ! If we are in serial
      !
      else

        ! There is no MPI in any form
        MPID%mpi = .False.
        MPID%lf = .False.

        ! Number of frequencies assigned
        allocate(MPID%nf(0:0))
        allocate(MPID%inf(0:0))
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:0))
        allocate(MPID%iif0(0:0))
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:0))
        allocate(MPID%iif1(0:0))

        ! The only processor sees everything
        MPID%inf(0) = nfreq
        MPID%iif0(0) = 1
        MPID%iif1(0) = nfreq

        ! If not limiting the ranges
        if (Input%lim_stk%nran.le.0) then

          MPID%nf(0) = nfreq
          MPID%if0(0) = 1
          MPID%if1(0) = nfreq

        ! Limiting ranges
        else

          MPID%inf(0) = Input%lim_stk%nn
          MPID%iif0(0) = 1
          MPID%iif1(0) = Input%lim_stk%nn

        end if

        !
        ! Define maximums of sizes
        !

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = nfreq

      end if ! Many, 2, or 1 CPU

      ! Reset the requests to null
      call reset_mpirequest(MPID)

      ! Write information about the distribution of tasks
      if (gpid.eq.0.and.Input%keep_mpil) then

        ! Open file
        open (200,file=trim(Input%folder)//'/MPI_info', &
              status='unknown', iostat=ios, action='write')

        write(200,'("Running with",1x,i4,1x,"threads and",'// &
                   '1x,i4,1x,"processes")') &
              nthread,nproc

        ! Write MPI info
        do iproc=1,nproc-1

          write(200,'("Process",1x,i4,1x,"received",1x,i6,1x,'// &
                    '"frequencies, with a total weight of",'// &
                    '1x,i7,1x," with maximum of",1x,i7)')&
                iproc,MPID%inf(iproc), &
                sum(IW_freq_in(MPID%iif0(iproc):MPID%iif1(iproc))), &
                maxval(IW_freq_in(MPID%iif0(iproc):MPID%iif1(iproc)))

        end do

        ! Close file
        close(200)

      end if ! Master

      ! Check if everything is fine
      call control

      return

      end subroutine setmpi_CLE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Readjust weights based on information from a previous run.\n
      !!    Input(Input_class): Structure with settings data\n
      !!   IW_freq(integer(:)): Vector with work weights for each
      !!                        frequency node
      subroutine adjust_IW(Input,IW_freq)

      ! I/O

      type(Input_class), intent(in):: Input
      integer, dimension(:), intent(inout):: IW_freq

      ! Local

      integer:: ios,infreq,inproc,iter,cpu
      integer, dimension(:), allocatable:: if0,if1,nf
      integer, dimension(:), allocatable:: iIW_freq
      double precision:: time,ff
      double precision, dimension(:), allocatable:: tot,tim,coun

      ! Skip if none
      if (trim(Input%MPIdetail).eq.'N'.or. &
          trim(Input%operform).eq.'N') return

      ! If master
      if (pid.eq.0) then

        !
        ! Read detailed processor data
        !

        ! Open file
        open (200,file=trim(Input%MPIdetail), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', form='unformatted')

        ! Read head
        read(200,err=1100) infreq
        allocate(iIW_freq(infreq))
        read(200,err=1100) iIW_freq

        ! Only if same number of frequencies
        if (infreq.eq.nfreq) then

          ! Continue reading
          read(200,err=1100) inproc
          allocate(if0(0:inproc-1),if1(0:inproc-1),nf(0:inproc-1))
          read(200,err=1100) if0
          read(200,err=1100) if1
          read(200,err=1100) nf

          ! Close file
          close(200)

          !
          ! Allocate and initialize performance data
          !
          allocate(coun(inproc-1),tot(inproc-1),tim(inproc-1))
          coun = 0d0
          tot = 0d0
          tim = 0d0

          !
          ! Read performance data
          !
          open (200,file=trim(Input%operform), &
                iostat=ios,err=2000)

          ! Until EoF
          do while (.True.)

            ! Read line of data
            read(200,*,iostat=ios,end=1111) iter,cpu,time
            if (ios.ne.0) exit

            ! Ignore master
            if (cpu.eq.0) cycle

            ! If skipping first
            if (Input%IWskip.and.iter.eq.1) cycle

            ! Add to CPU
            coun(cpu) = coun(cpu) + 1
            tim(cpu) = tim(cpu) + time - tot(cpu)
            tot(cpu) = time

          end do ! Until EoF

1111      close(200)

          ! Check if data for everyone
          if (minval(coun).le.0d0) then

            ! Message (failure)
            urou = 'adjust_IW'
            write(umsg,'(A,i5)') &
                ' # There is no data for every CPU in the '// &
                'provided performance file.'
                call abortedS(umsg,urou,-1,.True.,.True.)

          else

            ! Average times
            tim = tim/coun

            ! Norm
            ff = dble(inproc-1)/sum(tim)

            ! Normalize sections
            do cpu=1,inproc-1

              ! Scale IW_freq
              IW_freq(if0(cpu):if1(cpu)) = &
                 nint(dble(iIW_freq(if0(cpu):if1(cpu)))*tim(cpu)*ff)

            end do


            ! Message (success)
            write(umsg,'(A,i5)') &
                ' - Adjusted the weigths based in a previous '// &
                'run performance.'
                call verbose

            call control

          end if

        ! Message (failure)
        else

          urou = 'adjust_IW'
          write(umsg,'(A,i5)') &
              ' # The number of frequencies in this problem '// &
              'differs with the one in the supplied MPI data. '// &
              'It cannot be used to improve estimations.'
              call abortedS(umsg,urou,-1,.True.,.True.)

        end if ! Same number of frequencies

      ! Slave
      else

        call control

      end if

      ! If aborting
      if (laborted) return

      ! Share new weights
      call MPI_BCAST(IW_freq(1), nfreq, MPI_INTEGER, 0, &
                     MPI_COMM_RT, ios)

      return

1000  umsg = 'Error opening MPI details file'
      urou = 'adjust_IW'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error reading MPI details file'
      urou = 'adjust_IW'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
2000  umsg = 'Error opening MPI performance file'
      urou = 'adjust_IW'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine adjust_IW

!#####################################################################
!#####################################################################
!#####################################################################

      !> Assigns the sizes of MPI messages.\n
      !!       MPID(MPI_class): Structure with MPI data\n
      !! GeomI(Geometry_class): Structure with geometry data for the
      !!                        intensity problem\n
      !!  Geom(Geometry_class): Structure with geometry data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          lio(logical): If solving for intensity\n
      !!           lp(logical): If solving for polarization\n
      !!           lJ(logical): If solving for continuum\n
      !!         lGen(logical): If correcting for multi-term\n
      !!        reval(logical): In here to reevaluate the sizes
      subroutine setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,lJ,lGen, &
                              reval)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Geometry_class), intent(in):: Geom,GeomI
      type(Frequency_class), intent(in):: Frec
      logical, intent(in):: lio,lp,lJ,lGen,reval

      ! Local

      integer:: iproc,istep

      double precision:: b1,b2,b3,b4,maxbuffer


      !
      ! Determine the maximum size of buffer depending on the number
      ! of CPU
      !
      if (nproc.gt.5) then
        maxbuffer = maxbuffer_supercom
      else
        maxbuffer = maxbuffer_personal
      end if

      ! If reevaluating
      if (reval) then

        ! If buffer did not change, just return
        if (MPID%size3(0).eq.4*MPID%nf(0)*Rnz) return


        !
        ! Compute new buffer sizes
        !
        do iproc=0,nproc-1

          ! Contribution function to master
          MPID%size3(iproc) = 4*MPID%nf(iproc)*Rnz
          ! Stokes (full dependence) to master
          MPID%size4(iproc) = 4*MPID%nf(iproc)*Geom%nTh*Geom%nPh*Rnz
          ! Profile to master
          MPID%size5(iproc) = 2*Frec%Mntfreq(iproc)*Geom%nTh* &
                             Geom%nPh*Rnz
          ! JKQ/JKQS to broadcast
          MPID%size6(iproc) = 5*3*nxtran*Rnz
          ! JKQ(nu) to broadcast
          MPID%size7(iproc) = 5*3*nfreq*Rnz
          ! Stokes (full dependence) to broadcast
          MPID%size8(iproc) = 4*nfreq*Geom%nTh*Geom%nPh*Rnz

          ! Lambda_photo to master
          MPID%sizei0(iproc) = nxb*Frec%Mnpfreq(iproc)* &
                               GeomI%nTh*GeomI%nPh*Rnz
          ! Lambda_line to broadcast
          MPID%sizei2(iproc) = nxb*2*nxphot*Rnz
          ! J00P to broadcast
          MPID%sizei3(iproc) = 2*nxphot*Rnz
          ! Stokes (full dependence) in intensity to master
          MPID%sizei4(iproc) = MPID%nf(iproc)*GeomI%nTh*GeomI%nPh*Rnz
          ! Duplicate
          MPID%sizei4b(iproc) = MPID%sizei4(iproc)
          ! Profile in intensity to master
          MPID%sizei5(iproc) = 2*Frec%Mntfreqi(iproc)*GeomI%nTh* &
                              GeomI%nPh*Rnz
          ! J00/J00S to broadcast
          MPID%sizei6(iproc) = nxt*Rnz
          ! J00(nu) to broadcast
          MPID%sizei7(iproc) = nfreq*Rnz
          ! Stokes (full dependence) in intensity to
          ! broadcast
          MPID%sizei8(iproc) = nfreq*GeomI%nTh*GeomI%nPh*Rnz
          ! Lambda_line to master
          MPID%sizei9(iproc) = nxb*Frec%Mntfreqi(iproc)*GeomI%nTh* &
                              GeomI%nPh*Rnz
          ! Lambda_photo to broadcast
          MPID%sizei10(iproc) = nxb*nxt*Rnz
          ! Contribution function in intensity to master
          MPID%sizei14(iproc) = MPID%nf(iproc)*Rnz

        end do


        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the continuum problem
        !

        ! If doing intensity
        if (lio.and.lJ) then

          ! Compute expected buffer size
          b1 = 8d-6*MPID%nxfreq*GeomI%nTh*GeomI%nPh*Rnz

          ! If buffers are going to be too big
          if (b1.gt.maxbuffer) then

            ! We need to use the alternative
            MPID%alternJ = .True.

            ! For each processor
            do iproc=0,nproc-1
              MPID%sizei4b(iproc) = MPID%nf(iproc)*Rnz
            end do

          end if ! Buffer too big
        end if ! Going through solveJ

        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the intensity problem
        !

        ! If doing intensity
        if (lio) then

          ! Compute expected buffer size
          b1 = 8d-6*MPID%nxfreq*GeomI%nTh*GeomI%nPh*Rnz
          b2 = 8d-6*MPID%nxtfreqi*2*GeomI%nTh*GeomI%nPh*Rnz
          b3 = 8d-6*MPID%nxtfreqi*nxb*GeomI%nTh*GeomI%nPh*Rnz
          b4 = 8d-6*MPID%nxpfreq*nxb*GeomI%nTh*GeomI%nPh*Rnz

          ! If buffers are going to be too big
          if ((b1+b2+b3+b4).gt.maxbuffer) then

            ! We need to use the alternative
            MPID%alternI = .True.

            ! For each processor
            do iproc=0,nproc-1
              MPID%sizei4(iproc) = MPID%nf(iproc)*Rnz
              MPID%sizei5(iproc) = 2*Frec%Mntfreqi(iproc)*Rnz
              MPID%sizei9(iproc) = nxb*Frec%Mntfreqi(iproc)*Rnz
              MPID%sizei0(iproc) = nxb*Frec%Mnpfreq(iproc)*Rnz
            end do

          end if ! Buffer to big
        end if ! Going through solveri

        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the polarization correction problem
        !

        ! If correcting intensity
        if (lGen) then

          if (PRD.and..not.AV) then

            ! Compute expected buffer size
            b1 = 8d-6*MPID%nxtfreq*2*Geom%nTh*Geom%nPh*Rnz

            ! If buffers are going to be too big
            if (b1.gt.maxbuffer) then

              ! We need to use the alternative
              MPID%alternJgen = .True.

            end if ! Big buffer
          end if ! PRD and AD
        end if ! Doing JKQ_gen


        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the polarization problem
        !

        ! If iterating polarization
        if (lp) then

          ! Compute expected buffer size
          b1 = 8d-6*MPID%nxfreq*4*Geom%nTh*Geom%nPh*Rnz
          b2 = 8d-6*MPID%nxtfreq*2*Geom%nTh*Geom%nPh*Rnz

          ! If buffers are going to be too big
          if ((b1+b2).gt.maxbuffer) then

            ! We need to use the alternative
            MPID%alternP = .True.

            do iproc=0,nproc-1
              MPID%size4(iproc) = 4*MPID%nf(iproc)*Rnz
              MPID%size5(iproc) = 2*Frec%Mntfreq(iproc)*Rnz
            end do

          end if
        end if

      !
      ! Not reevaluating
      !
      else

        !
        ! Share ntfreq information
        !

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive sizes
            call MPI_RECV(Frec%Mntfreq(0), nproc, MPI_INTEGER, &
                          MPID%recv, pid, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            ! Receive sizes
            call MPI_RECV(Frec%Mntfreqi(0), nproc, MPI_INTEGER, &
                          MPID%recv, 1000000+pid, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            ! Receive sizes
            call MPI_RECV(Frec%Mnpfreq(0), nproc, MPI_INTEGER, &
                          MPID%recv, 2000000+pid, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)

          end if ! No Master

          ! For each send
          do istep=1,MPID%nsend

            ! Send sizes
            call MPI_ISEND(Frec%Mntfreq(0), nproc, MPI_INTEGER, &
                           MPID%lsend(istep), MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,1), &
                           ierr)
            call MPI_ISEND(Frec%Mntfreqi(0), nproc, MPI_INTEGER, &
                           MPID%lsend(istep), &
                           1000000+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,2), &
                           ierr)
            call MPI_ISEND(Frec%Mnpfreq(0), nproc, MPI_INTEGER, &
                           MPID%lsend(istep), &
                           2000000+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,3), &
                           ierr)
          end do ! sends

        ! Normal bcast
        else

          call MPI_BCAST(Frec%Mntfreq(0), nproc, MPI_INTEGER, &
                         0, MPI_COMM_RT, ierr)
          call MPI_BCAST(Frec%Mntfreqi(0), nproc, MPI_INTEGER, &
                         0, MPI_COMM_RT, ierr)
          call MPI_BCAST(Frec%Mnpfreq(0), nproc, MPI_INTEGER, &
                         0, MPI_COMM_RT, ierr)

        end if ! Type of bcast

        ! Maximum size of profile
        MPID%nxtfreq = maxval(Frec%Mntfreq(1:nproc-1))
        MPID%nxtfreqi = maxval(Frec%Mntfreqi(1:nproc-1))
        MPID%nxpfreq = maxval(Frec%Mnpfreq(1:nproc-1))

        ! Allocate size arrays

        ! Size of RT information sent between domains
        allocate(MPID%size1(0:nproc-1))
        ! Size of profile information sent between domains
        allocate(MPID%size2(0:nproc-1))
        ! Size of location indexes vector sent to master
        allocate(MPID%size3(0:nproc-1))
        ! Size of Stokes chunk sent to master
        allocate(MPID%size4(0:nproc-1))
        ! Size of profile chunk sent to master
        allocate(MPID%size5(0:nproc-1))
        ! Size of JKQ to broadcast
        allocate(MPID%size6(0:nproc-1))
        ! Size of JKQ frequency dependent to broadcast
        allocate(MPID%size7(0:nproc-1))
        ! Size of Stokes to broadcast
        allocate(MPID%size8(0:nproc-1))
        ! Size of stokes in emergence sent to master
        allocate(MPID%size10(0:nproc-1))
        ! Size of RT information sent between domains in intensity
        allocate(MPID%sizei1(0:nproc-1))
        ! Size of Lambda_line to broadcast
        allocate(MPID%sizei2(0:nproc-1))
        ! Size of J00P to broadcast
        allocate(MPID%sizei3(0:nproc-1))
        ! Size of Stokes chunk sent to master in intensity
        allocate(MPID%sizei4(0:nproc-1))
        ! Size of Stokes chunk sent to master in intensity
        allocate(MPID%sizei4b(0:nproc-1))
        ! Size of profile chunk sent to master in intensity
        allocate(MPID%sizei5(0:nproc-1))
        ! Size of J00 to broadcast
        allocate(MPID%sizei6(0:nproc-1))
        ! Size of J00 frequency dependent to broadcast
        allocate(MPID%sizei7(0:nproc-1))
        ! Size of Stokes to broadcast in intensity
        allocate(MPID%sizei8(0:nproc-1))
        ! Size of Lambda_line sent to master
        allocate(MPID%sizei9(0:nproc-1))
        ! Size of Lamda_photo sent to master
        allocate(MPID%sizei0(0:nproc-1))
        ! Size of Lambda_photo to broadcast
        allocate(MPID%sizei10(0:nproc-1))
        ! Size of profiles between slaves
        allocate(MPID%sizei11(0:nproc-1))
        ! Size of b-b r between slaves
        allocate(MPID%sizei12(0:nproc-1))
        ! Size of b-f r between slaves
        allocate(MPID%sizei13(0:nproc-1))
        ! Size of coontribution function intensity
        allocate(MPID%sizei14(0:nproc-1))

        !
        ! Compute buffer sizes
        !

        do iproc=0,nproc-1

          ! RT data1 between domains
          MPID%size1(iproc) = MPID%nf(iproc)*4*6
          ! Profile between domains
          MPID%size2(iproc) = 2*Frec%Mntfreq(iproc)
          ! Contribution function to master
          MPID%size3(iproc) = 4*MPID%nf(iproc)*nz
          ! Stokes (full dependence) to master
          MPID%size4(iproc) = 4*MPID%nf(iproc)*Geom%nTh*Geom%nPh*nz
          ! Profile to master
          MPID%size5(iproc) = 2*Frec%Mntfreq(iproc)*Geom%nTh* &
                             Geom%nPh*nz
          ! JKQ/JKQS to broadcast
          MPID%size6(iproc) = 5*3*nxtran*nz
          ! JKQ(nu) to broadcast
          MPID%size7(iproc) = 5*3*nfreq*nz
          ! Stokes (full dependence) to broadcast
          MPID%size8(iproc) = 4*nfreq*Geom%nTh*Geom%nPh*nz
          ! Stokes (emergent) to master
          MPID%size10(iproc) = 4*MPID%nf(iproc)


          ! Lambda_photo to master
          MPID%sizei0(iproc) = nxb*Frec%Mnpfreq(iproc)* &
                               GeomI%nTh*GeomI%nPh*nz
          ! RT data1 between domains in intensity
          MPID%sizei1(iproc) = MPID%nf(iproc)*3
          ! Lambda_line to broadcast
          MPID%sizei2(iproc) = nxb*2*nxphot*nz
          ! J00P to broadcast
          MPID%sizei3(iproc) = 2*nxphot*nz
          ! Stokes (full dependence) in intensity to master
          MPID%sizei4(iproc) = MPID%nf(iproc)*GeomI%nTh*GeomI%nPh*nz
          ! Duplicate
          MPID%sizei4b(iproc) = MPID%sizei4(iproc)
          ! Profile in intensity to master
          MPID%sizei5(iproc) = 2*Frec%Mntfreqi(iproc)*GeomI%nTh* &
                              GeomI%nPh*nz
          ! J00/J00S to broadcast
          MPID%sizei6(iproc) = nxt*nz
          ! J00(nu) to broadcast
          MPID%sizei7(iproc) = nfreq*nz
          ! Stokes (full dependence) in intensity to
          ! broadcast
          MPID%sizei8(iproc) = nfreq*GeomI%nTh*GeomI%nPh*nz
          ! Lambda_line to master
          MPID%sizei9(iproc) = nxb*Frec%Mntfreqi(iproc)*GeomI%nTh* &
                              GeomI%nPh*nz
          ! Lambda_photo to broadcast
          MPID%sizei10(iproc) = nxb*nxt*nz
          ! Profiles in intensity between domains
          MPID%sizei11(iproc) = 2*Frec%Mntfreqi(iproc)
          ! b-b r between domains
          MPID%sizei12(iproc) = nxb*Frec%Mntfreqi(iproc)
          ! b-f r between domains
          MPID%sizei13(iproc) = nxb*Frec%Mnpfreq(iproc)
          ! Contribution function in intensity to master
          MPID%sizei14(iproc) = MPID%nf(iproc)*nz

        end do


        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the continuum problem
        !

        ! If doing intensity
        if (lio.and.lJ) then

          ! Compute expected buffer size
          b1 = 8d-6*MPID%nxfreq*GeomI%nTh*GeomI%nPh*nz

          ! If buffers are going to be too big
          if (b1.gt.maxbuffer) then

            ! We need to use the alternative
            MPID%alternJ = .True.

            ! Compute expected buffer size in the alternative scheme
            b1 = 8d-6*MPID%nxfreq*nz

            ! If buffers are going to be too big
            if (b1.gt.maxbuffer) then
              umsg = 'The buffer is too big even with '// &
                     'the alternative routines. Reduce the '// &
                     'number of frequencies or use more CPU'
              call aborted
              return
            end if

            ! Send warning
            if (gpid.eq.0) then
              umsg = ' # The buffers for MPI are too big. The '// &
                     'code will use the alternative '// &
                     'communication scheme for continuum.'
              call verbose
            end if

            ! For each processor
            do iproc=0,nproc-1
              MPID%sizei4b(iproc) = MPID%nf(iproc)*nz
            end do

          end if ! Buffer too big
        end if ! Going through solveJ

        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the intensity problem
        !

        ! If doing intensity
        if (lio) then

          ! Compute expected buffer size
          b1 = 8d-6*MPID%nxfreq*GeomI%nTh*GeomI%nPh*nz
          b2 = 8d-6*MPID%nxtfreqi*2*GeomI%nTh*GeomI%nPh*nz
          b3 = 8d-6*MPID%nxtfreqi*nxb*GeomI%nTh*GeomI%nPh*nz
          b4 = 8d-6*MPID%nxpfreq*nxb*GeomI%nTh*GeomI%nPh*nz

          ! If buffers are going to be too big
          if ((b1+b2+b3+b4).gt.maxbuffer) then

            ! We need to use the alternative
            MPID%alternI = .True.

            ! Compute expected buffer size in the alternative scheme
            b1 = 8d-6*MPID%nxfreq*nz
            b2 = 8d-6*MPID%nxtfreqi*2*nz
            b3 = 8d-6*MPID%nxtfreqi*nxb*nz
            b4 = 8d-6*MPID%nxpfreq*nxb*nz

            ! If buffers are going to be too big
            if ((b1+b2+b3+b4).gt.maxbuffer) then
              umsg = 'The buffer is too big even with '// &
                     'the alternative routines. Reduce the '// &
                     'number of frequencies or use more CPU'
              call aborted
              return
            end if

            ! Send warning
            if (gpid.eq.0) then
              umsg = ' # The buffers for MPI are too big. The '// &
                     'code will use the alternative '// &
                     'communication scheme for intensity.'
              call verbose
            end if

            ! For each processor
            do iproc=0,nproc-1
              MPID%sizei4(iproc) = MPID%nf(iproc)*nz
              MPID%sizei5(iproc) = 2*Frec%Mntfreqi(iproc)*nz
              MPID%sizei9(iproc) = nxb*Frec%Mntfreqi(iproc)*nz
              MPID%sizei0(iproc) = nxb*Frec%Mnpfreq(iproc)*nz
            end do

          end if ! Buffer to big
        end if ! Going through solveri


        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the polarization correction problem
        !

        ! If correcting intensity
        if (lGen) then

          ! If CRD or AV
          if (.not.PRD.or.AV) then

            ! Compute expected buffer size
            b1 = 8d-6*MPID%nxfreq*2*nz

            ! If buffers are going to be too big
            if (b1.gt.maxbuffer) then
              umsg = 'The buffer is too big. Reduce the '// &
                     'number of frequencies or use more CPU'
              call aborted
              return
            end if
          end if

          if (PRD.and..not.AV) then

            ! Compute expected buffer size
            b1 = 8d-6*MPID%nxtfreq*2*Geom%nTh*Geom%nPh*nz

            ! If buffers are going to be too big
            if (b1.gt.maxbuffer) then

              ! We need to use the alternative
              MPID%alternJgen = .True.

              ! Compute expected buffer size in the alternative scheme
              b1 = 8d-6*MPID%nxtfreq*2*nz

              ! If buffers are going to be too big
              if (b1.gt.maxbuffer) then
                umsg = 'The buffer is too big even with '// &
                       'the alternative routines. Reduce '// &
                       'the number of frequencies or use '// &
                       'more CPU'
                call aborted
                return
              end if

              ! Send warning
              if (gpid.eq.0) then
                umsg = ' # The buffers for MPI are too big. '// &
                       'The code will use the alternative '// &
                       'communication scheme for the '// &
                       'polarization correction.'
                call verbose
              end if ! Master
            end if ! Big buffer
          end if ! PRD and AD
        end if ! Doing JKQ_gen


        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the polarization problem
        !

        ! If iterating polarization
        if (lp) then

          ! Compute expected buffer size
          b1 = 8d-6*MPID%nxfreq*4*Geom%nTh*Geom%nPh*nz
          b2 = 8d-6*MPID%nxtfreq*2*Geom%nTh*Geom%nPh*nz

          ! If buffers are going to be too big
          if ((b1+b2).gt.maxbuffer) then

            ! We need to use the alternative
            MPID%alternP = .True.

            ! Compute expected buffer size in the alternative scheme
            b1 = 8d-6*MPID%nxfreq*4*nz
            b2 = 8d-6*MPID%nxtfreq*2*nz

            ! If buffers are going to be too big
            if ((b1+b2).gt.maxbuffer) then
              umsg = 'The buffer is too big even with '// &
                     'the alternative routines. Reduce the '// &
                     'number of frequencies or use more CPU'
              call aborted
              return
            end if

            ! Send warning
            if (gpid.eq.0) then
              umsg = ' # The buffers for MPI are too big. The '// &
                     'code will use the alternative '// &
                     'communication scheme for polarization.'
              call verbose
            end if

            do iproc=0,nproc-1
              MPID%size4(iproc) = 4*MPID%nf(iproc)*nz
              MPID%size5(iproc) = 2*Frec%Mntfreq(iproc)*nz
            end do

          end if
        end if

        ! Alternative bcast, wait until done
        if (MPID%altbcast) then

          ! For each slave to send
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing
            call MPI_WAIT(MPID%requestA(iproc,1), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,2), &
                          MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%requestA(iproc,3), &
                          MPI_STATUS_IGNORE,ierr)

          end do
        end if
      end if

      ! Control
      call control

      return

      end subroutine setmpi_sizes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resets request MPI variables.\n
      !!   MPID(MPI_class): Structure with MPI data
      subroutine reset_mpirequest(MPID)

      ! I/O
      type(MPI_class), intent(inout):: MPID

      ! Local
      integer:: iproc,ii,nreq

      ! Reset numbered requests
      MPID%request1 = MPI_REQUEST_NULL
      MPID%request2 = MPI_REQUEST_NULL
      MPID%request3 = MPI_REQUEST_NULL
      MPID%request4 = MPI_REQUEST_NULL
      MPID%request5 = MPI_REQUEST_NULL
      MPID%request6 = MPI_REQUEST_NULL
      MPID%request7 = MPI_REQUEST_NULL
      MPID%request8 = MPI_REQUEST_NULL
      MPID%request9 = MPI_REQUEST_NULL
      MPID%request0 = MPI_REQUEST_NULL
      MPID%request11 = MPI_REQUEST_NULL

      if (.not.MPID%mpi) return

      ! Reset array request
      if (MPID%altbcast) then
        nreq = MPID%nsend
      else
        nreq = 0
      end if
      do iproc=1,nreq
        do ii=1,8

          MPID%requestA(iproc,ii) = MPI_REQUEST_NULL

        end do
      end do

      return

      end subroutine reset_mpirequest

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a report on RAM use\n
      !!      MPID(MPI_class): Structure with MPI data\n
      !! folder(character(:)): Path to the output folder\n
      !!         pol(integer): Indicates if doing polarization\n
      !!      formal(integer): Indicates if doing formal solution
      subroutine RAMreport(MPID,folder,pol,formal)

      ! I/O
      type(MPI_class), intent(in):: MPID
      character(len=500), intent(in):: folder
      integer, intent(in):: pol,formal

      ! Local
      character(len=11):: fname
      integer:: iproc,ios
      double precision, dimension(:), allocatable:: lRAM,lPRAM,lVRAM
      double precision, dimension(:), allocatable:: lWRAM,lRRAM,lBRAM

      ! Allocate processors
      allocate(lRAM(nproc))
      allocate(lPRAM(nproc))
      allocate(lVRAM(nproc))
      allocate(lWRAM(nproc))
      allocate(lRRAM(nproc))
      allocate(lBRAM(nproc))

      ! If MPI
      if (MPID%mpi) then

        ! Tell everything to master
        call MPI_GATHER(MPID%RAM, 1, MPI_DOUBLE_PRECISION, &
                        lRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, ierr)
        call MPI_GATHER(MPID%PRAM, 1, MPI_DOUBLE_PRECISION, &
                        lPRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, ierr)
        call MPI_GATHER(MPID%VRAM, 1, MPI_DOUBLE_PRECISION, &
                        lVRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, ierr)
        call MPI_GATHER(MPID%WRAM, 1, MPI_DOUBLE_PRECISION, &
                        lWRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, ierr)
        call MPI_GATHER(MPID%RRAM, 1, MPI_DOUBLE_PRECISION, &
                        lRRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, ierr)
        call MPI_GATHER(MPID%BRAM, 1, MPI_DOUBLE_PRECISION, &
                        lBRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_WORLD, ierr)

      ! No MPI
      else

        lRAM(1) = MPID%RAM
        lPRAM(1) = MPID%PRAM
        lVRAM(1) = MPID%VRAM
        lWRAM(1) = MPID%WRAM
        lRRAM(1) = MPID%RRAM
        lBRAM(1) = MPID%BRAM

      end if ! MPI

      ! Only master writes
      if (pid.eq.0) then

        ! Choose filename
        if (pol.eq.1) then
          if (formal.eq.1) then
            fname = 'RAM-P-F.dat'
          else
            fname = 'RAM-P-E.dat'
          endif
        else
          if (formal.eq.1) then
            fname = 'RAM-I-F.dat'
          else
            fname = 'RAM-I-E.dat'
          endif
        end if

        ! Open file
        open (200,file=trim(folder)//'/'//fname, &
              status='unknown', iostat=ios, action='write')

        ! Line
        write(200,'(A)') '------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--'

        ! Header
        write(200,'(A)') '#  CPU'// &
                         ' | Background '// &
                         ' | Photoioni. '// &
                         ' |  Radiation '// &
                         ' | Voigt prf. '// &
                         ' | Redistrib. '// &
                         ' |    Total   '// &
                         ' |'

        ! Line
        write(200,'(A)') '------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--'

        ! For each CPU
        do iproc=1,nproc
          write(200,'(2x,i4,6(" | ",2x,f9.3)," |")') &
            iproc-1, &
            lBRAM(iproc),lPRAM(iproc),lRRAM(iproc), &
            lVRAM(iproc),lWRAM(iproc),lRAM(iproc)
        end do

        ! Line
        write(200,'(A)') '------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--'

        ! Close file
        close(200)

      end if ! Master

      ! Control
      call control

      return

      end subroutine RAMreport

!#####################################################################
!#####################################################################
!#####################################################################

#ifdef _OPENMP

      !> Activate OpenMP sections and work split for some sections.
      subroutine setomp

      ! Check if multiple threads
      if (nthread.gt.1) then
        omp = .True.
      else
        omp = .False.
      end if

      end subroutine setomp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Decides how the work is going to be split with OpenMP in
      !! the magnetic polarization first order RT coefficients.\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!  Bstrength(dfloat(:)): Magnetic field strength
      subroutine setomp_magn(Atom,Bstrength)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: Nfield, anyyes

      integer:: iz,ia,if0,if1,fpt,iaux,tid,iterml,itermu
      integer:: jtran,nt
      integer:: iL,iMl,iMu,iU,ncom,nMl,nMu
      integer, dimension(:), allocatable:: vnf

      double precision:: S,rLu,rLl,rJumax,rJlmax,rMl,rMu


      ! Routine name
      urou = 'setomp_magn'

      ! Initialize No field
      Nfield = .True.

      ! If MPI, the master does not care
      if (nproc.gt.1.and.pid.eq.0) return

      ! Check if there is magnetic field
      do iz=Rz0,Rz1

        ! There is magnetic field somewhere
        if (Bstrength(iz).gt.TINYB) then
          Nfield = .False.
          exit
        end if

      end do ! Checking magnetic field at every height

      ! If already did
      if (allocated(Atom(1)%omp_comp_1ord)) then

        ! All true?

        ! If no field
        if (Nfield) then

          ! Make all falses
          do ia=1,nA
            anyyes = any(Atom(ia)%omp_comp_1ord)
            if (anyyes) Atom(ia)%omp_comp_1ord = .False.
          end do

          ! And return
          return

        ! If yes field
        else

          ! Is allocated?
          if (allocated(Atom(1)%omp_1c)) then

            ! Already did, skip
            return

          ! Not allocated
          else

            ! Must have been non-magnetic, reset
            do ia=1,nA
              deallocate(Atom(ia)%omp_comp_1ord)
            end do

          end if ! Allocated?
        end if ! Field?
      end if ! Already went through here

      ! If there is no magnetic field or a single
      ! thread
      if (Nfield.or.nthread.eq.1) then

        ! Set bools to false for every atom and return
        do ia=1,nA
          allocate(Atom(ia)%omp_comp_1ord(Atom(ia)%ntran))
          Atom(ia)%omp_comp_1ord = .False.
        end do

        return

      end if ! Only one thread

      ! For each atom
      do ia=1,nA

        ! Allocate bools
        allocate(Atom(ia)%omp_comp_1ord(Atom(ia)%ntran))
        allocate(Atom(ia)%omp_1c(Atom(ia)%ntran))

        ! For each transition
        do jtran=1,Atom(ia)%ntran

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          iterml = Atom(ia)%fst(jtran)%iterml

          ! Get frequency limits
          if0 = Atom(ia)%if0(jtran)
          if1 = Atom(ia)%if1(jtran)
          nt = if1 - if0 + 1

          ! Hypothetical frequency per thread
          fpt = nt/nthread

          ! Spin
          S = Atom(ia)%Sval(itermu)

          ! Orbital angular momentum
          rLu = Atom(ia)%rLval(itermu)
          rLl = Atom(ia)%rLval(iterml)

          ! Determine the maximum angular momentum and the number
          ! of magnetic sublevels for that maximum momentum
          rJumax = rLu+S
          nMu = nint(2d0*rJumax+1d0)
          rJlmax = rLl+S
          nMl = nint(2d0*rJlmax+1d0)

          !
          ! Count components
          !

          ! Initialize counter
          ncom = 0

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! For each mu_u
            do iU=1,Atom(ia)%nblk(iMu,itermu)

              ! For each Ml
              do iMl=1,nMl

                ! Value of Ml
                rMl = -rJlmax + dble(iMl-1)

                ! If not pi nor sigma, skip
                if (nint(abs(rMu-rMl)).gt.1) cycle

                ! For each mu_l
                do iL=1,Atom(ia)%nblk(iMl,iterml)

                  ! Add count
                  ncom = ncom + 1

                end do ! mu_l
              end do ! Ml
            end do ! mu_u
          end do ! Mu

          ! If more threads than components or less frequencies than
          ! threads, split in components
          if (nthread.le.ncom.or.fpt.lt.1) then

            ! Allocate nf if not allocated
            if (.not.allocated(vnf)) allocate(vnf(nthread))

            ! Set to true
            Atom(ia)%omp_comp_1ord(jtran) = .True.

            ! Allocate limits
            allocate(Atom(ia)%omp_1c(jtran)%if0(nthread))
            allocate(Atom(ia)%omp_1c(jtran)%if1(nthread))

            !
            ! Simple split
            !

            ! Split equally first
            iaux = ncom/nthread

            ! Initialize nf
            vnf = iaux

            ! If remaining
            if (iaux*nthread.ne.ncom) then

              ! Remainder
              iaux = ncom - iaux*nthread

              ! Add to firsts processes
              do tid=1,iaux
                vnf(tid) = vnf(tid) + 1
              end do

            end if

            ! Determine component limits
            Atom(ia)%omp_1c(jtran)%if0(1) = 1
            Atom(ia)%omp_1c(jtran)%if1(1) = vnf(1)

            ! For each other thread
            do tid=2,nthread
              Atom(ia)%omp_1c(jtran)%if0(tid) = 1 + &
                                   Atom(ia)%omp_1c(jtran)%if1(tid-1)
              Atom(ia)%omp_1c(jtran)%if1(tid) = vnf(tid) - 1 + &
                                     Atom(ia)%omp_1c(jtran)%if0(tid)
            end do

          ! If not, split in frequency
          else

            ! Set to false
            Atom(ia)%omp_comp_1ord(jtran) = .False.

          end if ! Splitting in frequencies

        end do ! Output transition
      end do ! Atom

      ! Back
      return

      end subroutine setomp_magn

!#####################################################################
!#####################################################################
!#####################################################################

      !> Decides how the work is going to be split with OpenMP in
      !! the emiss2ord routine.\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !! Frec(Frequency_class): Structure with frequency data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!  Bstrength(dfloat(:)): Magnetic field strength
      subroutine setomp_magn_2ord(Atom,Flgsg,Bstrength)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: Nfield,integrate,done

      integer:: jtran,mina,maxa,minto,maxto,nt,ktran,itran
      integer:: itermu,iterml,itermf,nMu,nMf,nMl,ncom,iti,iz,ia
      integer:: iMf,mF,indF,iMu,iq,iU,indU,iMu1,iq1,iQQ,iU1,indU1
      integer:: iMl,ip,iL,indL,iMl1,ip1,iPP,iL1,indL1,icom,K1,K
      integer:: kLb,iJlb,kLb1,iJlb1,kL,iJl,kU2,iJu2
      integer:: kL1,iJl1,kU3,iJu3,kF,iJf,kU,iJu
      integer:: kF1,iJf1,kU1,iJu1,tid
      integer, dimension(:), allocatable:: nf,nnf

      double precision:: S,rLu,rLf,rLl,rJumax,rJfmax,rJlmax
      double precision:: rMf,rMu,q,rMu1,q1,QQ,rMl,p,rMl1,p1,PP,rK1,rK
      double precision:: ftmp,f1tmp,Clb,rJlb,Clb1,rJlb1,rhoc,Cl,rJl
      double precision:: Cu2,rJu2,CC2,Cl1,rJl1,Cu3,rJu3,CC3,Cf,rJf
      double precision:: Cu,rJu,CC,Cf1,rJf1,Cu1,rJu1,CC1

      complex(kind=8):: tmpK


      ! Routine name
      urou = 'setomp_magn_2ord'

      ! Initialize No field
      Nfield = .True.

      ! If MPI, the master does not care
      if (nproc.gt.1.and.pid.eq.0) return

      ! If there is no PRD
      if (.not.PRD) return

      ! If a single thread, skip
      if (nthread.eq.1) return

      ! Check if there is magnetic field
      do iz=Rz0,Rz1

        ! There is magnetic field somewhere
        if (Bstrength(iz).gt.TINYB) then
          Nfield = .False.
          exit
        end if

      end do ! Checking magnetic field at every height

      ! If there is no magnetic field, just return
      if (Nfield) return

      ! Check we did not do this already
      done = .False.

      ! For each atom
      do ia=1,nA

        ! Check allocated
        if (allocated(Atom(ia)%omp_2c)) then
          done = .True.
          exit
        end if

      end do

      ! If done already, go out
      if (done) return

      ! Initialize indexes
      mina = 10000
      maxa = 0

      ! For each atom
      do ia=1,nA
        ! For all transitions
        do jtran=1,Atom(ia)%ntran

          ! If PRD line
          if (Atom(ia)%lemiss2(jtran)) then
            ! If not absent in this CPU
            if (.not.Atom(ia)%fflag(jtran)%absent) then

              ! Update limits
              if (ia.lt.mina) mina = ia
              if (ia.gt.maxa) maxa = ia

            end if ! Presence of line
          end if ! PRD line

        end do ! Transitions
      end do ! Atoms

      ! For each atom
      do ia=mina,maxa

        !
        ! Count the transition combinations
        !

        ! Reset index
        minto = Atom(ia)%ntran + 1
        maxto = 0
        nt = 0

        ! For each upper term
        do jtran=1,Atom(ia)%ntran

          ! PRD line
          if (.not.Atom(ia)%lemiss2(jtran)) cycle
          ! Present line
          if (Atom(ia)%fflag(jtran)%absent) cycle

          ! Update size
          nt = nt + 1

          ! Update limits
          if (jtran.lt.minto) minto = jtran
          if (jtran.gt.maxto) maxto = jtran

        end do ! Transitions

        ! If none, skip
        if (nt.lt.1) cycle

        ! Allocate split data
        allocate(Atom(ia)%omp_2c(nt))

        ! For each transition
        do jtran=minto,maxto

          ! Check if skipping
          if (.not.Atom(ia)%lemiss2(jtran).or. &
              Atom(ia)%fflag(jtran)%absent) cycle

          ! Get continuous index
          ktran = Atom(ia)%itrano(jtran)

          ! Get terms and quantum numbers
          itermu = Atom(ia)%fst(jtran)%itermu
          itermf = Atom(ia)%fst(jtran)%iterml
          S = Atom(ia)%Sval(itermu)
          rLu = Atom(ia)%rLval(itermu)
          rLf = Atom(ia)%rLval(itermf)

          ! Determine the maximum angular momentum and the number
          ! of magnetic sublevels for that maximum momentum
          rJumax = rLu + S
          nMu = nint(2d0*rJumax+1d0)
          rJfmax = rLf + S
          nMf = nint(2d0*rJfmax+1d0)

          ! Get number of input transitions
          nt = Atom(ia)%trano(ktran)%nt

          ! Allocate limits
          allocate(Atom(ia)%omp_2c(ktran)%if0(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%if1(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mnU(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mxU(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mnU1(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mxU1(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%nif0(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%nif1(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mnnU(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mxnU(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mnnU1(nthread,nt))
          allocate(Atom(ia)%omp_2c(ktran)%mxnU1(nthread,nt))

          ! Initialize
          Atom(ia)%omp_2c(ktran)%mnU = 1000000
          Atom(ia)%omp_2c(ktran)%mnU1 = 1000000
          Atom(ia)%omp_2c(ktran)%mnnU = 1000000
          Atom(ia)%omp_2c(ktran)%mnnU1 = 1000000
          Atom(ia)%omp_2c(ktran)%mxU = 0
          Atom(ia)%omp_2c(ktran)%mxU1 = 0
          Atom(ia)%omp_2c(ktran)%mxnU = 0
          Atom(ia)%omp_2c(ktran)%mxnU1 = 0

          ! For each input transition
          do iti=1,nt

            ! Get index real transition
            itran = Atom(ia)%trano(ktran)%ind(iti)

            ! Get term and quantum number
            iterml = Atom(ia)%fst(itran)%iterml
            rLl = Atom(ia)%rLval(iterml)

            ! Determine maximum value of J and number of magnetic
            ! sublevels for this maximum J
            rJlmax = rLl + S
            nMl = nint(2d0*rJlmax+1d0)

            ! Get number of components
            ncom = maxval(Atom(ia)%trano(ktran)%trani(iti)%Wind)

            ! Allocate weights with and wighout coherent lower
            ! term
            allocate(nf(ncom),nnf(ncom))
            nf = 0
            nnf = 0

            !
            ! Calculate weight of each component assuming coherent
            !

      ! Reset indent

      ! For each Mf
      do iMf=1,nMf

        ! Value of Mf
        rMf = -rJfmax + dble(iMf-1)

        ! For each mu_f
        do mF=1,Atom(ia)%nblk(iMf,itermf)

          ! Get indexes
          indF = Atom(ia)%i_Wind(itermf)%ind(mF,iMf)

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! Difference between M momentums, done integer
            q = rMu - rMf
            iq = nint(q)

            ! If not pi nor sigma, skip
            if(abs(iq).gt.1) cycle

            ! For each mu_u
            do iU=1,Atom(ia)%nblk(iMu,itermu)

              ! Get indexes
              indU = Atom(ia)%i_Wind(itermu)%ind(iU,iMu)

              ! For each Mu'
              do iMu1=1,nMu

                ! Value of Mu'
                rMu1 = -rJumax + dble(iMu1-1)

                ! Difference between M momentums
                q1 = rMu1-rMf
                QQ = q1-q

                ! Convert to integers
                iq1 = nint(q1)
                iQQ = nint(QQ)

                ! If not pi or sigma, skip
                if(abs(iq1).gt.1) cycle

                ! For each mu_u'
                do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                  ! Get indexes
                  indU1 = Atom(ia)%i_Wind(itermu)%ind(iU1,iMu1)

                  ! For each Ml
                  do iMl=1,nMl

                    ! Value of Ml
                    rMl = -rJlmax + dble(iMl-1)

                    ! Difference between M momentums, in integer
                    p = rMu-rMl
                    ip = nint(p)

                    ! If not pi nor sigma, skip
                    if(abs(ip).gt.1) cycle

                    ! For each mu_l
                    do iL=1,Atom(ia)%nblk(iMl,iterml)

                      ! Get indexes
                      indL = Atom(ia)%i_Wind(iterml)%ind(iL,iMl)

                      ! For each Ml'
                      do iMl1=1,nMl

                        ! Value of Ml'
                        rMl1 = -rJlmax + dble(iMl1-1)

                        ! Difference between M momentums
                        p1 = rMu1-rMl1
                        PP = p1-p

                        ! Convert to integer
                        ip1 = nint(p1)
                        iPP = nint(PP)

                        ! If not pi nor sigma, skip
                        if(abs(ip1).gt.1) cycle

                        ! For each mu_l'
                        do iL1 = 1,Atom(ia)%nblk(iMl1,iterml)

                          ! Get indexes
                          indL1 = Atom(ia)%i_Wind(iterml)%ind(iL1,iMl1)

      !
      ! Reset indenting
      !

      ! Get the component index
      icom = Atom(ia)%trano(ktran)%trani(iti)% &
                  Wind(indL1,indL,indF,indU1,indU)

      ! For each K'
      do K1=abs(iPP),2

        ! Get real value
        rK1 = dble(K1)

        ! Integrate
        integrate = .True.

        ! For each K
        do K=abs(iQQ),2

          ! Get real value
          rK = dble(K)

          ! Racah algebra
          ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

          ! If forbidden (3j-sym=0), skip
          if(abs(ftmp).lt.TINYJS) cycle

          ! Racah algebra
          f1tmp = fun3j(1d0,1d0,rK1,-p,p1,-PP,Flgsg)

          ! If forbidden (3j-sym=0), skip
          if(abs(f1tmp).lt.TINYJS) cycle

          ! Add the rest of the factor
          ftmp = ftmp*Flgsg%sg(iq)*sqrt(2d0*rK+1d0)
          f1tmp = f1tmp*Flgsg%sg(ip1)*sqrt(2d0*rK1+1d0)

          ! Reset temporal variable
          tmpK = cZero

          ! For each Jlb
          do kLb=1,Atom(ia)%nblk(iMl,iterml)

            ! Get eigenvector for lower b level
            Clb = maxval(abs(Atom(ia)%evec(kLb,iL,iMl,iterml,:)))

            ! If coefficient too small, skip
            if(abs(Clb).lt.TINYEV) cycle

            ! Get J level index
            iJlb = Atom(ia)%iJval(kLb,iMl,iterml)

            ! Get angular momentum
            rJlb = Atom(ia)%rJval(iJlb,iterml)

            ! For each Jlb'
            do kLb1=1,Atom(ia)%nblk(iMl1,iterml)

              ! Get eigenvector for lower b level
              Clb1 = &
                    maxval(abs(Atom(ia)%evec(kLb1,iL1,iMl1,iterml,:)))

              ! If coefficient too small, skip
              if(abs(Clb1).lt.TINYEV) cycle

              ! Get J level index
              iJlb1 = Atom(ia)%iJval(kLb1,iMl1,iterml)

              ! Get angular momentum
              rJlb1 = Atom(ia)%rJval(iJlb1,iterml)

              ! Add coefficients to rhoKQ (assummed 1d0)
              rhoc = Flgsg%sg(nint(rJlb-rMl))*Clb*Clb1

              ! For each Jl
              do kL=1,Atom(ia)%nblk(iMl,iterml)

                ! Get eigenvector for lower level
                Cl = maxval(abs(Atom(ia)%evec(kL,iL,iMl,iterml,:)))

                ! If coefficient too small, skip
                if(abs(Cl).lt.TINYEV) cycle

                ! Get J level index
                iJl = Atom(ia)%iJval(kL,iMl,iterml)

                ! Get angular momentum
                rJl = Atom(ia)%rJval(iJl,iterml)

                ! For each Ju''
                do kU2=1,Atom(ia)%nblk(iMu,itermu)

                  ! Get eigenvector for upper'' level
                  Cu2 = &
                       maxval(abs(Atom(ia)%evec(kU2,iU,iMu,itermu,:)))

                  ! If coefficient too small, skip
                  if(abs(Cu2).lt.TINYEV) cycle

                  ! Get J level index
                  iJu2 = Atom(ia)%iJval(kU2,iMu,itermu)

                  ! Get angular momentum
                  rJu2 = Atom(ia)%rJval(iJu2,itermu)

                  ! Add coefficient to dipole strength
                  CC2 = Cl*Cu2*Atom(ia)%rdip(itran)% &
                               rdip(ip,iMu,iMl,iJu2,iJl)

                  ! If coefficient too small, skip
                  if(abs(CC2).lt.TINYCO) cycle

                  ! For each Jl'
                  do kL1=1,Atom(ia)%nblk(iMl1,iterml)

                    ! Get eigenvector for lower' level
                    Cl1 = &
                     maxval(abs(Atom(ia)%evec(kL1,iL1,iMl1,iterml,:)))

                    ! If coefficient too small, skip
                    if(abs(Cl1).lt.TINYEV) cycle

                    ! Get J level index
                    iJl1 = Atom(ia)%iJval(kL1,iMl1,iterml)

                    ! Get angular momentum
                    rJl1 = Atom(ia)%rJval(iJl1,iterml)

                    ! For each Ju'''
                    do kU3=1,Atom(ia)%nblk(iMu1,itermu)

                      ! Get eigenvector for upper''' level
                      Cu3 = &
                     maxval(abs(Atom(ia)%evec(kU3,iU1,iMu1,itermu,:)))

                      ! If coefficient too small, skip
                      if(abs(Cu3).lt.TINYEV) cycle

                      ! Get J level index
                      iJu3 = Atom(ia)%iJval(kU3,iMu1,itermu)

                      ! Get angular momentum
                      rJu3 = Atom(ia)%rJval(iJu3,itermu)

                      ! Add coefficients to dipole strength
                      CC3 = Cl1*Cu3* &
                            Atom(ia)%rdip(itran)% &
                                     rdip(ip1,iMu1,iMl1,iJu3,iJl1)

                      ! If coefficient too small, skip
                      if(abs(CC3).lt.TINYCO) cycle

      !
      ! Reset identation
      !

      ! For each Jf
      do kF=1,Atom(ia)%nblk(iMf,itermf)

        ! Get eigenvector for final level
        Cf = maxval(abs(Atom(ia)%evec(kF,mF,iMf,itermf,:)))

        ! If coefficient too small, skip
        if(abs(Cf).lt.TINYEV) cycle

        ! Get J level index
        iJf = Atom(ia)%iJval(kF,iMf,itermf)

        ! Get angular momentum
        rJf = Atom(ia)%rJval(iJf,itermf)

        ! For each Ju
        do kU=1,Atom(ia)%nblk(iMu,itermu)

          ! Get eigenvector for upper level
          Cu = maxval(abs(Atom(ia)%evec(kU,iU,iMu,itermu,:)))

          ! If coefficient too small, skip
          if(abs(Cu).lt.TINYEV) cycle

          ! Get J level index
          iJu = Atom(ia)%iJval(kU,iMu,itermu)

          ! Get angular momentum
          rJu = Atom(ia)%rJval(iJu,itermu)

          ! Add coefficients to dipole strength
          CC = Cf*Cu*Atom(ia)%rdip(jtran)% &
                              rdip(iq,iMu,iMf,iJu,iJf)

          ! If coefficient too small, skip
          if(abs(CC).lt.TINYCO) cycle

          ! For each Jf'
          do kF1=1,Atom(ia)%nblk(iMf,itermf)

            ! Get eigenvector for final' level
            Cf1 = maxval(abs(Atom(ia)%evec(kF1,mF,iMf,itermf,:)))

            ! If coefficient too small, skip
            if(abs(Cf1).lt.TINYEV) cycle

            ! Get J level index
            iJf1 = Atom(ia)%iJval(kF1,iMf,itermf)

            ! Get angular momentum
            rJf1 = Atom(ia)%rJval(iJf1,itermf)

            ! For each Ju'
            do kU1=1,Atom(ia)%nblk(iMu1,itermu)

              ! Get eigenvector for upper' level
              Cu1 = maxval(abs(Atom(ia)%evec(kU1,iU1,iMu1,itermu,:)))

              ! If coefficient too small, skip
              if(abs(Cu1).lt.TINYEV) cycle

              ! Get J level index
              iJu1 = Atom(ia)%iJval(kU1,iMu1,itermu)

              ! Get angular momentum
              rJu1 = Atom(ia)%rJval(iJu1,itermu)

              ! Add coefficients to dipole strength
              CC1 = Cf1*Cu1*Atom(ia)%rdip(itran)% &
                                     rdip(iq1,iMu1,iMf,iJu1,iJf1)

              ! If coefficient big enough, add contribution to
              ! temporal variable
              if(abs(CC1).gt.TINYCO) &
                tmpK = f1tmp*ftmp*CC*CC1*CC2*CC3*rhoc + tmpK

            end do ! kU1
          end do ! kF1
        end do ! kU
      end do ! kF
                    end do ! kU3
                  end do ! kL1
                end do ! kU2
              end do ! kL
            end do ! kLb1
          end do ! kLb

          ! Integration here
          if (integrate) then

            ! Add weight to always active
            nf(icom) = nf(icom) + 1

            ! Add weight to NCHLT
            if (iMl.eq.iMl1.and.iL.eq.iL1) &
              nnf(icom) = nnf(icom) + 1

          end if

        end do ! K1
      end do ! K

                          !
                          ! Restore indent
                          !

                        end do ! For each mu_l'
                      end do ! For each Ml'
                    end do ! For each mu_l
                  end do ! For each Ml
                end do ! For each mu_u'
              end do ! For each Mu'
            end do ! For each mu_u
          end do ! For each Mu
        end do ! For each mu_f
      end do ! For each Mf

      !
      !
      ! Now distribute indexes
      !
      !

      ! Set indexes CHLT
      call set_index(ncom,nf,sum(nf), &
                     Atom(ia)%omp_2c(ktran)%if0(:,iti), &
                     Atom(ia)%omp_2c(ktran)%if1(:,iti))

      ! Set indexes NCHLT
      call set_index(ncom,nnf,sum(nnf), &
                     Atom(ia)%omp_2c(ktran)%nif0(:,iti), &
                     Atom(ia)%omp_2c(ktran)%nif1(:,iti))

      !
      ! Set maximum and minimum U
      !

      ! For each Mf
      do iMf=1,nMf

        ! Value of Mf
        rMf = -rJfmax + dble(iMf-1)

        ! For each mu_f
        do mF=1,Atom(ia)%nblk(iMf,itermf)

          ! Get indexes
          indF = Atom(ia)%i_Wind(itermf)%ind(mF,iMf)

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! Difference between M momentums, done integer
            q = rMu - rMf
            iq = nint(q)

            ! If not pi nor sigma, skip
            if(abs(iq).gt.1) cycle

            ! For each mu_u
            do iU=1,Atom(ia)%nblk(iMu,itermu)

              ! Get indexes
              indU = Atom(ia)%i_Wind(itermu)%ind(iU,iMu)

              ! For each Mu'
              do iMu1=1,nMu

                ! Value of Mu'
                rMu1 = -rJumax + dble(iMu1-1)

                ! Difference between M momentums
                q1 = rMu1-rMf
                QQ = q1-q

                ! Convert to integers
                iq1 = nint(q1)
                iQQ = nint(QQ)

                ! If not pi or sigma, skip
                if(abs(iq1).gt.1) cycle

                ! For each mu_u'
                do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                  ! Get indexes
                  indU1 = Atom(ia)%i_Wind(itermu)%ind(iU1,iMu1)

                  ! For each Ml
                  do iMl=1,nMl

                    ! Value of Ml
                    rMl = -rJlmax + dble(iMl-1)

                    ! Difference between M momentums, in integer
                    p = rMu-rMl
                    ip = nint(p)

                    ! If not pi nor sigma, skip
                    if(abs(ip).gt.1) cycle

                    ! For each mu_l
                    do iL=1,Atom(ia)%nblk(iMl,iterml)

                      ! Get indexes
                      indL = Atom(ia)%i_Wind(iterml)%ind(iL,iMl)

                      ! For each Ml'
                      do iMl1=1,nMl

                        ! Value of Ml'
                        rMl1 = -rJlmax + dble(iMl1-1)

                        ! Difference between M momentums
                        p1 = rMu1-rMl1
                        PP = p1-p

                        ! Convert to integer
                        ip1 = nint(p1)
                        iPP = nint(PP)

                        ! If not pi nor sigma, skip
                        if(abs(ip1).gt.1) cycle

                        ! For each mu_l'
                        do iL1 = 1,Atom(ia)%nblk(iMl1,iterml)

                          ! Get indexes
                          indL1 = Atom(ia)%i_Wind(iterml)% &
                                           ind(iL1,iMl1)

      !
      ! Reset indenting
      !

      ! Get the component index
      icom = Atom(ia)%trano(ktran)%trani(iti)% &
                  Wind(indL1,indL,indF,indU1,indU)

      ! For each thread
      do tid=1,nthread

        ! Check if in range CHLT
        if (icom.ge.Atom(ia)%omp_2c(ktran)%if0(tid,iti).and. &
            icom.le.Atom(ia)%omp_2c(ktran)%if1(tid,iti)) then

          ! Update limits
          if (iU.lt.Atom(ia)%omp_2c(ktran)%mnU(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mnU(tid,iti) = iU
          if (iU.gt.Atom(ia)%omp_2c(ktran)%mxU(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mxU(tid,iti) = iU
          if (iU1.lt.Atom(ia)%omp_2c(ktran)%mnU1(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mnU1(tid,iti) = iU1
          if (iU1.gt.Atom(ia)%omp_2c(ktran)%mxU1(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mxU1(tid,iti) = iU1
        end if

        ! Check if in range NCHLT
        if (icom.ge.Atom(ia)%omp_2c(ktran)%nif0(tid,iti).and. &
            icom.le.Atom(ia)%omp_2c(ktran)%nif1(tid,iti)) then

          ! Update limits
          if (iU.lt.Atom(ia)%omp_2c(ktran)%mnnU(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mnnU(tid,iti) = iU
          if (iU.gt.Atom(ia)%omp_2c(ktran)%mxnU(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mxnU(tid,iti) = iU
          if (iU1.lt.Atom(ia)%omp_2c(ktran)%mnnU1(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mnnU1(tid,iti) = iU1
          if (iU1.gt.Atom(ia)%omp_2c(ktran)%mxnU1(tid,iti)) &
            Atom(ia)%omp_2c(ktran)%mxnU1(tid,iti) = iU1
        end if

      end do ! Threads

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

            !
            ! Restore indent
            !

            ! Deallocate weights
            deallocate(nf,nnf)

          end do ! Input transitions
        end do ! Output transitions
      end do ! Atom

      ! Back
      return

      end subroutine setomp_magn_2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Decides how the work is going to be split with OpenMP in
      !! the second order emissivity.\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !! Frec(Frequency_class): Structure with frequency data\n
      !!  Geom(Geometry_class): Structure with geometry data\n
      !!  Bstrength(dfloat(:)): Magnetic field strength\n
      !! polarization(logical): If organizing for polarization\n
      !!    emergence(logical): If emergent calculation
      subroutine setomp_2ord(Atom,Frec,Geom,Bstrength, &
                             polarization,emergence)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(inout):: Frec
      type(Geometry_class), intent(in):: Geom
      logical, intent(in):: polarization,emergence
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: lAV
      logical, dimension(:), allocatable:: FAset
      logical, dimension(nthread):: set

      integer:: indx,jndx,jj,ii,kk,tid,nn,MW,TW,lif1,rif0
      integer:: jdir,ia,jtran,iz,mina,maxa,minto,maxto
      integer, dimension(:), allocatable:: nf,vmfreq

      double precision:: dMW,dW

      ! Pointers
      type(Frequencyc2_class), pointer:: pd

      ! Routine name
      urou = 'setomp_2ord'

      if (polarization) then
        lAV = AV
      else
        lAV = AVI
      end if

      ! If doing LOS and no need to recompute
      if (lAV.and..not.dyn.and.emergence) then
        call control
        return
      end if

      ! If a single thread, skip
      if (nthread.eq.1) then
        call control
        return
      end if

      ! Only if Frequency exists
      if (.not.associated(Frec%dzao)) then
        call control
        return
      end if

      ! Allocate Fset
      allocate(FAset(Frec%ndzao))
      FAset = .False.

      ! If polarization
      if (polarization) then

        ! Initialize indexes
        mina = 10000
        maxa = 0
        minto = 10000
        maxto = 0

        ! For each atom
        do ia=1,nA
          ! For all transitions
          do jtran=1,Atom(ia)%ntran

            ! If PRD line
            if (Atom(ia)%lemiss2(jtran)) then
              ! If not absent in this CPU
              if (.not.Atom(ia)%fflag(jtran)%absent) then

                ! Update limits
                if (jtran.lt.minto) minto = jtran
                if (jtran.gt.maxto) maxto = jtran
                if (ia.lt.mina) mina = ia
                if (ia.gt.maxa) maxa = ia

              end if ! Presence of line
            end if ! PRD line

          end do ! Transitions
        end do ! Atoms

        ! Run indexes
        do jdir=1,Frec%ndir
        do iz=Rz0,Rz1

          ! If magnetic field, skip
          if (Bstrength(iz).gt.TINYB) cycle

        do ia=mina,maxa
        do jtran=minto,maxto

          ! Get index
          indx = Frec%indx(jtran,ia,iz,jdir)

          ! If not valid, skip
          if (indx.lt.1) cycle

          ! If already set, skip
          if (FAset(indx)) cycle
          FAset(indx) = .True.

          ! Pointer
          pd => Frec%dzao(indx)

          ! Initialize maximum weight
          MW = -1
          ii = 0
          nn = 1
          allocate(vmfreq(nn))

          ! For each input transition
          do jndx=1,size(pd%trani)

            ! If allocated mfreq
            if (.not.allocated(pd%trani(jndx)%mfreq)) cycle

            ! Total weight
            TW = sum(pd%trani(jndx)%mfreq)

            ! Larger weight
            if (TW.gt.MW) then
              MW = TW

              ! If angle-averaged
              if (lAV) then
                ii = size(pd%trani(jndx)%mfreq)
              ! If angle dependent
              else
                ii = size(pd%trani(jndx)%mfreq)*Geom%nTh*Geom%nPh2
              end if

              ! Different size
              if (ii.ne.nn) then
                nn = ii
                deallocate(vmfreq)
                allocate(vmfreq(nn))
              end if

              ! If angle-averaged
              if (lAV) then
                vmfreq = pd%trani(jndx)%mfreq
              ! If angle dependent
              else
                ! Get size of output frequency
                ii = size(pd%trani(jndx)%mfreq)
                ! Get angular dimension
                tid = Geom%nTh*Geom%nPh2
                ! For each angle
                do jj=1,tid
                  kk = (jj-1)*ii
                  vmfreq(kk+1:kk+ii) = pd%trani(jndx)%mfreq
                end do
              end if
            end if

          end do ! For each input transition

          ! Allocate given threads
          allocate(pd%oif0(nthread))
          allocate(pd%oif1(nthread))

          ! Allocate size
          if (.not.allocated(nf)) allocate(nf(nthread))

          ! Modify weight if AD
          if (.not.lAV) MW = MW*Geom%nTh*Geom%nph2

          ! Average weight per CPU
          dMW = dble(MW)/dble(nthread)
          dW = dble(MW)/dble(nn)

          ! Modify average
          dMW = dMW + dW

          ! Initialize set
          set = .False.

          ! First thread
          tid = 1
          jj = 1
          pd%oif0(tid) = jj
          pd%oif1(tid) = jj
          nf(tid) = 1
          TW = vmfreq(jj)
          set(tid) = .True.

          ! Check if adding more
          if ((TW-dMW).lt.-TINYSP) then
            do while (.True.)
              if ((TW+vmfreq(jj+1)-dMW).ge.TINYSP) exit
              jj = jj + 1
              nf(tid) = nf(tid) + 1
              pd%oif1(tid) = pd%oif1(tid) + 1
              TW = TW + vmfreq(jj)
              if ((TW-dMW).ge.TINYSP) exit
              if (jj.ge.nn) exit
            end do
          end if
          lif1 = jj

          ! If more than two
          if (nthread.gt.2) then

            ! Last thread
            tid = nthread
            jj = nn
            pd%oif0(tid) = jj
            pd%oif1(tid) = jj
            nf(tid) = 1
            TW = vmfreq(jj)
            set(tid) = .True.

            ! Check if adding more
            if ((TW-dMW).lt.-TINYSP) then
              do while (.True.)
                if ((TW+vmfreq(jj-1)-dMW).ge.TINYSP) exit
                jj = jj - 1
                nf(tid) = nf(tid) + 1
                pd%oif0(tid) = pd%oif0(tid) - 1
                TW = TW + vmfreq(jj)
                if ((TW-dMW).ge.TINYSP) exit
                if (jj.ge.nn) exit
              end do
            end if
            rif0 = jj

            ! Go alternating advance
            do ii=1,nthread-1

              ! For each direction
              do kk=-1,1,2

                ! Left
                if (kk.gt.0) then

                  tid = 1 + ii

                  if (set(tid)) cycle
                  set(tid) = .True.

                  pd%oif0(tid) = pd%oif1(tid-1) + 1

                  if (set(tid+1)) then
                    pd%oif1(tid) = pd%oif0(tid+1) - 1
                    lif1 = pd%oif1(tid)
                    cycle
                  end if

                  jj = pd%oif0(tid)
                  pd%oif1(tid) = jj
                  nf(tid) = 1
                  TW = vmfreq(jj)

                  ! Check if adding more
                  if ((TW-dMW).lt.-TINYSP) then
                    do while (.True.)
                      if ((TW+vmfreq(jj+1)-dMW).ge.TINYSP) exit
                      jj = jj + 1
                      nf(tid) = nf(tid) + 1
                      pd%oif1(tid) = pd%oif1(tid) + 1
                      TW = TW + vmfreq(jj)
                      if ((TW-dMW).ge.TINYSP) exit
                      if (jj.ge.rif0-1) exit
                    end do
                  end if
                  lif1 = jj

                else

                  tid = nthread - ii

                  if (set(tid)) cycle

                  pd%oif1(tid) = pd%oif0(tid+1) - 1
                  set(tid) = .True.

                  if (set(tid-1)) then
                    pd%oif0(tid) = pd%oif1(tid-1) + 1
                    rif0 = pd%oif0(tid)
                    cycle
                  end if

                  jj = pd%oif1(tid)
                  pd%oif0(tid) = jj
                  nf(tid) = 1
                  TW = vmfreq(jj)

                  ! Check if adding more
                  if ((TW-dMW).lt.-TINYSP) then
                    do while (.True.)
                      if ((TW+vmfreq(jj-1)-dMW).ge.TINYSP) exit
                      jj = jj - 1
                      nf(tid) = nf(tid) + 1
                      pd%oif0(tid) = pd%oif0(tid) - 1
                      TW = TW + vmfreq(jj)
                      if ((TW-dMW).ge.TINYSP) exit
                      if (jj.le.lif1+1) exit
                    end do
                  end if
                  rif0 = jj
                end if

              end do
            end do

          else

            ! Last thread
            tid = nthread
            jj = nn
            pd%oif0(tid) = pd%oif1(1) + 1
            pd%oif1(tid) = jj

          end if

          deallocate(vmfreq)

        end do ! Frec%dzao indexes
        end do
        end do
        end do

        ! Allocate given threads
        if (allocated(nf)) deallocate(nf)

      ! Intensity case
      else

        ! Look for largest weights
        do indx=1,Frec%ndzao

          ! If already set, skip
          if (FAset(indx)) cycle
          FAset(indx) = .True.

          ! Pointer
          pd => Frec%dzao(indx)

          MW = -1
          ii = 0
          nn = 1
          allocate(vmfreq(nn))

          ! For each input transition
          do jndx=1,size(pd%trani)

            ! If allocated mfreq
            if (.not.allocated(pd%trani(jndx)%mfreq)) cycle

            ! Total weight
            TW = sum(pd%trani(jndx)%mfreq)

            ! Larger weight
            if (TW.gt.MW) then
              MW = TW

              ! If angle-averaged
              if (lAV) then
                ii = size(pd%trani(jndx)%mfreq)
              ! If angle dependent
              else
                ii = size(pd%trani(jndx)%mfreq)*Geom%nTh*Geom%nPh2
              end if

              ! Different size
              if (ii.ne.nn) then
                nn = ii
                deallocate(vmfreq)
                allocate(vmfreq(nn))
              end if

              ! If angle-averaged
              if (lAV) then
                vmfreq = pd%trani(jndx)%mfreq
              ! If angle dependent
              else
                ! Get size of output frequency
                ii = size(pd%trani(jndx)%mfreq)
                ! Get angular dimension
                tid = Geom%nTh*Geom%nPh2
                ! For each angle
                do jj=1,tid
                  kk = (jj-1)*ii
                  vmfreq(kk+1:kk+ii) = pd%trani(jndx)%mfreq
                end do
              end if
            end if

          end do ! For each input transition

          ! Allocate given threads
          allocate(pd%oif0(nthread))
          allocate(pd%oif1(nthread))

          ! Allocate size
          if (.not.allocated(nf)) allocate(nf(nthread))

          ! Modify weight if AD
          if (.not.lAV) MW = MW*Geom%nTh*Geom%nph2

          ! Distribute indexes
          call set_index(nn,vmfreq,MW,pd%oif0,pd%oif1)

          ! Deallocate weight vector
          deallocate(vmfreq)

        end do ! Frec%dzao indexes

        ! Allocate given threads
        if (allocated(nf)) deallocate(nf)

      end if ! Polarization vs Intensity

      ! Control and back
      call control

      ! Nullify pointer
      if (associated(pd)) nullify(pd)

      return

      end subroutine setomp_2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Distribute weights via indexes.\n
      !!      nn(integer): Size of weight vector\n
      !!      MW(integer): Maximum weight to compute averages\n
      !!   WW(integer(:)): Weight vector\n
      !!  if0(integer(:)): Initial indexes to compute\n
      !!  if1(integer(:)): Final indexes to compute
      subroutine set_index(nn,WW,MW,if0,if1)

      ! I/O
      integer, intent(in):: nn,MW
      integer, dimension(:), intent(in):: WW
      integer, dimension(:), intent(out):: if0,if1

      ! Local

      logical, dimension(nthread):: set

      integer:: tid,jj,TW,lif1,rif0,ii,kk
      integer, dimension(nthread):: nf

      double precision:: dMW,dW


      ! Average weight per CPU
      dMW = dble(MW)/dble(nthread)
      dW = dble(MW)/dble(nn)

      ! Modify average
      dMW = dMW + dW

      ! Initialize set
      set = .False.

      ! First thread
      tid = 1
      jj = 1
      if0(tid) = jj
      if1(tid) = jj
      nf(tid) = 1
      TW = WW(jj)
      set(tid) = .True.

      ! Check if adding more
      if ((TW-dMW).lt.-TINYSP) then
        do while (.True.)
          if ((TW+WW(jj+1)-dMW).ge.TINYSP) exit
          jj = jj + 1
          nf(tid) = nf(tid) + 1
          if1(tid) = if1(tid) + 1
          TW = TW + WW(jj)
          if ((TW-dMW).ge.TINYSP) exit
          if (jj.ge.nn) exit
        end do
      end if
      lif1 = jj

      ! If more than two
      if (nthread.gt.2) then

        ! Last thread
        tid = nthread
        jj = nn
        if0(tid) = jj
        if1(tid) = jj
        nf(tid) = 1
        TW = WW(jj)
        set(tid) = .True.

        ! Check if adding more
        if ((TW-dMW).lt.-TINYSP) then
          do while (.True.)
            if ((TW+WW(jj-1)-dMW).ge.TINYSP) exit
            jj = jj - 1
            nf(tid) = nf(tid) + 1
            if0(tid) = if0(tid) - 1
            TW = TW + WW(jj)
            if ((TW-dMW).ge.TINYSP) exit
            if (jj.ge.nn) exit
          end do
        end if
        rif0 = jj

        ! Go alternating advance
        do ii=1,nthread-1
          ! For each direction
          do kk=-1,1,2

            ! Left
            if (kk.gt.0) then

              tid = 1 + ii

              if (set(tid)) cycle
              set(tid) = .True.

              if0(tid) = if1(tid-1) + 1

              if (set(tid+1)) then
                if1(tid) = if0(tid+1) - 1
                lif1 = if1(tid)
                cycle
              end if

              jj = if0(tid)
              if1(tid) = jj
              nf(tid) = 1
              TW = WW(jj)

              ! Check if adding more
              if ((TW-dMW).lt.-TINYSP) then
                do while (.True.)
                  if ((TW+WW(jj+1)-dMW).ge.TINYSP) exit
                  jj = jj + 1
                  nf(tid) = nf(tid) + 1
                  if1(tid) = if1(tid) + 1
                  TW = TW + WW(jj)
                  if ((TW-dMW).ge.TINYSP) exit
                  if (jj.ge.rif0-1) exit
                end do
              end if
              lif1 = jj

            ! Right
            else

              tid = nthread - ii

              if (set(tid)) cycle
              set(tid) = .True.

              if1(tid) = if0(tid+1) - 1

              if (set(tid-1)) then
                if0(tid) = if1(tid-1) + 1
                rif0 = if0(tid)
                cycle
              end if

              jj = if1(tid)
              if0(tid) = jj
              nf(tid) = 1
              TW = WW(jj)

              ! Check if adding more
              if ((TW-dMW).lt.-TINYSP) then
                do while (.True.)
                  if ((TW+WW(jj-1)-dMW).ge.TINYSP) exit
                  jj = jj - 1
                  nf(tid) = nf(tid) + 1
                  if0(tid) = if0(tid) - 1
                  TW = TW + WW(jj)
                  if ((TW-dMW).ge.TINYSP) exit
                  if (jj.le.lif1+1) exit
                end do
              end if
              rif0 = jj
            end if

          end do
        end do

      else

        ! Last thread
        tid = nthread
        jj = nn
        if0(tid) = if1(1) + 1
        if1(tid) = jj

      end if

      end subroutine set_index
#endif

!#####################################################################
!#####################################################################
!#####################################################################

      end module setmpi_mod
