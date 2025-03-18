      !> Distribution of MPI tasks
      module setmpi_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     18/03/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/03/2025:    V4.0.1 - Implemented the option to skip ALI for
!                             photoionizations, changing the sizes
!                             for the relevant buffers (TdPA)
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
!  setmpi15D
!    Split tasks in groups of CPU, with each group taking care of
!  independent solutions of the radiation transfer problem
!
!  setmpi
!    Split tasks in the radiation transfer problem
!
!  setmpi_CLE
!    Split tasks in the radiation transfer problem in the CLE
!  synthesis
!
!  weighted_split
!    Split the ranges of a vector of positive numbers to minimize
!  the maximum sum in the individual ranges
!
!  adjust_IW
!    Read weigths and times from a previous run to try readjust the
!  expected weights
!
!  setmpi_sizes
!    Calculate the size of the buffers in MPI messages
!
!  reset_mpirequest
!    Set all request variables to MPI_REQUEST_NULL
!
!  RAMreport
!    Write the amomunt of RAM accounted for in an ASCII file
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use cram_mod
      use funnj_mod
      use parameters_mod , only : TINYSP , TINYB , cZero , TINYCO , &
                                  TINYEV , TINYJS
      use types_mod

      ! Parameters

      ! Maximum size of solver buffer in MB depending on the
      ! type of computer
      double precision, parameter:: maxbuffer_personal = 5000
      double precision, parameter:: maxbuffer_supercom = 900

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Split tasks in groups of CPU, with each group taking care of
      !! independent solutions of the radiation transfer problem\n
      !!     MPID(MPI_class): Structure with MPI data\n
      !!  Input(Input_class): Structure with configuration data
      subroutine setmpi15D(MPID,Input)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(inout):: Input

      ! Local

      integer:: nslave,ngroup,indx,jndx,lid,icolor
      integer, dimension(:), allocatable:: nelem
      integer, dimension(0:gnproc-1):: colors

#ifdef DEBUGSYN
      ! Only one group doing synthesis
      Input%rt_group_n = gnproc - 1
#endif


      ! If only 1 CPU
      if (gnproc.eq.1) then

        ! Write message depending on the run mode
        if (run_mode.eq.-1) then
          umsg = ' # It is recommended to use MPI for the 1.5D '// &
                 'inversion'
        else
          umsg = ' # It is recommended to use MPI for the 1.5D '// &
                 'synthesis'
        end if

        ! Warn this may not be a good idea
        call verbose

        ! Trivial
        pid = gpid
        nproc = gnproc
        MPID%pid = pid
        MPID%nproc = nproc
        MPID%mpi15d = .False.
        MPID%ngroup = 1

        ! Define communicator
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, gpid, &
                            MPI_COMM_RT, ierr)
        call MPI_COMM_RANK(MPI_COMM_RT,pid,ierr)
        call MPI_COMM_SIZE(MPI_COMM_RT,nproc,ierr)

        ! We can do MPI
      else

        ! Number of slaves
        nslave = gnproc - 1

        ! Master by itself
        colors(0) = 0

        ! Control expectations
        if (Input%rt_group_n.gt.nslave) then

          ! Limit RT group to the number of slaves
          Input%rt_group_n = nslave

          ! Global master
          if (gpid.eq.0) then

            ! Verbose
            write(umsg,'(A,i5)') &
              ' # You choose too many RT slaves for the '// &
              'number of processes, set to ',nslave
            call verbose

          end if ! Global master
        end if ! Asking for too many RT processes

        ! Number of groups
        ngroup = nslave/Input%rt_group_n
        MPID%ngroup = ngroup

        ! Allocate number of elements per group
        allocate(nelem(ngroup))
        nelem = Input%rt_group_n

        ! Add elements until everything is distributed
        indx = ngroup
        do while (sum(nelem).lt.nslave)
          if (indx.lt.1) indx = ngroup
          nelem(indx) = nelem(indx) + 1
          indx = indx-1
        end do

        ! Global master
        if (gpid.eq.0) then

          ! Allocate leaders
          allocate(MPID%ltslave(ngroup))
          MRAMc = MRAMc + 1d-6*sizeof(MPID%ltslave)

        end if ! Global master

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

          ! Global master save leader
          if (gpid.eq.0) MPID%ltslave(indx) = lid

          ! Rest of members of the group
          do jndx=2,nelem(indx)

            ! Advance who and give color
            lid = lid + 1
            colors(lid) = icolor

          end do ! Members
        end do ! Groups

        ! Define new communicator
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

      !> Split tasks in the radiation transfer problem\n
      !!      MPID(MPI_class): Structure with MPI data\n
      !!   Input(Input_class): Structure with configuration data\n
      !!  IW_freq(integer(:)): Expected work load for each frequency
      !!                       node
      subroutine setmpi(MPID,Input,IW_freq)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(in):: Input
      integer, dimension(:), allocatable, intent(inout):: IW_freq

      ! Local

      integer:: ios,nnd,iaux,iproc,maxcoun,ifbl,Mwg


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


        !
        ! Allocate variables
        !

        ! Number of frequencies assigned
        allocate(MPID%nf(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%nf)
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if0)
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if1)


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

        ! Only the master
        if (pid.eq.0) then

          ! Now optimize how they are distributed
          call weighted_split(nproc-1,maxcoun,IW_freq, &
                              MPID%if0(1:nproc-1), &
                              MPID%if1(1:nproc-1), &
                              MPID%nf(1:nproc-1))

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

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = maxval(MPID%nf(1:nproc-1))

      !
      ! If there is only one slave
      !
      else if (nproc.eq.2) then

        ! We are doing MPI
        MPID%mpi = .True.

        ! The number of processors
        nnd = 1
        MPID%nnd = 1

        !
        ! Allocate variables
        !

        ! Number of frequencies assigned
        allocate(MPID%nf(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%nf)
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if0)
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if1)

        ! Determine number of frequencies per CPU
        MPID%nf(0) = nfreq
        MPID%nf(1) = nfreq

        ! Maximum number of frequencies that a processor have
        MPID%nxfreq = nfreq

        ! Set the limits of each processor
        MPID%if0(0) = 1
        MPID%if0(1) = 1
        MPID%if1(0) = nfreq
        MPID%if1(1) = nfreq

      !
      ! If we are in serial
      !
      else

        ! There is no MPI in any form
        MPID%mpi = .False.

        ! Number of frequencies assigned
        allocate(MPID%nf(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%nf)
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if0)
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if1)

        ! The only processor sees everything
        MPID%nf(0) = nfreq
        MPID%if0(0) = 1
        MPID%if1(0) = nfreq

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

        ! Write number of processes
        write(200,'("Running with",1x,i4,1x,"processes")') nproc

        ! For each slave
        do iproc=1,nproc-1

          ! Write split info
          write(200,'("Process",1x,i4,1x,"received",1x,i6,1x,'// &
                    '"frequencies, with a total weight of",'// &
                    '1x,i7,1x," with maximum of",1x,i7)')&
                iproc,MPID%nf(iproc), &
                sum(IW_freq(MPID%if0(iproc):MPID%if1(iproc))), &
                maxval(IW_freq(MPID%if0(iproc):MPID%if1(iproc)))

        end do ! Slaves

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

      ! Check if everything is fine
      call control

      return

      end subroutine setmpi

!#####################################################################
!#####################################################################
!#####################################################################

      !> Split tasks in the radiation transfer problem in the CLE
      !! synthesis\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!     Input(Input_class): Structure with configuration data\n
      !! IW_freq_ou(integer(:)): Expected work load in the output
      !!                         frequency axis\n
      !! IW_freq_in(integer(:)): Expected work load in the input
      !!                         frequency axis
      subroutine setmpi_CLE(MPID,Input,IW_freq_ou,IW_freq_in)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(in):: Input
      integer, dimension(:), intent(inout), target:: IW_freq_in
      integer, dimension(:), intent(inout), target:: IW_freq_ou

      ! Local

      integer:: ios,iaux,iproc,maxcoun,Mwg
      integer:: ifbl,iaxis,mfreq,iran,jfreq,ifreq
      integer, dimension(0:nproc-1):: nf,if0,if1
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

        !
        ! Allocate variables
        !

        ! Number of frequencies assigned
        allocate(MPID%nf(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%nf)
        allocate(MPID%inf(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%inf)
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if0)
        allocate(MPID%iif0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%iif0)
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if1)
        allocate(MPID%iif1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%iif1)


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

            ! Issue error
            umsg = 'Too many processors for this frequency grid'
            call aborted
            return

          end if ! Too many frequencies

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

          ! Only the master
          if (pid.eq.0) then

            ! First CPU
            if0(0) = 1
            if1(0) = nf(0)

            ! For each group
            do ifbl=1,nproc-1
              if0(ifbl) = if1(ifbl-1) + 1
              if1(ifbl) = if0(ifbl) + nf(ifbl) - 1
            end do

            ! Now optimize how they are distributed
            call weighted_split(nproc,maxcoun, &
                                IW_freq(if0(0):if1(nproc-1)), &
                                if0,if1,nf)
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
      ! If we are in serial
      !
      else

        ! There is no MPI in any form
        MPID%mpi = .False.

        ! Number of frequencies assigned
        allocate(MPID%nf(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%nf)
        allocate(MPID%inf(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%inf)
        ! First index of the frequencies assigned
        allocate(MPID%if0(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if0)
        allocate(MPID%iif0(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%iif0)
        ! Last index of the frequencies assigned
        allocate(MPID%if1(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%if1)
        allocate(MPID%iif1(0:0))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%iif1)

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

          MPID%nf(0) = Input%lim_stk%nn
          MPID%if0(0) = 1
          MPID%if1(0) = Input%lim_stk%nn

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

        ! Write number of processes
        write(200,'("Running with",1x,i4,1x,"processes")') nproc

        ! For each slave
        do iproc=1,nproc-1

          ! Write MPI info
          write(200,'("Process",1x,i4,1x,"received",1x,i6,1x,'// &
                    '"frequencies, with a total weight of",'// &
                    '1x,i7,1x," with maximum of",1x,i7)')&
                iproc,MPID%inf(iproc), &
                sum(IW_freq_in(MPID%iif0(iproc):MPID%iif1(iproc))), &
                maxval(IW_freq_in(MPID%iif0(iproc):MPID%iif1(iproc)))

        end do ! Slaves

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

      !> Read weigths and times from a previous run to try readjust
      !! the expected weights\n
      !!      lnproc(integer): Number of processes to split in\n
      !!     maxcoun(integer): Maximum number of tries for the
      !!                       split\n
      !!  IW_freq(integer(:)): Expected work load array
      !!      if0(integer(:)): Initial index for ramge\n
      !!      if0(integer(:)): Final index for range\n
      !!       nf(integer(:)): Number of frequencies for each range
      subroutine weighted_split(lnproc,maxcoun,IW_freq,if0,if1,nf)

      ! I/O

      integer, intent(in):: lnproc,maxcoun
      integer, dimension(:), intent(in):: IW_freq
      integer, dimension(:), intent(inout):: if0,if1,nf

      ! Local

      logical:: change

      integer:: coun,ifbl,ifblb,ii,MTW,LW,TW,Mwg1,TWb,Mwg2,mMwg
      integer:: TTW,TTWb,TWM


      ! Set logical variable
      change = .True.
      coun = 0

      ! While there is change, repeat
      do while (change)

        ! Reset flag
        change = .False.

        ! Get maximum summed weight
        MTW = sum(IW_freq(if0(1):if1(1)))
        do ifbl=2,lnproc
          LW = sum(IW_freq(if0(ifbl):if1(ifbl)))
          if (LW.gt.MTW) MTW = LW
        end do

        !
        ! First, check if there are very big differences and try
        ! big shifts
        !

        ! For each slave
        do ifbl=1,lnproc-1

          ! Get local maximum and sum
          TW = sum(IW_freq(if0(ifbl):if1(ifbl)))
          Mwg1 = maxval(IW_freq(if0(ifbl):if1(ifbl)))

          ! For each other processor
          do ifblb=ifbl+1,lnproc

            ! Get local maximum and sum
            TWb = sum(IW_freq(if0(ifblb):if1(ifblb)))
            Mwg2 = maxval(IW_freq(if0(ifblb): &
                                  if1(ifblb)))


            ! Minimum of the maximum weights
            mMwg = min(Mwg1,Mwg2)

            ! If the difference between the sums if larger than the
            ! minimum of local maxima
            if (abs(TW - TWb).gt.mMwg) then

              ! If the first has larger sum
              if (TW.gt.TWb) then

                ! If the first has more than one node
                if (nf(ifbl).gt.1) then

                  ! Save current sum
                  TTWb = TWb
                  TWM = TW

                  ! For all ranges, get maximum sum
                  TTW = TW
                  do ii=ifbl,ifblb
                    LW = sum(IW_freq(if0(ii):if1(ii)))
                    if (LW.gt.TTW) TTW = LW
                  end do

                  ! Get what would be the sum by removing one node
                  ! from the first
                  LW = sum(IW_freq(if0(ifbl): &
                                   if1(ifbl)-1))

                  ! If the sum would be larger than the second sum,
                  ! assign this to the second
                  if (LW.gt.TWb) TWb = LW

                  ! If the sum would be larger than the maximum sum
                  ! in between, assign this
                  if (LW.gt.TTWb) TTWb = LW

                  ! For all intermediate processes, get maximum sum
                  ! after shifting everything to the left
                  do ii=ifbl+1,ifblb-1
                    LW = sum(IW_freq(if0(ii)-1: &
                                     if1(ii)-1))
                    if (LW.gt.TTWb) TTWb = LW
                  end do

                  ! Get sum of second after adding one to the left
                  LW = sum(IW_freq(if0(ifblb)-1: &
                                   if1(ifblb)))

                  ! If sum in the last is larger than the original,
                  ! assign this
                  if (LW.gt.TWb) TWb = LW

                  ! Is sum in the last is larger than intermediates,
                  ! assing this
                  if (LW.gt.TTWb) TTWb = LW

                  ! If the weight of the second range is smaller than
                  ! the initial first range and the weight of any
                  ! intermediate range or the second range is smaller
                  ! or equal than the previous maximum size among the
                  ! intermediate ranges plus the change in weight
                  if (TWb.lt.TWM.and. &
                      TTWb.le.TTW+abs(TWb-TWM)) then

                    ! Remove one point from the first range
                    if1(ifbl) = if1(ifbl)-1
                    nf(ifbl) = nf(ifbl)-1

                    ! Shift all intermediate ranges
                    do ii=ifbl+1,ifblb-1
                      if0(ii) = if0(ii)-1
                      if1(ii) = if1(ii)-1
                    end do

                    ! Add one point to the second range
                    if0(ifblb) = if0(ifblb)-1
                    nf(ifblb) = nf(ifblb)+1

                    ! Update sum and maximum
                    TW = sum(IW_freq(if0(ifbl): &
                                     if1(ifbl)))
                    Mwg1 = maxval(IW_freq(if0(ifbl): &
                                          if1(ifbl)))

                    ! Flag there was a change
                    change = .True.

                  end if ! It is worth to interchange nodes
                end if ! First range can give node

              ! If the second has larger sum
              else

                ! If the second has more than one node
                if (nf(ifblb).gt.1) then

                  ! Save current sum
                  TTWb = TWb
                  TWM = TWb

                  ! For all ranges, get maximum sum
                  TTW = TW
                  do ii=ifbl,ifblb
                    LW = sum(IW_freq(if0(ii):if1(ii)))
                    if (LW.gt.TTW) TTW = LW
                  end do

                  ! Get sum of first after adding one to the right
                  TWb = sum(IW_freq(if0(ifbl): &
                                          if1(ifbl)+1))

                  ! For all intermediate processes, get maximum sum
                  ! after shifting everything to the right
                  do ii=ifbl+1,ifblb-1
                    LW = sum(IW_freq(if0(ii)+1: &
                                     if1(ii)+1))
                    if (LW.gt.TWb) TWb = LW
                  end do

                  ! Get what would be the sum by removing one node
                  ! from the second
                  LW = sum(IW_freq(if0(ifblb)+1: &
                                   if1(ifblb)))

                  ! If the sum would be larger than the second sum,
                  ! assign this to the second
                  if (LW.gt.TWb) TWb = LW

                  ! If the weight of the second range is smaller than
                  ! the initial first range and the weight of any
                  ! intermediate range or the second range is smaller
                  ! or equal than the previous maximum size among the
                  ! intermediate ranges plus the change in weight
                  if (TWb.lt.TWM.and. &
                      TTWb.le.TTW+abs(TWb-TWM)) then

                    ! Add one point to the first range
                    if1(ifbl) = if1(ifbl)+1
                    nf(ifbl) = nf(ifbl)+1

                    ! Shift all intermediate ranges
                    do ii=ifbl+1,ifblb-1
                      if0(ii) = if0(ii)+1
                      if1(ii) = if1(ii)+1
                    end do

                    ! Remove one point from the second range
                    if0(ifblb) = if0(ifblb)+1
                    nf(ifblb) = nf(ifblb)-1
                    TW = sum(IW_freq(if0(ifbl):if1(ifbl)))
                    Mwg1 = maxval(IW_freq(if0(ifbl): &
                                          if1(ifbl)))

                    ! Flag there was a change
                    change = .True.

                  end if ! It is worth to interchange nodes
                end if ! Second range can give node
              end if ! What range has larger weight
            end if ! It is worth to interchange nodes

          end do ! Each other processor
        end do ! Each processor

        ! Cycle of changes
        coun = coun + 1

        ! If we have changed more than the limit
        if (coun.gt.maxcoun) then

          ! Issue warning
          change = .False.
          umsg = ' # Limit of iterations in task distributer'
          call verbose

        end if ! Changed more than the limit

      end do ! while changing

      end subroutine weighted_split

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the size of the buffers in MPI messages\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!  GeomI(Geometry_class): Structure with geometric data for the
      !!                         intensity problem\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!           lio(logical): If solving for intensity\n
      !!            lp(logical): If solving for polarization\n
      !!            lJ(logical): If solving for continuum\n
      !!          lGen(logical): If correcting for multi-term\n
      !!     ALI_photo(logical): If doing ALI in photoionization\n
      !!         reval(logical): In here to recalculate the sizes
      subroutine setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,lJ,lGen, &
                              ALI_photo,reval)

      ! I/O

      type(MPI_class), intent(inout):: MPID
      type(Geometry_class), intent(in):: Geom,GeomI
      type(Frequency_class), intent(in):: Frec
      logical, intent(in):: lio,lp,lJ,lGen,reval,ALI_photo

      ! Local

      integer:: iproc

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

        ! For each CPU
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

        end do ! For each CPU

        ! Correct photo
        if (.not.ALI_photo) then
          MPID%sizei0(iproc) = 0
          MPID%sizei10(iproc) = 0
        end if

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
          if (.not.ALI_photo) b4 = 0

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
            if (.not.ALI_photo) MPID%sizei0 = 0

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

        ! bcast
        call MPI_BCAST(Frec%Mntfreq(0), nproc, MPI_INTEGER, &
                       0, MPI_COMM_RT, ierr)
        call MPI_BCAST(Frec%Mntfreqi(0), nproc, MPI_INTEGER, &
                       0, MPI_COMM_RT, ierr)
        call MPI_BCAST(Frec%Mnpfreq(0), nproc, MPI_INTEGER, &
                       0, MPI_COMM_RT, ierr)

        ! Maximum size of profile
        MPID%nxtfreq = maxval(Frec%Mntfreq(1:nproc-1))
        MPID%nxtfreqi = maxval(Frec%Mntfreqi(1:nproc-1))
        MPID%nxpfreq = maxval(Frec%Mnpfreq(1:nproc-1))

        ! Allocate size arrays

        ! Size of RT information sent between domains
        allocate(MPID%size1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size1)
        ! Size of profile information sent between domains
        allocate(MPID%size2(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size2)
        ! Size of location indexes vector sent to master
        allocate(MPID%size3(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size3)
        ! Size of Stokes chunk sent to master
        allocate(MPID%size4(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size4)
        ! Size of profile chunk sent to master
        allocate(MPID%size5(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size5)
        ! Size of JKQ to broadcast
        allocate(MPID%size6(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size6)
        ! Size of JKQ frequency dependent to broadcast
        allocate(MPID%size7(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size7)
        ! Size of Stokes to broadcast
        allocate(MPID%size8(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size8)
        ! Size of stokes in emergence sent to master
        allocate(MPID%size10(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%size10)
        ! Size of RT information sent between domains in intensity
        allocate(MPID%sizei1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei1)
        ! Size of Lambda_line to broadcast
        allocate(MPID%sizei2(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei2)
        ! Size of J00P to broadcast
        allocate(MPID%sizei3(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei3)
        ! Size of Stokes chunk sent to master in intensity
        allocate(MPID%sizei4(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei4)
        ! Size of Stokes chunk sent to master in intensity
        allocate(MPID%sizei4b(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei4b)
        ! Size of profile chunk sent to master in intensity
        allocate(MPID%sizei5(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei5)
        ! Size of J00 to broadcast
        allocate(MPID%sizei6(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei6)
        ! Size of J00 frequency dependent to broadcast
        allocate(MPID%sizei7(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei7)
        ! Size of Stokes to broadcast in intensity
        allocate(MPID%sizei8(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei8)
        ! Size of Lambda_line sent to master
        allocate(MPID%sizei9(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei9)
        ! Size of Lamda_photo sent to master
        allocate(MPID%sizei0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei0)
        ! Size of Lambda_photo to broadcast
        allocate(MPID%sizei10(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei10)
        ! Size of profiles between slaves
        allocate(MPID%sizei11(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei11)
        ! Size of b-b r between slaves
        allocate(MPID%sizei12(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei12)
        ! Size of b-f r between slaves
        allocate(MPID%sizei13(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei13)
        ! Size of coontribution function intensity
        allocate(MPID%sizei14(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(MPID%sizei14)

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

        end do ! Processes

        ! Correct photo
        if (.not.ALI_photo) then
          MPID%sizei0(iproc) = 0
          MPID%sizei10(iproc) = 0
        end if


        !
        ! Check sizes of specific buffers than can go beyond limits
        ! for the continuum problem
        !

        ! If doing intensity and continuum
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
              print*,gpid,b1,maxbuffer
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
          if (.not.ALI_photo) b4 = 0

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
              print*,gpid,b1,b2,b3,b4,b1+b2+b3+b4,maxbuffer
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
            if (.not.ALI_photo) MPID%sizei0 = 0

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

          ! If PRD and angle-dependent
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

            ! Get new sizes
            do iproc=0,nproc-1
              MPID%size4(iproc) = 4*MPID%nf(iproc)*nz
              MPID%size5(iproc) = 2*Frec%Mntfreq(iproc)*nz
            end do

          end if ! Big buffer
        end if ! Polarization
      end if ! Recalculating buffer sizes or first calculation

      ! Control
      call control

      return

      end subroutine setmpi_sizes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set all request variables to MPI_REQUEST_NULL\n
      !!  MPID(MPI_class): Structure with MPI data
      subroutine reset_mpirequest(MPID)

      ! I/O

      type(MPI_class), intent(inout):: MPID


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

      end subroutine reset_mpirequest

!#####################################################################
!#####################################################################
!#####################################################################

      !! Write the amomunt of RAM accounted for in an ASCII file\n
      !!  folder(character(:)): Path to the output folder\n
      !!          pol(integer): Indicates if doing polarization\n
      !!       formal(integer): Indicates if doing formal solution
      subroutine RAMreport(folder,pol,formal)

      ! I/O

      character(len=500), intent(in):: folder
      integer, intent(in):: pol,formal

      ! Local

      character(len=11):: fname

      integer:: iproc,ios

      double precision, dimension(:), allocatable:: lRAM,lPRAM,lVRAM
      double precision, dimension(:), allocatable:: lWRAM,lRRAM,lSRAM
      double precision, dimension(:), allocatable:: lBRAM,lMRAM,lTRAM
      double precision, dimension(:), allocatable:: lERAM,lORAM,lFRAM

      ! Allocate processors
      allocate(lRAM(nproc))
      allocate(lPRAM(nproc))
      allocate(lVRAM(nproc))
      allocate(lWRAM(nproc))
      allocate(lRRAM(nproc))
      allocate(lSRAM(nproc))
      allocate(lBRAM(nproc))
      allocate(lMRAM(nproc))
      allocate(lTRAM(nproc))
      allocate(lERAM(nproc))
      allocate(lORAM(nproc))
      allocate(lFRAM(nproc))


      ! If MPI
      if (nproc.gt.1) then

        ! Tell everything to master
        call MPI_GATHER(PRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lPRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(VRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lVRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(WRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lWRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(RRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lRRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(SRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lSRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(BRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lBRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(MRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lMRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(TRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lTRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(ERAMc, 1, MPI_DOUBLE_PRECISION, &
                        lERAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(ORAMc, 1, MPI_DOUBLE_PRECISION, &
                        lORAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
        call MPI_GATHER(FRAMc, 1, MPI_DOUBLE_PRECISION, &
                        lFRAM, 1, MPI_DOUBLE_PRECISION, 0, &
                        MPI_COMM_RT, ierr)
      ! No MPI
      else

        ! Copy from commons
        lPRAM(1) = PRAMc
        lVRAM(1) = VRAMc
        lWRAM(1) = WRAMc
        lRRAM(1) = RRAMc
        lSRAM(1) = SRAMc
        lBRAM(1) = BRAMc
        lMRAM(1) = MRAMc
        lTRAM(1) = TRAMc
        lERAM(1) = ERAMc
        lORAM(1) = ORAMc
        lFRAM(1) = FRAMc

      end if ! MPI

      ! Get sum
      lRAM = lPRAM + lVRAM + lWRAM + LRRAM + lSRAM + &
             LBRAM + lMRAM + lTRAM + lERAM + lORAM + &
             lFRAM


      ! Only master writes
      if (pid.eq.0) then

        !
        ! Choose filename
        !

        ! Polarization
        if (pol.eq.1) then

          ! Formal solution
          if (formal.eq.1) then
            fname = 'RAM-P-F.dat'

          ! Self-consistent solution
          else
            fname = 'RAM-P-E.dat'

          endif ! Type of solution

        ! Intensity
        else

          ! Formal solution
          if (formal.eq.1) then
            fname = 'RAM-I-F.dat'

          ! Self-consistent solution
          else
            fname = 'RAM-I-E.dat'

          endif ! Type of solution
        end if ! Polarization

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
                         ' |   Solution '// &
                         ' | Memoizati. '// &
                         ' | Voigt prf. '// &
                         ' | Redistrib. '// &
                         ' | 1st o. PRD '// &
                         ' | Rad. Tran. '// &
                         ' |      Misc. '// &
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
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--'

        ! For each CPU
        do iproc=1,nproc

          ! Write RAM count
          write(200,'(2x,i4,11(" | ",2x,f9.3)," |")') &
            iproc-1, &
            lBRAM(iproc),lPRAM(iproc),lRRAM(iproc), &
            lSRAM(iproc),lERAM(iproc),lVRAM(iproc), &
            lWRAM(iproc),lORAM(iproc),lTRAM(iproc), &
            lMRAM(iproc),lRAM(iproc)

        end do ! CPUs

        ! Line
        write(200,'(A)') '------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
                         '--------------'// &
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

      end module setmpi_mod
