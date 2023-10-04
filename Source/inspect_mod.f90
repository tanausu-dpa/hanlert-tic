      !> Manages input spectrum for CLE
      module inspect_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     11/14/2022
!  Last version:
!     11/24/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/24/2022:    V3.0.0 - First version (TdPA)
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
!    Manages the input radiation for CLE
!
!    rinspect:
!      Read and store input spectra
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : isfrac , bbdis, bfdis , c
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the input spectrum file.\n
      !!    Input(Input_class): Structure with settings data\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!  Geom(Geometry_class): Structure with quadrature data\n
      !!    spect(spect_class): Structure with the input spectra
      !!                        data\n
      !!      omega(double(:)): Frequency axis
      subroutine rinspect(Input,Atom,Geom,spect,omega)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Input_class):: Input
      type(Geometry_class), intent(in):: Geom
      type(spect_class), intent(out):: spect
      double precision, dimension(:), intent(in):: omega

      ! Local

      logical:: found
      logical, dimension(:), allocatable:: bspecin,fspecin

      integer:: ios,ia,i0,i1,di,ji,jj,j0,j1,dj,ki,kk,k0,k1,dk
      integer:: jtran,ifreq,jfreq,iproc
      integer, dimension(5):: sizes
      integer, dimension(:), allocatable:: add

      double precision:: w0,w1,lim
      double precision, dimension(:), allocatable:: v1,v2,v3
      double precision, dimension(:,:,:,:), allocatable:: v1234


      ! Slaves need to initialize some quantities initialize
      if (gpid.gt.0) then

        ! For each atom
        do ia=1,nA

          ! If b-b, initialize
          if (Atom(ia)%ntran.gt.0) then
            allocate(Atom(ia)%bbspecin(Atom(ia)%ntran))
            Atom(ia)%bbspecin = .False.
          end if

          ! If b-f, initialize
          if (Atom(ia)%nphot.gt.0) then
            allocate(Atom(ia)%bfspecin(Atom(ia)%nphot))
            Atom(ia)%bfspecin = .False.
          end if

        end do ! For each atom

      end if ! Slaves only

      ! Check if input file
      if (trim(Input%spect_input).eq.'NONE') then
        Input%lspect_input = .False.
        spect%valid = .False.
        return
      end if

      ! There is a file
      Input%lspect_input = .True.


      !
      ! Master
      !
      if (gpid.eq.0) then

        ! Open file with spectral data and read it
        open(200,file=trim(Input%spect_input),status='old', &
             iostat=ios,err=1000,access='stream', &
             form='unformatted')

        ! Read file dimensions
        read(200,err=1100) sizes(1:4)
        sizes(5) = sizes(1)*sizes(2)*sizes(3)*sizes(4)

        ! Allocate buffers
        allocate(v1(sizes(1)),v2(sizes(2)),v3(sizes(3)))
        allocate(v1234(sizes(4),sizes(3),sizes(2),sizes(1)))

        ! Read buffers
        read(200,err=1100) v1
        read(200,err=1100) v2
        read(200,err=1100) v3
        read(200,err=1100) v1234

        close(200)

        ! control reading
        call control

        ! Broadcast
        call MPI_BCAST(sizes(1), 5, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v1(1), sizes(1), MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v2(1), sizes(2), MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v3(1), sizes(3), MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v1234(1,1,1,1), sizes(5), &
                       MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)

      !
      ! Slave
      !
      else

        ! control reading
        call control

        ! Broadcast
        call MPI_BCAST(sizes(1), 5, MPI_INTEGER, 0, &
                       MPI_COMM_WORLD, ios)

        ! Allocate buffers
        allocate(v1(sizes(1)),v2(sizes(2)),v3(sizes(3)))
        allocate(v1234(sizes(4),sizes(3),sizes(2),sizes(1)))
        allocate(add(sizes(3)))

        ! Receive by bcast
        call MPI_BCAST(v1(1), sizes(1), MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v2(1), sizes(2), MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v3(1), sizes(3), MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)
        call MPI_BCAST(v1234(1,1,1,1), sizes(5), &
                       MPI_DOUBLE_PRECISION, 0, &
                       MPI_COMM_WORLD, ios)

        ! Initialize included frequencies
        add = 0

      end if ! Master or slave

      ! Sanity check the stokes size
      if (sizes(4).ne.1.and.sizes(4).ne.4) then

        ! Abort
        umsg = 'The Stokes size in the spectrum file must be '// &
               '1 or 4'
        urou = 'rinspect'
        call gaborted

      end if

      ! Check if polarized
      if (sizes(4).eq.4) then
        spect%pol = .True.
        spect%nstk = 3
      else
        spect%pol = .False.
        spect%nstk = 0
      end if

      ! Check if axial
      spect%axial = sizes(2).eq.1

      !
      ! Master can leave now
      !
      if (gpid.eq.0) return

      ! Change units
      ! m -> 10^5 cm^-1
      v3 = 1d-7/v3

      ! Choose order
      if(v3(2).gt.v3(1))then
        i0 = 1
        i1 = sizes(3)
        di = 1
      else
        i0 = sizes(3)
        i1 = 1
        di = -1
      end if

      ! For each atom
      do ia=1,nA

        ! Allocate and initialize limits
        if (Atom(ia)%ntran.gt.0) then
          allocate(Atom(ia)%sbif0(Atom(ia)%ntran))
          allocate(Atom(ia)%sbif1(Atom(ia)%ntran))
          Atom(ia)%sbif0 = nfreq + 1
          Atom(ia)%sbif1 = 0
        end if
        if (Atom(ia)%nphot.gt.0) then
          allocate(Atom(ia)%sfif0(Atom(ia)%nphot))
          allocate(Atom(ia)%sfif1(Atom(ia)%nphot))
          Atom(ia)%sfif0 = nfreq + 1
          Atom(ia)%sfif1 = 0
        end if

        ! For each b-b transition
        do jtran=1,Atom(ia)%ntran

          ! If this CPU does not have frequencies in this line, skip
          if (Atom(ia)%fflag(jtran)%absent) cycle

          ! Get boundaries
          w0 = omega(Atom(ia)%if0(jtran))
          w1 = omega(Atom(ia)%if1(jtran))

          ! Look for frequencies in the input axis
          do ifreq=i0,i1,di

            ! Ordered frequency
            jfreq = abs(ifreq - i0) + 1

            ! Check if frequency in range
            if (w0.le.v3(ifreq).and.w1.ge.v3(ifreq)) then

              ! If this line not initialized
              if (Atom(ia)%sbif0(jtran).gt.nfreq) then

                ! If not first point
                if (jfreq.gt.1) then

                  ! Left limit, fraction of the line range size
                  lim = w0 - (w1 - w0)*isfrac

                  ! If previous within this limit
                  if (v3(ifreq-di).ge.lim) then
                    Atom(ia)%sbif0(jtran) = jfreq - 1
                    Atom(ia)%sbif1(jtran) = jfreq - 1
                    add(ifreq-di) = 1
                  end if ! Previous point within limit
                end if ! Not first point in input array
              end if ! Line not initialized

              ! Update limits
              if (Atom(ia)%sbif0(jtran).gt.jfreq) &
                Atom(ia)%sbif0(jtran) = jfreq
              if (Atom(ia)%sbif1(jtran).lt.jfreq) &
                Atom(ia)%sbif1(jtran) = jfreq
              add(ifreq) = 1

            ! If above
            else if (v3(ifreq).gt.w1) then

              ! Right limit, fraction of the line range size
              lim = w1 + (w1 - w0)*isfrac

              ! If previous within this limit
              if (v3(ifreq).le.lim) then

                ! If this line not initialized
                if (Atom(ia)%sbif0(jtran).gt.nfreq) &
                  Atom(ia)%sbif0(jtran) = jfreq
                Atom(ia)%sbif1(jtran) = jfreq
                add(ifreq) = 1
              end if

              exit

            end if ! Above line range

          end do ! Input frequencies
        end do ! b-b transitions

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! Get boundaries
          w0 = Atom(ia)%phot(jtran)%infreq(1)
          w1 = Atom(ia)%phot(jtran)%infreq(Atom(ia)%phot(jtran)%nfreq)

          ! Look for frequencies in the input axis
          do ifreq=i0,i1,di

            ! Ordered frequency
            jfreq = abs(ifreq - i0) + 1

            ! Check if frequency in range
            if (w0.le.v3(ifreq).and.w1.ge.v3(ifreq)) then

              ! If this line not initialized
              if (Atom(ia)%sfif0(jtran).gt.nfreq) then

                ! If not first point
                if (jfreq.gt.1) then

                  ! Left limit, fraction of the line range size
                  lim = w0 - (w1 - w0)*isfrac

                  ! If previous within this limit
                  if (v3(ifreq-di).ge.lim) then
                    Atom(ia)%sfif0(jtran) = jfreq - 1
                    Atom(ia)%sfif1(jtran) = jfreq - 1
                    add(ifreq-di) = 1
                  end if ! Previous point within limit
                end if ! Not first point in input array
              end if ! Line not initialized

              ! Update limits
              if (Atom(ia)%sfif0(jtran).gt.jfreq) &
                Atom(ia)%sfif0(jtran) = jfreq
              if (Atom(ia)%sfif1(jtran).lt.jfreq) &
                Atom(ia)%sfif1(jtran) = jfreq
              add(ifreq) = 1

            end if ! In frequency range

            ! If above
            if (v3(ifreq).gt.w1) then

              ! Right limit, fraction of the line range size
              lim = w1 + (w1 - w0)*isfrac

              ! If previous within this limit
              if (v3(ifreq).le.lim) then
                ! If this line not initialized
                if (Atom(ia)%sfif0(jtran).gt.nfreq) &
                  Atom(ia)%sfif0(jtran) = jfreq
                Atom(ia)%sfif1(jtran) = jfreq
                add(ifreq) = 1
              end if

              exit

            end if ! Above line range

          end do ! Input frequencies
        end do ! b-f transitions
      end do ! Atoms

      ! Now only keep the flagged frequencies
      spect%nmu = sizes(1)
      spect%nphi = sizes(2)
      spect%nfreq = sum(add)

      ! Not valid dimensions
      if (spect%nmu.lt.1.or.spect%nphi.lt.1.or. &
          spect%nfreq.lt.1) then
        spect%valid = .False.
        call control
        return
      end if

      ! Check direction of mu
      if (v1(2).gt.v1(1)) then
        j0 = 1
        j1 = sizes(1)
        dj = 1
      else
        j0 = sizes(1)
        j1 = 1
        dj = -1
      end if

      ! Check direction of phi
      if (sizes(2).gt.1) then
        if (v2(2).gt.v2(1)) then
          k0 = 1
          k1 = sizes(2)
          dk = 1
        else
          k0 = sizes(2)
          k1 = 1
          dk = -1
        end if
      else
        k0 = 1
        k1 = 1
        dk = 1
      end if

      ! Allocate and copy
      allocate(spect%mu(spect%nmu))
      allocate(spect%phi(spect%nphi))
      allocate(spect%omega(spect%nfreq))
     !allocate(spect%stokes(0:spect%nstk,spect%nfreq, &
     !                      spect%nphi,spect%nmu))
      allocate(spect%stokes(spect%nmu,spect%nphi,spect%nfreq, &
                            0:spect%nstk))
      spect%mu = v1(j0:j1:dj)
      spect%phi = v2(k0:k1:dk)

      ! Initialize index
      jj = 0

      ! Mu directions
      do ji=j0,j1,dj

        ! Advance index
        jj = jj + 1

        ! Initialize index
        kk = 0

        ! Axial directions
        do ki=k0,k1,dk

          ! Advance index
          kk = kk + 1

          ! Initialize freq. index
          jfreq = 0

          ! Go by the whole input axis
          do ifreq=i0,i1,di

            ! If to keep, add
            if (add(ifreq).gt.0) then

              ! Advance index
              jfreq = jfreq + 1

              ! Store frequency and intensity (expected in SI)
              spect%omega(jfreq) = v3(ifreq)
             !spect%stokes(:,jfreq,kk,jj) = v1234(:,ifreq,ki,ji)* &
              spect%stokes(jj,kk,jfreq,:) = v1234(:,ifreq,ki,ji)* &
                                            c*1d13

            end if ! Frequency to keep

          end do ! Check the whole frequency axis
        end do ! Axial direction
      end do ! Mu directions


      !
      ! Check which lines can be interpolated from
      ! the spectra
      !

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do jtran=1,Atom(ia)%ntran

          ! Get line limits
          w0 = omega(Atom(ia)%if0(jtran))
          w1 = omega(Atom(ia)%if1(jtran))

          ! Initialize
          found = .False.
          jfreq = spect%nfreq

          ! Check red limit
          do ifreq=1,spect%nfreq
            if (abs(spect%omega(ifreq)-w0).lt.bbdis) then
              found = .True.
              jfreq = ifreq
              exit
            end if
            if (spect%omega(ifreq).gt.w0) exit
          end do

          ! If found red, check blue
          if (found) then
            found = .False.
            do ifreq=jfreq,spect%nfreq
              if (abs(spect%omega(ifreq)-w1).lt.bbdis) then
                found = .True.
                exit
              end if
              if (spect%omega(ifreq).gt.w1) exit
            end do
          end if

          Atom(ia)%bbspecin(jtran) = found

        end do ! b-b transitions

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! Get photo limits
          w0 = Atom(ia)%phot(jtran)%infreq(1)
          w1 = Atom(ia)%phot(jtran)%infreq(Atom(ia)%phot(jtran)%nfreq)

          ! Initialize
          found = .False.
          jfreq = spect%nfreq

          ! Check red limit
          do ifreq=1,spect%nfreq
            if (abs(spect%omega(ifreq)-w0).lt.bfdis) then
              found = .True.
              jfreq = ifreq
              exit
            end if
            if (spect%omega(ifreq).gt.w0) exit
          end do

          ! If found red, check blue
          if (found) then
            found = .False.
            do ifreq=jfreq,spect%nfreq
              if (abs(spect%omega(ifreq)-w1).lt.bfdis) then
                found = .True.
                exit
              end if
              if (spect%omega(ifreq).gt.w1) exit
            end do
          end if

          Atom(ia)%bfspecin(jtran) = found

        end do ! b-f transitions

        ! If splitting frequency, we need to communicate
        if (nproc.gt.1) then

          ! The leader
          if (pid.eq.0) then

            ! Allocate buffers
            allocate(bspecin(Atom(ia)%ntran))
            if (Atom(ia)%nphot.gt.0) &
              allocate(fspecin(Atom(ia)%nphot))

            ! For all processes
            do iproc=1,nproc-1

              ! Receive
              call MPI_RECV(bspecin(1), Atom(ia)%ntran, &
                            MPI_LOGICAL, iproc, iproc, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ios)
              if (Atom(ia)%nphot.gt.0) &
                call MPI_RECV(fspecin(1), Atom(ia)%nphot, &
                              MPI_LOGICAL, iproc, 1000000+iproc, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ios)

              ! Update own info
              do jtran=1,Atom(ia)%ntran
                Atom(ia)%bbspecin(jtran) = bspecin(jtran).and. &
                                           Atom(ia)%bbspecin(jtran)
              end do ! Transitions
              do jtran=1,Atom(ia)%nphot
                Atom(ia)%bfspecin(jtran) = fspecin(jtran).and. &
                                           Atom(ia)%bfspecin(jtran)
              end do ! Photoionizations
            end do ! Rest of members

            ! Send to the rest of members
            call MPI_BCAST(Atom(ia)%bbspecin(1), Atom(ia)%ntran, &
                           MPI_LOGICAL,0,MPI_COMM_RT,ios)
            if (Atom(ia)%nphot.gt.0) &
              call MPI_BCAST(Atom(ia)%bfspecin(1), Atom(ia)%nphot, &
                             MPI_LOGICAL,0,MPI_COMM_RT,ios)

            ! Deallocate buffers
            deallocate(bspecin)
            if (allocated(fspecin)) deallocate(fspecin)

          ! Everyone but the leader
          else

            ! Send to leader
            call MPI_SEND(Atom(ia)%bbspecin(1), Atom(ia)%ntran, &
                          MPI_LOGICAL, 0, pid, MPI_COMM_RT, ios)
            if (Atom(ia)%nphot.gt.0) &
              call MPI_SEND(Atom(ia)%bfspecin(1), Atom(ia)%nphot, &
                            MPI_LOGICAL, 0, 1000000+pid, &
                            MPI_COMM_RT, ios)

            ! Receive from leader
            call MPI_BCAST(Atom(ia)%bbspecin(1), Atom(ia)%ntran, &
                          MPI_LOGICAL, 0, MPI_COMM_RT, ios)
            if (Atom(ia)%nphot.gt.0) &
              call MPI_BCAST(Atom(ia)%bfspecin(1), Atom(ia)%nphot, &
                             MPI_LOGICAL, 0, MPI_COMM_RT, ios)

          end if ! Leader
        end if ! Splitting in frequency

      end do ! Atoms

      ! Check that something remains
      spect%valid = .False.
      do ia=1,nA

        spect%valid = spect%valid.OR.ANY(Atom(ia)%bbspecin)

        if (Atom(ia)%nphot.gt.0) &
          spect%valid = spect%valid.OR.ANY(Atom(ia)%bfspecin)

      end do

      ! If there is data
      if (spect%valid) then

        ! If axial
        if (spect%axial) then

          allocate(spect%mustokes(0:spect%nstk,spect%nfreq, &
                                  1,Geom%nTh))

        ! Not axial
        else

          allocate(spect%mustokes(0:spect%nstk,spect%nfreq, &
                                  Geom%nPh,Geom%nTh))

        end if

      ! If nothing, remove data
      else

        deallocate(spect%mu,spect%phi,spect%omega,spect%stokes)

      end if

      ! Control
      call control

      return

1000  umsg = 'Error opening spectrum file'
      urou = 'specin'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1100  umsg = 'Error reading spectrum file'
      urou = 'specin'
      close(100)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine rinspect

!#####################################################################
!#####################################################################
!#####################################################################

      end module inspect_mod
