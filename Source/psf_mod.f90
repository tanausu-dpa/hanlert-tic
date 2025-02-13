      !> Postprocessing out emerging Stokes profiles
      module psf_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     23/02/2023
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Removed references to threads in the
!                             calls to abortedS (TdPA)
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
!  PSF_Stokes
!    Apply spectral convolution to the Stokes profiles
!
!  Profiles_Out
!    Process the emergent Stokes profiles for the inversion
!
!  cut_stk
!    Slice the Stokes parameters to only include the requested
!  grid points
!
!  set_psf_ranges
!    Prepare the ranges to apply the PSF
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use inter_mod
      use parameters_mod, only: sqrt2, fw2sg, TINYO
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Apply spectral convolution to the Stokes profiles\n
      !!  FWHM(FWHM_helper_class): Structure with the data on the
      !!                           spectral PSF\n
      !!               n(integer): Dimension of the wavelength axis\n
      !!            nstk(integer): Last index in Stokes dimension\n
      !!      Stokes(double(:,:)): Stokes parameters\n
      !!        lambda(double(:)): Wavelength axis
      subroutine PSF_Stokes(FWHM,n,nstk,Stokes,lambda)

      ! I/O

      type(FWHM_helper_class), dimension(:), intent(inout):: FWHM
      integer, intent(in):: n,nstk
      double precision, dimension(n), intent(in):: lambda
      double precision, dimension(0:nstk,n), intent(inout):: Stokes

      ! Local

      integer:: i,j,k,iran,il,ir,inl,inr,ifreq0,lfirst,last

      double precision:: Sigma,Sigma2,Sigma20,dt,dt2
      double precision, dimension(0:nstk,n):: Stokes_tmp
      double precision, dimension(:), allocatable:: twave
      double precision, dimension(:,:), allocatable:: aStokes


      ! Copy Stokes parameters
      Stokes_tmp = Stokes

      ! For each FWHM range
      do iran=1,FWHM(1)%nn

        ! Skip if not to apply
        if (FWHM(iran)%indx(2).lt.1) cycle

        ! If Gaussian
        if (FWHM(iran)%gaussian) then

          ! Skip if not to apply
          if (FWHM(iran)%doub(3).le.0d0) cycle

          ! Get sigma of the gaussian
          Sigma = FWHM(iran)%doub(3)*fw2sg
          Sigma2 = sqrt2*Sigma
          Sigma20 = 20d0*Sigma

          ! Get indexes
          il = FWHM(iran)%indx(1)
          ir = FWHM(iran)%indx(2)
          inl = FWHM(iran)%indx(3)
          inr = FWHM(iran)%indx(4)

          ! Initialize Stokes output in relevant range
          Stokes(:,il:ir) = 0d0

          ! For each output wavelength
          do i=il,ir

            ! If closer than 20 sigmas to the lower limit
            if (lambda(i) - lambda(inl).lt.Sigma20) then

              ! Step
              dt = 0.5d0*(lambda(inl+1) - lambda(inl))

              ! Output
              Stokes(:,i) = 0.5d0*Stokes_tmp(:,inl)* &
                           (erf((lambda(inl)-lambda(i)+dt)/Sigma2) + &
                            1d0)

            ! If closer than 20 sigmas to the upper limit
            else if (lambda(inr) - lambda(i).lt.Sigma20) then

              ! Step
              dt = 0.5d0*(lambda(inr) - lambda(inr-1))

              ! Output
              Stokes(:,i) = 0.5d0*Stokes_tmp(:,inr)* &
                            (1d0 - &
                             erf((lambda(inr)-lambda(i)-dt)/Sigma2))

            end if ! 20 sigmas to any limit

            ! For all wavelengths (except extremes)
            do j=inl+1,inr-1

              ! If closer than 20 sigmas to the output
              if (abs(lambda(j)-lambda(i)).lt.Sigma20) then

                ! Step to the next
                dt = 0.5d0*(lambda(j+1) - lambda(j))
                dt2 = 0.5d0*(lambda(j) - lambda(j-1))

                ! Add contribution
                Stokes(:,i) = Stokes(:,i) + &
                              0.5d0*Stokes_tmp(:,j)* &
                 (erf((lambda(j)-lambda(i) + dt)/Sigma2) - &
                  erf((lambda(j)-lambda(i) - dt2)/Sigma2))

              end if ! Distance between wavelengths

            end do ! Running integral
          end do ! Output frequency

        ! NOT gaussian
        else

          ! Get indexes
          il = FWHM(iran)%indx(1)
          ir = FWHM(iran)%indx(2)
          inl = FWHM(iran)%indx(3)
          inr = FWHM(iran)%indx(4)
          ifreq0 = inl

          ! Allocate true wavelength axis and auxiliar Stokes
          allocate(twave(FWHM(iran)%nfreq))
          allocate(aStokes(0:nstk,FWHM(iran)%nfreq))

          ! Initialize Stokes output
          Stokes(:,il:ir) = 0d0

          ! If need to initialize
          if (FWHM(iran)%toinit) then

            ! Allocate space for auxiliars for interpolation
            allocate(FWHM(iran)%indx1(FWHM(iran)%nfreq,il:ir))
            allocate(FWHM(iran)%indx2(FWHM(iran)%nfreq,il:ir))
            allocate(FWHM(iran)%idx(FWHM(iran)%nfreq,il:ir))
            MRAMc = MRAMc + 1d-6*sizeof(FWHM(iran)%indx1)
            MRAMc = MRAMc + 1d-6*sizeof(FWHM(iran)%indx2)
            MRAMc = MRAMc + 1d-6*sizeof(FWHM(iran)%idx)

            ! Flag initialized
            FWHM(iran)%toinit = .False.

            ! Initialize last first
            lfirst = inl

            ! For each output wavelength
            do i=il,ir

              ! Get new axis
              twave = FWHM(iran)%wave + lambda(i)

              ! Initialize last
              last = lfirst

              ! For each input wavelength
              do j=1,FWHM(iran)%nfreq

                ! If out of range, take boundary
                if (twave(j).le.lambda(1)+TINYO) then

                  ! Set last
                  last = 1
                  if (j.eq.1) lfirst = 1

                  ! Set RAM
                  FWHM(iran)%indx1(j,i) = 1
                  FWHM(iran)%indx2(j,i) = 1
                  FWHM(iran)%idx(j,i) = 0d0

                ! Out of range, take boundary
                else if (twave(j).ge.lambda(inr)) then

                  ! Set last
                  last = inr
                  if (j.eq.1) lfirst = inr

                  ! Set RAM
                  FWHM(iran)%indx1(j,i) = inr
                  FWHM(iran)%indx2(j,i) = inr
                  FWHM(iran)%idx(j,i) = 0d0

                ! Within boundaries
                else

                  ! Search between last found and all
                  ! but boundary
                  do k=last,inr-1

                    ! Exact
                    if (abs(twave(j)-lambda(k)).lt.TINYO) then

                      ! Set last
                      last = k
                      if (j.eq.1) lfirst = k

                      ! Set RAM
                      FWHM(iran)%indx1(j,i) = k
                      FWHM(iran)%indx2(j,i) = k
                      FWHM(iran)%idx(j,i) = 0d0

                      ! Done
                      exit

                    ! Between this and the next
                    else if (twave(j).ge.lambda(k).and. &
                             twave(j).lt.lambda(k+1)) then

                      ! Set last
                      last = k
                      if (j.eq.1) lfirst = k

                      ! Set RAM
                      FWHM(iran)%indx1(j,i) = k
                      FWHM(iran)%indx2(j,i) = k+1
                      FWHM(iran)%idx(j,i) = (twave(j) - lambda(k))/ &
                                            (lambda(k+1) - lambda(k))

                      ! Done
                      exit

                    end if ! Exact or between

                  end do ! Search in frequencies

                end if ! Check boundaries

              end do ! Input wavelengths
            end do ! Output wavelengths

          end if ! Initialize interpolation

          ! For each output wavelength
          do i=il,ir

            ! For each input
            do j=1,FWHM(iran)%nfreq

              ! Initialize auxiliar Stokes
              aStokes(:,j) = Stokes_tmp(:,FWHM(iran)%indx1(j,i))

              ! If interpolating
              if (FWHM(iran)%idx(j,i).gt.0d0) then

                ! Complete
                aStokes(:,j) = aStokes(:,j) + &
                              (Stokes_tmp(:,FWHM(iran)%indx2(j,i)) - &
                               aStokes(:,j))*FWHM(iran)%idx(j,i)
              end if

            end do ! Input frequencies

            ! Compute convolution
            Stokes(0,i) = sum(aStokes(0,:)*FWHM(iran)%kernel)

            ! Pol
            if (nstk.gt.0) then
              Stokes(1,i) = sum(aStokes(1,:)*FWHM(iran)%kernel)
              Stokes(2,i) = sum(aStokes(2,:)*FWHM(iran)%kernel)
              Stokes(3,i) = sum(aStokes(3,:)*FWHM(iran)%kernel)
            end if

          end do ! Output frequencies

          ! Free
          deallocate(twave,aStokes)

        end if ! Type of PSF

      end do ! FWHM ranges

      return

      end subroutine PSF_Stokes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Process the emergent Stokes profiles for the inversion\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!    e_stk(double(:,:,:,:)): Emergent Stokes parameters\n
      !!              lpe(logical): If polarization\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable\n
      !!   FWHM(FWHM_helper_class): Structure with the data on the
      !!                            spectral PSF
      subroutine Profiles_Out(Frec,Sol,e_stk,lpe,buff,FWHM)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      type(FWHM_helper_class), dimension(:), &
                               allocatable, intent(inout):: FWHM
      type(Frequency_class), intent(in):: Frec
      type(Solution_class), intent(inout):: Sol
      logical, intent(in):: lpe
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: e_stk

      ! Local

      logical:: cut,PSF

      integer:: nstk

      double precision, dimension(:), allocatable:: lambda
      double precision, dimension(:,:,:,:), allocatable:: c_stk


      ! Only the Master does this
      if (pid.eq.0) then

        ! Check if we need to cut
        cut = buff%nran.gt.0
        PSF = allocated(FWHM)

        ! Initialize
        Sol%Stokes_out = 0d0

        ! If polarization
        if (lpe) then

          ! Last index is 3 for V
          nstk = 3

        ! Intensity
        else

          ! Last index is 0 for I
          nstk = 0

        end if

        ! Get wavelength axis
        allocate(lambda(nfreq))
        lambda = 1d2/Frec%omega(nfreq:1:-1)

        ! Copy Stokes
        allocate(c_stk(0:nstk,nfreq,1,1))
        c_stk = e_stk(0:nstk,nfreq:1:-1,1:1,1:1)

        ! If FWHM, apply PSF
        if (PSF) call PSF_Stokes(FWHM,nfreq,nstk,c_stk,lambda)

        ! Interpolate into observations axis
        call Intpol_Lin_stk(lambda,c_stk,4,nstk+1, &
                            nfreq,Sol%omega_input, &
                            Sol%Stokes_out,Sol%Num_Wavelength)

        ! Cut Stokes
        if (cut) then

          ! Copy again
          c_stk = e_stk(0:nstk,:,1:1,1:1)

          ! Cut
          call cut_stk(e_stk,c_stk,buff,nstk)

        end if ! Cut Stokes

        ! Free
        deallocate(c_stk,lambda)

      end if ! Master

      ! Share result
      call MPI_BCAST(Sol%Stokes_out(0,1), 4*Sol%Num_Wavelength, &
                     MPI_DOUBLE_PRECISION, 0, MPI_COMM_RT, ierr)

      return

      end subroutine Profiles_Out

!#####################################################################
!#####################################################################
!#####################################################################

      !> Slice the Stokes parameters to only include the requested
      !! grid points\n
      !!    o_stk(double(:,:,:,:)): Cut emergent Stokes parameters\n
      !!    i_stk(double(:,:,:,:)): Full emergent Stokes parameters\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable\n
      !!             nstk(integer): Polarization last index
      subroutine cut_stk(o_stk,i_stk,buff,nstk)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      integer, intent(in):: nstk
      double precision, dimension(:,:,:,:), intent(in):: i_stk
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: o_stk

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Reallocate
      deallocate(o_stk)
      allocate(o_stk(0:nstk,buff%nn,1,1))

      ! Initialize buffer
      ii = 0

      ! For each entry to write
      do iran=1,buff%nran

        ! Range and size
        i0 = buff%indx(1,iran)
        i1 = buff%indx(2,iran)
        nn = i1-i0+1

        ! Fill buffer
        o_stk(0:nstk,ii+1:ii+nn,1,1) = i_stk(1:nstk+1,i0:i1,1,1)

        ! Advance index
        ii = ii + nn

      end do ! Ranges

      end subroutine cut_stk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the ranges to apply the PSF\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!   buff(FWHM_helper_class): Structure with the data on the
      !!                            spectral PSF\n
      !!    files(strarr_class(:)): PSF kernel filenames\n
      !!       omega_in(double(:)): Frequency axis from data\n
      !!          omega(double(:)): Frequency axis for radiation
      !!                            transfer
      subroutine set_psf_ranges(Inf_Stokes,buff,files,omega_in,omega)

      ! I/O

      type(Stokes_class), intent(in):: Inf_Stokes
      type(FWHM_helper_class), dimension(:), &
                               allocatable, intent(inout):: buff
      type(strarr_class), dimension(:), &
                          allocatable, intent(inout):: files
      double precision, dimension(:), intent(in):: omega_in,omega

      ! Local

      logical:: left,right

      integer:: ir,ir1,irl,irr,ileft,iright,ifreq,ios

      double precision:: sig20,wleft,wright,wileft,wiright
      double precision, dimension(:), allocatable:: lambda

      ! If no PSF, skip
      if (.not.allocated(buff)) return

      ! For each range
      do ir=1,buff(1)%nn

        ! Allocate indexes (1 and 2 for in, 3 and 4 for out)
        allocate(buff(ir)%indx(4))
        MRAMc = MRAMc + 1d-6*sizeof(buff(ir)%indx)
        buff(ir)%indx(1) = nfreq + 1
        buff(ir)%indx(3) = buff(ir)%indx(1)
        buff(ir)%indx(2) = 0
        buff(ir)%indx(4) = 0

      end do ! Ranges

      ! Get radiation transfer wavelength axis
      allocate(lambda(nfreq))
      lambda = 1d2/omega(nfreq:1:-1)

      ! For each PSF range
      do ir=1,buff(1)%nn

        ! If gaussian
        if (buff(ir)%gaussian) then

          ! Get 20 sigma of the gaussian
          sig20 = buff(ir)%doub(3)
          sig20 = sqrt2*sig20*20d0

          ! No sigma, skip
          if (sig20.le.0d0) cycle

          ! Wavelength limits
          wleft = buff(ir)%doub(1)
          wright = buff(ir)%doub(2)

          ! Look over wavelength ranges in data
          do ir1=1,Inf_Stokes%Num_Range

            ! Shorter limits
            irl = Inf_Stokes%Range(ir1,1)
            irr = Inf_Stokes%Range(ir1,2)
            wileft = omega_in(irl)
            wiright = omega_in(irr)

            ! Out of bounds comparisons
            left = wright.lt.wileft
            right = wleft.gt.wiright

            ! Skip if out of bounds
            if (left.or.right) cycle

            !
            ! Prepare output
            !

            ! Left index
            ileft = minloc(abs(lambda - wileft),1)

            ! Move one?
            if (lambda(ileft).gt.wileft) then
              if (ileft.gt.irl.and.ileft.gt.1) ileft = ileft - 1
            end if

            ! Right index
            iright = minloc(abs(lambda - wiright),1)

            ! Move one?
            if (lambda(iright).lt.wiright) then
              if (iright.lt.irr.and.iright.lt.nfreq) &
                iright = iright + 1
            end if

            ! Store
            buff(ir)%indx(1) = min(buff(ir)%indx(1),ileft)
            buff(ir)%indx(2) = max(buff(ir)%indx(2),iright)

          end do ! Ranges in data

          ! Skip if not found
          if (buff(ir)%indx(1).gt.buff(ir)%indx(2)) cycle

          !
          ! Prepare input
          !

          ! Limits
          wleft = lambda(buff(ir)%indx(1)) - sig20
          wright = lambda(buff(ir)%indx(2)) + sig20

          ! Initialize
          buff(ir)%indx(3:4) = buff(ir)%indx(1:2)

          ! Left index
          ileft = minloc(abs(lambda - wleft),1)

          ! Move one?
          if (lambda(ileft).gt.wleft) then
            if (ileft.gt.1) ileft = ileft - 1
          end if

          ! Right index
          iright = minloc(abs(lambda - wright),1)

          ! Move one?
          if (lambda(iright).lt.wright) then
           if (iright.lt.nfreq) iright = iright + 1
          end if

          ! Store
          buff(ir)%indx(3) = ileft
          buff(ir)%indx(4) = iright

        ! Not gaussian
        else

          ! Open file
          open (200,file=trim(files(ir)%str), &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='read', form='unformatted')

          ! Read size of array
          read(200,err=1100) buff(ir)%nfreq

          ! Allocate space for wavelength and kernel
          allocate(buff(ir)%wave(buff(ir)%nfreq))
          allocate(buff(ir)%kernel(buff(ir)%nfreq))
          MRAMc = MRAMc + 1d-6*sizeof(buff(ir)%wave)
          MRAMc = MRAMc + 1d-6*sizeof(buff(ir)%kernel)

          ! Read wavelengths and kernel
          read(200,err=1100) buff(ir)%wave
          read(200,err=1100) buff(ir)%kernel

          ! Close file
          close(200)

          ! No sigma, skip
          if (maxval(buff(ir)%kernel).le.0d0) cycle

          ! First
          buff(ir)%kernel(1) = buff(ir)%kernel(1)*0.5d0* &
                               (buff(ir)%wave(2) - buff(ir)%wave(1))

          ! Last
          buff(ir)%kernel(buff(ir)%nfreq) = &
                              buff(ir)%kernel(buff(ir)%nfreq)*0.5d0* &
                              (buff(ir)%wave(buff(ir)%nfreq) - &
                               buff(ir)%wave(buff(ir)%nfreq-1))
          ! Middle
          do ifreq=2,buff(ir)%nfreq-1

            buff(ir)%kernel(ifreq) = buff(ir)%kernel(ifreq)*0.5d0* &
                                     (buff(ir)%wave(ifreq+1) - &
                                      buff(ir)%wave(ifreq-1))

          end do

          ! Normalize
          buff(ir)%kernel = buff(ir)%kernel/sum(buff(ir)%kernel)

          ! Wavelength limits
          wleft = buff(ir)%doub(1)
          wright = buff(ir)%doub(2)

          ! Look over wavelength ranges in data
          do ir1=1,Inf_Stokes%Num_Range

            ! Shorter limits
            irl = Inf_Stokes%Range(ir1,1)
            irr = Inf_Stokes%Range(ir1,2)
            wileft = omega_in(irl)
            wiright = omega_in(irr)

            ! Out of bounds comparisons
            left = wright.lt.wileft
            right = wleft.gt.wiright

            ! Skip if out of bounds
            if (left.or.right) cycle

            !
            ! Prepare output
            !

            ! Left index
            ileft = minloc(abs(lambda - wileft),1)

            ! Move one?
            if (lambda(ileft).gt.wileft) then
              if (ileft.gt.irl.and.ileft.gt.1) ileft = ileft - 1
            end if

            ! Right index
            iright = minloc(abs(lambda - wiright),1)

            ! Move one?
            if (lambda(iright).lt.wiright) then
              if (iright.lt.irr.and.iright.lt.nfreq) &
                iright = iright + 1
            end if

            ! Store
            buff(ir)%indx(1) = min(buff(ir)%indx(1),ileft)
            buff(ir)%indx(2) = max(buff(ir)%indx(2),iright)

          end do ! Ranges in data

          ! If no outputs, skip
          if (buff(ir)%indx(1).gt.buff(ir)%indx(2)) then

            ! Free memory and continue with next
            MRAMc = MRAMc - 1d-6*sizeof(buff(ir)%wave)
            MRAMc = MRAMc - 1d-6*sizeof(buff(ir)%kernel)
            deallocate(buff(ir)%wave,buff(ir)%kernel)
            cycle

          end if

          !
          ! Prepare input
          !

          ! Limits
          wleft = lambda(buff(ir)%indx(1)) + minval(buff(ir)%wave)
          wright = lambda(buff(ir)%indx(2)) + maxval(buff(ir)%wave)

          ! Initialize
          buff(ir)%indx(3:4) = buff(ir)%indx(1:2)

          ! Left index
          ileft = minloc(abs(lambda - wleft),1)

          ! Move one?
          if (lambda(ileft).gt.wleft) then
            if (ileft.gt.1) ileft = ileft - 1
          end if

          ! Right index
          iright = minloc(abs(lambda - wright),1)

          ! Move one?
          if (lambda(iright).lt.wright) then
           if (iright.lt.nfreq) iright = iright + 1
          end if

          ! Store
          buff(ir)%indx(3) = ileft
          buff(ir)%indx(4) = iright

          ! Flag no initialized
          buff(ir)%toinit = .True.

        end if ! Type of kernel

      end do ! PSF ranges

      ! Free
      deallocate(lambda)
      if (allocated(files)) deallocate(files)

      return

1000  write(umsg,'(A,1x,i2)') 'Error opening PSF file',ir
      urou = 'set_psf_ranges'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  write(umsg,'(A,1x,i2)') 'Error reading PSF file',ir
      urou = 'set_psf_ranges'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine set_psf_ranges

!#####################################################################
!#####################################################################
!#####################################################################

      end module psf_mod
