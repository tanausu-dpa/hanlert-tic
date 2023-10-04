      !> Postprocessing out emerging Stokes profiles
      module psf_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/23/2023
!  Last version:
!     10/04/2023 V3.1.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/04/2023:    V3.1.4 - Bugfix: the deallocations of unused
!                             arrays in the non-gaussian PSF case
!                             were being done in the gaussian case,
!                             resulting in the deallcoation of not
!                             allocated arrays (TdPA)
!
!     09/19/2023:    V3.1.3 - Bugfix: Wrong index for the Gaussian
!                             FWHM (HL)
!
!     08/25/2023:    V3.1.2 - Bugfix: Wrong last index when
!                             normalizing PSF kernel (TdPA)
!
!     08/24/2023:    V3.1.1 - Added a skip if there are no output
!                             frequencies for a PSF range (TdPA)
!
!     08/24/2023:    V3.1.0 - Added the possibility of non-gaussian
!                             PSF which are to be read from a
!                             provided file path (TdPA)
!
!     07/03/2023:    V3.0.2 - Removed PSF_StokesI and generalized
!                             PSF_Stokes instead (TdPA)
!                           - Changed the approach to the range
!                             limitation for the PSF application
!                             by defining input and output index
!                             limits (TdPA)
!                           - Added cut_stk in order to crop the
!                             emergent Stokes profiles to store
!                             only what is being inverted (TdPA)
!                           - Added set_psf_ranges to define the
!                             limits in index for the application
!                             of the PSF convolution (TdPA)
!
!     06/13/2023:    V3.0.1 - Update to reduce the calling of error
!                             function (HL)
!                           - Update for multi waveleng ranges (HL)
!                           - The wavelength in Sol changed to nm (HL)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/23/2023:    V0.0.0 - Started from 12/05/2020
!                             TIC@pdf_mod.f90 revision (TdPA)
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
!    PSF_Stokes:
!      Apply spectral convolution with a gaussian to Stokes profiles
!
!    Profiles_Out:
!      Postprocessing of the emergent Stokes profiles
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

      !> Apply spectral convolution with a gaussian to Stokes
      !! profiles\n
      !!  FWHM(FWHM_helper_class): Structure with ranges and PSF\n
      !!               n(integer): Dimension of wavelength axis\n
      !!            nstk(integer): Last index Stokes dimension\n
      !! Stokes(double(0:nstk,n)): Stokes parameters\n
      !!        lambda(double(n)): Wavelength axis
      subroutine PSF_Stokes(FWHM,n,nstk,Stokes,lambda)

      ! IO
      type(FWHM_helper_class), dimension(:), intent(inout):: FWHM
      integer, intent(in):: n,nstk
      double precision, dimension(n), intent(in):: lambda
      double precision, dimension(0:nstk,n), intent(inout):: Stokes

      ! Local
      integer:: i, j, k, iran, il, ir, inl, inr, ifreq0
      integer:: lfirst, last
      double precision:: Sigma, Sigma2, Sigma20, dt, dt2
      double precision, dimension(0:nstk,n):: Stokes_tmp
      double precision, dimension(:), allocatable:: twave
      double precision, dimension(:,:), allocatable:: aStokes

      ! Get wavelength axis and reversed copy of Stokes
      Stokes_tmp = Stokes

      ! For each FWHM range
      do iran=1,FWHM(1)%nn

        ! Skip if not to apply
        if (FWHM(iran)%indx(2).lt.1) cycle

        ! If gaussian
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

          ! Initialize Stokes output
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
            end if

            ! For all wavelengths (except last)
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

            ! Allocate space for auxiliars to interpolation
            allocate(FWHM(iran)%indx1(FWHM(iran)%nfreq,il:ir))
            allocate(FWHM(iran)%indx2(FWHM(iran)%nfreq,il:ir))
            allocate(FWHM(iran)%idx(FWHM(iran)%nfreq,il:ir))

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

      !> Process the emergent Stokes profiles\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!      Sol(Solution_class): Class with the data of the RT
      !!                           solution\n
      !!   e_stk(double(:,:,:,:)): Emergent Stokes parameters\n
      !!             lpe(logical): If polarization\n
      !! buff(IO_helper_class(:)): Info about what to store\n
      !!  FWHM(FWHM_helper_class): Structure with ranges and PSF
      subroutine Profiles_Out(Frec,Sol,e_stk,lpe,buff,FWHM)

      ! IO
      type(IO_helper_class), intent(in):: buff
      type(FWHM_helper_class), dimension(:), allocatable, &
                                                  intent(inout):: FWHM
      type(Frequency_class), intent(in):: Frec
      type(Solution_class), intent(inout):: Sol
      logical, intent(in):: lpe
      double precision, dimension(:,:,:,:), allocatable:: e_stk

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

          nstk = 3

        ! Intensity
        else

          nstk = 0

        end if

        ! Get wavelength axis
        allocate(lambda(nfreq))
        lambda = 1d2/Frec%omega(nfreq:1:-1)

        ! Copy Stokes
        allocate(c_stk(0:nstk,nfreq,1,1))
        c_stk = e_stk(0:nstk,nfreq:1:-1,1:1,1:1)

        ! If FWHM
        if (PSF) call PSF_Stokes(FWHM,nfreq,nstk,c_stk,lambda)

        ! Interpolate into observations axis
        call Intpol_Lin_stk(lambda, c_stk, 4, nstk+1, &
                            nfreq, Sol%omega_input, &
                            Sol%Stokes_out, Sol%Num_Wavelength)

        ! Cut Stokes
        if (cut) then

          ! Copy again
          c_stk = e_stk(0:nstk,:,1:1,1:1)

          ! Cut
          call cut_stk(e_stk,c_stk,buff,nstk)

        end if

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

      !> Process the emergent Stokes profiles\n
      !!   o_stk(double(:,:,:,:)): Cut emergent Stokes parameters\n
      !!   i_stk(double(:,:,:,:)): Full emergent Stokes parameters\n
      !! buff(IO_helper_class(:)): Info about what to store\n
      !!            nstk(integer): Polarization last index
      subroutine cut_stk(o_stk,i_stk,buff,nstk)

      ! IO
      type(IO_helper_class), intent(in):: buff
      integer, intent(in):: nstk
      double precision, dimension(:,:,:,:), intent(in):: i_stk
      double precision, dimension(:,:,:,:), allocatable:: o_stk

      ! Local
      integer:: iran,ii,i0,i1,nn


      ! Reallocate
      deallocate(o_stk)
      allocate(o_stk(0:nstk,buff%nn,1,1))

      ! Initialize buffer
      ii = 0

      ! For each entry to write
      do iran=1,buff%nran

        ! Atom and transition
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

      !> Prepare the ranges to apply the PSF later\n
      !!   Inf_Stokes(Stokes_class): Structure with the Stokes data\n
      !! buff(FWHM_helper_class(:)): Info about ranges\n
      !!     files(strarr_class(:)): Kernel filenames\n
      !!        omega_in(double(:)): Frequency axis from data\n
      !!           omega(double(:)): Frequency axis for RT
      subroutine set_psf_ranges(Inf_Stokes,buff,files,omega_in,omega)

      ! I/O
      type(Stokes_class), intent(in):: Inf_Stokes
      type(FWHM_helper_class), dimension(:), allocatable:: buff
      type(strarr_class), dimension(:), allocatable:: files
      double precision, dimension(:), intent(in):: omega_in,omega

      ! Local
      logical:: left,right

      integer:: ir,ir1,irl,irr,ileft,iright,ifreq,ios

      double precision:: sig20,wleft,wright,wileft,wiright
      double precision, dimension(:), allocatable:: lambda

      ! If no PSF, skip
      if (.not.allocated(buff)) return

      ! Allocate indexes (1 and 2 for in, 3 and 4 for out)
      do ir=1,buff(1)%nn
        allocate(buff(ir)%indx(4))
        buff(ir)%indx(1) = nfreq + 1
        buff(ir)%indx(3) = buff(ir)%indx(1)
        buff(ir)%indx(2) = 0
        buff(ir)%indx(4) = 0
      end do

      ! Get wavelength axis
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

          ! Allocate
          allocate(buff(ir)%wave(buff(ir)%nfreq))
          allocate(buff(ir)%kernel(buff(ir)%nfreq))

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
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  write(umsg,'(A,1x,i2)') 'Error reading PSF file',ir
      urou = 'set_psf_ranges'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine set_psf_ranges

!#####################################################################
!#####################################################################
!#####################################################################

      end module psf_mod
