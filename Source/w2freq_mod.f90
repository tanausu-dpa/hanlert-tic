      !> Unit conversions for wavelength and profiles
      module w2freq_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/17/2023
!  Last version:
!     05/20/2024 V3.0.12
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     05/20/2024:   V3.0.12 - Added an additional way of controling
!                             the weights in the inversion by
!                             introducing enhancing factors over the
!                             existing inputs (TdPA)
!
!     05/17/2024:   V3.0.11 - The number of degrees of freedom is
!                             now counted when defining the weights,
!                             although it is not used anymore (TdPA)
!                           - New normalization for the weights when
!                             creating them to keep the order of
!                             the merit function in the expected
!                             ranges regardless of the weights in the
!                             user input (TdPA)
!
!     10/04/2023:   V3.0.10 - Verbosity update (TdPA)
!
!     09/18/2023:    V3.0.9 - Verbosity update (TdPA)
!
!     08/11/2023:    V3.0.8 - Bugfix: Verbosity directed to the
!                             wrong files (TdPA)
!
!     08/11/2023:    V3.0.7 - Added inversion_weights to manage the
!                             definition of the inversion weights
!                             accounting for specified ranges or
!                             for a weights file (TdPA)
!
!     07/03/2023:    V3.0.6 - Update errors to parallel (TdPA)
!                           - In Profile_conversion the wavelengths
!                             ranges are copied into the
!                             Solution_class structure (TdPA)
!                           - Added back the allocation of the
!                             scales (TdPA)
!                           - Updated verbosity of ranges (TdPA)
!                           - Added management of diffuse light
!                             profile (TdPA)
!
!     06/12/2023:    V3.0.5 - Removed the allocation of scales and
!                             wavelength range (HL)
!                           - Rename the Wavelength_Conversion since
!                             we will not reverse the wavelength (HL)
!
!     05/16/2023:    V3.0.4 - Bugfix: Sigma must be always scaled if
!                             it is flagged, and not only if the
!                             input was wavelength dependent (TdPA)
!
!     04/26/2023:    V3.0.3 - Made the step to consider different
!                             ranges a parameter in this module (TdPA)
!                           - Bugfix: There was a typo in the
!                             automatic definition of weights for
!                             polarization (TdPA)
!
!     04/11/2023:    V3.0.2 - Update for multi-wavelength ranges (HL)
!
!     03/15/2023:    V3.0.1 - Removed unused Check_Profile (TdPA)
!                           - Removed some commented lines (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/17/2023:    V0.0.0 - Started from 05/12/2020
!                             TIC@w2freq_mod.f90 revision from
!                             Hao (TdPA)
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
!    Range_Check:
!      Check the wavelength ranges
!
!    inversion_weights:
!      Prepare the wavelength dependent weights for the inversion from
!      the inputs
!
!    Profile_Conversion
!      Convert the Stokes profiles from HanleRT units to SI
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use commons_mod
      use parameters_mod, only: c
      use types_mod

      double precision, parameter:: lstep = 50d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate the wavelength ranges from the input\n
      !!         Lambda(double(:)): Wavelength axis\n
      !!  Inf_Stokes(Stokes_class): Structure with the Stokes data
      subroutine Range_Check(Lambda,Inf_Stokes)

      ! IO
      type(Stokes_class), intent(inout):: Inf_Stokes
      double precision, dimension(:), allocatable, intent(inout):: &
                                                                Lambda

      ! Local
      integer:: i, j
      integer, dimension(20):: tmp

      ! Initialize range
      Inf_Stokes%Num_Range = 1

      !
      ! Find wavelength ranges
      !

      ! For each wavelength
      do i=1,Inf_Stokes%Num_Wavelength-1

        ! Check if separation larger than step
        if (abs(lambda(i+1)-lambda(i)).gt.lstep) then

          ! Store the split in temporal
          tmp(Inf_Stokes%Num_Range) = i

          ! Add a range
          Inf_Stokes%Num_Range = Inf_Stokes%Num_Range+1

        end if ! If larger jump that prefixed step

      end do ! Wavelengths

      ! If ranges not allocated, do it
      if (.not.allocated(Inf_Stokes%Range)) &
        allocate(Inf_Stokes%Range(Inf_Stokes%Num_Range,2))

      ! If more than one range found
      if (Inf_Stokes%Num_Range.gt.1) then

        ! Set final split (the end)
        tmp(Inf_Stokes%Num_Range) = Inf_Stokes%Num_Wavelength

        ! Initialize counter
        j=1

        ! For each range
        do i=1,Inf_Stokes%Num_Range

          ! Set limits
          Inf_Stokes%Range(i,1) = j
          Inf_Stokes%Range(i,2) = tmp(i)

          ! Update running index
          j = tmp(i) + 1

        end do ! Ranges

      ! Just one range
      else

        ! Full range
        Inf_Stokes%Range(1,1) = 1
        Inf_Stokes%Range(1,2) = Inf_Stokes%Num_Wavelength

      end if ! Number of ranges

      return

      end subroutine Range_Check

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate the wavelength ranges from the input\n
      !!             Input(Input): Structure with settings data\n
      !! Inf_Stokes(Stokes_class): Structure with the Stokes data\n
      !!      omega_in(double(:)): Frequency axis from data\n
      subroutine inversion_weights(Input,Inf_Stokes,omega_in)

      ! IO
      type(Input_class):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      double precision, dimension(:), intent(in):: omega_in

      ! Local
      logical:: left, right

      integer:: i, j, il, ir, nl, nstk, ios
      integer:: ileft, iright

      double precision:: wileft, wiright, wleft, wright


      !
      ! Sanity
      if (Inf_Stokes%auto_weight) return

      !
      ! If factors
      if (allocated(Input%Weight_Factor)) then

        ! Run over all entries
        do i=1,size(Input%Weight_Factor,2)-1

          ! Flag that limits are to be check
          left = .True.

          ! Run over all other entry
          do j=i+1,size(Input%Weight_Factor,2)

            ! If Stokes is not the same, skip
            if (nint(Input%Weight_Factor(1,i)).ne. &
                nint(Input%Weight_Factor(1,j))) cycle

            ! Get index limits for first if not yet
            if (left) then

              ! Get indexes
              il = minloc(abs(omega_in - Input%Weight_Factor(2,i)),1)
              ir = minloc(abs(omega_in - Input%Weight_Factor(3,i)),1)

              ! Flag that found
              left = .False.

            end if ! Need to get indexes for first

            ! Get indexes for second
            nl = minloc(abs(omega_in - Input%Weight_Factor(2,j)),1)
            nstk = minloc(abs(omega_in - Input%Weight_Factor(3,j)),1)

            ! Check superposition
            if ((nl.ge.il.and.nl.le.ir).or. &
                (nstk.ge.il.and.nstk.le.ir)) then

              ! Abort if superposition
              select case (nint(Input%Weight_Factor(1,i)))

                case (0)

                  ! Abort
                  write(umsg, &
                        '(A,2(A,i2,"(",es15.8,",",es15.8,")"))') &
                         'There is superposition of at least two '// &
                         'ranges in WEIGHT_FACTOR for Stokes I: ', &
                         'range ',i,omega_in(il),omega_in(ir), &
                         ' and range ',j,omega_in(nl),omega_in(nstk)

                case (1)

                  ! Abort
                  write(umsg, &
                        '(A,2(A,i2,"(",es15.8,",",es15.8,")"))') &
                         'There is superposition of at least two '// &
                         'ranges in WEIGHT_FACTOR for Stokes Q: ', &
                         'range ',i,omega_in(il),omega_in(ir), &
                         ' and range ',j,omega_in(nl),omega_in(nstk)

                case (2)

                  ! Abort
                  write(umsg, &
                        '(A,2(A,i2,"(",es15.8,",",es15.8,")"))') &
                         'There is superposition of at least two '// &
                         'ranges in WEIGHT_FACTOR for Stokes U: ', &
                         'range ',i,omega_in(il),omega_in(ir), &
                         ' and range ',j,omega_in(nl),omega_in(nstk)

                case (3)

                  ! Abort
                  write(umsg, &
                        '(A,2(A,i2,"(",es15.8,",",es15.8,")"))') &
                         'There is superposition of at least two '// &
                         'ranges in WEIGHT_FACTOR for Stokes V: ', &
                         'range ',i,omega_in(il),omega_in(ir), &
                         ' and range ',j,omega_in(nl),omega_in(nstk)

                case default

                  ! Abort
                  write(umsg, &
                        '(A,2(A,i2,"(",es15.8,",",es15.8,")"))') &
                         'There is superposition of at least two '// &
                         'ranges in WEIGHT_FACTOR: ', &
                         'range ',i,omega_in(il),omega_in(ir), &
                         ' and range ',j,omega_in(nl),omega_in(nstk)
              end select

              ! Finish the error
              urou = 'inversion_weights'
              call abortedS(umsg,urou,-1,.True.,.True.)
              call control
              return

            end if

          end do ! All other entries
        end do ! All entries

      end if ! If weight factors

      !
      ! If weight from file
      !
      if (Input%linv_weight) then

        ! Routine name
        urou = 'inversion_weights'

        ! Prepare variable
        allocate(Inf_Stokes%weight(0:3,Inf_Stokes%Num_wavelength))

        ! If thermal, initialize to 0
        if (Input%Type_Inversion.eq.0) &
          Inf_Stokes%weight(1:3,:) = 0d0

        !
        ! Master reads
        !
        if (gpid.eq.0) then

          ! Open
          open (200,file=trim(Input%inv_weight), &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='read', form='unformatted')

          ! Read wavelengths
          read(200,err=1100) nl

          ! Sanity
          if (Inf_Stokes%Num_wavelength.ne.nl) then

            ! Abort
            umsg = 'number of wavelength in weight file do not '// &
                   'match the data'
            urou = 'inversion_weights'
            call abortedS(umsg,urou,-1,.True.,.True.)
            call control
            return

          end if

          ! Read Stokes
          read(200,err=1100) nstk

          ! Sanity
          if (nstk.ne.1.and.nstk.ne.4) then

            ! Abort
            write(umsg,'(A,i4)') &
                   'the number of Stokes parameters in the '// &
                   'weight file must be 1 or 4, but it is ',nstk
            urou = 'inversion_weights'
            call abortedS(umsg,urou,-1,.True.,.True.)
            call control
            return

          end if

          ! Sanity
          if (nstk.eq.1.and.Input%Type_Inversion.gt.0) then

            ! Abort
            umsg = 'you are loading weights for only intensity '// &
                   'for an inversion with polarization'
            urou = 'inversion_weights'
            call abortedS(umsg,urou,-1,.True.,.True.)
            call control
            return

          end if

          ! Read weights for each Stokes parameter
          do i=0,nstk-1
            read(200,err=1100) Inf_Stokes%weight(i,:)
          end do

          ! Close
          close(200)

        end if

        ! Control
        call control
        if (laborted) return

        ! Share
        call MPI_BCAST(Inf_Stokes%weight(0,1), &
                       4*Inf_Stokes%num_wavelength, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      !
      ! If weight from numbers in input
      !
      else

        ! Allocate and  initialize
        allocate(Inf_Stokes%weight(0:3,Inf_Stokes%Num_wavelength))
        Inf_Stokes%weight = 0d0

        ! For each weight range
        do i=1,Input%Num_weight

          ! Wavelength limits
          wleft  = Input%weight(4,i)
          wright = Input%weight(5,i)

          ! Look over wavelength ranges in data
          do j=1,Inf_Stokes%Num_Range

            ! Shorter limits
            il = Inf_Stokes%Range(j,1)
            ir = Inf_Stokes%Range(j,2)
            wileft = omega_in(il)
            wiright = omega_in(ir)

            ! Out of bounds comparisons
            left = wright.lt.wileft
            right = wleft.gt.wiright

            ! Skip if out of bounds
            if (left.or.right) cycle

            !
            ! Prepare actual range
            !

            ! Left index
            ileft = minloc(abs(omega_in - wileft),1)

            ! Move one?
            if (omega_in(ileft).gt.wileft) then
              if (ileft.gt.il.and.ileft.gt.1) ileft = ileft - 1
            end if

            ! Right index
            iright = minloc(abs(omega_in - wiright),1)

            ! Move one?
            if (omega_in(iright).lt.wiright) then
              if (iright.lt.ir.and.iright.lt.nfreq) &
                iright = iright + 1
            end if

            ! Store
            do nstk=0,3
              Inf_Stokes%weight(nstk,ileft:iright) = &
                                                  Input%weight(nstk,i)
            end do ! Stokes

          end do ! Ranges in data
        end do ! Weight ranges

      end if ! Type of input


      !
      ! Factors
      if (allocated(Input%Weight_Factor)) then

        ! Run over all entries
        do i=1,size(Input%Weight_Factor,2)

          ! Get indexes
          nstk = nint(Input%Weight_Factor(1,i))
          il = minloc(abs(omega_in - Input%Weight_Factor(2,i)),1)
          ir = minloc(abs(omega_in - Input%Weight_Factor(3,i)),1)

          ! Apply weight
          Inf_Stokes%weight(nstk,il:ir) = Input%Weight_Factor(4,i)* &
                                         Inf_Stokes%weight(nstk,il:ir)
        end do ! Factor entries

        ! Free
        deallocate(Input%Weight_Factor)

      end if ! There are input factors


      !
      ! Sanity check
      !
      if (minval(Inf_Stokes%weight).lt.0d0) then

        ! Abort
        umsg = 'Found a negative weight for the inversion'
        urou = 'inversion_weights'
        call aborted
        return

      end if

      !
      ! Normalization for the weights
      Inf_Stokes%Weight = Inf_Stokes%Weight/ &
                          sqrt(sum(Inf_Stokes%Weight* &
                                   Inf_Stokes%Weight))

      ! Initialize number of data points
      Inf_Stokes%Num_freedomI = 0
      Inf_Stokes%Num_freedom = 0

      ! Wavelengths
      do j=1,Inf_Stokes%Num_wavelength

        ! If non-zero intensity weight, add to count
        if (Inf_Stokes%Weight(0,j).gt.0d0) then

          ! Total
          Inf_Stokes%Num_freedom = 1 + Inf_Stokes%Num_freedom
          Inf_Stokes%Num_freedomI = 1 + Inf_Stokes%Num_freedomI

        end if ! Non-zero weight

        ! For each Stokes parameter
        do i=1,3

          ! If non-zero weight, add to count
          if (Inf_Stokes%Weight(i,j).gt.0d0) &
            Inf_Stokes%Num_freedom = 1 + Inf_Stokes%Num_freedom

        end do ! Stokes parameters
      end do ! Wavelengths

      ! Free
      if (allocated(Input%weight)) deallocate(Input%weight)

      return

1000  write(umsg,'(A,1x,i2)') 'Error opening weight file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  write(umsg,'(A,1x,i2)') 'Error reading weight file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine inversion_weights

!#####################################################################
!#####################################################################
!#####################################################################

      !> Convert the Stokes profiles from HanleRT units to SI\n
      !! Inf_Stokes(Stokes_class): Structure with the Stokes data\n
      !!      Sol(Solution_class): Structure with the solution data
      subroutine Profile_Conversion(Inf_Stokes,Sol)

      ! IO
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Solution_class), intent(inout):: Sol

      ! Local
      integer:: i, j

      ! Copy ranges
      Sol%Num_Range = Inf_Stokes%Num_Range
      if (.not.allocated(Sol%Range)) &
        allocate(Sol%Range(Inf_Stokes%Num_Range,2))
      Sol%Range = Inf_Stokes%Range

      ! Unit convesion
      Inf_Stokes%Stokes_Ob = Inf_Stokes%Stokes_Ob*c*1d14

      ! Q (rotate 90º)
      Inf_Stokes%Stokes_Ob(1,:) = -Inf_Stokes%Stokes_Ob(1,:)

      ! U (rotate 90º)
      Inf_Stokes%Stokes_Ob(2,:) = -Inf_Stokes%Stokes_Ob(2,:)

      ! Allocate scales if not already there
      if (.not.allocated(Inf_Stokes%Scales)) &
        allocate(Inf_Stokes%Scales(Inf_Stokes%Num_Range,0:3))
      if (.not.allocated(Sol%Scal_Stokes)) &
        allocate(Sol%Scal_Stokes(Inf_Stokes%Num_Range))

      ! For each wavelength range
      do i=1,Inf_Stokes%Num_Range
        ! For each Stokes parameter
        do j = 0, 3

          ! Set scales
          Inf_Stokes%Scales(i,j) = &
            sum(abs(Inf_Stokes%Stokes_Ob(j,Sol%Range(i,1): &
                                           Sol%Range(i,2))))/ &
              dble(Sol%Range(i,2)-Sol%Range(i,1)+1)*1d-2

        end do
      end do

      ! Store intensity scale
      Sol%Scal_Stokes = Inf_Stokes%Scales(:,0)

      ! If more than one range
      if (Sol%Num_Range.gt.1) then

        ! Master
        if (pid.eq.0) then

          ! Header
          write(umsg,'(A)') ' * Setting-up weights:'
          call verboseI(3)

          ! Run over ranges
          do i=1,Sol%Num_Range

            ! Verbose wavelength ranges
            write(umsg,'(A,i4,3x,2F15.4)')  &
              '   Wavelength range ', &
              i,Sol%omega_input(Sol%Range(i,2)), &
                Sol%omega_input(Sol%Range(i,1))
            call verboseI(3)

            ! Verbose the weight
            write(umsg,'(A, 4F10.4,3x, A, es13.5)')  &
              '     Weights at first wavelength in range = ', &
              Inf_Stokes%weight(0:3,Sol%Range(i,1)), &
              ' scaling factor = ',Sol%Scal_Stokes(i)
            call verboseI(3)

          end do ! Ranges

        end if ! Master
      end if ! More than one range

      ! If weights are automatic
      if (Inf_Stokes%auto_weight) then

        ! Initialize to 1
        Inf_Stokes%Weight(1,0) = 1d0

        ! If more than one range
        if (Inf_Stokes%Num_Range.gt.1) then

          ! For each range
          do i=2,Inf_Stokes%Num_Range

            ! Set weight
            Inf_Stokes%Weight(i,0) = Inf_Stokes%Scales(1,0)/ &
                                     Inf_Stokes%Scales(i,0)
          end do

        end if ! Multiple ranges

        ! For each range
        do i=Inf_Stokes%Num_Range,1,-1

          ! Check reductions in Q of a factor 2d3
          if (Inf_Stokes%Scales(i,0).gt.Inf_Stokes%Scales(i,1)*2d3) &
            Inf_Stokes%Scales(i,1) = Inf_Stokes%Scales(i,0)/2d3

          ! Check reductions in U of a factor 2d3
          if (Inf_Stokes%Scales(i,0).gt.Inf_Stokes%Scales(i,2)*2d3) &
            Inf_Stokes%Scales(i,2) = Inf_Stokes%Scales(i,0)/2d3

          ! Check reductions in V of a factor 2d3
          if (Inf_Stokes%Scales(i,0).gt.Inf_Stokes%Scales(i,3)*2d3) &
            Inf_Stokes%Scales(i,3) = Inf_Stokes%Scales(i,0)/2d3

          !
          ! Compute automatic weights
          !

          ! Q
          Inf_Stokes%Weight(i, 1) = 0.5d0*Inf_Stokes%Scales(i,0)* &
                                    (1d0/Inf_Stokes%Scales(i,1) + &
                                     1d0/Inf_Stokes%Scales(i,2))
          ! U
          Inf_Stokes%Weight(i, 2) = Inf_Stokes%Weight(i, 1)
          ! V
          Inf_Stokes%Weight(i, 3) = Inf_Stokes%Scales(i,0)/ &
                                    Inf_Stokes%Scales(i,3)

          ! Master notify
          if (pid.eq.0) then
            write(umsg,'(A,i4,4(1x,es20.3))') &
              '   o Automatic weights for range',i, &
              Inf_Stokes%weight(i,:)
            call verbosev
          end if ! Master

        end do

      end if ! Automatic weights

      ! If fractional Stokes
      if (Sol%Fractional) then

        ! Q/I
        Inf_Stokes%Stokes_Ob(1,:) = Inf_Stokes%Stokes_Ob(1,:)/ &
                                    Inf_Stokes%Stokes_Ob(0,:)*1d2
        ! U/I
        Inf_Stokes%Stokes_Ob(2,:) = Inf_Stokes%Stokes_Ob(2,:)/ &
                                    Inf_Stokes%Stokes_Ob(0,:)*1d2
        ! V/I
        Inf_Stokes%Stokes_Ob(3,:) = Inf_Stokes%Stokes_Ob(3,:)/ &
                                    Inf_Stokes%Stokes_Ob(0,:)*1d2

        ! For each range
        do i = 1,Sol%Num_Range

          ! I/I_scale
          Inf_Stokes%Stokes_Ob(0,Sol%Range(i,1):Sol%Range(i,2)) = &
              Inf_Stokes%Stokes_Ob(0,Sol%Range(i,1):Sol%Range(i,2))/ &
              Sol%Scal_Stokes(i)

        end do

      ! If not fractional
      else

        ! For each range
        do i = 1,Sol%Num_Range

          ! Just scale
          Inf_Stokes%Stokes_Ob(0:3,Sol%Range(i,1):Sol%Range(i,2)) = &
            Inf_Stokes%Stokes_Ob(0:3,Sol%Range(i,1):Sol%Range(i,2))/ &
            Sol%Scal_Stokes(i)

        end do

      end if ! Fractional

      ! If there are wavelength dependent sigmas
      if (Inf_Stokes%Sigma_Flag) then

        ! If fractional Stokes
        if (Sol%Fractional) then

          ! For each range
          do i=1,Sol%Num_Range

            ! Scale for intensity
            Inf_Stokes%Sigma_W(0,Sol%Range(i,1):Sol%Range(i,2)) = &
              Inf_Stokes%Sigma_W(0,Sol%Range(i,1):Sol%Range(i,2))* &
              c*1d14/Sol%Scal_Stokes(i)

          end do

          ! Convert to consistent scale for Stokes
          do i=1,3
            Inf_Stokes%Sigma_W(i,:) = Inf_Stokes%Sigma_W(i,:)*1d2
          end do

        ! No fractional
        else

          ! For each range
          do i=1,Sol%Num_Range

            ! Scale
            Inf_Stokes%Sigma_W(0:3,Sol%Range(i,1):Sol%Range(i,2)) =&
             Inf_Stokes%Sigma_W(0:3,Sol%Range(i,1):Sol%Range(i,2))*&
             c*1d14/Sol%Scal_Stokes(i)

          end do

        end if ! Fractional

      end if

      ! If there is diffuse light
      if (Inf_Stokes%Diff_flag) then

        ! Change units
        Sol%Stokes_diff = Sol%Stokes_diff*c*1d14

      end if

      return

      end subroutine Profile_Conversion

!#####################################################################
!#####################################################################
!#####################################################################

      end module w2freq_mod
