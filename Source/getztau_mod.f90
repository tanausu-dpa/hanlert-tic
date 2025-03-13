      !> Background opacities
      module getztau_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     03/12/2019
!  Last version:
!     12/03/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     12/03/2025:    V4.0.1 - Explicitly account for changes in the
!                             memory allocated in Atmo_class (TdPA)
!                           - Only consider the limitation for PRD
!                             if there is PRD (TdPA)
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
!  getztau
!    Computes height or tau scale given the other
!
!  find_tau
!    Find the height indexes that restrict the calculation to the
!  specified optical depth range
!
!  find_z
!    Find the height indexes that restrict the calculation to the
!  specified height range
!
!  restrict_zaxis
!    Find the indexes to limit the height axis in the RT problem from
!  the specified input in heights and/or continuum optical depth
!
!  restrict_LTE_lines
!    Find the index to limit the height axis in the RT problem for a
!  LTE line
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes height or tau scale given the other\n
      !!    Atmo(Atmo_class): Structure with atmospheric data\n
      !!  can_leave(logical): Indicate if slaves can leave
      subroutine getztau(Atmo,can_leave)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      logical, intent(in):: can_leave

      ! Local

      integer:: iz,iz0

      double precision:: dz,dt,mint,chi

      ! If MPI, master needs to receive chi500
      if (nproc.gt.1) then

        ! Master
        if (pid.eq.0) then

          ! Allocate chi500
          if (.not.allocated(Atmo%chi500)) then
            allocate(Atmo%chi500(nZ))
            MRAMc = MRAMc + 1d-6*sizeof(Atmo%chi500)
          end if

          call MPI_RECV(Atmo%chi500,nZ,MPI_DOUBLE_PRECISION, &
                        1,1,MPI_COMM_RT,MPI_STATUS_IGNORE,ierr)

        ! Slaves
        else

          ! It is enough with the first process
          if (pid.eq.1) then

            ! Send the chi500 data
            call MPI_SEND(Atmo%chi500,nZ,MPI_DOUBLE_PRECISION, &
                          0,1,MPI_COMM_RT,ierr)

          end if

          ! If no tau scale, chi500 is not needed
          if (.not.ztau.and.can_leave) then

            ! Free
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%chi500)
            deallocate(Atmo%chi500)

          end if ! Keeping Atmo%chi500
        end if ! Master or slave
      end if ! MPI

      ! Rest of CPU can leave now
      if (pid.gt.0.and.can_leave) return

      ! If zalt is already allocated
      if (allocated(Atmo%zalt)) then

        ! If not correct size
        if (nZ.ne.size(Atmo%zalt)) then

          ! Clean and allocate again
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%zalt)
          deallocate(Atmo%zalt)
          allocate(Atmo%zalt(nZ))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%zalt)

        end if ! Incorrect size

      ! zalt not allocated
      else

        ! Allocate
        allocate(Atmo%zalt(nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%zalt)

      end if ! zalt status

      ! If tau was given
      if (ztau) then

        ! Initialize mint
        mint = 1d99
        iz0 = -1

        !
        ! Look for node closest to tau=1
        !

        ! For every height starting from the top
        do iz=1,nZ

          ! If this optical depth is closer to 1
          if (abs(Atmo%z(iz)-1d0).lt.abs(mint-1d0)) then

            ! Save
            iz0 = iz
            mint = Atmo%z(iz)

          end if ! Closer to one

          ! Leave if already went beyond one
          if (Atmo%z(iz).gt.1d0) exit

        end do ! Height nodes

        ! Initialize height at the bottom of the atmosphere
        Atmo%zalt(nZ) = 0d0

        ! For each height starting from the second to last at the
        ! bottom boundary
        do iz=nZ-1,1,-1

          ! Optical path to the next point
          dt = Atmo%z(iz+1) - Atmo%z(iz)

          ! Opacity sum
          chi = Atmo%chi500(iz+1) + Atmo%chi500(iz)

          ! Get new height
          ! 10^-5 cm -> km
          Atmo%zalt(iz) = Atmo%zalt(iz+1) + 2d-5*dt/chi

        end do ! Height nodes

        ! Shift height so z=0 ~ tau=1
        Atmo%zalt = Atmo%zalt - Atmo%zalt(iz0)

      ! If z was given
      else

        ! Initialize top at 0 optical depth
        Atmo%zalt(1) = 0d0

        ! For every height going down, add optical depth at reference
        ! frequency
        do iz=2,nZ

          ! Height step
          dz = Atmo%z(iz-1) - Atmo%z(iz)

          ! Optical depth step
          dt = 0.5*dz*(Atmo%chi500(iz) + Atmo%chi500(iz-1))

          ! Add to optical depth
          Atmo%zalt(iz) = Atmo%zalt(iz-1) + dt

        end do ! Height nodes

      end if ! Given tau or height axis

      return

      end subroutine getztau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Find the height indexes that restrict the calculation to the
      !! specified optical depth range\n
      !!  z(double(:)): Optical depth axis\n
      !!  str(logical): If the restriction is strict\n
      !!    r0(double): Minimum value of optical depth to consider\n
      !!    r1(double): Maximum value of optical depth to consider\n
      !!  rz0(integer): Minimum index of the height range to
      !!                consider\n
      !!  rz0(integer): Maximum index of the height range to
      !!                consider
      subroutine find_tau(z,str,r0,r1,rz0,rz1)

      ! I/O

      logical, intent(in):: str
      integer, intent(out):: rz0,rz1
      double precision, intent(in):: r0,r1
      double precision, dimension(:), intent(in):: z

      ! Local

      integer:: iz


      ! Initialize indexes
      rz0 = -1
      rz1 = -1

      ! If the top is already larger than the lower limit,
      ! just specify it
      if (z(1).ge.r0) rz0 = 1

      ! If the bottom is already smaller than the upper limit,
      ! just specify it
      if (z(nz).le.r1) rz1 = nz

      ! For all non-boundary heights going down
      do iz=2,nz-1

        ! If the lower limit is between the current point
        ! and the previous one
        if (z(iz-1).lt.r0.and.z(iz).ge.r0) then

          ! If strict restriction
          if (str) then

            ! Cut at this point
            rz0 = iz

          ! If lax restriction
          else

            ! Cut at previous point
            rz0 = iz - 1

          end if ! Strict or lax restriction
        end if ! If lower limit between current and previous

        ! If the upper limit is between the current point
        ! and the next one
        if (z(iz).le.r1.and.z(iz+1).gt.r1) then

          ! If strict restriction
          if (str) then

            ! Cut at this point
            rz1 = iz

          ! If lax restriction
          else

            ! Cut at previous point
            rz1 = iz + 1

          end if ! Strict or lax restriction
        end if ! If upper limit between current and next

        ! If we have found both limits, leave
        if (rz1.gt.0.and.rz0.gt.0) exit

      end do ! Intermediate heights

      end subroutine find_tau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Find the height indexes that restrict the calculation to the
      !! specified height range\n
      !!  z(double(:)): Height axis\n
      !!  str(logical): If the restriction is strict\n
      !!    r0(double): Minimum value of optical depth to consider\n
      !!    r1(double): Maximum value of optical depth to consider\n
      !!  rz0(integer): Minimum index of the height range to
      !!                consider\n
      !!  rz0(integer): Maximum index of the height range to
      !!                consider
      subroutine find_z(z,str,r0,r1,rz0,rz1)

      ! I/O

      logical, intent(in):: str
      integer, intent(out):: rz0,rz1
      double precision, intent(in):: r0,r1
      double precision, dimension(:), intent(in):: z

      ! Local

      integer:: iz


      ! Initialize indexes
      rz0 = -1
      rz1 = -1

      ! If the bottom is already larger than the lower limit,
      ! just specify it
      if (z(nz).ge.r0) rz1 = nz

      ! If the top is already larger than the lower limit,
      ! just specify it
      if (z(1).le.r1) rz0 = 1

      ! For all non-boundary heights going down
      do iz=2,nz-1

        ! If the upper limit is between the current point
        ! and the previous one
        if (z(iz-1).gt.r1.and.z(iz).le.r1) then

          ! If strict restriction
          if (str) then

            ! Cut at this point
            rz0 = iz

          ! If lax restriction
          else

            ! Cut at previous point
            rz0 = iz - 1

          end if ! Strict or lax restriction
        end if ! If lower limit between current and previous

        ! If the lower limit is between the current point
        ! and the next one
        if (z(iz).ge.r0.and.z(iz+1).lt.r0) then

          ! If strict restriction
          if (str) then

            ! Cut at this point
            rz1 = iz

          ! If lax restriction
          else

            ! Cut at previous point
            rz1 = iz + 1

          end if ! Strict or lax restriction
        end if ! If lower limit between current and previous

        ! If we have found both limits, leave
        if (rz1.gt.0.and.rz0.gt.0) exit

      end do ! Intermediate heights

      end subroutine find_z

!#####################################################################
!#####################################################################
!#####################################################################

      !> Find the indexes to limit the height axis in the RT problem
      !! from the specified input in heights and/or continuum optical
      !! depth\n
      !!    Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Input(Input_class): Structure with configuration data\n
      subroutine restrict_zaxis(Atmo,Input)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ierr,tz0,tz1,zz0,zz1

      double precision:: r0,r1


      ! Initialize full range
      Rz0 = 1
      Rz1 = nz
      Rnz = nz
      Rz1_PRD = nz


      ! If not restricting, leave
      if (.not.(Input%rest_tau.or.Input%rest_z.or. &
                Input%rest_tau_red.or. &
                Input%rest_z_red)) return

      ! If Master
      if (pid.eq.0) then

        ! If restricting by tau
        if (Input%rest_tau) then

          ! Shorten inputs
          r0 = Input%r0tc
          r1 = Input%r1tc

          ! Tau is the main axis
          if (ztau) then

            ! Find indexes of the height range to consider
            call find_tau(Atmo%z,Input%rest_tau_strc,r0,r1,tz0,tz1)

          ! Tau is the secondary axis
          else

            ! Find indexes of the height range to consider
            call find_tau(Atmo%zalt,Input%rest_tau_strc,r0,r1,tz0,tz1)

          end if ! Tau main

          ! If not restricting in Z
          if (.not.Input%rest_z) then

            ! We are done
            Rz0 = tz0
            Rz1 = tz1

          end if ! Not restrict in z
        end if ! Restrict tau

        ! Restrict Height
        if (Input%rest_z) then

          ! Shorten inputs
          r0 = Input%r0z*1d5
          r1 = Input%r1z*1d5

          ! Height is the secondary axis
          if (ztau) then

            ! Find indexes of the height range to consider
            call find_z(Atmo%zalt,Input%rest_z_strc,r0,r1,zz0,zz1)

          ! Height is the main axis
          else

            ! Find indexes of the height range to consider
            call find_z(Atmo%z,Input%rest_z_strc,r0,r1,zz0,zz1)

          end if ! Tau main

          ! If not restricting in tau
          if (.not.Input%rest_tau) then

            ! We are done
            Rz0 = zz0
            Rz1 = zz1

          end if ! Not restrict in tau
        end if ! Restrict Height

        ! If restricting both tau and height
        if (Input%rest_z.and.Input%rest_tau) then

          ! Get the more conservative numbers
          Rz0 = max(zz0,tz0)
          Rz1 = min(zz1,tz1)

        end if ! Restricting both tau and z

        ! Correct PRD lim.
        if (PRD.and.(Input%rest_z_red.or.Input%rest_tau_red)) then

          ! Initialize
          zz1 = Rz1
          tz1 = Rz1

          ! Restricting tau
          if (Input%rest_tau_red) then

            ! Shorten
            r1 = Input%r1tc_prd

            ! Tau is the main axis
            if (ztau) then

              ! Shorten
              r0 = Atmo%z(1)

              ! Find indexes of height range to consider
              call find_tau(Atmo%z,.True.,r0,r1,tz0,tz1)

            ! Tau is the secondary axis
            else

              ! Shorten
              r0 = Atmo%zalt(1)

              ! Find indexes of height range to consider
              call find_tau(Atmo%zalt,.True.,r0,r1,tz0,tz1)

            end if ! Tau main

            ! If not restricting in z, done
            if (.not.Input%rest_z_red) Rz1_PRD = tz1

          end if ! Restrict tau

          ! Restricting tau
          if (Input%rest_z_red) then

            ! Shorten
            r0 = Input%r1z_prd*1d5

            ! Height is the secondary axis
            if (ztau) then

              ! Shorten
              r1 = Atmo%zalt(1)

              ! Find indexes of height range to consider
              call find_z(Atmo%zalt,.True.,r0,r1,zz0,zz1)

            ! Height is the main axis
            else

              ! Shorten
              r1 = Atmo%z(1)

              ! Find indexes of height range to consider
              call find_z(Atmo%z,.True.,r0,r1,zz0,zz1)

            end if ! Tau main

            ! If not restricting in tau, done
            if (.not.Input%rest_tau_red) Rz1_PRD = zz1

          end if ! Restrict Height

          ! If restricting both tau and height
          if (Input%rest_z_red.and.Input%rest_tau_red) then

            ! Get the more conservative numbers
            Rz1_PRD = min(zz1,tz1)

          end if ! Restricting both tau and z
        end if ! Restricting PRD range
      end if ! Master

      ! If MPI
      if (nproc.gt.1) then

        ! Broadcast limits
        call MPI_BCAST(Rz0,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)
        call MPI_BCAST(Rz1,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)
        call MPI_BCAST(Rz1_PRD,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)

      end if ! MPI

      ! Dimension
      Rnz = Rz1 - Rz0 + 1

      ! If smaller than 3, there is a problem
      if (Rnz.lt.3) then

        ! Error message
        urou = 'restrict_zaxis'
        umsg = ' # The resulting atmosphere after the height '// &
               'restriction has less than three nodes'
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,.True.,.False.)

      end if ! Too few nodes

      ! If non-sense
      if (Rz0.lt.1) then

        ! Error message
        urou = 'restrict_zaxis'
        write(umsg,'(A,i3)') ' # The resulting lower limit is '// &
                             'below 1 :',Rz0
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,.True.,.False.)

      end if ! Non-sense index

      ! If non-sense
      if (Rz1.gt.nz) then

        ! Error message
        urou = 'restrict_zaxis'
        write(umsg,'(A,i3,A,i3)') ' # The resulting upper '// &
                                  'limit is above ',nz,' :',Rz1
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,.True.,.False.)

      end if ! Non-sense index

      ! If non-sense
      if (Rz1_PRD.gt.nz) then

        ! Error message
        urou = 'restrict_zaxis'
        write(umsg,'(A,i3,A,i3)') ' # The resulting PRD upper '// &
                                  'limit is above ',nz,' :',Rz1_PRD
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,.True.,.False.)

      end if ! Non-sense index

      ! If non-sense
      if (Rz1_PRD.gt.Rz1) then

        ! Error message
        urou = 'restrict_zaxis'
        write(umsg,'(A,i3,A,i3)') ' # The resulting PRD upper '// &
                                  'limit is above the restricted ', &
                                  Rz1,' :',Rz1_PRD
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,.True.,.False.)

      end if ! Non-sense index

      ! If non-sense
      if (Rz1_PRD.lt.Rz0) then

        ! Error message
        urou = 'restrict_zaxis'
        write(umsg,'(A,i3,A,i3)') ' # The resulting PRD upper '// &
                                  'limit is below the beginning ', &
                                  Rz0,' :',Rz1_PRD
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,.True.,.False.)

      end if ! Non-sense index

      return

      end subroutine restrict_zaxis

!#####################################################################
!#####################################################################
!#####################################################################

      !> Find the index to limit the height axis in the RT problem for
      !! a LTE line
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!  lines(LTEline_class(:)): Structures with LTE line data
      subroutine restrict_LTE_lines(Atmo,lines)

      ! I/O
      type(LTEline_class), dimension(:), intent(inout):: lines
      type(Atmo_class), intent(in):: Atmo

      ! Local

      integer:: ierr,ia,iz,tz0,zz0,mTiz

      double precision:: r0


      ! For each LTE line
      do ia=1,nLTEl

        ! Initialize to current limit
        lines(ia)%Rz0 = Rz0

        ! If not limiting, skip the determination
        if (.not.lines(ia)%taulim_l.and..not.lines(ia)%Tlim_l) cycle

        ! If Master
        if (pid.eq.0) then

          ! If limiting tau
          if (lines(ia)%taulim_l) then

            ! Initialize index
            tz0 = -1

            ! Shorten input
            r0 = lines(ia)%taulim

            ! If tau is in main axis
            if (ztau) then

              ! If already larger than the limit, take it
              if (Atmo%z(Rz0).ge.r0) tz0 = Rz0

              ! For the rest of heights going down
              do iz=Rz0+1,Rz1

                ! If limit between this and previous node
                if (Atmo%z(iz-1).lt.r0.and.Atmo%z(iz).ge.r0) then

                  ! Cut at previous index and leave
                  tz0 = iz - 1
                  exit

                end if ! Limit between this and previous node

              end do ! Heights

            ! If tau is in secondary axis
            else

              ! If already larger than the limit, take it
              if (Atmo%zalt(Rz0).ge.r0) tz0 = Rz0

              ! For the rest of heights going down
              do iz=Rz0+1,Rz1

                ! If limit between this and previous node
                if (Atmo%zalt(iz-1).lt.r0.and. &
                    Atmo%zalt(iz).ge.r0) then

                  ! Cut at previous index and leave
                  tz0 = iz - 1

                end if ! Limit between this and previous node

              end do ! Heights

            end if ! Tau in main or secondary axis

            ! If not restricting in temperature, done
            if (.not.lines(ia)%Tlim_l) lines(ia)%Rz0 = tz0

          end if ! If restricting in optical depth

          ! If restricting in temperature
          if (lines(ia)%Tlim_l) then

            ! Initialize index
            zz0 = -1

            ! Get location of minimum of temperature
            mTiz = minloc(Atmo%T(Rz0:Rz1),1) + Rz0 - 1

            ! Go up from the location of the minimum of temperature
            ! til the top
            do iz=mTiz,Rz0,-1

              ! If the temperature is larger than the limit
              if (Atmo%T(iz).gt.lines(ia)%Tlim) then

                ! Save index and leave
                zz0 = iz
                exit

              end if ! Temperature larger than the limit

            end do ! Heights

            ! If not found
            if (zz0.lt.0) zz0 = 1

            ! If not restricting in optical depth, done
            if (.not.lines(ia)%taulim_l) lines(ia)%Rz0 = zz0

          end if ! If restricting in temperature

          ! If restricting both tau and temperature, take the more
          ! conservative
          if (lines(ia)%taulim_l.and.lines(ia)%Tlim_l) &
            lines(ia)%Rz0 = max(zz0,tz0)

        end if ! Master

        ! If MPI, broadcast the limit
        if (nproc.gt.1) &
          call MPI_BCAST(lines(ia)%Rz0,1,MPI_INTEGER,0, &
                         MPI_COMM_RT,ierr)

        ! If non-sense
        if (lines(ia)%Rz0.lt.Rz1) then

          ! Error message
          urou = 'restrict_LTE_lines'
          write(umsg,'(3(A,i3))') &
            ' # The resulting lower for LTE line ',ia, &
            'is the global depth limit ',lines(ia)%Rz0,' > ',Rz1
          if (pid.eq.0) call verbose

        end if ! Non-sensical index

      end do ! LTE lines

      return

      end subroutine restrict_LTE_lines

!#####################################################################
!#####################################################################
!#####################################################################

      end module getztau_mod
