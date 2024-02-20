      !> Background opacities
      module getztau_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     12/03/2019
!  Last version:
!     02/16/2024 V3.0.8
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     02/16/2024:    V3.0.8 - Added additional mode to restrict nodes
!                             in RT (TdPA)
!                           - Bugfix: wrong initialization of z
!                             indexes when restricting height and
!                             the upper limit does not apply (TdPA)
!
!     10/16/2023:    V3.0.7 - Ensure Atmo%zalt is allocated with the
!                             correct size in getztau (TdPA)
!
!     08/07/2023:    V3.0.6 - Added exit condition on the lower
!                             index in restrict_zaxis (TdPA)
!                           - Added restrict_LTE_lines (TdPA)
!
!     07/03/2023:    V3.0.5 - In getztau, added argument to avoid
!                             slaves releasing data needed in the
!                             inversion (TdPA)
!                           - chi500 is now part of the atmospheric
!                             model, so the continuum structure is
!                             not needed anymore (TdPA)
!
!     02/10/2022:    V3.0.4 - Bugfix: It is necessary to add a check
!                             on the lower limit in restrict_zaxis
!                             because the upper limit could be
!                             predetermined and the algorithm then
!                             skips at the first step without
!                             finding the lower index (TdPA)
!                           - Added a non-sense check for the output
!                             of restrict_zaxis (TdPA)
!
!     10/27/2022:    V3.0.3 - Removed debugging messages (TdPA)
!
!     10/25/2022:    V3.0.2 - Added restrict_zaxis routine (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case changed the
!                             MPI communicator from MPI_COMM_WORLD
!                             to MPI_COMM_RT (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     05/11/2020:    V1.0.1 - Bugfix: Wrote continue instead of cycle
!                             to advance loop when domain there is
!                             domain decomposition, making the master
!                             listen for the wrong process (TdPA)
!
!     12/03/2019:    V1.0.0 - Started coding (TdPA)
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
!  getztau:
!    Computes height or tau scale given the other
!
!  restrict_zaxis:
!    Find the indexes to limit the height axis in RT from the
!  specified input in heights and/or continuum optical depth
!
!  restrict_LTE_lines:
!    Find the index to limit the height axis in RT for a LTE line
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

      !> Computes height or tau scale given the other.\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        can_leave(logical): Indicate if slaves can leave
      subroutine getztau(Atmo,MPID,can_leave)

      ! I/O
      type(Atmo_class), intent(inout):: Atmo
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: can_leave

      ! Local

      integer:: iz,iz0
      double precision:: dz,dt,mint,chi

      ! If MPI, master needs to receive chi500
      if (MPID%mpi) then

        ! Master
        if (pid.eq.0) then

          ! Allocate chi500
          if (.not.allocated(Atmo%chi500)) &
            allocate(Atmo%chi500(nZ))

          call MPI_RECV(Atmo%chi500,nZ,MPI_DOUBLE_PRECISION, &
                        1,1,MPI_COMM_RT,MPI_STATUS_IGNORE,ierr)

        ! Slaves
        else

          ! Is enough with the first process
          if (pid.eq.1) then

            ! Send the chi500 data
            call MPI_SEND(Atmo%chi500,nZ,MPI_DOUBLE_PRECISION, &
                          0,1,MPI_COMM_RT,ierr)

          end if

          ! If no tau scale, chi500 is not needed
          if (.not.ztau.and.can_leave) deallocate(Atmo%chi500)

        end if ! Master or slave
      end if ! MPI

      ! Rest of CPU can leave now
      if (pid.gt.0.and.can_leave) return

      ! Allocate zalt scale
      if (allocated(Atmo%zalt)) then
        if (nZ.ne.size(Atmo%zalt)) then
          deallocate(Atmo%zalt)
          allocate(Atmo%zalt(nZ))
        end if
      else
        allocate(Atmo%zalt(nZ))
      end if

      ! If tau was given
      if (ztau) then

        ! Initialize mint
        mint = 1d99
        iz0 = -1

        ! Look for node closest to tau=1
        do iz=1,nZ

          ! Check tau closer to 1
          if (abs(Atmo%z(iz)-1d0).lt.abs(mint-1d0)) then
            iz0 = iz
            mint = Atmo%z(iz)
          end if

          ! Abort if larger than one
          if (Atmo%z(iz).gt.1d0) exit

        end do ! Height nodes

        ! Initialize bottom of the atmosphere
        Atmo%zalt(nZ) = 0d0

        ! For each height
        do iz=nZ-1,1,-1

          ! Optical path
          dt = Atmo%z(iz+1) - Atmo%z(iz)

          ! Opacity
          chi = Atmo%chi500(iz+1) + Atmo%chi500(iz)

          ! 10^-5 cm -> km
          Atmo%zalt(iz) = Atmo%zalt(iz+1) + 2d-5*dt/chi

        end do ! Height nodes

        ! Shift height to tau=1
        Atmo%zalt = Atmo%zalt - Atmo%zalt(iz0)

      ! If z was given
      else

        ! Initialize top at 0
        Atmo%zalt(1) = 0d0

        ! For every height, add optical depth at reference freq
        do iz=2,NZ

          ! Height step
          dz = Atmo%z(iz-1) - Atmo%z(iz)

          ! Optical depth step
          dt = 0.5*dz*(Atmo%chi500(iz) + Atmo%chi500(iz-1))

          ! Add to optical depth
          Atmo%zalt(iz) = Atmo%zalt(iz-1) + dt

        end do ! Height nodes

      end if ! Given tau or height

      return

      end subroutine getztau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Restrict the vertical domain to compute in from a height or
      !! optical depth constrain
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Input(Input_class): Structure with settings data\n
      !!      MPID(MPI_class): Structure with MPI data
      subroutine restrict_zaxis(Atmo,Input,MPID)

      ! I/O
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID

      ! Local

      integer:: ierr
      integer:: iz,tz0,tz1,zz0,zz1
      double precision:: r0,r1

      ! If not restricting
      if (.not.(Input%rest_tau.or.Input%rest_z)) then
        Rz0 = 1
        Rz1 = nz
        Rnz = nz
        return
      end if

      ! If Master
      if (pid.eq.0) then

        ! Restrict tau
        if (Input%rest_tau) then

          ! Initialize
          tz0 = -1
          tz1 = -1

          ! Shortener
          r0 = Input%r0tc
          r1 = Input%r1tc

          ! Tau is main
          if (ztau) then

            ! If already larger
            if (Atmo%z(1).ge.r0) tz0 = 1
            ! If already smaller
            if (Atmo%z(nz).le.r1) tz1 = nz

            ! Find
            do iz=2,nz-1
              if (Atmo%z(iz-1).lt.r0.and.Atmo%z(iz).ge.r0) then
                if (Input%rest_tau_strc) then
                  tz0 = iz
                else
                  tz0 = iz - 1
                end if
              end if
              if (Atmo%z(iz).le.r1.and.Atmo%z(iz+1).gt.r1) then
                if (Input%rest_tau_strc) then
                  tz1 = iz
                else
                  tz1 = iz + 1
                end if
              end if
              if (tz1.gt.0.and.tz0.gt.0) exit
            end do

          ! Tau is secondary
          else

            ! If already larger
            if (Atmo%zalt(1).ge.r0) tz0 = 1
            ! If already smaller
            if (Atmo%zalt(nz).le.r1) tz1 = nz

            ! Find
            do iz=2,nz-1
              if (Atmo%zalt(iz-1).lt.r0.and.Atmo%zalt(iz).ge.r0) then
                if (Input%rest_tau_strc) then
                  tz0 = iz
                else
                  tz0 = iz - 1
                end if
              end if
              if (Atmo%zalt(iz).le.r1.and.Atmo%zalt(iz+1).gt.r1) then
                if (Input%rest_tau_strc) then
                  tz1 = iz
                else
                  tz1 = iz + 1
                end if
              end if
              if (tz1.gt.0.and.tz0.gt.0) exit
            end do

          end if ! Tau main

          ! If not z, done
          if (.not.Input%rest_z) then
            Rz0 = tz0
            Rz1 = tz1
          end if ! Not restrict in z
        end if ! Restrict tau

        ! Restrict Height
        if (Input%rest_z) then

          ! Initialize
          zz0 = -1
          zz1 = -1

          ! Shortener
          r0 = Input%r0z*1d5
          r1 = Input%r1z*1d5

          ! Tau is main
          if (ztau) then

            ! If already larger
            if (Atmo%zalt(nz).ge.r0) zz1 = nz
            ! If already smaller
            if (Atmo%zalt(1).le.r1) zz0 = 1

            ! Find
            do iz=2,nz-1
              if (Atmo%zalt(iz-1).gt.r1.and.Atmo%zalt(iz).le.r1) then
                if (Input%rest_z_strc) then
                  zz0 = iz
                else
                  zz0 = iz - 1
                end if
              end if
              if (Atmo%zalt(iz).ge.r0.and.Atmo%zalt(iz+1).lt.r0) then
                if (Input%rest_z_strc) then
                  zz1 = iz
                else
                  zz1 = iz + 1
                end if
              end if
              if (zz1.gt.0.and.zz0.gt.0) exit
            end do

          ! Tau is secondary
          else

            ! If already larger
            if (Atmo%z(nz).ge.r0) zz1 = nz
            ! If already smaller
            if (Atmo%z(1).le.r1) zz0 = 1

            ! Find
            do iz=2,nz-1
              if (Atmo%z(iz-1).gt.r1.and.Atmo%z(iz).le.r1) then
                if (Input%rest_z_strc) then
                  zz0 = iz
                else
                  zz0 = iz - 1
                end if
              end if
              if (Atmo%z(iz).ge.r0.and.Atmo%z(iz+1).lt.r0) then
                if (Input%rest_z_strc) then
                  zz1 = iz
                else
                  zz1 = iz + 1
                end if
              end if
              if (zz1.gt.0.and.zz0.gt.0) exit
            end do

          end if ! Tau main

          ! If not tau, done
          if (.not.Input%rest_tau) then
            Rz0 = zz0
            Rz1 = zz1
          end if

        end if ! Restrict Height

        ! Restrict both tau and height
        if (Input%rest_z.and.Input%rest_tau) then
          Rz0 = max(zz0,tz0)
          Rz1 = min(zz1,tz1)
        end if
      end if ! Master

      ! If MPI
      if (MPID%mpi) then

        ! Broadcast limits
        call MPI_BCAST(Rz0,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)
        call MPI_BCAST(Rz1,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)

      end if

      ! Dimension
      Rnz = Rz1 - Rz0 + 1

      ! If smaller than 3, problem
      if (Rnz.lt.3) then

        urou = 'restrict_zaxis'
        umsg = ' # The resulting atmosphere after the height '// &
               'restriction has less than three nodes'
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,0,.True.,.False.)

      end if

      ! If non-sense
      if (Rz0.lt.1) then

        urou = 'restrict_zaxis'
        write(umsg,'(A,i3)') ' # The resulting lower limit is '// &
                             'below 1 :',Rz0
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,0,.True.,.False.)

      end if

      ! If non-sense
      if (Rz1.gt.nz) then

        urou = 'restrict_zaxis'
        write(umsg,'(A,i3,A,i3)') ' # The resulting upper '// &
                                  'limit is above ',nz,' :',Rz1
        if (pid.eq.0) call verbose
        call abortedS(umsg,urou,0,.True.,.False.)

      end if

      return

      end subroutine restrict_zaxis

!#####################################################################
!#####################################################################
!#####################################################################

      !> Restrict the region of validity of LTE lines depending on
      !! its input\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !! lines(LTEline_class): Structure with the LTE line data\n
      !!      MPID(MPI_class): Structure with MPI data
      subroutine restrict_LTE_lines(Atmo,lines,MPID)

      ! I/O
      type(LTEline_class), dimension(:):: lines
      type(Atmo_class), intent(in):: Atmo
      type(MPI_class), intent(in):: MPID

      ! Local

      integer:: ierr,ia
      integer:: iz,tz0,zz0,mTiz
      double precision:: r0

      ! For each LTE line
      do ia=1,nLTEl

        ! Initialize
        lines(ia)%Rz0 = Rz0

        ! If no limit, skip
        if (.not.lines(ia)%taulim_l.and..not.lines(ia)%Tlim_l) cycle

        ! If Master
        if (pid.eq.0) then

          ! taulim
          if (lines(ia)%taulim_l) then

            ! Initialize
            tz0 = -1

            ! Shortener
            r0 = lines(ia)%taulim

            ! Tau is main
            if (ztau) then

              ! If already larger
              if (Atmo%z(Rz0).ge.r0) tz0 = Rz0

              ! Find
              do iz=Rz0+1,Rz1
                if (Atmo%z(iz-1).lt.r0.and.Atmo%z(iz).ge.r0) then
                  tz0 = iz - 1
                  exit
                end if
              end do

            ! Tau is secondary
            else

              ! If already larger
              if (Atmo%zalt(Rz0).ge.r0) tz0 = Rz0

              ! Find
              do iz=Rz0+1,Rz1
                if (Atmo%zalt(iz-1).lt.r0.and. &
                    Atmo%zalt(iz).ge.r0) then
                  tz0 = iz - 1
                  exit
                end if
              end do

            end if ! Tau main

            ! If not T, done
            if (.not.lines(ia)%Tlim_l) lines(ia)%Rz0 = tz0

          end if ! Restrict tau

          ! Restrict Temperature
          if (lines(ia)%Tlim_l) then

            ! Initialize
            zz0 = -1

            ! Get location minimum of T
            mTiz = minloc(Atmo%T(Rz0:Rz1),1) + Rz0 - 1

            ! Find
            do iz=mTiz,Rz0,-1
              if (Atmo%T(iz).gt.lines(ia)%Tlim) then
                zz0 = iz
                exit
              end if
            end do

            ! If not found
            if (zz0.lt.0) zz0 = 1

            ! If not tau, done
            if (.not.lines(ia)%taulim_l) lines(ia)%Rz0 = zz0

          end if ! Restrict Height

          ! Restrict both tau and height
          if (lines(ia)%taulim_l.and.lines(ia)%Tlim_l) &
            lines(ia)%Rz0 = max(zz0,tz0)

        end if ! Master

        ! If MPI, broadcast the limit
        if (MPID%mpi) &
          call MPI_BCAST(lines(ia)%Rz0,1,MPI_INTEGER,0, &
                         MPI_COMM_RT,ierr)

        ! If non-sense
        if (lines(ia)%Rz0.lt.Rz1) then

          urou = 'restrict_LTE_lines'
          write(umsg,'(3(A,i3))') &
            ' # The resulting lower for LTE line ',ia, &
            'is the global depth limit ',lines(ia)%Rz0,' > ',Rz1
          if (pid.eq.0) call verbose

        end if

      end do ! LTE lines

      return

      end subroutine restrict_LTE_lines

!#####################################################################
!#####################################################################
!#####################################################################

      end module getztau_mod
