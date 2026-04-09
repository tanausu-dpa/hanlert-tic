      !> Reading magnetic data
      module rbfield_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     17/04/2017
!  Last version:
!     09/04/2026 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/04/2026:    V4.0.1 - Changed the call to python to python3
!                             as some systems do not have just python
!                             anymore (TdPA)
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
!  rBField
!    Read the magnetic field stratification from the specified file or
!  input
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : RAD, TINYB
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read the magnetic field stratification from the specified
      !! file or input\n
      !!  filename(character(:)): Name of the file to read\n
      !!    source(character(:)): Path to the source code folder\n
      !!        ID(character(:)): ID of this run\n
      !!    Bfield(Bfield_class): Structure with magnetic field data\n
      !!            lnz(integer): Expected height nodes for the file\n
      !!      Input(Input_class): Structure with configuration data
      subroutine rBField(filename,source,ID,Bfield,lnz,Input)

      ! I/O

      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(in):: Input
      integer, intent(in):: lnz

      ! Local

      character(len=3):: cdump
      character(len=6):: formtBs
      character(len=4):: formtBt,formtBp
      character(len=65):: formt

      integer:: iz,nZ_bfield,ios

      double precision:: B_l(lnZ),T_l(lnZ),F_l(lnZ)


      ! If numeric values
      if (Input%bfieldn) then

        ! Allocate arrays
        allocate(Bfield%Bstrength(lnZ))
        MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bstrength)
        allocate(Bfield%Btheta(lnZ))
        MRAMc = MRAMc + 1d-6*sizeof(Bfield%Btheta)
        allocate(Bfield%Bphi(lnZ))
        MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bphi)

        ! If no field
        if (Input%bfieldv(1).lt.0d0) then

          ! Set to zero
          Bfield%Bstrength = 0d0
          Bfield%Btheta = 0d0
          Bfield%Bphi = 0d0

          ! If master
          if (gpid.eq.0) then

            ! Verbose zero field
            umsg = ' - No magnetic field'
            if (run_mode.eq.-1) then
              call verbosev
            else
              call verbose
            end if ! Run mode
          end if ! Master

        ! If yes field
        else

          ! Save in array (homogeneous)
          Bfield%Bstrength = Input%bfieldv(1)
          Bfield%Btheta = Input%bfieldv(2)
          Bfield%Bphi = Input%bfieldv(3)

          ! Adjust angles (they are still in degrees)
          do while (Bfield%Btheta(1).lt.0d0)
            Bfield%Btheta = Bfield%Btheta + 360d0
          end do
          do while (Bfield%Btheta(1).gt.360d0)
            Bfield%Btheta = Bfield%Btheta - 360d0
          end do
          do while (Bfield%Bphi(1).lt.0d0)
            Bfield%Bphi = Bfield%Bphi + 360d0
          end do
          do while (Bfield%Bphi(1).gt.360d0)
            Bfield%Bphi = Bfield%Bphi - 360d0
          end do

          ! If master (global)
          if (gpid.eq.0) then

            ! Determine correct format for |B|
            if (Bfield%Bstrength(1).gt.9999.9d0.or. &
                Bfield%Bstrength(1).lt.0.1d0) then
              formtBs = 'es7.1'
            else
              if (Bfield%Bstrength(1).lt.10.0) then
                formtBs = 'f3.1'
              else if (Bfield%Bstrength(1).lt.100.0) then
                formtBs = 'f4.1'
              else if (Bfield%Bstrength(1).lt.1000.0) then
                formtBs = 'f5.1'
              else
                formtBs = 'f6.1'
              end if
            end if

            ! Determine correct format for \Theta_B
            if (Bfield%Btheta(1).lt.10d0) then
              formtBt = 'f3.1'
            else if (Bfield%Btheta(1).lt.100d0) then
              formtBt = 'f4.1'
            else
              formtBt = 'f5.1'
            end if

            ! Determine correct format for \Phi_B
            if (Bfield%Bphi(1).lt.10d0) then
              formtBp = 'f3.1'
            else if (Bfield%Bphi(1).lt.100d0) then
              formtBp = 'f4.1'
            else
              formtBp = 'f5.1'
            end if

            ! Verbose magnetic field constant value
            formt = "(' - Constant magnetic field: ',"// &
                    trim(formtBs)//",'G ',"// &
                    formtBt//",'º ',"//formtBp//",'º')"
            write(umsg,trim(formt)) Bfield%Bstrength(1), &
                                    Bfield%Btheta(1), &
                                    Bfield%Bphi(1)
            if (run_mode.eq.-1) then
              call verbosev
            else
              call verbose
            end if ! Type of run
          end if ! Global master

          ! Transform angles to radians
          Bfield%Btheta = Bfield%Btheta/RAD
          Bfield%Bphi = Bfield%Bphi/RAD

        end if ! If actual magnetic field value

        ! And leave
        return

      end if ! Numeric constant input for magnetic field


      ! If no file
      if(trim(filename).eq.'NONE')then

        ! No Magnetic field
        do iz=1,lnZ
          B_l(iz) = 0d0
          T_l(iz) = 0d0
          F_l(iz) = 0d0
        end do

        ! If global Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = ' - No magnetic field'
          if (run_mode.eq.-1) then
            call verbosev
          else
            call verbose
          end if ! Type of run
        end if ! Global Master

      ! If there is a file, read the magnetic field
      else

        ! Translate the magnetic field file in python
        if(gpid.eq.0) call system('python3 '//trim(source)// &
                                 'rbfield.py '//trim(filename)// &
                                 ' '//ID//' '//verbosef)

        ! Wait for the master to finish
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ! Open translated file
        open(100, file='tmp_bfield_'//ID, status='old', &
             iostat=ios,err=1000)

        ! Success
        read (100,*,err=1100) ios

        ! If no correct file
        if (ios.lt.0) then

          ! Issue error
          umsg = 'Problem translating the magnetic file'
          goto 1200

        end if ! Wrong file

        ! Read dimension in input model
        read(100,*,err=1100) nZ_bfield

        ! If height dependent and equal to model atmosphere
        if(nZ_bfield.eq.lnZ)then

          ! Read angular units
          read(100,*,err=1100) cdump

          ! For each height
          do iz=1,lnZ

            ! Read magnetic field in polar coordinates
            read(100,*,err=1100) B_l(iz), T_l(iz), F_l(iz)

          end do ! Heights

          ! If input is in degrees
          if(trim(cdump).eq.'DEG')then

            ! Transform into radians
            T_l = T_l/RAD
            F_l = F_l/RAD

          end if ! Input in degrees

        ! If homogeneous
        else if(nZ_bfield.eq.1)then

          ! Read angular units
          read(100,*,err=1100) cdump

          ! Read magnetic field in polar coordinates
          read(100,*,err=1100) B_l(1), T_l(1), F_l(1)

          ! If input is in degrees
          if(trim(cdump).eq.'DEG')then

            ! Transform into radians
            T_l(1) = T_l(1)/RAD
            F_l(1) = F_l(1)/RAD

          end if ! Input in degrees

          ! For the rest of heights
          do iz=2,lnZ

            ! Copy value
            B_l(iz) = B_l(1)
            T_l(iz) = T_l(1)
            F_l(iz) = F_l(1)

          end do ! Heights

        ! If none of the others
        else

          ! Wrong input, issue error
          umsg = 'The field provided does not have the'// &
                 ' same length than the atmosphere'
          goto 1200

        end if ! Height dimension of input

        ! Close file
        close(100)

        ! Control that everything went fine
        call control

        ! Global Master
        if (gpid.eq.0) then

          ! Delete temporal input file
          call system('rm tmp_bfield_'//ID)

          ! Verbose
          umsg = ' - Magnetic field '//trim(filename)//' read'
          if (run_mode.eq.-1) then
            call verbosev
          else
            call verbose
          end if ! Type of run
        end if ! Global Master
      end if ! Specified file or not

      ! Allocate structure arrays
      allocate(Bfield%Bstrength(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bstrength)
      allocate(Bfield%Btheta(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Btheta)
      allocate(Bfield%Bphi(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bphi)

      ! For each height
      do iz=1,lnZ

        ! Copy |B|
        Bfield%Bstrength(iz) = B_l(iz)

        ! If |B| is too small
        if(Bfield%Bstrength(iz).le.TINYB)then

          ! No angles
          Bfield%Btheta(iz) = 0d0
          Bfield%Bphi(iz) = 0d0

        ! Significant |B|
        else

          ! Store angles
          Bfield%Btheta(iz) = T_l(iz)
          Bfield%Bphi(iz) = F_l(iz)

        end if ! |B| value

      enddo ! Heights

      ! Control that everything went fine
      call control

      return

1000  umsg = ' - Error opening field file'
      call aborted
1100  umsg = ' - Error reading field file'
1200  close(100)
      call aborted

      end subroutine rBfield

!#####################################################################
!#####################################################################
!#####################################################################

      end module rbfield_mod
