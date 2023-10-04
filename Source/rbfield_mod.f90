      !> Reading magnetic data
      module rbfield_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/17/2017
!  Last version:
!     09/21/2023 V3.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/21/2023:    V3.0.3 - Bugfix: Need to distinguish the run mode
!                             to choose the correct verbose routine
!                             to call (TdPA)
!                           - For safety reasons, rbfield must be
!                             called by all processes in
!                             MPI_COMM_WORLD (TdPA)
!
!     08/11/2023:    V3.0.2 - Bugfix: verbosity for reading directed
!                             to the wrong file and not limited to the
!                             global master (TdPA)
!
!     03/15/2023:    V3.0.1 - Added the atmosphere size as argument
!                             for the inversion branch (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     01/13/2021:    V1.1.3 - Removed capitalization (TdPA)
!
!     11/19/2019:    V1.1.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     06/12/2019:    V1.1.1 - Fixed constant magnetic field numeric
!                             input (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!                           - Checks for success of python routine
!                             and unit is now 100 (TdPA)
!                           - Now it can interpret numerical values
!                             in the input for homogeneous fields
!                             (TdPA)
!
!     08/08/2018:    V1.0.2 - Abort if the nodes in the magnetic field
!                             file do not fit with the atmospheric
!                             model (TdPA)
!
!     09/14/2017:    V1.0.1 - Added a path and ID to the file (TdPA)
!
!     04/17/2016:    V1.0.0 - First version (TdPA)
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
!    This subroutine reads the magnetic field data
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

      !> Reads a file with magnetic field data.\n
      !!  filename(character(:)): Name of the file to read\n
      !!    source(character(:)): Path to the source code\n
      !!        ID(character(:)): ID of this run\n
      !!    Bfield(Bfield_class): Structure with magnetic field data\n
      !!            lnz(integer): Expected height nodes for the file\n
      !!      Input(Input_class): Structure with input data
      subroutine rBField(filename,source,ID,Bfield,lnz,Input)

      ! I/O

      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(in):: Input
      integer, intent(in):: lnz

      ! Local

      character(len=65):: formt
      character(len=3):: cdump
      character(len=6):: formtBs
      character(len=4):: formtBt,formtBp

      integer:: iz, nZ_bfield, ios

      double precision:: B_l(lnZ), T_l(lnZ), F_l(lnZ)


      ! If numeric values
      if (Input%bfieldn) then

        ! Store in the structure
        allocate(Bfield%Bstrength(lnZ))
        allocate(Bfield%Btheta(lnZ))
        allocate(Bfield%Bphi(lnZ))

        ! If no field
        if (Input%bfieldv(1).lt.0d0) then

          Bfield%Bstrength = 0d0
          Bfield%Btheta = 0d0
          Bfield%Bphi = 0d0

          ! If master
          if (gpid.eq.0) then
            umsg = ' - No magnetic field'
            if (run_mode.eq.-1) then
              call verbosev
            else
              call verbose
            end if
          end if

        ! If yes field
        else

          Bfield%Bstrength = Input%bfieldv(1)
          Bfield%Btheta = Input%bfieldv(2)
          Bfield%Bphi = Input%bfieldv(3)

          ! Adjust angles
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
            if (Bfield%Btheta(1).lt.10d0) then
              formtBt = 'f3.1'
            else if (Bfield%Btheta(1).lt.100d0) then
              formtBt = 'f4.1'
            else
              formtBt = 'f5.1'
            end if
            if (Bfield%Bphi(1).lt.10d0) then
              formtBp = 'f3.1'
            else if (Bfield%Bphi(1).lt.100d0) then
              formtBp = 'f4.1'
            else
              formtBp = 'f5.1'
            end if
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
            end if
          end if

          Bfield%Btheta = Bfield%Btheta/RAD
          Bfield%Bphi = Bfield%Bphi/RAD

        end if

        return

      end if ! Numeric field


      ! If no file, no magnetic field
      if(trim(filename).eq.'NONE')then

        ! No Magnetic field
        do iz=1,lnZ
          B_l(iz) = 0d0
          T_l(iz) = 0d0
          F_l(iz) = 0d0
        end do

        ! If Master
        if (gpid.eq.0) then
          umsg = ' - No magnetic field'
          if (run_mode.eq.-1) then
            call verbosev
          else
            call verbose
          end if
        end if

      ! If there is a file, read the magnetic field
      else

        ! Yes magnetic field
        if(gpid.eq.0) call system('python '//trim(source)// &
                                 'rbfield.py '//trim(filename)// &
                                 ' '//ID//' '//verbosef)
        ! Wait for the master to finish
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        open(100, file='tmp_bfield_'//ID, status='old', &
             iostat=ios,err=1000)

        ! Success
        read (100,*,err=1100) ios

        ! If no correct file, abort
        if (ios.lt.0) then

          umsg = 'Problem translating the magnetic file'
          goto 1200

        end if

        read(100,*,err=1100) nZ_bfield

        ! If height dependent
        if(nZ_bfield.eq.lnZ)then

          ! Correct dimensions
          read(100,*,err=1100) cdump

          do iz=1,lnZ

            read(100,*,err=1100) B_l(iz), T_l(iz), F_l(iz)

            if(trim(cdump).eq.'DEG')then

              T_l(iz) = T_l(iz)/RAD
              F_l(iz) = F_l(iz)/RAD

            end if

          end do

        ! If homogeneous
        else if(nZ_bfield.eq.1)then

          ! Simplified input
          read(100,*,err=1100) cdump
          read(100,*,err=1100) B_l(1), T_l(1), F_l(1)

          if(trim(cdump).eq.'DEG')then

            T_l(1) = T_l(1)/RAD
            F_l(1) = F_l(1)/RAD

          end if

          do iz=2,lnZ

            B_l(iz) = B_l(1)
            T_l(iz) = T_l(1)
            F_l(iz) = F_l(1)

          end do

        ! If none of the others
        else

          ! Wrong input
          umsg = 'The field provided does not have the'// &
                 ' same length than the atmosphere'
          goto 1200

        end if

        close(100)

        ! Control that everything went fine
        call control

        ! Delete temporal input file
        if (gpid.eq.0) then
          call system('rm tmp_bfield_'//ID)
          umsg = ' - Magnetic field '//trim(filename)//' read'
          if (run_mode.eq.-1) then
            call verbosev
          else
            call verbose
          end if
        end if

      end if

      ! Store in the structure
      allocate(Bfield%Bstrength(lnZ))
      allocate(Bfield%Btheta(lnZ))
      allocate(Bfield%Bphi(lnZ))

      do iz=1,lnZ

        Bfield%Bstrength(iz) = B_l(iz)

        ! Make sure that is no field no rotation if there is no B
        if(Bfield%Bstrength(iz).le.TINYB)then
          Bfield%Btheta(iz) = 0d0
          Bfield%Bphi(iz) = 0d0
        else
          Bfield%Btheta(iz) = T_l(iz)
          Bfield%Bphi(iz) = F_l(iz)
        end if

      enddo

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
