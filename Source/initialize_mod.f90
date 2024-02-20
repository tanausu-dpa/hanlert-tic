      !> Initialization of radiation field tensors
      module initialize_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
!  Last version:
!     02/14/2024 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     02/14/2024:    V3.0.5 - RAM variable for radiation does not
!                             assume not initialized (TdPA)
!
!     08/07/2023:    V3.0.4 - Added initialize_failread and
!                             initializeI_failread (TdPA)
!
!     11/10/2022:    V3.0.3 - When doing PRD with velocities, the
!                             slaves also need JKQa to add it in
!                             emiss2ord/NB (TdPA)
!
!     10/25/2022:    V3.0.2 - Implemented the height range limitation
!                             in the allocation of the radiation
!                             variables (TdPA)
!                           - Implemented the 1.5D case for the ad-hoc
!                             radiation field tensors (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Allocate dummy asymmetry JKQ if not
!                             input (TdPA)
!
!     01/13/2021:    V1.2.6 - Added initialize_asym subroutine (TdPA)
!
!     01/12/2021:    V1.2.5 - Use KSTK to determine how to initialize
!                             Stokes (TdPA)
!
!     09/11/2020:    V1.2.4 - Now the RAM used by the radiation
!                             variables is counted and kept (TdPA)
!
!     03/18/2020:    V1.2.3 - Initialized radiation quantities when
!                             allocating (TdPA)
!
!     11/19/2019:    V1.2.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     09/13/2019:    V1.2.1 - Added condition to store only 2 heights
!                             of Stokes to only the static case when
!                             there is angle-averaged PRD (TdPA)
!
!     05/08/2019:    V1.2.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     05/16/2018:    V1.0.3 - Stokes has dimension nPh (and not nPh2)
!                             for the azimuth. If there is axial
!                             symmetry, that dimension is not
!                             necessary (TdPA)
!
!     06/23/2017:    V1.0.2 - Only initializing JKQC (TdPA)
!                           - Cleaned unused variables (TdPA)
!
!     06/22/2017:    V1.0.1 - Removed initializePtoI and
!                             initializeItoP (TdPA)
!                           - Not initializing to zero (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!    This subroutine initializes the radiation field tensors with a
!  Planckian function and allocates memory for everything
!
!  initialize
!    Initialize Stokes(polarization) and JKQ(terms)
!
!  initializeI
!    Initialize Stokes(intensity) and J00(levels)
!
!  initialize_failread
!    Initialize Stokes(polarization) and JKQC when cannot be read from
!  solution file
!
!  initializeI_failread
!    Initialize Stokes(intensity) and J00C when cannot be read from
!  solution file
!
!  initialize_asym
!    Initialize factor for ad-hoc asymmetries in JKQ
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : cZero
      use planck_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes Stokes parameters and radiation
      !! field tensors (initializes when not reading a previous
      !! solution).\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!           mode(character): Mode of operation (solving,
      !!                            reading or both)
      subroutine initialize(Frec,Geom,Atmo,MPID,Stokes,JKQ, &
                            JKQS,JKQC,J00P,mode)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      type(MPI_class), intent(inout):: MPID
      character(len=1), intent(in)::  mode
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:,:), allocatable:: J00P
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      integer:: iz,ifreq

      double precision:: omg


      !
      ! Allocations
      !

      ! If we are doing angle averaged, we only need one height,
      ! we allocate two to store the emergence in the quadrature
      if (KSTK) then
        allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
        giz0 = Rz0
        giz1 = Rz1
        MPID%RRAM = MPID%RRAM + &
                    8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*Rnz)
      else
        allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
        giz0 = Rz0
        giz1 = Rz0+1
        MPID%RRAM = MPID%RRAM + &
                    8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*2)
      end if
      Stokes = 0d0

      ! JKQ for absorptivity
      allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
      JKQ = cZero

      ! JKQ for stimulated emission
      allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))
      JKQS = cZero

      ! JKQ frequency dependent
      allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))

      ! J00 for photoionizations
      allocate(J00P(nxphot,2,Rz0:Rz1))
      J00P = 0d0

      ! Compute allocated memory
      MPID%RRAM = MPID%RRAM + 8d-6*dble(Rnz*(nxphot*2 + &
                                        2*5*3*(2*nxtran + nfreq)))

      ! If we are reading a solution, do not initialize to LTE
      if (mode.eq.'R'.or.mode.eq.'B') return

      !
      ! Make J00 frequency dependent planckian
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each frequency
        do ifreq=1,nfreq

          ! Current frequency
          omg = Frec%omega(ifreq)

          ! Initialize with planckian
          JKQC(0,0,ifreq,iz) = dcmplx(planck(omg,Atmo%T(iz)),0d0)

        end do ! frequencies
      end do ! heights

      ! Control
      call control

      return

      end subroutine initialize

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes intensity and radiation mean
      !! intensity (initializes when not reading a previous
      !! solution)\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!           mode(character): Mode of operation (solving,
      !!                            reading or both)
      subroutine initializeI(Frec,Geom,Atmo,MPID,Stokes,J00, &
                            J00S,J00C,J00P,mode)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      type(MPI_class), intent(inout):: MPID
      character(len=1), intent(in)::  mode
      double precision, dimension(:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P

      ! Local

      integer:: iz,ifreq

      double precision:: omg


      !
      ! Allocations
      !

      ! If we are doing angle averaged, we only need one height,
      ! we allocate two to store the emergence in the quadrature
      if (KSTK) then
        allocate(Stokes(nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
        giz0 = Rz0
        giz1 = Rz1
        MPID%RRAM = MPID%RRAM + &
                    8d-6*dble(nfreq*Geom%nph*Geom%nTh*Rnz)
      else
        allocate(Stokes(nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
        giz0 = Rz0
        giz1 = Rz0+1
        MPID%RRAM = MPID%RRAM + &
                    8d-6*dble(nfreq*Geom%nph*Geom%nTh*2)
      end if
      Stokes = 0d0

      ! J00 for absorptivity
      allocate(J00(nxt,Rz0:Rz1))
      J00 = 0d0

      ! J00 for stimulated emission
      allocate(J00S(nxt,Rz0:Rz1))
      J00S = 0d0

      ! J00 frequency dependent
      allocate(J00C(nfreq,Rz0:Rz1))

      ! If we are reading a solution, do not initialize to LTE
      if (mode.eq.'R'.or.mode.eq.'B') return

      ! J00 for photoionizations
      allocate(J00P(nxphot,2,Rz0:Rz1))
      J00P = 0d0

      ! Compute allocated memory
      MPID%RRAM = MPID%RRAM + 8d-6*dble(Rnz*(nxphot*2 + &
                                        2*nxt + nfreq))


      !
      ! Make J00 frequency dependent planckian
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each frequency
        do ifreq=1,nfreq

          ! Current frequency
          omg = Frec%omega(ifreq)

          ! Initialize with planckian
          J00C(ifreq,iz) = planck(omg,Atmo%T(iz))

        end do ! Frequencies
      end do ! heights

      ! Control
      call control

      return

      end subroutine initializeI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initializes JKQC for K=Q=0 when failing to read file\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      subroutine initialize_failread(Frec,Atmo,Stokes,JKQC)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      integer:: iz,ifreq

      double precision:: omg


      !
      ! Make J00 frequency dependent planckian
      !

      ! Zero out multipoles
      JKQC(:,1:2,:,:) = cZero

      ! For each height
      do iz=Rz0,Rz1

        ! For each frequency
        do ifreq=1,nfreq

          ! Current frequency
          omg = Frec%omega(ifreq)

          ! Initialize with planckian
          JKQC(0,0,ifreq,iz) = dcmplx(planck(omg,Atmo%T(iz)),0d0)

        end do ! frequencies
      end do ! heights

      ! Need Stokes
      if (KSTK) then

        ! For each height
        do iz=Rz0,Rz1

          ! For each frequency
          do ifreq=1,nfreq

            ! Initialize with planckian
            Stokes(0,ifreq,:,:,iz) = dble(JKQC(0,0,ifreq,iz))

          end do ! frequencies
        end do ! heights

      end if

      ! Control
      call control

      return

      end subroutine initialize_failread

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes intensity and radiation mean
      !! intensity (initializes when not reading a previous
      !! solution)\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence
      subroutine initializeI_failread(Frec,Atmo,Stokes,J00C)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      double precision, dimension(:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:), allocatable:: J00C

      ! Local

      integer:: iz,ifreq

      double precision:: omg


      !
      ! Make J00 frequency dependent planckian
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each frequency
        do ifreq=1,nfreq

          ! Current frequency
          omg = Frec%omega(ifreq)

          ! Initialize with planckian
          J00C(ifreq,iz) = planck(omg,Atmo%T(iz))

        end do ! Frequencies
      end do ! heights

      ! Need Stokes
      if (KSTK) then

        ! For each height
        do iz=Rz0,Rz1

          ! For each frequency
          do ifreq=1,nfreq

            ! Initialize with planckian
            Stokes(ifreq,:,:,iz) = dble(J00C(ifreq,iz))

          end do ! frequencies
        end do ! heights

      end if

      ! Control
      call control

      return

      end subroutine initializeI_failread

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes anisotropy inputs.\n
      !!      Input(Input_class): Structure with input data
      !!         MPID(MPI_class): Structure with MPI data\n
      !!      Flgsg(Fctsg_class): Structure with factorials and
      !!                          signs\n
      !!        JKQin(double(:)): Radiation field tensors extra
      !!                          asymmetries, 1.5D input\n
      !!   JKQa(dcomplex(:,:,:)): Radiation field tensors extra
      !!                          asymmetries\n
      subroutine initialize_asym(Input,MPID,Flgsg,JKQin,JKQa)

      ! I/O

      type(Input_class), intent(inout):: Input
      type(MPI_class), intent(inout):: MPID
      type(Fctsg_class):: Flgsg
      double precision, dimension(:), allocatable:: JKQin
      complex(kind=8), dimension(:,:,:), allocatable:: JKQa

      ! Local

      logical, dimension(0:2,1:2):: setable
      integer:: iz,ientry,K,iQ,ifentry,nfentry,fnz,ios
      double precision:: daux1,daux2


      !
      ! Only the master needs this (although leave if no
      ! inputs)
      !
      if (pid.gt.0.or.Input%nasym.le.0) then

        ! If there is input, it is slave, and it is dynamic
        if (Input%nasym.ge.1.and.dyn.and.PRD) then

          !
          ! Allocate asymmetry input and add to RAM
          !
          allocate(JKQa(-2:2,1:2,Rz0:Rz1))
          MPID%RRAM = MPID%RRAM + 16d-6*dble(nz)*10d0

          ! Size
          iz = 5*2*Rnz

          ! Share
          call MPI_BCAST(JKQa(-2,1,Rz0),iz,MPI_DOUBLE_COMPLEX,0, &
                         MPI_COMM_RT,ios)

        ! There is no input or it is static
        else

          ! Allocate dummy
          allocate(JKQa(1,1,1))

        end if ! Need to share

        ! Control and leave
        call control
        return

      end if ! There is something to read

      !
      ! Allocate asymmetry input and add to RAM
      !
      allocate(JKQa(-2:2,1:2,Rz0:Rz1))
      MPID%RRAM = MPID%RRAM + 16d-6*dble(nz)*10d0

      !
      ! Initialize to zero
      !
      JKQa = cZero

      !
      ! 1.5D case (input buffer)
      !
      if (allocated(JKQin)) then

        ! Convert to JKQa array
        do iz=Rz0,Rz1
          JKQa(0,1,iz) = dcmplx(JKQin(iz),0d0)
          JKQa(1,1,iz) = dcmplx(JKQin(nz+iz),JKQin(2*nz+iz))
          JKQa(0,2,iz) = dcmplx(JKQin(3*nz+iz),0d0)
          JKQa(1,2,iz) = dcmplx(JKQin(4*nz+iz),JKQin(5*nz+iz))
          JKQa(2,2,iz) = dcmplx(JKQin(6*nz+iz),JKQin(7*nz+iz))
        end do

        ! Get negatives
        JKQa(-1,1,:) = -conjg(JKQa(1,1,:))
        JKQa(-1,2,:) = -conjg(JKQa(1,2,:))
        JKQa(-2,2,:) =  conjg(JKQa(2,2,:))

      !
      ! 1D case
      !
      else

        !
        ! Initialize set
        !
        setable = .True.

        !
        ! For each numberic entry
        !
        do ientry=1,Input%nasym_num

          ! Get K, Q
          K = nint(dble(Input%asym_num(1,ientry)))
          iQ = nint(dimag(Input%asym_num(1,ientry)))


          ! Check setable
          if (setable(iQ,K)) then
            setable(iQ,K) = .False.
          else
            write(umsg,'(A,1x,i1,1x,i1)') &
              'You have specified more than one anysotropy '// &
              'factor for the multipole',K,iQ
            urou = 'initialize_asym'
            call abortedS(umsg,urou,-1,.True.,.True.)
          end if

          !
          ! Set factor

          ! Get other Q
          if (iQ.gt.0) then

            JKQa(iQ,K,:) = Input%asym_num(2,ientry)
            JKQa(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JKQa(iQ,K,:))

          ! Q = 0, ensure real
          else

            JKQa(0,K,:) = dcmplx(dble(Input%asym_num(2,ientry)), 0d0)

          end if ! Q > 0 or Q = 0

        end do ! Numerical entries

        !
        ! For each file
        !
        do ientry=1,Input%nasym_fil

          ! Deal with the file in python
           call system('python '//trim(Input%source)//'rasym.py '// &
                        trim(Input%asym_fil(ientry)%str)// &
                       ' '//Input%ID//' '//verbosef)

          ! Open the file
           open(100, file='tmp_asym_'//Input%ID, status='old', &
                iostat=ios, err=1000)

          ! Success
          read (100,*,err=1100) ios

          ! If no correct file, abort
          if (ios.lt.0) then

            umsg = 'Problem translating the asymmetry file'
            urou = 'initialize_asym'
            call abortedS(umsg,urou,-1,.True.,.True.)

          end if

          ! Get number of entries in this file
          read(100,*,err=1100) nfentry

          ! For each entry
          do ifentry=1,nfentry

            ! Get K, Q, and input nz
            read(100,*,err=1100) K, iQ, fnz

            ! Check setable
            if (setable(iQ,K)) then
              setable(iQ,K) = .False.
            else
              write(umsg,'(A,1x,i1,1x,i1)') &
                'You have specified more than one anysotropy '// &
                'factor for the multipole',K,iQ
              urou = 'initialize_asym'
              call abortedS(umsg,urou,-1,.True.,.True.)
            end if

            ! If only one height
            if (fnz.eq.1) then

              ! Read factor
              read(100,*,err=1100) daux1,daux2

              !
              ! Set factor

              ! Get other Q
              if (iQ.gt.0) then

                JKQa(iQ,K,:) = dcmplx(daux1, daux2)
                JKQa(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JKQa(iQ,K,:))

              ! Q = 0, ensure real
              else

                JKQa(0,K,:) = dcmplx(daux1, 0d0)

              end if ! Q > 0 or Q = 0

            ! If the correct number of heights
            else if (fnz.eq.nz) then

              ! For each height
              do iz=1,nz

                ! Read factor
                read(100,*,err=1100) daux1,daux2

                ! Check in limits
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                !
                ! Set factor

                ! Get other Q
                if (iQ.gt.0) then

                  JKQa(iQ,K,iz) = dcmplx(daux1, daux2)
                  JKQa(-iQ,K,iz) = Flgsg%sg(iQ)*conjg(JKQa(iQ,K,iz))

                ! Q = 0, ensure real
                else

                  JKQa(0,K,iz) = dcmplx(daux1, 0d0)

                end if ! Q > 0 or Q = 0

              end do ! Height nodes

            ! None of them
            else

              write(umsg,'(A,1x,i3,A,1x,i3)') &
                'You have specified the wrong number of '// &
                'height nodes:',fnz,'. Must be:',nz
              urou = 'initialize_asym'
              call abortedS(umsg,urou,-1,.True.,.True.)

            end if

          end do ! Entry in the current file

          ! Close file
          close(100)

          ! Delete temporal input file
          if (pid.eq.0) then
            call system('rm tmp_asym_'//Input%ID)
            umsg = ' - Asymmetry input '// &
                   trim(Input%asym_fil(ientry)%str)//' read'
            call verbose
          end if

        end do ! File entries

      end if ! Type of run

      ! If MPI and dynamic, share
      if (dyn.and.PRD.and.nproc.gt.1) then

        ! Size
        iz = 5*2*Rnz

        ! Share
        call MPI_BCAST(JKQa(-2,1,Rz0),iz,MPI_DOUBLE_COMPLEX,0, &
                       MPI_COMM_RT,ios)

      end if

      ! Control
      call control

      return

1000  umsg = ' - Error opening asymmetry file'
      urou = 'initialize_asym'
      call abortedS(umsg,urou,-1,.True.,.True.)
1100  umsg = ' - Error reading asymmetry file'
      urou = 'initialize_asym'
      close(100)
      call abortedS(umsg,urou,-1,.True.,.True.)

      end subroutine initialize_asym

!#####################################################################
!#####################################################################
!#####################################################################

      end module initialize_mod
