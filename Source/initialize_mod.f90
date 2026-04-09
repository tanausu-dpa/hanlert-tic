      !> Initialization of radiation field tensors
      module initialize_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     09/04/2026 V4.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/04/2026:    V4.0.4 - Changed the call to python to python3
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
!  initialize
!    Allocate and initialize the Stokes parameters and radiation field
!  tensors for the polarization problem
!
!  initializeI
!    Allocate and initialize the intensity parameters and mean
!  intensity
!
!  initialize_failread
!    Initialize the Stokes parameters and JKQC radiation field tensors
!  if they could not be read from the solution file
!
!  initializeI_failread
!    Initialize the intensity and J00C mean intensity if they could
!  not be read from the solution file
!
!  initialize_asym
!    Initialize factor for ad-hoc asymmetries in radiation field
!  tensors
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

      !> Allocate and initialize the Stokes parameters and radiation
      !! field tensors for the polarization problem\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates\n
      !!            mode(character): Mode of operation (solving,
      !!                             reading, or both)
      subroutine initialize(Frec,Geom,Atmo,Stokes,JKQ,JKQS, &
                            JKQC,J00P,mode)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      character(len=1), intent(in)::  mode
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Stokes
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: J00P
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQC

      ! Local

      integer:: iz,ifreq

      double precision:: omg


      !
      ! Allocations
      !

      ! If for any reason the Stokes parameters are needed at
      ! all heights
      if (KSTK) then

        ! Allocate and save limits
        allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
        giz0 = Rz0
        giz1 = Rz1

      ! If only needed at one height
      else

        ! Allocate and save limits
        allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
        giz0 = Rz0
        giz1 = Rz0+1

      end if ! If Stokes is needed at all heights

      ! Add memory and initialize
      RRAMc = RRAMc + 1d-6*sizeof(Stokes)
      Stokes = 0d0

      ! Atoms
      if (nA.gt.0) then

        ! JKQ integrated over absorption profiles
        allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(JKQ)
        JKQ = cZero

        ! JKQ integrated over emission profiles
        allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(JKQS)
        JKQS = cZero

        ! J00 for photoionizations
        allocate(J00P(nxphot,2,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00P)
        J00P = 0d0

      end if ! Atoms

      ! JKQ frequency dependent
      allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
      RRAMc = RRAMc + 1d-6*sizeof(JKQC)

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

        end do ! Frequencies
      end do ! Heights

      ! Control
      call control

      return

      end subroutine initialize

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate and initialize the intensity parameters and mean
      !! intensity\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Stokes(double(:,:,:,:)): Intensity\n
      !!         J00(double(:,:)): Mean intensity integrated over the
      !!                           absorption profile\n
      !!        J00S(double(:,:)): Mean intensity integrated over the
      !!                           emission profile\n
      !!        J00C(double(:,:)): Mean intensity with frequency
      !!                           dependence\n
      !!      J00P(double(:,:,:)): Intensity integrals in the
      !!                           photoionization rates\n
      !!          mode(character): Mode of operation (solving,
      !!                           reading, or both)
      subroutine initializeI(Frec,Geom,Atmo,Stokes,J00,J00S, &
                             J00C,J00P,mode)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      character(len=1), intent(in)::  mode
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(out):: Stokes
      double precision, dimension(:,:), allocatable, intent(out):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: J00P

      ! Local

      integer:: iz,ifreq

      double precision:: omg


      !
      ! Allocations
      !

      ! If for any reason the intensity is needed at all heights
      if (KSTK) then

        ! Allocate and save limits
        allocate(Stokes(nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
        giz0 = Rz0
        giz1 = Rz1

      ! If only needed at one height
      else

        ! Allocate and save limits
        allocate(Stokes(nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
        giz0 = Rz0
        giz1 = Rz0+1

      end if ! If intensity is needed at all heights

      ! Add memory and initialize
      RRAMc = RRAMc + 1d-6*sizeof(Stokes)
      Stokes = 0d0

      ! Atoms
      if (nA.gt.0) then

        ! J00 integrated over absorption profiles
        allocate(J00(nxt,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00)
        J00 = 0d0

        ! J00 integrated over emission profiles
        allocate(J00S(nxt,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00S)
        J00S = 0d0

        ! J00 for photoionizations
        allocate(J00P(nxphot,2,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00P)
        J00P = 0d0

      end if ! Atoms

      ! J00 frequency dependent
      allocate(J00C(nfreq,Rz0:Rz1))
      RRAMc = RRAMc + 1d-6*sizeof(J00C)

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
          J00C(ifreq,iz) = planck(omg,Atmo%T(iz))

        end do ! Frequencies
      end do ! Heights

      ! Control
      call control

      return

      end subroutine initializeI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the Stokes parameters and JKQC radiation field
      !! tensors if they could not be read from the solution file\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence
      subroutine initialize_failread(Frec,Atmo,Stokes,JKQC)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Stokes
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQC

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

        end do ! Frequencies
      end do ! Heights

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

      end if ! Need Stokes

      ! Control
      call control

      return

      end subroutine initialize_failread

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the intensity and J00C mean intensity if they
      !! could not be read from the solution file\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Stokes(double(:,:,:,:)): Intensity\n
      !!        J00C(double(:,:)): Mean intensity with frequency
      !!                           dependence
      subroutine initializeI_failread(Frec,Atmo,Stokes,J00C)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: Stokes
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C

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
      end do ! Heights

      ! Need Stokes
      if (KSTK) then

        ! For each height
        do iz=Rz0,Rz1

          ! For each frequency
          do ifreq=1,nfreq

            ! Initialize with planckian
            Stokes(ifreq,:,:,iz) = dble(J00C(ifreq,iz))

          end do ! Frequencies
        end do ! Heights

      end if ! Need Stokes

      ! Control
      call control

      return

      end subroutine initializeI_failread

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize factor for ad-hoc asymmetries in radiation field
      !! tensors\n
      !!     Input(Input_class): Structure with configuration data\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                         J-symbols\n
      !!       JKQin(double(:)): Data with ad-hoc JKQ tensors\n
      !!  JKQa(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                         field tensors\n
      !!          asym(logical): If the are ad-hoc JKQ tensors
      subroutine initialize_asym(Input,Flgsg,JKQin,JKQa,asym)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(out):: asym
      double precision, dimension(:), allocatable, intent(in):: JKQin
      complex(kind=8), dimension(:,:,:), &
                       allocatable, intent(out):: JKQa

      ! Local

      logical, dimension(0:2,1:2):: setable

      integer:: iz,ientry,K,iQ,ifentry,nfentry,fnz,ios

      double precision:: daux1,daux2


      ! Initialize
      asym = Input%nasym.gt.0

      ! If a slave CPU or no inputs
      if (pid.gt.0.or.Input%nasym.le.0) then

        ! If there are inputs (so this is a slave) and there is PRD
        ! PRD and either dynamics or angle-dependent
        if (Input%nasym.ge.1.and.(dyn.or..not.AV).and.PRD) then

          ! Allocate asymmetry input
          allocate(JKQa(-2:2,1:2,Rz0:Rz1))

          ! Size
          iz = 5*2*Rnz

          ! Master is going to share
          call MPI_BCAST(JKQa(-2,1,Rz0),iz,MPI_DOUBLE_COMPLEX,0, &
                         MPI_COMM_RT,ios)

        ! There is no input, it is static or angle-averaged, or
        ! not even PRD
        else

          ! Allocate dummy variable
          allocate(JKQa(1,1,1))

        end if ! Need to share

        ! Memory count
        if (allocated(JKQa)) RRAMc = RRAMc + 1d-6*sizeof(JKQa)

        ! Control and leave
        call control
        return

      end if ! There is something to read or it is a slave

      !
      ! Only the master with available data should reach this
      !

      ! Allocate asymmetry input and add to RAM
      allocate(JKQa(-2:2,1:2,Rz0:Rz1))
      RRAMc = RRAMc + 1d-6*sizeof(JKQa)

      ! Initialize to zero
      JKQa = cZero

      !
      ! 1.5D case (input buffer)
      !
      if (allocated(JKQin)) then

        ! For every considered height
        do iz=Rz0,Rz1

          ! Convert to JKQa array
          JKQa(0,1,iz) = dcmplx(JKQin(iz),0d0)
          JKQa(1,1,iz) = dcmplx(JKQin(nz+iz),JKQin(2*nz+iz))
          JKQa(0,2,iz) = dcmplx(JKQin(3*nz+iz),0d0)
          JKQa(1,2,iz) = dcmplx(JKQin(4*nz+iz),JKQin(5*nz+iz))
          JKQa(2,2,iz) = dcmplx(JKQin(6*nz+iz),JKQin(7*nz+iz))

        end do ! Heights

        ! Get negatives from conjugation properties
        JKQa(-1,1,:) = -conjg(JKQa(1,1,:))
        JKQa(-1,2,:) = -conjg(JKQa(1,2,:))
        JKQa(-2,2,:) =  conjg(JKQa(2,2,:))

      !
      ! 1D case
      !
      else

        ! Initialize set
        setable = .True.

        ! For each numberic entry
        do ientry=1,Input%nasym_num

          ! Get K and Q
          K = nint(dble(Input%asym_num(1,ientry)))
          iQ = nint(dimag(Input%asym_num(1,ientry)))

          ! Check setable
          if (setable(iQ,K)) then

            ! Flag as not setable anymore
            setable(iQ,K) = .False.

          ! No setable
          else

            ! Error
            write(umsg,'(A,1x,i1,1x,i1)') &
              'You have specified more than one anysotropy '// &
              'factor for the multipole',K,iQ
            urou = 'initialize_asym'
            call abortedS(umsg,urou,.True.,.True.)

          end if ! Check setable

          !
          ! Set factor
          !

          ! If asymmetry
          if (iQ.gt.0) then

            ! Get from input
            JKQa(iQ,K,:) = Input%asym_num(2,ientry)

            ! Get other Q from conjugation properties
            JKQa(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JKQa(iQ,K,:))

          ! Q = 0
          else

            ! Ensure real
            JKQa(0,K,:) = dcmplx(dble(Input%asym_num(2,ientry)), 0d0)

          end if ! Q > 0 or Q = 0

        end do ! Numerical entries

        !
        ! For each file
        !
        do ientry=1,Input%nasym_fil

          ! Deal with the file in python
          call system('python3 '//trim(Input%source)//'rasym.py '// &
                       trim(Input%asym_fil(ientry)%str)// &
                      ' '//Input%ID//' '//verbosef)

          ! Open the file
          open(100, file='tmp_asym_'//Input%ID, status='old', &
               iostat=ios, err=1000)

          ! Read if success
          read (100,*,err=1100) ios

          ! If no correct file
          if (ios.lt.0) then

            ! Abort
            umsg = 'Problem translating the asymmetry file'
            urou = 'initialize_asym'
            call abortedS(umsg,urou,.True.,.True.)

          end if ! Reading issue

          ! Get number of entries in this file
          read(100,*,err=1100) nfentry

          ! For each entry
          do ifentry=1,nfentry

            ! Get K, Q, and input nz
            read(100,*,err=1100) K, iQ, fnz

            ! Check setable
            if (setable(iQ,K)) then

              ! Flag as not setable
              setable(iQ,K) = .False.

            ! If was already set
            else

              ! Issue error
              write(umsg,'(A,1x,i1,1x,i1)') &
                'You have specified more than one anysotropy '// &
                'factor for the multipole',K,iQ
              urou = 'initialize_asym'
              call abortedS(umsg,urou,.True.,.True.)

            end if ! If can be set

            ! If only one height
            if (fnz.eq.1) then

              ! Read factor
              read(100,*,err=1100) daux1,daux2

              ! If Q > 0
              if (iQ.gt.0) then

                ! Save JKQ
                JKQa(iQ,K,:) = dcmplx(daux1, daux2)

                ! Get other sign from conjugation properties
                JKQa(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JKQa(iQ,K,:))

              ! Q = 0
              else

                ! Ensure real
                JKQa(0,K,:) = dcmplx(daux1, 0d0)

              end if ! Q > 0 or Q = 0

            ! If the correct number of heights
            else if (fnz.eq.nz) then

              ! For each height
              do iz=1,nz

                ! Read factor
                read(100,*,err=1100) daux1,daux2

                ! Check within considered limits
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                ! If Q > 0
                if (iQ.gt.0) then

                  ! Set JKQ
                  JKQa(iQ,K,iz) = dcmplx(daux1, daux2)

                  ! Compute other sign from conjugation properties
                  JKQa(-iQ,K,iz) = Flgsg%sg(iQ)*conjg(JKQa(iQ,K,iz))

                ! Q = 0
                else

                  ! Ensure real
                  JKQa(0,K,iz) = dcmplx(daux1, 0d0)

                end if ! Q > 0 or Q = 0

              end do ! Heights

            ! Wrong number of heights
            else

              ! Issue error
              write(umsg,'(A,1x,i3,A,1x,i3)') &
                'You have specified the wrong number of '// &
                'height nodes:',fnz,'. Must be:',nz
              urou = 'initialize_asym'
              call abortedS(umsg,urou,.True.,.True.)

            end if ! Number of heights in input file

          end do ! Entry in the current file

          ! Close file
          close(100)

          ! If master
          if (pid.eq.0) then

            ! Delete temporal input file
            call system('rm tmp_asym_'//Input%ID)
            umsg = ' - Asymmetry input '// &
                   trim(Input%asym_fil(ientry)%str)//' read'
            call verbose

          end if ! Master

        end do ! File entries

      end if ! Type of run

      ! If MPI, PRD, and dynamic or angle-dependent
      if ((dyn.or..not.AV).and.PRD.and.nproc.gt.1) then

        ! Size
        iz = 5*2*Rnz

        ! Share with rest of CPU in RT calculation
        call MPI_BCAST(JKQa(-2,1,Rz0),iz,MPI_DOUBLE_COMPLEX,0, &
                       MPI_COMM_RT,ios)

      end if ! MPI, PRD, and dynamic

      ! Control
      call control

      return

1000  umsg = ' - Error opening asymmetry file'
      urou = 'initialize_asym'
      call abortedS(umsg,urou,.True.,.True.)
1100  umsg = ' - Error reading asymmetry file'
      urou = 'initialize_asym'
      close(100)
      call abortedS(umsg,urou,.True.,.True.)

      end subroutine initialize_asym

!#####################################################################
!#####################################################################
!#####################################################################

      end module initialize_mod
