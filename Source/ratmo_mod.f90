      !> Reading of atmospheric data
      module ratmo_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!  Start:
!     17/04/2017
!  Last version:
!     29/01/2026 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     29/01/2026:    V4.0.3 - Added FALA, FALF and MCO models to the
!                             hard-coded models in gAtmo (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    It allows for helium density input, but the code does nothing
!  with those
!
!    It allows as an input log g, but the it is not used in synthesis
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
!  iAtmo_p
!    Nullify pointers in Atmo_class structure
!
!  cBfield
!    Create a copy of a the magnetic field stratification
!
!  cAtmo
!    Create a copy of a model atmosphere
!
!  dAtmo
!    Create a dummy copy of a model atmosphere with dimension 1
!
!  rAtmo
!    Read a 1D atmospheric model from an ASCII file
!
!  gAtmo
!    Setup a model atmosphere with the FALC or the FALP stratification
!
!  gAtmo_Strat
!    Generate an optical depth stratification from parameters
!
!  rAtmo_frombuffer
!    Interpret the atmospheric model data in a vectorial buffer
!
!  rAtmo_cle_prep
!    Prepare the Atmo structure in the main loop of the CLE synthesis
!
!  rAtmo_cle_init
!    Interpret the position in space of the current CLE node from the
!  vectorial buffer and setup the location of variables depending
!  on the type of model
!
!  rAtmo_cle
!    Interpret the atmospheric model data in a vectorial buffer in
!  the CLE synthesis
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : c , me , TINYB , kb , TINYVEL
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Nullify pointers in Atmo_class structure\n
      !!  Atmo(Atmo_class): Structure with atmospheric data
      subroutine iAtmo_p(Atmo)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo

      ! Nullify pointers
      nullify(Atmo%z,Atmo%T,Atmo%vmi,Atmo%vx,Atmo%vy,Atmo%vz)
      nullify(Atmo%Bx,Atmo%By,Atmo%Bz,Atmo%vxa,Atmo%vya,Atmo%vza)
      nullify(Atmo%zeros)

      end subroutine iAtmo_p

!#####################################################################
!#####################################################################
!#####################################################################

      !> Create a copy of a the magnetic field stratification\n
      !!  Bin(Bfield_class): Structure to copy from\n
      !!  Bou(Bfield_class): Structure to copy to
      subroutine cBfield(Bin,Bou)

      ! I/O

      type(Bfield_class), intent(in):: Bin
      type(Bfield_class), intent(out):: Bou


      ! Routine name
      urou = 'cBfield'

      ! Hard copy
      Bou = Bin

      ! Memory count
      if (allocated(Bou%Bstrength)) &
        MRAMc = MRAMc + 3d-6*sizeof(Bou%Bstrength)
      if (allocated(Bou%Blos)) &
        MRAMc = MRAMc + 3d-6*sizeof(Bou%Blos)

      end subroutine cBfield

!#####################################################################
!#####################################################################
!#####################################################################

      !> Create a copy of a model atmosphere\n
      !!  Ain(Atmo_class): Structure to copy from\n
      !!  Aou(Atmo_class): Structure to copy to
      subroutine cAtmo(Ain,Aou)

      ! I/O

      type(Atmo_class), intent(in):: Ain
      type(Atmo_class), intent(out):: Aou

      ! Local

      integer:: lnz


      ! Keep global
      lnz = nz
      nz = Ain%nz

      ! Routine name
      urou = 'cAtmo'

      ! Hard copy
      Aou = Ain

      ! Copy pointers
      if (associated(Ain%z)) then
        nullify(Aou%z)
        allocate(Aou%z(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%z)
        Aou%z = Ain%z
      end if
      if (associated(Ain%T)) then
        nullify(Aou%T)
        allocate(Aou%T(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%T)
        Aou%T = Ain%T
      end if
      if (associated(Ain%vmi)) then
        nullify(Aou%vmi)
        allocate(Aou%vmi(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vmi)
        Aou%vmi = Ain%vmi
      end if
      if (associated(Ain%vx)) then
        nullify(Aou%vx)
        allocate(Aou%vx(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vx)
        Aou%vx = Ain%vx
      end if
      if (associated(Ain%vy)) then
        nullify(Aou%vy)
        allocate(Aou%vy(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vy)
        Aou%vy = Ain%vy
      end if
      if (associated(Ain%vz)) then
        nullify(Aou%vz)
        allocate(Aou%vz(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vz)
        Aou%vz = Ain%vz
      end if
      if (associated(Ain%Bx)) then
        nullify(Aou%Bx)
        allocate(Aou%Bx(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%Bx)
        Aou%Bx = Ain%Bx
      end if
      if (associated(Ain%By)) then
        nullify(Aou%By)
        allocate(Aou%By(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%By)
        Aou%By = Ain%By
      end if
      if (associated(Ain%Bz)) then
        nullify(Aou%Bz)
        allocate(Aou%Bz(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%Bz)
        Aou%Bz = Ain%Bz
      end if
      if (associated(Ain%zeros)) then
        nullify(Aou%zeros)
        allocate(Aou%zeros(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Aou%zeros)
        Aou%zeros = Ain%zeros
      end if
      nullify(Aou%vxa)
      nullify(Aou%vya)
      nullify(Aou%vza)

      ! Allocs
      Aou%alloc_a = .True.
      Aou%alloc_b = .True.

      ! Return global
      nz = lnz

      ! Count memory in ele
      if (allocated(Aou%ele)) then
        do lnz=lbound(Aou%ele,1),ubound(Aou%ele,1)
          if (allocated(Aou%ele(lnz)%Ei)) &
            MRAMc = MRAMc + 1d-6*sizeof(Aou%ele(lnz)%Ei)
          if (allocated(Aou%ele(lnz)%pf)) &
            MRAMc = MRAMc + 1d-6*sizeof(Aou%ele(lnz)%pf)
          MRAMc = MRAMc + 1d-6*sizeof(Aou%ele(lnz))
        end do
      end if

      ! Count memory in arrays
      if (allocated(Aou%nht)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%nht)
      if (allocated(Aou%nhm)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%nhm)
      if (allocated(Aou%Pg)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%Pg)
      if (allocated(Aou%rho)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%rho)
      if (allocated(Aou%Pe)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%Pe)
      if (allocated(Aou%zalt)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%zalt)
      if (allocated(Aou%nHa)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%nHa)
      if (allocated(Aou%pT)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%pT)
      if (allocated(Aou%abund)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%abund)
      if (allocated(Aou%JKQin)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%JKQin)
      if (allocated(Aou%ne)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%ne)
      if (allocated(Aou%vlos)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vlos)
      if (allocated(Aou%vpos)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vpos)
      if (allocated(Aou%vphi)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%vphi)
      if (allocated(Aou%chi500)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%chi500)
      if (allocated(Aou%nh)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%nh)
      if (allocated(Aou%nhe)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%nhe)

      return

      end subroutine cAtmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Create a dummy copy of a model atmosphere with dimension 1\n
      !!  Ain(Atmo_class): Structure to copy from\n
      !!  Aou(Atmo_class): Structure to copy to
      subroutine dAtmo(Ain,Aou)

      ! I/O

      type(Atmo_class), intent(in):: Ain
      type(Atmo_class), intent(out):: Aou

      ! Local

      integer:: lnz


      ! Routine name
      urou = 'dAtmo'

      ! Dimension
      lnz = 1

      ! Generate pointers
      nullify(Aou%z)
      nullify(Aou%T)
      allocate(Aou%T(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%T)
      nullify(Aou%vmi)
      nullify(Aou%vx)
      nullify(Aou%vy)
      nullify(Aou%vz)
      nullify(Aou%Bx)
      nullify(Aou%By)
      nullify(Aou%Bz)
      nullify(Aou%zeros)
      nullify(Aou%vxa)
      nullify(Aou%vya)
      nullify(Aou%vza)

      ! Allocs
      Aou%alloc_a = .True.
      Aou%alloc_b = .False.

      ! Count memory in arrays
      allocate(Aou%nht(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%nht)
      allocate(Aou%nhm(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%nhm)
      allocate(Aou%Pg(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%Pg)
      allocate(Aou%rho(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%rho)
      allocate(Aou%Pe(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%Pe)
      allocate(Aou%nHa(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%nHa)
      Aou%NT = Ain%NT
      Aou%pT = Ain%pT
      if (allocated(Aou%pT)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%pT)
      Aou%abund = Ain%abund
      if (allocated(Aou%abund)) &
        MRAMc = MRAMc + 1d-6*sizeof(Aou%abund)
      allocate(Aou%ne(lnz))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%ne)
      allocate(Aou%nh(lnz,6))
      MRAMc = MRAMc + 1d-6*sizeof(Aou%nh)

      ! Copy ele
      Aou%ele = Ain%ele
      Aou%nele = Ain%nele

      ! Count memory in ele
      if (allocated(Aou%ele)) then
        do lnz=lbound(Aou%ele,1),ubound(Aou%ele,1)
          if (allocated(Aou%ele(lnz)%Ei)) &
            MRAMc = MRAMc + 1d-6*sizeof(Aou%ele(lnz)%Ei)
          if (allocated(Aou%ele(lnz)%pf)) &
            MRAMc = MRAMc + 1d-6*sizeof(Aou%ele(lnz)%pf)
          MRAMc = MRAMc + 1d-6*sizeof(Aou%ele(lnz))
        end do
      end if

      return

      end subroutine dAtmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read a 1D atmospheric model from an ASCII file\n
      !!  filename(character(:)): Name of the file to read\n
      !!    source(character(:)): Path to the source code folder\n
      !!        ID(character(:)): ID of this run\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!            fvmi(double): Forced microturbulence value
      subroutine rAtmo(filename,source,ID,Atmo,fvmi)

      ! I/O

      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      type(Atmo_class), intent(inout):: Atmo
      double precision, intent(in):: fvmi

      ! Local

      integer:: iz, ios

      double precision:: ikbcgs
      double precision, dimension(6):: ddump6
      double precision, dimension(4):: ddump4


      ! Routine name
      urou = 'rAtmo'

      ! Constant
      ikbcgs = 1d-7/kb

      ! Python translation of the model
      if(gpid.eq.0) call system('python '//trim(source)// &
                               'ratmo.py '//trim(filename)//' '// &
                               ID//' '//verbosef)

      ! Wait for the master to finish
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! Open the translated temporal file
      open (100,file='tmp_atmo_'//ID,status='old',iostat=ios,err=1000)

      ! Success
      read (100,*,err=1100) ios

      ! If no correct reading of file
      if (ios.lt.0) then

        ! Abort
        umsg = 'Problem translating the atmospheric file'
        goto 1200

      end if

      ! Nullify pointers
      nullify(Atmo%zeros,Atmo%Bx,Atmo%By,Atmo%Bz)

      ! Type of scale, reference frequency, log(g), and number of
      ! nodes
      read (100,*,err=1100) Atmo%scal
      read (100,*,err=1100) Atmo%tfreq
      read (100,*,err=1100) Atmo%logg
      read (100,*,err=1100) Atmo%nZ

      ! If optical depth scale
      ztau = Atmo%scal.eq.'T'

      !
      ! Allocations
      !

      ! Zeros
      allocate(Atmo%zeros(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%zeros)
      Atmo%zeros = 0d0

      ! Height/tau axis
      allocate(Atmo%z(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%z)

      ! Temperature
      allocate(Atmo%T(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%T)

      ! Microturbulence
      allocate(Atmo%vmi(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vmi)

      ! Electron density
      allocate(Atmo%ne(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%ne)

      ! Total hydrogen density
      allocate(Atmo%nht(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nht)
      Atmo%nht = 0d0

      ! Atomic hydrogen density
      allocate(Atmo%nha(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nha)
      Atmo%nha = 0d0

      ! H^- density
      allocate(Atmo%nhm(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhm)
      Atmo%nhm = 0d0

      ! Velocity vector (x,y,z)
      allocate(Atmo%vx(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vx)
      allocate(Atmo%vy(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vy)
      allocate(Atmo%vz(Atmo%nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vz)

      ! H density (5 HI levels + p+ density)
      allocate(Atmo%nh(Atmo%nZ,6))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nh)

      ! Read type of density
      read (100,*,err=1100) Atmo%typo

      !
      ! Read themodynamics
      !

      ! For each height
      do iz=1,Atmo%nZ

        ! Read height, temperature, electron density, z component of
        ! velocity, microturbulence, horizontal components of the
        ! velocity
        read (100,*,err=1100) Atmo%z(iz), Atmo%T(iz), Atmo%ne(iz), &
                              Atmo%vz(iz), Atmo%vmi(iz), &
                              Atmo%vx(iz), Atmo%vy(iz)

      end do ! heights

      ! If forced micro, apply it
      if (fvmi.ge.0d0) Atmo%vmi = fvmi

      !
      ! If number densities
      !
      if (Atmo%typo.eq.0) then

        !
        ! Read hydrogen densityes
        !

        ! For each height
        do iz=1,Atmo%nZ

          ! Read in a dump 6 elements vector
          read (100,*,err=1100) ddump6

          ! Store the individual densities
          Atmo%nh(iz,:) = ddump6

          ! And the total in atomic
          Atmo%nha(iz) = sum(ddump6)

          ! And whole
          Atmo%nht(iz) = Atmo%nha(iz)

        end do ! For each height


        !
        ! Read helium densityes
        !

        ! Check if there is helium quantities
        read (100,*,err=1100) ios

        ! If there is no helium input
        if (ios.eq.0) then

          ! Just make a flag to know that there is no helium input
          allocate(Atmo%nhe(1,1))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
          Atmo%nhe(1,1) = -1

        ! If there is helium input
        else

          ! Allocate the density array
          allocate(Atmo%nhe(Atmo%nZ,4))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)

          ! For each height
          do iz=1,Atmo%nZ

            ! Read into a dump 4 elements vector
            read (100,*,err=1100) ddump4

            ! Store the helium density
            Atmo%nhe(iz,:) = ddump4

          end do ! heights

        end if ! there is no helium input

      ! No number densities
      else

        ! Initialize H to zero
        Atmo%nh = 0d0
        Atmo%nht = 0d0
        Atmo%nha = 0d0

        ! Just make a flag to know that there is no helium input
        allocate(Atmo%nhe(1,1))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
        Atmo%nhe(1,1) = -1

        ! Electron density does not require changes
        ! Atmo%typo.eq.1 no source needed

        !
        ! Electron pressure
        !
        if (Atmo%typo.eq.2) then

          ! Convert to electron number density
          Atmo%ne = Atmo%ne*ikbcgs/Atmo%T

          ! Typo is now 1
          Atmo%typo = 1

        !
        ! Electron mass density
        !
        else if (Atmo%typo.eq.3) then

          ! Convert to electron number density
          Atmo%ne = Atmo%ne*1d-3/me

          ! Typo is now 1
          Atmo%typo = 1

        !
        ! Gas pressure
        !
        else if (Atmo%typo.eq.4) then

          ! Move into gas pressure
          allocate(Atmo%Pg(Atmo%nz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)
          Atmo%Pg = Atmo%ne
          Atmo%ne = 0d0

        !
        ! Mass density
        !
        else if (Atmo%typo.eq.5) then

          ! Move into mass density
          allocate(Atmo%rho(Atmo%nz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%rho)
          Atmo%rho = Atmo%ne
          Atmo%ne = 0d0

        end if ! No electron number density
      end if ! Type of input densities

      !
      ! Unit conversions
      !

      ! Transform height to cgs
      if (.not.ztau) Atmo%z = Atmo%z*1d5

      ! Divide velocities by c (1d5*1d-11/cbar)
      Atmo%vx = Atmo%vx*1d-6/c
      Atmo%vy = Atmo%vy*1d-6/c
      Atmo%vz = Atmo%vz*1d-6/c
      Atmo%vmi = Atmo%vmi*1d-6/c

      ! Check if dynamic (yes if > 1m/s)
      dyn = maxval(Atmo%vx*Atmo%vx + Atmo%vy*Atmo%vy + &
                   Atmo%vz*Atmo%vz).gt.TINYVEL

      ! Allocs
      Atmo%alloc_a = .True.
      Atmo%alloc_b = .True.

      ! Close the temporal file
      close (100)

      ! Control that everything went fine
      call control

      ! Global master
      if(gpid.eq.0) then

        ! Delete temporal input atmospheric file
        call system('rm tmp_atmo_'//ID)

        ! Verbose
        umsg = ' - Atmosphere read: '//trim(filename)
        if (run_mode.eq.-1) then
          call verbosev
        else
          call verbose
        end if
      end if ! Global Master

      return

1000  umsg = 'Error opening atmospheric file'
      call aborted
      return
1100  umsg = 'Error reading atmospheric file'
1200  close(100)
      call aborted
      return

      end subroutine rAtmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Setup a model atmosphere with the FALC or the FALP
      !! stratification\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!     init(integer): What hard-coded model to setup
      subroutine gAtmo(Atmo,init)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      integer, intent(in):: init

      ! Local

      integer:: iz
      integer, target:: nzC = 70
      integer, target:: nzP = 66
      integer, target:: nzA = 70
      integer, target:: nzX = 80
      integer, target:: nzF = 70

      ! Pointers

      integer, pointer:: lnz

      double precision, dimension(:), pointer:: z,tau,T,ne,vmi
      double precision, dimension(:,:), pointer:: nh


      !
      !************** FAL-C atmosphere model
      !Fontenla et al. (1993), ApJ, 406, 319
      !

      ! Optical depth (log10(tau))
     !double precision, dimension(70), target:: TAUC= (/ &
     ! -10.00d0, -9.162d0, -8.867d0, -8.698d0, -8.508d0, &
     ! -8.367d0, -8.305d0, -8.247d0, -8.192d0, -8.138d0, &
     ! -8.111d0, -8.086d0, -8.060d0, -8.033d0, -8.006d0, &
     ! -7.979d0, -7.952d0, -7.921d0, -7.894d0, -7.875d0, &
     ! -7.846d0, -7.824d0, -7.774d0, -7.709d0, -7.383d0, &
     ! -7.020d0, -6.772d0, -6.598d0, -6.498d0, -6.453d0, &
     ! -6.409d0, -6.351d0, -6.282d0, -6.200d0, -6.082d0, &
     ! -6.001d0, -5.893d0, -5.777d0, -5.688d0, -5.593d0, &
     ! -5.509d0, -5.421d0, -5.331d0, -5.206d0, -5.092d0, &
     ! -4.969d0, -4.872d0, -4.766d0, -4.650d0, -4.518d0, &
     ! -4.330d0, -4.090d0, -3.842d0, -3.597d0, -3.340d0, &
     ! -3.045d0, -2.681d0, -2.328d0, -1.985d0, -1.650d0, &
     ! -1.322d0, -1.003d0, -0.694d0, -0.377d0,  0.037d0, &
     !  0.232d0,  0.479d0,  0.755d0,  1.030d0,  1.283d0 /)

      ! Geometrical height (km)
      double precision, dimension(70), target:: HGHC= (/ &
      2219.43d0, 2217.88d0, 2216.43d0, 2215.08d0, 2212.92d0, &
      2210.75d0, 2209.65d0, 2208.53d0, 2207.40d0, 2206.27d0, &
      2205.70d0, 2205.18d0, 2204.66d0, 2204.13d0, 2203.63d0, &
      2203.16d0, 2202.70d0, 2202.22d0, 2201.83d0, 2201.56d0, &
      2201.16d0, 2200.84d0, 2200.10d0, 2199.00d0, 2190.00d0, &
      2168.00d0, 2140.00d0, 2110.00d0, 2087.00d0, 2075.00d0, &
      2062.00d0, 2043.00d0, 2017.00d0, 1980.00d0, 1915.00d0, &
      1860.00d0, 1775.00d0, 1670.00d0, 1580.00d0, 1475.00d0, &
      1378.00d0, 1278.00d0, 1180.00d0, 1065.00d0,  980.00d0, &
       905.00d0,  855.00d0,  805.00d0,  755.00d0,  705.00d0, &
       650.00d0,  600.00d0,  560.00d0,  525.00d0,  490.00d0, &
       450.00d0,  400.00d0,  350.00d0,  300.00d0,  250.00d0, &
       200.00d0,  150.00d0,  100.00d0,   50.00d0,    0.00d0, &
       -20.00d0,  -40.00d0,  -60.00d0,  -80.00d0, -100.00d0 /)

      ! Temperature
      double precision, dimension(70), target:: TEMPC= (/ &
     102770d0, 98790d0, 94800d0, 90816d0, 83891d0, 75934d0, 71336d0, &
      66145d0, 60170d0, 53284d0, 49385d0, 45416d0, 41178d0, 36594d0, &
      32145d0, 27972d0, 27972d0, 20416d0, 17925d0, 16500d0, 15000d0, &
      14250d0, 13500d0, 13000d0, 12000d0, 11150d0, 10550d0,  9900d0, &
       9450d0,  9200d0,  8950d0,  8700d0,  8400d0,  8050d0,  7650d0, &
       7450d0,  7250d0,  7050d0,  6900d0,  6720d0,  6560d0,  6390d0, &
       6230d0,  6040d0,  5900d0,  5755d0,  5650d0,  5490d0,  5280d0, &
       5030d0,  4750d0,  4550d0,  4430d0,  4400d0,  4410d0,  4460d0, &
       4560d0,  4660d0,  4770d0,  4880d0,  4990d0,  5150d0,  5410d0, &
       5790d0,  6520d0,  6980d0,  7590d0,  8220d0,  8860d0,  9400d0 /)

      ! Electron density
      double precision, dimension(70), target:: NeC= (/ &
          6.560d09,  6.810d09,  7.070d09,  7.360d09,  7.910d09,  &
          8.660d09,  9.170d09,  9.820d09,  1.070d10,  1.190d10,  &
          1.280d10,  1.380d10,  1.510d10,  1.670d10,  1.880d10,  &
          2.120d10,  2.410d10,  2.770d10,  3.080d10,  3.290d10,  &
          3.530d10,  3.650d10,  3.730d10,  3.720d10,  3.560d10,  &
          3.810d10,  4.060d10,  4.240d10,  4.310d10,  4.340d10,  &
          4.360d10,  4.370d10,  4.390d10,  4.420d10,  4.540d10,  &
          4.700d10,  5.080d10,  5.740d10,  6.320d10,  7.040d10,  &
          7.860d10,  8.760d10,  9.800d10,  1.100d11,  1.180d11,  &
          1.230d11,  1.240d11,  1.080d11,  9.410d10,  8.280d10,  &
          9.030d10,  1.260d11,  1.770d11,  2.410d11,  3.280d11,  &
          4.660d11,  7.240d11,  1.120d12,  1.710d12,  2.600d12,  &
          3.930d12,  6.040d12,  9.890d12,  1.970d13,  7.680d13,  &
          1.730d14,  4.480d14,  1.050d15,  2.210d15,  3.870d15  /)

      ! Hydrogen density
      double precision, dimension(6,70), target:: NHC= reshape((/ &
          1.110d05,3.700d-02, 0d0, 0d0, 0d0, 5.480d09, &
          1.790d05,5.900d-02, 0d0, 0d0, 0d0, 5.700d09, &
          2.460d05,8.000d-02, 0d0, 0d0, 0d0, 5.940d09, &
          3.440d05,1.110d-01, 0d0, 0d0, 0d0, 6.210d09, &
          7.000d05,2.200d-01, 0d0, 0d0, 0d0, 6.730d09, &
          2.050d06,6.290d-01, 0d0, 0d0, 0d0, 7.450d09, &
          4.390d06, 1.330d00, 0d0, 0d0, 0d0, 7.940d09, &
          9.810d06, 2.970d00, 0d0, 0d0, 0d0, 8.560d09, &
          2.310d07, 7.090d00, 0d0, 0d0, 0d0, 9.410d09, &
          5.700d07, 1.830d01, 0d0, 0d0, 0d0, 1.060d10, &
          9.340d07, 3.120d01, 0d0, 0d0, 0d0, 1.140d10, &
          1.470d08, 5.190d01, 0d0, 0d0, 0d0, 1.230d10, &
          2.320d08, 8.800d01, 0d0, 0d0, 0d0, 1.350d10, &
          3.690d08, 1.530d02, 0d0, 0d0, 0d0, 1.500d10, &
          5.700d08, 2.610d02, 0d0, 0d0, 0d0, 1.690d10, &
          8.500d08, 4.230d02, 0d0, 0d0, 0d0, 1.910d10, &
          1.240d09, 6.620d02, 0d0, 0d0, 0d0, 2.180d10, &
          1.780d09, 1.020d03, 0d0, 0d0, 0d0, 2.510d10, &
          2.350d09, 1.400d03, 0d0, 0d0, 0d0, 2.810d10, &
          2.820d09, 1.720d03, 0d0, 0d0, 0d0, 3.010d10, &
          3.600d09, 2.260d03, 0d0, 0d0, 0d0, 3.250d10, &
          4.280d09, 2.750d03, 0d0, 0d0, 0d0, 3.380d10, &
          5.890d09, 4.000d03, 0d0, 0d0, 0d0, 3.470d10, &
          8.240d09, 5.990d03, 0d0, 0d0, 0d0, 3.480d10, &
          1.710d10, 1.780d04, 0d0, 0d0, 0d0, 3.300d10, &
          2.120d10, 2.590d04, 0d0, 0d0, 0d0, 3.540d10, &
          2.630d10, 2.960d04, 0d0, 0d0, 0d0, 3.770d10, &
          3.540d10, 3.150d04, 0d0, 0d0, 0d0, 3.920d10, &
          4.450d10, 3.230d04, 0d0, 0d0, 0d0, 3.980d10, &
          5.030d10, 3.270d04, 0d0, 0d0, 0d0, 4.000d10, &
          5.700d10, 3.320d04, 0d0, 0d0, 0d0, 4.010d10, &
          6.690d10, 3.400d04, 0d0, 0d0, 0d0, 4.010d10, &
          8.180d10, 3.540d04, 0d0, 0d0, 0d0, 4.010d10, &
          1.060d11, 3.800d04, 0d0, 0d0, 0d0, 4.030d10, &
          1.570d11, 4.370d04, 0d0, 0d0, 0d0, 4.190d10, &
          2.120d11, 4.960d04, 0d0, 0d0, 0d0, 4.410d10, &
          3.310d11, 6.130d04, 0d0, 0d0, 0d0, 4.870d10, &
          5.740d11, 8.010d04, 0d0, 0d0, 0d0, 5.580d10, &
          9.280d11, 9.690d04, 0d0, 0d0, 0d0, 6.160d10, &
          1.660d12, 1.180d05, 0d0, 0d0, 0d0, 6.890d10, &
          2.890d12, 1.420d05, 0d0, 0d0, 0d0, 7.660d10, &
          5.300d12, 1.690d05, 0d0, 0d0, 0d0, 8.530d10, &
          9.920d12, 2.000d05, 0d0, 0d0, 0d0, 9.500d10, &
          2.150d13, 2.350d05, 0d0, 0d0, 0d0, 1.050d11, &
          3.920d13, 2.630d05, 0d0, 0d0, 0d0, 1.100d11, &
          6.790d13, 2.740d05, 0d0, 0d0, 0d0, 1.110d11, &
          9.920d13, 2.700d05, 0d0, 0d0, 0d0, 1.090d11, &
          1.480d14, 2.120d05, 0d0, 0d0, 0d0, 8.540d10, &
          2.270d14, 1.380d05, 0d0, 0d0, 0d0, 6.180d10, &
          3.560d14, 7.100d04, 0d0, 0d0, 0d0, 3.490d10, &
          6.030d14, 3.160d04, 0d0, 0d0, 0d0, 1.350d10, &
          9.890d14, 1.860d04, 0d0, 0d0, 0d0, 5.460d09, &
          1.480d15, 1.490d04, 0d0, 0d0, 0d0, 2.970d09, &
          2.080d15, 1.800d04, 0d0, 0d0, 0d0, 2.570d09, &
          2.900d15, 2.670d04, 0d0, 0d0, 0d0, 2.760d09, &
          4.190d15, 5.150d04, 0d0, 0d0, 0d0, 3.730d09, &
          6.540d15, 1.420d05, 0d0, 0d0, 0d0, 6.700d09, &
          1.010d16, 3.800d05, 0d0, 0d0, 0d0, 1.210d10, &
          1.540d16, 1.040d06, 0d0, 0d0, 0d0, 2.380d10, &
          2.330d16, 2.730d06, 0d0, 0d0, 0d0, 4.940d10, &
          3.470d16, 6.950d06, 0d0, 0d0, 0d0, 1.100d11, &
          5.050d16, 2.120d07, 0d0, 0d0, 0d0, 3.220d11, &
          7.100d16, 8.960d07, 0d0, 0d0, 0d0, 1.400d12, &
          9.550d16, 5.070d08, 0d0, 0d0, 0d0, 7.490d12, &
          1.180d17, 6.180d09, 0d0, 0d0, 0d0, 6.000d13, &
          1.240d17, 2.150d10, 0d0, 0d0, 0d0, 1.530d14, &
          1.280d17, 8.630d10, 0d0, 0d0, 0d0, 4.250d14, &
          1.300d17, 2.900d11, 0d0, 0d0, 0d0, 1.020d15, &
          1.300d17, 8.230d11, 0d0, 0d0, 0d0, 2.180d15, &
          1.310d17, 1.790d12, 0d0, 0d0, 0d0, 3.830d15/), (/6,70/))

      ! Microturbulent velocity (km/s)
      double precision, dimension(70), target:: VmiC= (/ &
        11.90d0, 11.80d0, 11.70d0, 11.60d0, 11.50d0, 11.30d0, &
        11.20d0, 11.00d0, 10.80d0, 10.60d0, 10.50d0, 10.30d0, &
        10.10d0, 9.930d0, 9.700d0, 9.470d0, 9.220d0, 8.950d0, &
        8.750d0, 8.620d0, 8.470d0, 8.380d0, 8.280d0, 8.190d0, &
        7.960d0, 7.780d0, 7.600d0, 7.380d0, 7.220d0, 7.130d0, &
        7.040d0, 6.920d0, 6.760d0, 6.550d0, 6.220d0, 5.970d0, &
        5.610d0, 5.210d0, 4.880d0, 4.520d0, 4.200d0, 3.860d0, &
        3.500d0, 3.030d0, 2.640d0, 2.300d0, 2.060d0, 1.820d0, &
        1.570d0, 1.330d0, 1.070d0, 0.865d0, 0.725d0, 0.628d0, &
        0.555d0, 0.501d0, 0.477d0, 0.501d0, 0.627d0, 0.867d0, &
        1.170d0, 1.480d0, 1.750d0, 1.970d0, 2.000d0, 2.000d0, &
        2.000d0, 2.000d0, 2.000d0, 2.000d0 /)
      ! Microturbulent velocity (scaled)
     !double precision, dimension(70), target:: VmiC= (/ &
     !  3.969d-05, 3.936d-05, 3.903d-05, 3.869d-05, 3.836d-05, &
     !  3.769d-05, 3.736d-05, 3.669d-05, 3.602d-05, 3.536d-05, &
     !  3.502d-05, 3.436d-05, 3.369d-05, 3.312d-05, 3.236d-05, &
     !  3.159d-05, 3.075d-05, 2.985d-05, 2.919d-05, 2.875d-05, &
     !  2.825d-05, 2.795d-05, 2.762d-05, 2.732d-05, 2.655d-05, &
     !  2.595d-05, 2.535d-05, 2.462d-05, 2.408d-05, 2.378d-05, &
     !  2.348d-05, 2.308d-05, 2.255d-05, 2.185d-05, 2.075d-05, &
     !  1.991d-05, 1.871d-05, 1.738d-05, 1.628d-05, 1.508d-05, &
     !  1.401d-05, 1.288d-05, 1.167d-05, 1.011d-05, 8.806d-06, &
     !  7.672d-06, 6.871d-06, 6.071d-06, 5.237d-06, 4.436d-06, &
     !  3.569d-06, 2.885d-06, 2.418d-06, 2.095d-06, 1.851d-06, &
     !  1.671d-06, 1.591d-06, 1.671d-06, 2.091d-06, 2.892d-06, &
     !  3.903d-06, 4.937d-06, 5.837d-06, 6.571d-06, 6.671d-06, &
     !  6.671d-06, 6.671d-06, 6.671d-06, 6.671d-06, 6.671d-06 /)
      !
      !Fontenla et al. (1993), ApJ, 406, 319
      !************** FAL-C atmosphere model
      !



      !
      !************** FAL-P atmosphere model
      !Fontenla et al. (1993), ApJ, 406, 319
      !

      ! Optical depth
     !double precision, dimension(66), target:: TAUP= (/ &
     !    -10.00d0, -8.776d0, -8.483d0, -8.323d0, -8.210d0, &
     !    -8.133d0, -8.070d0, -8.017d0, -7.996d0, -7.974d0, &
     !    -7.951d0, -7.930d0, -7.903d0, -7.891d0, -7.872d0, &
     !    -7.858d0, -7.837d0, -7.816d0, -7.778d0, -7.755d0, &
     !    -7.733d0, -7.711d0, -7.683d0, -7.644d0, -7.619d0, &
     !    -7.590d0, -7.554d0, -7.503d0, -7.475d0, -7.390d0, &
     !    -7.268d0, -7.046d0, -6.836d0, -6.411d0, -6.061d0, &
     !    -5.752d0, -5.396d0, -5.146d0, -4.978d0, -4.839d0, &
     !    -4.753d0, -4.574d0, -4.459d0, -4.340d0, -4.249d0, &
     !    -4.018d0, -3.901d0, -3.775d0, -3.624d0, -3.432d0, &
     !    -3.190d0, -2.906d0, -2.597d0, -2.279d0, -1.956d0, &
     !    -1.636d0, -1.321d0, -1.010d0, -0.701d0, -0.375d0, &
     !     0.031d0,  0.223d0,  0.466d0,  0.740d0,  1.015d0, &
     !     1.267d0  /)

      ! Geometrical height (km)
      double precision, dimension(66), target:: HGHP= (/ &
      1741.970d0, 1741.570d0, 1741.220d0, 1740.930d0, 1740.680d0, &
      1740.490d0, 1740.330d0, 1740.200d0, 1740.150d0, 1740.100d0, &
      1740.050d0, 1740.010d0, 1739.960d0, 1739.940d0, 1739.910d0, &
      1739.890d0, 1739.860d0, 1739.830d0, 1739.780d0, 1739.750d0, &
      1739.720d0, 1739.690d0, 1739.650d0, 1739.590d0, 1739.550d0, &
      1739.500d0, 1739.430d0, 1739.320d0, 1739.250d0, 1739.000d0, &
      1738.500d0, 1737.000d0, 1734.500d0, 1723.000d0, 1700.000d0, &
      1660.000d0, 1575.030d0, 1475.030d0, 1380.000d0, 1280.000d0, &
      1210.000d0, 1065.000d0,  980.000d0,  905.000d0,  855.000d0, &
       750.000d0,  700.000d0,  650.000d0,  600.000d0,  550.000d0, &
       500.000d0,  450.000d0,  400.000d0,  350.000d0,  300.000d0, &
       250.000d0,  200.000d0,  150.000d0,  100.000d0,   50.000d0, &
         0.000d0,  -20.000d0,  -40.000d0,  -60.000d0,  -80.000d0, &
      -100.000d0 /)

      ! Temperature
      double precision, dimension(66), target:: TEMPP= (/ &
          121000d0, 111000d0, 101000d0, 91000d0, 81000d0, 71000d0, &
           61000d0,  51000d0,  46000d0, 41000d0, 36000d0, 31000d0, &
           26000d0,  23500d0,  21000d0, 19000d0, 17000d0, 15000d0, &
           13000d0,  12000d0,  11500d0, 11000d0, 10500d0, 10000d0, &
            9800d0,   9600d0,   9400d0,  9200d0,  9100d0,  8900d0, &
            8800d0,   8700d0,   8600d0,  8500d0,  8400d0,  8250d0, &
            8000d0,   7700d0,   7420d0,  7150d0,  6980d0,  6600d0, &
            6390d0,   6220d0,   6090d0,  5740d0,  5480d0,  5220d0, &
            5070d0,   4960d0,   4910d0,  4900d0,  4940d0,  5000d0, &
            5070d0,   5170d0,   5280d0,  5400d0,  5550d0,  5880d0, &
            6520d0,   6980d0,   7590d0,  8220d0,  8860d0,  9400d0  /)

      ! Electron density
      double precision, dimension(66), target:: NeP= (/ &
          6.01d10, 6.59d10, 7.23d10, 7.99d10, 8.91d10, 1.01d11, &
          1.16d11, 1.36d11, 1.50d11, 1.67d11, 1.87d11, 2.14d11, &
          2.49d11, 2.73d11, 3.00d11, 3.27d11, 3.56d11, 3.91d11, &
          4.29d11, 4.51d11, 4.62d11, 4.73d11, 4.84d11, 4.86d11, &
          4.84d11, 4.79d11, 4.68d11, 4.50d11, 4.39d11, 4.07d11, &
          3.66d11, 3.36d11, 3.17d11, 2.98d11, 3.09d11, 3.44d11, &
          4.11d11, 4.57d11, 4.79d11, 4.84d11, 4.91d11, 4.85d11, &
          4.73d11, 4.75d11, 4.64d11, 3.84d11, 3.02d11, 2.50d11, &
          2.67d11, 3.27d11, 4.48d11, 6.44d11, 9.51d11, 1.41d12, &
          2.10d12, 3.13d12, 4.70d12, 7.11d12, 1.11d13, 2.20d13, &
          7.58d13, 1.71d14, 4.42d14, 1.04d15, 2.19d15, 3.82d15  /)

      ! Hydrogen density
      double precision, dimension(6,66), target:: NHP= reshape((/ &
          9.190d05, 2.980d00, 0d0, 0d0, 0d0, 5.020d10, &
          1.420d06, 4.620d00, 0d0, 0d0, 0d0, 5.530d10, &
          2.470d06, 7.980d00, 0d0, 0d0, 0d0, 6.110d10, &
          6.700d06, 2.140d01, 0d0, 0d0, 0d0, 6.810d10, &
          2.730d07, 8.600d01, 0d0, 0d0, 0d0, 7.700d10, &
          1.190d08, 3.780d02, 0d0, 0d0, 0d0, 8.820d10, &
          4.600d08, 1.550d03, 0d0, 0d0, 0d0, 1.030d11, &
          1.500d09, 5.740d03, 0d0, 0d0, 0d0, 1.220d11, &
          2.600d09, 1.090d04, 0d0, 0d0, 0d0, 1.350d11, &
          4.290d09, 1.990d04, 0d0, 0d0, 0d0, 1.500d11, &
          6.880d09, 3.550d04, 0d0, 0d0, 0d0, 1.680d11, &
          1.080d10, 6.180d04, 0d0, 0d0, 0d0, 1.920d11, &
          1.720d10, 1.070d05, 0d0, 0d0, 0d0, 2.250d11, &
          2.170d10, 1.390d05, 0d0, 0d0, 0d0, 2.460d11, &
          2.790d10, 1.830d05, 0d0, 0d0, 0d0, 2.710d11, &
          3.450d10, 2.280d05, 0d0, 0d0, 0d0, 2.960d11, &
          4.380d10, 2.880d05, 0d0, 0d0, 0d0, 3.260d11, &
          5.740d10, 3.720d05, 0d0, 0d0, 0d0, 3.610d11, &
          8.040d10, 5.090d05, 0d0, 0d0, 0d0, 4.040d11, &
          1.000d11, 6.240d05, 0d0, 0d0, 0d0, 4.270d11, &
          1.150d11, 7.120d05, 0d0, 0d0, 0d0, 4.420d11, &
          1.360d11, 8.320d05, 0d0, 0d0, 0d0, 4.550d11, &
          1.650d11, 9.910d05, 0d0, 0d0, 0d0, 4.690d11, &
          2.080d11, 1.220d06, 0d0, 0d0, 0d0, 4.740d11, &
          2.340d11, 1.350d06, 0d0, 0d0, 0d0, 4.720d11, &
          2.670d11, 1.510d06, 0d0, 0d0, 0d0, 4.670d11, &
          3.090d11, 1.690d06, 0d0, 0d0, 0d0, 4.560d11, &
          3.660d11, 1.890d06, 0d0, 0d0, 0d0, 4.370d11, &
          4.010d11, 2.000d06, 0d0, 0d0, 0d0, 4.260d11, &
          4.870d11, 2.150d06, 0d0, 0d0, 0d0, 3.920d11, &
          5.730d11, 2.080d06, 0d0, 0d0, 0d0, 3.490d11, &
          6.430d11, 1.810d06, 0d0, 0d0, 0d0, 3.170d11, &
          6.970d11, 1.590d06, 0d0, 0d0, 0d0, 2.980d11, &
          7.830d11, 1.450d06, 0d0, 0d0, 0d0, 2.840d11, &
          8.770d11, 1.530d06, 0d0, 0d0, 0d0, 3.030d11, &
          1.060d12, 1.730d06, 0d0, 0d0, 0d0, 3.400d11, &
          1.650d12, 2.080d06, 0d0, 0d0, 0d0, 4.080d11, &
          2.910d12, 2.270d06, 0d0, 0d0, 0d0, 4.530d11, &
          5.030d12, 2.350d06, 0d0, 0d0, 0d0, 4.740d11, &
          9.060d12, 2.330d06, 0d0, 0d0, 0d0, 4.780d11, &
          1.390d13, 2.360d06, 0d0, 0d0, 0d0, 4.830d11, &
          3.540d13, 2.310d06, 0d0, 0d0, 0d0, 4.740d11, &
          6.350d13, 2.220d06, 0d0, 0d0, 0d0, 4.580d11, &
          1.070d14, 2.240d06, 0d0, 0d0, 0d0, 4.530d11, &
          1.540d14, 2.140d06, 0d0, 0d0, 0d0, 4.350d11, &
          3.440d14, 1.440d06, 0d0, 0d0, 0d0, 3.300d11, &
          5.240d14, 8.220d05, 0d0, 0d0, 0d0, 2.240d11, &
          8.150d14, 4.370d05, 0d0, 0d0, 0d0, 1.360d11, &
          1.260d15, 3.480d05, 0d0, 0d0, 0d0, 9.830d10, &
          1.950d15, 3.280d05, 0d0, 0d0, 0d0, 7.390d10, &
          3.000d15, 4.020d05, 0d0, 0d0, 0d0, 6.570d10, &
          4.590d15, 5.940d05, 0d0, 0d0, 0d0, 6.780d10, &
          6.960d15, 1.100d06, 0d0, 0d0, 0d0, 8.730d10, &
          1.040d16, 2.200d06, 0d0, 0d0, 0d0, 1.240d11, &
          1.560d16, 4.550d06, 0d0, 0d0, 0d0, 1.870d11, &
          2.300d16, 1.050d07, 0d0, 0d0, 0d0, 3.260d11, &
          3.340d16, 2.460d07, 0d0, 0d0, 0d0, 6.100d11, &
          4.800d16, 5.830d07, 0d0, 0d0, 0d0, 1.210d12, &
          6.800d16, 1.490d08, 0d0, 0d0, 0d0, 2.690d12, &
          9.200d16, 6.680d08, 0d0, 0d0, 0d0, 1.000d13, &
          1.150d17, 6.030d09, 0d0, 0d0, 0d0, 5.940d13, &
          1.210d17, 2.100d10, 0d0, 0d0, 0d0, 1.520d14, &
          1.240d17, 8.410d10, 0d0, 0d0, 0d0, 4.190d14, &
          1.260d17, 2.820d11, 0d0, 0d0, 0d0, 1.010d15, &
          1.270d17, 8.020d11, 0d0, 0d0, 0d0, 2.150d15, &
          1.280d17, 1.740d12, 0d0, 0d0, 0d0, 3.780d15/), (/6,66/))

      ! Microturbulent velocity (km/s)
      double precision, dimension(66), target:: VmiP= (/ &
         7.954d0,  7.809d0,  7.664d0,  7.510d0,  7.340d0,  7.158d0, &
         6.957d0,  6.744d0,  6.619d0,  6.490d0,  6.355d0,  6.200d0, &
         6.023d0,  5.927d0,  5.825d0,  5.734d0,  5.635d0,  5.530d0, &
         5.410d0,  5.343d0,  5.299d0,  5.254d0,  5.200d0,  5.146d0, &
         5.120d0,  5.092d0,  5.062d0,  5.027d0,  5.007d0,  4.964d0, &
         4.931d0,  4.903d0,  4.878d0,  4.831d0,  4.765d0,  4.654d0, &
         4.416d0,  4.129d0,  3.847d0,  3.528d0,  3.279d0,  2.703d0, &
         2.335d0,  2.011d0,  1.793d0,  1.346d0,  1.137d0,  0.942d0, &
        0.7776d0, 0.6448d0,  0.549d0, 0.4927d0, 0.4778d0, 0.5079d0, &
        0.6325d0, 0.8584d0,  1.137d0,  1.435d0,  1.720d0,  1.944d0, &
         2.000d0,  2.000d0,  2.000d0,  2.000d0,  2.000d0,  2.000d0  /)
      ! Microturbulent velocity (scaled)
     !double precision, dimension(66), target:: VmiP= (/ &
     !  2.653d-05, 2.605d-05, 2.556d-05, 2.505d-05, 2.448d-05, &
     !  2.388d-05, 2.321d-05, 2.250d-05, 2.208d-05, 2.165d-05, &
     !  2.120d-05, 2.068d-05, 2.009d-05, 1.977d-05, 1.943d-05, &
     !  1.913d-05, 1.880d-05, 1.845d-05, 1.805d-05, 1.782d-05, &
     !  1.768d-05, 1.753d-05, 1.735d-05, 1.717d-05, 1.708d-05, &
     !  1.699d-05, 1.689d-05, 1.677d-05, 1.670d-05, 1.656d-05, &
     !  1.645d-05, 1.635d-05, 1.627d-05, 1.611d-05, 1.589d-05, &
     !  1.552d-05, 1.473d-05, 1.377d-05, 1.283d-05, 1.177d-05, &
     !  1.094d-05, 9.016d-06, 7.789d-06, 6.708d-06, 5.981d-06, &
     !  4.490d-06, 3.793d-06, 3.142d-06, 2.594d-06, 2.151d-06, &
     !  1.831d-06, 1.643d-06, 1.594d-06, 1.694d-06, 2.110d-06, &
     !  2.863d-06, 3.793d-06, 4.787d-06, 5.737d-06, 6.484d-06, &
     !  6.671d-06, 6.671d-06, 6.671d-06, 6.671d-06, 6.671d-06, &
     !  6.671d-06  /)
      !
      !Fontenla et al. (1993), ApJ, 406, 319
      !************** FAL-P atmosphere model
      !

      !
      !************** FAL-A atmosphere model
      !Fontenla et al. (1993), ApJ, 406, 319
      !

      ! Geometrical height (km)
      double precision, dimension(70), target:: HGHA= (/ &
      2276.73d0, 2273.96d0, 2271.36d0, 2268.95d0, 2265.13d0, &
      2261.33d0, 2259.42d0, 2257.48d0, 2255.52d0, 2253.56d0, &
      2252.58d0, 2251.65d0, 2250.73d0, 2249.79d0, 2248.90d0, &
      2248.05d0, 2247.21d0, 2246.33d0, 2245.61d0, 2244.89d0, &
      2243.80d0, 2242.54d0, 2239.42d0, 2235.00d0, 2200.00d0, &
      2167.00d0, 2149.00d0, 2126.00d0, 2108.00d0, 2074.00d0, &
      2050.00d0, 2022.00d0, 1988.00d0, 1960.00d0, 1910.00d0, &
      1870.00d0, 1785.00d0, 1675.00d0, 1575.00d0, 1475.00d0, &
      1375.00d0, 1275.00d0, 1175.00d0, 1075.00d0,  975.00d0, &
       900.00d0,  850.00d0,  800.00d0,  750.00d0,  700.00d0, &
       650.00d0,  600.00d0,  550.00d0,  525.00d0,  500.00d0, &
       450.00d0,  400.00d0,  350.00d0,  300.00d0,  250.00d0, &
       200.00d0,  150.00d0,  100.00d0,   50.00d0,    0.00d0, &
       -20.00d0,  -40.00d0,  -60.00d0,  -80.00d0, -100.00d0 /)

      ! Temperature
      double precision, dimension(70), target:: TEMPA= (/ &
      102770d0, 98790d0, 94800d0, 90816d0, 83891d0, 75934d0, 71336d0, &
       66145d0, 60170d0, 53284d0, 49385d0, 45416d0, 41178d0, 36594d0, &
       32145d0, 27972d0, 24056d0, 20416d0, 17925d0, 15941d0, 13897d0, &
       12650d0, 11650d0, 11100d0, 10250d0,  9750d0,  9500d0,  9200d0, &
        9000d0,  8600d0,  8300d0,  7950d0,  7650d0,  7450d0,  7180d0, &
        7060d0,  6850d0,  6650d0,  6480d0,  6320d0,  6150d0,  5990d0, &
        5820d0,  5660d0,  5480d0,  5300d0,  5150d0,  4910d0,  4680d0, &
        4510d0,  4390d0,  4280d0,  4240d0,  4240d0,  4260d0,  4320d0, &
        4420d0,  4560d0,  4680d0,  4820d0,  4970d0,  5150d0,  5410d0, &
        5790d0,  6520d0,  6980d0,  7590d0,  8220d0,  8860d0,  9400d0 /)


      ! Electron density
      double precision, dimension(70), target:: NeA= (/ &
          3.320d09, 3.440d09, 3.580d09, 3.720d09, 4.000d09, &
          4.370d09, 4.620d09, 4.940d09, 5.380d09, 6.000d09, &
          6.430d09, 6.920d09, 7.550d09, 8.380d09, 9.370d09, &
          1.050d10, 1.200d10, 1.370d10, 1.510d10, 1.650d10, &
          1.810d10, 1.900d10, 1.880d10, 1.760d10, 1.650d10, &
          1.690d10, 1.710d10, 1.720d10, 1.730d10, 1.730d10, &
          1.730d10, 1.730d10, 1.740d10, 1.740d10, 1.780d10, &
          1.830d10, 1.970d10, 2.250d10, 2.550d10, 2.840d10, &
          3.100d10, 3.290d10, 3.470d10, 3.620d10, 3.750d10, &
          3.640d10, 3.470d10, 3.200d10, 3.440d10, 4.680d10, &
          7.140d10, 1.130d11, 1.800d11, 2.260d11, 2.840d11, &
          4.470d11, 7.000d11, 1.090d12, 1.690d12, 2.600d12, &
          3.970d12, 6.130d12, 9.960d12, 1.970d13, 7.670d13, &
          1.730d14, 4.470d14, 1.050d15, 2.210d15, 3.860d15 /)

      ! Hydrogen density
      double precision, dimension(6,70), target:: NHA= reshape((/ &
          6.050d04,9.150d-03, 0d0, 0d0, 0d0, 2.770d09, &
          1.080d05,1.600d-02, 0d0, 0d0, 0d0, 2.890d09, &
          1.680d05,2.460d-02, 0d0, 0d0, 0d0, 3.010d09, &
          2.700d05,3.930d-02, 0d0, 0d0, 0d0, 3.140d09, &
          6.650d05,9.570d-02, 0d0, 0d0, 0d0, 3.410d09, &
          2.110d06,3.000d-01, 0d0, 0d0, 0d0, 3.770d09, &
          4.380d06,6.240d-01, 0d0, 0d0, 0d0, 4.020d09, &
          9.240d06, 1.330d00, 0d0, 0d0, 0d0, 4.330d09, &
          2.020d07, 3.010d00, 0d0, 0d0, 0d0, 4.750d09, &
          4.550d07, 7.280d00, 0d0, 0d0, 0d0, 5.330d09, &
          7.040d07, 1.190d01, 0d0, 0d0, 0d0, 5.720d09, &
          1.070d08, 1.920d01, 0d0, 0d0, 0d0, 6.170d09, &
          1.620d08, 3.130d01, 0d0, 0d0, 0d0, 6.740d09, &
          2.470d08, 5.220d01, 0d0, 0d0, 0d0, 7.480d09, &
          3.690d08, 8.480d01, 0d0, 0d0, 0d0, 8.380d09, &
          5.370d08, 1.320d02, 0d0, 0d0, 0d0, 9.450d09, &
          7.670d08, 1.990d02, 0d0, 0d0, 0d0, 1.070d10, &
          1.090d09, 2.950d02, 0d0, 0d0, 0d0, 1.230d10, &
          1.430d09, 3.930d02, 0d0, 0d0, 0d0, 1.370d10, &
          1.840d09, 5.090d02, 0d0, 0d0, 0d0, 1.510d10, &
          2.590d09, 7.170d02, 0d0, 0d0, 0d0, 1.680d10, &
          3.590d09, 9.970d02, 0d0, 0d0, 0d0, 1.780d10, &
          6.180d09, 1.740d03, 0d0, 0d0, 0d0, 1.780d10, &
          9.660d09, 2.740d03, 0d0, 0d0, 0d0, 1.680d10, &
          1.710d10, 4.700d03, 0d0, 0d0, 0d0, 1.550d10, &
          2.170d10, 5.180d03, 0d0, 0d0, 0d0, 1.580d10, &
          2.460d10, 5.270d03, 0d0, 0d0, 0d0, 1.590d10, &
          2.890d10, 5.320d03, 0d0, 0d0, 0d0, 1.600d10, &
          3.260d10, 5.330d03, 0d0, 0d0, 0d0, 1.600d10, &
          4.100d10, 5.400d03, 0d0, 0d0, 0d0, 1.590d10, &
          4.820d10, 5.480d03, 0d0, 0d0, 0d0, 1.590d10, &
          5.820d10, 5.680d03, 0d0, 0d0, 0d0, 1.590d10, &
          7.160d10, 5.980d03, 0d0, 0d0, 0d0, 1.580d10, &
          8.400d10, 6.300d03, 0d0, 0d0, 0d0, 1.590d10, &
          1.100d11, 7.000d03, 0d0, 0d0, 0d0, 1.630d10, &
          1.350d11, 7.670d03, 0d0, 0d0, 0d0, 1.680d10, &
          2.080d11, 9.570d03, 0d0, 0d0, 0d0, 1.840d10, &
          3.660d11, 1.320d04, 0d0, 0d0, 0d0, 2.130d10, &
          6.250d11, 1.750d04, 0d0, 0d0, 0d0, 2.440d10, &
          1.090d12, 2.220d04, 0d0, 0d0, 0d0, 2.740d10, &
          1.960d12, 2.680d04, 0d0, 0d0, 0d0, 2.990d10, &
          3.640d12, 3.090d04, 0d0, 0d0, 0d0, 3.170d10, &
          7.060d12, 3.440d04, 0d0, 0d0, 0d0, 3.300d10, &
          1.430d13, 3.690d04, 0d0, 0d0, 0d0, 3.340d10, &
          3.010d13, 3.780d04, 0d0, 0d0, 0d0, 3.260d10, &
          5.460d13, 3.240d04, 0d0, 0d0, 0d0, 2.830d10, &
          8.300d13, 2.540d04, 0d0, 0d0, 0d0, 2.290d10, &
          1.310d14, 1.480d04, 0d0, 0d0, 0d0, 1.390d10, &
          2.120d14, 7.840d03, 0d0, 0d0, 0d0, 6.620d09, &
          3.450d14, 5.100d03, 0d0, 0d0, 0d0, 3.070d09, &
          5.660d14, 4.390d03, 0d0, 0d0, 0d0, 1.670d09, &
          9.430d14, 4.140d03, 0d0, 0d0, 0d0, 9.510d08, &
          1.560d15, 5.560d03, 0d0, 0d0, 0d0, 7.700d08, &
          2.000d15, 7.140d03, 0d0, 0d0, 0d0, 7.700d08, &
          2.550d15, 1.010d04, 0d0, 0d0, 0d0, 8.550d08, &
          4.110d15, 2.280d04, 0d0, 0d0, 0d0, 1.190d09, &
          6.510d15, 6.410d04, 0d0, 0d0, 0d0, 2.130d09, &
          1.010d16, 2.190d05, 0d0, 0d0, 0d0, 4.890d09, &
          1.550d16, 6.500d05, 0d0, 0d0, 0d0, 1.050d10, &
          2.340d16, 2.030d06, 0d0, 0d0, 0d0, 2.760d10, &
          3.470d16, 6.320d06, 0d0, 0d0, 0d0, 8.470d10, &
          5.030d16, 2.110d07, 0d0, 0d0, 0d0, 3.040d11, &
          7.070d16, 8.930d07, 0d0, 0d0, 0d0, 1.390d12, &
          9.510d16, 5.050d08, 0d0, 0d0, 0d0, 7.470d12, &
          1.180d17, 6.160d09, 0d0, 0d0, 0d0, 5.990d13, &
          1.240d17, 2.140d10, 0d0, 0d0, 0d0, 1.530d14, &
          1.270d17, 8.590d10, 0d0, 0d0, 0d0, 4.240d14, &
          1.290d17, 2.890d11, 0d0, 0d0, 0d0, 1.020d15, &
          1.300d17, 8.200d11, 0d0, 0d0, 0d0, 2.180d15, &
          1.310d17, 1.780d12, 0d0, 0d0, 0d0, 3.820d15 /), (/6,70/))

      ! Microturbulent velocity (km/s)
      double precision, dimension(70), target:: VmiA= (/ &
        13.10d0, 13.10d0, 13.00d0, 12.90d0, 12.80d0, 12.60d0, &
        12.40d0, 12.30d0, 12.10d0, 11.90d0, 11.80d0, 11.60d0, &
        11.40d0, 11.20d0, 11.00d0, 10.70d0, 10.50d0, 10.20d0, &
        9.960d0, 9.760d0, 9.520d0, 9.350d0, 9.150d0, 8.980d0, &
        8.630d0, 8.410d0, 8.290d0, 8.120d0, 8.000d0, 7.770d0, &
        7.600d0, 7.390d0, 7.170d0, 7.000d0, 6.720d0, 6.510d0, &
        6.090d0, 5.600d0, 5.180d0, 4.800d0, 4.440d0, 4.080d0, &
        3.700d0, 3.280d0, 2.810d0, 2.430d0, 2.170d0, 1.890d0, &
        1.610d0, 1.350d0, 1.100d0, 0.884d0, 0.708d0, 0.638d0, &
        0.581d0, 0.504d0, 0.477d0, 0.501d0, 0.630d0, 0.870d0, &
        1.170d0, 1.470d0, 1.750d0, 1.970d0, 2.000d0, 2.000d0, &
        2.000d0, 2.000d0, 2.000d0, 2.000d0 /)

      !
      !Fontenla et al. (1993), ApJ, 406, 319
      !************** FAL-A atmosphere model
      !

      !
      !************** FAL-X (MCO) atmosphere model
      !Ayres (1986), Highlights of Astronomy, 7, 425
      !

      ! Geometrical height (km)
      double precision, dimension(80), target:: HGHX= (/ &
      2167.302d0, 2165.582d0, 2163.952d0, 2161.808d0, 2159.664d0, &
      2158.580d0, 2157.479d0, 2156.360d0, 2155.245d0, 2154.690d0, &
      2154.168d0, 2153.648d0, 2153.118d0, 2152.620d0, 2152.151d0, &
      2151.686d0, 2151.200d0, 2150.802d0, 2150.521d0, 2150.109d0, &
      2149.770d0, 2149.012d0, 2147.895d0, 2138.818d0, 2116.624d0, &
      2088.354d0, 2058.080d0, 2034.852d0, 2022.723d0, 2009.583d0, &
      1990.388d0, 1964.138d0, 1926.844d0, 1861.745d0, 1806.527d0, &
      1721.992d0, 1618.715d0, 1529.763d0, 1427.505d0, 1333.292d0, &
      1237.594d0, 1145.563d0, 1045.505d0,  967.258d0,  896.738d0, &
       847.862d0,  799.010d0,  750.104d0,  701.129d0,  647.515d0, &
       598.423d0,  558.732d0,  523.850d0,  488.955d0,  449.201d0, &
       399.889d0,  350.530d0,  301.127d0,  251.680d0,  202.190d0, &
       177.067d0,  151.940d0,  126.803d0,  101.655d0,   76.485d0, &
        51.298d0,   36.166d0,   21.035d0,   10.945d0,    0.855d0, &
        -9.235d0,  -19.325d0,  -29.415d0,  -39.505d0,  -49.594d0, &
       -59.683d0,  -69.772d0,  -79.861d0,  -89.950d0, -100.039d0 /)

      ! Temperature
      double precision, dimension(80), target:: TEMPX= (/ &
     100000d0, 95600d0, 90816d0, 83891d0, 75934d0, 71336d0, 66145d0, &
      60170d0, 53284d0, 49385d0, 45416d0, 41178d0, 36594d0, 32145d0, &
      27972d0, 24056d0, 20416d0, 17925d0, 16500d0, 15000d0, 14250d0, &
      13500d0, 13000d0, 12000d0, 11150d0, 10550d0,  9900d0,  9450d0, &
       9200d0,  8950d0,  8700d0,  8400d0,  8050d0,  7650d0,  7450d0, &
       7250d0,  7050d0,  6900d0,  6550d0,  6100d0,  5500d0,  4800d0, &
       4000d0,  3850d0,  3890d0,  3930d0,  3970d0,  4040d0,  4120d0, &
       4218d0,  4307d0,  4379d0,  4442d0,  4505d0,  4577d0,  4666d0, &
       4756d0,  4846d0,  4936d0,  5025d0,  5070d0,  5115d0,  5180d0, &
       5280d0,  5480d0,  5790d0,  5980d0,  6180d0,  6340d0,  6520d0, &
       6720d0,  6980d0,  7280d0,  7590d0,  7900d0,  8220d0,  8540d0, &
       8860d0,  9140d0,  9400d0 /)

      ! Electron density
      double precision, dimension(80), target:: NeX= (/ &
          0.2689d10, 0.2803d10, 0.2939d10, 0.3162d10, 0.3465d10, &
          0.3669d10, 0.3932d10, 0.4288d10, 0.4790d10, 0.5133d10, &
          0.5538d10, 0.6051d10, 0.6726d10, 0.7548d10, 0.8529d10, &
          0.9718d10, 0.1116d11, 0.1240d11, 0.1320d11, 0.1408d11, &
          0.1458d11, 0.1506d11, 0.1524d11, 0.1464d11, 0.1466d11, &
          0.1544d11, 0.1633d11, 0.1695d11, 0.1731d11, 0.1769d11, &
          0.1813d11, 0.1877d11, 0.1967d11, 0.2122d11, 0.2245d11, &
          0.2409d11, 0.2575d11, 0.2704d11, 0.2759d11, 0.2723d11, &
          0.2611d11, 0.2462d11, 0.2224d11, 0.2018d11, 0.2229d11, &
          0.2124d11, 0.2402d11, 0.3285d11, 0.5061d11, 0.8350d11, &
          0.1304d12, 0.1856d12, 0.2522d12, 0.3419d12, 0.4831d12, &
          0.7407d12, 0.1128d13, 0.1710d13, 0.2581d13, 0.3889d13, &
          0.4781d13, 0.5880d13, 0.7285d13, 0.9160d13, 0.1230d14, &
          0.1964d14, 0.2770d14, 0.4058d14, 0.5498d14, 0.7700d14, &
          0.1109d15, 0.1734d15, 0.2815d15, 0.4494d15, 0.6948d15, &
          0.1054d16, 0.1553d16, 0.2225d16, 0.2994d16, 0.3886d16 /)

      ! Hydrogen density
      double precision, dimension(6,80), target:: NHX= reshape((/ &
   0.4090d6 ,0.4311d-1,0.2039d-1,0.2295d-1,0.2784d-1,0.2243d10, &
   0.9762d6 , 0.1022d0,0.4693d-1,0.5186d-1,0.6232d-1,0.2346d10, &
   0.1893d7 , 0.1975d0,0.8798d-1,0.9578d-1, 0.1144d0,0.2468d10, &
   0.4264d7 , 0.4446d0, 0.1877d0, 0.2000d0, 0.2374d0,0.2667d10, &
   0.9777d7 , 0.1025d1, 0.3975d0, 0.4124d0, 0.4858d0,0.2938d10, &
   0.1532d8 , 0.1619d1, 0.5885d0, 0.6000d0, 0.7034d0,0.3118d10, &
   0.2424d8 , 0.2592d1, 0.8616d0, 0.8596d0, 0.1001d1,0.3349d10, &
   0.3891d8 , 0.4228d1, 0.1238d1, 0.1200d1, 0.1386d1,0.3658d10, &
   0.6315d8 , 0.6987d1, 0.1698d1, 0.1584d1, 0.1807d1,0.4091d10, &
   0.8099d8 , 0.9026d1, 0.1930d1, 0.1754d1, 0.1984d1,0.4384d10, &
   0.1028d9 , 0.1151d2, 0.2114d1, 0.1864d1, 0.2086d1,0.4731d10, &
   0.1310d9 , 0.1467d2, 0.2225d1, 0.1887d1, 0.2081d1,0.5169d10, &
   0.1688d9 , 0.1883d2, 0.2219d1, 0.1784d1, 0.1928d1,0.5748d10, &
   0.2155d9 , 0.2390d2, 0.2096d1, 0.1563d1, 0.1645d1,0.6455d10, &
   0.2729d9 , 0.3003d2, 0.1888d1, 0.1283d1, 0.1304d1,0.7305d10, &
   0.3465d9 , 0.3788d2, 0.1669d1, 0.1019d1, 0.9976d0,0.8346d10, &
   0.4466d9 , 0.4872d2, 0.1528d1, 0.8532d0, 0.8193d0,0.9634d10, &
   0.5510d9 , 0.6028d2, 0.1508d1, 0.8302d0, 0.8123d0,0.1080d11, &
   0.6385d9 , 0.7021d2, 0.1549d1, 0.8675d0, 0.8712d0,0.1162d11, &
   0.7846d9 , 0.8727d2, 0.1660d1, 0.9667d0, 0.1006d1,0.1266d11, &
   0.9157d9 , 0.1031d3, 0.1768d1, 0.1058d1, 0.1121d1,0.1320d11, &
   0.1221d10, 0.1420d3, 0.1988d1, 0.1230d1, 0.1320d1,0.1369d11, &
   0.1684d10, 0.2064d3, 0.2300d1, 0.1436d1, 0.1549d1,0.1392d11, &
   0.4786d10, 0.9001d3, 0.6016d1, 0.2803d1, 0.2776d1,0.1328d11, &
   0.7671d10, 0.2339d4, 0.1422d2, 0.5814d1, 0.5138d1,0.1313d11, &
   0.9535d10, 0.2996d4, 0.1793d2, 0.7253d1, 0.6379d1,0.1378d11, &
   0.1206d11, 0.3348d4, 0.1982d2, 0.8044d1, 0.7149d1,0.1450d11, &
   0.1444d11, 0.3585d4, 0.2107d2, 0.8579d1, 0.7690d1,0.1499d11, &
   0.1588d11, 0.3721d4, 0.2178d2, 0.8888d1, 0.8008d1,0.1527d11, &
   0.1756d11, 0.3874d4, 0.2259d2, 0.9236d1, 0.8365d1,0.1556d11, &
   0.2005d11, 0.4090d4, 0.2378d2, 0.9723d1, 0.8845d1,0.1589d11, &
   0.2378d11, 0.4435d4, 0.2570d2, 0.1050d2, 0.9595d1,0.1635d11, &
   0.2991d11, 0.5025d4, 0.2903d2, 0.1180d2, 0.1083d2,0.1701d11, &
   0.4353d11, 0.6263d4, 0.3622d2, 0.1451d2, 0.1331d2,0.1814d11, &
   0.5874d11, 0.7452d4, 0.4338d2, 0.1710d2, 0.1562d2,0.1906d11, &
   0.9098d11, 0.9460d4, 0.5616d2, 0.2150d2, 0.1948d2,0.2048d11, &
   0.1553d12, 0.1226d5, 0.7543d2, 0.2777d2, 0.2484d2,0.2236d11, &
   0.2465d12, 0.1487d5, 0.9488d2, 0.3378d2, 0.2992d2,0.2405d11, &
   0.4408d12, 0.1767d5, 0.1184d3, 0.4053d2, 0.3562d2,0.2553d11, &
   0.7931d12, 0.1940d5, 0.1363d3, 0.4509d2, 0.3944d2,0.2592d11, &
   0.1563d13, 0.2044d5, 0.1502d3, 0.4823d2, 0.4186d2,0.2521d11, &
   0.3409d13, 0.2112d5, 0.1608d3, 0.5011d2, 0.4339d2,0.2367d11, &
   0.1047d14, 0.1995d5, 0.1595d3, 0.4803d2, 0.4146d2,0.2033d11, &
   0.2504d14, 0.1508d5, 0.1237d3, 0.3630d2, 0.3141d2,0.1656d11, &
   0.5359d14, 0.1488d5, 0.1251d3, 0.3647d2, 0.3143d2,0.1510d11, &
   0.8892d14, 0.8806d4, 0.7538d2, 0.2181d2, 0.1873d2,0.9495d10, &
   0.1481d15, 0.5318d4, 0.4563d2, 0.1321d2, 0.1133d2,0.5046d10, &
   0.2443d15, 0.2955d4, 0.2552d2, 0.7400d1, 0.6348d1,0.2051d10, &
   0.3986d15, 0.2262d4, 0.1956d2, 0.5703d1, 0.4909d1,0.1016d10, &
   0.6759d15, 0.2790d4, 0.2417d2, 0.7069d1, 0.6096d1, 0.7410d9, &
   0.1081d16, 0.5728d4, 0.4985d2, 0.1463d2, 0.1263d2, 0.9565d9, &
   0.1560d16, 0.1214d5, 0.1064d3, 0.3137d2, 0.2710d2,0.1420d10, &
   0.2139d16, 0.2378d5, 0.2112d3, 0.6248d2, 0.5388d2,0.2052d10, &
   0.2920d16, 0.4648d5, 0.4231d3, 0.1254d3, 0.1074d3,0.2983d10, &
   0.4146d16, 0.9903d5, 0.9510d3, 0.2811d3, 0.2378d3,0.4622d10, &
   0.6379d16, 0.2485d6, 0.2702d4, 0.7913d3, 0.6494d3,0.8141d10, &
   0.9734d16, 0.6115d6, 0.8141d4, 0.2376d4, 0.1867d4,0.1523d11, &
   0.1474d17, 0.1467d7, 0.2542d5, 0.7656d4, 0.5785d4,0.3120d11, &
   0.2215d17, 0.3436d7, 0.7697d5, 0.2515d5, 0.1870d5,0.6805d11, &
   0.3305d17, 0.7833d7, 0.2119d6, 0.7670d5, 0.5782d5,0.1451d12, &
   0.4026d17, 0.1176d8, 0.3458d6, 0.1305d6, 0.9963d5,0.2082d12, &
   0.4897d17, 0.1757d8, 0.5552d6, 0.2167d6, 0.1676d6,0.2920d12, &
   0.5921d17, 0.2840d8, 0.9659d6, 0.3888d6, 0.3048d6,0.4443d12, &
   0.7092d17, 0.5242d8, 0.1949d7, 0.8121d6, 0.6472d6,0.7907d12, &
   0.8302d17, 0.1389d9, 0.5913d7, 0.2596d7, 0.2121d7,0.2131d13, &
   0.9463d17, 0.5026d9, 0.2592d8, 0.1223d8, 0.1034d8,0.7510d13, &
   0.1019d18,0.1036d10, 0.6002d8, 0.2955d8, 0.2548d8,0.1426d14, &
   0.1094d18,0.2109d10, 0.1373d9, 0.7050d8, 0.6198d8,0.2576d14, &
   0.1140d18,0.3563d10, 0.2532d9, 0.1342d9, 0.1197d9,0.3917d14, &
   0.1183d18,0.6189d10, 0.4833d9, 0.2648d9, 0.2400d9,0.6010d14, &
   0.1222d18,0.1098d11, 0.9470d9, 0.5375d9, 0.4951d9,0.9275d14, &
   0.1250d18,0.2164d11,0.2106d10,0.1247d10,0.1172d10,0.1539d15, &
   0.1270d18,0.4419d11,0.4896d10,0.3034d10,0.2911d10,0.2603d15, &
   0.1286d18,0.8692d11,0.1089d11,0.7049d10,0.6899d10,0.4262d15, &
   0.1300d18,0.1620d12,0.2274d11,0.1531d11,0.1526d11,0.6693d15, &
   0.1309d18,0.2924d12,0.4573d11,0.3198d11,0.3244d11,0.1026d16, &
   0.1315d18,0.5039d12,0.8709d11,0.6306d11,0.6502d11,0.1522d16, &
   0.1317d18,0.8327d12,0.1579d12,0.1181d12,0.1236d12,0.2192d16, &
   0.1323d18,0.1259d13,0.2576d12,0.1978d12,0.2096d12,0.2958d16, &
   0.1329d18,0.1810d13,0.3955d12,0.3109d12,0.3330d12,0.3848d16 /), &
                                                             (/6,80/))

      ! Microturbulent velocity (km/s)
      double precision, dimension(80), target:: VmiX= (/ &
        13.459d0, 13.372d0, 13.274d0, 13.123d0, 12.934d0, &
        12.817d0, 12.675d0, 12.499d0, 12.275d0, 12.136d0, &
        11.984d0, 11.821d0, 11.635d0, 11.433d0, 11.218d0, &
        10.986d0, 10.737d0, 10.555d0, 10.436d0, 10.293d0, &
        10.217d0, 10.133d0, 10.063d0,  9.837d0,  9.627d0, &
         9.473d0,  9.296d0,  9.157d0,  9.080d0,  8.997d0, &
         8.887d0,  8.738d0,  8.526d0,  8.185d0,  7.901d0, &
         7.464d0,  6.853d0,  6.438d0,  5.952d0,  5.479d0, &
         4.993d0,  4.408d0,  3.499d0,  2.846d0,  2.347d0, &
         2.053d0,  1.774d0,  1.512d0,  1.337d0,  1.138d0, &
         0.976d0,  0.874d0,  0.795d0,  0.717d0,  0.653d0, &
         0.557d0,  0.523d0,  0.568d0,  0.622d0,  0.785d0, &
         0.876d0,  0.980d0,  1.087d0,  1.199d0,  1.297d0, &
         1.395d0,  1.451d0,  1.513d0,  1.550d0,  1.596d0, &
         1.641d0,  1.681d0,  1.709d0,  1.738d0,  1.766d0, &
         1.792d0,  1.807d0,  1.816d0,  1.830d0,  1.830d0 /)

      !
      !Ayres (1986), Highlights of Astronomy, 7, 425
      !************** FAL-X (MCO) atmosphere model
      !

      !
      !************** FAL-F atmosphere model
      !Fontenla et al. (1993), ApJ, 406, 319
      !

      ! Geometrical height (km)
      double precision, dimension(70), target:: HGHF= (/ &
      2021.54d0, 2021.02d0, 2020.54d0, 2020.09d0, 2019.37d0, &
      2018.66d0, 2018.31d0, 2017.94d0, 2017.58d0, 2017.21d0, &
      2017.03d0, 2016.86d0, 2016.69d0, 2016.52d0, 2016.36d0, &
      2016.21d0, 2016.06d0, 2015.90d0, 2015.78d0, 2015.66d0, &
      2015.47d0, 2015.26d0, 2014.73d0, 2014.00d0, 2008.00d0, &
      1995.00d0, 1977.00d0, 1955.00d0, 1935.00d0, 1923.00d0, &
      1900.00d0, 1880.00d0, 1854.00d0, 1818.00d0, 1790.00d0, &
      1760.00d0, 1730.00d0, 1680.00d0, 1590.00d0, 1475.00d0, &
      1375.00d0, 1275.00d0, 1175.00d0, 1075.00d0,  975.00d0, &
       900.00d0,  850.00d0,  800.00d0,  750.00d0,  700.00d0, &
       650.00d0,  600.00d0,  550.00d0,  525.00d0,  500.00d0, &
       450.00d0,  400.00d0,  350.00d0,  300.00d0,  250.00d0, &
       200.00d0,  150.00d0,  100.00d0,   50.00d0,    0.00d0, &
       -20.00d0,  -40.00d0,  -60.00d0,  -80.00d0, -100.00d0 /)

      ! Temperature
      double precision, dimension(70), target:: TEMPF= (/ &
     102770d0, 98790d0, 94800d0, 90816d0, 83891d0, 75934d0, 71336d0, &
      66145d0, 60170d0, 53284d0, 49385d0, 45416d0, 41178d0, 36594d0, &
      32145d0, 27972d0, 24056d0, 20416d0, 17925d0, 15941d0, 13897d0, &
      12650d0, 11650d0, 11100d0, 10500d0, 10000d0,  9700d0,  9400d0, &
       9200d0,  9080d0,  8900d0,  8750d0,  8550d0,  8300d0,  8100d0, &
       7900d0,  7750d0,  7550d0,  7340d0,  7120d0,  6930d0,  6750d0, &
       6550d0,  6360d0,  6180d0,  6000d0,  5870d0,  5730d0,  5540d0, &
       5300d0,  5020d0,  4870d0,  4760d0,  4720d0,  4700d0,  4700d0, &
       4740d0,  4800d0,  4880d0,  4970d0,  5100d0,  5250d0,  5450d0, &
       5790d0,  6520d0,  6980d0,  7590d0,  8220d0,  8860d0,  9400d0 /)

      ! Electron density
      double precision, dimension(70), target:: NeF= (/ &
          1.880d10, 1.950d10, 2.020d10, 2.100d10, 2.260d10, &
          2.480d10, 2.620d10, 2.810d10, 3.060d10, 3.420d10, &
          3.670d10, 3.960d10, 4.330d10, 4.820d10, 5.420d10, &
          6.140d10, 7.010d10, 8.090d10, 9.030d10, 9.940d10, &
          1.100d11, 1.170d11, 1.180d11, 1.130d11, 1.100d11, &
          1.150d11, 1.190d11, 1.220d11, 1.240d11, 1.250d11, &
          1.270d11, 1.290d11, 1.310d11, 1.340d11, 1.370d11, &
          1.390d11, 1.430d11, 1.490d11, 1.590d11, 1.750d11, &
          1.900d11, 2.080d11, 2.240d11, 2.390d11, 2.600d11, &
          2.600d11, 2.540d11, 2.430d11, 2.070d11, 1.720d11, &
          1.500d11, 1.780d11, 2.430d11, 2.930d11, 3.570d11, &
          5.390d11, 8.170d11, 1.240d12, 1.860d12, 2.800d12, &
          4.250d12, 6.480d12, 1.030d13, 1.970d13, 7.670d13, &
          1.730d14, 4.470d14, 1.050d15, 2.210d15, 3.860d15 /)

      ! Hydrogen density
      double precision, dimension(6,70), target:: NHF= reshape((/ &
          3.270d05,3.220d-01, 0d0, 0d0, 0d0, 1.570d10, &
          5.540d05,5.370d-01, 0d0, 0d0, 0d0, 1.630d10, &
          8.100d05,7.770d-01, 0d0, 0d0, 0d0, 1.700d10, &
          1.230d06, 1.170d00, 0d0, 0d0, 0d0, 1.780d10, &
          2.880d06, 2.680d00, 0d0, 0d0, 0d0, 1.930d10, &
          9.110d06, 8.340d00, 0d0, 0d0, 0d0, 2.140d10, &
          1.920d07, 1.750d01, 0d0, 0d0, 0d0, 2.280d10, &
          4.240d07, 3.900d01, 0d0, 0d0, 0d0, 2.460d10, &
          9.540d07, 9.010d01, 0d0, 0d0, 0d0, 2.700d10, &
          2.270d08, 2.290d02, 0d0, 0d0, 0d0, 3.040d10, &
          3.590d08, 3.820d02, 0d0, 0d0, 0d0, 3.270d10, &
          5.540d08, 6.310d02, 0d0, 0d0, 0d0, 3.540d10, &
          8.590d08, 1.060d03, 0d0, 0d0, 0d0, 3.870d10, &
          1.340d09, 1.820d03, 0d0, 0d0, 0d0, 4.320d10, &
          2.030d09, 3.040d03, 0d0, 0d0, 0d0, 4.860d10, &
          2.970d09, 4.840d03, 0d0, 0d0, 0d0, 5.510d10, &
          4.300d09, 7.520d03, 0d0, 0d0, 0d0, 6.310d10, &
          6.260d09, 1.150d04, 0d0, 0d0, 0d0, 7.310d10, &
          8.160d09, 1.550d04, 0d0, 0d0, 0d0, 8.200d10, &
          1.050d10, 2.020d04, 0d0, 0d0, 0d0, 9.090d10, &
          1.490d10, 2.920d04, 0d0, 0d0, 0d0, 1.020d11, &
          2.050d10, 4.090d04, 0d0, 0d0, 0d0, 1.090d11, &
          3.490d10, 7.260d04, 0d0, 0d0, 0d0, 1.100d11, &
          5.300d10, 1.140d05, 0d0, 0d0, 0d0, 1.060d11, &
          7.470d10, 1.880d05, 0d0, 0d0, 0d0, 1.020d11, &
          8.470d10, 2.120d05, 0d0, 0d0, 0d0, 1.070d11, &
          9.810d10, 2.240d05, 0d0, 0d0, 0d0, 1.110d11, &
          1.190d11, 2.370d05, 0d0, 0d0, 0d0, 1.140d11, &
          1.390d11, 2.500d05, 0d0, 0d0, 0d0, 1.160d11, &
          1.520d11, 2.570d05, 0d0, 0d0, 0d0, 1.170d11, &
          1.790d11, 2.720d05, 0d0, 0d0, 0d0, 1.200d11, &
          2.050d11, 2.850d05, 0d0, 0d0, 0d0, 1.230d11, &
          2.450d11, 3.020d05, 0d0, 0d0, 0d0, 1.260d11, &
          3.090d11, 3.250d05, 0d0, 0d0, 0d0, 1.300d11, &
          3.710d11, 3.420d05, 0d0, 0d0, 0d0, 1.330d11, &
          4.500d11, 3.600d05, 0d0, 0d0, 0d0, 1.370d11, &
          5.390d11, 3.780d05, 0d0, 0d0, 0d0, 1.400d11, &
          7.220d11, 4.080d05, 0d0, 0d0, 0d0, 1.460d11, &
          1.190d12, 4.570d05, 0d0, 0d0, 0d0, 1.570d11, &
          2.270d12, 5.290d05, 0d0, 0d0, 0d0, 1.720d11, &
          4.040d12, 6.010d05, 0d0, 0d0, 0d0, 1.860d11, &
          7.350d12, 6.930d05, 0d0, 0d0, 0d0, 2.030d11, &
          1.380d13, 7.680d05, 0d0, 0d0, 0d0, 2.180d11, &
          2.650d13, 8.450d05, 0d0, 0d0, 0d0, 2.320d11, &
          5.240d13, 9.520d05, 0d0, 0d0, 0d0, 2.490d11, &
          8.950d13, 9.170d05, 0d0, 0d0, 0d0, 2.430d11, &
          1.300d14, 8.560d05, 0d0, 0d0, 0d0, 2.310d11, &
          1.900d14, 7.620d05, 0d0, 0d0, 0d0, 2.120d11, &
          2.850d14, 5.500d05, 0d0, 0d0, 0d0, 1.630d11, &
          4.380d14, 3.200d05, 0d0, 0d0, 0d0, 1.090d11, &
          6.940d14, 1.490d05, 0d0, 0d0, 0d0, 5.580d10, &
          1.090d15, 1.160d05, 0d0, 0d0, 0d0, 3.540d10, &
          1.730d15, 1.080d05, 0d0, 0d0, 0d0, 2.360d10, &
          2.180d15, 1.120d05, 0d0, 0d0, 0d0, 2.000d10, &
          2.740d15, 1.280d05, 0d0, 0d0, 0d0, 1.850d10, &
          4.290d15, 2.010d05, 0d0, 0d0, 0d0, 1.920d10, &
          6.650d15, 3.850d05, 0d0, 0d0, 0d0, 2.460d10, &
          1.020d16, 8.070d05, 0d0, 0d0, 0d0, 3.540d10, &
          1.550d16, 1.830d06, 0d0, 0d0, 0d0, 5.740d10, &
          2.330d16, 4.260d06, 0d0, 0d0, 0d0, 1.010d11, &
          3.440d16, 1.150d07, 0d0, 0d0, 0d0, 2.240d11, &
          4.980d16, 3.230d07, 0d0, 0d0, 0d0, 5.640d11, &
          7.030d16, 1.040d08, 0d0, 0d0, 0d0, 1.720d12, &
          9.510d16, 5.050d08, 0d0, 0d0, 0d0, 7.530d12, &
          1.180d17, 6.160d09, 0d0, 0d0, 0d0, 5.990d13, &
          1.240d17, 2.150d10, 0d0, 0d0, 0d0, 1.530d14, &
          1.270d17, 8.600d10, 0d0, 0d0, 0d0, 4.240d14, &
          1.290d17, 2.890d11, 0d0, 0d0, 0d0, 1.020d15, &
          1.300d17, 8.210d11, 0d0, 0d0, 0d0, 2.180d15, &
          1.310d17, 1.780d12, 0d0, 0d0, 0d0, 3.820d15 /), (/6,70/))

      ! Microturbulent velocity (km/s)
      double precision, dimension(70), target:: VmiF= (/ &
        9.890d0, 9.830d0, 9.750d0, 9.670d0, 9.530d0, &
        9.340d0, 9.230d0, 9.100d0, 8.940d0, 8.740d0, &
        8.610d0, 8.480d0, 8.320d0, 8.140d0, 7.940d0, &
        7.740d0, 7.520d0, 7.300d0, 7.130d0, 6.980d0, &
        6.810d0, 6.690d0, 6.560d0, 6.460d0, 6.340d0, &
        6.260d0, 6.170d0, 6.060d0, 5.970d0, 5.920d0, &
        5.820d0, 5.740d0, 5.630d0, 5.490d0, 5.380d0, &
        5.260d0, 5.150d0, 4.970d0, 4.680d0, 4.320d0, &
        4.000d0, 3.670d0, 3.300d0, 2.890d0, 2.460d0, &
        2.120d0, 1.900d0, 1.670d0, 1.450d0, 1.220d0, &
        1.010d0, 0.829d0, 0.678d0, 0.617d0, 0.566d0, &
        0.499d0, 0.477d0, 0.502d0, 0.630d0, 0.867d0, &
        1.160d0, 1.470d0, 1.750d0, 1.970d0, 2.000d0, &
        2.000d0, 2.000d0, 2.000d0, 2.000d0, 2.000d0 /)

      !
      !Fontenla et al. (1993), ApJ, 406, 319
      !************** FAL-F atmosphere model
      !


      ! Routine name
      urou = 'gAtmo'

      ! If FALC
      if (init.eq.0) then

        ! Point to FALC model
        lnz => nzC
       !tau => TAUC
        z => HGHC
        T => TEMPC
        ne => NeC
        nh => NHC
        vmi => VmiC

      ! If FALP
      else if (init.eq.1) then

        ! Point to FALP model
        lnz => nzP
       !tau => TAUP
        z => HGHP
        T => TEMPP
        ne => NeP
        nh => NHP
        vmi => VmiP

      ! If FALA
      else if (init.eq.2) then

        ! Point to FALP model
        lnz => nzA
        z => HGHA
        T => TEMPA
        ne => NeA
        nh => NHA
        vmi => VmiA

      ! If FALX
      else if (init.eq.3) then

        ! Point to FALX (MCO) model
        lnz => nzX
        z => HGHX
        T => TEMPX
        ne => NeX
        nh => NHX
        vmi => VmiX

      ! If FALF
      else if (init.eq.4) then

        ! Point to FALF model
        lnz => nzF
        z => HGHF
        T => TEMPF
        ne => NeF
        nh => NHF
        vmi => VmiF

      end if ! FALC, FALP, FALA, FALX, or FALF

      ! Nullify pointers
      nullify(Atmo%zeros,Atmo%Bx,Atmo%By,Atmo%Bz)

      ! Type of scale, reference frequency, log(g), and number of
      ! nodes
      Atmo%scal = 'H'
      Atmo%tfreq = 0.2d0
      Atmo%logg = 4.44d0
      Atmo%nZ = lnz

      ! Set to geometrical scale
      ztau = .False.

      !
      ! Allocations
      !

      ! Zeros
      allocate(Atmo%zeros(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%zeros)
      Atmo%zeros = 0d0

      ! Height axis
      allocate(Atmo%z(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%z)

      ! Temperature
      allocate(Atmo%T(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%T)

      ! Microturbulence
      allocate(Atmo%vmi(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vmi)

      ! Electron density
      allocate(Atmo%ne(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%ne)

      ! Total hydrogen density
      allocate(Atmo%nht(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nht)
      Atmo%nht = 0d0

      ! Atomic hydrogen density
      allocate(Atmo%nha(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nha)
      Atmo%nha = 0d0

      ! H^- density
      allocate(Atmo%nhm(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhm)
      Atmo%nhm = 0d0

      ! Velocity vector (x,y,z)
      allocate(Atmo%vx(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vx)
      allocate(Atmo%vy(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vy)
      allocate(Atmo%vz(lnZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%vz)

      ! H density (5 HI levels + p+ density)
      allocate(Atmo%nh(lnZ,6))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nh)

      !
      ! Read themodynamics
      !

      ! Copy from pointers
      Atmo%z = z
      Atmo%T = T
      Atmo%ne = ne
      Atmo%vz = 0d0
      Atmo%vmi = vmi
      Atmo%vx = 0d0
      Atmo%vy = 0d0

      ! Set type of density
      Atmo%typo = 0

      !
      ! Read hydrogen densityes
      !

      ! For each height
      do iz=1,lnZ

        ! Store the individual densities
        Atmo%nh(iz,:) = nh(:,iz)

        ! And the total in atomic
        Atmo%nha(iz) = sum(nh(:,iz))

        ! And whole
        Atmo%nht(iz) = Atmo%nha(iz)

      end do ! For each height

      ! Just make a flag to know that there is no helium input
      allocate(Atmo%nhe(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
      Atmo%nhe(1,1) = -1

      ! Free pointers
      nullify(lnz,tau,T,ne,nh,vmi)

      ! Convert to cgs
      Atmo%z = Atmo%z*1d5

      ! Scale micro
      Atmo%vmi = Atmo%vmi*1d-6/c

      ! Allocs
      Atmo%alloc_a = .True.
      Atmo%alloc_b = .True.

      return

      end subroutine gAtmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate an optical depth stratification from parameters\n
      !!     Tau_min(double): Minimum optical depth\n
      !!     Tau_max(double): Maximum optical depth\n
      !!  Strat(double(:,:)): Stratification information to refine the
      !!                      grid\n
      !!         nz(integer): Final number of nodes\n
      !!        z(double(:)): Resulting stratification
      subroutine gAtmo_Strat(Tau_min,Tau_max,Strat,lnz,z)

      ! I/O

      integer, intent(in):: lnz
      double precision, intent(in):: Tau_min,Tau_max
      double precision, dimension(:,:), &
                        allocatable, intent(in):: Strat
      double precision, dimension(:), allocatable, intent(out):: z

      ! Local

      integer:: i,j,k,n,M,nn

      double precision:: Delt,iDelt,d,Delt0,iDelt0


      ! Allocate stratification
      allocate(z(lnz))

      ! If no stratification input
      if (.not.allocated(Strat)) then

        ! Compute the step
        Delt = (Tau_max - Tau_min)/dble(lnz-1)

        ! Initialize first node to minimum
        z(1) = Tau_min

        ! For the rest of nodes
        do i=2,lnz-1

          ! Advance step
          z(i) = z(i-1) + Delt

        end do ! For nodes after first

        ! Last node
        z(lnz) = Tau_max

      ! Stratification input present
      else

        ! Get size of input
        n = size(Strat, dim=2)

        ! First input
        d = (Strat(2,1) - Strat(1,1))*Strat(3,1)

        ! If an independent first part
        if (Tau_min.lt.Strat(1,1)) d = d + Strat(1,1) - Tau_min

        ! Middle
        do i=2,n

          ! Add contribution
          d = d + (Strat(2,i) - Strat(1,i))*Strat(3,i)

          ! If an independent first part
          if (Strat(2,i-1).lt.Strat(1,i)) d = d + Strat(1,i) - &
                                                  Strat(2,i-1)

        end do ! Middle

        ! Last
        if (Strat(2,n).lt.Tau_max) d = d + Tau_max - Strat(2,n)

        ! Get Step
        Delt = d/dble(lnz - 1)
        iDelt = 1d0/Delt

        ! Save step in case we fail
        Delt0 = Delt
        iDelt0 = iDelt

        !
        ! Adjust step
        !
        do j=1,10

          ! Initialize count
          M = 1

          ! If natural first range, add points
          if (Tau_min.lt.Strat(1,1)) M = M + max(1, &
                                                 int((Strat(1,1) - &
                                                      Tau_min)*iDelt))

          ! First range proper, add points
          M = M + max(1,int((Strat(2,1) - &
                             Strat(1,1))*Strat(3,1)*iDelt))

          ! For each range besides the first
          do k=2,n

            ! If natural range before, add points
            if (Strat(2,k-1).lt.Strat(1,k)) M = M + &
                         max(1,int((Strat(1,k) - Strat(2,k-1))*iDelt))

            ! Range proper, add points
            M = M + max(1,int((Strat(2,k) - &
                               Strat(1,k))*Strat(3,k)*iDelt))

          end do ! Ranges but the first

          ! Last natural range
          if (Strat(2,n).lt.Tau_max) M = M + &
                              max(1,int((Tau_max - Strat(2,n))*iDelt))

          ! If correct number of points, leave
          if (M.eq.lnz) exit

          ! If not, adjust Step
          Delt = Delt*dble(M)/dble(lnz)

          iDelt = 1d0/Delt

          ! If exhausted iterations
          if (j.eq.10) then

            ! Recover
            Delt = Delt0
            iDelt = iDelt0

            ! Set extremes
            z(1) = Tau_min
            z(lnz) = Tau_max

            ! Initialize counter
            nn = 1

            ! For intermadiate points
            do i=2,lnz

              ! Advance
              z(i) = z(i-1) + Delt

              ! Beginning of a range
              if ((z(i-1).lt.Strat(1,nn)).and. &
                  (z(i).gt.Strat(1,nn))) then

                ! Adjust
                z(i) = z(i) - &
                            (z(i)-Strat(1,nn))* &
                            (Strat(3,nn)-1)/Strat(3,nn)

                ! Adjust the delta
                Delt = Delt/Strat(3,nn)

              ! End of a range
              else if ((z(i-1).lt.Strat(2,nn)).and. &
                       (z(i).gt.Strat(2,nn))) then

                ! The last range
                if (nn.eq.n) then

                  ! Fix
                  z(i) = z(i) + &
                              (z(i)-Strat(2,nn))* &
                              (Strat(3,nn)-1)

                  ! Adjust the delta
                  Delt = Delt*Strat(3,nn)

                ! Any other range
                else

                  ! Beginning of a range
                  if ((z(i-1).lt.Strat(1,nn+1)).and. &
                      (z(i).gt.Strat(1,nn+1))) then

                    ! Adjust
                    z(i) = z(i) + &
                                (z(i)-Strat(2,nn))* &
                                (Strat(3,nn)-1)
                    z(i) = z(i) - &
                                (z(i)-Strat(1,nn+1))* &
                                (Strat(3,nn+1)-1)/Strat(3,nn+1)

                    ! Adjust the delta
                    Delt = Delt*Strat(3,nn)/Strat(3,nn+1)

                  ! Not Beginning
                  else

                    ! Adjust
                    z(i) = z(i) + &
                                (z(i)-Strat(2,nn))* &
                                (Strat(3,nn)-1)

                    ! Adjust the delta
                    Delt = Delt*Strat(3,nn)

                  end if ! Beginning of range

                  ! move to next range
                  nn = nn + 1

                end if ! Last range?
              end if ! Beginning/end of range

            end do ! Intermediate points

            ! Extreme
            z(lnz) = Tau_max

            ! And return
            return

          end if ! Failed

        end do ! Adjust step

        !
        ! Now build the stratification
        !

        ! First index
        i = 1

        ! Initialize first node to minimum
        z(i) = Tau_min

        ! If natural first range
        if (Tau_min.lt.Strat(1,1)) then

          ! Number of points in first range
          M = int((Strat(1,1) - Tau_min)*iDelt)

          ! If just one node
          if (M.le.1) then

            ! Full step
            i = i + 1
            z(i) = Strat(1,1)

          else

            ! Real step in first natural range
            d = (Strat(1,1) - Tau_min)/dble(M)

            ! For point in first natural range
            do j=1,M-1

              ! Advance
              i = i + 1
              z(i) = z(i-1) + d

            end do

            ! Last
            i = i + 1
            z(i) = Strat(1,1)

          end if ! Valid nodes
        end if ! First natural range

        !
        ! Now the step within the first range proper
        !

        ! Number of points in first range
        M = int((Strat(2,1) - Strat(1,1))*iDelt*Strat(3,1))

        ! If just one node
        if (M.le.1) then

          ! Full step
          i = i + 1
          z(i) = Strat(2,1)

        else

          ! Real step in first range
          d = (Strat(2,1) - Strat(1,1))/dble(M)

          ! For points in first range
          do j=1,M-1

            ! Advance
            i = i + 1
            z(i) = z(i-1) + d

          end do

          ! Last
          i = i + 1
          z(i) = Strat(2,1)

        end if ! Valid nodes

        ! For each range in input
        do k=2,n

          ! If natural range before
          if (Strat(2,k-1).lt.Strat(1,k)) then

            ! Number of points in first range
            M = int((Strat(1,k) - Strat(2,k-1))*iDelt)

            ! If just one node
            if (M.le.1) then

              ! Full step
              i = i + 1
              z(i) = Strat(1,k)

            else

              ! Real step in first range
              d = (Strat(1,k) - Strat(2,k-1))/dble(M)

              ! For points in natural range
              do j=1,M-1

                ! Advance
                i = i + 1
                z(i) = z(i-1) + d

              end do

              ! Last
              i = i + 1
              z(i) = Strat(1,k)

            end if ! Valid nodes
          end if ! First natural range

          !
          ! Now the step within the range proper
          !

          ! Number of points in range
          M = int((Strat(2,k) - Strat(1,k))*iDelt*Strat(3,k))

          ! If just one node
          if (M.le.1) then

            ! Full step
            i = i + 1
            z(i) = Strat(2,k)

          else

            ! Real step in range
            d = (Strat(2,k) - Strat(1,k))/dble(M)

            ! For points in range
            do j=1,M-1

              ! Advance
              i = i + 1
              z(i) = z(i-1) + d

            end do

            ! Last
            i = i + 1
            z(i) = Strat(2,k)

          end if ! Valid nodes

        end do ! Input ranges

        ! If natural last range
        if (Strat(2,n).lt.Tau_max) then

          ! Number of points in last natural range
          M = int((Tau_max - Strat(2,n))*iDelt)

          ! If just one node
          if (M.le.1) then

            ! Full step
            i = i + 1
            z(i) = Tau_max

          else

            ! Real step in first range
            d = (Tau_max - Strat(2,n))/dble(M)

            ! For point in first range
            do j=1,M-1

              ! Advance
              i = i + 1
              z(i) = z(i-1) + d

            end do

            ! Last
            i = i + 1
            z(i) = Tau_max

          end if ! Valid nodes
        end if ! Last natural range
      end if ! Stratification input

      return

      end subroutine gAtmo_Strat

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpret the atmospheric model data in a vectorial buffer\n
      !!     buffer(double(:)): Vector with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data\n
      !!      dims(integer(:)): Grid dimensions (X,Y,Z)
      subroutine rAtmo_frombuffer(buffer,Input,Atmo,Bfield,dims)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(out):: Bfield
      integer, dimension(:), intent(in):: dims
      double precision, dimension(:), target, intent(inout):: buffer

      ! Local

      integer:: iz

      double precision:: ikbcgs


      ! Constant
      ikbcgs = 1d-7/kb

      ! Initialize zero if not allocated
      allocate(Atmo%zeros(nz))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%zeros)
      Atmo%zeros = 0d0

      ! Scale, reference frequency, and nodes
      Atmo%scal = Input%atm_scale
      Atmo%tfreq = Input%omega_ref
      Atmo%nz = dims(3)

      !
      ! Identify type of scale and point to proper height
      !

      ! Tau scale
      if (Atmo%scal.eq.'T') then

        ! Set to tau scale and point to its position
        ztau = .True.
        Atmo%z => buffer(nz+1:2*nz)

      ! Geometrical scale
      else

        ! Set to geometrical scale and point to its position
        ztau = .False.
        Atmo%z => buffer(1:nz)

      end if ! Vertical scale

      ! If respecting the alternative scale
      if (Input%respect_zalt) then

        ! Allocate array
        allocate(Atmo%zalt(nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%zalt)

        ! Tau scale
        if (Atmo%scal.eq.'T') then

          ! Read height
          Atmo%zalt = buffer(1:nz)

        ! Geometrical scale
        else

          ! Read tau scale
          Atmo%zalt = buffer(nz+1:2*nz)

        end if ! Vertical scale
      end if ! Respecting alternative scale

      ! Temperature
      Atmo%T => buffer(3*nz+1:4*nz)

      ! Microturbulence
      Atmo%vmi => buffer(12*nz+1:13*nz)

      ! If forced microturbulence
      if (Input%fvmicro.ge.0d0) Atmo%vmi = Input%fvmicro

      ! If static
      if (Input%static) then

        ! Zero velocity
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros
        Atmo%vz => Atmo%zeros

      ! Dynamic
      else

        ! Point to positions
        Atmo%vx => buffer(9*nz+1:10*nz)
        Atmo%vy => buffer(10*nz+1:11*nz)
        Atmo%vz => buffer(11*nz+1:12*nz)

      end if

      ! Type of atmosphere
      Atmo%typo = Input%atmo_char

      ! Allocate nH and ne
      allocate(Atmo%nH(nz,6),Atmo%ne(nz))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nH)
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%ne)

      ! H^- density
      allocate(Atmo%nhm(nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhm)

      !
      ! Depending on type of scale
      !

      ! Density (full)
      if (Atmo%typo.eq.0) then

        ! Copy from buffer
        Atmo%ne = buffer(14*nz+1:15*nz)
        Atmo%nH = reshape(buffer(18*nz+1:24*nz), (/ nz, 6 /))

        ! H^-
        Atmo%nhm = buffer(17*nz+1:18*nz)

      else

        ! To zero
        Atmo%nH = 0d0
        Atmo%nhm = 0d0

        ! Only electron density
        if (Atmo%typo.eq.1) then

          ! Copy from buffer
          Atmo%ne = buffer(14*nz+1:15*nz)

        ! Electron pressure
        else if (Atmo%typo.eq.2) then

          ! Copy from buffer
          Atmo%ne = buffer(13*nz+1:14*nz)

          ! Transform to electron number density
          Atmo%ne = Atmo%ne*ikbcgs/Atmo%T
          Atmo%typo = 1

        ! Gas pressure
        else if (Atmo%typo.eq.4) then

          ! Allocate
          allocate(Atmo%Pg(nz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)

          ! Copy from buffer
          Atmo%Pg = buffer(4*nz+1:5*nz)
          Atmo%ne = 0d0

        ! Mass density
        else if (Atmo%typo.eq.5) then

          ! Allocate
          allocate(Atmo%rho(nz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%rho)

          ! Copy from buffer
          Atmo%rho = buffer(5*nz+1:6*nz)
          Atmo%ne = 0d0

        end if ! Type of atmospheric model
      end if ! No full number density atmospheric model

      ! Total hydrogen density
      allocate(Atmo%nht(nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nht)

      ! Atomic hydrogen density
      allocate(Atmo%nha(nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nha)

      !
      ! Initialize
      !

      ! If number densities
      if (Atmo%typo.eq.0) then

        !
        ! Read hydrogen densityes
        !

        ! For each height
        do iz=1,nZ

          ! And the total in atomic
          Atmo%nha(iz) = sum(Atmo%nH(iz,:)) + Atmo%nhm(iz)

          ! And whole
          Atmo%nht(iz) = Atmo%nha(iz)

        end do ! For each height

      ! No number densities
      else

        ! Initialize H to zero
        Atmo%nh = 0d0
        Atmo%nht = 0d0
        Atmo%nha = 0d0

      end if ! Type of input densities

      ! Just make a flag to know that there is no helium input
      allocate(Atmo%nhe(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
      Atmo%nhe(1,1) = -1

      !
      ! Unit conversions
      !

      ! Transform height to cgs
      if (.not.ztau) Atmo%z = Atmo%z*1d5

      ! Divide velocities by c (1d5*1d-11/cbar)
      Atmo%vx = Atmo%vx*1d-6/c
      Atmo%vy = Atmo%vy*1d-6/c
      Atmo%vz = Atmo%vz*1d-6/c
      Atmo%vmi = Atmo%vmi*1d-6/c

      ! Check if dynamic (yes if > 1m/s)
      dyn = maxval(Atmo%vx(:)*Atmo%vx(:) + Atmo%vy(:)*Atmo%vy(:) + &
                   Atmo%vz(:)*Atmo%vz(:)).gt.TINYVEL

      !
      ! Magnetic field
      !

      ! Allocate if not already
      allocate(Bfield%Bstrength(nz))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bstrength)
      allocate(Bfield%Btheta(nz))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Btheta)
      allocate(Bfield%Bphi(nz))
      MRAMc = MRAMc + 1d-6*sizeof(Bfield%Bphi)

      ! No magnetic field
      if (Input%unmagnetized) then

        Atmo%Bx => Atmo%zeros
        Atmo%By => Atmo%zeros
        Atmo%Bz => Atmo%zeros
        Bfield%Bstrength = 0d0
        Bfield%Btheta = 0d0
        Bfield%Bphi = 0d0

      ! Yes magnetic field
      else

        ! Get from buffer
        Atmo%Bx => buffer(6*nz+1:7*nz)
        Atmo%By => buffer(7*nz+1:8*nz)
        Atmo%Bz => buffer(8*nz+1:9*nz)

        ! Compute module and angles for B at each height
        do iz=1,nz

          ! Module
          Bfield%Bstrength(iz) = sqrt(Atmo%Bx(iz)*Atmo%Bx(iz) + &
                                      Atmo%By(iz)*Atmo%By(iz) + &
                                      Atmo%Bz(iz)*Atmo%Bz(iz))

          ! There is a magnetic field
          if (Bfield%Bstrength(iz).gt.TINYB) then

            ! Get angles
            Bfield%Btheta(iz) = acos(Atmo%Bz(iz)/Bfield%Bstrength(iz))
            Bfield%Bphi(iz) = atan2(Atmo%By(iz),Atmo%Bx(iz))

          ! There is no field
          else

            ! Trivial angles
            Bfield%Bstrength(iz) = 0d0
            Bfield%Btheta(iz) = 0d0
            Bfield%Bphi(iz) = 0d0

          end if ! Magnetic field

        end do ! Heights

      end if ! There is a magnetic field

      ! Type of allocation in this model
      Atmo%alloc_a = .False.
      Atmo%alloc_b = .False.

      return

      end subroutine rAtmo_frombuffer

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the Atmo structure in the main loop of the CLE
      !! synthesis\n
      !!  Input(Input_class): Structure with configuration data\n
      !!    Atmo(Atmo_class): Structure with atmospheric data\n
      !!       mode(integer): Type of model atmosphere in CLE
      !!                      synthesis\n
      !!       norm(integer): If the geometrical axes are normalized
      subroutine rAtmo_cle_prep(Input,Atmo,mode,norm)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Atmo_class), intent(inout):: Atmo
      integer, intent(in):: mode,norm


      ! Store mode and type of normalization in structure
      Atmo%mode = mode
      Atmo%norm = norm

      ! Allocate zeros
      allocate(Atmo%zeros(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%zeros)
      Atmo%zeros = 0d0

      ! Initialize velocity is static
      if (Input%static) then
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros
        Atmo%vz => Atmo%zeros
        dyn = .False.
      end if

      ! Type of atmosphere
      Atmo%typo = Input%atmo_char

      ! Allocations
      allocate(Atmo%nH(1,6))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nH)
      allocate(Atmo%nht(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nht)
      allocate(Atmo%nha(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nha)
      allocate(Atmo%nhm(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhm)
      allocate(Atmo%nhe(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%nhe)
      allocate(Atmo%ne(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%ne)

      return

      end subroutine rAtmo_cle_prep

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpret the position in space of the current CLE node from
      !! the vectorial buffer and setup the location of variables
      !! depending on the type of model\n
      !!   buffer(double(:)): Vector with atmospheric data\n
      !!  Input(Input_class): Structure with configuration data\n
      !!        x(double(:)): X axis\n
      !!        y(double(:)): Y axis\n
      !!        z(double(:)): Z axis\n
      !!    Atmo(Atmo_class): Structure with atmospheric data\n
      !!    dims(integer(:)): Grid dimensions (X,Y,Z)
      subroutine rAtmo_cle_init(buffer,Input,x,y,z,Atmo,dims)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Input_class), intent(in):: Input
      integer, dimension(:), intent(in):: dims
      double precision, intent(in):: y,z
      double precision, dimension(:), intent(in), target:: x
      double precision, dimension(:), intent(inout), target:: buffer


      ! IMPORTANT: The X dimension, which is the LOS axis, is stored
      ! in z variables for compatibility reasons
      Atmo%nz = dims(1)
      nz = dims(1)

      ! Cartesian mode
      if (Atmo%mode.eq.0) then

        ! Y,Z position in PoS
        Atmo%ypos = y
        Atmo%zpos = z

        ! LOS axis
        Atmo%z => x

        ! Indexes for buffer call
        Atmo%d0 = nz
        Atmo%di0 = 1
        Atmo%di1 = nz
        Atmo%it = 1
        Atmo%imi = 10
        Atmo%iv = 7
        Atmo%ine = 12
        Atmo%inh = 16
        Atmo%ipe = 11
        Atmo%ipg = 2
        Atmo%irh = 3
        Atmo%ib = 4

      ! Slab
      elseif (Atmo%mode.eq.1) then

        ! Y has theta angle
        Atmo%ypos = buffer(2)

        ! Z has optical depth
        Atmo%zpos = buffer(3)

        ! Z axis has height
        Atmo%z => buffer(1:1)
        if (Atmo%norm.le.0) Atmo%z = Atmo%z/Input%R_star

        ! Indexes for buffer call
        Atmo%d0 = 1
        Atmo%di0 = 0
        Atmo%di1 = 0
        Atmo%it = 5
        Atmo%imi = 14
        Atmo%iv = 11
        Atmo%ine = 16
        Atmo%inh = 20
        Atmo%ipe = 15
        Atmo%ipg = 6
        Atmo%irh = 7
        Atmo%ib = 8

      ! non-cartesian
      else if (Atmo%mode.eq.2) then

        ! LOS axis
        Atmo%z => buffer(1:nz)

        !
        ! Y,Z position in PoS and units
        !

        ! If not normalized
        if (Atmo%norm.le.0) then

          ! Normalize now
          Atmo%ypos = y/Input%R_star
          Atmo%zpos = z/Input%R_star
          Atmo%z = Atmo%z/Input%R_star

        ! Already normalized
        else

          ! Just copy
          Atmo%ypos = y
          Atmo%zpos = z

        end if

        ! Indexes for buffer call
        Atmo%d0 = nz
        Atmo%di0 = 1
        Atmo%di1 = nz
        Atmo%it = 2
        Atmo%imi = 11
        Atmo%iv = 8
        Atmo%ine = 13
        Atmo%inh = 17
        Atmo%ipe = 12
        Atmo%ipg = 3
        Atmo%irh = 4
        Atmo%ib = 5

      end if

      ! Reset values
      Atmo%nht = 0d0
      Atmo%nha = 0d0
      Atmo%nhm = 0d0
      Atmo%nhe(1,1) = -1

      return

      end subroutine rAtmo_cle_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpret the atmospheric model data in a vectorial buffer in
      !! the CLE synthesis\n
      !!   buffer(double(:)): Column with atmospheric data\n
      !!  Input(Input_class): Structure with configuration data\n
      !!    Atmo(Atmo_class): Structure with atmospheric data\n
      !!         ix(integer): Position index along the LOS axis
      subroutine rAtmo_cle(buffer,Input,Atmo,ix)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Atmo_class), intent(inout):: Atmo
      double precision, dimension(:), intent(inout), target:: buffer
      integer, intent(in):: ix

      ! Local

      integer:: iz,i0
      double precision:: ikbcgs
      double precision, dimension(Atmo%nz,6):: daux


      ! Constant
      ikbcgs = 1d-7/kb

      ! Coordinate
      iz = Atmo%di0*ix

      ! Temperature
      i0 = Atmo%it*Atmo%d0
      Atmo%T => buffer(i0+iz:i0+iz)

      ! Micro
      i0 = Atmo%imi*Atmo%d0
      Atmo%vmi => buffer(i0+iz:i0+iz)

      ! If forced micro
      if (Input%fvmicro.ge.0d0) Atmo%vmi = Input%fvmicro

      ! No velocity
      if (Input%static) then

        ! Point to zero
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros
        Atmo%vz => Atmo%zeros
        dyn = .False.

      ! Not necessarily static
      else

        ! Point to components
        i0 = Atmo%iv*Atmo%d0
        Atmo%vx => buffer(i0+iz:i0+iz)
        i0 = i0 + Atmo%d0
        Atmo%vy => buffer(i0+iz:i0+iz)
        i0 = i0 + Atmo%d0
        Atmo%vz => buffer(i0+iz:i0+iz)

        ! Slab
        if (Atmo%mode.eq.1) then

          ! Unit conversions
          ! Divide velocities by c (1d5*1d-11/cbar)
          Atmo%vx = Atmo%vx*1d-6/c

          ! Check if dynamic (yes if > 1m/s)
          if (maxval(Atmo%vx).gt.TINYVEL) then

            dyn = .True.

          ! Not dynamic
          else

            dyn = .False.

          end if ! dynamics

        ! No slab
        else

          ! Unit conversions
          ! Divide velocities by c (1d5*1d-11/cbar)
          Atmo%vx = Atmo%vx*1d-6/c
          Atmo%vy = Atmo%vy*1d-6/c
          Atmo%vz = Atmo%vz*1d-6/c

          ! Check if dynamic (yes if > 1m/s)
          if (maxval(Atmo%vx*Atmo%vx + Atmo%vy*Atmo%vy + &
                     Atmo%vz*Atmo%vz).gt.TINYVEL) then

            dyn = .True.

          ! Not dynamic
          else

            dyn = .False.

          end if ! dynamics
        end if ! Model type
      end if ! Forced static

      !
      ! Depending on type of scale
      !

      ! Density (full)
      if (Atmo%typo.eq.0) then

        ! Get electron number density
        i0 = Atmo%ine*Atmo%d0
        Atmo%ne = buffer(i0+iz:i0+iz)
        i0 = Atmo%inh*Atmo%d0

        ! If slab
        if (Atmo%mode.eq.1) then

          ! Get hydrogen number density
          Atmo%nH = reshape(buffer(Atmo%inh:Atmo%inh+5), (/ nz, 6 /))

        ! Others
        else

          ! Reshape from buffer
          daux = reshape(buffer(Atmo%inh*Atmo%d0+Atmo%di0: &
                                Atmo%inh*Atmo%d0+Atmo%di1*6), &
                         (/ Atmo%nz, 6 /))

          ! Copy in array
          Atmo%nH = daux(ix:ix,:)

        end if ! Slab or other model

      ! No full number density
      else

        ! Only electron number density
        if (Atmo%typo.eq.1) then

          ! Copy electron number density
          i0 = Atmo%ine*Atmo%d0
          Atmo%ne = buffer(i0+iz:i0+iz)

        ! Electron pressure
        else if (Atmo%typo.eq.2) then

          ! Copy electron pressure
          i0 = Atmo%ipe*Atmo%d0
          Atmo%ne = buffer(i0+iz:i0+iz)

          ! Transform into electron number density
          Atmo%ne = Atmo%ne*ikbcgs/Atmo%T

        ! Gas pressure
        else if (Atmo%typo.eq.4) then

          ! Allocate space
          i0 = Atmo%ipg*Atmo%d0
          allocate(Atmo%Pg(1))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)

          ! Get gas pressure and reset electron number density
          Atmo%Pg = buffer(i0+iz:i0+iz)
          Atmo%ne = 0d0

        ! Mass density
        else if (Atmo%typo.eq.5) then

          ! Allocate space
          i0 = Atmo%irh*Atmo%d0
          allocate(Atmo%rho(1))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo%rho)

          ! Get mass density and reset electron number density
          Atmo%rho = buffer(i0+iz:i0+iz)
          Atmo%ne = 0d0

        end if ! Type of model if not full number density
      end if ! Full number density

      ! Reset H^-
      Atmo%nhm = 0d0

      !
      ! If number densities
      !
      if (Atmo%typo.eq.0) then

        ! Add the total in atomic hydrogen number density
        Atmo%nha = sum(Atmo%nH(1,:))

        ! Add whole hydrogen number density
        Atmo%nht = Atmo%nha

      ! No number densities
      else

        ! Initialize H to zero
        Atmo%nh = 0d0
        Atmo%nht = 0d0
        Atmo%nha = 0d0

      end if ! Type of input densities

      ! Just make a flag to know that there is no helium input
      Atmo%nhe(1,1) = -1

      !
      ! Unit conversions
      !

      ! Divide microturbulent velocity by c (1d5*1d-11/cbar)
      Atmo%vmi = Atmo%vmi*1d-6/c

      !
      ! Magnetic field
      !

      ! No magnetic field
      if (Input%unmagnetized) then

        ! Point to zero
        Atmo%Bx => Atmo%zeros
        Atmo%By => Atmo%zeros
        Atmo%Bz => Atmo%zeros

      ! Yes magnetic field
      else

        ! Get from buffer
        i0 = Atmo%ib*Atmo%d0
        Atmo%Bx => buffer(i0+iz:i0+iz)
        i0 = i0 + Atmo%d0
        Atmo%By => buffer(i0+iz:i0+iz)
        i0 = i0 + Atmo%d0
        Atmo%Bz => buffer(i0+iz:i0+iz)

      end if ! There is a magnetic field

      return

      end subroutine rAtmo_cle

!#####################################################################
!#####################################################################
!#####################################################################

      end module ratmo_mod
