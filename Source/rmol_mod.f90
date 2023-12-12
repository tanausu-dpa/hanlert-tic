      !> Reading of molecular data
      module rmol_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     11/24/2023 V3.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/24/2023:    V3.0.4 - Added Moleq_T (TdPA)
!
!     09/26/2023:    V3.0.3 - Bugfix: Needed absolute value for
!                             energy, not mavxal (TdPA)
!
!     09/25/2023:    V3.0.2 - Receive frequencies in cm^-1 and
!                             transform after reading (TdPA)
!
!     07/03/2023:    V3.0.1 - Now preparemol only allocates and
!                             initializes the number density (TdPA)
!                           - Added setupmol_eq and setupmol_pf to
!                             allocate and prepare Mol%eq and Mol%pf
!                             variables (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Variables with height dimension are
!                                now allocated elsewhere.
!                              o Atmo is no longer an input for
!                                rMol.
!                              o Changed aborted calls to gaborted
!                                in rMol.
!                              o Added routine preparemol.
!                             (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     11/19/2019:    V1.1.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!                           - Checks for success of python routine
!                             and unit is now 100 (TdPA)
!                           - Added allocatemol (TdPA)
!
!     11/26/2018:    V1.0.3 - Removed an empty loop (TdPA)
!
!     09/14/2017:    V1.0.2 - Added a path and ID to the file (TdPA)
!
!     07/20/2017:    V1.0.1 - Changed the declaration of T in
!                             Molpf and Moleq to fixed size (TdPA)
!
!     04/18/2017:    V1.0.0 - Started coding (TdPA)
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
!  This subroutine reads the molecular data
!
!  rMol:
!    The proper reader
!
!  Molpf:
!    Interpolates the partition function given the type of fit
!
!  Moleq:
!    Interpolates the equilibrium constant given the type of fit
!
!  Moleq_T:
!    Interpolates the equilibrium constant given the type of fit
!  for a given temperature
!
!  allocatemol:
!    Allocates array of Mol_class
!
!  preparemol:
!    Allocates variables the molecule number density and initializes
!  them to zero
!
!  setupmol_eq:
!    Prepares molecule %eq data for current model
!
!  setupmol_pf:
!    Prepares molecule %pf data for current model
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : dopp , ktoev , kb , fktoJ , fktoev
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Reads a file with molecular data.\n
      !!   filename(character(:)): Name of the file to read\n
      !!     source(character(:)): Path to the source code\n
      !!         ID(character(:)): ID of this run\n
      !!           Mol(Mol_class): Structure with the molecule data
      subroutine rMol(filename,source,ID,Mol)

      ! I/O

      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      type(Mol_class), intent(inout):: Mol

      ! Local

      character(len=2):: cdump

      integer:: ios, ii


      ! Routine name
      urou = 'rMol'


      ! Read the atomic data
      if(pid.eq.0) call system('python '//trim(source)// &
                               'rmol.py '//trim(filename)//' '// &
                               ID//' '//verbosef)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      open(100,file='tmp_mol_'//ID,status='old',iostat=ios,err=1000)

      ! Success
      read (100,*,err=1100) ios

      ! If no correct file, abort
      if (ios.lt.0) then

        umsg = 'Problem translating the molecular file '// &
               trim(filename)
        goto 1200

      end if

      ! Molecule name
      ! Read the number of characters and allocate the string
      read(100,*,err=1100) ii
      allocate(character(len=ii) :: Mol%Molecule)

      ! Read the proper name
      read(100,*,err=1100) Mol%Molecule

      ! Mass
      read(100,*,err=1100) Mol%rmass

      ! Charge
      read(100,*,err=1100) Mol%Charge

      ! Components
      read(100,*,err=1100) Mol%nA
      allocate(Mol%natom(Mol%nA))
      allocate(Mol%catom(Mol%nA))
      allocate(Mol%iatom(Mol%nA))
      Mol%nAT = 0
      do ii=1,Mol%nA
        read(100,*,err=1100) Mol%natom(ii)
        read(100,*,err=1100) cdump
        if (cdump(2:2).eq.' ') then
          cdump(2:2) = cdump(1:1)
          cdump(1:1) = ' '
        end if
        Mol%catom(ii)%s = cdump
        Mol%nAT = Mol%nAT + Mol%natom(ii)
      end do

      ! Dissociation energy
      read(100,*,err=1100) Mol%Den
      Mol%Den = Mol%Den*1d-5

      ! Sanity check
      if (abs(Mol%Den).lt.1d-4) then

        write(umsg,'(A)') ' # Dissociation energy in '// &
             trim(Mol%Molecule)//' molecule is smaller than '// &
             '10 cm^-1, you may be using an old atomic model format'
        call verbose

      end if

      ! Type of partition
      read(100,*,err=1100) Mol%pffit

      ! Minimum and maximum temperatures
      read(100,*,err=1100) Mol%Tmin, Mol%Tmax

      ! Partition function coefficients
      read(100,*,err=1100) Mol%npfcoeff
      if (Mol%npfcoeff.gt.0) then
        allocate(Mol%pfcoeff(Mol%npfcoeff))
        read(100,*,err=1100) Mol%pfcoeff(Mol%npfcoeff:1:-1)
      end if

      ! Equilibrium constant coefficients
      read(100,*,err=1100) Mol%neqcoeff
      if (Mol%neqcoeff.gt.0) then
        allocate(Mol%eqcoeff(Mol%neqcoeff))
        read(100,*,err=1100) Mol%eqcoeff(Mol%neqcoeff:1:-1)
      end if

      ! Calculate the molecular part of the Doppler width
      Mol%cDopp = dopp/sqrt(Mol%rmass)

      close (100)

      ! Control that everything went fine
      call control

      ! Remove temporal file
      if (pid.eq.0) then
        call system('rm tmp_mol_'//ID)
        umsg = ' - Molecule '//trim(filename)//' read'
        call verbose
      end if

      return

1000  umsg = 'Error opening molecule file'
      urou = 'rMol'
      call gaborted
1100  umsg = 'Error reading molecule file'
1200  close(100)
      urou = 'rMol'
      call gaborted

      end subroutine rMol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the molecular partition function for the atmospheric
      !! temperature.\n
      !!    Mol(Mol_class): Structure with the molecule data\n
      !!      T(dfloat(:)): Temperature
      subroutine Molpf(Mol,T)

      ! I/O

      type(Mol_class), intent(inout):: Mol
      double precision, dimension(nz), intent(in):: T

      ! Local

      integer:: ii,iz

      double precision:: x

      ! Routine name
      urou = 'Molpf'

      ! Initialize
      Mol%pf = 0d0

      ! If there are no coefficients, return
      if (Mol%npfcoeff.eq.0) return

      ! For each height node temperature
      do iz=1,NZ

        ! If temperature out of limits of molecule existence, skip
        if (T(iz).lt.Mol%Tmin.or.T(iz).gt.Mol%Tmax) cycle

        Mol%pf(iz) = Mol%pfcoeff(1)

        if (Mol%pffit.eq.0) then

          x = T(iz)

          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          Mol%pf(iz) = exp(Mol%pf(iz))

        else if (Mol%pffit.eq.1) then

          x = T(iz)*1d-4

          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          Mol%pf(iz) = exp(Mol%pf(iz))

        else if (Mol%pffit.eq.2) then

          x = log10(ktoev/T(iz))

          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          Mol%pf(iz) = 1d1**(Mol%pf(iz))

        else if (Mol%pffit.eq.3) then

          x = log(T(iz))

          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          Mol%pf(iz) = exp(Mol%pf(iz))

        else if (Mol%pffit.eq.4) then

          return

        else

          umsg = 'Unknown fit for molecular partition '// &
                 'function'
          call aborted

        end if

      end do

      end subroutine Molpf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the molecular equilibrium parameter for the
      !! atmospheric temperature\n
      !!    Mol(Mol_class): Structure with the molecule data\n
      !!      T(dfloat(:)): Temperature
      subroutine Moleq(Mol,T)

      ! I/O

      type(Mol_class), intent(inout):: Mol
      double precision, dimension(nz), intent(in):: T

      ! Local

      integer:: ii,iz

      double precision:: x1, x2, x3, C0, fact


      ! Routine name
      urou = 'Moleq'


      ! Initialize
      Mol%eq = 0d0

      ! If there are no coefficients, return
      if (Mol%neqcoeff.eq.0) return

      ! For each height node temperature
      do iz=1,NZ

        ! If temperature out of limits of molecule existence, skip
        if (T(iz).lt.Mol%Tmin.or.T(iz).gt.Mol%Tmax) cycle

        Mol%eq(iz) = Mol%eqcoeff(1)

        if (Mol%pffit.eq.0) then

          x1 = T(iz)
          x2 = fktoJ/kb/T(iz)
          x3 = log(T(iz))
          C0 = Mol%nAT - 1 - Mol%charge

          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x1 + Mol%eqcoeff(ii)
          end do

          Mol%eq(iz) = exp(Mol%Den*x2 + Mol%eq(iz) - 1.5d0*C0*x3)

        else if (Mol%pffit.eq.1) then

          x1 = T(iz)*1d-4
          x2 = fktoJ/kb/T(iz)
          x3 = log(T(iz))
          C0 = Mol%nAT - 1 - Mol%charge

          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x1 + Mol%eqcoeff(ii)
          end do

          Mol%eq(iz) = exp(Mol%Den*x2 + Mol%eq(iz) - 1.5d0*C0*x3)

        else if (Mol%pffit.eq.2.OR.Mol%pffit.eq.3) then

          x1 = ktoev/T(iz)
          x2 = kb*T(iz)
          x3 = log10(x1)
          C0 = Mol%nAT - 1 - Mol%charge

          fact = 1d4/(1d-2**C0) !J->erg and m->cm / 10^3 ad-hoc

          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x3 + Mol%eqcoeff(ii)
          end do

          Mol%eq(iz) = 1d1**(Mol%Den*fktoev*x1 - Mol%eq(iz))*x2*fact

        else if (Mol%pffit.eq.4) then

          x1 = ktoev/T(iz)
          x2 = kb*T(iz)

          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x1 + Mol%eqcoeff(ii)
          end do

          Mol%eq(iz) = 1d1**(-Mol%eq(iz))*x2*x2

        else

          umsg = ' # Unknown fit for equilibrium constant '// &
                 'function'
          call aborted

        end if

      end do

      end subroutine Moleq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the molecular equilibrium parameter for a given
      !! temperature\n
      !!    Mol(Mol_class): Structure with the molecule data\n
      !!         T(dfloat): Temperature
      double precision function Moleq_T(Mol,T)

      ! I/O

      type(Mol_class), intent(inout):: Mol
      double precision, intent(in):: T

      ! Local

      integer:: ii

      double precision:: x1, x2, x3, C0, fact


      ! Routine name
      urou = 'Moleq_T'


      ! Initialize
      Moleq_T = 0d0

      ! If there are no coefficients, return
      if (Mol%neqcoeff.eq.0) return

      ! If temperature out of limits of molecule existence, skip
      if (T.lt.Mol%Tmin.or.T.gt.Mol%Tmax) return

      Moleq_T = Mol%eqcoeff(1)

      if (Mol%pffit.eq.0) then

        x1 = T
        x2 = fktoJ/kb/T
        x3 = log(T)
        C0 = Mol%nAT - 1 - Mol%charge

        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x1 + Mol%eqcoeff(ii)
        end do

        Moleq_T = exp(Mol%Den*x2 + Moleq_T - 1.5d0*C0*x3)

      else if (Mol%pffit.eq.1) then

        x1 = T*1d-4
        x2 = fktoJ/kb/T
        x3 = log(T)
        C0 = Mol%nAT - 1 - Mol%charge

        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x1 + Mol%eqcoeff(ii)
        end do

        Moleq_T = exp(Mol%Den*x2 + Moleq_T - 1.5d0*C0*x3)

      else if (Mol%pffit.eq.2.or.Mol%pffit.eq.3) then

        x1 = ktoev/T
        x2 = kb*T
        x3 = log10(x1)
        C0 = Mol%nAT - 1 - Mol%charge

        fact = 1d4/(1d-2**C0) !J->erg and m->cm / 10^3 ad-hoc

        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x3 + Mol%eqcoeff(ii)
        end do

        Moleq_T = 1d1**(Mol%Den*fktoev*x1 - Moleq_T)*x2*fact

      else if (Mol%pffit.eq.4) then

        x1 = ktoev/T
        x2 = kb*T

        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x1 + Mol%eqcoeff(ii)
        end do

        Moleq_T = 1d1**(-Moleq_T)*x2*x2

      else

        umsg = ' # Unknown fit for equilibrium constant '// &
               'function'
        call aborted

      end if

      end function Moleq_T

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates array of Mol_class.\n
      !!  Mol(Mol_class): Structure to allocate\n
      !!     nn(integer): Size to allocate
      subroutine allocatemol(Mol,nn)

      ! I/O
      type(Mol_class), dimension(:), allocatable:: Mol
      integer, intent(in):: nn

      ! Local
      integer:: ios


      ! Routine name
      urou = 'allocatemol'

      ! Allocate 1 even if no molecules
      if (nn.lt.1) then

        allocate(Mol(1), stat=ios)

      ! Allocate molecules to read
      else

        allocate(Mol(nn), stat=ios)

      end if

      ! Control that everything went fine
      call control

      return

      end subroutine allocatemol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate molecule data\n
      !!    Mol(Mol_class): Structure to allocate\n
      !!       nn(integer): Size of molecule array
      subroutine preparemol(Mol,nn)

      ! I/O
      type(Mol_class), dimension(:), intent(inout):: Mol
      integer, intent(in):: nn

      ! Local
      integer:: im

      ! For each molecule
      do im=1,nn

        ! Allocate space por populations and parameters
        allocate(Mol(im)%n(NZ))
        Mol(im)%n = 0d0

      end do

      end subroutine preparemol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepares molecule eq data for current model\n
      !!    Mol(Mol_class): Structure with molecules\n
      !!       nn(integer): Size of molecule array\n
      !!  Atmo(Atmo_class): Structure with atmospheric data
      subroutine setupmol_eq(Mol,nn,Atmo)

      ! I/O
      type(Mol_class), dimension(:), intent(inout):: Mol
      type(Atmo_class), intent(in):: Atmo
      integer, intent(in):: nn

      ! Local
      logical:: free
      integer:: im

      !
      ! Check if need to free input data
      !
      free = run_mode.eq.0

      ! For each molecule
      do im=1,nn

        ! Make sure allocated
        if (allocated(Mol(im)%eq)) then
          if (size(Mol(im)%eq).ne.NZ) then
            deallocate(Mol(im)%eq)
            allocate(Mol(im)%eq(NZ))
          end if
        else
          allocate(Mol(im)%eq(NZ))
        end if

        ! Interpolate into T grid
        call Moleq(Mol(im),Atmo%T)
        if (Mol(im)%neqcoeff.gt.0.and.free) &
          deallocate(Mol(im)%eqcoeff)

      end do

      end subroutine setupmol_eq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepares molecule pf data for current model\n
      !!    Mol(Mol_class): Structure with molecules\n
      !!       nn(integer): Size of molecule array\n
      !!  Atmo(Atmo_class): Structure with atmospheric data
      subroutine setupmol_pf(Mol,nn,Atmo)

      ! I/O
      type(Mol_class), dimension(:), intent(inout):: Mol
      type(Atmo_class), intent(in):: Atmo
      integer, intent(in):: nn

      ! Local
      logical:: free
      integer:: im

      !
      ! Check if need to free input data
      !
      free = run_mode.eq.0

      ! For each molecule
      do im=1,nn

        ! Make sure allocated
        if (allocated(Mol(im)%pf)) then
          if (size(Mol(im)%pf).ne.NZ) then
            deallocate(Mol(im)%pf)
            allocate(Mol(im)%pf(NZ))
          end if
        else
          allocate(Mol(im)%pf(NZ))
        end if

        ! Interpolate into T grid
        call Molpf(Mol(im),Atmo%T)
        if (Mol(im)%npfcoeff.gt.0.and.free) &
          deallocate(Mol(im)%pfcoeff)

      end do

      end subroutine setupmol_pf

!#####################################################################
!#####################################################################
!#####################################################################

      end module rmol_mod
