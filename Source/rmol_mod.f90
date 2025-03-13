      !> Reading of molecular data
      module rmol_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     18/04/2017
!  Last version:
!     17/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     17/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  rMol
!    Read a model molecule from the specified file
!
!  Molpf
!    Calculate the partition function given a temperature
!  stratification
!
!  Moleq
!    Calculate the equilibrium constant given a temperature
!  stratification
!
!  Moleq_T
!    Calculate the equilibrium constant given the type of fit
!  for a given temperature
!
!  allocatemol
!    Allocate an array of Mol_class
!
!  preparemol
!    Allocate the molecules number density and initialize them to zero
!
!  setupmol_eq
!    Calculate the equilibrium constant of all molecules in a given
!  model atmosphere
!
!  setupmol_pf
!    Calculate the partition function of all molecules in a given
!  model atmosphere
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

      !> Read a model molecule from the specified file\n
      !!   filename(character(:)): Name of the file to read from\n
      !!     source(character(:)): Path to the source code\n
      !!         ID(character(:)): ID of this run\n
      !         Mol(Mol_class(:)): Structures with molecular data
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

      ! Master translate the molecule model with python
      if(pid.eq.0) call system('python '//trim(source)// &
                               'rmol.py '//trim(filename)//' '// &
                               ID//' '//verbosef)

      ! Wait till master is done
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! Open translated file
      open(100,file='tmp_mol_'//ID,status='old',iostat=ios,err=1000)

      ! Read success
      read (100,*,err=1100) ios

      ! If no correct file
      if (ios.lt.0) then

        ! Issue error
        umsg = 'Problem translating the molecular file '// &
               trim(filename)
        goto 1200

      end if ! Wrong file

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Mol)

      !
      ! Molecule name
      !

      ! Read the number of characters and allocate the string
      read(100,*,err=1100) ii

      ! Allocate string for molecule name
      allocate(character(len=ii) :: Mol%Molecule)
      MRAMc = MRAMc + 1d-6*sizeof(Mol%molecule)

      ! Read the proper name
      read(100,*,err=1100) Mol%Molecule

      ! Mass
      read(100,*,err=1100) Mol%rmass

      ! Charge
      read(100,*,err=1100) Mol%Charge

      !
      ! Components
      !

      ! Number of components
      read(100,*,err=1100) Mol%nA

      ! Allocate components
      allocate(Mol%natom(Mol%nA))
      MRAMc = MRAMc + 1d-6*sizeof(Mol%natom)
      allocate(Mol%catom(Mol%nA))
      allocate(Mol%iatom(Mol%nA))
      MRAMc = MRAMc + 1d-6*sizeof(Mol%iatom)

      ! Initialize number of atoms per species
      Mol%nAT = 0

      ! For each component
      do ii=1,Mol%nA

        ! Memory count
        MRAMc = MRAMc + 1d-6*sizeof(Mol%catom(ii))

        ! Read number of atoms and name
        read(100,*,err=1100) Mol%natom(ii)
        read(100,*,err=1100) cdump

        ! Correct single letter names
        if (cdump(2:2).eq.' ') then
          cdump(2:2) = cdump(1:1)
          cdump(1:1) = ' '
        end if

        ! Save name and number of atoms
        Mol%catom(ii)%s = cdump
        Mol%nAT = Mol%nAT + Mol%natom(ii)

      end do ! Components

      ! Dissociation energy
      read(100,*,err=1100) Mol%Den

      ! Convert units to kaiser
      Mol%Den = Mol%Den*1d-5

      ! Sanity check
      if (abs(Mol%Den).lt.1d-4) then

        ! Issue error
        write(umsg,'(A)') ' # Dissociation energy in '// &
             trim(Mol%Molecule)//' molecule is smaller than '// &
             '10 cm^-1, you may be using an old atomic model format'
        call verbose

      end if ! Sanity check

      ! Type of partition
      read(100,*,err=1100) Mol%pffit

      ! Minimum and maximum temperatures
      read(100,*,err=1100) Mol%Tmin, Mol%Tmax

      ! Partition function coefficients
      read(100,*,err=1100) Mol%npfcoeff
      if (Mol%npfcoeff.gt.0) then
        allocate(Mol%pfcoeff(Mol%npfcoeff))
        MRAMc = MRAMc + 1d-6*sizeof(Mol%pfcoeff)
        read(100,*,err=1100) Mol%pfcoeff(Mol%npfcoeff:1:-1)
      end if

      ! Equilibrium constant coefficients
      read(100,*,err=1100) Mol%neqcoeff
      if (Mol%neqcoeff.gt.0) then
        allocate(Mol%eqcoeff(Mol%neqcoeff))
        MRAMc = MRAMc + 1d-6*sizeof(Mol%neqcoeff)
        read(100,*,err=1100) Mol%eqcoeff(Mol%neqcoeff:1:-1)
      end if

      ! Calculate the molecular part of the Doppler width
      Mol%cDopp = dopp/sqrt(Mol%rmass)

      ! Close file
      close (100)

      ! Control that everything went fine
      call control

      ! Master
      if (pid.eq.0) then

        ! Remove temporal file
        call system('rm tmp_mol_'//ID)

        ! Verbose
        umsg = ' - Molecule '//trim(filename)//' read'
        call verbose

      end if ! Master

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

      !> Calculate the partition function given a temperature
      !! stratification\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!       T(double(:)): Temperature
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

        ! Initialize constant
        Mol%pf(iz) = Mol%pfcoeff(1)

        ! Kurucz 1970
        if (Mol%pffit.eq.0) then

          ! Get temperature
          x = T(iz)

          ! Complete calculation
          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          ! Take exponential
          Mol%pf(iz) = exp(Mol%pf(iz))

        ! Kurucz 1985
        else if (Mol%pffit.eq.1) then

          ! Scale temperature
          x = T(iz)*1d-4

          ! Complete calculation
          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          ! Take exponential
          Mol%pf(iz) = exp(Mol%pf(iz))

        ! Sauval and Tatum 1984
        else if (Mol%pffit.eq.2) then

          ! Get logarithm of temperature in eV
          x = log10(ktoev/T(iz))

          ! Complete calculation
          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          ! Take power of 10
          Mol%pf(iz) = 1d1**(Mol%pf(iz))

        ! Irwin 1981
        else if (Mol%pffit.eq.3) then

          ! Get temperature logarithm
          x = log(T(iz))

          ! Complete calculation
          do ii=2,Mol%npfcoeff
            Mol%pf(iz) = Mol%pf(iz)*x + Mol%pfcoeff(ii)
          end do

          ! Take exponential
          Mol%pf(iz) = exp(Mol%pf(iz))

        ! Tsuji 1973
        else if (Mol%pffit.eq.4) then

          ! Constant
          return

        ! Unknown
        else

          ! Issue error
          umsg = 'Unknown fit for molecular partition '// &
                 'function'
          call aborted

        end if ! Type of fit

      end do ! Heights

      end subroutine Molpf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the equilibrium constant given a temperature
      !! stratification\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!       T(double(:)): Temperature
      subroutine Moleq(Mol,T)

      ! I/O

      type(Mol_class), intent(inout):: Mol
      double precision, dimension(nz), intent(in):: T

      ! Local

      integer:: ii,iz

      double precision:: x1,x2,x3,C0,fact


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

        ! Initialize equilibrium constant
        Mol%eq(iz) = Mol%eqcoeff(1)

        ! Kurucz 1970
        if (Mol%pffit.eq.0) then

          ! Get inverse kb*T and its logarithm
          x1 = T(iz)
          x2 = fktoJ/kb/T(iz)
          x3 = log(T(iz))
          C0 = Mol%nAT - 1 - Mol%charge

          ! Add contributions
          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x1 + Mol%eqcoeff(ii)
          end do

          ! Get exponential
          Mol%eq(iz) = exp(Mol%Den*x2 + Mol%eq(iz) - 1.5d0*C0*x3)

        ! Kurucz 1985
        else if (Mol%pffit.eq.1) then

          ! Get the inverse of the scaled kb*T and its logarithm
          x1 = T(iz)*1d-4
          x2 = fktoJ/kb/T(iz)
          x3 = log(T(iz))
          C0 = Mol%nAT - 1 - Mol%charge

          ! Add contributions
          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x1 + Mol%eqcoeff(ii)
          end do

          ! Get exponential
          Mol%eq(iz) = exp(Mol%Den*x2 + Mol%eq(iz) - 1.5d0*C0*x3)

        ! Sauval and Tatum 1984 or Irwin 1981
        else if (Mol%pffit.eq.2.OR.Mol%pffit.eq.3) then

          ! Get inverse temperature in eV and its logarithm
          x1 = ktoev/T(iz)
          x2 = kb*T(iz)
          x3 = log10(x1)
          C0 = Mol%nAT - 1 - Mol%charge

          ! Factor
          ! J->erg and m->cm / 10^3 ad-hoc
          fact = 1d4/(1d-2**C0)

          ! Add contributions
          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x3 + Mol%eqcoeff(ii)
          end do

          ! Get power of ten
          Mol%eq(iz) = 1d1**(Mol%Den*fktoev*x1 - Mol%eq(iz))*x2*fact

        ! Tsuji 1973
        else if (Mol%pffit.eq.4) then

          ! Get T in eV and kb*T
          x1 = ktoev/T(iz)
          x2 = kb*T(iz)

          ! Add contributions
          do ii=2,Mol%neqcoeff
            Mol%eq(iz) = Mol%eq(iz)*x1 + Mol%eqcoeff(ii)
          end do

          ! Get power of ten
          Mol%eq(iz) = 1d1**(-Mol%eq(iz))*x2*x2

        ! Unknown fit method
        else

          ! Issue error
          umsg = ' # Unknown fit for equilibrium constant '// &
                 'function'
          call aborted

        end if ! Fit method

      end do ! Heights

      end subroutine Moleq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the equilibrium constant given the type of fit for
      !! a given temperature\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!          T(double): Temperature
      double precision function Moleq_T(Mol,T)

      ! I/O

      type(Mol_class), intent(inout):: Mol
      double precision, intent(in):: T

      ! Local

      integer:: ii

      double precision:: x1,x2,x3,C0,fact


      ! Routine name
      urou = 'Moleq_T'

      ! Initialize
      Moleq_T = 0d0

      ! If there are no coefficients, return
      if (Mol%neqcoeff.eq.0) return

      ! If temperature out of limits of molecule existence, skip
      if (T.lt.Mol%Tmin.or.T.gt.Mol%Tmax) return

      ! Initialize equilibrium constant
      Moleq_T = Mol%eqcoeff(1)

      ! Kurucz 1970
      if (Mol%pffit.eq.0) then

        ! Get inverse kb*T and its logarithm
        x1 = T
        x2 = fktoJ/kb/T
        x3 = log(T)
        C0 = Mol%nAT - 1 - Mol%charge

        ! Add contributions
        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x1 + Mol%eqcoeff(ii)
        end do

        ! Get exponential
        Moleq_T = exp(Mol%Den*x2 + Moleq_T - 1.5d0*C0*x3)

      ! Kurucz 1985
      else if (Mol%pffit.eq.1) then

        ! Get the inverse of the scaled kb*T and its logarithm
        x1 = T*1d-4
        x2 = fktoJ/kb/T
        x3 = log(T)
        C0 = Mol%nAT - 1 - Mol%charge

        ! Add contributions
        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x1 + Mol%eqcoeff(ii)
        end do

        ! Get exponential
        Moleq_T = exp(Mol%Den*x2 + Moleq_T - 1.5d0*C0*x3)

      ! Sauval and Tatum 1984 or Irwin 1981
      else if (Mol%pffit.eq.2.or.Mol%pffit.eq.3) then

        ! Get inverse temperature in eV and its logarithm
        x1 = ktoev/T
        x2 = kb*T
        x3 = log10(x1)
        C0 = Mol%nAT - 1 - Mol%charge

        ! Factor
        ! J->erg and m->cm / 10^3 ad-hoc
        fact = 1d4/(1d-2**C0)

        ! Add contributions
        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x3 + Mol%eqcoeff(ii)
        end do

        ! Get power of ten
        Moleq_T = 1d1**(Mol%Den*fktoev*x1 - Moleq_T)*x2*fact

      ! Tsuji 1973
      else if (Mol%pffit.eq.4) then

        ! Get T in eV and kb*T
        x1 = ktoev/T
        x2 = kb*T

        ! Add contributions
        do ii=2,Mol%neqcoeff
          Moleq_T = Moleq_T*x1 + Mol%eqcoeff(ii)
        end do

        ! Get power of ten
        Moleq_T = 1d1**(-Moleq_T)*x2*x2

      ! Unknown fit method
      else

        ! Issue error
        umsg = ' # Unknown fit for equilibrium constant '// &
               'function'
        call aborted

      end if ! Fit method

      end function Moleq_T

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate an array of Mol_class\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!        nn(integer): Size to allocate
      subroutine allocatemol(Mol,nn)

      ! I/O

      type(Mol_class), dimension(:), allocatable, intent(out):: Mol
      integer, intent(in):: nn

      ! Local

      integer:: ios


      ! Routine name
      urou = 'allocatemol'

      ! If no molecules
      if (nn.lt.1) then

        ! Allocate 1 at least
        allocate(Mol(1), stat=ios)
        MRAMc = MRAMc + 1d-6*sizeof(Mol(1))

      ! Proper molecules
      else

        ! Allocate molecules to read
        allocate(Mol(nn), stat=ios)

      end if

      ! Control that everything went fine
      call control

      return

      end subroutine allocatemol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate the molecules number density and initialize them to
      !! zero\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!        nn(integer): Size of molecule array
      subroutine preparemol(Mol,nn)

      ! I/O

      type(Mol_class), dimension(:), intent(inout):: Mol
      integer, intent(in):: nn

      ! Local

      integer:: im


      ! For each molecule
      do im=1,nn

        ! Allocate space por populations
        allocate(Mol(im)%n(NZ))
        MRAMc = MRAMc + 1d-6*sizeof(Mol(im)%n)

        ! Initialize populations to zero
        Mol(im)%n = 0d0

      end do ! Molecules

      end subroutine preparemol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the equilibrium constant of all molecules in a
      !! given model atmosphere\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!        nn(integer): Size of molecule array\n
      !!   Atmo(Atmo_class): Structure with atmospheric data
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

        ! If allocated
        if (allocated(Mol(im)%eq)) then

          ! If wrong size
          if (size(Mol(im)%eq).ne.NZ) then

            ! Reallocate
            MRAMc = MRAMc - 1d-6*sizeof(Mol(im)%eq)
            deallocate(Mol(im)%eq)
            allocate(Mol(im)%eq(NZ))
            MRAMc = MRAMc + 1d-6*sizeof(Mol(im)%eq)

          end if ! Wrong size

        ! Not allocated
        else

          ! Allocate
          allocate(Mol(im)%eq(NZ))

        end if ! Allocated equilibrium constant array

        ! Calculate equilibrium constant for temperature grid
        call Moleq(Mol(im),Atmo%T)

        ! If there are saved coefficients and we can free them
        if (Mol(im)%neqcoeff.gt.0.and.free) then

          ! Free tabulation
          MRAMc = MRAMc - 1d-6*sizeof(Mol(im)%eqcoeff)
          deallocate(Mol(im)%eqcoeff)

        end if ! Free tabulation

      end do ! Molecules

      end subroutine setupmol_eq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the partition function of all molecules in a given
      !! model atmosphere\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!        nn(integer): Size of molecule array\n
      !!   Atmo(Atmo_class): Structure with atmospheric data
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

        ! If allocated
        if (allocated(Mol(im)%pf)) then

          ! If wrong size
          if (size(Mol(im)%pf).ne.NZ) then

            ! Reallocate
            MRAMc = MRAMc - 1d-6*sizeof(Mol(im)%pf)
            deallocate(Mol(im)%pf)
            allocate(Mol(im)%pf(NZ))
            MRAMc = MRAMc + 1d-6*sizeof(Mol(im)%pf)

          end if ! Wrong size

        ! Not allocated
        else

          ! Allocate
          allocate(Mol(im)%pf(NZ))

        end if ! Allocated equilibrium constant array

        ! Calculate partition function for temperature grid
        call Molpf(Mol(im),Atmo%T)

        ! If there are saved coefficients and we can free them
        if (Mol(im)%npfcoeff.gt.0.and.free) then

          ! Free tabulation
          MRAMc = MRAMc - 1d-6*sizeof(Mol(im)%pfcoeff)
          deallocate(Mol(im)%pfcoeff)

        end if ! Free tabulation

      end do ! Molecules

      end subroutine setupmol_pf

!#####################################################################
!#####################################################################
!#####################################################################

      end module rmol_mod
