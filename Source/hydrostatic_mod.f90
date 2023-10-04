      !> Hydrostatic equilibrium calculations
      module hydrostatic_mod
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
!     08/07/2023 V3.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.3 - Added dummy argument for compatibility
!                             with module (TdPA)
!
!     07/03/2023:    V3.0.2 - The model atmosphere now really works
!                             with T and Pgas (TdPA)
!                           - Use the new routines to allocate and
!                             deallocate atomic and molecular
!                             quantities (TdPA)
!                           - Use the new routines to initialize
!                             populations and call the eos and
!                             chemical equilibrium (TdPA)
!                           - Use the new chi_freq instead of the
!                             pseudo-copy of the last fragment
!                             of the old background routine which
!                             was in Compute_beta, removing the
!                             latter (TdPA)
!
!     03/15/2023:    V3.0.1 - Removed Flgsg as an argument from
!                             several routines (TdPA)
!                           - Removed unecessary broadcasts (TdPA)
!                           - Removed some commented blocks remaining
!                             from the original TIC (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/17/2023:    V0.0.0 - Started from 05/12/2020
!                             TIC@hydrostatic_mod.f90 revision from
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
!    Compute_Pressure_all:
!      Compute pressurefrom hydrostatic equilibrium
!
!    Fill_Atmo:
!      Put calculated node into the proper atmospheric model
!
!    Compute_Opacity:
!      Compute opacity at given point
!
!    Compute_population
!      Startup populations
!
!#####################################################################
!#####################################################################
!#####################################################################

      use background_mod
      use chemic_mod
      use commons_mod
      use free_mod
      use initmodel_mod
      use initpopu_mod
      use parameters_mod, only: vacuum
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the pressure asymming hydrostatic equilibrium\n
      !!      Atmo(Atmo_class): Structure with the model\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!     Atomb(Atom_class): Structure with the atomic data for
      !!                        background opacities\n
      !!        Mol(Mol_class): Structure with the molecule data\n
      !!    Input(Input_class): Structure with settings data\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!      Pg_input(double): Boundary pressure\n
      !! Reference:\n
      !! de la Cruz Rodríguez (2019), Mihalas (1970).
      subroutine Compute_Pressure_all(Atmo,Atom,Atomb,Mol,Input, &
                                      fudge,Pg_input)

      ! IO
      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(fudge_class):: fudge
      double precision, intent(in):: Pg_input

      ! Local
      type(Atmo_class):: Atmo_tmp

      integer:: local_nZ, ia, iter, maxiter

      double precision:: Pg, Pg_new, Beta_old, Beta_new, dtao, dif
      double precision:: limit, gravity


      ! Save true nz
      local_nZ = nZ

      ! Lie with nz
      nZ = 1

      ! Allocate variables in temporal atmosphere
      allocate(Atmo_tmp%ne(1),Atmo_tmp%T(1),Atmo_tmp%nh(1,6))
      allocate(Atmo_tmp%nHT(1),Atmo_tmp%nHm(1),Atmo_tmp%nHa(1))
      allocate(Atmo_tmp%rho(1),Atmo_tmp%Pg(1))
      nullify(Atmo_tmp%z,Atmo_tmp%vmi,Atmo_tmp%zeros)
      nullify(Atmo_tmp%Bx,Atmo_tmp%By,Atmo_tmp%Bz)
      nullify(Atmo_tmp%vx,Atmo_tmp%vy,Atmo_tmp%vz)
      nullify(Atmo_tmp%vxa,Atmo_tmp%vya,Atmo_tmp%vza)

      ! Setup temporal atmosphere defined by the gas pressure
      ! typo = 4
      Atmo_tmp%alloc_a = .True.
      Atmo_tmp%alloc_b = .False.
      Atmo_tmp%tfreq = Atmo%tfreq
      Atmo_tmp%typo = 4
      Atmo_tmp%Pg(1) = Pg_input
      Pg = Pg_input
      Atmo_tmp%T(1) = Atmo%T(1)
      gravity = 10d0**(Atmo%logg)
      Atmo_tmp%nh = 0d0
      Atmo_tmp%nHT = 0d0
      Atmo_tmp%nHm = 0d0
      Atmo_tmp%nHa = 0d0
      Atmo_tmp%ele = Atmo%ele
      Atmo_tmp%nele = Atmo%nele
      Atmo_tmp%NT = Atmo%NT
      Atmo_tmp%pT = Atmo%pT
      Atmo_tmp%logg = Atmo%logg
      Atmo_tmp%tfreq = Atmo%tfreq
      Atmo_tmp%ele = Atmo%ele
      Atmo_tmp%scal = Atmo%scal
      Atmo_tmp%abund = Atmo%abund
      Atmo_tmp%nZ = 1

      ! Fake total population size in active atoms
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! Compute opacity at current point
      call Compute_Opacity(Atmo_tmp, Atom, Atomb, Mol, &
                           Input, fudge, Beta_old, Pg)

      ! Dump data into boundary node
      call Fill_Atmo(Atmo, Atmo_tmp, 1)

      ! Initialize parameters
      limit = 1d-2
      maxiter = 20

      ! For the rest of heights
      do ia=2,local_nZ

        ! Step in vertical coordinate
        dtao = Atmo%z(ia) - Atmo%z(ia-1)

        ! Initialize diff
        dif = 1d0

        ! Set-up new auxiliar node
        Atmo_tmp%Pg(1) = Pg + gravity*dtao/Beta_old
        Atmo_tmp%T(1) = Atmo%T(ia)
        Atmo_tmp%nh = 0d0
        Atmo_tmp%nHT = 0d0
        Atmo_tmp%nHm = 0d0
        Atmo_tmp%nHa = 0d0

        ! Compute opacity at current point
        call Compute_Opacity(Atmo_tmp, Atom, Atomb, Mol, &
                             Input, fudge, Beta_new, Pg_new)

        ! Iterate
        do iter=1,maxiter

          ! Initialize diff (again?)
          dif = Pg_new

          ! Iterate density
          Atmo_tmp%Pg(1) = Pg + &
                           gravity*dtao/(Beta_new-Beta_old)* &
                           log10(Beta_new/Beta_old)

          ! Clean atmosphere
          Atmo_tmp%nh = 0d0
          Atmo_tmp%nHT = 0d0
          Atmo_tmp%nHm = 0d0
          Atmo_tmp%nHa = 0d0

          ! Compute opacity
          call Compute_Opacity(Atmo_tmp, Atom, Atomb, Mol, &
                               Input, fudge, Beta_new, Pg_new)

          ! Compute difference
          dif = abs((dif-Pg_new)/(dif+Pg_new))

          ! If converged, leave
          if (dif.le.limit) exit

        end do ! Iterations

        ! Shift results
        Pg = Pg_new
        Beta_old = Beta_new

        ! Store in atmospheric model
        call Fill_Atmo(Atmo,Atmo_tmp,ia)

      end do ! Heights

      ! Clean temporal model
      call free_Atmo(Atmo_tmp,.True.)

      ! Recover actual stratification
      nZ = local_nZ

      ! Free global populations
      call free_gpop(Atom,Atomb,Mol)

      ! Control
      call control

      return

      end subroutine Compute_Pressure_all

!#####################################################################
!#####################################################################
!#####################################################################

      !> Fill actual model with auxiliar node\n
      !!      Atmo(Atmo_class): Structure with the model\n
      !!  Atmo_tmp(Atmo_class): Structure with the auxiliar model\n
      !!         indx(integer): Height index where to store
      subroutine Fill_Atmo(Atmo,Atmo_tmp,indx)

      ! IO
      type(Atmo_class), intent(inout):: Atmo
      type(Atmo_class), intent(inout):: Atmo_tmp
      integer:: indx

      Atmo%ne(indx) = Atmo_tmp%ne(1)
      Atmo%Pg(indx) = Atmo_tmp%Pg(1)
      Atmo%nHT(indx) = Atmo_tmp%nHT(1)
      Atmo%nHm(indx) = Atmo_tmp%nHm(1)
      Atmo%nh(indx,:) = Atmo_tmp%nh(1,:)
      Atmo%nHa(indx) = Atmo_tmp%nHa(1)

      end subroutine Fill_Atmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the opacity at current point setup in Atmo\n
      !!    Atmo(Atmo_class): Structure with the model\n
      !!    Atom(Atom_class): Structure with the atomic data\n
      !!   Atomb(Atom_class): Structure with the atomic data for
      !!                      background opacities\n
      !!      Mol(Mol_class): Structure with the molecule data\n
      !!  Input(Input_class): Structure with settings data\n
      !!  fudge(fudge_class): Structure with fudge data\n
      !!        beta(double): Beta factor\n
      !!          Pg(double): Gas pressure
      subroutine Compute_Opacity(Atmo,Atom,Atomb,Mol,Input,fudge, &
                                 Beta,Pg)

      ! IO
      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(fudge_class):: fudge
      double precision, intent(inout):: Beta, Pg

      ! Local
      integer, dimension(Atmo%nele):: nlte, depar
      double precision, dimension(1):: aBeta

      ! Initialize beta
      Beta = 0d0

      ! Initialize nlte
      nlte = 0
      depar = 0

      ! Solve equation of state
      call eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! Compute populations
      call Compute_population(Atmo,Atom,Atomb,Mol,Input)

      ! Compute beta factor
      call chi_freq(Atom,Atomb,Mol,Atmo,fudge,Input, &
                    Atmo%tfreq,aBeta,1,1,.False.)

      ! Iterate
      Beta = aBeta(1)/Atmo%rho(1)
      Pg = Atmo%Pg(1)

      ! Free level populations
      call free_lpop(Atom,Atomb)

      return

      end subroutine Compute_Opacity

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute LTE populations\n
      !!      Atmo(Atmo_class): Structure with the model\n
      !!      Atom(Atom_class): Structure with the atomic data\n
      !!     Atomb(Atom_class): Structure with the atomic data for
      !!                        background opacities\n
      !!        Mol(Mol_class): Structure with the molecule data\n
      !!    Input(Input_class): Structure with settings data
      subroutine Compute_population(Atmo,Atom,Atomb,Mol,Input)

      ! IO
      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input

      ! Local
      type(LTEline_class), dimension(:), allocatable:: dummy
      integer:: ia

      !
      ! Set LTE
      call setlte(Atom,Atomb,Atmo,Input)

      !
      ! Set populations
      call setuppopu(Atom,Atomb,Atmo)

      ! Initialize molecule populations
      do ia=1,nM
        Mol(ia)%n = 0d0
      end do

      ! Compute chemical equilibrium
      call chemeq(Atom,Atomb,dummy,Mol,Atmo)

      return

      end subroutine Compute_population

!#####################################################################
!#####################################################################
!#####################################################################

      end module hydrostatic_mod
