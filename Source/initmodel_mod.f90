      !> Initialization of atoms, molecules, chemical eq., etc.
      module initmodel_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     06/16/2023
!  Last version:
!     07/18/2024 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     07/18/2024:    V3.0.5 - Added call to write elastic rates (TdPA)
!
!     10/16/2023:    V3.0.4 - The damping coefficients are now
!                             calculated in the new broadening
!                             routines (TdPA)
!                           - Added broadening and broadening_line
!                             subroutines (TdPA)
!                           - The density matrix is initialized
!                             via calls from the new Initrhoes
!                             routines (TdPA)
!                           - Ensure memory deallocation when an
!                             error happens in prepare_syn (TdPA)
!                           - Updated arguments of free_local_Atom
!                             call (TdPA)
!
!     09/29/2023:    V3.0.3 - Added arguments to Initcols (TdPA)
!
!     08/22/2023:    V3.0.2 - Changed the optical depth correction
!                             when checking the input atmospheric
!                             model when the top boundary has 0
!                             optical depth (TdPA)
!
!     08/07/2023:    V3.0.1 - Added setlte_lines (TdPA)
!                           - Added arguments for LTE lines (TdPA)
!
!     07/03/2023:    V3.0.0 - First version (TdPA)
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
!    Routines to initialize atoms, molecules (total densities) and
!  to apply the eq. of state and chemical equilibrium, setting up
!  the atmospheric model
!
!  prepareatomol:
!    Call for the allocation of the atom and molecule structures.
!
!  setpopufiles:
!    Call for the reading of population files and the set-up of
!  nlte populations or departure coefficients.
!
!  reviseH_init:
!    Revise the H populations in the atmosphere or in the model atom
!  depending on the inputs.
!
!  reviseH_out:
!    Revise the H populations in the atmosphere or in the model atom
!  depending on the inputs, but after having running the eos.
!
!  setlte:
!    Initialize LTE populations in atoms.
!
!  setlte_lines:
!    Initialize populations of the LTE lines.
!
!  setcols:
!    Initialize collisional rates in atoms and call for its output.
!
!  setuppopu:
!    Initialize populations and density matrices in atoms.
!
!  broadening:
!    Calculate the line broadening damping coefficient.
!
!  broadening_line:
!    Calculate the line broadening damping coefficient for LTE lines.
!
!  chemical:
!    Manage the call to chemical equilibrium and revise model
!  atmosphere Hydrogen number density.
!
!  Initrhoes:
!    Manage the initialization of the density matrix.
!
!  updateatmo:
!    Recalculate electron number density or write atmospheric model
!  into file.
!
!  prepare_syn:
!    Prepare an input model atmosphere for a forward solution.
!
!  setup_Atmo_ininv:
!    Prepare the initial atmospheric model in an inversion to put it
!  in the common expected form.
!
!  setup_Atmo_ouinv:
!    Prepare the final atmospheric model in an inversion to write into
!  the Result file.
!
!#####################################################################
!#####################################################################
!#####################################################################

      use background_mod
      use broad_mod
      use chemic_mod
      use commons_mod
      use free_mod
      use getztau_mod
      use initcols_mod
      use initpopu_mod
      use iosolution_mod
      use ratmo_mod
      use ratom_mod
      use rmol_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate space for atoms and molecules\n
      !!    Atom(Atom_class): Structure with the atomic data\n
      !!   Atomb(Atom_class): Structure with the atomic data for
      !!                      background opacities\n
      !!      Mol(Mol_class): Structure with the molecule data\n
      !!         nn(integer): Number of molecules
      subroutine prepareatomol(Atom,Atomb,Mol,nn)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      integer, intent(in):: nn

      ! Prepare active atoms
      call prepareatom(Atom,nA)
      ! Prepare passive atoms
      if (nAb.gt.0) call prepareatom(Atomb,nAb)

      ! Prepare molecules
      if (nn.gt.0) call preparemol(Mol,nn)

      end subroutine prepareatomol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere dealing with the equation of
      !! state\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Input(Input_class): Structure with settings data\n
      !!     nlte(integer(:)): Array with information about
      !!                       loaded populations\n
      !!    depar(integer(:)): Array with information about
      !!                       loaded departure coefficients
      subroutine setpopufiles(Atom,Atomb,Atmo,Input,nlte,depar)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo
      type(Input_class):: Input
      integer, dimension(:), allocatable, intent(out):: nlte,depar

      ! Local
      integer:: ia

      !
      ! Read NLTE populations
      !

      ! For every active atom
      do ia=1,nA
        call Initpopu_file(Atom(ia),Input%popu(ia)%str,Atmo)
      end do

      ! Control
      if (laborted) return

      ! For every passive atom
      do ia=1,nAb
        call Initpopu_file(Atomb(ia),Input%popuback(ia)%str,Atmo)
      end do

      ! Control
      if (laborted) return

      !
      ! Initialize nlte
      !
      call initializenlte(nlte,depar,Atom,Atomb,Atmo)

      end subroutine setpopufiles

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if populations of H atom are read\n
      !! state\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data
      subroutine reviseH_init(Atom,Atomb,Atmo)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo

      ! Local
      logical:: l1

      integer:: ia


      !
      ! Check the hydrogen populations in the model atmosphere
      ! if H populations are given
      !

      ! Check if there are hydrogen populations already set

      ! Initialize logical control variables
      l1 = .False.

      ! For every atom
      do ia=1,nA

        ! If it is hydrogen
        if (Atom(ia)%element.eq.' H') then

          ! Flag hydrogen found
          l1 = .True.

          ! Check if NLTE populations
          if (allocated(Atom(ia)%popu)) then

            ! Revise populations and exit
            call ReviseHatmo(Atom(ia),Atmo)
            exit

          end if ! If input file
        end if ! If Hydrogen
      end do

      ! If hydrogen not active
      if (.not.l1) then

        ! For every passive atom
        do ia=1,nAb

          ! If it is hydrogen
          if (Atomb(ia)%element.eq.' H') then

            ! Flag hydrogen found
            l1 = .True.

            ! Check if NLTE populations
            if (allocated(Atomb(ia)%popu)) then

              ! Revise populations and exit
              call ReviseHatmo(Atomb(ia),Atmo)
              exit

            end if ! If input file
          end if ! If Hydrogen

        end do ! passive atoms

      end if ! Hydrogen not found between active atoms

      end subroutine reviseH_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere dealing with the equation of
      !! state\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Input(Input_class): Structure with settings data
      subroutine reviseH_out(Atom,Atomb,Atmo,Input)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo
      type(Input_class):: Input

      ! Local
      integer:: ia

      ! Check active atoms
      do ia=1,nA
        if (Atom(ia)%element.eq.' H') then
          Atom(ia)%mol_protect = Input%protect_H
          if (allocated(Atom(ia)%popu)) then
            if (Input%protect_H) then
              call ReviseHatmo(Atom(ia),Atmo)
            else
              call ReviseHatom(Atom(ia),Atmo)
            end if
          end if
          exit
        end if
      end do

      ! Check passive atoms
      do ia=1,nAb
        if (Atomb(ia)%element.eq.' H') then
          Atomb(ia)%mol_protect = Input%protect_H
          if (allocated(Atomb(ia)%popu)) then
            if (Input%protect_H) then
              call ReviseHatmo(Atomb(ia),Atmo)
            else
              call ReviseHatom(Atomb(ia),Atmo)
            end if
          end if
          exit
        end if
      end do

      end subroutine reviseH_out

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere dealing with the equation of
      !! state\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Input(Input_class): Structure with settings data
      subroutine setlte(Atom,Atomb,Atmo,Input)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo
      type(Input_class):: Input

      ! Local
      integer:: ia

      ! Initialize LTE populations

      ! Active atoms
      do ia=1,nA
        call Initpopu_LTE(Atom(ia),Atmo,.True.,0)
      end do

      ! Control
      if (laborted) return

      ! Passive atoms
      do ia=1,nAb
        call Initpopu_LTE(Atomb(ia),Atmo,.False.,0)
      end do

      end subroutine setlte

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere dealing with the equation of
      !! state\n
      !!     Atomb(Atom_class): Structure with the atomic data for
      !!                        background opacities\n
      !! LTElines(Input_class): Structure with LTE line data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with settings data
      subroutine setlte_lines(Atomb,LTElines,Atmo,Input)

      ! I/O
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo
      type(Input_class):: Input

      ! Local
      integer:: ia,jdir,iz

      ! Copy from input
      LTElines = Input%LTEline

      ! For each LTE line
      do ia=1,nLTEl

        ! Nullify pointer
        nullify(LTElines(ia)%prof)

        ! Check Normp is not allocated
        if (associated(LTElines(ia)%prof)) then
          do jdir=1,size(LTElines(ia)%prof,2)
            do iz=LTElines(ia)%Rz0, &
                  size(LTElines(ia)%prof,1)+LTElines(ia)%Rz0
              if (allocated(LTElines(ia)%prof(iz,jdir)%p)) &
                deallocate(LTElines(ia)%prof(iz,jdir)%p)
              if (allocated(LTElines(ia)%prof(iz,jdir)%comp)) &
                deallocate(LTElines(ia)%prof(iz,jdir)%comp)
            end do
          end do
          deallocate(LTElines(ia)%prof)
        end if

        ! Nullify pointer
        nullify(LTElines(ia)%prof)

        ! Get LTE populations
        call Initpopu_LTE_line(LTElines(ia),Atmo)

        ! If passive
        if (LTElines(ia)%is_passive) then

          ! Get abundance from passive
          LTElines(ia)%n = LTElines(ia)%n*Atomb(LTElines(ia)%ia)%abun

        ! Not passive
        else

          ! Get abundance from atmosphere
          LTElines(ia)%n = LTElines(ia)%n*Atmo%abund(LTElines(ia)%ele)

        end if

      end do

      end subroutine setlte_lines

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere dealing with the equation of
      !! state\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Input(Input_class): Structure with settings data\n
      !!   Flgsg(Fctsg_class): Structure with factorials and
      !!                       signs
      subroutine setcols(Atom,Atomb,Atmo,Input,Flgsg)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input

      ! Local
      integer:: ia


      ! Active atoms
      do ia=1,nA
        call Initcols(Atom(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.True.)
      end do

      ! Control
      if (laborted) return

      ! If keeping collisions and master, call writer
      if (Input%keep_cols) then
        call writecols(Atom,Input%folder,&
                       Input%lim_cols_tt,Input%lim_cols_ll)
      end if

      ! Control
      if (laborted) return

      ! Passive atoms
      do ia=1,nAb
        call Initcols(Atomb(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.False.)
      end do

      end subroutine setcols

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize populations\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data
      subroutine setuppopu(Atom,Atomb,Atmo)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo

      ! Local
      integer:: ia

      !
      ! Initialize populations properly
      !

      ! Active atoms
      do ia=1,nA
        call Initpopu(Atom(ia),Atmo,.True.)
      end do

      ! Control
      if (laborted) return

      ! Passive atoms
      do ia=1,nAb
        call Initpopu(Atomb(ia),Atmo,.False.)
      end do

      end subroutine setuppopu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate line broadening\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!       Atomb(Atom_class): Structure with the atomic data for
      !!                          background opacities\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Input(Input_class): Structure with settings data
      subroutine broadening(Atom,Atomb,Atmo,Input)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo
      type(Input_class):: Input

      ! Local
      integer:: ia

      ! For active atoms
      do ia=1,nA
        call broad(Atom(ia),Atmo,Input%folder,Input%keep_aparam)
      end do

      ! Control
      if (laborted) return

      ! If keeping damping and master, call writer
      if (Input%keep_damp) call writedamp(Atom,Atmo,Input%folder, &
                                          Input%lim_damp)

      ! If keeping elastic rates and master, call writer
      if (Input%keep_qel) call writeqel(Atom,Atmo,Input%folder, &
                                        Input%lim_qel)

      ! Control
      if (laborted) return

      ! And for background atoms
      if (Input%nAb.gt.0) then
        do ia=1,Input%nAb
          call broad(Atomb(ia),Atmo,Input%folder,.False.)
        end do
      end if

      end subroutine broadening

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate line broadening\n
      !! LTElines(LTEline_class): Structure with the LTE line data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data
      subroutine broadening_line(LTElines,Atmo)

      ! I/O
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class):: Atmo

      ! Local
      integer:: ia

      ! If LTE lines
      do ia=1,nLTEl
        call broad_line(LTElines(ia),Atmo)
      end do

      end subroutine broadening_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere dealing with the equation of
      !! state\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!       Atomb(Atom_class): Structure with the atomic data for
      !!                          background opacities\n
      !! LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Mol(Mol_class): Structure with the molecule data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Input(Input_class): Structure with settings data
      subroutine chemical(Atom,Atomb,LTEline,Mol,Atmo,Input)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(LTEline_class), dimension(:), allocatable:: LTEline
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(Input_class):: Input

      ! Local
      integer:: ia

      !
      ! Protect atoms from chemical equilibrium
      !
      if (Input%chem_protect_all) then
        ! Active
        do ia=1,nA
          Atom(ia)%mol_protect = .True.
        end do
        ! Passive
        do ia=1,nAb
          Atomb(ia)%mol_protect = .True.
        end do
      end if

      !
      ! Calculate chemical equilibrium
      !
      call chemeq(Atom,Atomb,LTEline,Mol,Atmo)

      ! Control
      if (laborted) return

      !
      ! Check again atmospheric hydrogen number density
      !

      ! Check active atoms
      do ia=1,nA
        if (Atom(ia)%element.eq.' H') then
          if (allocated(Atom(ia)%popu)) then
            call ReviseHatmo(Atom(ia),Atmo)
          end if
          exit
        end if
      end do

      ! Check passive atoms
      do ia=1,nAb
        if (Atomb(ia)%element.eq.' H') then
          if (allocated(Atomb(ia)%popu)) then
            call ReviseHatmo(Atomb(ia),Atmo)
          end if
          exit
        end if
      end do

      end subroutine chemical

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the density matrices\n
      !!     Atom(Atom_class): Structure with the atomic data
      subroutine Initrhoes(Atom)

      ! I/O
      type(Atom_class), dimension(:):: Atom

      ! Local
      integer:: ia


      ! Active atoms
      do ia=1,nA
        call Initcrho(Atom(ia))
      end do

      end subroutine Initrhoes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Update model atmosphere\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !! Bfield(Bfield_class): Structure with magnetic field
      !!                       data\n
      !!   Input(Input_class): Structure with settings data
      subroutine updateatmo(Atom,Atomb,Atmo,Bfield,Input)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Input_class):: Input

      ! Local
      integer:: ia
      integer, dimension(:), allocatable:: nlte,depar


      !
      ! Check if updating the model
      !

      if (Input%update_atmos.lt.0) return

      !
      ! Revise nlte
      !
      call activenlte(nlte,depar,Atom)

      !
      ! Revise the Hydrogen if it was active
      !

      ! For every atom
      do ia=1,nA

        ! Correct populations
        call correctpop(Atom(ia),1)

        ! If it is hydrogen
        if (Atom(ia)%element.EQ.' H') then

          ! Revise populations and exit
          call ReviseHatmo(Atom(ia),Atmo)

        end if ! If Hydrogen

      end do ! Go through all active atoms

      !
      ! Redo the electrons if specified in input
      !
      if (Input%redo_ne.eq.1.or.Input%redo_ne.eq.11) &
        call redo_ne(Atom,Atomb,nlte,depar,Atmo)

      if (Input%keep_atmo) call writeatmo(Atmo,Bfield, &
                                          Input%folder, &
                                          Input%lim_atmo)

      !
      ! Output the new atmosphere if 1D
      !
      if (run_mode.eq.0) &
        call wAtmo(Atmo,Input%update_atmos,Input%folder, &
                   Input%atmo)

      return

      end subroutine updateatmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare for the synthesis\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!            Mol(Mol_class): Structure with the molecule data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Input(Input_class): Structure with settings data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      subroutine prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(Fctsg_class):: Flgsg
      type(Input_class):: Input

      ! Local
      integer, dimension(:), allocatable:: nlte,depar


      ! Prepare atoms and molecules
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! If error, skip
      if (laborted) goto 1000

      ! Check populations from files
      call setpopufiles(Atom,Atomb,Atmo,Input,nlte,depar)

      ! If error, skip
      if (laborted) goto 1000

      ! Revise H population
      call reviseH_init(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Eq of state
      if (Atmo%typo.gt.0.or.Input%keep_atmo.or. &
          Input%redo_ne.gt.0) then

        ! Call equation of state
        call eqstate(Atmo,Atom,Atomb,nlte,depar)

        ! Control
        if (laborted) goto 1000

        ! Revise H population
        if (Atmo%typo.gt.0) &
          call reviseH_out(Atom,Atomb,Atmo,Input)

      end if

      ! If error, skip
      if (laborted) goto 1000

      !
      ! Recalculate electron density
      !
      if (Input%redo_ne.ge.10) &
        call redo_ne(Atom,Atomb,nlte,depar,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! LTE populations
      call setlte(Atom,Atomb,Atmo,Input)

      ! If error, skip
      if (laborted) goto 1000

      ! Collisions
      call setcols(Atom,Atomb,Atmo,Input,Flgsg)

      ! If error, skip
      if (laborted) goto 1000

      ! LTE lines
      if (nLTEl.gt.0) then

        ! Prepare LTE lines
        call setlte_lines(Atomb,LTElines,Atmo,Input)

      end if

      ! Initialize populations
      call setuppopu(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Broadening
      call broadening(Atom,Atomb,Atmo,Input)

      ! If error, skip
      if (laborted) goto 1000

      ! LTE lines
      if (nLTEl.gt.0) call broadening_line(LTElines,Atmo)

      ! Chemical equilibrium
      call chemical(Atom,Atomb,LTElines,Mol,Atmo,Input)

      ! If error, free atom memory
      if (laborted) goto 1000

      ! Initialize density matrix
      call Initrhoes(Atom)

      return

      ! Free memory if aborting
1000  call free_lpop(Atom,Atomb)
      call free_gpop(Atom,Atomb,Mol)

      return

      end subroutine prepare_syn

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up input model atmosphere for inversion\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!       Mol(Mol_class): Structure with the molecule data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!      MPID(MPI_class): Structure with MPI data\n
      !!   Input(Input_class): Structure with settings data\n
      !!   Flgsg(Fctsg_class): Structure with factorials and
      !!                       signs\n
      !!   fudge(fudge_class): Structure with fudge data\n
      !!      zalt(double(:)): Alternative height axis\n
      !!       alloc(logical): If tau axis must be allocated
      subroutine setup_Atmo_ininv(Atom,Atomb,Mol,Atmo,MPID,Input, &
                                  fudge,zalt,alloc)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(fudge_class):: fudge
      type(Input_class):: Input
      type(MPI_class):: MPID
      logical, intent(in):: alloc
      double precision, dimension(:), allocatable, target:: zalt

      ! Local
      type(Atmo_class):: Atmo_tmp
      type(LTEline_class), dimension(:), allocatable:: dummy

      integer:: ii
      integer, dimension(Atmo%nele):: nlte,depar


      ! Keep current nz
      ii = nz
      nz = Atmo%nz

      ! No data
      nlte = 0
      depar = 0

      ! Prepare atoms and molecules
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! If error, skip
      if (laborted) return

      ! Revise H populations
      call reviseH_init(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) return

      ! Eq. of state
      call eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! If error, skip
      if (laborted) return

      ! If not tau scale
      if (Atmo%scal.ne.'T') then

        ! Get a copy of Atmo
        call cAtmo(Atmo,Atmo_tmp)

        ! Initialize LTE populations
        call setlte(Atom,Atomb,Atmo,Input)

        ! Initialize populations
        call setuppopu(Atom,Atomb,Atmo_tmp)

        ! Chemical equilibrium
        call chemical(Atom,Atomb,dummy,Mol,Atmo_tmp,Input)
        if (laborted) return

        ! Allocate chi500
        if (.not.allocated(Atmo_tmp%chi500)) &
          allocate(Atmo_tmp%chi500(nz))


        !
        ! Calculate continuum opacity at reference frequency
        !
        call chi_freq(Atom,Atomb,Mol,Atmo_tmp,fudge,Input, &
                      Atmo_tmp%tfreq,Atmo_tmp%chi500,1,nz, &
                      MPID%mpi)

        ! Free
        call free_lpop(Atom,Atomb)
        call free_mol(Mol)

        !
        ! Compute missing height or tau
        !
        call getztau(Atmo_tmp,MPID,.False.)


        !
        ! Get the tau scale
        !

        ! If alloc
        if (alloc) then

          Atmo%z = Atmo_tmp%zalt

        ! Not alloc
        else

          if (.not.allocated(zalt)) allocate(zalt(nz))
          zalt = Atmo_tmp%zalt
          Atmo%z => zalt

        end if

        !
        ! Wipe the copy clean
        call free_Atmo(Atmo_tmp,.True.)

      end if ! No tau scale

      ! Atmo typo is for sure 4 now, and scale tau
      Atmo%typo = 4
      Atmo%scal = 'T'
      ztau = .True.

      ! Reset variables
      Atmo%ne = 0d0
      Atmo%nh = 0d0
      Atmo%nht = 0d0
      Atmo%nha = 0d0
      Atmo%nhm = 0d0
      if (allocated(Atmo%rho)) deallocate(Atmo%rho)
      if (allocated(Atmo%Pe)) deallocate(Atmo%Pe)

      ! Restore nz
      nz = ii

      ! Sanity check tau
      if (minval(Atmo%z).le.0d0) then

        ! If really negative
        if (minval(Atmo%z).lt.0d0) then

          umsg = 'The input model atmosphere has '// &
                 'negative optical depths'
          urou = 'setup_Atmo'
          call aborted
          return

        end if

        ! Try to follow some progression
        Atmo%z(1) = Atmo%z(2)*Atmo%z(2)/Atmo%z(3)

        ! If still negative, try with the largest step
        if (Atmo%z(1).le.0d0) &
          Atmo%z(1) = 2d0*Atmo%z(2) - Atmo%z(3)

        ! If still zero or negative, put the upper
        ! boundary one orders of magnitude smaller
        if (Atmo%z(1).le.0d0) Atmo%z(1) = Atmo%z(2)*1d-1

        ! For consistency, add to all the depths
        Atmo%z(2:Atmo%nz) = Atmo%z(2:Atmo%nz) + Atmo%z(1)

        ! If still 0, give up
        if (minval(Atmo%z).le.0d0) then

          umsg = 'The input model atmosphere has '// &
                 'a non valid optical depth stratification'
          urou = 'setup_Atmo'
          call aborted

        end if ! Still 0 optical depth
      end if ! 0 or negative optical depth

      ! Free
      call free_local_Atom(Atom)
      call free_gpop(Atom,Atomb,Mol)


      end subroutine setup_Atmo_ininv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up output model atmosphere for inversion\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atomb(Atom_class): Structure with the atomic data for
      !!                       background opacities\n
      !!       Mol(Mol_class): Structure with the molecule data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      !!      MPID(MPI_class): Structure with MPI data\n
      !!   Input(Input_class): Structure with settings data\n
      !!   Flgsg(Fctsg_class): Structure with factorials and
      !!                       signs\n
      !!   fudge(fudge_class): Structure with fudge data
      subroutine setup_Atmo_ouinv(Atom,Atomb,Mol,Atmo,MPID,Input, &
                                  fudge)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(fudge_class):: fudge
      type(Input_class):: Input
      type(MPI_class):: MPID

      ! Local
      type(LTEline_class), dimension(:), allocatable:: dummy

      integer, dimension(Atmo%nele):: nlte,depar

      ! No data
      nlte = 0
      depar = 0

      ! Prepare atoms and molecules
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! If error, skip
      if (laborted) return

      ! Revise H populations
      call reviseH_init(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) return

      ! Eq. of state
      call eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! If error, skip
      if (laborted) return

      ! Initialize LTE populations
      call setlte(Atom,Atomb,Atmo,Input)

      ! Initialize populations
      call setuppopu(Atom,Atomb,Atmo)

      ! Chemical equilibrium
      call chemical(Atom,Atomb,dummy,Mol,Atmo,Input)
      if (laborted) return

      ! Allocate chi500
      if (.not.allocated(Atmo%chi500)) &
        allocate(Atmo%chi500(nz))

      !
      ! Calculate continuum opacity at reference frequency
      !
      call chi_freq(Atom,Atomb,Mol,Atmo,fudge,Input, &
                    Atmo%tfreq,Atmo%chi500,1,nz,MPID%mpi)

      ! Free
      call free_lpop(Atom,Atomb)
      call free_mol(Mol)
      call free_local_Atom(Atom)
      call free_gpop(Atom,Atomb,Mol)


      end subroutine setup_Atmo_ouinv

!#####################################################################
!#####################################################################
!#####################################################################

      end module initmodel_mod
