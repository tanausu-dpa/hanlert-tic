      !> Initialization of atoms, molecules, chemical eq., etc.
      module initmodel_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     16/06/2023
!  Last version:
!     04/11/2025 V4.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     04/11/2025:    V4.0.5 - Bugfix: in setup_Atmo_ininv, an error
!                             in optical depth setup could make the
!                             routine to return prematurely (not
!                             freeing all memory) and thus producing
!                             an even worse error (TdPA)
!                           - Changed atmosphere in argument for
!                             setlte in the initialization of the
!                             model atmosphere for inversions (TdPA)
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
!  prepareatomol
!    Allocate space needed in atoms and molecules to prepare the
!  model atmophere for the synthesis
!
!  setpopufiles
!    Try to initialize the populations of atoms from existing files
!  and indicate for which atoms these exist and in which format, i.e.,
!  populations or departure coefficients
!
!  reviseH_init
!    Revise the hydrogen populations in the model atmosphere if NLTE
!  populations were read from a file
!
!  reviseH_out:
!    Revise the hydrogen  populations in the model atmosphere or in
!  the model atom depending on the inputs once the EoS has been solved
!
!  setlte:
!    Initialize LTE populations in all atoms
!
!  setlte_lines
!    Initialize LTE population of the ion of the LTE lines
!
!  setcols
!    Initialize the collisional rates in atoms and write the output
!
!  setuppopu
!    Initialize populations and density matrices for all atoms
!
!  broadening:
!    Calculate the line broadening damping coefficient
!
!  broadening_line
!    Calculate the line broadening damping coefficient for LTE lines
!
!  chemical
!    Solve the chemical equilibrium and revise model the hydrogen
!  number density in the model atmosphere
!
!  Initrhoes
!    Initialize the density matrix of all atoms
!
!  updateatmo
!    Update the model atmosphere, recalculating number densities if
!  necessary
!
!  prepare_syn
!    Prepare a model atmosphere for the hanle routine to solve the
!  self-consistent problem and/or calculate formal solutions
!
!  setup_Atmo_ininv
!    Prepare a model atmosphere in an inversion to ensure that it is
!  in the expected format
!
!  setup_Atmo_ouinv
!    Prepare the atmospheric model resulting from a inversion to write
!  in the result file
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

      !> Allocate space needed in atoms and molecules to prepare the
      !! model atmophere for the synthesis\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!           nn(integer): Number of molecules
      subroutine prepareatomol(Atom,Atomb,Mol,nn)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      integer, intent(in):: nn

      ! Prepare active atoms
      if (nA.gt.0) call prepareatom(Atom,nA)

      ! Prepare background atoms
      if (nAb.gt.0) call prepareatom(Atomb,nAb)

      ! Prepare molecules
      if (nn.gt.0) call preparemol(Mol,nn)

      end subroutine prepareatomol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Try to initialize the populations of atoms from existing
      !! files and indicate for which atoms these exist and in which
      !! format, i.e., populations or departure coefficients\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!      nlte(integer(:)): Array with information about
      !!                        loaded populations\n
      !!     depar(integer(:)): Array with information about
      !!                        loaded departure coefficients
      subroutine setpopufiles(Atom,Atomb,Atmo,Input,nlte,depar)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input
      integer, dimension(:), allocatable, intent(inout):: nlte,depar

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

      ! For every background atom
      do ia=1,nAb
        call Initpopu_file(Atomb(ia),Input%popuback(ia)%str,Atmo)
      end do

      ! Control
      if (laborted) return

      ! Initialize nlte and depar variables
      call initializenlte(nlte,depar,Atom,Atomb,Atmo)

      end subroutine setpopufiles

!#####################################################################
!#####################################################################
!#####################################################################

      !> Revise the hydrogen populations in the model atmosphere if
      !! NLTE populations were read from a file\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data
      subroutine reviseH_init(Atom,Atomb,Atmo)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(Atmo_class), intent(inout):: Atmo

      ! Local

      logical:: l1

      integer:: ia


      ! Initialize logical control variables
      l1 = .False.

      ! For every active atom
      do ia=1,nA

        ! If it is hydrogen
        if (Atom(ia)%element.eq.' H') then

          ! Flag hydrogen found
          l1 = .True.

          ! Check if NLTE populations
          if (allocated(Atom(ia)%popu)) then

            ! Revise populations in model atmosphere and exit
            call ReviseHatmo(Atom(ia),Atmo)
            exit

          end if ! If input file
        end if ! If Hydrogen

      end do ! Active atoms

      ! If hydrogen not found among active atoms
      if (.not.l1) then

        ! For every background atom
        do ia=1,nAb

          ! If it is hydrogen
          if (Atomb(ia)%element.eq.' H') then

            ! Flag hydrogen found
            l1 = .True.

            ! Check if NLTE populations
            if (allocated(Atomb(ia)%popu)) then

              ! Revise populations in model atmosphere and exit
              call ReviseHatmo(Atomb(ia),Atmo)
              exit

            end if ! If input file
          end if ! If Hydrogen atom

        end do ! Background atoms

      end if ! Hydrogen not found between active atoms

      end subroutine reviseH_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Revise the hydrogen populations in the model atmosphere or
      !! in the model atom depending on the inputs once the EoS has
      !! been solved\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data
      subroutine reviseH_out(Atom,Atomb,Atmo,Input)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Atmo_class), intent(inout):: Atmo
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia


      ! For each active atom
      do ia=1,nA

        ! If it is hydrogen
        if (Atom(ia)%element.eq.' H') then

          ! Save if protecting
          Atom(ia)%mol_protect = Input%protect_H

          ! If allocated NLTE populations
          if (allocated(Atom(ia)%popu)) then

            ! If protected
            if (Input%protect_H) then

              ! Revise populations in model atmosphere
              call ReviseHatmo(Atom(ia),Atmo)

            ! If not protected
            else

              ! Revice populations in model atom
              call ReviseHatom(Atom(ia),Atmo)

            end if ! Protected
          end if ! Allocated populations

          ! Leave
          exit

        end if ! Is hydrogen

      end do ! Active atoms

      ! For each background atom
      do ia=1,nAb

        ! If it is hydrogen
        if (Atomb(ia)%element.eq.' H') then

          ! Save if protecting
          Atomb(ia)%mol_protect = Input%protect_H

          ! If allocated NLTE populations
          if (allocated(Atomb(ia)%popu)) then

            ! If protected
            if (Input%protect_H) then

              ! Revise populations in model atmosphere
              call ReviseHatmo(Atomb(ia),Atmo)

            ! If not protected
            else

              ! Revice populations in model atom
              call ReviseHatom(Atomb(ia),Atmo)

            end if ! Protected
          end if ! Allocated populations

          ! Leave
          exit

        end if ! Is hydrogen

      end do ! Background atoms

      end subroutine reviseH_out

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize LTE populations in all atoms\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data\n
      subroutine setlte(Atom,Atomb,Atmo,Input)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia


      ! Active atoms
      do ia=1,nA

        ! Get LTE populations
        call Initpopu_LTE(Atom(ia),Atmo,.True.,0)

      end do ! Active atoms

      ! Control
      if (laborted) return

      ! Background atoms
      do ia=1,nAb

        ! Get LTE populations
        call Initpopu_LTE(Atomb(ia),Atmo,.False.,0)

      end do ! Background atoms

      return

      ! Deveice compiler
      ia = Input%iter_min

      end subroutine setlte

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize LTE population of the ion of the LTE lines\n
      !!        Atomb(Atom_class(:)): Structures with atomic data for
      !!                              background atoms\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      subroutine setlte_lines(Atomb,LTElines,Atmo,Input)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia


      ! Copy line data from input
      LTElines = Input%LTEline

      ! For each LTE line
      do ia=1,nLTEl

        ! Memory counter
        MRAMc = MRAMc + 1d-6*sizeof(LTElines(ia))
        MRAMc = MRAMc + 1d-6*sizeof(LTElines(ia)%broad_args)

        ! Get LTE populations
        call Initpopu_LTE_line(LTElines(ia),Atmo)

        ! If there is a background model atom
        if (LTElines(ia)%is_passive) then

          ! Get abundance from background model atom
          LTElines(ia)%n = LTElines(ia)%n*Atomb(LTElines(ia)%ia)%abun

        ! Does not have a passive model atom
        else

          ! Get abundance from atmosphere
          LTElines(ia)%n = LTElines(ia)%n*Atmo%abund(LTElines(ia)%ele)

        end if ! Element has a background model atom

      end do ! LTE lines

      end subroutine setlte_lines

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the collisional rates in atoms and write the
      !! output\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols
      subroutine setcols(Atom,Atomb,Atmo,Input,Flgsg)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Atmo_class), intent(in):: Atmo
      type(Fctsg_class), intent(inout):: Flgsg
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia


      ! For each active atom
      do ia=1,nA

        ! Calculate collisional rates
        call Initcols(Atom(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.True.)

      end do ! Active atoms

      ! Control
      if (laborted) return

      ! If keeping collisions, call writer
      if (Input%keep_cols.and.nA.gt.0) &
        call writecols(Atom,Input%folder,&
                       Input%lim_cols_tt,Input%lim_cols_ll)

      ! Control
      if (laborted) return

      ! For each background atom
      do ia=1,nAb

        ! Calculate collisional rates
        call Initcols(Atomb(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.False.)

      end do ! Background atoms

      end subroutine setcols

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize populations and density matrices for all atoms\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data
      subroutine setuppopu(Atom,Atomb,Atmo)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Atmo_class), intent(in):: Atmo

      ! Local

      integer:: ia


      ! For each active atom
      do ia=1,nA

        ! Initialize populations
        call Initpopu(Atom(ia),Atmo,.True.)

      end do ! Active atoms

      ! Control
      if (laborted) return

      ! For each background atom
      do ia=1,nAb

        ! Initialize populations
        call Initpopu(Atomb(ia),Atmo,.False.)

      end do ! Background atoms

      end subroutine setuppopu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the line broadening damping coefficient\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data
      subroutine broadening(Atom,Atomb,Atmo,Input)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia


      ! For each active atom
      do ia=1,nA

        ! Compute line broadening damping parameter
        call broad(Atom(ia),Atmo,Input%folder,Input%keep_aparam)

      end do ! Active atoms

      ! Control
      if (laborted) return

      ! If keeping damping, call writer
      if (Input%keep_damp.and.nA.gt.0) &
        call writedamp(Atom,Atmo,Input%folder, &
                       Input%lim_damp)

      ! If keeping elastic rates, call writer
      if (Input%keep_qel.and.nA.gt.0) &
        call writeqel(Atom,Input%folder,Input%lim_qel)

      ! Control
      if (laborted) return

      ! If there are background atoms
      if (Input%nAb.gt.0) then

        ! For each background atom
        do ia=1,Input%nAb

          ! Compute line broadening damping parameter
          call broad(Atomb(ia),Atmo,Input%folder,.False.)

        end do ! Background atoms

      end if ! There are background atoms

      end subroutine broadening

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the line broadening damping coefficient for LTE
      !! lines\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric data
      subroutine broadening_line(LTElines,Atmo)

      ! I/O

      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTElines
      type(Atmo_class), intent(in):: Atmo

      ! Local

      integer:: ia


      ! For every LTE line
      do ia=1,nLTEl

        ! Compute line broadening damping parameter
        call broad_line(LTElines(ia),Atmo)

      end do ! LTE line

      end subroutine broadening_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the chemical equilibrium and revise model the hydrogen
      !! number density in the model atmosphere\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!        Atomb(Atom_class(:)): Structures with atomic data for
      !!                              background atoms\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!           Mol(Mol_class(:)): Structures with molecular data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data
      subroutine chemical(Atom,Atomb,LTEline,Mol,Atmo,Input)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTEline
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia


      !
      ! Protect atoms from chemical equilibrium
      !

      ! If protecting all atoms
      if (Input%chem_protect_all) then

        ! Protect every active atom
        do ia=1,nA
          Atom(ia)%mol_protect = .True.
        end do

        ! Protect every background atom
        do ia=1,nAb
          Atomb(ia)%mol_protect = .True.
        end do

      end if ! Protecting all atoms

      !
      ! Calculate chemical equilibrium
      call chemeq(Atom,Atomb,LTEline,Mol,Atmo,Input%protect_Hm)

      ! Control
      if (laborted) return

      !
      ! Check again atmospheric hydrogen number density
      !

      ! For each active atom
      do ia=1,nA

        ! If it is hydrogen
        if (Atom(ia)%element.eq.' H') then

          ! If it had NLTE populations, revise atmosphere
          if (allocated(Atom(ia)%popu)) &
            call ReviseHatmo(Atom(ia),Atmo)

          ! Leave
          exit

        end if ! If it hydrogen

      end do ! Active atoms

      ! For each background atom
      do ia=1,nAb

        ! If it is hydrogen
        if (Atomb(ia)%element.eq.' H') then

          ! If it had NLTE populations, revise atmosphere
          if (allocated(Atomb(ia)%popu)) &
            call ReviseHatmo(Atomb(ia),Atmo)

          ! Leave
          exit

        end if ! If it hydrogen

      end do ! Background atoms

      end subroutine chemical

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the density matrix of all atoms\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine Initrhoes(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ia


      ! For each active atom
      do ia=1,nA

        ! Initialize density matrix
        call Initcrho(Atom(ia))

      end do ! Active atoms

      end subroutine Initrhoes

!#####################################################################
!#####################################################################
!#####################################################################

      !> Update the model atmosphere, recalculating number densities
      !! if necessary\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data\n
      !!    Input(Input_class): Structure with configuration data
      subroutine updateatmo(Atom,Atomb,Atmo,Bfield,Input)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Input_class), intent(in):: Input

      ! Local

      integer:: ia
      integer, dimension(:), allocatable:: nlte,depar


      !
      ! Check if updating the model
      !

      if (Input%update_atmos.lt.0) return

      ! Activate all nlte and deactivate all depar flags
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

          ! Revise populations in model atmosphere and exit
          call ReviseHatmo(Atom(ia),Atmo)

        end if ! If Hydrogen

      end do ! Active atoms

      !
      ! Redo the electrons if specified in input
      if (Input%redo_ne.eq.1.or.Input%redo_ne.eq.11) &
        call redo_ne(Atom,Atomb,nlte,depar,Atmo)

      ! Write model atmosphere if indicated in input
      if (Input%keep_atmo) call writeatmo(Atmo,Bfield, &
                                          Input%folder, &
                                          Input%lim_atmo)

      ! Write new model in ASCII if 1D
      if (run_mode.eq.0) &
        call wAtmo(Atmo,Input%update_atmos,Input%folder, &
                   Input%atmo)

      return

      end subroutine updateatmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare a model atmosphere for the hanle routine to solve the
      !! self-consistent problem and/or calculate formal solutions\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!        Atomb(Atom_class(:)): Structures with atomic data for
      !!                              background atoms\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!           Mol(Mol_class(:)): Structures with molecular data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols
      subroutine prepare_syn(Atom,Atomb,LTElines,Mol,Atmo,Input,Flgsg)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTElines
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(Fctsg_class), intent(inout):: Flgsg
      type(Input_class), intent(in):: Input

      ! Local

      integer, dimension(:), allocatable:: nlte,depar


      ! Prepare space for the total population of atoms and molecules
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! If error, skip
      if (laborted) goto 1000

      ! Check populations from files
      call setpopufiles(Atom,Atomb,Atmo,Input,nlte,depar)

      ! If error, skip
      if (laborted) goto 1000

      ! Revise H population in model atmosphere
      call reviseH_init(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Eq. of state
      if (Atmo%typo.gt.0.or.Input%keep_atmo.or. &
          Input%redo_ne.gt.0) then

        ! Call equation of state
        call eqstate(Atmo,Atom,Atomb,nlte,depar)

        ! Control
        if (laborted) goto 1000

        ! Revise H population
        if (Atmo%typo.gt.0) &
          call reviseH_out(Atom,Atomb,Atmo,Input)

      end if ! Eq. of state

      ! If error, skip
      if (laborted) goto 1000

      ! Recalculate electron density
      if (Input%redo_ne.ge.10) &
        call redo_ne(Atom,Atomb,nlte,depar,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Calculate LTE populations
      call setlte(Atom,Atomb,Atmo,Input)

      ! If error, skip
      if (laborted) goto 1000

      ! Calculate collisions
      call setcols(Atom,Atomb,Atmo,Input,Flgsg)

      ! If error, skip
      if (laborted) goto 1000

      ! If there are LTE lines
      if (nLTEl.gt.0) then

        ! Prepare LTE lines
        call setlte_lines(Atomb,LTElines,Atmo,Input)

      end if

      ! Initialize populations
      call setuppopu(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Calculate line broadening damping parameter
      call broadening(Atom,Atomb,Atmo,Input)

      ! If error, skip
      if (laborted) goto 1000

      ! If there are LTE lines, calculate their line broadening
      ! damping parameter
      if (nLTEl.gt.0) call broadening_line(LTElines,Atmo)

      ! Solve chemical equilibrium
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

      !> Prepare a model atmosphere in an inversion to ensure that it
      !! is in the expected format\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!       zalt(double(:)): Alternative height axis\n
      !!        alloc(logical): If the optical depth axis is allocated
      subroutine setup_Atmo_ininv(Atom,Atomb,Mol,Atmo,Input,fudge, &
                                  zalt,alloc)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(fudge_class), intent(in):: fudge
      type(Input_class), intent(in):: Input
      logical, intent(in):: alloc
      double precision, dimension(:), allocatable, &
                        target, intent(inout):: zalt

      ! Local

      type(Atmo_class):: Atmo_tmp
      type(LTEline_class), dimension(:), allocatable:: dummy

      integer:: ii
      integer, dimension(Atmo%nele):: nlte,depar


      ! Keep current number of height nodes, nz
      ii = nz
      nz = Atmo%nz

      ! No data in inversion
      nlte = 0
      depar = 0

      ! Prepare total populations of atoms and molecules
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! If error, skip
      if (laborted) goto 1000

      ! Revise H populations
      call reviseH_init(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Eq. of state
      call eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! If error, skip
      if (laborted) goto 1000

      ! If model not in tau scale
      if (Atmo%scal.ne.'T') then

        ! Get a copy of Atmo
        call cAtmo(Atmo,Atmo_tmp)

        ! Initialize LTE populations
        call setlte(Atom,Atomb,Atmo_tmp,Input)

        ! Initialize populations
        call setuppopu(Atom,Atomb,Atmo_tmp)

        ! Solve chemical equilibrium
        call chemical(Atom,Atomb,dummy,Mol,Atmo_tmp,Input)

        ! If error, skip
        if (laborted) goto 1000

        ! Allocate chi500
        if (.not.allocated(Atmo_tmp%chi500)) then
          allocate(Atmo_tmp%chi500(nz))
          MRAMc = MRAMc + 1d-6*sizeof(Atmo_tmp%chi500)
        end if

        ! Calculate continuum opacity at reference frequency
        call chi_freq(Atom,Atomb,Mol,Atmo_tmp,fudge,Input, &
                      Atmo_tmp%tfreq,Atmo_tmp%chi500,1,nz, &
                      nproc.gt.1)

        ! Free memory
        call free_lpop(Atom,Atomb)
        call free_mol(Mol)

        !
        ! Compute missing height or tau
        !
        call getztau(Atmo_tmp,.False.,.False.)


        !
        ! Get the tau scale
        !

        ! If allocated
        if (alloc) then

          ! Copy axis from temporal model
          Atmo%z = Atmo_tmp%zalt

        ! Not allocated, but pointed
        else

          ! If zalt is not allocated, do it
          if (.not.allocated(zalt)) then
            allocate(zalt(nz))
            MRAMc = MRAMc + 1d-6*sizeof(zalt)
          end if

          ! Copy from model
          zalt = Atmo_tmp%zalt

          ! And point it
          Atmo%z => zalt

        end if ! Allocated or pointed

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

      ! Restore nz to the original value
      nz = ii

      ! If minimum tau is smaller or equal to zero
      if (minval(Atmo%z).le.0d0) then

        ! If really negative
        if (minval(Atmo%z).lt.0d0) then

          ! Issue error
          umsg = 'The input model atmosphere has '// &
                 'negative optical depths'
          urou = 'setup_Atmo'
          call aborted
          goto 1000

        end if ! Minimum optical depth is negative

        ! Try to follow some progression to compute the top
        ! optical depth value
        Atmo%z(1) = Atmo%z(2)*Atmo%z(2)/Atmo%z(3)

        ! If still negative, try with the latest step
        if (Atmo%z(1).le.0d0) &
          Atmo%z(1) = 2d0*Atmo%z(2) - Atmo%z(3)

        ! If still zero or negative, put the upper boundary at onei
        ! order of magnitude smaller than second to last
        if (Atmo%z(1).le.0d0) Atmo%z(1) = Atmo%z(2)*1d-1

        ! For consistency, add to all the depths
        Atmo%z(2:Atmo%nz) = Atmo%z(2:Atmo%nz) + Atmo%z(1)

        ! If it is still 0
        if (minval(Atmo%z).le.0d0) then

          ! Issue error and give up
          umsg = 'The input model atmosphere has '// &
                 'a non valid optical depth stratification'
          urou = 'setup_Atmo'
          call aborted

        end if ! Still 0 optical depth
      end if ! 0 or negative optical depth

      ! Free mass density if allocated
1000  if (allocated(Atmo%rho)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%rho)
        deallocate(Atmo%rho)
      end if

      ! Free electron pressure if allocated
      if (allocated(Atmo%Pe)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pe)
        deallocate(Atmo%Pe)
      end if

      ! Free local variables
      call free_local_Atom(Atom)
      call free_gpop(Atom,Atomb,Mol)

      end subroutine setup_Atmo_ininv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the atmospheric model resulting from a inversion to
      !! write in the result file\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!    fudge(fudge_class): Structure with fudge data
      subroutine setup_Atmo_ouinv(Atom,Atomb,Mol,Atmo,Input,fudge)


      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(fudge_class), intent(in):: fudge
      type(Input_class), intent(in):: Input

      ! Local

      type(LTEline_class), dimension(:), allocatable:: dummy

      integer, dimension(Atmo%nele):: nlte,depar


      ! No data in inversion
      nlte = 0
      depar = 0

      ! Prepare total populations of atoms and molecules
      call prepareatomol(Atom,Atomb,Mol,Input%nM)

      ! If error, skip
      if (laborted) goto 1000

      ! Revise H populations
      call reviseH_init(Atom,Atomb,Atmo)

      ! If error, skip
      if (laborted) goto 1000

      ! Eq. of state
      call eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! If error, skip
      if (laborted) goto 1000

      ! Initialize LTE populations
      call setlte(Atom,Atomb,Atmo,Input)

      ! Initialize populations
      call setuppopu(Atom,Atomb,Atmo)

      ! Solve chemical equilibrium
      call chemical(Atom,Atomb,dummy,Mol,Atmo,Input)

      ! If error, skip
      if (laborted) goto 1000

      ! Allocate chi500
      if (.not.allocated(Atmo%chi500)) then
        allocate(Atmo%chi500(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%chi500)
      end if

      ! Calculate continuum opacity at reference frequency
      call chi_freq(Atom,Atomb,Mol,Atmo,fudge,Input, &
                    Atmo%tfreq,Atmo%chi500,1,nz,nproc.gt.1)

      ! Free memory
1000  call free_lpop(Atom,Atomb)
      call free_mol(Mol)
      call free_local_Atom(Atom)
      call free_gpop(Atom,Atomb,Mol)

      end subroutine setup_Atmo_ouinv

!#####################################################################
!#####################################################################
!#####################################################################

      end module initmodel_mod
