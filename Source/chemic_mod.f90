      !> Chemical equilibrium
      module chemic_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/19/2017
!  Last version:
!     08/07/2023 V3.0.7
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.7 - The chemical equilibrium must reduce the
!                             populations of LTE lines as well, if
!                             involved in molecules (TdPA)
!
!     07/11/2023:    V3.0.6 - Some combinations of temperature and
!                             number densities can lead to the Newton
!                             method in chemical equilibrium to
!                             produce negative populations in atoms
!                             and molecules that cannot be
!                             self-corrected, leading to exhausting
!                             the iterations without hope to converge.
!                             I have added a "shake-up" snippet of
!                             code to tackle that. At certain
!                             iteration steps, the code checks for
!                             negative populations. It changes the
!                             populations of molecules with negative
!                             populations and reduces the populations
!                             of the rest of molecules until there are
!                             no negative atomic populations. This
!                             cannot always solve the issue, so the
!                             safest bet is to avoid temperatures
!                             below 2000 K (TdPA)
!
!     07/03/2023:    V3.0.5 - The molecular quantities are initialized
!                             when computing the chemical equilibrium,
!                             as they are not used elsewhere (TdPA)
!                           - The electron density is no longer
!                             recycled, and the equation of state
!                             expects the actual input variable to
!                             be defined (TdPA)
!
!     03/08/2023:    V3.0.4 - In chemeq, do not signal exhausting the
!                             iterations if doing inversions, as it
!                             can lead to significant spam of warning
!                             lines as it is called more times than in
!                             synthesis (TdPA)
!
!     02/14/2023:    V3.0.3 - Limited chemeq messages to the 1D
!                             synthesis case. Suggested by Hao (TdPA)
!
!     07/13/2022:    V3.0.2 - The resource argument in no longer
!                             needed in chemeq, eqstate, and
!                             redo_ne (TdPA)
!                           - Atmo is now a required argument for
!                             initializenlte, because it contains
!                             the number of elements (TdPA)
!                           - The resource argument in no longer
!                             neeeded when calling getpf (TdPA)
!                           - Because getpf does not read files
!                             anymore, it is not necessary to check
!                             for failure after calling it (TdPA)
!                           - Abundances are now taken from the Atmo
!                             structure and not from the hard-coded
!                             table in the chemicaux_mod module (TdPA)
!                           - Number of elements is now taken from
!                             the Atmo structure and not from the
!                             chemicaux_mod module (TdPA)
!
!     07/08/2022:    V3.0.1 - The message about modified populations
!                             in the chemical equilibrium is only
!                             written if the global master is
!                             participating, that is, in the single
!                             1D synthesis (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case return
!                             clauses have been added after every
!                             call to aborted or control, as well
!                             as checks to determine is the
!                             return is necessary (TdPA)
!                           - Fixed a typo in an error message (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     11/12/2020:    V1.3.2 - The chemical equilibrium now issues a
!                             warning if the maximum number of
!                             iterations is reached and if
!                             negative populations are obtained (TdPA)
!
!     09/11/2020:    V1.3.1 - eqstate and redo_ne now take depar as
!                             argument and pass it along to other
!                             subroutines (TdPA)
!                           - redo_ne takes into account departure
!                             coefficients to compute electron
!                             contributions if provided (TdPA)
!                           - initializenlte and activenlte take
!                             into account the possibility to input
!                             departure coefficients instead of atomic
!                             populations (TdPA)
!
!     03/05/2020:    V1.3.0 - Chemeq now takes into account the
!                             different ways of working with Hydrogen,
!                             as well as the option to protect it
!                             from reducing its population by
!                             molecules (TdPA)
!                           - Limited the repetition of the message
!                             about atoms without neutral stage in
!                             chemeq (TdPA)
!                           - Added possibility to protect atomic
!                             populations from changing properly, that
!                             is, consistently with the equations to
!                             be solved in chemeq (TdPA)
!                           - Removed correction for charged molecules
!                             in chemeq because I am not sure it
!                             should be there (TdPA)
!                           - Bugfix: In chemeq, one should check the
!                             absolute change, not the change itself,
!                             in order to decide if converged (TdPA)
!                           - Severely changed eqstate, now it takes
!                             into account your inputs for nlte
!                             populations (TdPA)
!                           - Added redo_ne routine (TdPA)
!                           - Added initializenlte routine (TdPA)
!                           - Added activatenlte routine (TdPA)
!
!     01/14/2020:    V1.2.2 - In eqstate, forgot that recallabund
!                             already gives abundances in fractions
!                             and not in 12 + log10 (TdPA)
!
!     11/19/2019:    V1.2.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     09/26/2019:    V1.2.0 - Added eqstate subroutine, that computes
!                             the hydrogen and electron number
!                             density from electron number density,
!                             electron pressure, electron density,
!                             gas pressure, or gas density (TdPA)
!
!     04/08/2019:    V1.1.1 - Change for compatibility with new
!                             getpf arguments (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     12/18/2017:    V1.0.4 - Added explicit dependency of
!                             parameters_mod (TdPA)
!
!     09/15/2017:    V1.0.3 - Receiving Input%resource (TdPA)
!
!     07/20/2017:    V1.0.2 - Moved some parameters to simple
!                             variables (TdPA)
!
!     05/01/2017:    V1.0.1 - Limiting calculations to the height
!                             range of each processor (TdPA)
!
!     04/19/2017:    V1.0.0 - First version (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!#####################################################################
!#####################################################################
!
!  Notes:
!
!    The method for chemeq and for redo_ne were originally the same
!  than in the RH code, except for the added functionalities.
!
!    The method for eqstate was originally the same than in the SIR
!  code, except for the added functionalities.
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  chemeq
!    This subroutine calculates chemical equilibrium with molecules
!
!  eqstate
!    Computes number densities of electrons and hydrogen from
!  pressure/densities (gas or electron)
!
!  redo_ne
!    Recalculate the electron density from temperature and other
!  atom densities (LTE if not in atom list)
!
!  initializenlte
!    Allocates and initializes the nlte array, with information of
!  what atomic populations are available
!
!  actiavenlte
!    Set all active atoms as available in terms of populations
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use parameters_mod , only : hplanck , PI , me , kb, Avog
      use initpopuaux_mod
      use rmol_mod
      use types_mod

      ! Parameters

      ! Maximum iterations when input is density in eqstate
      integer, parameter:: maxiter = 3
      integer, parameter:: maxneiter = 15
      double precision, parameter:: epsne = 1d-2

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the chemical equilibrium to get the molecular
      !! populations, and modifies the atomic populations if
      !! necessary\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!       Atomb(Atom_class): Structure with the atomic data for
      !!                          background opacities\n
      !! LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Mol(Mol_class): Structure with the molecule data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data
      subroutine chemeq(Atom,Atomb,LTElines,Mol,Atmo)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), intent(inout):: Atomb
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Mol_class), dimension(:), intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo

      ! Local

      type(catm_class), dimension(nM*2):: atoms
      type(Atom_class):: Atom_l

      character(len=2):: sname

      logical:: nfound,warning,nwarning

      integer:: iz, ia, iterm, iJ, ilevel, iter, info, nres
      integer:: iatom, imol, ieq, iatom1, natom, nmol, neq
      integer:: max_it_chem, check_it_chem, check_it_res, nres_opt
      integer, dimension(:), allocatable:: indx
      integer, dimension(:,:), allocatable:: atom_index

      double precision:: phiHm, aaHm, C0, phil, frac
      double precision:: mrc_chem, minv, maxv
      double precision, dimension(:), allocatable:: aa, xx, bb, phi
      double precision, dimension(:,:), allocatable:: aap, frc


      ! Prepare molecular data
      call setupmol_eq(Mol,nM,Atmo)

      ! Parameters
      max_it_chem = 500
      check_it_chem = 20
      check_it_res = 3
      nres_opt = 4
      mrc_chem = 1d-5

      ! Initialize
      warning = .True.
      nwarning = .True.

      ! Constant needed below
      C0 = hplanck*hplanck/2d0/PI/me/kb

      ! If there are molecules, solve the chemical equilibrium
      if (nM.gt.0) then

        ! Local variable
        nmol = nM

        ! First atom in the sistem of equations is always hydrogen
        natom = 1
        atoms(1)%s = ' H'
        atoms(1)%inmod = .True.
        atoms(1)%abun = 1d0

        ! Initialize the atoms structure
        allocate(atoms(1)%imol(nM))
        atoms(1)%imol = 0
        allocate(atoms(1)%nmol(nM))
        atoms(1)%nmol = 0

        ! Get the partition function for the H atom
        call getpf(' H',atoms(1)%nstg,atoms(1)%pf, &
                   atoms(1)%Eion,Atmo)


        !
        ! Go through all the molecules and add the elements that
        ! constitute them to the list
        !

        ! For each molecule
        do imol=1,nmol

          ! For each component of the molecule
          do iatom=1,Mol(imol)%nA

            ! Reset the not found flag
            nfound = .True.

            ! For each atom already listed
            do iatom1=1,natom

              ! If the atom was in the list, at the molecular
              ! information
              if (Mol(imol)%catom(iatom)%s.eq.atoms(iatom1)%s) then

                atoms(iatom1)%pnmol = atoms(iatom1)%pnmol + 1
                atoms(iatom1)%imol(atoms(iatom1)%pnmol) = imol
                atoms(iatom1)%nmol(atoms(iatom1)%pnmol) = &
                                                Mol(imol)%natom(iatom)
                nfound = .False.

                ! And tell the molecule which atom it was
                Mol(imol)%iatom(iatom) = iatom1
                exit

              end if

            end do ! atoms in the list

            ! If not in the list, then add it
            if (nfound) then

              ! Run the index and fill its information
              natom = natom + 1
              allocate(atoms(natom)%imol(nM))
              atoms(natom)%imol(nM) = 0
              allocate(atoms(natom)%nmol(nM))
              atoms(natom)%nmol(nM) = 0
              atoms(natom)%s = Mol(imol)%catom(iatom)%s
              atoms(natom)%pnmol = 1
              atoms(natom)%imol(1) = imol
              atoms(natom)%nmol(1) = Mol(imol)%natom(iatom)

              ! Get the partition function of this atom
              call getpf(atoms(natom)%s,atoms(natom)%nstg, &
                         atoms(natom)%pf,atoms(natom)%Eion,Atmo)

              ! Look for the abundance in the list of active atoms
              nfound = .True.
              do ia=1,nA

                if (atoms(natom)%s.eq.Atom(ia)%Element) then
                  nfound = .False.
                  atoms(natom)%inmod = .True.
                  atoms(natom)%abun = Atom(ia)%abun
                  exit
                end if

              end do

              ! If not found, look in the background atoms
              if (nfound) then
                do ia=1,nAb

                  if (atoms(natom)%s.eq.Atomb(ia)%Element) then
                    nfound = .False.
                    atoms(natom)%inmod = .True.
                    atoms(natom)%abun = Atomb(ia)%abun
                    exit
                  end if

                end do
              end if

              ! If not found, look in the table
              if (nfound) then

                atoms(natom)%inmod = .False.
                atoms(natom)%abun = &
                    Atmo%abund(atom_char2index(atoms(natom)%s))

              end if

              ! And tell the molecule which atom it was
              Mol(imol)%iatom(iatom) = natom

            end if

          end do ! Components of the molecule
        end do ! Molecule

        !
        ! Construct the system of equations
        !

        ! Number of equations is number of atoms in molecules, plus
        ! Hydrogen (if not in molecule) and Molecules
        neq = natom + nM

        ! Allocate variables
        ! Distance to total dissociation
        allocate(aa(neq))
        ! Solution
        allocate(xx(neq))
        ! Independent term
        allocate(bb(neq))
        ! Molecular coefficient
        allocate(phi(neq))
        ! System of equations (SoE)
        allocate(aap(neq,neq))
        ! Variable needed by BLAS
        allocate(indx(neq))
        ! Fraction of atoms in the first two stages
        allocate(frc(2,neq))
        ! Indexing of the atoms in the system
        allocate(atom_index(natom,2))
        atom_index = -1

        !
        ! Index the variables in the system to retrieve later
        !

        ! For each atom involved in the chemical equilibrium
        do iatom=1,natom

          ! Look for it between the active atoms
          do ia=1,nA

            if (Atom(ia)%Element.eq.atoms(iatom)%s) then
              atom_index(iatom,1) = ia
              atom_index(iatom,2) = 1
              exit
            end if

          end do

          ! Look for it between the background atoms
          do ia=1,nAb

            if (Atomb(ia)%Element.eq.atoms(iatom)%s) then
              atom_index(iatom,1) = ia
              atom_index(iatom,2) = 2
              exit
            end if

          end do

        end do ! Atoms in the SoE


        !
        ! Solve the chemical equilibrium
        !

        ! For each height
        do iz=1,nz

          ! Re_initialize
          frc = 0d0
          bb = 0d0
          nres = 0

          ! For each atom in the SoE
          do iatom=1,natom

            ! Independent element
            ! If Hydrogen
            if (iatom.eq.1) then

              ! If there is a model
              if (atoms(iatom)%inmod) then

                ! Get model
                if (atom_index(iatom,2).eq.1) then
                  Atom_l = Atom(atom_index(iatom,1))
                else if (atom_index(iatom,2).eq.2) then
                  Atom_l = Atomb(atom_index(iatom,1))
                else
                  umsg = 'Something went wrong when '// &
                         'indexing the atoms in chemical '// &
                         'equilibrium'
                  urou = 'chemeq'
                  call aborted
                  return
                end if

                ! If protected
                if (Atom_l%mol_protect) then
                  bb(iatom) = Atom_l%n(iz)
                ! If not protected
                else
                  bb(iatom) = Atmo%nHt(iz)
                end if

              ! If no model
              else
                bb(iatom) = Atmo%nHt(iz)
              end if

            ! Any other element
            else
              bb(iatom) = atoms(iatom)%abun*Atmo%nHT(iz)
            end if


            !
            ! Get fraction in neutral and first ion
            !

            ! Initialize the not found flag
            nfound = .True.

            ! If the atom has a model, store it in the local variable
            if (atoms(iatom)%inmod) then

              if (atom_index(iatom,2).eq.1) then
                Atom_l = Atom(atom_index(iatom,1))
              else if (atom_index(iatom,2).eq.2) then
                Atom_l = Atomb(atom_index(iatom,1))
              else
                umsg = 'Something went wrong when '// &
                       'indexing the atoms in chemical '// &
                       'equilibrium'
                urou = 'chemeq'
                call aborted
                return
              end if

              ! Check that the model have neutral stage
              if (Atom_l%stage(1).eq.1) then
                nfound = .False.
              else
                if (pid.eq.0.and.iz.eq.1) then
                  umsg = ' - The atomic model for atom '// &
                         atoms(iatom)%s//' does not have'// &
                         ' neutral stage, using '// &
                         'partition functions for the '// &
                         'chemical equilibrium'
                  call verbose
                end if ! Master in first height
              end if ! Model has neutral
            end if ! There is model

            ! If we don't have the necessary populations, use
            ! partition functions to calculate the fractions
            if (nfound) then

              call getfrc(atoms(iatom)%nstg,atoms(iatom)%pf(:,iz), &
                          atoms(iatom)%Eion,Atmo%T(iz),Atmo%ne(iz), &
                          -1,frc(:,iatom))

            ! If we have the populations, use them
            else

              ilevel = 0
              do iterm=1,Atom_l%nMulti
                if (Atom_l%stage(iterm).gt.2) exit
                do iJ=1,Atom_l%nJ(iterm)

                  ilevel = ilevel + 1

                  frc(Atom_l%stage(iterm),iatom) = &
                                  frc(Atom_l%stage(iterm),iatom) + &
                                  Atom_l%popu(ilevel,iz)

                end do
              end do

              frc(:,iatom) = frc(:,iatom)/Atom_l%n(iz)

            end if

          end do ! atoms in molecules

          ! Add Hminus
          phiHm = .25d6*((C0/Atmo%T(iz))**(1.5d0))* &
                  exp(8.74980963338d3/Atmo%T(iz))
          aaHm = Atmo%ne(iz)*frc(1,1)*phiHm

          ! Get equilibrium constant for every molecule
          do imol=1,nM
            phi(imol) = Mol(imol)%eq(iz)
          end do

          ! Initial solution, complete dissociation
          xx = bb


          !
          ! Force protections
          !

          ! Atoms in SoE
          do iatom=1,natom

            ! If has model
            if (atoms(iatom)%inmod) then
              if (atom_index(iatom,2).eq.1) then
                Atom_l = Atom(atom_index(iatom,1))
              else if (atom_index(iatom,2).eq.2) then
                Atom_l = Atomb(atom_index(iatom,1))
              else
                umsg = 'Something went wrong when '// &
                       'indexing the atoms in chemical '// &
                       'equilibrium'
                urou = 'chemeq'
                call aborted
                return
              end if
            end if ! Has model

            ! Check if protected
            if (Atom_l%mol_protect) xx(iatom) = Atom_l%n(iz)

          end do ! Atoms in SoE


          !
          ! Iterate the SoE for this height
          !
          do iter=1,max_it_chem

            ! Reset the SoE matrix
            aap = 0d0

            ! Compute the difference between solution and independent
            ! term of total dissociation
            aa = xx - bb

            ! Diagonal is identity
            do ieq=1,neq

              aap(ieq,ieq) = 1d0

            end do

            ! Add Hminus to the Hydrogen row
            aa(1) = aa(1) + aaHm*xx(1)
            aap(1,1) = aap(1,1) + phiHm

            ! Run over molecules
            do imol=1,nM

              ! Determine the row for this molecule and get the
              ! equilibrum constant
              ieq = imol + natom
              phil = phi(imol)

              ! Run over components
              do ia=1,Mol(imol)%nA

                ! Ask which atom it is in the list
                iatom = Mol(imol)%iatom(ia)

                ! Update the coefficient using the fraction of
                ! populations
                phil = phil*((frc(1,iatom)*xx(iatom))** &
                             (Mol(imol)%natom(ia)))

               !! Correction for charged molecules
               !if (Mol(imol)%Charge.gt.0) &
               !  phil = phil*frc(2,iatom)/frc(1,iatom)

                ! Add the contribution to the independent term of the
                ! row of the atomic component
                aa(iatom) = aa(iatom) + &
                            Mol(imol)%natom(ia)*xx(ieq)

              end do ! Atomic components of molecule

              ! Get the new independent term
              phil = phil/((Atmo%ne(iz))**Mol(imol)%Charge)
              aa(ieq) = aa(ieq) - phil

              ! Derivative
              do ia=1,Mol(imol)%nA

                ! Ask which atom it is in the list
                iatom = Mol(imol)%iatom(ia)

                ! Contribution due to molecule in the atomic row
                aap(iatom,ieq) = aap(iatom,ieq) + &
                                 Mol(imol)%natom(ia)
                ! Contribution due to atom in the molecule row
                aap(ieq,iatom) = -phil*Mol(imol)%natom(ia)/xx(iatom)

              end do

            end do ! Molecules

            !
            ! Force protections
            !

            ! Atoms in SoE
            do iatom=1,natom

              ! If has model
              if (atoms(iatom)%inmod) then
                if (atom_index(iatom,2).eq.1) then
                  Atom_l = Atom(atom_index(iatom,1))
                else if (atom_index(iatom,2).eq.2) then
                  Atom_l = Atomb(atom_index(iatom,1))
                else
                  umsg = 'Something went wrong when '// &
                         'indexing the atoms in chemical '// &
                         'equilibrium'
                  urou = 'chemeq'
                  call aborted
                  return
                end if
              end if ! Has model

              ! Check if protected
              if (Atom_l%mol_protect) then
                aap(iatom,:) = 0d0
                aap(iatom,iatom) = 1d0
                aa(iatom) = 0d0
              end if

            end do ! Atoms in SoE

            ! Solve system of equations
            call DGESV(neq,1,aap,neq,indx,aa,neq,info)

            ! Get the new solution
            xx = xx - aa

            ! And calculate the MRC
            do ieq=1,neq

              if (abs(xx(ieq)).gt.1d-30) aa(ieq) = aa(ieq)/xx(ieq)

            end do

            ! If we reached convergence, go out
            if (maxval(abs(aa)).lt.mrc_chem.and.minval(xx).ge.0d0) &
              exit

            ! Warning iterations
            if (iter.eq.max_it_chem.and.nwarning.and.pid.eq.0) then
              if (minval(xx).lt.0d0) then
                umsg = 'Negative fractions obtained in '// &
                       'chemical equilibrium'
                urou = 'chemical'
                call abortedS(umsg,urou,-1,.True.,.True.)
                nwarning = .False.
                warning = .False.
                cycle
              else if (warning.and.run_mode.ne.-1) then
                umsg = 'Maximum number of iterations reached '//&
                       'in the chemical equilibrium'
                urou = 'chemical'
                call abortedS(umsg,urou,-1,.False.,.True.)
                warning = .False.
                cycle
              end if ! Negative solution

            ! Last iteration
            else if (iter.eq.max_it_chem) then

              cycle

            end if ! Can issue more warnings and last iteration


            ! If we have iterated for a while times and populations
            ! are out of control
            if (mod(iter,check_it_chem).eq.0.and. &
                minval(xx).lt.0d0) then

              ! Reset atomic populations to dissociation
              xx(1:natom) = bb(1:natom)

              ! If we are beyond the number of resets
              if (nres.gt.check_it_res) then

                ! Try with minimum
                if (mod(nres,nres_opt).eq.0) then

                  !
                  ! Get minimum non-zero
                  minv = minval(xx(natom+1:neq), &
                                xx(natom+1:neq).gt.0d0)

                  ! For each molecule
                  do imol=1,nM

                    ! Get equation index
                    ieq = imol + natom

                    ! Remove negatives
                    if (xx(ieq).lt.0d0) xx(ieq) = minv

                    ! Run over components
                    do ia=1,Mol(imol)%nA

                      ! Ask which atom it is in the list
                      iatom = Mol(imol)%iatom(ia)

                      ! Remove molecules
                      xx(iatom) = xx(iatom) - &
                                  Mol(imol)%natom(ia)*xx(ieq)


                    end do ! Components
                  end do ! Molecules

                ! Try with maximum
                else if (mod(nres,nres_opt).eq.1) then

                  !
                  ! Get maximum
                  maxv = maxval(xx(natom+1:neq))

                  ! For each molecule
                  do imol=1,nM

                    ! Get equation index
                    ieq = imol + natom

                    ! Remove negatives
                    if (xx(ieq).lt.0d0) xx(ieq) = maxv

                    ! Run over components
                    do ia=1,Mol(imol)%nA

                      ! Ask which atom it is in the list
                      iatom = Mol(imol)%iatom(ia)

                      ! Remove molecules
                      xx(iatom) = xx(iatom) - &
                                  Mol(imol)%natom(ia)*xx(ieq)


                    end do ! Components
                  end do ! Molecules

                ! Try average
                else if (mod(nres,nres_opt).eq.2) then

                  !
                  ! Get minimum non-zero
                  minv = minval(xx(natom+1:neq),xx(natom+1:neq).gt.0d0)
                  ! Get maximum non-zero
                  maxv = maxval(xx(natom+1:neq))

                  ! For each molecule
                  do imol=1,nM

                    ! Get equation index
                    ieq = imol + natom

                    ! Remove negatives
                    if (xx(ieq).lt.0d0) xx(ieq) = 0.5d0*(minv+maxv)

                    ! Run over components
                    do ia=1,Mol(imol)%nA

                      ! Ask which atom it is in the list
                      iatom = Mol(imol)%iatom(ia)

                      ! Remove molecules
                      xx(iatom) = xx(iatom) - &
                                  Mol(imol)%natom(ia)*xx(ieq)


                    end do ! Components
                  end do ! Molecules

                ! Try abs
                else if (mod(nres,nres_opt).eq.3) then

                  ! For each molecule
                  do imol=1,nM

                    ! Get equation index
                    ieq = imol + natom

                    ! Remove negatives
                    if (xx(ieq).lt.0d0) xx(ieq) = abs(xx(ieq))

                    ! Run over components
                    do ia=1,Mol(imol)%nA

                      ! Ask which atom it is in the list
                      iatom = Mol(imol)%iatom(ia)

                      ! Remove molecules
                      xx(iatom) = xx(iatom) - &
                                  Mol(imol)%natom(ia)*xx(ieq)


                    end do ! Components
                  end do ! Molecules

                end if ! Type of try

              ! First resets
              else

                ! For each molecule
                do imol=1,nM

                  ! Get equation index
                  ieq = imol + natom

                  ! Remove negatives
                  if (xx(ieq).lt.0d0) xx(ieq) = 0d0

                  ! Run over components
                  do ia=1,Mol(imol)%nA

                    ! Ask which atom it is in the list
                    iatom = Mol(imol)%iatom(ia)

                    ! Remove molecules
                    xx(iatom) = xx(iatom) - &
                                Mol(imol)%natom(ia)*xx(ieq)


                  end do ! Components
                end do ! Molecules

              end if

              ! While there are negative atomic populations
              do while (minval(xx(1:natom)).lt.0d0)

                ! Reduce molecules by 10
                xx(natom+1:neq) = xx(natom+1:neq)*0.1d0

                ! Reset to dissociation
                xx(1:natom) = bb(1:natom)

                ! For each molecule
                do imol=1,nM

                  ! Get equation index
                  ieq = imol + natom

                  ! Run over components
                  do ia=1,Mol(imol)%nA

                    ! Ask which atom it is in the list
                    iatom = Mol(imol)%iatom(ia)

                    ! Remove molecules
                    xx(iatom) = xx(iatom) - &
                                Mol(imol)%natom(ia)*xx(ieq)


                  end do ! Components
                end do ! Molecules
              end do ! While there are negative populations

              ! Add reset to counter
              nres = nres + 1

            end if ! Out of control populations

          end do ! Iterations


          !
          ! Use the solution to affect the populations of atoms
          ! and molecules
          !

          ! For each atom, if it is loaded, modify its population
          do iatom=1,natom

            ! Check if it has model
            if (atom_index(iatom,2).gt.0) then

              ! Get its index
              ia = atom_index(iatom,1)

              ! If it is active
              if (atom_index(iatom,2).eq.1) then

                ! If protected, skip
                if (Atom(ia)%mol_protect) then
                  if (iz.eq.1.and.pid.eq.0.and.run_mode.eq.0) then
                    umsg = ' - '//atoms(iatom)%s// &
                           ' is in molecules, but its '// &
                           ' read population was kept'
                    call verbose
                  end if
                  cycle
                end if

                frac = xx(iatom)/Atom(ia)%n(iz)
                Atom(ia)%popu(:,iz) = Atom(ia)%popu(:,iz)*frac
                Atom(ia)%populte(:,iz) = Atom(ia)%populte(:,iz)*frac
                Atom(ia)%n(iz) = xx(iatom)

              ! If it is background
              else if (atom_index(iatom,2).eq.2) then

                ! If protected, skip
                if (Atomb(ia)%mol_protect) then
                  if (iz.eq.1.and.pid.eq.0.and.run_mode.eq.0) then
                    umsg = ' - Atom '//atoms(iatom)%s// &
                           ' is in molecules, but the '// &
                           'read populations were kept'
                    call verbose
                  end if
                  cycle
                end if

                frac = xx(iatom)/Atomb(ia)%n(iz)
                Atomb(ia)%popu(:,iz) = Atomb(ia)%popu(:,iz)*frac
                Atomb(ia)%populte(:,iz) = Atomb(ia)%populte(:,iz)*frac
                Atomb(ia)%n(iz) = xx(iatom)

              end if

              ! Notify the change
              if (iz.eq.1.and.gpid.eq.0.and.run_mode.eq.0) then
                umsg = ' - Population of '//atoms(iatom)%s// &
                       ' modified by molecules'
                call verbose
              end if

            end if

          end do

          ! If LTE lines and allocated
          if (nLTEl.gt.0.and.allocated(LTElines)) then

            ! For each LTE line
            do ia=1,nLTEl

              sname = atom_index2char(LTElines(ia)%ele)

              ! For each atom in chemical eq.
              do iatom=1,natom

                ! If same atom
                if (atoms(iatom)%s.eq.sname) then

                  ! Get population
                  LTElines(ia)%n(iz) = xx(iatom)

                end if

              end do ! Atoms in chemical eq.
            end do ! LTE lines

          end if ! LTE lines


          !
          ! Store Hm density
          !
          Atmo%nHm(iz) = Atmo%ne(iz)*xx(1)*phiHm


          !
          ! Store molecular density
          !
          do imol=1,nM
            Mol(imol)%n(iz) = xx(imol + natom)
          end do

        end do ! Heights

      ! If no molecules, calculate just Hminus
      else

        do iz=1,nz

          phiHm = .25d6*((C0/Atmo%T(iz))**(1.5d0))* &
                  exp(8.74980963338d3/Atmo%T(iz))
          Atmo%nHm(iz) = Atmo%ne(iz)*Atmo%nHT(iz)*phiHm

        end do

      end if

      ! Control
      call control

      return

      end subroutine chemeq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the equation of state with the Wittmann method
      !! as in the SIR code, to go from densities/pressures to number
      !! densities of H and electrons\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!       Atomb(Atom_class): Structure with the atomic data for
      !!                          background opacities\n
      !!        nlte(integer(:)): Array with information about
      !!                          loaded populations\n
      !!       depar(integer(:)): Array with information about
      !!                          loaded departure coefficients
      subroutine eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom, Atomb
      type(Atmo_class), intent(inout):: Atmo
      integer, dimension(:), intent(in):: nlte, depar

      ! Local

      type(catm_class), dimension(:), allocatable:: atoms

      integer:: iz,ia,natom,iter,istg

      double precision:: kbcgs,ikbcgs,Watom,Aatom,daux

      ! Inverse of Boltzmann in cgs
      kbcgs = kb*1d7
      ikbcgs = 1d-7/kb

      ! Compute sum of atomic weights and abundances
      Watom = 0d0
      Aatom = 0d0
      natom = Atmo%nele

      ! Allocate partition function
      allocate(atoms(natom))

      ! Get partition functions
      do ia=1,natom
        atoms(ia)%s = atom_index2char(ia)
        atoms(ia)%abun = Atmo%abund(ia)
        call getpf(atoms(ia)%s, atoms(ia)%nstg, atoms(ia)%pf, &
                   atoms(ia)%Eion, Atmo)

        ! For each height and stage, convert to linear partition
        do iz=1,nz
          do istg=1,atoms(ia)%nstg
            atoms(ia)%pf(istg,iz) = exp(atoms(ia)%pf(istg,iz))
          end do
        end do

        ! Convert energy to eV
        atoms(ia)%Eion = atoms(ia)%Eion*fktoev

      end do

      ! Allocate

      ! Electron Pressure
      if (.not.allocated(Atmo%Pe)) then
        allocate(Atmo%Pe(nZ))
        Atmo%Pe = 0d0
      end if

      ! Gas Pressure
      if (.not.allocated(Atmo%Pg)) then
        allocate(Atmo%Pg(nZ))
        Atmo%Pg = 0d0
      end if

      ! Density
      if (.not.allocated(Atmo%rho)) then
        allocate(Atmo%rho(nZ))
        Atmo%rho = 0d0
      end if

      ! For each contribution
      do ia=1,natom
        daux = Atmo%abund(ia)
        Watom = Watom + recallmass_ind(ia)*daux
        Aatom = Aatom + daux
      end do

      ! If only here because you want to output an atmospheric file
      if (Atmo%typo.eq.0) then

        ! Compute pressure
        Atmo%Pe = Atmo%ne*kbcgs*Atmo%T

        ! Compute partial pressures
        call eqstate_known(Atmo,Atom,Atomb,nlte,depar,atoms)

        ! Preliminar rho
        Atmo%rho = Atmo%Pg/Atmo%T
        Atmo%rho = Atmo%rho*ikbcgs/Avog

        ! For each height, get density
        do iz=1,nz
          daux = Watom/(Aatom + Atmo%Pe(iz)/Atmo%Pg(iz))
          Atmo%rho(iz) = Atmo%rho(iz)*daux
        end do

      ! If an electronic quantity if given
      else if (Atmo%typo.le.3) then

        ! If number density
        if (Atmo%typo.eq.1) then

          ! Compute pressure
          Atmo%Pe = Atmo%ne*kbcgs*Atmo%T

        ! If no number density
        else

          ! Error
          urou = 'eqstate'
          umsg = 'Eq. of state called with Pe or rhoe, not '// &
                 'supposed to happen'
          call abortedS(umsg,urou,-1,.True.,.True.)
          call control
          return

        end if ! Type of electron input

        call eqstate_ele(Atmo,Atom,Atomb,nlte,depar,atoms)

        ! Preliminar rho
        Atmo%rho = Atmo%Pg/Atmo%T
        Atmo%rho = Atmo%rho*ikbcgs/Avog

        ! For each height, get density
        do iz=1,nz
          daux = Watom/(Aatom + Atmo%Pe(iz)/Atmo%Pg(iz))
          Atmo%rho(iz) = Atmo%rho(iz)*daux
        end do

      ! If gas pressure
      else if (Atmo%typo.eq.4) then

        call eqstate_gas(Atmo,Atom,Atomb,nlte,depar,atoms)

        ! Preliminar rho
        Atmo%rho = Atmo%Pg/Atmo%T
        Atmo%rho = Atmo%rho*ikbcgs/Avog

        ! For each height, compute density
        do iz=1,nz
          daux = Watom/(Aatom + Atmo%Pe(iz)/Atmo%Pg(iz))
          Atmo%rho(iz) = Atmo%rho(iz)*daux
        end do

      ! If density
      else if (Atmo%typo.eq.5) then

        ! Stimate density
        daux = Avog*kbcgs*Aatom/Watom
        Atmo%Pg = Atmo%rho*Atmo%T*daux

        ! Iterate
        do iter=1,maxiter

          ! Compute Pe
          call eqstate_gas(Atmo,Atom,Atomb,nlte,depar,atoms)

          ! Preliminar Pg
          Atmo%Pg = Atmo%rho*Atmo%T
          Atmo%Pg = Atmo%Pg*kbcgs*Avog

          ! For each height, compute density
          do iz=1,nz
            daux = Watom/(Aatom + Atmo%Pe(iz)/Atmo%Pg(iz))
            Atmo%Pg(iz) = Atmo%Pg(iz)/daux
          end do ! Heights
        end do ! Iterations

      end if ! Type of input

      ! Deallocate partition function
      deallocate(atoms)

      call control

      return

      end subroutine eqstate

!#####################################################################
!#####################################################################
!#####################################################################

      !> Recalculate the electron density from temperature and
      !! other atom densities (LTE if not specified in input)\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!       Atomb(Atom_class): Structure with the atomic data for
      !!                          background opacities\n
      !!        nlte(integer(:)): Array with information about
      !!                          loaded populations\n
      !!       depar(integer(:)): Array with information about
      !!                          loaded departure coefficients\n
      !!        Atmo(Atmo_class): Structure with atmospheric data
      subroutine redo_ne(Atom,Atomb,nlte,depar,Atmo)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom, Atomb
      type(Atmo_class), intent(inout):: Atmo
      integer, dimension(:), intent(in):: nlte, depar

      ! Local

      type(catm_class), dimension(:), allocatable:: atoms

      logical:: active
      integer:: natom,iatom,ia,istg,iterm,iJ,ilevel,iter,iz,Mstg

      double precision:: C0,ikb,T,iT,ikT,np,nea,ne,ne0,dne,nht,nha
      double precision:: PhiH,arg,exu,err,SS,U0,ctr
      double precision:: kbcgs,ikbcgs,daux,Watom,Aatom
      double precision, dimension(:), allocatable:: frc, dfrc, popu


      ! Constants needed below
      C0 = hplanck*hplanck/2d0/PI/me/kb
      kbcgs = 1d7*kb
      ikbcgs = 1d-7/kb
      ikb = fktoJ/kb

      ! Get number of atoms in database
      natom = Atmo%nele

      !
      ! Compute sum of atomic weights and abundances
      !

      ! Initialize
      Watom = 0d0
      Aatom = 0d0

      ! For each contribution
      do ia=1,natom
        daux = Atmo%abund(ia)
        Watom = Watom + recallmass_ind(ia)*daux
        Aatom = Aatom + daux
      end do

      ! Allocate partition function
      allocate(atoms(natom))

      !
      ! Prepare atoms
      !

      ! Get partition function and compute ionization fraction
      do ia=1,natom

        atoms(ia)%s = atom_index2char(ia)
        atoms(ia)%abun = Atmo%abund(ia)
        call getpf(atom_index2char(ia),atoms(ia)%nstg,atoms(ia)%pf, &
                   atoms(ia)%Eion,Atmo)

      end do

      ! For each height
      do iz=1,nz

        ! Local temperature and proton density
        T = Atmo%T(iz)
        iT = 1d0/T
        ikT = ikb*iT
        np = Atmo%nH(iz,6)
        nht = Atmo%nht(iz)
        nha = Atmo%nha(iz)
        nea = Atmo%ne(iz)

        ! If no file and the atmospheric file did not have
        ! Hydrogen in it
        if (nlte(1).eq.0.and.Atmo%typo.ne.0) then

          ! Get partition
          U0 = atoms(1)%pf(1,iz)
          arg = U0 + atoms(1)%Eion(1)*ikT
          if (arg.lt.0d0) then
            arg = -arg
            exu = diexp(arg)
          else
            exu = ddexp(arg)
          end if
          PhiH = 0.5d6*((C0*iT)**(1.5d0))*exu
          ne0 = 0.5d0*(sqrt(1d0 + 4d0*nht*PhiH) - 1d0)/PhiH

        ! If there was either a file or data in atmosphere
        else

          ! Initialize with current electrons
          ne0 = nea

        end if ! Type of atmospheric input regarding densities

        ! Iterate
        do iter=1,maxneiter

          ! Update error and initialize sum
          err = ne0/nht
          SS = 0d0

          ! For every atomic index in the database
          do iatom=1,natom

            ! Get maximum stage, from atoms or from
            ! the atomic model
            ! LTE
            if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

              Mstg = atoms(iatom)%nstg

            ! NLTE
            else

              !
              ! Atom index

              ! Departure
              if (depar(ia).ne.0) then
                ia = abs(depar(iatom))
              else
                ia = abs(nlte(iatom))
              end if

              ! Active
              if (nlte(iatom).gt.0.or.depar(iatom).gt.0) then

                Mstg = maxval(Atom(ia)%stage)
                active = .True.

              ! Passive
              else

                Mstg = maxval(Atomb(ia)%stage)
                active = .False.

              end if ! Active/passive
            end if ! LTE/NLTE

            ! Sanity check
            if (Mstg.lt.1) cycle

            ! Check size for ionization fraction and
            ! derivative
            if (allocated(frc)) then
              if (size(frc).lt.Mstg) then
                deallocate(frc, dfrc)
                allocate(frc(Mstg))
                allocate(dfrc(Mstg))
              end if
            else
              allocate(frc(Mstg))
              allocate(dfrc(Mstg))
            end if

            ! Initialize
            frc(1:Mstg) = 0d0
            dfrc(1:Mstg) = 0d0

            ! If LTE
            if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

              ! Get ionization fraction and derivative
              call getdfrc(atoms(iatom)%nstg,atoms(iatom)%pf(:,iz), &
                           atoms(iatom)%Eion,T,ne0, &
                           frc(1:Mstg),dfrc(1:Mstg))

            ! If NLTE
            else

              ! If departure coefficients active atom
              if (depar(iatom).gt.0) then

                ! Atom index
                ia = depar(iatom)

                ! Allocate popu
                allocate(popu(Atom(ia)%nlevel))

                ! LTE and departure
                call LTEiz(Atom(ia),Atmo,iz,popu)
                popu = popu*Atom(ia)%depar(:,iz)

              ! If departure coefficients passive atom
              else if (depar(iatom).lt.0) then

                ! Atom index
                ia = depar(iatom)

                ! Allocate popu
                allocate(popu(Atomb(ia)%nlevel))

                ! LTE and departure
                call LTEiz(Atomb(ia),Atmo,iz,popu)
                popu = popu*Atomb(ia)%depar(:,iz)

              ! Active population
              else if (nlte(iatom).gt.0) then

                ! Allocate popu and copy
                allocate(popu(Atom(ia)%nlevel))
                popu = Atom(ia)%popu(:,iz)

              ! Passive population
              else

                ! Allocate popu and copy
                allocate(popu(Atomb(ia)%nlevel))
                popu = Atomb(ia)%popu(:,iz)

              end if ! Type of NLTE input and atom

              ! Total population
              U0 = 0d0

              ! Reset level index
              ilevel = 0

              ! Active
              if (active) then

                ! For each level
                do iterm=1,Atom(ia)%nmulti
                  do iJ=1,Atom(ia)%nJ(iterm)

                    ! Advance index
                    ilevel = ilevel + 1

                    ! Stage
                    istg = Atom(ia)%stage(iterm)

                    ! Add to total
                    U0 = U0 + popu(ilevel)

                    ! Add to fraction
                    frc(istg) = frc(istg) + popu(ilevel)

                  end do ! FS level
                end do ! Term

              ! Passive
              else

                ! For each level
                do iterm=1,Atomb(ia)%nmulti
                  do iJ=1,Atomb(ia)%nJ(iterm)

                    ! Advance index
                    ilevel = ilevel + 1

                    ! Stage
                    istg = Atomb(ia)%stage(iterm)

                    ! Add to total
                    U0 = U0 + popu(ilevel)

                    ! Add to fraction
                    frc(istg) = frc(istg) + popu(ilevel)

                  end do ! FS level
                end do ! Term

              end if ! Active/passive

              ! Deallocate popu
              deallocate(popu)

              ! Normalize fractions
              frc(1:Mstg) = frc(1:Mstg)/U0

            end if ! N/LTE

            ! If Hydrogen remove H minus
            if (iatom.eq.1) then

              PhiH = 0.25d6*((C0*iT)**(1.5d0))* &
                     exp(8.74980963338d3*iT)
              err = err + ne0*frc(1)*PhiH
              SS = SS - (frc(1) + ne0*dfrc(1))*PhiH

            end if ! Hydrogen

            ! For every stage but neutral
            do istg=2,Mstg

              ! Electron contribution
              ctr = dble(istg-1)*Atmo%abund(iatom)
              err = err - ctr*frc(istg)
              SS = SS + ctr*dfrc(istg)

            end do ! All ionization stages
          end do ! atom in the database

          ! Update electrons
          ne = ne0 - nht*err/(1d0 - nht*SS)

          ! Check change
          dne = abs((ne - ne0)/ne0)

          ! Move to old
          ne0 = ne

          ! Check change
          if (dne.le.epsne) exit

          ! Last iteration
          if (iter.eq.maxneiter) then
            ! If master, announce no convergence
            if (pid.eq.0) then
              write(umsg,'(A,1x,i3,1x,A,1x,es13.6,1x,A,1x,i3)') &
                ' # Warning: electron density did not converge '// &
                'after',iter,'iterations, relative change was',dne, &
                'at height',iz
              call verbose
            end if ! Master
          end if ! Last iteration
        end do ! Iterations

        Atmo%ne(iz) = ne
        Atmo%Pe(iz) = Atmo%ne(iz)*kbcgs*Atmo%T(iz)

        ! Rewrite rho
        Atmo%rho(iz) = Atmo%Pg(iz)/Atmo%T(iz)
        Atmo%rho(iz) = Atmo%rho(iz)*ikbcgs/Avog
        Atmo%rho(iz) = Atmo%rho(iz)*Watom/ &
                       (Aatom + Atmo%Pe(iz)/Atmo%Pg(iz))

      end do ! Every height

      ! Deallocate arrays
      if (allocated(atoms)) deallocate(atoms)
      if (allocated(frc)) deallocate(frc)
      if (allocated(dfrc)) deallocate(dfrc)

      call control
      return

      return

      end subroutine redo_ne

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates the array nlte and initializes it\n
      !!  nlte(integer(:)): Array that will contain information about
      !!                    population files\n
      !! depar(integer(:)): Array that will contain information about
      !!                    departure coefficient files\n
      !!  Atom(Atom_class): Structure with the atomic data\n
      !! Atomb(Atom_class): Structure with the atomic data for
      !!                    background opacities\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      subroutine initializenlte(nlte,depar,Atom,Atomb,Atmo)

      ! I/O

      integer, dimension(:), allocatable:: nlte, depar
      type(Atom_class), dimension(:), intent(inout):: Atom, Atomb
      type(Atmo_class), intent(inout):: Atmo

      ! Local

      integer:: ia

      if (.not.allocated(nlte)) then
        allocate(nlte(Atmo%nele))
        allocate(depar(Atmo%nele))
      end if
      nlte = 0
      depar = 0

      ! Active
      do ia=1,nA
        if (allocated(Atom(ia)%popu)) &
          nlte(atom_char2index(Atom(ia)%element)) = ia
        if (allocated(Atom(ia)%depar)) &
          depar(atom_char2index(Atom(ia)%element)) = ia
      end do

      ! Passive
      do ia=1,nAb
        if (allocated(Atomb(ia)%popu)) &
          nlte(atom_char2index(Atomb(ia)%element)) = -ia
        if (allocated(Atomb(ia)%depar)) &
          depar(atom_char2index(Atomb(ia)%element)) = -ia
      end do

      return

      end subroutine initializenlte

!#####################################################################
!#####################################################################
!#####################################################################

      !> Activate all active atom indexes in nlte array\n
      !!  nlte(integer(:)): Array that will contain information about
      !!                    population files\n
      !! depar(integer(:)): Array that will contain information about
      !!                    departure coefficient files\n
      !!  Atom(Atom_class): Structure with the atomic data
      subroutine activenlte(nlte,depar,Atom)

      ! I/O

      integer, dimension(:), allocatable:: nlte, depar
      type(Atom_class), dimension(:), intent(inout):: Atom

      ! Local

      integer:: ia

      ! Allocate arrays
      allocate(nlte(na),depar(na))

      ! Active
      do ia=1,nA
        nlte(atom_char2index(Atom(ia)%element)) = ia
        depar(atom_char2index(Atom(ia)%element)) = 0
      end do

      return

      end subroutine activenlte

!#####################################################################
!#####################################################################
!#####################################################################

      end module chemic_mod
