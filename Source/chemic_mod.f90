      !> Chemical equilibrium
      module chemic_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     25/03/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     25/03/2025:    V4.0.1 - The mass density was being computed
!                             following what is done in the NICOLE
!                             code, but the expression did not seem
!                             right. The determination of the mass
!                             density now follows Eq.(3-52) in Mihalas
!                             1970 (TdPA)
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
!    Calculate chemical equilibrium, set molecular populations, and
!  modify atomic populations
!
!  chemeq_T
!    Solve the system of equations for the chemical equilibrium for
!  a given temperature
!
!  eqstate
!    Solve the equation of state with the Wittmann method as in the
!  SIR code
!
!  redo_ne
!    Recalculate the electron number density from temperature and
!  other other atomic number densities (LTE if not specified in input)
!
!  initializenlte
!    Allocates and initializes the nlte and depar arrays, with
!  information on what atomic populations or departure coefficients
!  are available
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
      use parameters_mod , only : hplanck , PI , me , kb, Avog, mhm
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

      !> Calculate chemical equilibrium, set molecular populations,
      !! and modify atomic populations\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atomb(Atom_class(:)): Structures with atomic data for
      !!                             background atoms\n
      !! LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!          Mol(Mol_class(:)): Structures with molecular data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data
      subroutine chemeq(Atom,Atomb,LTElines,Mol,Atmo)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), intent(inout):: Atomb
      type(LTEline_class), dimension(:), &
                           intent(inout), allocatable:: LTElines
      type(Mol_class), dimension(:), intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo

      ! Local

      type(catm_class), dimension(nM*2):: atoms

      character(len=2):: sname

      logical:: nfound,warning,nwarning

      integer:: iz,ia
      integer:: iatom,imol,iatom1,natom,nmol,neq
      integer:: max_it_chem,check_it_chem,check_it_res,nres_opt
      integer, dimension(:,:), allocatable:: atom_index

      double precision:: phiHm,C0,frac
      double precision:: mrc_chem,minT,minT2
      double precision, dimension(:), allocatable:: xx,xxT


      ! Parameters
      max_it_chem = 500
      check_it_chem = 20
      check_it_res = 3
      nres_opt = 4
      mrc_chem = 1d-5
      minT = 2.5d3

      ! Minimum T for stable chemeq
      minT = 3.0d3
      minT2 = 2.0d3

      ! Initialize warning flags
      warning = .True.
      nwarning = .True.

      ! Constant needed below
      C0 = hplanck*hplanck/2d0/PI/me/kb

      ! If there are molecules, solve the chemical equilibrium
      if (nM.gt.0) then

        ! Local variable
        nmol = nM

        ! First atom in the system of equations is always hydrogen
        natom = 1
        atoms(1)%s = ' H'
        atoms(1)%inmod = .True.
        atoms(1)%abun = 1d0

        ! Initialize the atoms structure
        allocate(atoms(1)%imol(nM))
        atoms(1)%imol = 0
        allocate(atoms(1)%nmol(nM))
        atoms(1)%nmol = 0

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

              ! If the atom was in the list
              if (Mol(imol)%catom(iatom)%s.eq.atoms(iatom1)%s) then

                ! Add molecular information to the atom
                atoms(iatom1)%pnmol = atoms(iatom1)%pnmol + 1
                atoms(iatom1)%imol(atoms(iatom1)%pnmol) = imol
                atoms(iatom1)%nmol(atoms(iatom1)%pnmol) = &
                                                Mol(imol)%natom(iatom)

                ! Unflag not-found
                nfound = .False.

                ! And tell the molecule which atom it was
                Mol(imol)%iatom(iatom) = iatom1

                ! And stop searching
                exit

              end if ! Atom in the list

            end do ! atoms already listed

            ! If not in the list
            if (nfound) then

              !
              ! Add it
              !

              ! Advance index
              natom = natom + 1

              ! Fill information
              allocate(atoms(natom)%imol(nM))
              atoms(natom)%imol(nM) = 0
              allocate(atoms(natom)%nmol(nM))
              atoms(natom)%nmol(nM) = 0
              atoms(natom)%s = Mol(imol)%catom(iatom)%s
              atoms(natom)%pnmol = 1
              atoms(natom)%imol(1) = imol
              atoms(natom)%nmol(1) = Mol(imol)%natom(iatom)

              ! Initialize not found
              nfound = .True.

              ! In the list of active atoms
              do ia=1,nA

                ! If same atom
                if (atoms(natom)%s.eq.Atom(ia)%Element) then

                  ! Unflag not-found
                  nfound = .False.

                  ! Take abundance from model
                  atoms(natom)%inmod = .True.
                  atoms(natom)%abun = Atom(ia)%abun

                  ! Stop searching
                  exit

                end if ! Same atom

              end do ! Active atoms

              ! If not found yet
              if (nfound) then

                ! In the list of passive atoms
                do ia=1,nAb

                  ! If same atom
                  if (atoms(natom)%s.eq.Atomb(ia)%Element) then

                    ! Unflag not found
                    nfound = .False.

                    ! Take abundance from model
                    atoms(natom)%inmod = .True.
                    atoms(natom)%abun = Atomb(ia)%abun

                    ! Stop searching
                    exit

                  end if ! Same atom

                end do ! Passive atoms

              end if ! Not found yet

              ! If not found
              if (nfound) then

                ! Flag no model and look in the tabulation
                atoms(natom)%inmod = .False.
                atoms(natom)%abun = &
                    Atmo%abund(atom_char2index(atoms(natom)%s))

              end if ! Atom not found

              ! Tell the molecule which atom it was
              Mol(imol)%iatom(iatom) = natom

            end if ! The atom was in the list of atoms

          end do ! Components of the molecule
        end do ! Molecule

        !
        ! Construct the system of equations
        !

        ! Number of equations is number of atoms in molecules, plus
        ! Hydrogen (if not in molecule) and Molecules
        neq = natom + nM

        ! Allocate and initialize solution arrays
        ! Solution
        allocate(xx(neq),xxT(neq))
        xxT = 0d0

        ! Allocate and initialize indexing of the atoms in the system
        allocate(atom_index(natom,2))
        atom_index = -1

        !
        ! Index the variables in the system to retrieve later
        !

        ! For each atom involved in the chemical equilibrium
        do iatom=1,natom

          ! Look for it between the active atoms
          do ia=1,nA

            ! If element involved
            if (Atom(ia)%Element.eq.atoms(iatom)%s) then

              ! Save index and type
              atom_index(iatom,1) = ia
              atom_index(iatom,2) = 1
              exit

            end if ! Element involved

          end do ! Active atoms

          ! Look for it between the passive atoms
          do ia=1,nAb

            ! If element involved
            if (Atomb(ia)%Element.eq.atoms(iatom)%s) then

              ! Save index and type
              atom_index(iatom,1) = ia
              atom_index(iatom,2) = 2
              exit

            end if ! Element involved

          end do ! Passive atoms
        end do ! Atoms in the system of equations


        !
        ! Solve the chemical equilibrium
        !

        ! For each height
        do iz=1,nz

          ! If temperature too small
          if (Atmo%T(iz).lt.minT) then

            ! Solve for minT
            call chemeq_T(Atom,Atomb,Mol,Atmo,natom,nmol,neq, &
                          atoms,atom_index,warning,nwarning, &
                          max_it_chem,check_it_chem,check_it_res, &
                          nres_opt,mrc_chem,c0,iz,.True., &
                          minT,xx,xxT)

            ! If temperature even smaller
            if (Atmo%T(iz).lt.minT2) then

              ! Solve for T
              call chemeq_T(Atom,Atomb,Mol,Atmo,natom,nmol,neq, &
                            atoms,atom_index,warning,nwarning, &
                            max_it_chem,check_it_chem,check_it_res, &
                            nres_opt,mrc_chem,c0,iz,.False., &
                            minT2,xxT,xx)

              ! Update initial
              xxT = xx

            end if

            ! Solve for actual T
            call chemeq_T(Atom,Atomb,Mol,Atmo,natom,nmol,neq, &
                          atoms,atom_index,warning,nwarning, &
                          max_it_chem,check_it_chem,check_it_res, &
                          nres_opt,mrc_chem,c0,iz,.False., &
                          Atmo%T(iz),xxT,xx)

          ! Normal temperature
          else

            ! Solve for T
            call chemeq_T(Atom,Atomb,Mol,Atmo,natom,nmol,neq, &
                          atoms,atom_index,warning,nwarning, &
                          max_it_chem,check_it_chem,check_it_res, &
                          nres_opt,mrc_chem,c0,iz,.True., &
                          Atmo%T(iz),xxT,xx)

          end if


          !
          ! Use the solution to affect the populations of atoms
          ! and molecules
          !

          ! For each atom
          do iatom=1,natom

            ! Check if it has model
            if (atom_index(iatom,2).gt.0) then

              ! Get its index
              ia = atom_index(iatom,1)

              ! If it is active
              if (atom_index(iatom,2).eq.1) then

                ! If protected
                if (Atom(ia)%mol_protect) then

                  ! Verbose if master in 1D synthesis and first time
                  if (iz.eq.1.and.pid.eq.0.and.run_mode.eq.0) then
                    umsg = ' - '//atoms(iatom)%s// &
                           ' is in molecules, but its '// &
                           ' read population was kept'
                    call verbose
                  end if

                  ! Skip correction
                  cycle

                end if ! Protected

                ! Calculate and apply fraction of population
                frac = xx(iatom)/Atom(ia)%n(iz)
                Atom(ia)%popu(:,iz) = Atom(ia)%popu(:,iz)*frac
                Atom(ia)%populte(:,iz) = Atom(ia)%populte(:,iz)*frac
                Atom(ia)%n(iz) = xx(iatom)

              ! If it is a background atom
              else if (atom_index(iatom,2).eq.2) then

                ! If protected
                if (Atomb(ia)%mol_protect) then

                  ! Verbose if master in 1D synthesis and first time
                  if (iz.eq.1.and.pid.eq.0.and.run_mode.eq.0) then
                    umsg = ' - Atom '//atoms(iatom)%s// &
                           ' is in molecules, but the '// &
                           'read populations were kept'
                    call verbose
                  end if

                  ! Skip correction
                  cycle

                end if ! Protected

                ! Calculate and apply fraction of population
                frac = xx(iatom)/Atomb(ia)%n(iz)
                Atomb(ia)%popu(:,iz) = Atomb(ia)%popu(:,iz)*frac
                Atomb(ia)%populte(:,iz) = Atomb(ia)%populte(:,iz)*frac
                Atomb(ia)%n(iz) = xx(iatom)

              end if ! Active or background atom

              ! If master in 1D synthesis and first time, verbose
              if (iz.eq.1.and.gpid.eq.0.and.run_mode.eq.0) then
                umsg = ' - Population of '//atoms(iatom)%s// &
                       ' modified by molecules'
                call verbose
              end if

            end if ! If the atom has a model

          end do ! For all involved atoms

          ! If LTE lines and allocated
          if (nLTEl.gt.0.and.allocated(LTElines)) then

            ! For each LTE line
            do ia=1,nLTEl

              ! Get name of element
              sname = atom_index2char(LTElines(ia)%ele)

              ! For each atom in chemical eq.
              do iatom=1,natom

                ! If same atom
                if (atoms(iatom)%s.eq.sname) then

                  ! Substitute population
                  LTElines(ia)%n(iz) = xx(iatom)

                end if

              end do ! Atoms in chemical eq.
            end do ! LTE lines

          end if ! LTE lines


          !
          ! Store molecular density
          !
          do imol=1,nM
            Mol(imol)%n(iz) = xx(imol + natom)
          end do

        end do ! Heights

      ! If no molecules, calculate just Hminus
      else

        ! For every height
        do iz=1,nz

          ! H- constant
          phiHm = .25d6*((C0/Atmo%T(iz))**(1.5d0))* &
                  exp(8.74980963338d3/Atmo%T(iz))

          ! Get number density
          Atmo%nHm(iz) = Atmo%ne(iz)*Atmo%nHT(iz)*phiHm

        end do ! Heights

      end if ! If there are molecules

      ! If 1D synthesis
      if (run_mode.eq.0) then

        ! For each molecule
        do imol=1,nM

          ! Free molecular data
          deallocate(Mol(imol)%eqcoeff)

        end do

      end if ! 1D

      ! Control
      call control

      return

      end subroutine chemeq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the system of equations for the chemical equilibrium
      !! for a given temperature\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!            natom(integer): Number of different atoms\n
      !!  atom_index(integer(:,:)): Indexing of atoms\n
      !!             nmol(integer): Number of different molecules\n
      !!              neq(integer): Number of equations\n
      !!      atoms(catm_class(:)): Atomic data for equations\n
      !!          warning(logical): If to issue warning for all\n
      !!         nwarning(logical): If to issue warning for solution
      !!                            issue\n
      !!      max_it_chem(integer): Maximum number of iterations\n
      !!    check_it_chem(integer): Iterations to check physics\n
      !!     check_it_res(integer): Number of resets allowed\n
      !!         nres_opt(integer): Type of reset to try first\n
      !!          mrc_chem(dfloat): Maximum relative change to
      !!                            achieve\n
      !!                C0(dfloat): Constant needed in equations\n
      !!               iz(integer): Height index\n
      !!             diss(logical): If to initialize full
      !!                            dissociation\n
      !!                 T(double): Temperature\n
      !!               xx0(double): Initial solution\n
      !!                xx(double): Solution
      subroutine chemeq_T(Atom,Atomb,Mol,Atmo,natom,nmol, &
                          neq,atoms,atom_index,warning,nwarning, &
                          max_it_chem,check_it_chem,check_it_res, &
                          nres_opt,mrc_chem,c0,iz,diss,T,xx0,xx)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), intent(inout):: Atomb
      type(Mol_class), dimension(:), intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(catm_class), dimension(:), intent(inout):: atoms
      logical, intent(in):: diss
      logical, intent(inout):: warning,nwarning
      integer, intent(in):: natom,nmol,neq,iz
      integer, intent(in):: max_it_chem,check_it_chem
      integer, intent(in):: check_it_res,nres_opt
      integer, dimension(:,:), intent(in):: atom_index
      double precision, intent(in):: mrc_chem,C0,T
      double precision, dimension(:), intent(in):: xx0
      double precision, dimension(:), intent(out):: xx

      ! Local

      type(Atom_class):: Atom_l

      logical:: nfound

      integer:: ia,iterm,iJ,ilevel,iter,info,nres
      integer:: iatom,imol,ieq
      integer, dimension(neq):: indx

      double precision:: phiHm,aaHm,phil
      double precision:: minv,maxv
      double precision, dimension(neq):: aa,bb,phi
      double precision, dimension(nmol):: MEq
      double precision, dimension(2,neq):: frc
      double precision, dimension(neq,neq):: aap


      ! For each molecule
      do imol=1,nmol

        ! Get Eq. constant
        MEq(imol) = Moleq_T(Mol(imol),T)

      end do ! Molecules

      ! For all involved atoms atoms
      do iatom=1,natom

        ! Get the partition function of this atom
        call getpf_T(atoms(iatom)%s,atoms(iatom)%nstg, &
                     atoms(iatom)%pfs,atoms(iatom)%Eion, &
                     Atmo,T)

      end do ! Involved atoms


      !
      ! Construct the system of equations
      !

      ! Initialize
      frc = 0d0
      bb = 0d0
      nres = 0

      ! For each atom in the SoE
      do iatom=1,natom

        !
        ! Independent element
        !

        ! If Hydrogen
        if (iatom.eq.1) then

          ! If there is a model
          if (atoms(iatom)%inmod) then

            ! If active atom
            if (atom_index(iatom,2).eq.1) then

              ! Get model
              Atom_l = Atom(atom_index(iatom,1))

            ! If passive atom
            else if (atom_index(iatom,2).eq.2) then

              ! Get model
              Atom_l = Atomb(atom_index(iatom,1))

            ! No model
            else

              ! Critical error
              umsg = 'Something went wrong when '// &
                     'indexing the atoms in chemical '// &
                     'equilibrium'
              urou = 'chemeq_T'
              call aborted
              return

            end if ! Type of model

            ! If protected
            if (Atom_l%mol_protect) then

              ! Take population from model
              bb(iatom) = Atom_l%n(iz)

            ! If not protected
            else

              ! Take population from atmosphere
              bb(iatom) = Atmo%nHt(iz)

            end if ! Protected

          ! If no model
          else

            ! Get population from atmosphere
            bb(iatom) = Atmo%nHt(iz)

          end if ! If there is model

        ! Any other element
        else

          ! Get population from abundance and atmospheric hydrogen
          bb(iatom) = atoms(iatom)%abun*Atmo%nHT(iz)

        end if ! Hydrogen or other element


        !
        ! Get fraction in neutral and first ion
        !

        ! Initialize the not found flag
        nfound = .True.

        ! If the atom has a model
        if (atoms(iatom)%inmod) then

          ! If active atom
          if (atom_index(iatom,2).eq.1) then

            ! Get model
            Atom_l = Atom(atom_index(iatom,1))

          ! If passive atom
          else if (atom_index(iatom,2).eq.2) then

            ! Get model
            Atom_l = Atomb(atom_index(iatom,1))

          ! No model
          else

            ! Critical error
            umsg = 'Something went wrong when '// &
                   'indexing the atoms in chemical '// &
                   'equilibrium'
            urou = 'chemeq_T'
            call aborted
            return

          end if ! Type of model

          ! If atom has neutral stage
          if (Atom_l%stage(1).eq.1) then

            ! Unflag not found
            nfound = .False.

          ! The atom does not have neutral stage
          else

            ! The master issues a warning the first time
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

        ! Missing the necessary populations
        if (nfound) then

          ! Use partition functions to calculate the fractions
          call getfrc(atoms(iatom)%nstg,atoms(iatom)%pfs(:), &
                      atoms(iatom)%Eion,T,Atmo%ne(iz), &
                      -1,frc(:,iatom))

        ! If we have the populations, use them
        else

          ! Initialize level index
          ilevel = 0

          ! For each term
          do iterm=1,Atom_l%nMulti

            ! Stop if beyond the single ionized species
            if (Atom_l%stage(iterm).gt.2) exit

            ! For each sublevel
            do iJ=1,Atom_l%nJ(iterm)

              ! Advance index
              ilevel = ilevel + 1

              ! Add population to fraction (numerator)
              frc(Atom_l%stage(iterm),iatom) = &
                              frc(Atom_l%stage(iterm),iatom) + &
                              Atom_l%popu(ilevel,iz)

            end do ! Sublevels
          end do ! Term

          ! Compute fraction
          frc(:,iatom) = frc(:,iatom)/Atom_l%n(iz)

        end if ! If populations available or not

      end do ! atoms in molecules

      ! Add Hminus
      phiHm = .25d6*((C0/T)**(1.5d0))*exp(8.74980963338d3/T)
      aaHm = Atmo%ne(iz)*frc(1,1)*phiHm

      ! Get equilibrium constant for every molecule
      do imol=1,nmol
        phi(imol) = MEq(imol)
      end do

      ! If initializing dissociated
      if (diss) then

        ! Initial solution, complete dissociation
        xx = bb

      ! Not initializing dissociated
      else

        ! Take provided solution
        xx = xx0

      end if


      !
      ! Force protections
      !

      ! Atoms in SoE
      do iatom=1,natom

        ! If has model
        if (atoms(iatom)%inmod) then

          ! Active atom
          if (atom_index(iatom,2).eq.1) then

            ! Get model
            Atom_l = Atom(atom_index(iatom,1))

          ! Passive atom
          else if (atom_index(iatom,2).eq.2) then

            ! Get model
            Atom_l = Atomb(atom_index(iatom,1))

          ! No model
          else

            ! Issue error
            umsg = 'Something went wrong when '// &
                   'indexing the atoms in chemical '// &
                   'equilibrium'
            urou = 'chemeq_T'
            call aborted
            return

          end if ! Type of model
        end if ! Has model

        ! If protected, get population from model
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
        do imol=1,nmol

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

          !
          ! Derivative

          ! Atoms in molecule
          do ia=1,Mol(imol)%nA

            ! Ask which atom it is in the list
            iatom = Mol(imol)%iatom(ia)

            ! Contribution due to molecule in the atomic row
            aap(iatom,ieq) = aap(iatom,ieq) + &
                             Mol(imol)%natom(ia)
            ! Contribution due to atom in the molecule row
            aap(ieq,iatom) = -phil*Mol(imol)%natom(ia)/xx(iatom)

          end do ! Atom in molecule

        end do ! Molecules

        !
        ! Force protections
        !

        ! Atoms in SoE
        do iatom=1,natom

          ! If has model
          if (atoms(iatom)%inmod) then

            ! Active atom
            if (atom_index(iatom,2).eq.1) then

              ! Get model
              Atom_l = Atom(atom_index(iatom,1))

            ! Passive atom
            else if (atom_index(iatom,2).eq.2) then

              ! Get model
              Atom_l = Atomb(atom_index(iatom,1))

            ! No model
            else

              ! Issue warning
              umsg = 'Something went wrong when '// &
                     'indexing the atoms in chemical '// &
                     'equilibrium'
              urou = 'chemeq_T'
              call aborted
              return

            end if ! Type of model
          end if ! Has model

          ! If protected
          if (Atom_l%mol_protect) then

            ! Trivial row in system of equations
            aap(iatom,:) = 0d0
            aap(iatom,iatom) = 1d0
            aa(iatom) = 0d0

          end if ! Protected

        end do ! Atoms in SoE

        ! Solve system of equations
        call DGESV(neq,1,aap,neq,indx,aa,neq,info)

        ! Get the new solution
        xx = xx - aa

        ! And calculate the MRC
        do ieq=1,neq
          if (abs(xx(ieq)).gt.1d-30) aa(ieq) = aa(ieq)/xx(ieq)
        end do

        ! If we reached convergence, leave
        if (maxval(abs(aa)).lt.mrc_chem.and.minval(xx).ge.0d0) &
          exit

        ! If maximum iterations, can warn, and master
        if (iter.eq.max_it_chem.and.nwarning.and.pid.eq.0) then

          ! If negative solution
          if (minval(xx).lt.0d0) then

            ! Issue error
            umsg = 'Negative fractions obtained in '// &
                   'chemical equilibrium'
            urou = 'chemeq_T'
            call abortedS(umsg,urou,.True.,.True.)
            nwarning = .False.
            warning = .False.
            cycle

          ! If not in inversion and can warn
          else if (warning.and.run_mode.ne.-1) then

            ! Notify warning
            umsg = 'Maximum number of iterations reached '//&
                   'in the chemical equilibrium'
            urou = 'chemeq_T'
            call abortedS(umsg,urou,.False.,.True.)
            warning = .False.
            cycle

          end if ! Negative solution or limit iterations

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
              do imol=1,nmol

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
              do imol=1,nmol

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
              minv = minval(xx(natom+1:neq), &
                            xx(natom+1:neq).gt.0d0)
              ! Get maximum non-zero
              maxv = maxval(xx(natom+1:neq))

              ! For each molecule
              do imol=1,nmol

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
              do imol=1,nmol

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
            do imol=1,nmol

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
            do imol=1,nmol

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
      ! Store Hm density in atmospheric model
      !
      Atmo%nHm(iz) = Atmo%ne(iz)*xx(1)*phiHm

      ! Free auxiliar variables
      do iatom=1,natom
        deallocate(atoms(iatom)%pfs,atoms(iatom)%Eion)
      end do

      return

      end subroutine chemeq_T

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the equation of state with the Wittmann method as in
      !! the SIR code\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      nlte(integer(:)): Array with information about
      !!                        loaded populations\n
      !!     depar(integer(:)): Array with information about
      !!                        loaded departure coefficients
      subroutine eqstate(Atmo,Atom,Atomb,nlte,depar)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom,Atomb
      type(Atmo_class), intent(inout):: Atmo
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      type(catm_class), dimension(:), allocatable:: atoms

      integer:: iz,ia,natom,iter,istg

      double precision:: kbcgs,ikbcgs,Watom,C0,daux


      ! Inverse of Boltzmann constant in cgs
      kbcgs = kb*1d7
      ikbcgs = 1d-7/kb

      ! Density constant
      C0 = ikbcgs*mhm/Avog/recallmass_ind(1)

      ! Compute sum of atomic weights and abundances
      Watom = 0d0
      natom = Atmo%nele

      ! Allocate partition function
      allocate(atoms(natom))

      ! For every atom
      do ia=1,natom

        ! Get name, abundance, and partition functions
        atoms(ia)%s = atom_index2char(ia)
        atoms(ia)%abun = Atmo%abund(ia)
        call getpf(atoms(ia)%s, atoms(ia)%nstg, atoms(ia)%pf, &
                   atoms(ia)%Eion, Atmo)

        ! For each height
        do iz=1,nz

          ! For each stage
          do istg=1,atoms(ia)%nstg

            ! Make partition function linear
            atoms(ia)%pf(istg,iz) = exp(atoms(ia)%pf(istg,iz))

          end do ! Stages
        end do ! Heights

        ! Convert energy to eV
        atoms(ia)%Eion = atoms(ia)%Eion*fktoev

      end do ! Atoms


      !
      ! Allocations
      !

      ! Electron Pressure
      if (.not.allocated(Atmo%Pe)) then
        allocate(Atmo%Pe(nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pe)
        Atmo%Pe = 0d0
      end if

      ! Gas Pressure
      if (.not.allocated(Atmo%Pg)) then
        allocate(Atmo%Pg(nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%Pg)
        Atmo%Pg = 0d0
      end if

      ! Mass density
      if (.not.allocated(Atmo%rho)) then
        allocate(Atmo%rho(nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%rho)
        Atmo%rho = 0d0
      end if

      ! For each atomic contribution
      do ia=1,natom

        ! Get abundance
        daux = Atmo%abund(ia)

        ! Get mass times abundance
        Watom = Watom + recallmass_ind(ia)*daux

      end do ! Atomic contributions

      ! If only here because you want to output an atmospheric file
      if (Atmo%typo.eq.0) then

        ! Compute electron pressure
        Atmo%Pe = Atmo%ne*kbcgs*Atmo%T

        ! Solve equation of state
        call eqstate_known(Atmo,Atom,Atomb,nlte,depar,atoms)

        ! Preliminar rho
        Atmo%rho = Atmo%Pe*C0/Atmo%T

        ! For each height
        do iz=1,nz

          ! Get mass density
          daux = Watom*Atmo%nHa(iz)/Atmo%ne(iz)
          Atmo%rho(iz) = Atmo%rho(iz)*daux

        end do ! Heights

      ! If an electronic quantity is given
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
          call abortedS(umsg,urou,.True.,.True.)
          call control
          return

        end if ! Type of electron input

        ! Solve equation of state
        call eqstate_ele(Atmo,Atom,Atomb,nlte,depar,atoms)

        ! Preliminar rho
        Atmo%rho = Atmo%Pe*C0/Atmo%T

        ! For each height
        do iz=1,nz

          ! Get mass density
          daux = Watom*Atmo%nHa(iz)/Atmo%ne(iz)
          Atmo%rho(iz) = Atmo%rho(iz)*daux

        end do ! Heights

      ! If gas pressure
      else if (Atmo%typo.eq.4) then

        ! Solve equation of state
        call eqstate_gas(Atmo,Atom,Atomb,nlte,depar,atoms)

        ! Preliminar rho
        Atmo%rho = Atmo%Pe*C0/Atmo%T

        ! For each height
        do iz=1,nz

          ! Get mass density
          daux = Watom*Atmo%nHa(iz)/Atmo%ne(iz)
          Atmo%rho(iz) = Atmo%rho(iz)*daux

        end do ! Heights

      ! If mass density
      else if (Atmo%typo.eq.5) then

        ! C0 is used with the inverse
        C0 = 1d0/C0

        ! Stimate gass pressure
        Atmo%Pe = Atmo%rho*Atmo%T*C0

        ! Iterate
        do iter=1,maxiter

          ! Solve equation of state
          call eqstate_ele(Atmo,Atom,Atomb,nlte,depar,atoms)

          ! Preliminar Pe
          Atmo%Pe = Atmo%rho*Atmo%T*C0

          ! For each height
          do iz=1,nz

            ! Compute electron pressure
            daux = Watom*Atmo%nHa(iz)/Atmo%ne(iz)
            Atmo%Pe(iz) = Atmo%Pe(iz)/daux

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

      !> Recalculate the electron number density from temperature and
      !! other atomic number densities (LTE if not specified in
      !! input)\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!      nlte(integer(:)): Array with information about loaded
      !!                        populations\n
      !!     depar(integer(:)): Array with information about loaded
      !!                        departure coefficients\n
      !!      Atmo(Atmo_class): Structure with atmospheric data
      subroutine redo_ne(Atom,Atomb,nlte,depar,Atmo)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom,Atomb
      type(Atmo_class), intent(inout):: Atmo
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      type(catm_class), dimension(:), allocatable:: atoms

      logical:: active

      integer:: natom,iatom,ia,istg,iterm,iJ,ilevel,iter,iz,Mstg

      double precision:: C0,C1,ikb,T,iT,ikT,np,nea,ne,ne0,dne,nht,nha
      double precision:: PhiH,arg,exu,err,SS,U0,ctr
      double precision:: kbcgs,ikbcgs,daux,Watom
      double precision, dimension(:), allocatable:: frc, dfrc, popu


      ! Constants needed below
      C0 = hplanck*hplanck/2d0/PI/me/kb
      kbcgs = 1d7*kb
      ikbcgs = 1d-7/kb
      ikb = fktoJ/kb

      ! Density constant
      C1 = ikbcgs*mhm/Avog/recallmass_ind(1)

      ! Get number of atoms in database
      natom = Atmo%nele

      !
      ! Compute sum of atomic weights and abundances
      !

      ! Initialize
      Watom = 0d0

      ! For each contribution
      do ia=1,natom

        ! Sum mass and mass times abundance
        daux = Atmo%abund(ia)
        Watom = Watom + recallmass_ind(ia)*daux

      end do ! Atomic contributions

      ! Allocate partition function
      allocate(atoms(natom))

      !
      ! Prepare atoms
      !

      ! For each atom
      do ia=1,natom

        ! Get name, abundance, and partition function
        atoms(ia)%s = atom_index2char(ia)
        atoms(ia)%abun = Atmo%abund(ia)
        call getpf(atom_index2char(ia),atoms(ia)%nstg,atoms(ia)%pf, &
                   atoms(ia)%Eion,Atmo)

      end do ! Atoms

      ! For each height
      do iz=1,nz

        ! Local temperature and densities
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

          ! Get electron from hydrogen partition function
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

            !
            ! Get maximum stage, from atoms or from
            ! the atomic model

            ! Full LTE
            if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

              ! Number of stages is known
              Mstg = atoms(iatom)%nstg

            ! NLTE
            else

              !
              ! Atom index

              ! Departure
              if (depar(ia).ne.0) then

                ! Get index
                ia = abs(depar(iatom))

              ! NLTE
              else

                ! Get index
                ia = abs(nlte(iatom))

              end if ! Departure or NLTE populations

              ! Active
              if (nlte(iatom).gt.0.or.depar(iatom).gt.0) then

                ! Get maximum stage from atomic model
                Mstg = maxval(Atom(ia)%stage)
                active = .True.

              ! Passive
              else

                ! Get maximum stage from atomic model
                Mstg = maxval(Atomb(ia)%stage)
                active = .False.

              end if ! Active/passive
            end if ! LTE/NLTE

            ! Sanity check, skip if no stages
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

            ! If full LTE
            if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

              ! Get ionization fraction and derivative from
              ! partition functions
              call getdfrc(atoms(iatom)%nstg,atoms(iatom)%pf(:,iz), &
                           atoms(iatom)%Eion,T,ne0, &
                           frc(1:Mstg),dfrc(1:Mstg))

            ! If NLTE
            else

              ! If departure coefficients in active atom
              if (depar(iatom).gt.0) then

                ! Atom index
                ia = depar(iatom)

                ! Allocate popu
                allocate(popu(Atom(ia)%nlevel))

                ! Calculate LTE and apply departure
                call LTEiz(Atom(ia),Atmo,iz,popu)
                popu = popu*Atom(ia)%depar(:,iz)

              ! If departure coefficients passive atom
              else if (depar(iatom).lt.0) then

                ! Atom index
                ia = depar(iatom)

                ! Allocate popu
                allocate(popu(Atomb(ia)%nlevel))

                ! Calculate LTE and apply departure
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

            ! If Hydrogen
            if (iatom.eq.1) then

              ! Remove H-
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

          ! Check change is small enough to finish
          if (dne.le.epsne) exit

          ! Last iteration
          if (iter.eq.maxneiter) then

            ! If master
            if (pid.eq.0) then

              ! Announce no convergence
              write(umsg,'(A,1x,i3,1x,A,1x,es13.6,1x,A,1x,i3)') &
                ' # Warning: electron density did not converge '// &
                'after',iter,'iterations, relative change was',dne, &
                'at height',iz
              call verbose

            end if ! Master
          end if ! Last iteration
        end do ! Iterations

        ! Save electron density and pressure
        Atmo%ne(iz) = ne
        Atmo%Pe(iz) = Atmo%ne(iz)*kbcgs*Atmo%T(iz)

        ! Rewrite rho
        Atmo%rho(iz) = Atmo%Pe(iz)*C1/Atmo%T(iz)
        Atmo%rho(iz) = Atmo%rho(iz)*Watom*Atmo%NHa(iz)/Atmo%ne(iz)

      end do ! Every height

      ! Deallocate arrays
      if (allocated(atoms)) deallocate(atoms)
      if (allocated(frc)) deallocate(frc)
      if (allocated(dfrc)) deallocate(dfrc)

      ! Control
      call control

      return

      end subroutine redo_ne

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes the nlte and depar arrays, with
      !! information on what atomic populations or departure
      !! coefficients are available\n
      !!     nlte(integer(:)): Array that will contain information
      !!                       about population files\n
      !!    depar(integer(:)): Array that will contain information
      !!                       about departure coefficient files\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !! Atomb(Atom_class(:)): Structures with atomic data for
      !!                       background atoms\n
      !!     Atmo(Atmo_class): Structure with atmospheric data\n
      subroutine initializenlte(nlte,depar,Atom,Atomb,Atmo)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom,Atomb
      type(Atmo_class), intent(in):: Atmo
      integer, dimension(:), allocatable, intent(inout):: nlte,depar

      ! Local

      integer:: ia

      ! If arrays not allocated, do it
      ! Not added to memory count because they live in the
      ! prepare_syn() subroutine
      if (.not.allocated(nlte)) then
        allocate(nlte(Atmo%nele))
        allocate(depar(Atmo%nele))
      end if

      ! Initialize arrays
      nlte = 0
      depar = 0

      ! Run over active atoms
      do ia=1,nA

        ! If there are populations already, signal it in nlte
        if (allocated(Atom(ia)%popu)) &
          nlte(atom_char2index(Atom(ia)%element)) = ia

        ! If there are populations already, signal it in depar
        if (allocated(Atom(ia)%depar)) &
          depar(atom_char2index(Atom(ia)%element)) = ia

      end do ! Active atoms

      ! Run over background atoms
      do ia=1,nAb

        ! If there are populations already, signal it in nlte
        if (allocated(Atomb(ia)%popu)) &
          nlte(atom_char2index(Atomb(ia)%element)) = -ia

        ! If there are populations already, signal it in depar
        if (allocated(Atomb(ia)%depar)) &
          depar(atom_char2index(Atomb(ia)%element)) = -ia

      end do ! Passive atoms

      return

      end subroutine initializenlte

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set all active atoms as available in terms of populations\n
      !!     nlte(integer(:)): Array with information about
      !!                       loaded populations\n
      !!    depar(integer(:)): Array with information about
      !!                       loaded departure coefficients
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      subroutine activenlte(nlte,depar,Atom)

      ! I/O

      integer, dimension(:), allocatable, intent(inout):: nlte,depar
      type(Atom_class), dimension(:), intent(in):: Atom

      ! Local

      integer:: ia

      ! Allocate arrays
      allocate(nlte(na),depar(na))

      ! Run over active atoms
      do ia=1,nA

        ! Signal nlte and not departure
        nlte(atom_char2index(Atom(ia)%element)) = ia
        depar(atom_char2index(Atom(ia)%element)) = 0

      end do ! Active atoms

      return

      end subroutine activenlte

!#####################################################################
!#####################################################################
!#####################################################################

      end module chemic_mod
