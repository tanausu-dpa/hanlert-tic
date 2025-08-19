      !> Reading atomic data
      module ratom_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     18/04/2017
!  Last version:
!     18/08/2025 V4.0.7
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/08/2025:    V4.0.7 - The non-magnetic indexes need to be
!                             initialized even if all heights have
!                             a magnetic field if we are using the
!                             two step solution (TdPA)
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
!  ratom
!    Read a model atom from the specified file
!
!  AtomH
!    Generate a 6 levels hydrogen atomic model without fine structure
!
!  abund
!    Normalize the abundances of several atoms of the same species
!
!  setFScoeff
!    Calculate the fine-structure Einstein coefficients for the given
!  atom
!
!  setup_LTE_transition
!    Prepare the known quantities for a given LTE transition
!
!  remove_LTE_transition
!    Check if an LTE transition is part of an atomic model and remove
!  it from the model
!
!  allocateatom
!    Allocate an array of atomic models
!
!  set_atom_label
!    Provide each atomic model with a unique label
!
!  shiftatoms
!    Shift all atoms one position to the right in the array, as well
!  as their associated filenames
!
!  prepareatom
!    Allocate height dependent pupulation variables
!
!  prepareatomMPI
!    Allocate Master MPI variables for each model atom
!
!  initMblock
!    Initialize the size of the M blocks assuming there is a magnetic
!  field
!
!  set_atom_indexes
!    Setup all necessary atomic indexings for transitions components
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use funnj_mod
      use parameters_mod , only : c, c2, dopp, PI, ConvF , TINYSP
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read a model atom from the specified file\n
      !!    filename(character(:)): Name of the file to read\n
      !!      source(character(:)): Path to the source code folder\n
      !!          ID(character(:)): ID of this run\n
      !!        skip_wave(logical): If this atom does not contribute
      !!                            to the wavelength axis\n
      !!  Kcut_input(integer(:,:)): Term wise K cut input data\n
      !!          Atom(Atom_class): Structure with atomic data\n
      !!             indx(integer): Index of the currect atom in the
      !!                            list of atoms\n
      !!            isPRD(logical): If the current atom has at least
      !!                            one PRD line\n
      !!           active(logical): If the current atom is active
      subroutine rAtom(filename,source,ID,skip_wave,Kcut_input, &
                       Atom,indx,isPRD,active)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: filename,source
      character(len=9), intent(in):: ID
      logical, intent(in):: active,skip_wave
      logical, intent(inout):: isPRD
      integer, intent(in):: indx
      integer, dimension(:,:), allocatable, intent(in):: Kcut_input

      ! Local

      character(len=1):: cdump

      logical:: flin,changed

      integer:: ios,iterm,iterm1,ilevel,ilevel1
      integer:: ii,jj,i,itran,ftran,icol,stagl,stagu
      integer:: K,iQ,iJ,iJ1,iJl,iJu,up,low,Tindex,minK,maxK,maxfst
      integer:: l1,l2,l3,l4,nTmp,nclin,cind,col_type,nion

      double precision, parameter:: dLSJ = .25d0

      double precision:: d1,d2,d3,rL,S,rJ,rJ1,rJmin,rJmax
      double precision:: freq,ftmp,radAcoeff
      double precision, dimension(:), allocatable:: vaux1,vaux2

      ! Pointer

      type(Tbox_class), pointer:: p_T


      ! Initialize
      nullify(p_T,Atom%Ccoeff_special,Atom%Tbox)

      ! Master
      if(pid.eq.0) then

        ! Verbose
        umsg = ' - Read atom '//trim(filename)
        call verbose

        ! Translate the atomic model in python
        call system('python '//trim(source)//'ratom.py '// &
                    trim(filename)//' '//ID//' '//verbosef)

      end if ! Master

      ! Wait for the master to finish
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! Open translated file
      open(100,file='tmp_atom_'//ID,status='old',iostat=ios,err=1000)

      ! Success
      read (100,*,err=1100) ios

      ! If no correct file
      if (ios.lt.0) then

        ! Issue error
        umsg = 'Problem translating the atomic file '//trim(filename)
        goto 1200

      end if ! Wrong file

      ! Add to memory
      MRAMc = MRAMc + 1d-6*sizeof(Atom)

      ! Read element name
      read (100,*,err=1100) Atom%Element

      ! Fix element name spaces
      if (Atom%Element(2:2).eq.' ') then
        Atom%Element(2:2) = Atom%Element(1:1)
        Atom%Element(1:1) = ' '
      end if

      ! Atomic mass
      read (100,*,err=1100) Atom%rmass

      ! Atomic abundance (12 + log(Atom/H))
      read (100,*,err=1100) Atom%abun

      ! Transform to 12 + log10()
      Atom%abun = 1d1**(Atom%abun - 12d0)

      ! Multiplicative modifier to this abundance
      read (100,*,err=1100) Atom%abun_mod

      ! If to normalize abundance multiplier to 1 for same ID atoms
      read (100,*,err=1100) ios

      ! If atom is active, transform last entry to bool
      if (active) Atom%anorm = ios.eq.1

      ! Read type of model (multi-term or multi-level)
      read (100,*,err=1100) l1
      Atom%ML = l1.eq.1

      ! Read dimensions of the model: Number or terms, number of
      ! transitions, number of photoionization transitions, number
      ! of elastic collisional rates
      read (100,*,err=1100) Atom%nMulti,Atom%ntran,Atom%nphot,Atom%ngk

      ! Read maximum number of levels per term
      read (100,*,err=1100) Atom%nJmax

      ! Maximum vectorial dimension for levels
      Atom%NNN = Atom%nMulti*Atom%nJmax


      !!!!!!!!!!!!!!!!!!!!!
      !
      ! Level information
      !
      !!!!!!!!!!!!!!!!!!!!!

      !
      ! Allocations
      !

      ! Orbital angular momentum
      allocate(Atom%rLval(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rLval)
      ! Spin angular momentum
      allocate(Atom%Sval(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%Sval)
      ! Number of FS levels
      allocate(Atom%nJ(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%nJ)
      ! Ionization stage
      allocate(Atom%stage(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%stage)
      ! Angular momentum
      allocate(Atom%rJval(Atom%nJmax,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rJval)
      Atom%rJval = -1d0
      ! Frequency of FS levels
      allocate(Atom%FSfreq(Atom%nJmax,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%FSfreq)
      ! Frequency of term
      allocate(Atom%TRfreq(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%TRfreq)
      ! Term degeneracy
      allocate(Atom%deg(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%deg)
      ! If multi-level atom
      if (Atom%ML) then
        ! Lande factor
        allocate(Atom%gL(Atom%nMulti))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%gL)
      end if

      ! Initialize level index
      l3 = 0

      ! For each term
      do iterm=1,Atom%nMulti

        ! Read orbital and spin ang. momenta and ionization stage
        read (100,*,err=1100) rL,S,l1

        ! Update maximum value of quantum numbers
        if (nint(2d0*S).gt.nxS) nxS = nint(2d0*S)
        if (nint(2d0*rL).gt.nxL) nxL = nint(2d0*rL)

        ! Store in atomic structure
        Atom%rLval(iterm) = rL
        Atom%Sval(iterm) = S
        Atom%stage(iterm) = l1

        ! If multilevel
        if (Atom%ML) then

          ! Only one level per term is possible
          Atom%nJ(iterm) = 1

          ! Read angular momentum, energy, Landé factor, term index,
          ! and level index
          read (100,*,err=1100) rJ,freq,d1,l1,l2

          ! Advance level index
          l3 = l3 + 1

          ! Sanity check read
          if (l1.ne.iterm.or.l2.ne.l3) then

            ! If read term index is not current index
            if (l1.ne.iterm) then

              ! Issue term error
              umsg = 'The term index does not '// &
                     'correspond to the term'

            ! The error must be the level instead
            else

              ! Issue level error
              umsg = 'The level index does not '// &
                     'correspond to the level'

            end if ! Wrong term or level index

            ! Jump to abortion
            goto 1200

          end if ! Sanity check

          ! Energy units
          freq = freq*1d-5

          ! Store in atomic structure
          Atom%rJval(1,iterm) = rJ
          Atom%deg(iterm) = 2d0*rJ + 1d0
          Atom%FSfreq(1,iterm) = freq
          Atom%TRfreq(iterm) = freq

          ! Update maximum total angular momenta
          if (nint(2d0*rJ).gt.nxJ) nxJ = nint(2d0*rJ)

          ! If input Landé factor non-physical
          if (d1.lt.0d0) then

            ! If J is 0
            if (rJ.lt.dLSJ) then

              ! Landé factor must be zero as well
              Atom%gL(iterm) = 1d0

            ! If J is not 0
            else

              ! Assume LS coupling
              Atom%gL(iterm) = 1d0 + .5d0*(rJ*(rJ+1d0) + S*(S+1d0) - &
                                           rL*(rL+1d0))/rJ/(rJ+1d0)

            end if ! J value

          ! Physical Landé factor
          else

            ! Save in structure
            Atom%gl(iterm) = d1

          end if ! Input Landé factor

        ! If multiterm
        else

          ! Calculate degeneration of term
          Atom%deg(iterm) = (2d0*rL + 1d0)*(2d0*S + 1d0)

          ! Minimum and maximum angular momentums
          rJmin = abs(rL - S)
          rJmax = rL + S

          ! Number of FS levels for this term
          Atom%nJ(iterm) = nint(rJmax - rJmin) + 1

          ! Initialize sum
          ftmp = 0d0

          ! For each level within the term
          do iJ=1,Atom%nJ(iterm)

            ! Read angular momentum, energy, term index and level
            ! index
            read (100,*,err=1100) rJ,freq,l1,l2

            ! Advance level index
            l3 = l3 + 1

            ! Sanity check read
            if (l1.ne.iterm.or.l2.ne.l3) then

              ! If read term index is not current index
              if (l1.ne.iterm) then

                ! Issue term error
                umsg = 'The term index does not '// &
                       'correspond to the term'

              ! The error must be the level instead
              else

                ! Issue level error
                umsg = 'The level index does not '// &
                       'correspond to the level'

              end if ! Wrong term or level index

              ! Jump to abortion
              goto 1200

            end if ! Sanity check

            ! Energy units
            freq = freq*1d-5

            ! Store in atomic structure
            Atom%rJval(iJ,iterm) = rJ
            Atom%FSfreq(iJ,iterm) = freq

            ! Update maximum total angular momenta
            if (nint(2d0*rJ).gt.nxJ) nxJ = nint(2d0*rJ)

            ! Add weighted contribution to the term energy
            ftmp = (2d0*rJ + 1d0)*freq + ftmp

          end do ! Levels within term

          ! Calculate the multiplet energy as a weighted average of
          ! the F-levels energies
          Atom%TRfreq(iterm) = ftmp/Atom%deg(iterm)

        end if ! Multilevel or multiterm

      end do ! Term

      ! Sanity check energies
      if (maxval(Atom%TRfreq).lt.1d-4) then

        ! Issue warning
        write(umsg,'(A)') ' # Maximum energy in '// &
                          Atom%element// &
                          ' model atom is smaller than 10 cm^-1, '// &
                          'you may be using an old atomic model'// &
                          ' format'
        call verbose

      end if ! Suspicious energy values

      ! Determine the number of FS levels in this atom
      Atom%nlevel = sum(Atom%nJ)

      ! Calculate the maximum K multipole
      Atom%nKmax = nint(maxval(Atom%rJval)*2d0)

      ! Update max K value if active atom
      if (Atom%nKmax.gt.nkx.and.active) nkx = Atom%nKmax

      ! Calculate the maximum number of magnetic sublevels for this
      ! atom
      Atom%nMmax = nint(maxval(Atom%rJval)*2d0 + 1d0)

      ! Allocate atomic quantities to be used later
      allocate(Atom%nblk(Atom%nMmax,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%nblk)
      allocate(Atom%iJval(Atom%nJmax,Atom%nMmax,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%iJval)

      ! Term given the level
      allocate(Atom%term(Atom%nlevel))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%term)
      ! Sublevel in its term given the global level
      allocate(Atom%sublevel(Atom%nlevel))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%sublevel)
      ! Allocate indexing by term
      allocate(Atom%irho(Atom%nMulti))

      !
      ! Index levels and sublevels
      !

      ! Initialize level counter
      ii = 0

      ! For each term
      do iterm=1,Atom%nMulti

        ! Add memory in indexing
        MRAMc = MRAMc + 1d-6*sizeof(Atom%irho(iterm))

        ! Allocate rho indexing
        allocate(Atom%irho(iterm)%irho_ij(Atom%nJ(iterm)))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%irho(iterm)%irho_ij)

        ! For each level
        do iJ=1,Atom%nJ(iterm)

          ! Advance index
          ii = ii + 1

          ! Index term--sublevel--level
          Atom%irho(iterm)%irho_ij(iJ) = ii
          Atom%term(ii) = iterm
          Atom%sublevel(ii) = iJ

        end do ! Levels
      end do ! Terms

      ! If active atom
      if (active) then

        ! Add contribution to the number of transitions (b-b and b-f)
        nxtran = nxtran + Atom%ntran
        nxphot = nxphot + Atom%nphot

      end if ! Active atom

      ! Initialize continuous index shift
      Atom%tshift = 0
      Atom%tfshift = 0
      Atom%pshift = 0


      !!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Radiative transitions
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!

      !
      ! Allocations
      !

      ! Structure for FS information
      allocate(Atom%fst(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst)
      ! Energy of the term-term transition
      allocate(Atom%Dfreq(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%Dfreq)
      ! Type of elastic broadening
      allocate(Atom%broad_type(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%broad_type)
      ! Arguments to calculate the broadening
      allocate(Atom%broad_args(4,Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%broad_args)
      ! Argument for Stark broadening
      allocate(Atom%broad_stark(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%broad_stark)
      ! Number of frequencies to use in the line
      allocate(Atom%nfreqt(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%nfreqt)
      ! How many of those frequencies for the core region
      allocate(Atom%nfreqtc(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%nfreqtc)
      ! Doppler widths that the line holds
      allocate(Atom%dwvl(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%dwvl)
      ! Dopple widths that the core of the line holds
      allocate(Atom%dwvlc(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%dwvlc)
      ! Flag of second order emissivity
      allocate(Atom%lemiss2(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%lemiss2)
      ! Flag of Split components
      allocate(Atom%splitf(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%splitf)
      ! Flag of line presence
      allocate(Atom%fflag(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fflag)
      ! Lower frequency limits
      allocate(Atom%if0(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%if0)
      ! Upper frequency limits
      allocate(Atom%if1(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%if1)
      ! Lower absolute frequency limits
      allocate(Atom%tif0(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%tif0)
      ! Upper absolute frequency limits
      allocate(Atom%tif1(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%tif1)
      ! Lower absolute frequency limits
      allocate(Atom%rif0(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rif0)
      ! Upper absolute frequency limits
      allocate(Atom%rif1(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rif1)
      ! Lower frequency weight
      allocate(Atom%W0(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%W0)
      ! Upper frequency weight
      allocate(Atom%W1(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%W1)
      ! Indexing of radiative transitions
      allocate(Atom%irad(Atom%nMulti,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%irad)
      Atom%irad = 0
      ! Einstein coefficients Aul and Blu
      allocate(Atom%Ecoeff(Atom%nMulti,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%Ecoeff)
      Atom%Ecoeff = 0d0

      ! Reset counters
      Atom%nfreq = 0
      Atom%nftran = 0

      ! For each transition
      do itran=1,Atom%ntran

        ! Read the terms involved and Aul
        read (100,*,err=1100) iterm1,iterm,radAcoeff

        ! Rate units
        radAcoeff = radAcoeff*1d-8

        ! Ensure correct term ordering (iterm1 upper)
        if (iterm1.lt.iterm) then
          l1 = iterm1
          iterm1 = iterm
          iterm = l1
        end if

        ! Enumerate permitted transitions
        Atom%irad(iterm1,iterm) = itran
        Atom%irad(iterm,iterm1) = itran
        Atom%fst(itran)%iterml = iterm
        Atom%fst(itran)%itermu = iterm1

        ! Calculate the frequency of the transition
        Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                                Atom%TRfreq(iterm))

        ! Read collisional broadening parameters
        read (100,*,err=1100) l1
        Atom%broad_type(itran) = l1
        do l1=1,4
          read (100,*,err=1100) d1
          Atom%broad_args(l1,itran) = d1
        end do

        ! Calculate Aul and Blu
        Atom%Ecoeff(iterm1,iterm) = radAcoeff
        Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                    Atom%deg(iterm1)/ &
                                    Atom%deg(iterm)/ &
                                   (ConvF*Atom%Dfreq(itran)* &
                                    1d21*(2d0*c)* &
                                    Atom%Dfreq(itran)**2d0)


        ! Read the Stark broadening parameter, the number of
        ! frequencies for this transition and its core, the
        ! Doppler widths of this transition and its core, the
        ! flag of second order emissivity, split between components
        read (100,*,err=1100) d3,l1,l2,d1,d2,l3,l4

        ! Store Stark broadening information
        Atom%broad_stark(itran) = d3

        ! Ensure odd numbers
        if (mod(l1,2).eq.0) l1 = l1 + 1
        if (mod(l2,2).eq.0) l2 = l2 + 1

        ! Store frequency data in atomic structure
        Atom%nfreqt(itran) = l1
        Atom%nfreqtc(itran) = l2
        Atom%Dwvl(itran) = d1
        Atom%Dwvlc(itran) = d2

        ! If PRD flag
        if (l3.ne.0) then

          ! This transition is PRD
          Atom%lemiss2(itran) = .True.
          isPRD = .True.

        ! No PRD flag
        else

          ! This transition is CRD
          Atom%lemiss2(itran) = .False.

        end if ! PRD flag

        ! Store the split components flag
        Atom%splitf(itran) = l4.ne.0

        ! Allocate FS transition indexing
        allocate(Atom%fst(itran)%irad(Atom%nJ(iterm1), &
                                      Atom%nJ(iterm)))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
        Atom%fst(itran)%irad = 0
        ! Allocate Aul for FS transitions
        allocate(Atom%fst(itran)%Aul(Atom%nJ(iterm1),Atom%nJ(iterm)))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%Aul)
        Atom%fst(itran)%Aul = 0d0
        ! Allocate Blu for FS transitions
        allocate(Atom%fst(itran)%Blu(Atom%nJ(iterm),Atom%nJ(iterm1)))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%Blu)
        Atom%fst(itran)%Blu = 0d0

        ! If Multi-level
        if (Atom%ML) then

          ! Trivial values in FS indexing
          Atom%fst(itran)%irad(1,1) = 1
          Atom%fst(itran)%nt = 1
          Atom%nftran = Atom%nftran + 1
          maxfst = 1
          l3 = 1
          allocate(Atom%fst(itran)%ilevell(l3))
          allocate(Atom%fst(itran)%ilevelu(l3))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
          MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
          Atom%fst(itran)%ilevell = 1
          Atom%fst(itran)%ilevelu = 1

        ! If multi-term
        else

          ! Initialize maximum FS trans
          maxfst = 0

          !
          ! Count FS transitions
          !

          ! Initialize FS counter
          l3 = 0

          ! For each upper level
          do iJu=1,Atom%nJ(iterm1)

            ! For each lower level
            do iJl=1,Atom%nJ(iterm)

              ! Check electric dipole selection rules
              if(abs(Atom%rJval(iJu,iterm1) - &
                     Atom%rJval(iJl,iterm)).gt.1 .or. &
                 Atom%rJval(iJu,iterm1) + &
                 Atom%rJval(iJl,iterm) .lt. .4d0)cycle

              ! Advance index and store
              l3 = l3 + 1
              Atom%fst(itran)%irad(iJu,iJl) = l3

            end do ! Lower levels
          end do ! Upper levels

          ! Store number of transitions
          Atom%fst(itran)%nt = l3

          ! Update maximum number of FS transitions
          if (l3.gt.maxfst) maxfst = l3

          ! Add to total number of FS transition in this atom
          Atom%nftran = Atom%nftran + l3

          ! Allocate and initialize sublevel indexing
          allocate(Atom%fst(itran)%ilevell(l3))
          allocate(Atom%fst(itran)%ilevelu(l3))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
          MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
          Atom%fst(itran)%ilevell = -1
          Atom%fst(itran)%ilevelu = -1

          !
          ! Index
          !

          ! Initialize FS counter
          l3 = 0

          ! For each upper level
          do iJu=1,Atom%nJ(iterm1)

            ! For each lower level
            do iJl=1,Atom%nJ(iterm)

              ! Check electric dipole selection rules
              if(abs(Atom%rJval(iJu,iterm1) - &
                     Atom%rJval(iJl,iterm)).gt.1 .or. &
                 Atom%rJval(iJu,iterm1) + &
                 Atom%rJval(iJl,iterm) .lt. .4d0)cycle

              ! Advance index and store
              l3 = l3 + 1
              Atom%fst(itran)%ilevell(l3) = iJl
              Atom%fst(itran)%ilevelu(l3) = iJu

            end do ! Lower levels
          end do ! Upper levels

        end if ! Multi-level/Multi-term

        ! Background atoms can stop here
        if (.not.active) cycle

        ! Add the frequency nodes to the total count of this atom
        Atom%nfreq = Atom%nfreq + l1*l3

      end do ! Transitions


      ! If an active atom
      if (active) then

        ! Update the total maximum number of FS transitions
        nxt = nxt + Atom%nftran

        !
        ! Index the FS transitions given the term transition and the
        ! sub FS transition
        !

        ! Allocate indexes
        allocate(Atom%ifst_ij(maxfst,Atom%ntran))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%ifst_ij)
        Atom%ifst_ij = 0
        allocate(Atom%ifst(Atom%nftran))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%ifst)
        Atom%ifst = 0
        allocate(Atom%ifstj(Atom%nftran))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%ifstj)
        Atom%ifstj = 0

        ! Reset counter
        ftran = 0

        ! For each transition
        do itran=1,Atom%ntran

          ! For each FS transition
          do l3=1,Atom%fst(itran)%nt

            ! Advance index
            ftran = ftran + 1

            ! Store indexes
            Atom%ifst_ij(l3,itran) = ftran
            Atom%ifst(ftran) = itran
            Atom%ifstj(ftran) = l3

          end do ! FS transitions
        end do ! Term transition

        ! Add the frequency nodes of this atom to the total count
        ! if not neglecting its wavelengths
        if (.not.skip_wave) nfreq = nfreq + Atom%nfreq

        !
        ! Control K cuts
        !

        ! Term-wise K cut
        allocate(Atom%Kcut(Atom%nMulti))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%Kcut)
        Atom%Kcut = Kcut
        ! Transition-wise Krad cut
        allocate(Atom%Krad(Atom%ntran))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%Krad)
        Atom%Krad = Krad

        ! If there are inputs, we need to take them into account
        if (allocated(Kcut_input)) then

          ! Initialize flag
          changed = .False.

          ! For each entry
          do ii=1,size(Kcut_input,2)

            ! If atom is indexed
            if (indx.eq.Kcut_input(1,ii)) then

              ! Sanity check
              if (Kcut_input(2,ii).lt.0.or. &
                  Kcut_input(3,ii).gt.Atom%nMulti) then

                ! Abort
                umsg = 'Wrong index in Kcut for atom '// &
                       trim(Atom%Element)
                urou = 'rAtom'
                call gaborted

              end if ! Sanity check

              ! Flag changed
              changed = .True.

              ! For specified range
              do jj=Kcut_input(2,ii),Kcut_input(3,ii)

                ! For this term, apply K cut
                Atom%Kcut(jj) = min(Kcut_input(4,ii), &
                                   nint(2d0*maxval(Atom%rJval(:,jj))))

              end do ! Specified range

            end if ! Atom is indexed

          end do ! Input entries

          ! If there were changes in the Kcut, control Krad as well
          if (changed) then

            ! For each transition
            do itran=1,Atom%ntran

              ! Get terms
              iterm = Atom%fst(itran)%iterml
              iterm1 = Atom%fst(itran)%itermu

              ! If there is PRD
              if (isPRD) then

                ! Check if PRD transition
                changed = Atom%lemiss2(itran)

                ! If no PRD, check all
                if (.not.changed) then

                  ! Check all terms
                  do icol=1,Atom%nMulti

                    ! Check transition with same upper term
                    ftran = Atom%irad(icol,iterm1)

                    ! Cycle null
                    if (ftran.le.0) cycle

                    ! If PRD
                    if (Atom%lemiss2(ftran)) then
                      changed = .True.
                      exit
                    end if

                  end do ! All terms

                end if ! Current transition no PRD

                ! If PRD in upper term
                if (changed) then

                  ! Check cut
                  if (Atom%Kcut(iterm1).ne.Kcut) then

                    ! Abort
                    umsg = 'Please, do not add specific K '// &
                           'limits to terms involved in '// &
                           'PRD transitions'
                    urou = 'rAtom'
                    call gaborted

                  end if ! Edited K cut

                  ! Check all terms
                  do icol=1,Atom%nMulti

                    ! Check transition with same upper term
                    ftran = Atom%irad(icol,iterm1)

                    ! Cycle null
                    if (ftran.le.0) cycle

                    ! Check cut
                    if (Atom%Kcut(icol).ne.Kcut) then

                      ! Abort
                      umsg = 'Please, do not add specific K '// &
                             'limits to terms involved in '// &
                             'PRD transitions'
                      urou = 'rAtom'
                      call gaborted

                    end if ! Edited K cut

                  end do ! All terms

                end if ! Upper level in PRD transition
              end if ! There is PRD at all

              ! Set maximum Krad
              Atom%Krad(itran) = min(Atom%Kcut(iterm) + &
                                     Atom%Kcut(iterm1),2)

            end do ! Transitions

          end if ! Changes in term-wise K cut
        end if ! There are custom K cut


        !
        ! Define the index matrix
        !

        ! Initialize running index
        i = 0

        ! For each term
        do iterm=1,Atom%nMulti

          ! Allocate J-J index
          allocate(Atom%irho(iterm)%Jrho(Atom%nJ(iterm), &
                                         Atom%nJ(iterm)))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%irho(iterm)%Jrho)

          ! For each level
          do iJ=1,Atom%nJ(iterm)

            ! Get J
            rJ = Atom%rJval(iJ,iterm)

            ! For each other level
            do iJ1=1,Atom%nJ(iterm)
           !do iJ1=iJ,iJ ! no J,J'

              ! Get J
              rJ1 = Atom%rJval(iJ1,iterm)

              ! Get K limits
              minK = nint(abs(rJ-rJ1))
              maxK = min(nint(rJ+rJ1),Atom%Kcut(iterm))

              ! Allocate and initialize KQ tab
              allocate(Atom%irho(iterm)%Jrho(iJ1,iJ)% &
                            kq(-maxK:maxK,minK:maxK))
              MRAMc = MRAMc + &
                      1d-6*sizeof(Atom%irho(iterm)%Jrho(iJ1,iJ)%kq)
              Atom%irho(iterm)%Jrho(iJ1,iJ)%kq = 0

              ! For each K
              do K=minK,maxK

                ! For each Q
                do iQ=-K,K

                  ! Advance index
                  i = i + 1

                  ! Store index
                  Atom%irho(iterm)%Jrho(iJ1,iJ)%kq(iQ,K) = i

                end do ! Q
              end do ! K
            end do ! J'
          end do ! J
        end do ! Term

        ! Update maximum dimensionality of SEE
        Atom%ndim = i

        ! Verbose SEE dimension
        if (pid.eq.0) then
          write(umsg,'(A,i6)') '   S.E. dimension =',Atom%ndim
          call verbose
        end if

      end if ! Active atom

      ! If there are transitions
      if (Atom%ntran.gt.0) then

        ! Initialize maximum
        ftmp = 0d0

        ! For each term
        do iterm=1,Atom%nMulti-1

          ! Get maximum Einstein coefficient for this term
          ftmp = max(ftmp, &
                     maxval(Atom%Ecoeff(iterm+1:Atom%nMulti,iterm)))

        end do ! Terms

        ! If too small Einstein coefficient
        if (ftmp.lt.1d-9.and.ftmp.gt.0d0) then

          ! Issue warning due to wrong format
          write(umsg,'(A)') ' # Maximum Einstein '// &
              'coefficient in '//Atom%element// &
              ' model atom is smaller than 0.1 s^-1, '// &
              'you may be using an old atomic model format'
          call verbose

        end if ! Too small Esintein coefficient
      end if ! There are transitions


      !!!!!!!!!!!!!!!!!!!!!!
      !
      ! Elastic collisions
      !
      !!!!!!!!!!!!!!!!!!!!!!

      ! If there are inputs
      if (Atom%ngk.ge.1) then

        ! Allocate elastic data
        allocate(Atom%elas(Atom%ngk))

        ! For each entry
        do ii=1,Atom%ngk

          ! Count memory
          MRAMc = MRAMc + 1d-6*sizeof(Atom%elas(ii))

          ! Read level and number of entries
          read (100,*,err=1100) ilevel,l1

          ! Store in elastic structure
          Atom%elas(ii)%ilevel = ilevel
          Atom%elas(ii)%nentry = l1

          ! Allocate entries
          allocate(Atom%elas(ii)%datum(l1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%elas(ii)%datum)

          ! For each line in this entry
          do jj=1,l1

            ! Read the multipole, the type of input, and the
            ! dimensionality
            read (100,*,err=1100) K
            read (100,*,err=1100) cdump
            read (100,*,err=1100) l2

            ! Parse into database
            Atom%elas(ii)%datum(jj)%K = K
            Atom%elas(ii)%datum(jj)%nz = l2

            ! If it is a fit input
            if(cdump.eq.'f')then

              ! Read the coefficients
              read (100,*,err=1100) d1,d2,d3

              ! Put in database
              Atom%elas(ii)%datum(jj)%typo = 0
              Atom%elas(ii)%datum(jj)%a = d1
              Atom%elas(ii)%datum(jj)%b = d2
              Atom%elas(ii)%datum(jj)%c = d3

            ! No type of input recognized
            else

              ! Issue error
              umsg = 'Mode of elastic collisions '// &
                     'not recognized'
              goto 1200

            end if ! Type of elastic input

          end do ! Input sub entry
        end do ! Input entry

        ! If passive, forget everything
        if (.not.active) then

          ! For each entry
          do ii=1,Atom%ngk

            ! Deallocate subentry
            MRAMc = MRAMc - 1d-6*(sizeof(Atom%elas(ii)) + &
                                  sizeof(Atom%elas(ii)%datum))
            deallocate(Atom%elas(ii)%datum)

          end do ! Entries

          ! Deallocate structure
          deallocate(Atom%elas)

        end if ! Passive

      ! There is no input
      else

        ! If master and active atom
        if(pid.eq.0.and.active) then

          ! Verbose
          write(umsg,'(A)') '   No elastic collisions'
          call verbose

        end if ! Master and active atom
      end if ! Elastic collisional rates input


      !!!!!!!!!!!!!!!!!!!!
      !
      ! Photoionizations
      !
      !!!!!!!!!!!!!!!!!!!!

      !
      ! Allocations
      !

      ! Indexing of b-f transitions
      allocate(Atom%iphot(Atom%nlevel,Atom%nlevel))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%iphot)
      Atom%iphot = 0

      ! If there are inputs
      if (Atom%nphot.ge.1) then

        ! Allocate photoionization information
        allocate(Atom%phot(Atom%nphot))

        ! For each input entry
        do ii=1,Atom%nphot

          ! Memory count
          MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii))

          ! Read the levels involved
          read (100,*,err=1100) ilevel1,ilevel

          ! Ensure ilevel1 > ilevel
          if (ilevel1.lt.ilevel) then
            l1 = ilevel1
            ilevel1 = ilevel
            ilevel = l1
          end if

          ! Find the term and level indexes
          iterm = Atom%term(ilevel)
          iJ = Atom%sublevel(ilevel)
          iterm1 = Atom%term(ilevel1)
          iJ1 = Atom%sublevel(ilevel1)

          ! Check that they have actually different ionization stages
          if (Atom%stage(iterm1).eq.Atom%stage(iterm)) then

            ! Issue error
            umsg = 'Photoionization between terms in '// &
                   'the same stage.'
            goto 1200

          end if ! Same stage

          ! Index the b-f transition
          Atom%iphot(ilevel1,ilevel) = ii
          Atom%iphot(ilevel,ilevel1) = ii
          Atom%phot(ii)%ilevell = ilevel
          Atom%phot(ii)%ilevelu = ilevel1

          ! Calculate the edge of the photoionization
          Atom%phot(ii)%edge = Atom%FSfreq(iJ1,iterm1) - &
                               Atom%FSfreq(iJ,iterm)

          ! Store the ratio between the degeneration of the levels
          Atom%phot(ii)%glu = (2d0*Atom%rJval(iJ,iterm)+1d0)/ &
                              (2d0*Atom%rJval(iJ1,iterm1)+1d0)

          ! Read the type of input and the number of frequencies to
          ! use for this b-f transition
          read (100,*,err=1100) cdump
          read (100,*,err=1100) l1

          ! Save in structure
          Atom%phot(ii)%nfreq = l1

          ! If the input is explicit
          if (cdump.eq.'e') then

            !
            ! Allocate explicit data
            !

            ! Frequencies in the input
            allocate(Atom%phot(ii)%infreq(l1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii)%infreq)
            ! Cross-sections in the input
            allocate(Atom%phot(ii)%inalpha(l1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii)%inalpha)

            ! Flag explicit
            Atom%phot(ii)%mode = 0

            ! For each entry
            do jj=1,l1

              ! Read wavelength and cross-section
              read(100,*,err=1100) d1,d2

              ! Store frequency (kaiser) and cross-section
              Atom%phot(ii)%infreq(jj) = 1d2/d1
              Atom%phot(ii)%inalpha(jj) = d2

            end do ! Entries

            !
            ! Check order
            !

            ! If more than one entry
            if (l1.gt.1) then

              ! If wrong order, revert it
              if (Atom%phot(ii)%infreq(1).gt. &
                  Atom%phot(ii)%infreq(2)) then

                ! Check not already allocated
                if (.not.allocated(vaux1)) then
                  allocate(vaux1(l1))
                  allocate(vaux2(l1))
                else
                  if (l1.gt.size(vaux1)) then
                    deallocate(vaux1)
                    deallocate(vaux2)
                    allocate(vaux1(l1))
                    allocate(vaux2(l1))
                  end if
                end if

                ! Save auxiliar
                vaux1 = Atom%phot(ii)%infreq
                vaux2 = Atom%phot(ii)%inalpha

                ! For each element
                do jj=1,l1

                  ! Reverse order
                  Atom%phot(ii)%infreq(jj) = vaux1(l1 - jj + 1)
                  Atom%phot(ii)%inalpha(jj) = vaux2(l1 - jj + 1)

                end do ! Elements

              end if ! Wrong order
            end if ! More than one element

          ! If the input is hydrogenic
          else if (cdump.eq.'h') then

            !
            ! Allocate hydrogenic data
            !

            ! Maximum frequency to reach
            allocate(Atom%phot(ii)%infreq(1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii)%infreq)
            ! Cross-section at edge
            allocate(Atom%phot(ii)%inalpha(1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii)%inalpha)

            ! Flag hydrogenic
            Atom%phot(ii)%mode = 1

            ! Read the maximum frequency and cross-section at edge
            read(100,*,err=1100) d1,d2

            ! Save frequency (kaiser) and cross-section
            Atom%phot(ii)%infreq(1) = 1d2/d1
            Atom%phot(ii)%inalpha(1) = d2

          end if ! Type of input

          ! If active and not neglecting the wavelengths of this atom,
          ! add the nodes to the total count
          if (active.and..not.skip_wave) &
            nfreq = nfreq + Atom%phot(ii)%nfreq + 1

        end do ! Photoionization inputs

      ! If there is no input
      else

        ! If master and active
        if (pid.eq.0.and.active) then

          ! Verbose
          write(umsg,'(A)') '   No photoionizations'
          call verbose

        end if ! Master and active

        ! Allocate dummy structure
        allocate(Atom%phot(1))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(1))

      end if ! There is photoionization inputs


      !!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Inelastic collisions
      !
      !!!!!!!!!!!!!!!!!!!!!!!!

      ! Read number of collisional entries
      read (100,*,err=1100) Atom%ncol

      ! If there are inputs
      if(Atom%ncol.gt.0)then

        ! Allocate database
        allocate(Atom%inelas(Atom%ncol))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas)

        ! Flag for forbidden collisions
        allocate(Atom%fcflag(Atom%nlevel,Atom%nlevel))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%fcflag)
        Atom%fcflag = 0

        ! Nullify Tbox pointer
        nullify(Atom%Tbox)

        ! Read number of entries (not the same that number of
        ! collisions)
        read (100,*,err=1100) nclin

        ! Reset index
        icol = 0

        ! Index of temperature block
        Tindex = 0

        ! For each entry
        do ii=1,nclin

          ! Check type of entry
          read (100,*,err=1100) cind

          ! It is a temperature entry
          if (cind.lt.0) then

            ! Read the type of the collisions to come, the type of
            ! interpolation, and the number of temperatures in the
            ! input table
            read (100,*,err=1100) col_type
            read (100,*,err=1100) nion
            read (100,*,err=1100) nTmp

            ! Memory count
            MRAMc = MRAMc + 1d-6*dble(4*5 + 8*nTmp)

            ! Advance the temperature index
            Tindex = Tindex + 1

            ! If structure is initialized
            if (associated(Atom%Tbox)) then

              ! Add new box
              p_T => Atom%Tbox
              do while (associated(p_T%next))
                p_T => p_T%next
              end do
              allocate(p_T%next)
              p_T => p_T%next
              nullify(p_T%next)

             ! Structure is empty
            else

              ! Initialize structure
              allocate(Atom%Tbox)
              p_T => Atom%Tbox
              nullify(p_T%next)

            end if ! Structure initialized

            ! Put data in the box
            p_T%nTmp = nTmp
            p_T%ind = Tindex
            p_T%col_type = col_type
            p_T%nion = nion
            allocate(p_T%temp(nTmp))

            ! Read temperatures
            read (100,*,err=1100) p_T%temp

            ! Determine if interpolation
            flin = nion.eq.2.or.nion.eq.3.or.nion.eq.-2

            ! Complete the box
            p_T%flin = flin

            ! Nullify
            nullify(p_T)

          ! If it is a collisional entry
          else if (cind.ge.0) then

            ! Run up the index
            icol = icol + 1

            ! Memory count
            MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(icol))

            ! Read if forbidden
            read (100,*,err=1100) l3

            ! Read the terms involced in the collision
            read (100,*,err=1100) up,low

            ! Ensure up > low
            if (low.gt.up) then
              l2 = up
              up = low
              low = l2
            end if

            ! Store in database
            Atom%inelas(icol)%ind = Tindex
            Atom%inelas(icol)%col_type = cind
            Atom%inelas(icol)%up = up
            Atom%inelas(icol)%low = low
            Atom%inelas(icol)%forbid = l3

            ! Allocate space
            allocate(Atom%inelas(icol)%Cul(nTmp))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(icol)%Cul)

            ! Read the data in the table
            read (100,*,err=1100) Atom%inelas(icol)%Cul

            ! If between terms
            if (col_type.eq.0) then

              ! Stages
              stagl = Atom%stage(low)
              stagu = Atom%stage(up)

            ! If between levels
            else if (col_type.eq.1) then

              ! Stages
              stagl = Atom%stage(Atom%term(low))
              stagu = Atom%stage(Atom%term(up))

            end if ! If between terms or levels


            !
            ! If it is a symmetric b-b excitation
            !
            if (cind.eq.0.or.cind.eq.2.or.cind.eq.3) then

              ! If different stages
              if (stagu.ne.stagl) then

                ! Issue error
                umsg = 'Exciting collisions cannot '// &
                       'change ion'
                goto 1200

              end if ! Different stages


            !
            ! If it is a symmetric b-f excitation
            !
            else if (cind.eq.1.or.cind.eq.4) then

              ! If term-wise (can only be level-wise)
              if (col_type.eq.0) then

                ! Issue error
                umsg = 'Ionizing collisions only admit '// &
                       'level wise rates'
                goto 1200

              end if ! Term-wise

              ! If same stage
              if (stagu.eq.stagl) then

                ! Issue error
                umsg = 'Ionizing collisions must '// &
                       'change ion'
                goto 1200

              end if ! Same stage


            !
            ! Charge transfer
            !
            else if (cind.eq.5.or.cind.eq.6) then

              ! If term-wise (can only be level wise)
              if (col_type.eq.0) then

                ! Issue error
                umsg = 'Charge transfer collisions only admit '// &
                       'level wise rates'
                goto 1200

              end if ! Term-wise

              ! If same stage
              if (stagu.eq.stagl) then

                ! Issue error
                umsg = 'Charge transfer collisions must '// &
                       'change ion'
                goto 1200

              end if ! Same stage
            end if ! Type of collision
          end if ! Collisional rate input

        end do ! Collision input lines in text

      ! No collisions
      else

        ! If Master and active atom
        if (active.and.pid.eq.0) then

          ! Verbose
          write(umsg,'(A)') '   No collisions'
          call verbose

        end if ! Master and active atom
      end if ! There are collisional inputs

      ! Calculate the atomic part of the Doppler width
      Atom%cDopp = dopp/sqrt(Atom%rmass)

      ! Close file
      close (100)

      ! Control that everything went fine
      call control

      ! Master, remove temporal file
      if (pid.eq.0) call system('rm tmp_atom_'//ID)

      ! If multilevel, trick the quantum numbers
      if (Atom%ML) then

        ! Make spin equal to zero
        Atom%Sval = 0d0

        ! For each term
        do iterm=1,Atom%nMulti

          ! L = J
          Atom%rLval(iterm) = Atom%rJval(1,iterm)

        end do ! Terms

      end if ! Multi-level atom

      ! If background atom
      if (.not.active) then

        ! Deallocate unecessary data if not active atom
        deallocate(Atom%nfreqt)
        deallocate(Atom%nfreqtc)
        deallocate(Atom%lemiss2)
        deallocate(Atom%fflag)
        MRAMc = MRAMc - 1d-6*(sizeof(Atom%nfreqt) + &
                              sizeof(Atom%nfreqtc) + &
                              sizeof(Atom%lemiss2) + &
                              sizeof(Atom%fflag))

      end if ! Background atom

      return

1000  umsg = 'Error opening atom file'
      urou = 'rAtom'
      call gaborted
1100  umsg = 'Error reading atom file'
1200  close(100)
      urou = 'rAtom'
      call gaborted

      end subroutine rAtom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate a 6 levels hydrogen atomic model without fine
      !! structure\n
      !!    Atom(Atom_class): Structure with the atomic data
      subroutine AtomH(Atom)

      ! I/O

      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iterm,iterm1,ii,itran,iJ,nTmp

      double precision:: radAcoeff
      double precision, dimension(:), allocatable:: deg

      ! Pointer

      type(Tbox_class), pointer:: p_T


      ! Initialize pointer
      nullify(p_T,Atom%Ccoeff_special)

      ! If master
      if (pid.eq.0) then

        ! Informative message
        umsg = ' - Could not find hydrogen atom, loading '// &
               'default with 6-n levels for the background'
        call verbose

      end if ! Master

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Atom)

      ! Atomic ID
      Atom%Element = ' H'
      ! Identify this as a custom atom
      Atom%cust = .True.
      ! Atomic mass
      Atom%rmass = 1.00794d0
      ! Abundance modifier
      Atom%abun_mod = 1d0
      ! Abundance with respect to Hydrogen
      Atom%abun = 1d0

      ! Number of 'terms'
      Atom%nMulti = 6
      ! Number of b-b transitions
      Atom%ntran = 10
      ! Number of b-f transitions
      Atom%nphot = 5
      ! Number of elastic collisions
      Atom%ngk = 0
      ! Maximum of sublevels
      Atom%nJmax = 1
      ! Maximum of levels
      Atom%NNN = Atom%nMulti*Atom%nJmax

      !
      ! Allocations
      !

      ! Orbital angular momentum
      allocate(Atom%rLval(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rLval)
      ! Spin angular momentum
      allocate(Atom%Sval(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%Sval)
      ! Degeneration (2L+1)*(2S+1)
      allocate(deg(Atom%nMulti))
      ! Number of FS levels per term
      allocate(Atom%nJ(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%nJ)
      ! Ionization stage
      allocate(Atom%stage(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%stage)
      ! Angular momentum
      allocate(Atom%rJval(Atom%nJmax,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%rJval)
      Atom%rJval = -1d0
      ! Frequency of FS levels
      allocate(Atom%FSfreq(Atom%nJmax,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%FSfreq)
      ! Frequency of term
      allocate(Atom%TRfreq(Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%TRfreq)

      ! (g-1)/2, with g statistical weight of n level
      Atom%rLval = (/ .5d0, 3.5d0, 8.5d0, 15.5d0, 24.5d0, 0d0 /)
      ! The spin is artificially put to 0, this atom is given in n
      ! levels
      Atom%Sval = 0d0
      ! Ionization state
      Atom%stage = 1
      Atom%stage(Atom%nMulti) = 2
      ! Degeneration of each level
      deg = (/ 2d0, 8d0, 18d0, 32d0, 50d0, 1d0 /)
      ! Number of sublevels for each level
      Atom%nJ = 1
      ! J = L because S = 0
      Atom%rJval(1,:) = Atom%rLval
      ! Energies
      Atom%TRfreq = (/ 0d0, .82258211d0, .97491219d0, &
                       1.02822766d0, 1.05290508d0, &
                       1.09677617d0 /)
      ! Multilevel
      Atom%ML = .True.
      ! nJ = 1 ==> FSfreq = TRfreq
      Atom%FSfreq(1,:) = Atom%TRfreq

      ! Determine the number of FS levels in this atom
      Atom%nlevel = sum(Atom%nJ)

      ! Term given the level
      allocate(Atom%term(Atom%nlevel))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%term)
      ! Sublevel in its term given the global level
      allocate(Atom%sublevel(Atom%nlevel))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%sublevel)
      ! Degeneration of each term
      allocate(Atom%deg(Atom%nlevel))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%deg)
      ! Allocate indexing by term
      allocate(Atom%irho(Atom%nMulti))

      ! Save degeneracy
      Atom%deg = deg

      ! Initialize
      ii = 0

      ! For each term
      do iterm=1,Atom%nMulti

        ! Memory count
        MRAMc = MRAMc + 1d-6*sizeof(Atom%irho(iterm))

        ! Allocate rho indexing
        allocate(Atom%irho(iterm)%irho_ij(Atom%nJ(iterm)))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%irho(iterm)%irho_ij)

        ! For each level
        do iJ=1,Atom%nJ(iterm)

          ! Advance indexing
          ii = ii + 1

          ! Store indexing
          Atom%irho(iterm)%irho_ij(iJ) = ii
          Atom%term(ii) = iterm
          Atom%sublevel(ii) = iJ

        end do ! Levels
      end do ! Terms


      !!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Radiative transitions
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!

      !
      ! Allocations
      !

      allocate(Atom%fst(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst)
      ! Energy of the term-term transition
      allocate(Atom%Dfreq(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%Dfreq)
      ! Type of elastic broadening
      allocate(Atom%broad_type(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%broad_type)
      ! Arguments to calculate the broadening
      allocate(Atom%broad_args(4,Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%broad_args)
      ! Argument for Stark broadening
      allocate(Atom%broad_stark(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%broad_stark)
      ! Doppler widths that the line holds
      allocate(Atom%dwvl(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%dwvl)
      ! Doppler widths that the core of the line holds
      allocate(Atom%dwvlc(Atom%ntran))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%dwvlc)
      ! Indexing of radiative transitions
      allocate(Atom%irad(Atom%nMulti,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%irad)
      Atom%irad = 0
      ! Einstein coefficients Aul and Blu
      allocate(Atom%Ecoeff(Atom%nMulti,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%Ecoeff)
      Atom%Ecoeff = 0d0

      !
      ! Hardcoded transitions
      !

      ! 1
      itran = 1
      iterm1 = 2
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = 4.6986d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 600
      Atom%Dwvlc(itran) = 15

      ! 2
      itran = 2
      iterm1 = 3
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .55751d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 10

      ! 3
      itran = 3
      iterm1 = 4
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .12785d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 100
      Atom%Dwvlc(itran) = 3

      ! 4
      itran = 4
      iterm1 = 5
      iterm = 1
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .04125d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 100
      Atom%Dwvlc(itran) = 3

      ! 5
      itran = 5
      iterm1 = 3
      iterm = 2
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .44101d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 3

      ! 6
      itran = 6
      iterm1 = 4
      iterm = 2
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .084193d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 3

      ! 7
      itran = 7
      iterm1 = 5
      iterm = 2
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .025304d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 250
      Atom%Dwvlc(itran) = 3

      ! 8
      itran = 8
      iterm1 = 4
      iterm = 3
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .089860d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 30
      Atom%Dwvlc(itran) = 2

      ! 9
      itran = 9
      iterm1 = 5
      iterm = 3
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .022008d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 30
      Atom%Dwvlc(itran) = 2

      ! 10
      itran = 10
      iterm1 = 5
      iterm = 4
      Atom%irad(iterm1,iterm) = itran
      Atom%irad(iterm,iterm1) = itran
      Atom%fst(itran)%iterml = iterm
      Atom%fst(itran)%itermu = iterm1
      allocate(Atom%fst(itran)%irad(1,1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%irad)
      Atom%fst(itran)%irad(1,1) = 1
      Atom%fst(itran)%nt = 1
      allocate(Atom%fst(itran)%ilevell(1))
      allocate(Atom%fst(itran)%ilevelu(1))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevell)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%fst(itran)%ilevelu)
      Atom%fst(itran)%ilevell = 1
      Atom%fst(itran)%ilevelu = 1
      Atom%Dfreq(itran) = abs(Atom%TRfreq(iterm1) - &
                              Atom%TRfreq(iterm))
      radAcoeff = .026993d0
      Atom%Ecoeff(iterm1,iterm) = radAcoeff
      Atom%Ecoeff(iterm,iterm1) = radAcoeff* &
                                  deg(iterm1)/deg(iterm)/ &
                                  (ConvF*Atom%Dfreq(itran)* &
                                   1d21*(2d0*c)* &
                                   Atom%Dfreq(itran)**2d0)
      Atom%broad_stark(itran) = 1d0
      Atom%broad_type(itran) = 1
      Atom%broad_args(:,itran) = (/ 1d0, 0d0, 1d0, 0d0 /)
      Atom%Dwvl(itran) = 30
      Atom%Dwvlc(itran) = 1


      !!!!!!!!!!!!!!!!!!!!
      !
      ! Photoionizations
      !
      !!!!!!!!!!!!!!!!!!!!

      !
      ! Allocations
      !

      ! Indexing of b-f transitions
      allocate(Atom%iphot(Atom%nMulti,Atom%nMulti))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%iphot)
      Atom%iphot = 0
      ! Photoionization information
      allocate(Atom%phot(Atom%nphot))

      !
      ! Hard-coded transitions
      !

      ! 1
      ii = 1
      MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii))
      iterm1 = 6
      iterm = 1
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      MRAMc = MRAMc + 1d-6*(sizeof(Atom%phot(ii)%infreq) + &
                            sizeof(Atom%phot(ii)%inalpha))
                            
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/22.794d0
      Atom%phot(ii)%inalpha(1) = 6.152d-22
      Atom%phot(ii)%ilevell = iterm

      ! 2
      ii = 2
      MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii))
      iterm1 = 6
      iterm = 2
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      MRAMc = MRAMc + 1d-6*(sizeof(Atom%phot(ii)%infreq) + &
                            sizeof(Atom%phot(ii)%inalpha))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/91.176d0
      Atom%phot(ii)%inalpha(1) = 1.379d-21
      Atom%phot(ii)%ilevell = iterm

      ! 3
      ii = 3
      MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii))
      iterm1 = 6
      iterm = 3
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      MRAMc = MRAMc + 1d-6*(sizeof(Atom%phot(ii)%infreq) + &
                            sizeof(Atom%phot(ii)%inalpha))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/205.147d0
      Atom%phot(ii)%inalpha(1) = 2.149d-21
      Atom%phot(ii)%ilevell = iterm

      ! 4
      ii = 4
      MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii))
      iterm1 = 6
      iterm = 4
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      MRAMc = MRAMc + 1d-6*(sizeof(Atom%phot(ii)%infreq) + &
                            sizeof(Atom%phot(ii)%inalpha))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/364.705d0
      Atom%phot(ii)%inalpha(1) = 2.923d-21
      Atom%phot(ii)%ilevell = iterm

      ! 5
      ii = 5
      MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(ii))
      iterm1 = 6
      iterm = 5
      Atom%iphot(iterm1,iterm) = ii
      Atom%iphot(iterm,iterm1) = ii
      Atom%phot(ii)%edge = Atom%TRfreq(iterm1) - Atom%TRfreq(iterm)
      allocate(Atom%phot(ii)%infreq(1))
      allocate(Atom%phot(ii)%inalpha(1))
      MRAMc = MRAMc + 1d-6*(sizeof(Atom%phot(ii)%infreq) + &
                            sizeof(Atom%phot(ii)%inalpha))
      Atom%phot(ii)%mode = 1
      Atom%phot(ii)%nfreq = 1
      Atom%phot(ii)%infreq(1) = 1d2/569.852d0
      Atom%phot(ii)%inalpha(1) = 3.699d-21
      Atom%phot(ii)%ilevell = iterm


      !!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Inelastic collisions
      !
      !!!!!!!!!!!!!!!!!!!!!!!!

      ! Number of collisional transitions
      Atom%ncol = 10

      ! Allocate inelastic database
      allocate(Atom%inelas(Atom%ncol))

      ! Temperature data points
      nTmp = 6

      ! Allocate Tbox
      allocate(Atom%Tbox)
      p_T => Atom%Tbox
      nullify(p_T%next)

      ! Put data in the box
      p_T%nTmp = nTmp
      p_T%ind = 1
      p_T%col_type = 0
      p_T%nion = 0
      p_T%flin = .False.

      ! Allocate and set temperature
      allocate(p_T%temp(nTmp))
      p_T%temp = (/ 3d3, 5d3, 7d3, 1d4, 2d4, 3d4 /)
      MRAMc = MRAMc + 1d-6*dble(4*5 + 8*nTmp)

      ! 1
      ii = 1
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 2
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 1.3351d-14,1.0780d-14,9.4856d-15,8.4125d-15, &
               7.0994d-15,6.7550d-15 /)

      ! 2
      ii = 2
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 3
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 8.7453d-16,7.1253d-16,6.3196d-16,5.6633d-16, &
               4.8995d-16,4.7362d-16 /)

      ! 3
      ii = 3
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 4
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 1.6240d-16,1.3263d-16,1.1792d-16,1.0600d-16, &
               9.2277d-17,8.9644d-17 /)

      ! 4
      ii = 4
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 1
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 4.7192d-17,3.8580d-17,3.4337d-17,3.0892d-17, &
               2.6995d-17,2.6265d-17 /)

      ! 5
      ii = 5
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 2
      Atom%inelas(ii)%up = 3
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 2.7435d-13,2.5384d-13,2.4973d-13,2.5293d-13, &
               2.7775d-13,2.9945d-13 /)

      ! 6
      ii = 6
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 2
      Atom%inelas(ii)%up = 4
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 1.8623d-14,1.7872d-14,1.8024d-14,1.8705d-14, &
               2.1454d-14,2.3746d-14 /)

      ! 7
      ii = 7
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 2
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 3.5405d-15,3.4405d-15,3.4966d-15,3.6592d-15, &
               4.2698d-15,4.7832d-15 /)

      ! 8
      ii = 8
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 3
      Atom%inelas(ii)%up = 4
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 9.5940d-13,1.0457d-12,1.1455d-12,1.2881d-12, &
               1.6451d-12,1.8677d-12 /)

      ! 9
      ii = 9
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 3
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 6.15d-14,6.8731d-14,7.6113d-14,8.64d-14, &
               1.1348d-13,1.3281d-13 /)

      ! 10
      ii = 10
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii))
      Atom%inelas(ii)%ind = 1
      Atom%inelas(ii)%col_type = 0
      Atom%inelas(ii)%low = 4
      Atom%inelas(ii)%up = 5
      allocate(Atom%inelas(ii)%Cul(nTmp))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%inelas(ii)%Cul)
      Atom%inelas(ii)%Cul = &
            (/ 2.7090d-12,3.3113d-12,3.8548d-12,4.5498d-12, &
               6.1112d-12,6.9947d-12 /)

      ! Calculate the atomic part of the Doppler width
      Atom%cDopp = dopp/sqrt(Atom%rmass)

      ! Nullify local pointer
      if (associated(p_T)) nullify(p_T)

      ! Control that everything went fine
      call control

      return

      end subroutine AtomH

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalize the abundances of several atoms of the same
      !! species\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine abund(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      logical, dimension(NA):: renorm, checked

      integer:: i1,i2,nat,natr
      integer, dimension(NA):: list

      double precision:: abmod


      ! Control logical variables
      checked = .False.
      renorm = .False.

      ! Collect the flags into a vector for each atom
      do i1=1,NA
        renorm(i1) = Atom(i1)%anorm
      end do

      ! If any of them ask for normalization, normalize
      if (any(renorm)) then

        ! For each atom
        do i1=1,NA

          ! If we already checked this atom, skip
          if (checked(i1)) cycle

          ! Initialize auxiliar variables
          nat = 1
          natr = 0
          list = 0
          list(1) = i1
          abmod = Atom(i1)%abun_mod

          ! If renormalizing this atom, add to count
          if (renorm(i1)) natr = natr + 1

          ! For all other atoms
          do i2=1,NA

            ! If the atom is checked or if it is the one already
            ! selected, skip
            if (checked(i2).or.i1.eq.i2) cycle

            ! If they are the same element, add to the counter
            if (Atom(i1)%Element.eq.Atom(i2)%Element) then

              ! Add to the list of checked atoms
              nat = nat + 1
              list(nat) = i2

              ! Add this atom to the ones to modify
              if (renorm(i2)) natr = natr + 1

              ! Accumulate the modification factors
              abmod = abmod + Atom(i2)%abun_mod

            end if ! Same element

          end do ! All other atoms

          ! Set all checket atoms to checked
          do i2=1,nat
            checked(list(i2)) = .True.
          end do

          ! If there was an atom asking for renorm
          if (natr.gt.0) then

            ! Master
            if (pid.eq.0) then

              ! Verbose
              umsg = ' - Atom '//trim(Atom(i1)%Element)// &
                     ' abundance modifier normalized'
              call verbose

            end if ! Master

            ! Normalize the atoms on the list
            do i2=1,nat

              ! Recalculate abundance
              Atom(list(i2))%abun_mod = Atom(list(i2))%abun_mod/abmod

            end do ! Normalize atoms on the list

          end if ! There was an atom asking to renormalize

        end do ! Atoms

      end if ! At least one atom wants to renormalize


      end subroutine abund

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the fine-structure Einstein coefficients for the
      !! given atom\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      subroutine setFScoeff(Atom,Flgsg)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg

      ! Local

      integer:: iterm,iterm1,iJ,iJ1,ilevel,ilevel1,itran,ftran

      double precision:: W6,S,rL,rL1,rJ,rJ1,Aul,AJul,BJlu,Dfreq


      ! If the atom is multi-level, ignore quantum numbers
      if (Atom%ML) then

        ! For each lower term/level
        do ilevel=1,Atom%nMulti-1

          ! For each upper term/level
          do ilevel1=ilevel+1,Atom%nMulti

            ! Check for a transition
            itran = Atom%irad(ilevel,ilevel1)

            ! If no transition, stkip
            if (itran.lt.1) cycle

            ! Get quantum numbers
            rJ = Atom%rJval(1,ilevel)
            rJ1 = Atom%rJval(1,ilevel1)

            ! Get transition frequency
            Dfreq = Atom%FSfreq(1,ilevel1) - &
                    Atom%FSfreq(1,ilevel)

            ! Get Aul
            AJul = Atom%Ecoeff(ilevel1,ilevel)

            ! Calculate the corresponding Blu
            BJlu = AJul*(2d0*rJ1+1d0)/(2d0*rJ+1d0)/ &
                   (ConvF*1d21*(2d0*c)*Dfreq**3d0)

            ! Save in structure
            Atom%fst(itran)%Aul(1,1) = AJul
            Atom%fst(itran)%Blu(1,1) = BJlu

          end do ! Upper level
        end do ! Lower level

      ! If atom is multi-term, do the branching
      else

        ! For each lower term
        do iterm=1,Atom%nMulti-1

          ! For each upper term
          do iterm1=iterm+1,Atom%nMulti

            ! Check for a transition
            itran = Atom%irad(iterm,iterm1)

            ! Skip if no transition
            if (itran.lt.1) cycle

            ! Get quantum numbers
            rL = Atom%rLval(iterm)
            rL1 = Atom%rLval(iterm1)
            S = Atom%Sval(iterm)

            ! Get Einstein coefficient
            Aul = Atom%Ecoeff(iterm1,iterm)

            ! For each level in the lower term
            do iJ=1,Atom%nJ(iterm)

              ! For each level in the upper term
              do iJ1=1,Atom%nJ(iterm1)

                ! Check for a transition
                ftran = Atom%fst(itran)%irad(iJ1,iJ)

                ! Skip if no fine-structure transition
                if (ftran.lt.1) cycle

                ! Get total angular momenta
                rJ = Atom%rJval(iJ,iterm)
                rJ1 = Atom%rJval(iJ1,iterm1)

                ! Get frequency
                Dfreq = Atom%FSfreq(iJ1,iterm1) - &
                        Atom%FSfreq(iJ,iterm)

                ! Branch the Aul
                W6 = fun6j(rL1,rL,1d0,rJ,rJ1,S,Flgsg)
                AJul = (2d0*rL1+1d0)*(2d0*rJ+1d0)*W6*W6*Aul

                ! Calculate the corresponding FS-Blu
                BJlu = AJul*(2d0*rJ1+1d0)/(2d0*rJ+1d0)/ &
                       (ConvF*1d21*(2d0*c)*Dfreq**3d0)

                ! Store Einstein coefficients
                Atom%fst(itran)%Aul(iJ1,iJ) = AJul
                Atom%fst(itran)%Blu(iJ,iJ1) = BJlu

              end do ! Upper levels
            end do ! Lower levels
          end do ! Upper terms
        end do ! Lower terms

      end if ! Multi-level or Multi-term

      return

      end subroutine setFScoeff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the known quantities for a given LTE transition\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!  line(LTEline_class): Structure with LTE line data
      subroutine setup_LTE_transition(Atom,line)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(LTEline_class), intent(inout):: line

      ! Local

      integer:: iMu,iMf

      double precision:: rMu,rMf


      ! If the line is from an atom with an atomic model
      if (line%is_passive) then

        ! Get data from model
        line%rmass = Atom(line%ia)%rmass
        line%abund = Atom(line%ia)%abun
        line%cDopp = Atom(line%ia)%cDopp

      ! If the line is from an atom without an atomic model
      else

        ! Get data from tabulation
        line%rmass = recallmass_ind(line%ele)
        line%abund = recallabund_ind(line%ele)
        line%cDopp = dopp/sqrt(line%rmass)

      end if ! The LTE is from an atom with an atomic model

      ! Transition freq. and Blu Einstein coefficient
      line%Dfreq = line%Eu - line%El
      line%Blu = line%Aul*(2d0*line%Ju + 1d0)/ &
                          (2d0*line%Jl + 1d0)/ &
                 (ConvF*line%Dfreq*1d21*(2d0*c)* &
                  line%Dfreq**2d0)

      ! Number of M levels indexes for each level
      line%nMu = nint(2d0*line%Ju+1d0)
      line%nMl = nint(2d0*line%Jl+1d0)

      ! Initialize number magnetic components
      line%ncom = 0

      ! Run over Mu
      do iMu=1,line%nMu

        ! Magnetic number
        rMu = -line%Ju + dble(iMu-1)

        ! Run over Ml
        do iMf=1,line%nMl

          ! Magnetic number
          rMf = -line%Jl + dble(iMf-1)

          ! Selection rules
          if (nint(abs(rMu-rMf)).gt.1) cycle

          ! Add to count
          line%ncom = line%ncom + 1

        end do ! Mf
      end do ! Mu

      ! If valid number of frequencies
      if (line%nfreq.gt.0) then

        ! Ensure odd number
        if (mod(line%nfreq,2).eq.0) line%nfreq = line%nfreq + 1
        if (mod(line%nfreqc,2).eq.0) line%nfreqc = line%nfreqc + 1

      end if ! Valid number of frequencies

      ! Add frequencies to global count if not skipping
      if (.not.line%nowave) nfreq = nfreq + line%nfreq

      end subroutine setup_LTE_transition

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if an LTE transition is part of an atomic model and
      !! remove it from the model\n
      !!     Atom(Atom_class): Structure with atomic data\n
      !!  line(LTEline_class): Structure with LTE line data
      subroutine remove_LTE_transition(Atom,line)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(LTEline_class), intent(in):: line

      ! Local

      logical:: no_found_l, no_found_u

      integer:: iterm,iJ,Jl2,Ju2,iterml,iJl,iJu,itermu

      double precision:: El,Eu,factl,factu


      !
      ! Find levels
      !

      ! Initialize
      no_found_l = .True.
      no_found_u = .True.

      ! Get data from line
      Jl2 = nint(line%Jl*2d0)
      Ju2 = nint(line%Ju*2d0)
      El = line%El
      Eu = line%Eu

      ! Relative energy factors
      if (abs(El).gt.0d0) then
        factl = 1d0/abs(El)
      else
        factl = 1d0
      end if
      if (abs(Eu).gt.0d0) then
        factu = 1d0/abs(Eu)
      else
        factu = 1d0
      end if

      ! For each term
      do iterm=1,Atom%nMulti

        ! If not even same stage, skip
        if (Atom%stage(iterm).ne.line%stage) cycle

        ! For each level within the term
        do iJ=1,Atom%nJ(iterm)

          ! If not found lower yet
          if (no_found_l) then

            ! Check J
            if (nint(Atom%rJval(iJ,iterm)*2d0).eq.Jl2) then

              ! Compare energies
              if (abs(El - Atom%FSfreq(iJ,iterm))*factl.le. &
                  TINYSP) then

                ! Found
                no_found_l = .False.
                iterml = iterm
                iJl = iJ

              end if ! Same energy
            end if ! Same J
          end if ! Not found l

          ! If not found upper yet
          if (no_found_u) then

            ! Check J
            if (nint(Atom%rJval(iJ,iterm)*2d0).eq.Ju2) then

              ! Compare energies
              if (abs(Eu - Atom%FSfreq(iJ,iterm))*factu.le. &
                  TINYSP) then

                  ! Found
                  no_found_u = .False.
                  itermu = iterm
                  iJu = iJ

              end if ! Same energy
            end if ! Same J
          end if ! Not found u

        end do ! Levels

        ! If found, stop search
        if (.not.no_found_l.and..not.no_found_u) exit

      end do ! Terms

      ! If not found, the line is not here and we can leave
      if (no_found_l.or.no_found_u) return

      ! Multi-level
      if (Atom%ML) then

        ! Just remove the transition (from the calculations)
        Atom%irad(itermu,iterml) = 0
        Atom%irad(iterml,itermu) = 0

      ! Multi-term
      else

        ! If transition is not forbidden by electric dipole rule
        if (abs(Jl2-Ju2).le.1.and.Jl2+Ju2.gt.0) then

          ! Remove between terms
          Atom%irad(itermu,iterml) = 0
          Atom%irad(iterml,itermu) = 0

        end if ! Not forbidden
      end if ! Multi-level/Multi-term

      return

      end subroutine remove_LTE_transition

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate an array of atomic models\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!          nn(integer): Size to allocate
      subroutine allocateatom(Atom,nn)

      ! I/O

      type(Atom_class), dimension(:), allocatable, intent(out):: Atom
      integer, intent(in):: nn

      ! Allocate array
      if (nn.gt.0) allocate(Atom(nn))

      return

      end subroutine allocateatom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Provide each atomic model with a unique label\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !   !       nn(integer): Size to allocate
      subroutine set_atom_label(Atom,nn)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      integer, intent(in):: nn

      ! Local

      character(len=10):: label

      logical:: repeated
      logical, dimension(nn):: check

      integer:: ia,ja,ka,irep


      ! For each atom, create default label
      do ia=1,nn

        ! Initialize label label
        Atom(ia)%file_label = '          '

        ! Write atomic label
        if (Atom(ia)%Element(1:1).eq.' ') then
          Atom(ia)%file_label(1:1) = Atom(ia)%Element(2:2)
        else
          Atom(ia)%file_label = Atom(ia)%Element
        end if

      end do ! Atoms

      ! Initialize
      check = .False.

      ! For each model atom but the last
      do ia=1,nn-1

        ! If already checked, skip
        if (check(ia)) cycle

        ! Current label
        label = Atom(ia)%file_label

        ! Initialize
        repeated = .False.
        irep = 0

        ! For the rest of atoms
        do ja=ia+1,nn

          ! If label is the same
          if (trim(Atom(ja)%file_label).eq.trim(label)) then

            ! Flag repeated
            repeated = .True.
            exit

          end if ! Same label

        end do ! Rest of atoms

        ! If repeated
        if (repeated) then

          ! Change current and initialize index
          Atom(ia)%file_label = trim(label)//'-1'

          ! Initialize index
          ka = 1

          ! For the rest of atoms
          do ja=ia+1,nn

            ! Check if same label
            if (trim(Atom(ja)%file_label).eq.trim(label)) then

              ! Advance
              ka = ka + 1

              ! Get string with index number
              if (ka.lt.10) then
                write(Atom(ja)%file_label,'("-",i1)') ka
              else if (ka.lt.100) then
                write(Atom(ja)%file_label,'("-",i2)') ka
              else if (ka.lt.1000) then
                write(Atom(ja)%file_label,'("-",i3)') ka
              end if

              ! Get new label
              Atom(ja)%file_label = trim(label)// &
                                    trim(Atom(ja)%file_label)

              ! Flag as checked
              check(ja) = .True.

            end if ! Same label

          end do ! Rest of atoms

        end if ! Repeated

        ! Flag current atom as checked
        check(ia) = .True.

      end do ! Atoms

      return

      end subroutine set_atom_label

!#####################################################################
!#####################################################################
!#####################################################################

      !> Shift all atoms one position to the right in the array, as
      !! well as their associated filenames\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!  files(strarr_class(:)): List of file names
      subroutine shiftatoms(Atom,files)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(strarr_class), dimension(:), &
                          allocatable, intent(inout):: files

      ! Local

      type(Atom_class), dimension(:), allocatable:: Atomaux
      type(strarr_class), dimension(:), allocatable:: filesaux


      ! If there are no background atoms
      if (nAb.lt.1) then

          ! No need to shift
          nAb = 1
          allocate(Atom(nAb))
          allocate(files(nAb))
          files(1)%str = 'N'

      ! If there is data to shift
      else

        ! Copy data into auxiliar
        allocate(Atomaux(nAb))
        allocate(filesaux(nAb))
        Atomaux = Atom
        filesaux = files
        deallocate(Atom)
        deallocate(files)

        ! Increase the number of background atoms
        nAb = nAb + 1

        ! Create a larger array and copy back
        allocate(Atom(nAb))
        allocate(files(nAb))
        Atom(2:nAb) = Atomaux
        files(2:nAb) = filesaux
        files(1)%str = 'N'

        ! Deallocate auxiliar
        deallocate(Atomaux)
        deallocate(filesaux)

      end if ! Previous inputs

      ! Check if everything is fine
      call control

      return

      end subroutine shiftatoms

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate height dependent pupulation variables\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!          nn(integer): Size of the structure
      subroutine prepareatom(Atom,nn)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      integer, intent(in):: nn

      ! Local

      integer:: ia


      ! For each atom
      do ia=1,nn

        ! Allocate total population
        allocate(Atom(ia)%n(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%n)

      end do ! Atoms

      end subroutine prepareatom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate Master MPI variables for each model atom\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine prepareatomMPI(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ia,ii


      ! If Master and MPI
      if (pid.eq.0.and.nproc.gt.1) then

        ! For each atom
        do ia=1,na

          ! Lower frequency limits for Master
          allocate(Atom(ia)%Mif0(Atom(ia)%ntran,0:nproc-1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%Mif0)
          ! Upper frequency limits for Master
          allocate(Atom(ia)%Mif1(Atom(ia)%ntran,0:nproc-1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%Mif1)
          ! Lower frequency weight for Master
          allocate(Atom(ia)%MW0(Atom(ia)%ntran,0:nproc-1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%MW0)
          ! Upper frequency weight for Master
          allocate(Atom(ia)%MW1(Atom(ia)%ntran,0:nproc-1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%MW1)

          ! For each b-f transition
          do ii=1,Atom(ia)%nphot
            ! Allocate limits for master
            allocate(Atom(ia)%phot(ii)%Mif0(0:nproc-1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%phot(ii)%Mif0)
            allocate(Atom(ia)%phot(ii)%Mif1(0:nproc-1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%phot(ii)%Mif1)
            ! Allocate weights for master
            allocate(Atom(ia)%phot(ii)%MW0(0:nproc-1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%phot(ii)%MW0)
            allocate(Atom(ia)%phot(ii)%MW1(0:nproc-1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%phot(ii)%MW1)
          end do ! b-f transitions
        end do ! Atoms

      end if ! Master of MPI

      end subroutine prepareatomMPI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize the size of the M blocks assuming there is a
      !! magnetic field\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine initMblock(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ia,nM,iM,iJ,i,iterm

      double precision:: rL,S,rJ,rJmin,rJmax,pS,rM,rJm


      ! For each atom
      do ia=1,nA

        ! For each term
        do iterm=1,Atom(ia)%nMulti

          ! Get term quantities
          rL = Atom(ia)%rLval(iterm)
          S = Atom(ia)%Sval(iterm)

          ! Multilevel
          if (Atom(ia)%ML) then

            ! Get level quantities
            rJ = Atom(ia)%rJval(1,iterm)
            nM = nint(2d0*rJ + 1d0)

            ! For each M
            do iM=1,nM

              ! Size of the block of this M (just one)
              Atom(ia)%nblk(iM,iterm) = 1

            end do ! Magnetic sublevels

            ! Rest of empty space
            do iM=nM+1,Atom(ia)%nMmax

              ! Size of the block of this M (just one)
              Atom(ia)%nblk(iM,iterm) = 0

            end do ! Magnetic sublevels

          ! Multiterm
          else

            ! Minimum and maximum total angular momenta
            rJmin = abs(rL - S)
            rJmax = rL + S

            ! Spin factor
            pS = S*(S+1d0)*(2d0*S + 1d0)

            ! Number of magnetic sublevels
            nM = nint(2d0*rJmax + 1d0)

            ! Run over the magnetic block
            do iM=1,nM

              ! Magnetic quantum number
              rM = -rJmax + dble(iM-1)

              ! Minimun value of total angular momentum
              rJm = max(abs(rM),rJmin)

              ! initialize the column index
              i=0

              ! Column loop
              do iJ=1,Atom(ia)%nJ(iterm)

                ! Total angular momentum
                rJ = Atom(ia)%rJval(iJ,iterm)

                ! If total angular momentum is above minimum
                if (rJ.ge.rJm) then

                  ! Increase the column index
                  i = i + 1

                end if ! Test J compatibility with M

              end do ! iJ

              ! Block dimension
              Atom(ia)%nblk(iM,iterm) = i

            end do ! iM

            ! Rest of empty space
            do iM=nM+1,Atom(ia)%nMmax

              ! Size of the block of this M (just one)
              Atom(ia)%nblk(iM,iterm) = 0

            end do ! Rest of empty space

          end if ! Multilevel or multiterm

        end do ! Terms
      end do ! Atoms

      end subroutine initMblock

!#####################################################################
!#####################################################################
!#####################################################################

      !> Setup all necessary atomic indexings for transitions
      !! components\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!   Input(Input_class): Structure with configuration data\n
      !!          li(logical): If doing intensity only calculations\n
      !!          lp(logical): If doing polarization calculations\n
      !!      nfield(logical): If expecting points without magnetic
      !!                       field\n
      !!      yfield(logical): If expecting points with magnetic
      !!                       field
      subroutine set_atom_indexes(Atom,Input,li,lp,nfield,yfield)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Input_class), intent(in):: Input
      logical, intent(in):: li,lp,nfield,yfield

      ! Local

      logical:: nfield2

      integer:: ia,iterm,iM,jtran,iti,itran,minto,maxto
      integer:: fjtran,ffjtran,ffktran,fitran,ffitran
      integer:: iMf,mF,iMu,iU,iMu1,iU1,iMl,iL,iMl1,iL1
      integer:: iJu,iJu1,iJf,iJl,iJl1,itermu,itermf,iterml
      integer:: nMm,nMu,nMf,nMl,nblk,ii,i1,i2,nti
      integer:: iq,ip,iq1,ip1,iQQ,iPP
      integer:: MiindU,MiindF,MiindU1,MiindL,MiindL1,nnchlt
      integer:: MindU,MindF,MindU1,MindL,MindL1
      integer:: indF,indL,indL1,indU,indU1

      double precision:: rJumax,rJfmax,rJlmax,rMf,rMu,rMl,rMu1,rMl1
      double precision:: rJf,rJu,rJl,rJu1,rJl1,q,p,q1,p1,QQ,PP


      ! Correct no field accounting for two-step
      nfield2 = nfield.or.Input%two_step_pol

      ! For each atom
      do ia=1,nA

        ! If we are doing intensity
        if (li) then

          ! If PRD
          if (PRD) then

            !
            ! Count and index output transition
            !

            ! Initialize counters
            minto = Atom(ia)%nftran + 1
            maxto = 0
            nti = 0

            ! For each output transition (term)
            do jtran=1,Atom(ia)%ntran

              ! If not PRD, skip
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! For each output transition (FS)
              do fjtran=1,Atom(ia)%fst(jtran)%nt

                ! Get rolling index
                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

                ! Update limits and count
                if (ffjtran.lt.minto) minto = ffjtran
                if (ffjtran.gt.maxto) maxto = ffjtran
                nti = nti + 1

              end do ! FS transitions
            end do ! Output transition (term)

            ! Store number of output transition
            Atom(ia)%ntrano = nti

            ! Allocate output transition structure
            allocate(Atom(ia)%tranoI(nti))

            ! Allocate and initialize output PRD transition indexing
            allocate(Atom(ia)%itrano(minto:maxto))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%itrano)
            Atom(ia)%itrano = -1

            ! Initialize counter
            nti = 0

            ! For each output transition (term)
            do jtran=1,Atom(ia)%ntran

              ! If not PRD, skip
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! For each output transition (FS)
              do fjtran=1,Atom(ia)%fst(jtran)%nt

                ! Get rolling index
                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

                ! Advance and store
                nti = nti + 1
                Atom(ia)%itrano(ffjtran) = nti

              end do ! FS transitions
            end do ! Output transition (term)

            !
            ! Count and index input transitions
            !

            ! For each output transition (term)
            do jtran=1,Atom(ia)%ntran

              ! Get terms
              itermu = Atom(ia)%fst(jtran)%itermu
              itermf = Atom(ia)%fst(jtran)%iterml

              ! If not PRD, skip
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! For each output transition (FS)
              do fjtran=1,Atom(ia)%fst(jtran)%nt

                ! Get rolling index and level indexes
                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
                iJf = Atom(ia)%fst(jtran)%ilevell(fjtran)

                ! Get trano rolling index
                ffktran = Atom(ia)%itrano(ffjtran)

                ! Initialize counter
                nti = 0

                ! For each other lower term
                do iterml=1,itermu-1

                  ! Get transition
                  itran = Atom(ia)%irad(itermu,iterml)

                  ! Skip no transition
                  if (itran.le.0) cycle

                  ! Skip Raman if indicated
                  if (.not.Input%Raman.and.itran.ne.jtran) cycle

                  ! For every level
                  do iJl=1,Atom(ia)%nJ(iterml)

                    ! Get index
                    fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

                    ! Skip no transition
                    if (fitran.le.0) cycle

                    ! Advance index
                    nti = nti + 1

                  end do ! Lower level
                end do ! Lower term

                ! Memory count
                MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%tranoI(ffktran))

                ! Store number of transition
                Atom(ia)%tranoI(ffktran)%nt = nti

                ! Allocate and initialize indexing
                allocate(Atom(ia)%tranoI(ffktran)%indT(nti))
                MRAMc = MRAMc + &
                        1d-6*sizeof(Atom(ia)%tranoI(ffktran)%indT)
                Atom(ia)%tranoI(ffktran)%indT = -1

                ! Initialize counter
                nti = 0

                ! For each other lower term
                do iterml=1,itermu-1

                  ! Get transition
                  itran = Atom(ia)%irad(itermu,iterml)

                  ! Skip no transition
                  if (itran.le.0) cycle

                  ! Skip Raman if indicated
                  if (.not.Input%Raman.and.itran.ne.jtran) cycle

                  ! For every level
                  do iJl=1,Atom(ia)%nJ(iterml)

                    ! Get indexes
                    fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

                    ! Skip no transition
                    if (fitran.le.0) cycle

                    ! Get rolling index
                    ffitran = Atom(ia)%ifst_ij(fitran,itran)

                    ! Advance index and store
                    nti = nti + 1
                    Atom(ia)%tranoI(ffktran)%indT(nti) = ffitran

                  end do ! Lower level
                end do ! Lower term
              end do ! FS transitions
            end do ! Output transition (term)

          end if ! PRD
        end if ! Doing intensity

        !
        ! If we are doing polarization
        !
        if (lp) then

          !
          ! Common transition indexing
          !

          ! Allocate output transition structure
          allocate(Atom(ia)%trano(Atom(ia)%ntran))

          !
          ! Initialize component count
          !

          ! For each transition
          do jtran=1,Atom(ia)%ntran

            ! Initialize count
            Atom(ia)%trano(jtran)%ncomNB = 0
            Atom(ia)%trano(jtran)%ncomB = 0

          end do ! Transitions

          ! If PRD
          if (PRD) then

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Get terms
              itermf = Atom(ia)%fst(jtran)%iterml
              itermu = Atom(ia)%fst(jtran)%itermu

              ! Skip CRD
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! Count input transitions
              nti = 0

              ! For each other lower term
              do iterml=1,itermu-1

                ! Get transition index
                itran = Atom(ia)%irad(itermu,iterml)

                ! No transition, skip
                if (itran.le.0) cycle

                ! No Raman, skip
                if (.not.Input%Raman.and.itran.ne.jtran) cycle

                ! Add count
                nti = nti + 1

              end do ! Other lower term

              ! Memory count
              MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%trano(jtran))

              ! Save size
              Atom(ia)%trano(jtran)%nt = nti

              ! Allocate PRD structure and input transition
              ! indexing
              allocate(Atom(ia)%trano(jtran)%trani(nti))
              allocate(Atom(ia)%trano(jtran)%indT(nti))
              MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%trano(jtran)%indT)

              ! Initialize index
              Atom(ia)%trano(jtran)%indT = -1

              ! Index input transitions
              ii = 0

              ! For each other lower term
              do iterml=1,itermu-1

                ! Get transition index
                itran = Atom(ia)%irad(itermu,iterml)

                ! No transition, skip
                if (itran.le.0) cycle

                ! Skip Raman if indicated
                if (.not.Input%Raman.and.itran.ne.jtran) cycle

                ! Advance index
                ii = ii + 1

                ! Store index
                Atom(ia)%trano(jtran)%indT(ii) = itran

              end do ! Other lower term
            end do ! For each transition

          end if ! If PRD


          !
          ! If we expect points without magnetic fields
          !
          if (nfield2) then

            !
            ! Transiton indexing
            !

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Initialize count
              Atom(ia)%trano(jtran)%ncomNB = 0

              ! Get terms
              itermf = Atom(ia)%fst(jtran)%iterml
              itermu = Atom(ia)%fst(jtran)%itermu

              !
              ! First order indexing
              !

              ! Two rounds
              do i1=1,2

                ! First round
                if (i1.eq.1) then

                  ! Initialize
                  MiindU = 1000000
                  MiindF = 1000000
                  MindU = 0
                  MindF = 0

                ! Second round
                else

                  ! Allocate and initialize
                  allocate(Atom(ia)%trano(jtran)% &
                                    indNB(MiindF:MindF, &
                                          MiindU:MindU))
                  MRAMc = MRAMc + &
                          1d-6*sizeof(Atom(ia)%trano(jtran)%indNB)
                  Atom(ia)%trano(jtran)%indNB = 0

                end if

                ! Reset indexing
                ii = 0

                ! For each Jf
                do iJf=1,Atom(ia)%nJ(itermf)

                  ! Get indexes
                  indF = Atom(ia)%irho(itermf)%irho_ij(iJf)

                  ! Get Jf
                  rJf = Atom(ia)%rJval(iJf,itermf)

                  ! For each Ju
                  do iJu=1,Atom(ia)%nJ(itermu)

                    ! Get indexes
                    indU = Atom(ia)%irho(itermu)%irho_ij(iJu)

                    ! Get Ju
                    rJu = Atom(ia)%rJval(iJu,itermu)

                    ! Electric dipole
                    if (abs(rJu-rJf).gt.1d0.or. &
                        rJu+rJf.lt..25) cycle

                    ! First round
                    if (i1.eq.1) then

                      ! Update limits
                      if (indU.lt.MiindU) MiindU = indU
                      if (indF.lt.MiindF) MiindF = indF
                      if (indU.gt.MindU) MindU = indU
                      if (indF.gt.MindF) MindF = indF

                    ! Second round
                    else

                      ! Advance index and store
                      ii = ii + 1
                      Atom(ia)%trano(jtran)%ncomNB = ii
                      Atom(ia)%trano(jtran)%indNB(indF,indU) = ii

                    end if ! Round

                  end do ! Upper level in term Ju
                end do ! Final level in term Jf
              end do ! Round i1

              ! If PRD and storing
              if (PRD.and.PRAM.and.Atom(ia)%lemiss2(jtran)) then

                ! For each other transition
                do iti=1,Atom(ia)%trano(jtran)%nt

                  ! Get transition and term indexes
                  itran = Atom(ia)%trano(jtran)%indT(iti)
                  iterml = Atom(ia)%fst(itran)%iterml

                  !
                  ! Normal indexing
                  !

                  ! Run twice
                  do i1=1,2

                    ! First round
                    if (i1.eq.1) then

                      ! Initialize
                      MiindU = 1000000
                      MiindU1 = 1000000
                      MiindL = 1000000
                      MiindL1 = 1000000
                      MiindF = 1000000
                      MindU = 0
                      MindU1 = 0
                      MindL = 0
                      MindL1 = 0
                      MindF = 0

                    ! Second round
                    else

                      ! Memory count
                      MRAMc = MRAMc + &
                              1d-6*sizeof(Atom(ia)%trano(jtran)% &
                                                   trani(iti))

                      ! Allocate and initialize
                      allocate(Atom(ia)%trano(jtran)%trani(iti)% &
                               indNB(MiindL1:MindL1,MiindL:MindL, &
                                     MiindF:MindF,MiindU1:MindU1, &
                                     MiindU:MindU))
                      MRAMc = MRAMc + &
                              1d-6*sizeof(Atom(ia)%trano(jtran)% &
                                                   trani(iti)%indNB)
                      Atom(ia)%trano(jtran)%trani(iti)%indNB = 0

                    end if ! Round

                    ! Reset indexing
                    ii = 0

                    ! For each Jf
                    do iJf=1,Atom(ia)%nJ(itermf)

                      ! Get indexes
                      indF = Atom(ia)%irho(itermf)%irho_ij(iJf)

                      ! Get Jf
                      rJf = Atom(ia)%rJval(iJf,itermf)

                      ! For each Ju
                      do iJu=1,Atom(ia)%nJ(itermu)

                        ! Get indexes
                        indU = Atom(ia)%irho(itermu)%irho_ij(iJu)

                        ! Get Ju
                        rJu = Atom(ia)%rJval(iJu,itermu)

                        ! Electric dipole
                        if (abs(rJu-rJf).gt.1d0.or. &
                            rJu+rJf.lt..25) cycle

                        ! For each Ju'
                        do iJu1=1,Atom(ia)%nJ(itermu)

                          ! Get indexes
                          indU1 = Atom(ia)%irho(itermu)%irho_ij(iJu1)

                          ! Get Ju'
                          rJu1 = Atom(ia)%rJval(iJu1,itermu)

                          ! Electric dipole
                          if (abs(rJu1-rJf).gt.1d0.or. &
                              rJu1+rJf.lt..25) cycle

                          ! For each Jl
                          do iJl=1,Atom(ia)%nJ(iterml)

                            ! Get indexes
                            indL = Atom(ia)%irho(iterml)%irho_ij(iJl)

                            ! Get Jl
                            rJl = Atom(ia)%rJval(iJl,iterml)

                            ! Electric dipole
                            if (abs(rJu-rJl).gt.1d0.or. &
                                rJu+rJl.lt..25) cycle

                            ! For each Jl'
                            do iJl1=1,Atom(ia)%nJ(iterml)

                              ! Get indexes
                              indL1 = Atom(ia)%irho(iterml)% &
                                               irho_ij(iJl1)

                              ! Get Jl1
                              rJl1 = Atom(ia)%rJval(iJl1,iterml)

                              ! Electric dipole
                              if (abs(rJu1-rJl1).gt.1d0.or. &
                                  rJu1+rJl1.lt..25) cycle

                              ! First round
                              if (i1.eq.1) then

                                ! Update limits
                                if (indU.lt.MiindU) MiindU = indU
                                if (indU1.lt.MiindU1) MiindU1 = indU1
                                if (indL.lt.MiindL) MiindL = indL
                                if (indL1.lt.MiindL1) MiindL1 = indL1
                                if (indF.lt.MiindF) MiindF = indF
                                if (indU.gt.MindU) MindU = indU
                                if (indU1.gt.MindU1) MindU1 = indU1
                                if (indL.gt.MindL) MindL = indL
                                if (indL1.gt.MindL1) MindL1 = indL1
                                if (indF.gt.MindF) MindF = indF

                              ! Second round
                              else

                                ! Advance index and store
                                ii = ii + 1
                                Atom(ia)%trano(jtran)%trani(iti)% &
                                  indNB(indL1,indL,indF, &
                                        indU1,indU) = ii

                              end if ! Round

                            end do ! Jl'
                          end do ! Jl
                        end do ! Ju'
                      end do ! Ju
                    end do ! Jf
                  end do ! Round of counting
                end do ! Input transitions

              end if ! PRD transition

            end do ! Transitions

          end if ! We expect points without magnetic field


          !
          ! If we expect magnetic fields
          !
          if (yfield) then

            !
            ! Term indexing
            !

            ! For each term
            do iterm=1,Atom(ia)%nMulti

              ! Number of M values and maximum block
              nMm = nint(2d0*(Atom(ia)%rLval(iterm) + &
                              Atom(ia)%Sval(iterm)) + 1d0)
              nblk = 0

              ! Find the maximum block size
              do iM=1,nMm
                if (Atom(ia)%nblk(iM,iterm).gt.nblk) &
                  nblk = Atom(ia)%nblk(iM,iterm)
              end do

              ! Allocate index size
              allocate(Atom(ia)%irho(iterm)%jM(nblk,nMm))
              MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%irho(iterm)%jM)

              !
              ! Index the components
              !

              ! Initialize M
              ii = 0

              ! For each M
              do i1=1,nMm

                ! For each mu in the block
                do i2=1,Atom(ia)%nblk(i1,iterm)

                  ! Advance index and store
                  ii = ii + 1
                  Atom(ia)%irho(iterm)%jM(i2,i1) = ii

                end do ! mu in block
              end do ! M values
            end do ! Terms

            !
            ! Transiton indexing
            !

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Get terms
              itermf = Atom(ia)%fst(jtran)%iterml
              itermu = Atom(ia)%fst(jtran)%itermu

              ! Number of magnetic components for each term
              rJumax = Atom(ia)%rLval(itermu) + &
                       Atom(ia)%Sval(itermu)
              nMu = nint(2d0*rJumax + 1d0)
              rJfmax = Atom(ia)%rLval(itermf) + &
                       Atom(ia)%Sval(itermf)
              nMf = nint(2d0*rJfmax + 1d0)

              !
              ! First order indexing
              !

              ! Two rounds
              do i1=1,2

                ! First round
                if (i1.eq.1) then

                  ! Initialize
                  MiindU = 1000000
                  MiindF = 1000000
                  MindU = 0
                  MindF = 0

                ! Second round
                else

                  ! Allocate and initialize
                  allocate(Atom(ia)%trano(jtran)% &
                                    indB(MiindF:MindF, &
                                         MiindU:MindU))
                  MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)% &
                                              trano(jtran)%indB)
                  Atom(ia)%trano(jtran)%indB = 0

                end if ! Round

                ! Reset indexing
                ii = 0

                ! For each Mf
                do iMf=1,nMf

                  ! Value of Mf
                  rMf = -rJfmax + dble(iMf-1)

                  ! For each mu_f
                  do mF=1,Atom(ia)%nblk(iMf,itermf)

                    ! Get jM index
                    indF = Atom(ia)%irho(itermf)%jM(mF,iMf)

                    ! For each Mu
                    do iMu=1,nMu

                      ! Value of Mu
                      rMu = -rJumax + dble(iMu-1)

                      ! Difference between M momentums, done integer
                      q = rMu - rMf
                      iq = nint(q)

                      ! If not pi nor sigma, skip
                      if(abs(iq).gt.1) cycle

                      ! For each mu_u
                      do iU=1,Atom(ia)%nblk(iMu,itermu)

                        ! Get jM index
                        indU = Atom(ia)%irho(itermu)%jM(iU,iMu)

                        ! First round
                        if (i1.eq.1) then

                          ! Update maximum
                          if (indU.lt.MiindU) MiindU = indU
                          if (indF.lt.MiindF) MiindF = indF
                          if (indU.gt.MindU) MindU = indU
                          if (indF.gt.MindF) MindF = indF

                        ! Second round
                        else

                          ! Advance and store index
                          ii = ii + 1
                          Atom(ia)%trano(jtran)%ncomB = ii
                          Atom(ia)%trano(jtran)%indB(indF,indU) = ii

                        end if ! Round

                      end do ! iU
                    end do ! iMu
                  end do ! mF
                end do ! iMf
              end do ! Round i1

              ! If PRD
              if (PRD.and.Atom(ia)%lemiss2(jtran)) then

                ! For each other transition
                do iti=1,Atom(ia)%trano(jtran)%nt

                  ! Get transition and term indexes
                  itran = Atom(ia)%trano(jtran)%indT(iti)
                  iterml = Atom(ia)%fst(itran)%iterml

                  ! Get M size
                  rJlmax = Atom(ia)%rLval(iterml) + &
                           Atom(ia)%Sval(iterml)
                  nMl = nint(2d0*rJlmax + 1d0)

                  !
                  ! Normal indexing
                  !

                  ! Run twice
                  do i1=1,2

                    ! First round
                    if (i1.eq.1) then

                      ! Initialize
                      MiindU = 1000000
                      MiindU1 = 1000000
                      MiindL = 1000000
                      MiindL1 = 1000000
                      MiindF = 1000000
                      MindU = 0
                      MindU1 = 0
                      MindL = 0
                      MindL1 = 0
                      MindF = 0

                    ! Second round
                    else

                      ! Memory count
                      MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)% &
                                                  trano(jtran)% &
                                                  trani(iti))

                      ! Allocate and initialize
                      allocate(Atom(ia)%trano(jtran)%trani(iti)% &
                               indB(MiindL1:MindL1, &
                                    MiindL:MindL, &
                                    MiindF:MindF, &
                                    MiindU1:MindU1, &
                                    MiindU:MindU))
                      MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)% &
                                                  trano(jtran)% &
                                                  trani(iti)%indB)
                      Atom(ia)%trano(jtran)%trani(iti)%indB = 0
                      nnchlt = 0

                    end if ! Round

                    ! Reset indexing
                    ii = 0

                    ! For each Mf
                    do iMf=1,nMf

                      ! Value of Mf
                      rMf = -rJfmax + dble(iMf-1)

                      ! For each mu_f
                      do mF=1,Atom(ia)%nblk(iMf,itermf)

                        ! Get jM index
                        indF = Atom(ia)%irho(itermf)%jM(mF,iMf)

                        ! For each Mu
                        do iMu=1,nMu

                          ! Value of Mu
                          rMu = -rJumax + dble(iMu-1)

                          ! Difference between M momentums, done
                          ! integer
                          q = rMu - rMf
                          iq = nint(q)

                          ! If not pi nor sigma, skip
                          if(abs(iq).gt.1) cycle

                          ! For each mu_u
                          do iU=1,Atom(ia)%nblk(iMu,itermu)

                            ! Get jM index
                            indU = Atom(ia)%irho(itermu)%jM(iU,iMu)

                            ! For each Mu'
                            do iMu1=1,nMu

                              ! Value of Mu'
                              rMu1 = -rJumax + dble(iMu1-1)

                              ! Difference between M momentums
                              q1 = rMu1-rMf
                              QQ = q1-q

                              ! Convert to integers
                              iq1 = nint(q1)
                              iQQ = nint(QQ)

                              ! If not pi or sigma, skip
                              if(abs(iq1).gt.1) cycle

                              ! For each mu_u'
                              do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                                ! Get jM index
                                indU1 = Atom(ia)%irho(itermu)% &
                                                 jM(iU1,iMu1)

      !
      ! Reset indentation
      !

      ! For each Ml
      do iMl=1,nMl

        ! Value of Ml
        rMl = -rJlmax + dble(iMl-1)

        ! Difference between M momentums
        p = rMu-rMl
        ip = nint(p)

        ! If not pi nor sigma, skip
        if(abs(ip).gt.1) cycle

        ! For each mu_l
        do iL=1,Atom(ia)%nblk(iMl,iterml)

          ! Get jM index
          indL = Atom(ia)%irho(iterml)%jM(iL,iMl)

          ! For each Ml'
          do iMl1=1,nMl

            ! Value of Ml'
            rMl1 = -rJlmax + dble(iMl1-1)

            ! Difference between M momentums
            p1 = rMu1-rMl1
            PP = p1-p

            ! Convert to integer
            ip1 = nint(p1)
            iPP = nint(PP)

            ! If not pi nor sigma, skip
            if(abs(ip1).gt.1) cycle

            ! For each mu_l'
            do iL1=1,Atom(ia)%nblk(iMl1,iterml)

              ! Get jM index
              indL1 = Atom(ia)%irho(iterml)%jM(iL1,iMl1)

              ! First round
              if (i1.eq.1) then

                ! Update limits
                if (indU.lt.MiindU) MiindU = indU
                if (indU1.lt.MiindU1) MiindU1 = indU1
                if (indL.lt.MiindL) MiindL = indL
                if (indL1.lt.MiindL1) MiindL1 = indL1
                if (indF.lt.MiindF) MiindF = indF
                if (indU.gt.MindU) MindU = indU
                if (indU1.gt.MindU1) MindU1 = indU1
                if (indL.gt.MindL) MindL = indL
                if (indL1.gt.MindL1) MindL1 = indL1
                if (indF.gt.MindF) MindF = indF

              ! Second round
              else

                ! Advance counter and store
                ii = ii + 1
                Atom(ia)%trano(jtran)%trani(iti)% &
                   indB(indL1,indL,indF,indU1,indU) = ii

                ! Advance nchlt if proceeds
                if (iL.ne.iL1) &
                  nnchlt = nnchlt + 1

              end if ! Round

            end do ! iL1
          end do ! iMl1
        end do ! iL
      end do ! iMl
                                !
                                ! Recover indentation
                                !

                              end do ! iU1
                            end do ! iMu1
                          end do ! iU
                        end do ! iMu
                      end do ! mF
                    end do ! iMf
                  end do ! Round i1

                  ! If non-coherent lower term
                  if (NCHLT) then

                    ! Add count of coherent terms
                    Atom(ia)%trano(jtran)%trani(iti)%nchlt = nnchlt

                  ! All coherent
                  else

                    ! Set count to zero
                    Atom(ia)%trano(jtran)%trani(iti)%nchlt = 0

                  end if ! If non-coherent lower term

                end do ! Input transitions

              end if ! PRD transition and storing redistribution

            end do ! Transitions

          end if ! We expect points with magnetic field
        end if ! Doing polarization

      end do ! Atoms

      end subroutine set_atom_indexes

!#####################################################################
!#####################################################################
!#####################################################################

      end module ratom_mod
