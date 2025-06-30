      !> Frees allocated memory
      module free_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     28/06/2022
!  Last version:
!     26/06/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     26/06/2025:    V4.0.3 - Updated due to changes in the Red_class
!                             and Redb_class structures (TdPA)
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
!  free_local_Atom
!    Deallocate quantities in Atom_class generated in a hanle call
!
!  free_los_geom
!    Deallocate LOS quantities in Geometry_class generated in a hanle
!  call
!
!  free_local_geom
!    Deallocate quantities in Geometry_class generated in a hanle call
!
!  free_local_CGF
!    Deallocate quantities in Continuum_class, Geometry_class, and
!  Frequency_class generated during a call to hanle
!
!  free_local
!    Deallocate quantities in Atom_class, Continuum_class,
!  Geometry_class, and Frequency_class generated during a call to
!  hanle
!
!  free_Atmo
!    Deallocate the model atmosphere, partially or completely
!
!  free_gpop
!    Deallocate the total populations of atoms and molecules
!
!  free_lpop
!    Deallocate the level population and density matrices of atoms
!
!  free_mol
!    Deallocate populations and interpolated partition function and
!  chemical equilibrium quantities of molecules
!
!  free_cols
!    Deallocate all collisions for atoms
!
!  free_damp
!    Deallocate all damping quantities for atoms
!
!  free_B
!    Deallocate magnetic field quantities
!
!  free_hanle
!    Deallocate quantities allocated in the call to hanle
!
!  free_pixel
!    Deallocate pixel wise quantities, including atoms, molecules,
!  and magnetic field quantities
!
!  free_cle_node
!    Deallocate quantities allocated in a node of a CLE calculation,
!  including atoms, molecules, model atmosphere, magnetic field, and
!  geometry
!
!  free_atom_full
!    Deallocate completely the data in an Atom_class array
!
!  free_LTElines_full
!    Deallocate completely the data in an LTEline_class array
!
!  free_mol_full
!    Deallocate completely the data in a Mol_class array
!
!  free_inv_solution
!    Deallocate the solution data in Solution_F_class used in the
!  inversion
!
!  free_warr
!    Deallocate the redistribution function
!
!  free_e2ord
!    Deallocate arrays to store the emissivity in PRD lines
!
!  free_ifreq
!    Deallocate the input frequency axis for PRD calculations
!
!  free_norm
!    Deallocate the profile normalization data
!
!  free_1stord
!    Deallocate the profile normalization data used to compute the
!  PRD emissivity
!
!  free_red
!    Deallocate all data in the Red_class structure
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate quantities in Atom_class generated in a hanle
      !! call\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine free_local_Atom(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ii,jj,kk


      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Eigenvalues
        if (allocated(Atom(ii)%eval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%eval)
          deallocate(Atom(ii)%eval)
        end if

        ! Eigenvectors
        if (allocated(Atom(ii)%evec)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%evec)
          deallocate(Atom(ii)%evec)
        end if

        ! Photoionization (thermal part)
        do jj=1,Atom(ii)%nphot
          if (allocated(Atom(ii)%phot(jj)%TEI)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%phot(jj)%TEI)
            deallocate(Atom(ii)%phot(jj)%TEI)
          end if
        end do

        ! NCHLT
        if (allocated(Atom(ii)%NCHLT)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%NCHLT)
          deallocate(Atom(ii)%NCHLT)
        end if

        ! Dipole strength energy eigenbasis
        if (allocated(Atom(ii)%rdipev)) then

          ! For each element
          do jj=lbound(Atom(ii)%rdipev,1),ubound(Atom(ii)%rdipev,1)

            ! Energy representation dipole strength
            if (allocated(Atom(ii)%rdipev(jj)%rdipev)) then

              ! For each element
              do kk=lbound(Atom(ii)%rdipev(jj)%rdipev,1), &
                    ubound(Atom(ii)%rdipev(jj)%rdipev,1)

                ! Dipole strength
                if (allocated(Atom(ii)%rdipev(jj)% &
                                       rdipev(kk)%rdip)) then
                  MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj)% &
                                                  rdipev(kk)%rdip)
                  deallocate(Atom(ii)%rdipev(jj)%rdipev(kk)%rdip)
                end if

                ! Free element space
                MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj)% &
                                                     rdipev(kk))

              end do ! Elements

              ! Array
              deallocate(Atom(ii)%rdipev(jj)%rdipev)

            end if ! Dipole strength

            ! Free element space
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj))

          end do ! Elements

          ! Array
          deallocate(Atom(ii)%rdipev)

        end if

      end do ! Atoms

      return

      end subroutine free_local_Atom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate LOS quantities in Geometry_class generated in a
      !! hanle call\n
      !!  Geom(Geometry_class): Structure with geometric data
      subroutine free_los_geom(Geom)

      ! I/O

      type(Geometry_class), intent(inout):: Geom


      !
      ! Geometry
      !

      ! TKQ S frame for LOS
      MRAMc = MRAMc - 1d-6*sizeof(Geom%TSL)
      deallocate(Geom%TSL)
      nullify(Geom%TSL)

      ! TKQ B frame for LOS
      MRAMc = MRAMc - 1d-6*sizeof(Geom%TBL)
      deallocate(Geom%TBL)
      nullify(Geom%TBL)

      return

      end subroutine free_los_geom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate quantities in Geometry_class generated in a hanle
      !! call\n
      !!  Geom(Geometry_class): Structure with geometric data
      subroutine free_local_geom(Geom)

      ! I/O

      type(Geometry_class), intent(inout):: Geom


      !
      ! Geometry
      !

      ! TKQ B frame
      if (associated(Geom%TB)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TB)
        deallocate(Geom%TB)
        nullify(Geom%TB)
        if (PRD.and..not.AV.and.axial) then
          MRAMc = MRAMc - 1d-6*sizeof(Geom%TBo)
          deallocate(Geom%TBo)
        end if
        nullify(Geom%TBo)
      end if

      ! TKQ S frame for LOS
      if (associated(Geom%TSL)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TSL)
        deallocate(Geom%TSL)
        nullify(Geom%TSL)
      end if

      ! TKQ B frame for LOS
      if (associated(Geom%TBL)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TBL)
        deallocate(Geom%TBL)
        nullify(Geom%TBL)
      end if

      ! Indexes mapping
      if (allocated(Geom%i_geom)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%i_geom)
        deallocate(Geom%i_geom)
      end if

      ! Index inclination
      if (allocated(Geom%ithv)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%ithv)
        deallocate(Geom%ithv)
      end if

      ! Index azimuth
      if (allocated(Geom%iphv)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%iphv)
        deallocate(Geom%iphv)
      end if

      return

      end subroutine free_local_geom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate quantities on Continuum_class, Geometry_class, and
      !! Frequency_class generated during a call to hanle\n
      !!  Cont(Continuum_class): Structure with background opacity
      !!                         data\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine free_local_CGF(Cont,Geom,Frec)

      ! I/O

      type(Continuum_class), intent(inout):: Cont
      type(Geometry_class), intent(inout):: Geom
      type(Frequency_class), intent(inout):: Frec


      !
      ! Continuum
      !
      if (allocated(Cont%c)) then
        BRAMc = 0d0
        deallocate(Cont%c)
      end if

      !
      ! Geometry
      !
      call free_local_geom(Geom)

      !
      ! Frequency
      !

      ! Cubic frequency
      if (allocated(Frec%omega3)) then
        PRAMc = PRAMc - 1d-6*sizeof(Frec%omega3)
        deallocate(Frec%omega3)
      end if

      ! Frequency exponential
      if (associated(Frec%exu)) then
        PRAMc = PRAMc - 1d-6*sizeof(Frec%exu)
        deallocate(Frec%exu)
        nullify(Frec%exu)
      end if

      return

      end subroutine free_local_CGF

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate quantities in Atom_class, Continuum_class,
      !! Geometry_class, and Frequency_class generated during a call
      !! to hanle\n
      !!    Atom(Atom_class(:)): Structures with atomic data
      !!  Cont(Continuum_class): Structure with background opacity
      !!                         data\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!  Frec(Frequency_class): Structure with frequency data
      subroutine free_local(Atom,Cont,Geom,Frec)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Continuum_class), intent(inout):: Cont
      type(Geometry_class), intent(inout):: Geom
      type(Frequency_class), intent(inout):: Frec

      ! Atom
      call free_local_Atom(Atom)

      ! Continuum, geometry, and frequency
      call free_local_CGF(Cont,Geom,Frec)

      return

      end subroutine free_local

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the model atmosphere, partially or completely\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!     full(logical): If deallocating everything
      subroutine free_Atmo(Atmo,full)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      logical, intent(in):: full

      ! Local

      integer:: ii


      !
      ! Pointers
      !

      ! If allocated group a
      if (Atmo%alloc_a) then

        ! Height
        if (associated(Atmo%z)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%z)
          deallocate(Atmo%z)
          nullify(Atmo%z)
        end if

        ! Temperature
        if (associated(Atmo%T)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%T)
          deallocate(Atmo%T)
          nullify(Atmo%T)
        end if

        ! Microturbulence
        if (associated(Atmo%vmi)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%vmi)
          deallocate(Atmo%vmi)
          nullify(Atmo%vmi)
        end if

        ! Velocity along x
        if (associated(Atmo%vx)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%vx)
          deallocate(Atmo%vx)
          nullify(Atmo%vx)
        end if

        ! Velocity along y
        if (associated(Atmo%vy)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%vy)
          deallocate(Atmo%vy)
          nullify(Atmo%vy)
        end if

        ! Velocity along z
        if (associated(Atmo%vz)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%vz)
          deallocate(Atmo%vz)
          nullify(Atmo%vz)
        end if

      ! If pointing
      else

        ! Nullify pointers
        if (associated(Atmo%z)) nullify(Atmo%z)
        if (associated(Atmo%T)) nullify(Atmo%T)
        if (associated(Atmo%vmi)) nullify(Atmo%vmi)
        if (associated(Atmo%vx)) nullify(Atmo%vx)
        if (associated(Atmo%vy)) nullify(Atmo%vy)
        if (associated(Atmo%vz)) nullify(Atmo%vz)

      end if ! Real allocation or pointing group a

      !
      ! If allocated group b
      if (Atmo%alloc_b) then

        ! B along x
        if (associated(Atmo%Bx)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%Bx)
          deallocate(Atmo%Bx)
          nullify(Atmo%Bx)
        end if

        ! B along y
        if (associated(Atmo%By)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%By)
          deallocate(Atmo%By)
          nullify(Atmo%By)
        end if

        ! B along z
        if (associated(Atmo%Bz)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%Bz)
          deallocate(Atmo%Bz)
          nullify(Atmo%Bz)
        end if

      ! If pointing
      else

        ! Nullify pointers
        if (associated(Atmo%Bx)) nullify(Atmo%Bx)
        if (associated(Atmo%By)) nullify(Atmo%By)
        if (associated(Atmo%Bz)) nullify(Atmo%Bz)

      end if ! Real allocation or pointing group b

      !
      ! Always pointers

      ! Auxiliar velocity alon x, y, and z
      if (associated(Atmo%vxa)) nullify(Atmo%vxa)
      if (associated(Atmo%vya)) nullify(Atmo%vya)
      if (associated(Atmo%vza)) nullify(Atmo%vza)

      !
      ! Always arrays or allocated pointers

      ! Vector of zeros
      if (associated(Atmo%zeros)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%zeros)
        deallocate(Atmo%zeros)
        nullify(Atmo%zeros)
      end if

      !
      ! Arrays!

      ! Total hydrogen
      if (allocated(Atmo%nHT)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nHT)
        deallocate(Atmo%nHT)
      end if

      ! Hydrogen minus
      if (allocated(Atmo%nHm)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nHm)
        deallocate(Atmo%nHm)
      end if

      ! Gas pressure
      if (allocated(Atmo%Pg)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pg)
        deallocate(Atmo%Pg)
      end if

      ! Mass density
      if (allocated(Atmo%rho)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%rho)
        deallocate(Atmo%rho)
      end if

      ! Electron pressure
      if (allocated(Atmo%Pe)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pe)
        deallocate(Atmo%Pe)
      end if

      ! Alternative height scale
      if (allocated(Atmo%zalt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%zalt)
        deallocate(Atmo%zalt)
      end if

      ! Total atomic hydrogen
      if (allocated(Atmo%nHa)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nHa)
        deallocate(Atmo%nHa)
      end if

      ! Opacity at reference wavelength
      if (allocated(Atmo%chi500)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%chi500)
        deallocate(Atmo%chi500)
      end if

      ! Electron number density
      if (allocated(Atmo%ne)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%ne)
        deallocate(Atmo%ne)
      end if

      ! Hydrogen level density
      if (allocated(Atmo%nh)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nh)
        deallocate(Atmo%nh)
      end if

      ! Helium level density
      if (allocated(Atmo%nhe)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%nhe)
        deallocate(Atmo%nhe)
      end if

      ! Longitudinal component of the velocity
      if (allocated(Atmo%vlos)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%vlos)
        deallocate(Atmo%vlos)
      end if

      ! Transversal component of the velocity
      if (allocated(Atmo%vpos)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%vpos)
        deallocate(Atmo%vpos)
      end if

      ! Azimuth in the plane of the sky for velocity
      if (allocated(Atmo%vphi)) then
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%vphi)
        deallocate(Atmo%vphi)
      end if

      !
      ! If fully deallocating
      !
      if (full) then

        ! Temperature tabulation
        if (allocated(Atmo%pT)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%pT)
          deallocate(Atmo%pT)
        end if

        ! Element data
        if (allocated(Atmo%ele)) then

          ! For all elements
          do ii=lbound(Atmo%ele,1),ubound(Atmo%ele,1)
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%ele(ii))
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%ele(ii)%pf)
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%ele(ii)%Ei)
            deallocate(Atmo%ele(ii)%pf)
            deallocate(Atmo%ele(ii)%Ei)
          end do

          ! Remove data
          deallocate(Atmo%ele)

        end if ! Element data

        ! Abundances
        if (allocated(Atmo%abund)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%abund)
          deallocate(Atmo%abund)
        end if

        ! Input ad-hoc radiation field
        if (allocated(Atmo%JKQin)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%JKQin)
          deallocate(Atmo%JKQin)
        end if

      end if ! Full deallocation

      return

      end subroutine free_Atmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the total populations of atoms and molecules\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data
      subroutine free_gpop(Atom,Atomb,Mol)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol

      ! Local

      integer:: ii


      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Population
        if (allocated(Atom(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%n)
          deallocate(Atom(ii)%n)
        end if

      end do ! Atoms

      ! For each passive atom
      do ii=1,nAb

        ! Population
        if (allocated(Atomb(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%n)
          deallocate(Atomb(ii)%n)
        end if

      end do ! Passive atoms

      !
      ! Molecules
      !

      ! For each molecule
      do ii=1,nM

        ! Population
        if (allocated(Mol(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%n)
          deallocate(Mol(ii)%n)
        end if

      end do

      return

      end subroutine free_gpop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the level population and density matrices of
      !! atoms\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms
      subroutine free_lpop(Atom,Atomb)

      ! I/O
      type(Atom_class), dimension(:), allocatable:: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb

      ! Local

      integer:: ii


      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Population
        if (allocated(Atom(ii)%popu)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%popu)
          deallocate(Atom(ii)%popu)
        end if

        ! LTE population
        if (allocated(Atom(ii)%populte)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%populte)
          deallocate(Atom(ii)%populte)
        end if

        ! Density matrix
        if (allocated(Atom(ii)%crho)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%crho)
          deallocate(Atom(ii)%crho)
        end if

        ! Small density matrix flag
        if (allocated(Atom(ii)%rhonull)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rhonull)
          deallocate(Atom(ii)%rhonull)
        end if

      end do ! Active atoms

      ! For each background atom
      do ii=1,nAb

        ! Population
        if (allocated(Atomb(ii)%popu)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%popu)
          deallocate(Atomb(ii)%popu)
        end if

        ! LTE population
        if (allocated(Atomb(ii)%populte)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%populte)
          deallocate(Atomb(ii)%populte)
        end if

      end do ! Passive atoms

      return

      end subroutine free_lpop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate populations and interpolated partition function
      !! and chemical equilibrium quantities of molecules\n
      !!  Mol(Mol_class(:)): Structures with molecular data
      subroutine free_mol(Mol)

      ! I/O

      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol

      ! Local

      integer:: ii


      ! For each molecule
      do ii=1,nM

        ! Population
        if (allocated(Mol(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%n)
          deallocate(Mol(ii)%n)
        end if

        ! Partition function
        if (allocated(Mol(ii)%pf)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%pf)
          deallocate(Mol(ii)%pf)
        end if

        ! Equilibrium constant
        if (allocated(Mol(ii)%eq)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%eq)
          deallocate(Mol(ii)%eq)
        end if

      end do ! Molecules

      return

      end subroutine free_mol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate all collisions for atoms\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine free_cols(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local
      integer:: ii


      !
      ! Atom
      !

      ! If allocated
      if (allocated(Atom)) then

        ! For each active atom
        do ii=1,size(Atom)

          ! Collisions between terms
          if (allocated(Atom(ii)%Ccoeff)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%Ccoeff)
            deallocate(Atom(ii)%Ccoeff)
          end if

          ! Collisions between levels
          if (allocated(Atom(ii)%CcoeffJ)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%CcoeffJ)
            deallocate(Atom(ii)%CcoeffJ)
          end if

          ! Depolarizing elastic collisions
          if (allocated(Atom(ii)%gk)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%gk)
            deallocate(Atom(ii)%gk)
          end if

          ! Indexing of collisions
          if (allocated(Atom(ii)%icol)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%icol)
            deallocate(Atom(ii)%icol)
          end if

        end do ! Atoms

      end if ! Allocated

      end subroutine free_cols

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate all damping quantities for atoms\n
      !!  Atom(Atom_class): Structure with atomic data
      subroutine free_damp(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ii


      !
      ! Atom
      !

      ! If allocated
      if (allocated(Atom)) then

        ! For each atom
        do ii=1,size(Atom)

          ! Term inverse lifetime
          if (allocated(Atom(ii)%damp)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%damp)
            deallocate(Atom(ii)%damp)
          end if

          ! Line damping
          if (allocated(Atom(ii)%ldamp)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%ldamp)
            deallocate(Atom(ii)%ldamp)
          end if

          ! Line elastic collision rate
          if (allocated(Atom(ii)%qel)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%qel)
            deallocate(Atom(ii)%qel)
          end if

        end do ! Atoms

      end if ! Allocated

      end subroutine free_damp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate magnetic field quantities\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine free_B(Bfield)

      ! I/O

      type(Bfield_class), intent(inout):: Bfield


      ! Magnetic field in polar coordinates
      if (allocated(Bfield%Bstrength)) then
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Bstrength)
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Btheta)
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Bphi)
        deallocate(Bfield%Bstrength,Bfield%Btheta,Bfield%Bphi)
      end if

      ! Magnetic field in LOS coordinates
      if (allocated(Bfield%Blos)) then
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Blos)
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Bpos)
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Azimuth)
        deallocate(Bfield%Blos,Bfield%Bpos,Bfield%Azimuth)
      end if

      return

      end subroutine free_B

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate quantities allocated in the call to hanle\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!     Bfield0(Bfield_class): Structure with magnetic field
      !!                            data for the unmagnetized call\n
      !!       Rho_old(Rhoc_class): Structure to store the density
      !!                            matrix of the previous iteration\n
      !!  JKQ_asym(dcmplex(:,:,:)): Extra asymmetry for the radiation
      !!                            field tensors\n
      subroutine free_hanle(Atom,Cont,Geom,Frec,Bfield0,Rho_old, &
                            JKQ_asym)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Continuum_class), intent(inout):: Cont
      type(Geometry_class), intent(inout):: Geom
      type(Frequency_class), intent(inout):: Frec
      type(Bfield_class), intent(inout):: Bfield0
      type(Rhoc_class), dimension(:), &
                        allocatable, intent(inout):: Rho_old
      complex(kind=8), dimension(:,:,:), &
                       allocatable, intent(inout):: JKQ_asym

      ! Local
      integer:: ia

      ! Free atom, cointinuum, geometry, and frequency data
      call free_local(Atom,Cont,Geom,Frec)

      ! Remove count of non-allocatables, as the structure is going
      ! to be leave behind in context
      MRAMc = MRAMc - 1d-6*sizeof(Cont)

      ! Free dummy zero magnetic field
      call free_B(Bfield0)

      ! Remove count of non-allocatables, as the structure is going
      ! to be leave behind in context
      MRAMc = MRAMc - 1d-6*sizeof(Bfield0)

      ! Rho_old allocated
      if (allocated(Rho_old)) then

        ! For each atom
        do ia=1,size(Rho_old)

          ! Density matrix
          if (allocated(Rho_old(ia)%crho)) then
            MRAMc = MRAMc - 1d-6*sizeof(Rho_old(ia)%crho)
            deallocate(Rho_old(ia)%crho)
          end if

        end do ! Atoms

        ! Remove count of non-allocatables, as the structure is going
        ! to be leave behind in context
        MRAMc = MRAMc - 1d-6*sizeof(Rho_old)

      end if ! Allocated Rho_old

      ! If ad-hoc radiation data
      if (allocated(JKQ_asym)) then
        RRAMc = RRAMc - 1d-6*sizeof(JKQ_asym)
        deallocate(JKQ_asym)
      end if

      end subroutine free_hanle

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory that depends on the pixel\n
      !> Deallocate pixel wise quantities, including atoms, molecules,
      !! and magnetic field quantities\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine free_pix(Atom,Atomb,Mol,Bfield)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Bfield_class), intent(inout):: Bfield

      ! Free global populations
      call free_gpop(Atom,Atomb,Mol)

      ! Free local populations
      call free_lpop(Atom,Atomb)

      ! Free molecule
      call free_mol(Mol)

      ! Free collisions
      call free_cols(Atom)

      ! Free damping
      call free_damp(Atom)
      call free_damp(Atomb)

      ! Free Magnetic field
      call free_B(Bfield)

      return

      end subroutine free_pix

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate quantities allocated in a node of a CLE
      !! calculation, including atoms, molecules, model atmosphere,
      !! magnetic field, and geometry\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      subroutine free_cle_node(Atom,Atomb,Mol,Atmo,Bfield,Geom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Geometry_class), intent(inout):: Geom

      ! Local

      integer:: ii,jj,kk


      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Total population
        if (allocated(Atom(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%n)
          deallocate(Atom(ii)%n)
        end if

        ! Level population
        if (allocated(Atom(ii)%popu)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%popu)
          deallocate(Atom(ii)%popu)
        end if

        ! LTE population
        if (allocated(Atom(ii)%populte)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%populte)
          deallocate(Atom(ii)%populte)
        end if

        ! Eigenvalues
        if (allocated(Atom(ii)%eval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%eval)
          deallocate(Atom(ii)%eval)
        end if

        ! Eigenvectors
        if (allocated(Atom(ii)%evec)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%evec)
          deallocate(Atom(ii)%evec)
        end if

        ! Density matrix
        if (allocated(Atom(ii)%crho)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%crho)
          deallocate(Atom(ii)%crho)
        end if

        ! Density matrix smallness flag
        if (allocated(Atom(ii)%rhonull)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rhonull)
          deallocate(Atom(ii)%rhonull)
        end if

        ! Inverse lifetime
        if (allocated(Atom(ii)%damp)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%damp)
          deallocate(Atom(ii)%damp)
        end if

        ! Line damping
        if (allocated(Atom(ii)%ldamp)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%ldamp)
          deallocate(Atom(ii)%ldamp)
        end if

        ! Line elastic rates
        if (allocated(Atom(ii)%qel)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%qel)
          deallocate(Atom(ii)%qel)
        end if

        ! Collisions between terms
        if (allocated(Atom(ii)%Ccoeff)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%Ccoeff)
          deallocate(Atom(ii)%Ccoeff)
        end if

        ! Collisions between levels
        if (allocated(Atom(ii)%CcoeffJ)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%CcoeffJ)
          deallocate(Atom(ii)%CcoeffJ)
        end if

        ! Depolarizing elastic collisions
        if (allocated(Atom(ii)%gk)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%gk)
          deallocate(Atom(ii)%gk)
        end if

        ! Collisions indexing
        if (allocated(Atom(ii)%icol)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%icol)
          deallocate(Atom(ii)%icol)
        end if

        ! Photoionization (thermal part)
        do jj=1,Atom(ii)%nphot
          if (allocated(Atom(ii)%phot(jj)%TEI)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%phot(jj)%TEI)
            deallocate(Atom(ii)%phot(jj)%TEI)
          end if
        end do

        ! NCHLT
        if (allocated(Atom(ii)%NCHLT)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%NCHLT)
          deallocate(Atom(ii)%NCHLT)
        end if

        ! Dipole strength energy eigenbasis
        if (allocated(Atom(ii)%rdipev)) then

          ! For each element
          do jj=lbound(Atom(ii)%rdipev,1),ubound(Atom(ii)%rdipev,1)

            ! Energy representation dipole strength
            if (allocated(Atom(ii)%rdipev(jj)%rdipev)) then

              ! For each element
              do kk=lbound(Atom(ii)%rdipev(jj)%rdipev,1), &
                    ubound(Atom(ii)%rdipev(jj)%rdipev,1)

                ! Dipole strength
                if (allocated(Atom(ii)%rdipev(jj)% &
                                       rdipev(kk)%rdip)) then
                  MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj)% &
                                                  rdipev(kk)%rdip)
                  deallocate(Atom(ii)%rdipev(jj)%rdipev(kk)%rdip)
                end if

                ! Free element space
                MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj)% &
                                                     rdipev(kk))

              end do ! Elements

              ! Array
              deallocate(Atom(ii)%rdipev(jj)%rdipev)

            end if ! Dipole strength

            ! Free element space
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj))

          end do ! Elements

          ! Array
          deallocate(Atom(ii)%rdipev)

        end if ! Allocated rdipev

      end do ! Active atoms

      ! For each passive atom
      do ii=1,nAb

        ! Total population
        if (allocated(Atomb(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%n)
          deallocate(Atomb(ii)%n)
        end if

        ! Level population
        if (allocated(Atomb(ii)%popu)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%popu)
          deallocate(Atomb(ii)%popu)
        end if

        ! LTE population
        if (allocated(Atomb(ii)%populte)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%populte)
          deallocate(Atomb(ii)%populte)
        end if

        ! Inverse lifetime
        if (allocated(Atomb(ii)%damp)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%damp)
          deallocate(Atomb(ii)%damp)
        end if

        ! Line damping
        if (allocated(Atomb(ii)%ldamp)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%ldamp)
          deallocate(Atomb(ii)%ldamp)
        end if

        ! Elastic collision rate
        if (allocated(Atomb(ii)%qel)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atomb(ii)%qel)
          deallocate(Atomb(ii)%qel)
        end if

      end do ! Passive atoms

      !
      ! Molecules
      !

      ! For each molecule
      do ii=1,nM

        ! Total population
        if (allocated(Mol(ii)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%n)
          deallocate(Mol(ii)%n)
        end if

        ! Partition function
        if (allocated(Mol(ii)%pf)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%pf)
          deallocate(Mol(ii)%pf)
        end if

        ! Chemical equilibrium constant
        if (allocated(Mol(ii)%eq)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ii)%eq)
          deallocate(Mol(ii)%eq)
        end if

      end do ! Molecules


      !
      ! Atmosphere
      !

      ! Nullify pointers
      nullify(Atmo%T,Atmo%vmi,Atmo%vx,Atmo%vy,Atmo%vz)
      if (associated(Atmo%Bx)) nullify(Atmo%Bx,Atmo%By,Atmo%Bz)

      !
      ! Magnetic field
      !
      if (allocated(Bfield%Bstrength)) then
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Bstrength)
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Btheta)
        MRAMc = MRAMc - 1d-6*sizeof(Bfield%Bphi)
        deallocate(Bfield%Bstrength,Bfield%Btheta,Bfield%Bphi)
      end if

      !
      ! Geometry
      !

      ! TKQ in S frame for quadrature
      if (associated(Geom%TS)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TS)
        deallocate(Geom%TS)
        nullify(Geom%TS)
      end if

      ! TKQ in B frame for quadrature
      if (associated(Geom%TB)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TB)
        deallocate(Geom%TB)
        nullify(Geom%TB)
      end if

      ! TKQ in S frame for line of sight
      if (associated(Geom%TSL)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TSL)
        deallocate(Geom%TSL)
        nullify(Geom%TSL)
      end if

      ! TKQ in B frame for line of sight
      if (associated(Geom%TBL)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TBL)
        deallocate(Geom%TBL)
        nullify(Geom%TBL)
      end if

      return

      end subroutine free_cle_node

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate completely the data in an Atom_class array\n
      !!   Atom(Atom_class(:)): Structures with atomic data
      subroutine free_atom_full(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ia,ii,jj,kk

      ! Pointers

      type(tmp_col_box_class), pointer:: p_col,p_col_p
      type(Tbox_class), pointer:: p_T,p_T_p


      ! Nullify pointers
      nullify(p_col,p_col_p,p_T,p_T_p)

      ! Return if empty structure
      if (.not.allocated(Atom)) return

      ! For each atom
      do ia=1,size(Atom)

        !
        ! fflag
        !
        if (allocated(Atom(ia)%fflag)) then

          ! For each element
          do ii=lbound(Atom(ia)%fflag,1),ubound(Atom(ia)%fflag,1)

            ! Absence from master
            if (allocated(Atom(ia)%fflag(ii)%Mabsent)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fflag(ii)%Mabsent)
              deallocate(Atom(ia)%fflag(ii)%Mabsent)
            end if

            ! Wavelength absence
            if (allocated(Atom(ia)%fflag(ii)%Vabsent)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fflag(ii)%Vabsent)
              deallocate(Atom(ia)%fflag(ii)%Vabsent)
            end if

            ! Non-array memory
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fflag(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%fflag)

        end if ! fflag


        !
        ! fst
        !
        if (allocated(Atom(ia)%fst)) then

          ! For each element
          do ii=lbound(Atom(ia)%fst,1),ubound(Atom(ia)%fst,1)

            ! Lower level index
            if (allocated(Atom(ia)%fst(ii)%ilevell)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fst(ii)%ilevell)
              deallocate(Atom(ia)%fst(ii)%ilevell)
            end if

            ! Upper level index
            if (allocated(Atom(ia)%fst(ii)%ilevelu)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fst(ii)%ilevelu)
              deallocate(Atom(ia)%fst(ii)%ilevelu)
            end if

            ! FS transition indexing
            if (allocated(Atom(ia)%fst(ii)%irad)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fst(ii)%irad)
              deallocate(Atom(ia)%fst(ii)%irad)
            end if

            ! Aul Einstein coefficient
            if (allocated(Atom(ia)%fst(ii)%Aul)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fst(ii)%Aul)
              deallocate(Atom(ia)%fst(ii)%Aul)
            end if

            ! Blu Einstein coefficient
            if (allocated(Atom(ia)%fst(ii)%Blu)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fst(ii)%Blu)
              deallocate(Atom(ia)%fst(ii)%Blu)
            end if

            ! Non-array memory
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fst(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%fst)

        end if ! fst


        !
        ! phot
        !
        if (allocated(Atom(ia)%phot)) then

          ! For each element
          do ii=lbound(Atom(ia)%phot,1),ubound(Atom(ia)%phot,1)

            ! Absence from master
            if (allocated(Atom(ia)%phot(ii)%Mabsent)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%Mabsent)
              deallocate(Atom(ia)%phot(ii)%Mabsent)
            end if

            ! Initial frequency index from master
            if (allocated(Atom(ia)%phot(ii)%Mif0)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%Mif0)
              deallocate(Atom(ia)%phot(ii)%Mif0)
            end if

            ! Final frequency index from master
            if (allocated(Atom(ia)%phot(ii)%Mif1)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%Mif1)
              deallocate(Atom(ia)%phot(ii)%Mif1)
            end if

            ! Cross-section
            if (allocated(Atom(ia)%phot(ii)%alpha)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%alpha)
              deallocate(Atom(ia)%phot(ii)%alpha)
            end if

            ! Thermal rate
            if (allocated(Atom(ia)%phot(ii)%TEI)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%TEI)
              deallocate(Atom(ia)%phot(ii)%TEI)
            end if

            ! Input frequency array
            if (allocated(Atom(ia)%phot(ii)%infreq)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%infreq)
              deallocate(Atom(ia)%phot(ii)%infreq)
            end if

            ! Input cross-section array
            if (allocated(Atom(ia)%phot(ii)%inalpha)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%inalpha)
              deallocate(Atom(ia)%phot(ii)%inalpha)
            end if

            ! Initial index weight
            if (allocated(Atom(ia)%phot(ii)%MW0)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%MW0)
              deallocate(Atom(ia)%phot(ii)%MW0)
            end if

            ! Final index weight
            if (allocated(Atom(ia)%phot(ii)%MW1)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii)%MW1)
              deallocate(Atom(ia)%phot(ii)%MW1)
            end if

            ! Non-arrays
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%phot(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%phot)

        end if ! phot


        !
        ! trano
        !
        if (allocated(Atom(ia)%trano)) then

          ! For each element
          do ii=lbound(Atom(ia)%trano,1),ubound(Atom(ia)%trano,1)

            ! Indexing of input transitions
            if (allocated(Atom(ia)%trano(ii)%indT)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%trano(ii)%indT)
              deallocate(Atom(ia)%trano(ii)%indT)
            end if

            ! Indexing of magnetic components
            if (allocated(Atom(ia)%trano(ii)%indB)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%trano(ii)%indB)
              deallocate(Atom(ia)%trano(ii)%indB)
            end if

            ! Indexing of FS components
            if (allocated(Atom(ia)%trano(ii)%indNB)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%trano(ii)%indNB)
              deallocate(Atom(ia)%trano(ii)%indNB)
            end if

            ! If input transitions
            if (allocated(Atom(ia)%trano(ii)%trani)) then

              ! For each input transition
              do jj=lbound(Atom(ia)%trano(ii)%trani,1), &
                    ubound(Atom(ia)%trano(ii)%trani,1)

                ! Indexing of magnetic components
                if (allocated(Atom(ia)%trano(ii)%trani(jj)%indB)) then
                  MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)% &
                                              trano(ii)% &
                                              trani(jj)%indB)
                  deallocate(Atom(ia)%trano(ii)%trani(jj)%indB)
                end if

                ! Indexing of non-magnetic components
                if (allocated(Atom(ia)%trano(ii)% &
                                       trani(jj)%indNB)) then
                  MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)% &
                                              trano(ii)% &
                                              trani(jj)%indNB)
                  deallocate(Atom(ia)%trano(ii)%trani(jj)%indNB)
                end if

                ! Array
                MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%trano(ii)% &
                                                     trani(jj))
              end do ! Input transitions

              ! Array
              deallocate(Atom(ia)%trano(ii)%trani)

            end if ! Input transitions

            ! Non-arrays
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%trano(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%trano)

        end if ! trano


        !
        ! tranoI
        !
        if (allocated(Atom(ia)%tranoI)) then

          ! For each element
          do ii=lbound(Atom(ia)%tranoI,1),ubound(Atom(ia)%tranoI,1)

            ! Input transition indexing
            if (allocated(Atom(ia)%tranoI(ii)%indT)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%tranoI(ii)%indT)
              deallocate(Atom(ia)%tranoI(ii)%indT)
            end if

            ! Non-arrays
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%tranoI(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%tranoI)

        end if ! tranoI


        !
        ! Ccoeff_special
        !
        if (associated(Atom(ia)%Ccoeff_special)) then

          ! Do while there is data
          do while (associated(Atom(ia)%Ccoeff_special))

            ! Initialize
            p_col => Atom(ia)%Ccoeff_special
            nullify(p_col_p)

            ! Go to the last one
            do while (associated(p_col%next))
              p_col_p => p_col
              p_col => p_col%next
            end do ! Forward navigation

            ! Deallocate content
            MRAMc = MRAMc - 1d-6*sizeof(p_col%C)
            deallocate(p_col%C)

            ! Remove the last one
            if (associated(p_col_p)) then
              nullify(p_col_p%next)
              nullify(p_col_p)
            ! Remove the first one
            else
              deallocate(Atom(ia)%Ccoeff_special)
              nullify(Atom(ia)%Ccoeff_special)
              nullify(p_col)
            end if

            ! Remove the non-arrays
            MRAMc = MRAMc - 12d-6

          end do ! There is data

        end if ! Ccoeff_special


        !
        ! inelas
        !
        if (allocated(Atom(ia)%inelas)) then

          ! For each element
          do ii=lbound(Atom(ia)%inelas,1),ubound(Atom(ia)%inelas,1)

            ! Collisional rate
            if (allocated(Atom(ia)%inelas(ii)%Cul)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%inelas(ii)%Cul)
              deallocate(Atom(ia)%inelas(ii)%Cul)
            end if

            ! Non-arrays
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%inelas(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%inelas)

        end if ! inelas


        !
        ! elas
        !
        if (allocated(Atom(ia)%elas)) then

          ! For each element
          do ii=lbound(Atom(ia)%elas,1),ubound(Atom(ia)%elas,1)

            ! Elastic collision data
            if (allocated(Atom(ia)%elas(ii)%datum)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%elas(ii)%datum)
              deallocate(Atom(ia)%elas(ii)%datum)
            end if

            ! Non-arrays
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%elas(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%elas)

        end if ! elas


        !
        ! Tbox
        !
        if (associated(Atom(ia)%Tbox)) then

          ! For every collision
          do while (associated(Atom(ia)%Tbox))

            ! Initialize
            p_T => Atom(ia)%Tbox

            ! Go to the last one
            do while (associated(p_T%next))
              p_T_p => p_T
              p_T => p_T%next
            end do ! Forward navigation

            ! Deallocate temp
            MRAMc = MRAMc - 1d-6*sizeof(p_T%temp)
            deallocate(p_T%temp)

            ! Remove the last one
            if (associated(p_T_p)) then
              nullify(p_T)
              deallocate(p_T_p%next)
              nullify(p_T_p%next)
              nullify(p_T_p)
            else
              nullify(p_T)
              deallocate(Atom(ia)%Tbox)
              nullify(Atom(ia)%Tbox)
            end if

            ! Non-array data
            MRAMc = MRAMc - 20d-6

          end do ! Every collision

        end if ! Tbox


        !
        ! rdip
        !
        if (allocated(Atom(ia)%rdip)) then

          ! For each element
          do ii=lbound(Atom(ia)%rdip,1),ubound(Atom(ia)%rdip,1)

            ! Dipole strength
            if (allocated(Atom(ia)%rdip(ii)%rdip)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rdip(ii)%rdip)
              deallocate(Atom(ia)%rdip(ii)%rdip)
            end if

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%rdip)

        end if ! rdip

        ! dripev
        if (allocated(Atom(ia)%rdipev)) then


          ! For each element
          do ii=lbound(Atom(ia)%rdipev,1),ubound(Atom(ia)%rdipev,1)

            ! Energy representation dipole strength
            if (allocated(Atom(ia)%rdipev(ii)%rdipev)) then

              ! For each element
              do jj=lbound(Atom(ia)%rdipev(ii)%rdipev,1), &
                    ubound(Atom(ia)%rdipev(ii)%rdipev,1)

                ! Dipole strength
                if (allocated(Atom(ia)%rdipev(ii)% &
                                       rdipev(jj)%rdip)) then
                  MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rdipev(ii)% &
                                                  rdipev(jj)%rdip)
                  deallocate(Atom(ia)%rdipev(ii)%rdipev(jj)%rdip)
                end if

                ! Free element space
                MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj)% &
                                                     rdipev(kk))

              end do ! Elements

              ! Array
              deallocate(Atom(ia)%rdipev(ii)%rdipev)

            end if ! Dipole strength

            ! Free element space
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ii)%rdipev(jj))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%rdipev)

        end if ! rdipev


        !
        ! irho
        !
        if (allocated(Atom(ia)%irho)) then

          ! For each element
          do ii=lbound(Atom(ia)%irho,1),ubound(Atom(ia)%irho,1)

            ! Level indexing
            if (allocated(Atom(ia)%irho(ii)%irho_ij)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%irho(ii)%irho_ij)
              deallocate(Atom(ia)%irho(ii)%irho_ij)
            end if

            ! Magnetic sublevels indexing
            if (allocated(Atom(ia)%irho(ii)%jM)) then
              MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%irho(ii)%jM)
              deallocate(Atom(ia)%irho(ii)%jM)
            end if

            ! Density matrix indexing
            if (allocated(Atom(ia)%irho(ii)%Jrho)) then

              ! J' indexes
              do kk=lbound(Atom(ia)%irho(ii)%Jrho,2), &
                    ubound(Atom(ia)%irho(ii)%Jrho,2)

                ! J indexes
                do jj=lbound(Atom(ia)%irho(ii)%Jrho,1), &
                      ubound(Atom(ia)%irho(ii)%Jrho,1)

                  ! KQ indexing
                  if (allocated(Atom(ia)%irho(ii)% &
                                         Jrho(jj,kk)%kq)) then
                    MRAMc = MRAMc - &
                            1d-6*sizeof(Atom(ia)%irho(ii)% &
                                                 Jrho(jj,kk)%kq)
                    deallocate(Atom(ia)%irho(ii)%Jrho(jj,kk)%kq)
                  end if ! KQ indexing

                end do ! J indexes
              end do ! J' indexes

              ! Array
              deallocate(Atom(ia)%irho(ii)%Jrho)

            end if ! Density matrix indexing

            ! Non-array
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%irho(ii))

          end do ! Elements

          ! Array
          deallocate(Atom(ia)%irho)

        end if ! irho

        !
        ! Normal arrays
        !

        !
        ! bbspecin
        if (allocated(Atom(ia)%bbspecin)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%bbspecin)
          deallocate(Atom(ia)%bbspecin)
        end if

        !
        ! bfspecin
        if (allocated(Atom(ia)%bfspecin)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%bfspecin)
          deallocate(Atom(ia)%bfspecin)
        end if

        !
        ! lemiss2
        if (allocated(Atom(ia)%lemiss2)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%lemiss2)
          deallocate(Atom(ia)%lemiss2)
        end if

        !
        ! splitf
        if (allocated(Atom(ia)%splitf)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%splitf)
          deallocate(Atom(ia)%splitf)
        end if

        !
        ! rhonull
        if (allocated(Atom(ia)%rhonull)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rhonull)
          deallocate(Atom(ia)%rhonull)
        end if

        !
        ! NCHLT
        if (allocated(Atom(ia)%NCHLT)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%NCHLT)
          deallocate(Atom(ia)%NCHLT)
        end if

        !
        ! nfreqt
        if (allocated(Atom(ia)%nfreqt)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%nfreqt)
          deallocate(Atom(ia)%nfreqt)
        end if

        !
        ! nfreqtc
        if (allocated(Atom(ia)%nfreqtc)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%nfreqtc)
          deallocate(Atom(ia)%nfreqtc)
        end if

        !
        ! nJ
        if (allocated(Atom(ia)%nJ)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%nJ)
          deallocate(Atom(ia)%nJ)
        end if

        !
        ! stage
        if (allocated(Atom(ia)%stage)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%stage)
          deallocate(Atom(ia)%stage)
        end if

        !
        ! broad_type
        if (allocated(Atom(ia)%broad_type)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%broad_type)
          deallocate(Atom(ia)%broad_type)
        end if

        !
        ! term
        if (allocated(Atom(ia)%term)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%term)
          deallocate(Atom(ia)%term)
        end if

        !
        ! sublevel
        if (allocated(Atom(ia)%sublevel)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%sublevel)
          deallocate(Atom(ia)%sublevel)
        end if

        !
        ! nfreqph
        if (allocated(Atom(ia)%nfreqph)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%nfreqph)
          deallocate(Atom(ia)%nfreqph)
        end if

        !
        ! col_type
        if (allocated(Atom(ia)%col_type)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%col_type)
          deallocate(Atom(ia)%col_type)
        end if

        !
        ! if0
        if (allocated(Atom(ia)%if0)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%if0)
            deallocate(Atom(ia)%if0)
        end if

        !
        ! if1
        if (allocated(Atom(ia)%if1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%if1)
          deallocate(Atom(ia)%if1)
        end if

        !
        ! ifst
        if (allocated(Atom(ia)%ifst)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ifst)
          deallocate(Atom(ia)%ifst)
        end if

        !
        ! itrano
        if (allocated(Atom(ia)%itrano)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%itrano)
          deallocate(Atom(ia)%itrano)
        end if

        !
        ! rif0
        if (allocated(Atom(ia)%rif0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rif0)
          deallocate(Atom(ia)%rif0)
        end if

        !
        ! rif1
        if (allocated(Atom(ia)%rif1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rif1)
          deallocate(Atom(ia)%rif1)
        end if

        !
        ! rif20
        if (allocated(Atom(ia)%rif20)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rif20)
          deallocate(Atom(ia)%rif20)
        end if

        !
        ! rif21
        if (allocated(Atom(ia)%rif21)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rif21)
          deallocate(Atom(ia)%rif21)
        end if

        !
        ! sbif0
        if (allocated(Atom(ia)%sbif0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%sbif0)
          deallocate(Atom(ia)%sbif0)
        end if

        !
        ! sbif1
        if (allocated(Atom(ia)%sbif1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%sbif1)
          deallocate(Atom(ia)%sbif1)
        end if

        !
        ! sfif0
        if (allocated(Atom(ia)%sfif0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%sfif0)
          deallocate(Atom(ia)%sfif0)
        end if

        !
        ! sfif1
        if (allocated(Atom(ia)%sfif1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%sfif1)
          deallocate(Atom(ia)%sfif1)
        end if

        !
        ! ilf0
        if (allocated(Atom(ia)%ilf0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ilf0)
          deallocate(Atom(ia)%ilf0)
        end if

        !
        ! ilf1
        if (allocated(Atom(ia)%ilf1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ilf1)
          deallocate(Atom(ia)%ilf1)
        end if

        !
        ! ipf0
        if (allocated(Atom(ia)%ipf0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ipf0)
          deallocate(Atom(ia)%ipf0)
        end if

        !
        ! ipf1
        if (allocated(Atom(ia)%ipf1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ipf1)
          deallocate(Atom(ia)%ipf1)
        end if

        !
        ! Kcut
        if (allocated(Atom(ia)%Kcut)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Kcut)
          deallocate(Atom(ia)%Kcut)
        end if

        !
        ! Krad
        if (allocated(Atom(ia)%Krad)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Krad)
          deallocate(Atom(ia)%Krad)
        end if

        !
        ! tif0
        if (allocated(Atom(ia)%tif0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%tif0)
          deallocate(Atom(ia)%tif0)
        end if

        !
        ! tif1
        if (allocated(Atom(ia)%tif1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%tif1)
          deallocate(Atom(ia)%tif1)
        end if

        !
        ! irad
        if (allocated(Atom(ia)%irad)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%irad)
          deallocate(Atom(ia)%irad)
        end if

        !
        ! icol
        if (allocated(Atom(ia)%icol)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%icol)
          deallocate(Atom(ia)%icol)
        end if

        !
        ! iphot
        if (allocated(Atom(ia)%iphot)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%iphot)
          deallocate(Atom(ia)%iphot)
        end if

        !
        ! Mif0
        if (allocated(Atom(ia)%Mif0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Mif0)
          deallocate(Atom(ia)%Mif0)
        end if

        !
        ! Mif1
        if (allocated(Atom(ia)%Mif1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Mif1)
          deallocate(Atom(ia)%Mif1)
        end if

        !
        ! CMif0
        if (allocated(Atom(ia)%CMif0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%CMif0)
          deallocate(Atom(ia)%CMif0)
        end if

        !
        ! CMif1
        if (allocated(Atom(ia)%CMif1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%CMif1)
          deallocate(Atom(ia)%CMif1)
        end if

        !
        ! ifst_ij
        if (allocated(Atom(ia)%ifst_ij)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ifst_ij)
          deallocate(Atom(ia)%ifst_ij)
        end if

        !
        ! nblk
        if (allocated(Atom(ia)%nblk)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%nblk)
          deallocate(Atom(ia)%nblk)
        end if

        !
        ! fcflag
        if (allocated(Atom(ia)%fcflag)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%fcflag)
          deallocate(Atom(ia)%fcflag)
        end if

        !
        ! iJval
        if (allocated(Atom(ia)%iJval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%iJval)
          deallocate(Atom(ia)%iJval)
        end if

        !
        ! rLval
        if (allocated(Atom(ia)%rLval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rLval)
          deallocate(Atom(ia)%rLval)
        end if

        !
        ! Sval
        if (allocated(Atom(ia)%Sval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Sval)
          deallocate(Atom(ia)%Sval)
        end if

        !
        ! TRfreq
        if (allocated(Atom(ia)%TRfreq)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%TRfreq)
          deallocate(Atom(ia)%TRfreq)
        end if

        !
        ! Dwvl
        if (allocated(Atom(ia)%Dwvl)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Dwvl)
          deallocate(Atom(ia)%Dwvl)
        end if

        !
        ! Dwvlc
        if (allocated(Atom(ia)%Dwvlc)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Dwvlc)
          deallocate(Atom(ia)%Dwvlc)
        end if

        !
        ! Dfreq
        if (allocated(Atom(ia)%Dfreq)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Dfreq)
          deallocate(Atom(ia)%Dfreq)
        end if

        !
        ! n
        if (allocated(Atom(ia)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%n)
          deallocate(Atom(ia)%n)
        end if

        !
        ! broad_stark
        if (allocated(Atom(ia)%broad_stark)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%broad_stark)
          deallocate(Atom(ia)%broad_stark)
        end if

        !
        ! deg
        if (allocated(Atom(ia)%deg)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%deg)
          deallocate(Atom(ia)%deg)
        end if

        !
        ! W0
        if (allocated(Atom(ia)%W0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%W0)
          deallocate(Atom(ia)%W0)
        end if

        !
        ! W1
        if (allocated(Atom(ia)%W1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%W1)
          deallocate(Atom(ia)%W1)
        end if

        !
        ! gL
        if (allocated(Atom(ia)%gL)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%gL)
          deallocate(Atom(ia)%gL)
        end if

        !
        ! damp
        if (allocated(Atom(ia)%damp)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%damp)
          deallocate(Atom(ia)%damp)
        end if

        !
        ! FSfreq
        if (allocated(Atom(ia)%FSfreq)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%FSfreq)
          deallocate(Atom(ia)%FSfreq)
        end if

        !
        ! rJval
        if (allocated(Atom(ia)%rJval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%rJval)
          deallocate(Atom(ia)%rJval)
        end if

        !
        ! Ecoeff
        if (allocated(Atom(ia)%Ecoeff)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Ecoeff)
          deallocate(Atom(ia)%Ecoeff)
        end if

        !
        ! broad_args
        if (allocated(Atom(ia)%broad_args)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%broad_args)
          deallocate(Atom(ia)%broad_args)
        end if

        !
        ! ldamp
        if (allocated(Atom(ia)%ldamp)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%ldamp)
          deallocate(Atom(ia)%ldamp)
        end if

        !
        ! popu
        if (allocated(Atom(ia)%popu)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%popu)
          deallocate(Atom(ia)%popu)
        end if

        !
        ! populte
        if (allocated(Atom(ia)%populte)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%populte)
          deallocate(Atom(ia)%populte)
        end if

        !
        ! depar
        if (allocated(Atom(ia)%depar)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%depar)
          deallocate(Atom(ia)%depar)
        end if

        !
        ! MW0
        if (allocated(Atom(ia)%MW0)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%MW0)
          deallocate(Atom(ia)%MW0)
        end if

        !
        ! MW1
        if (allocated(Atom(ia)%MW1)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%MW1)
          deallocate(Atom(ia)%MW1)
        end if

        !
        ! qel
        if (allocated(Atom(ia)%qel)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%qel)
          deallocate(Atom(ia)%qel)
        end if

        !
        ! Ccoeff
        if (allocated(Atom(ia)%Ccoeff)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%Ccoeff)
          deallocate(Atom(ia)%Ccoeff)
        end if

        !
        ! CcoeffJ
        if (allocated(Atom(ia)%CcoeffJ)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%CcoeffJ)
          deallocate(Atom(ia)%CcoeffJ)
        end if

        !
        ! eval
        if (allocated(Atom(ia)%eval)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%eval)
          deallocate(Atom(ia)%eval)
        end if

        !
        ! gk
        if (allocated(Atom(ia)%gk)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%gk)
          deallocate(Atom(ia)%gk)
        end if

        !
        ! evec
        if (allocated(Atom(ia)%evec)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%evec)
          deallocate(Atom(ia)%evec)
        end if

        !
        ! crho
        if (allocated(Atom(ia)%crho)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom(ia)%crho)
          deallocate(Atom(ia)%crho)
        end if

        ! Non-arrays
        MRAMc = MRAMc - 1d-6*sizeof(Atom(ia))

      end do ! Atoms

      ! Deallocate array
      deallocate(Atom)

      end subroutine free_atom_full

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate completely the data in an Atom_class array\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      subroutine free_LTElines_full(LTElines)

      ! I/O

      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTElines

      ! Local

      integer:: ia


      ! Return if empty structure
      if (.not.allocated(LTElines)) return

      ! For each LTE line
      do ia=1,size(LTElines)

        ! damp
        if (allocated(LTElines(ia)%damp)) then
          MRAMc = MRAMc - 1d-6*sizeof(LTElines(ia)%damp)
          deallocate(LTElines(ia)%damp)
        end if

        ! broad_args
        if (allocated(LTElines(ia)%broad_args)) then
          MRAMc = MRAMc - 1d-6*sizeof(LTElines(ia)%broad_args)
          deallocate(LTElines(ia)%broad_args)
        end if

        ! nl
        if (allocated(LTElines(ia)%nl)) then
          MRAMc = MRAMc - 1d-6*sizeof(LTElines(ia)%nl)
          deallocate(LTElines(ia)%nl)
        end if

        ! nu
        if (allocated(LTElines(ia)%nu)) then
          MRAMc = MRAMc - 1d-6*sizeof(LTElines(ia)%nu)
          deallocate(LTElines(ia)%nu)
        end if

        ! n
        if (allocated(LTElines(ia)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(LTElines(ia)%n)
          deallocate(LTElines(ia)%n)
        end if

        ! Non-arrays
        MRAMc = MRAMc - 1d-6*sizeof(LTElines(ia))

      end do ! LTE lines

      ! Free
      deallocate(LTElines)

      end subroutine free_LTElines_full

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate completely the data in a Mol_class array\n
      !!  Mol(Mol_class(:)): Structures with molecular data
      subroutine free_mol_full(Mol)

      ! I/O

      type(Mol_class), dimension(:), allocatable, intent(inout):: Mol

      ! Local

      integer:: ia,ii


      ! Return if empty structure
      if (.not.allocated(Mol)) return

      ! For each molecule
      do ia=1,size(Mol)

        ! catom
        if (allocated(Mol(ia)%catom)) then

          ! For each element
          do ii=lbound(Mol(ia)%catom,1), &
                ubound(Mol(ia)%catom,1)

            ! String
            MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%catom(ii))

          end do ! Elements

          ! Array
          deallocate(Mol(ia)%catom)

        end if ! catom

        !
        ! Normal arrays
        !

        !
        ! Molecule
        if (allocated(Mol(ia)%Molecule)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%Molecule)
          deallocate(Mol(ia)%Molecule)
        end if

        !
        ! natom
        if (allocated(Mol(ia)%natom)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%natom)
          deallocate(Mol(ia)%natom)
        end if

        !
        ! iatom
        if (allocated(Mol(ia)%iatom)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%iatom)
          deallocate(Mol(ia)%iatom)
        end if

        !
        ! pfcoeff
        if (allocated(Mol(ia)%pfcoeff)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%pfcoeff)
          deallocate(Mol(ia)%pfcoeff)
        end if

        !
        ! pf
        if (allocated(Mol(ia)%pf)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%pf)
          deallocate(Mol(ia)%pf)
        end if

        !
        ! eqcoeff
        if (allocated(Mol(ia)%eqcoeff)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%eqcoeff)
          deallocate(Mol(ia)%eqcoeff)
        end if

        !
        ! eq
        if (allocated(Mol(ia)%eq)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%eq)
          deallocate(Mol(ia)%eq)
        end if

        !
        ! n
        if (allocated(Mol(ia)%n)) then
          MRAMc = MRAMc - 1d-6*sizeof(Mol(ia)%n)
          deallocate(Mol(ia)%n)
        end if

        ! Non-arrays
        MRAMc = MRAMc - 1d-6*sizeof(Mol(ia))

      end do ! Molecules

      ! Deallocate array
      deallocate(Mol)

      end subroutine free_mol_full

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the solution data in Solution_F_class used in the
      !! inversion\n
      !!   Sol(Solution_F_class): Structure with the solution of the
      !!                          self-consistent problem and the
      !!                          corresponding emergent profiles,
      !!                          contribution function, and height
      !!                          for optical depth equal to one
      subroutine free_inv_solution(Sol)

      ! I/O

      type(Solution_F_class), intent(inout):: Sol

      ! Local

      integer:: ii


      ! Slaves, leave, because you do not allocate these
      if (pid.gt.0) return

      !
      ! Free all allocated
      !

      !
      ! i_StkI
      if (allocated(Sol%i_StkI)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_StkI)
        deallocate(Sol%i_StkI)
      end if

      !
      ! i_StkI_b
      if (allocated(Sol%i_StkI_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_StkI_b)
        deallocate(Sol%i_StkI_b)
      end if

      !
      ! i_StkI_t
      if (allocated(Sol%i_StkI_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_StkI_t)
        deallocate(Sol%i_StkI_t)
      end if

      !
      ! i_J00
      if (allocated(Sol%i_J00)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00)
        deallocate(Sol%i_J00)
      end if

      !
      ! i_J00C
      if (allocated(Sol%i_J00C)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00C)
        deallocate(Sol%i_J00C)
      end if

      !
      ! i_J00P
      if (allocated(Sol%i_J00P)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00P)
        deallocate(Sol%i_J00P)
      end if

      !
      ! i_J00_b
      if (allocated(Sol%i_J00_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00_b)
        deallocate(Sol%i_J00_b)
      end if

      !
      ! i_J00C_b
      if (allocated(Sol%i_J00C_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00C_b)
        deallocate(Sol%i_J00C_b)
      end if

      !
      ! i_J00P_b
      if (allocated(Sol%i_J00P_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00P_b)
        deallocate(Sol%i_J00P_b)
      end if

      !
      ! i_J00_t
      if (allocated(Sol%i_J00_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00_t)
        deallocate(Sol%i_J00_t)
      end if

      !
      ! i_J00C_t
      if (allocated(Sol%i_J00C_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00C_t)
        deallocate(Sol%i_J00C_t)
      end if

      !
      ! i_J00P_t
      if (allocated(Sol%i_J00P_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_J00P_t)
        deallocate(Sol%i_J00P_t)
      end if

      !
      ! i_rhoes
      if (allocated(Sol%i_rhoes)) then
        do ii=lbound(Sol%i_rhoes,1),ubound(Sol%i_rhoes,1)
          if (allocated(Sol%i_rhoes(ii)%rho)) then
            SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes(ii)%rho)
            deallocate(Sol%i_rhoes(ii)%rho)
          end if
          if (allocated(Sol%i_rhoes(ii)%crho)) then
            SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes(ii)%crho)
            deallocate(Sol%i_rhoes(ii)%crho)
          end if
        end do
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes)
        deallocate(Sol%i_rhoes)
      end if

      !
      ! i_rhoes_b
      if (allocated(Sol%i_rhoes_b)) then
        do ii=lbound(Sol%i_rhoes_b,1),ubound(Sol%i_rhoes_b,1)
          if (allocated(Sol%i_rhoes_b(ii)%rho)) then
            SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes_b(ii)%rho)
            deallocate(Sol%i_rhoes_b(ii)%rho)
          end if
          if (allocated(Sol%i_rhoes_b(ii)%crho)) then
            SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes_b(ii)%crho)
            deallocate(Sol%i_rhoes_b(ii)%crho)
          end if
        end do
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes_b)
        deallocate(Sol%i_rhoes_b)
      end if

      !
      ! i_rhoes_t
      if (allocated(Sol%i_rhoes_t)) then
        do ii=lbound(Sol%i_rhoes_t,1),ubound(Sol%i_rhoes_t,1)
          if (allocated(Sol%i_rhoes_t(ii)%rho)) then
            SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes_t(ii)%rho)
            deallocate(Sol%i_rhoes_t(ii)%rho)
          end if
          if (allocated(Sol%i_rhoes_t(ii)%crho)) then
            SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes_t(ii)%crho)
            deallocate(Sol%i_rhoes_t(ii)%crho)
          end if
        end do
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_rhoes_t)
        deallocate(Sol%i_rhoes_t)
      end if

      !
      ! i_Stk
      if (allocated(Sol%i_Stk)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_Stk)
        deallocate(Sol%i_Stk)
      end if

      !
      ! i_Stk_b
      if (allocated(Sol%i_Stk_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_Stk_b)
        deallocate(Sol%i_Stk_b)
      end if

      !
      ! i_Stk_t
      if (allocated(Sol%i_Stk_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_Stk_t)
        deallocate(Sol%i_Stk_t)
      end if

      !
      ! i_JKQ
      if (allocated(Sol%i_JKQ)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQ)
        deallocate(Sol%i_JKQ)
      end if

      !
      ! i_JKQS
      if (allocated(Sol%i_JKQS)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQS)
        deallocate(Sol%i_JKQS)
      end if

      !
      ! i_JKQC
      if (allocated(Sol%i_JKQC)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQC)
        deallocate(Sol%i_JKQC)
      end if

      !
      ! i_JKQ_b
      if (allocated(Sol%i_JKQ_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQ_b)
        deallocate(Sol%i_JKQ_b)
      end if

      !
      ! i_JKQS_b
      if (allocated(Sol%i_JKQS_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQS_b)
        deallocate(Sol%i_JKQS_b)
      end if

      !
      ! i_JKQC_b
      if (allocated(Sol%i_JKQC_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQC_b)
        deallocate(Sol%i_JKQC_b)
      end if

      !
      ! i_JKQ_t
      if (allocated(Sol%i_JKQ_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQ_t)
        deallocate(Sol%i_JKQ_t)
      end if

      !
      ! i_JKQC_t
      if (allocated(Sol%i_JKQC_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQC_t)
        deallocate(Sol%i_JKQC_t)
      end if

      !
      ! i_JKQS_t
      if (allocated(Sol%i_JKQS_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%i_JKQS_t)
        deallocate(Sol%i_JKQS_t)
      end if

      !
      ! e_Stk
      if (allocated(Sol%e_Stk)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_Stk)
        deallocate(Sol%e_Stk)
      end if

      !
      ! e_tau1
      if (allocated(Sol%e_tau1)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_tau1)
        deallocate(Sol%e_tau1)
      end if

      !
      ! e_Ctr
      if (allocated(Sol%e_Ctr)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_Ctr)
        deallocate(Sol%e_Ctr)
      end if

      !
      ! e_Stk_b
      if (allocated(Sol%e_Stk_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_Stk_b)
        deallocate(Sol%e_Stk_b)
      end if

      !
      ! e_tau1_b
      if (allocated(Sol%e_tau1_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_tau1_b)
        deallocate(Sol%e_tau1_b)
      end if

      !
      ! e_Ctr_b
      if (allocated(Sol%e_Ctr_b)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_Ctr_b)
        deallocate(Sol%e_Ctr_b)
      end if

      !
      ! e_Stk_t
      if (allocated(Sol%e_Stk_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_Stk_t)
        deallocate(Sol%e_Stk_t)
      end if

      !
      ! e_tau1_t
      if (allocated(Sol%e_tau1_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_tau1_t)
        deallocate(Sol%e_tau1_t)
      end if

      !
      ! e_Ctr_t
      if (allocated(Sol%e_Ctr_t)) then
        SRAMc = SRAMc - 1d-6*sizeof(Sol%e_Ctr_t)
        deallocate(Sol%e_Ctr_t)
      end if

      end subroutine free_inv_solution

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the redistribution function\n
      !!  Red(Red_class): Structure with redistribution input
      !!                  frequency data, redistribution function
      !!                  data, and profile or normalization data
      subroutine free_warr(Red)

      ! I/O

      type(Red_class), intent(inout):: Red

      ! Local

      integer:: indx,jndx

      !
      ! If thete is data
      if (associated(Red%rzao)) then

        ! Run over indexes
        do indx=1,Red%nzao

          ! If input transition data
          if (associated(Red%rzao(indx)%trani)) then

            ! Run over input transitions
            do jndx=1,size(Red%rzao(indx)%trani)

              ! If allocated
              if (allocated(Red%rzao(indx)%trani(jndx)%iPPRD)) &
                deallocate(Red%rzao(indx)%trani(jndx)%iPPRD)
              if (allocated(Red%rzao(indx)%trani(jndx)%IWarr2)) &
                deallocate(Red%rzao(indx)%trani(jndx)%IWarr2)
              if (allocated(Red%rzao(indx)%trani(jndx)%PWarr2)) &
                deallocate(Red%rzao(indx)%trani(jndx)%PWarr2)

            end do ! Input transitions

            ! Free trani
            deallocate(Red%rzao(indx)%trani)
            nullify(Red%rzao(indx)%trani)

          end if ! There are input transitions

        end do ! Indexes

        ! Free rzao
        deallocate(Red%rzao)
        nullify(Red%rzao)

      end if ! There is data

      ! Reset RAM counter
      WRAMc = 0d0

      end subroutine free_warr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate arrays to store the emissivity in PRD lines\n
      !!  Red(Red_class): Structure with redistribution input
      !!                  frequency data, redistribution function
      !!                  data, and profile or normalization data
      subroutine free_e2ord(Red)

      ! I/O

      type(Red_class), intent(inout):: Red

      ! Local

      integer:: indx


      !
      ! If thete is data
      if (associated(Red%zao)) then

        ! Run over indexes
        do indx=1,Red%nzao

          ! Free arrays
          if (allocated(Red%zao(indx)%eps20)) &
            deallocate(Red%zao(indx)%eps20)
          if (allocated(Red%zao(indx)%eps21)) &
            deallocate(Red%zao(indx)%eps21)
          if (allocated(Red%zao(indx)%eps22)) &
            deallocate(Red%zao(indx)%eps22)
          if (allocated(Red%zao(indx)%eps23)) &
            deallocate(Red%zao(indx)%eps23)
          if (allocated(Red%zao(indx)%rpf)) &
            deallocate(Red%zao(indx)%rpf)

        end do ! Indexes

      end if ! There is data

      end subroutine free_e2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the input frequency axis for PRD calculations\n
      !!  Red(Red_class): Structure with redistribution input
      !!                  frequency data, redistribution function
      !!                  data, and profile or normalization data
      subroutine free_ifreq(Red)

      ! I/O

      type(Red_class), intent(inout):: Red

      ! Local

      integer:: indx,jndx


      !
      ! If thete is data
      if (associated(Red%zao)) then

        ! Run over indexes
        do indx=1,Red%nzao

          ! Free arrays
          if (allocated(Red%zao(indx)%eps20)) &
            deallocate(Red%zao(indx)%eps20)
          if (allocated(Red%zao(indx)%eps21)) &
            deallocate(Red%zao(indx)%eps21)
          if (allocated(Red%zao(indx)%eps22)) &
            deallocate(Red%zao(indx)%eps22)
          if (allocated(Red%zao(indx)%eps23)) &
            deallocate(Red%zao(indx)%eps23)
          if (allocated(Red%zao(indx)%rpf)) &
            deallocate(Red%zao(indx)%rpf)

          ! If input transition data
          if (associated(Red%zao(indx)%trani)) then

            ! Run over input transitions
            do jndx=1,size(Red%zao(indx)%trani)

              ! If allocated
              if (allocated(Red%zao(indx)%trani(jndx)%mfreq)) &
                deallocate(Red%zao(indx)%trani(jndx)%mfreq)
              if (allocated(Red%zao(indx)%trani(jndx)%omega)) &
                deallocate(Red%zao(indx)%trani(jndx)%omega)
              if (allocated(Red%zao(indx)%trani(jndx)%W_freq)) &
                deallocate(Red%zao(indx)%trani(jndx)%W_freq)

            end do ! Input transitions

            ! Free trani
            deallocate(Red%zao(indx)%trani)
            nullify(Red%zao(indx)%trani)

          end if ! There is input transition data

        end do ! Indexes

        ! Free zao
        deallocate(Red%zao)
        nullify(Red%zao)

      end if

      !
      ! If thete is data
      if (associated(Red%ao)) then

        ! Run over indexes
        do indx=1,Red%nao

          ! Free arrays
          if (allocated(Red%ao(indx)%nn)) &
            deallocate(Red%ao(indx)%nn)
          if (allocated(Red%ao(indx)%Mi0)) &
            deallocate(Red%ao(indx)%Mi0)
          if (allocated(Red%ao(indx)%Mi1)) &
            deallocate(Red%ao(indx)%Mi1)
          if (allocated(Red%ao(indx)%if0)) &
            deallocate(Red%ao(indx)%if0)
          if (allocated(Red%ao(indx)%if1)) &
            deallocate(Red%ao(indx)%if1)

        end do ! Indexes

        ! Free zao
        deallocate(Red%ao)
        nullify(Red%ao)

      end if

      ! Free izao
      if (allocated(Red%izao)) deallocate(Red%izao)

      ! Reset RAM counter
      FRAMc = 0d0

      end subroutine free_ifreq

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the profile normalization data\n
      !!  Red(Red_class): Structure with redistribution input
      !!                  frequency data, redistribution function
      !!                  data, and profile or normalization data\n
      !!   full(logical): Indicate if the whole structure needs to
      !!                  be freed
      subroutine free_norm(Red,full)

      ! I/O

      type(Red_class), intent(inout):: Red
      logical, intent(in):: full

      ! Local

      integer:: indx


      !
      ! If thete is data
      if (associated(Red%dzao)) then

        ! Run over indexes
        do indx=1,Red%ndzao

          ! If allocated
          if (allocated(Red%dzao(indx)%Norm)) &
            deallocate(Red%dzao(indx)%Norm)
          if (allocated(Red%dzao(indx)%p)) &
            deallocate(Red%dzao(indx)%p)
          if (allocated(Red%dzao(indx)%cp)) &
            deallocate(Red%dzao(indx)%cp)
        end do ! Indexes

        ! Free dzao
        if (full) then
          deallocate(Red%dzao)
          nullify(Red%dzao)
        end if

      end if

      ! Free idzao indexing
      if (full) then
        if (allocated(Red%idzao)) deallocate(Red%idzao)
      end if

      ! Reset RAM counter
      VRAMc = 0d0

      end subroutine free_norm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate the profile normalization data used to compute the
      !! PRD emissivity\n
      !!  Red(Red_class): Structure with redistribution input
      !!                  frequency data, redistribution function
      !!                  data, and profile or normalization data
      subroutine free_1stord(Red)

      ! I/O

      type(Red_class), intent(inout):: Red

      ! Local

      integer:: indx


      !
      ! If thete is data
      if (associated(Red%pzao)) then

        ! Run over indexes
        do indx=1,Red%nzao

          ! If allocated
          if (allocated(Red%pzao(indx)%Norm)) &
            deallocate(Red%pzao(indx)%Norm)
          if (allocated(Red%pzao(indx)%p)) &
            deallocate(Red%pzao(indx)%p)
          if (allocated(Red%pzao(indx)%cp)) &
            deallocate(Red%pzao(indx)%cp)
        end do ! Indexes

        ! Free dzao
        deallocate(Red%pzao)
        nullify(Red%pzao)

      end if ! There is data

      ! Reset RAM counter
      ORAMc = 0d0

      end subroutine free_1stord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate all data in the Red_class structure\n
      !!  Red(Red_class): Structure with redistribution input
      !!                  frequency data, redistribution function
      !!                  data, and profile or normalization data
      subroutine free_red(Red)

      ! I/O

      type(Red_class), intent(inout):: Red

      ! Free CRD normalization
      call free_norm(Red,.True.)

      ! Free PRD variables
      call free_warr(Red)
      call free_1stord(Red)
      call free_ifreq(Red)

      end subroutine free_red

!#####################################################################
!#####################################################################
!#####################################################################

      end module free_mod
