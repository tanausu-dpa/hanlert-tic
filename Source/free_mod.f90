      !> Frees allocated memory
      module free_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Contributors:
!     Hao Li (IAC)
!  Start:
!     06/28/2022
!  Last version:
!     10/16/2023 V3.1.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/16/2023:    V3.1.3 - Moved the deallocation of the damping
!                             from free_local_Atom to the new
!                             free_damp routine (TdPA)
!                           - Added free_local_geom routine (TdPA)
!                           - Added free_damp routine. This one is
!                             called from free_pix (TdPA)
!
!     10/04/2023:    V3.1.2 - When freeing the profiles for LTE lines,
!                             get the range limits from the array
!                             metadata, as the variable could be a
!                             dummy (TdPA)
!
!     08/28/2023:    V3.1.1 - Add deallocation of LTE profiles (TdPA)
!
!     07/03/2023:    V3.1.0 - Split deallocations in different
!                             subroutines to allow for more
!                             flexibility to call them (TdPA)
!                           - Improved the deallocation of the
!                             model atmosphere (TdPA)
!
!     02/14/2023:    V3.0.4 - Do not deallocate atmosphere and
!                             magnetic field if doing inversion (TdPA)
!
!     11/24/2022:    V3.0.3 - Finished free_cle_node first working
!                             version (TdPA)
!
!     10/25/2022:    V3.0.2 - Updated free_ram header (TdPA)
!                           - Bugfix: Correctly free the memory in
!                             Atom%Normp (TdPA)
!                           - Bugfix: Deallocate Atom%NCHLT is
!                             present (TdPA)
!                           - Bugfix: Ensure atmospheric variables
!                             (secondary ones) are always deallocated
!                             if still present (TdPA)
!                           - Bugfix: Check if Bfield is allocated
!                             before freeing it (TdPA)
!                           - Added free_cle_node subroutine for the
!                             CLE problem. It is not finished (TdPA)
!
!     07/13/2022:    V3.0.1 - Added deallocation of Atmo%nha and
!                             Atmo%zalt (TdPA)
!
!     06/29/2022:    V3.0.0 - First version (TdPA)
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
!  free_local_Atom:
!    Deallocate atomic quantities generated during a call to hanle
!
!  free_local_geom:
!    Deallocate geometry generated during a call to hanle
!
!  free_local_CGF:
!    Deallocate continuum, geometry, and frequency quantities
!  generated during a call to hanle
!
!  free_local_LTE:
!    Deallocate profiles in LTE lines
!
!  free_local:
!    Calls free_local_Atom and free_local_CGF
!
!  free_Atmo:
!    Deallocate an atmospheric model partially or completely
!
!  free_gpop:
!    Deallocate total populations of atoms and molecules
!
!  free_lpop:
!    Deallocate (LTE and NLTE) populations and density matrices
!  of the atoms
!
!  free_mol:
!    Deallocate molecular quantities and population arrays
!
!  free_cols:
!    Deallocate collisions for active atoms
!
!  free_damp:
!    Deallocate damping coefficients for list of atoms
!
!  free_B:
!    Deallocate magnetic field arrays
!
!  free_pixel:
!    Calls free_gpop, free_lpop, free_mol, free_cols, and free_B
!
!  free_cle_node:
!    Deallocate memory that is expected to be allocated again when
!  running the CLE case
!
!  free_inv_solution:
!    Release the memory for the synthesis solutions in the inversion
!  process
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

      !> Deallocate memory that depends on the synthesis run\n
      !!      Atom(Atom_class): Atom structures (active)
      subroutine free_local_Atom(Atom)

      ! I/O
      type(Atom_class), dimension(:):: Atom

      ! Local
      integer:: ii,jj,jdir,jtran,iz

      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Eigenvalues/vectors
        if (allocated(Atom(ii)%eval)) deallocate(Atom(ii)%eval)
        if (allocated(Atom(ii)%evec)) deallocate(Atom(ii)%evec)

        ! Damping

        ! Photoionization
        do jj=1,Atom(ii)%nphot
          if (allocated(Atom(ii)%phot(jj)%TEI)) &
            deallocate(Atom(ii)%phot(jj)%TEI)
        end do

        ! Norm
        if (associated(Atom(ii)%Normp)) then
          do jdir=1,size(Atom(ii)%Normp,3)
            do iz=Rz0,Rz1
              do jtran=1,Atom(ii)%ntran
                if (allocated(Atom(ii)% &
                              Normp(jtran,iz,jdir)%prof)) then
                  deallocate(Atom(ii)%Normp(jtran,iz,jdir)%prof)
                end if
                if (allocated(Atom(ii)% &
                                   Normp(jtran,iz,jdir)%Norm)) then
                  deallocate(Atom(ii)%Normp(jtran,iz,jdir)%Norm)
                end if
              end do
            end do
          end do
          deallocate(Atom(ii)%Normp)
          nullify(Atom(ii)%Normp)
        end if

        ! NCHLT
        if (allocated(Atom(ii)%NCHLT)) deallocate(Atom(ii)%NCHLT)

      end do ! Atoms

      return

      end subroutine free_local_Atom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory that depends on the synthesis run\n
      !!  Geom(Geometry_class): Structure with geometry data\n
      subroutine free_local_geom(Geom)

      ! I/O
      type(Geometry_class):: Geom


      !
      ! Geometry
      !
      if (allocated(Geom%TB)) deallocate(Geom%TB)
      if (allocated(Geom%TBL)) deallocate(Geom%TBL)

      return

      end subroutine free_local_geom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory that depends on the synthesis run\n
      !! Cont(Continuum_class): Structure with background opacity
      !!  Geom(Geometry_class): Structure with geometry data\n
      !! Frec(Frequency_class): Structure with frequency data
      subroutine free_local_CGF(Cont,Geom,Frec)

      ! I/O
      type(Continuum_class):: Cont
      type(Geometry_class):: Geom
      type(Frequency_class):: Frec


      !
      ! Continuum
      !
      if (allocated(Cont%c)) deallocate(Cont%c)

      !
      ! Geometry
      !
      if (allocated(Geom%TB)) deallocate(Geom%TB)
      if (allocated(Geom%TBL)) deallocate(Geom%TBL)

      !
      ! Frequency
      !
      if (allocated(Frec%stype)) deallocate(Frec%stype)
      if (allocated(Frec%omega3)) deallocate(Frec%omega3)
      if (associated(Frec%exu)) then
        deallocate(Frec%exu)
        nullify(Frec%exu)
      end if

      return

      end subroutine free_local_CGF

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory that depends on the synthesis run\n
      !!   lines(LTEline_class): Structure with the LTE line data
      subroutine free_local_LTE(lines)

      ! I/O
      type(LTEline_class), dimension(:), allocatable:: lines

      ! Local
      integer:: ii,jdir,iz


      ! Allocated?
      if (.not.allocated(lines)) return

      ! For each line
      do ii=1,size(lines)

        ! Check prof is allocated
        if (associated(lines(ii)%prof)) then
          do jdir=1,size(lines(ii)%prof,2)
            do iz=lbound(lines(ii)%prof,1),ubound(lines(ii)%prof,1)
              if (allocated(lines(ii)%prof(iz,jdir)%p)) &
                deallocate(lines(ii)%prof(iz,jdir)%p)
              if (allocated(lines(ii)%prof(iz,jdir)%comp)) &
                deallocate(lines(ii)%prof(iz,jdir)%comp)
            end do
          end do
          deallocate(lines(ii)%prof)
          nullify(lines(ii)%prof)
        end if

      end do

      end subroutine free_local_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory that depends on the synthesis run\n
      !!        Atom(Atom_class): Atom structures (active)\n
      !!   Cont(Continuum_class): Structure with background opacity
      !!    Geom(Geometry_class): Structure with geometry data\n
      !!   Frec(Frequency_class): Structure with frequency data\n
      !! LTElines(LTEline_class): Structure with the LTE line data
      subroutine free_local(Atom,Cont,Geom,Frec,LTElines)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Continuum_class):: Cont
      type(Geometry_class):: Geom
      type(Frequency_class):: Frec

      call free_local_Atom(Atom)
      call free_local_CGF(Cont,Geom,Frec)
      call free_local_LTE(LTElines)

      return

      end subroutine free_local

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate a copy of atmospheric model\n
      !!    Atmo(Atmo_class): Atmosphere structure\n
      !!       full(logical): If deallocate everything
      subroutine free_Atmo(Atmo,full)

      ! I/O
      type(Atmo_class):: Atmo
      logical, intent(in):: full

      !
      ! Pointers
      !

      ! If allocated, group a
      if (Atmo%alloc_a) then

        if (associated(Atmo%z)) then
          deallocate(Atmo%z)
          nullify(Atmo%z)
        end if
        if (associated(Atmo%T)) then
          deallocate(Atmo%T)
          nullify(Atmo%T)
        end if
        if (associated(Atmo%vmi)) then
          deallocate(Atmo%vmi)
          nullify(Atmo%vmi)
        end if
        if (associated(Atmo%vx)) then
          deallocate(Atmo%vx)
          nullify(Atmo%vx)
        end if
        if (associated(Atmo%vy)) then
          deallocate(Atmo%vy)
          nullify(Atmo%vy)
        end if
        if (associated(Atmo%vz)) then
          deallocate(Atmo%vz)
          nullify(Atmo%vz)
        end if

      ! If pointing
      else

        if (associated(Atmo%z)) nullify(Atmo%z)
        if (associated(Atmo%T)) nullify(Atmo%T)
        if (associated(Atmo%vmi)) nullify(Atmo%vmi)
        if (associated(Atmo%vx)) nullify(Atmo%vx)
        if (associated(Atmo%vy)) nullify(Atmo%vy)
        if (associated(Atmo%vz)) nullify(Atmo%vz)

      end if

      !
      ! If allocated, group b
      if (Atmo%alloc_b) then

        if (associated(Atmo%Bx)) then
          deallocate(Atmo%Bx)
          nullify(Atmo%Bx)
        end if
        if (associated(Atmo%By)) then
          deallocate(Atmo%By)
          nullify(Atmo%By)
        end if
        if (associated(Atmo%Bz)) then
          deallocate(Atmo%Bz)
          nullify(Atmo%Bz)
        end if

      else

        if (associated(Atmo%Bx)) nullify(Atmo%Bx)
        if (associated(Atmo%By)) nullify(Atmo%By)
        if (associated(Atmo%Bz)) nullify(Atmo%Bz)

      end if

      !
      ! Always pointers
      if (associated(Atmo%vxa)) then
        nullify(Atmo%vxa)
      end if
      if (associated(Atmo%vya)) then
        nullify(Atmo%vya)
      end if
      if (associated(Atmo%vza)) then
        nullify(Atmo%vza)
      end if

      ! Always arrays
      if (associated(Atmo%zeros)) then
        deallocate(Atmo%zeros)
        nullify(Atmo%zeros)
      end if

      !
      ! Arrays!
      if (allocated(Atmo%nHT)) deallocate(Atmo%nHT)
      if (allocated(Atmo%nHm)) deallocate(Atmo%nHm)
      if (allocated(Atmo%Pg)) deallocate(Atmo%Pg)
      if (allocated(Atmo%rho)) deallocate(Atmo%rho)
      if (allocated(Atmo%Pe)) deallocate(Atmo%Pe)
      if (allocated(Atmo%zalt)) deallocate(Atmo%zalt)
      if (allocated(Atmo%nHa)) deallocate(Atmo%nHa)
      if (allocated(Atmo%chi500)) deallocate(Atmo%chi500)
      if (allocated(Atmo%ne)) deallocate(Atmo%ne)
      if (allocated(Atmo%nh)) deallocate(Atmo%nh)
      if (allocated(Atmo%nhe)) deallocate(Atmo%nhe)
      if (allocated(Atmo%vlos)) deallocate(Atmo%vlos)
      if (allocated(Atmo%vpos)) deallocate(Atmo%vpos)
      if (allocated(Atmo%vphi)) deallocate(Atmo%vphi)

      ! Deallocate full
      if (full) then
        if (allocated(Atmo%pT)) deallocate(Atmo%pT)
        if (allocated(Atmo%ele)) deallocate(Atmo%ele)
        if (allocated(Atmo%abund)) deallocate(Atmo%abund)
        if (allocated(Atmo%JKQin)) deallocate(Atmo%JKQin)
      end if

      return

      end subroutine free_Atmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate global populations\n
      !!      Atom(Atom_class): Atom structures (active)\n
      !!     Atomb(Atom_class): Atom structures (passive)\n
      !!        Mol(Mol_class): Molecule structures
      subroutine free_gpop(Atom,Atomb,Mol)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:):: Atomb
      type(Mol_class), dimension(:):: Mol

      ! Local
      integer:: ii

      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Population
        if (allocated(Atom(ii)%n)) deallocate(Atom(ii)%n)

      end do ! Atoms

      ! For each passive atom
      do ii=1,nAb

        ! Population
        if (allocated(Atomb(ii)%n)) deallocate(Atomb(ii)%n)

      end do ! Passive atoms

      !
      ! Molecules
      !

      ! For each molecule
      do ii=1,nM

        ! Population
        if (allocated(Mol(ii)%n)) deallocate(Mol(ii)%n)

      end do

      return

      end subroutine free_gpop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate level populations\n
      !!      Atom(Atom_class): Atom structures (active)\n
      !!     Atomb(Atom_class): Atom structures (passive)
      subroutine free_lpop(Atom,Atomb)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:):: Atomb

      ! Local
      integer:: ii

      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Population
        if (allocated(Atom(ii)%popu)) deallocate(Atom(ii)%popu)
        if (allocated(Atom(ii)%populte)) deallocate(Atom(ii)%populte)

        ! Rho
        if (allocated(Atom(ii)%crho)) deallocate(Atom(ii)%crho)
        if (allocated(Atom(ii)%rhonull)) deallocate(Atom(ii)%rhonull)

      end do ! Atoms

      ! For each passive atom
      do ii=1,nAb

        ! Population
        if (allocated(Atomb(ii)%popu)) deallocate(Atomb(ii)%popu)
        if (allocated(Atomb(ii)%populte)) &
          deallocate(Atomb(ii)%populte)

      end do ! Passive atoms

      return

      end subroutine free_lpop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate molecule arrays\n
      !!    Mol(Mol_class): Molecule structures
      subroutine free_mol(Mol)

      ! I/O
      type(Mol_class), dimension(:):: Mol

      ! Local
      integer:: ii

      ! For each molecule
      do ii=1,nM

        ! Population
        if (allocated(Mol(ii)%n)) deallocate(Mol(ii)%n)
        if (allocated(Mol(ii)%pf)) deallocate(Mol(ii)%pf)
        if (allocated(Mol(ii)%eq)) deallocate(Mol(ii)%eq)

      end do

      return

      end subroutine free_mol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Free data on collisions\n
      !!      Atom(Atom_class): Atom structures (active)
      subroutine free_cols(Atom)

      ! I/O
      type(Atom_class), dimension(:):: Atom

      ! Local
      integer:: ii

      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Collisions
        if (allocated(Atom(ii)%Ccoeff)) deallocate(Atom(ii)%Ccoeff)
        if (allocated(Atom(ii)%CcoeffJ)) deallocate(Atom(ii)%CcoeffJ)
        if (allocated(Atom(ii)%gk)) deallocate(Atom(ii)%gk)
        if (allocated(Atom(ii)%icol)) deallocate(Atom(ii)%icol)

      end do ! Atoms

      end subroutine free_cols

!#####################################################################
!#####################################################################
!#####################################################################

      !> Free data on damping\n
      !!    Atom(Atom_class): Atom structures\n
      !!         nn(integer): Size of Atom array
      subroutine free_damp(Atom,nn)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      integer, intent(in):: nn

      ! Local
      integer:: ii

      !
      ! Atom
      !

      ! For each atom
      do ii=1,nn

        if (allocated(Atom(ii)%damp)) deallocate(Atom(ii)%damp)
        if (allocated(Atom(ii)%ldamp)) deallocate(Atom(ii)%ldamp)

      end do ! Atoms

      end subroutine free_damp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory for magnetic field\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine free_B(Bfield)

      ! I/O
      type(Bfield_class):: Bfield

      ! Magnetic field
      if (allocated(Bfield%Bstrength)) &
        deallocate(Bfield%Bstrength,Bfield%Btheta,Bfield%Bphi)
      if (allocated(Bfield%Blos)) &
        deallocate(Bfield%Blos,Bfield%Bpos,Bfield%Azimuth)

      return

      end subroutine free_B

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory that depends on the pixel\n
      !!      Atom(Atom_class): Atom structures (active)\n
      !!     Atomb(Atom_class): Atom structures (passive)\n
      !!        Mol(Mol_class): Molecule structures\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine free_pix(Atom,Atomb,Mol,Bfield)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:):: Atomb
      type(Mol_class), dimension(:):: Mol
      type(Bfield_class):: Bfield


      ! Free global populations
      call free_gpop(Atom,Atomb,Mol)

      ! Free local populations
      call free_lpop(Atom,Atomb)

      ! Free molecule
      call free_mol(Mol)

      ! Free collisions
      call free_cols(Atom)

      ! Free damping
      call free_damp(Atom,nA)
      if (nAb.gt.0) call free_damp(Atomb,nAb)

      ! Free Magnetic field
      call free_B(Bfield)

      return

      end subroutine free_pix

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate memory at node level in CLE\n
      !!      Atom(Atom_class): Atom structures (active)\n
      !!     Atomb(Atom_class): Atom structures (passive)\n
      !!        Mol(Mol_class): Molecule structures\n
      !!      Atmo(Atmo_class): Atmosphere structure\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine free_cle_node(Atom,Atomb,Mol,Atmo,Bfield,Geom)

      ! I/O
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:):: Atomb
      type(Mol_class), dimension(:):: Mol
      type(Atmo_class):: Atmo
      type(Bfield_class):: Bfield
      type(Geometry_class):: Geom

      ! Local
      integer:: ii,jj

      !
      ! Atom
      !

      ! For each active atom
      do ii=1,nA

        ! Population
        if (allocated(Atom(ii)%n)) deallocate(Atom(ii)%n)
        if (allocated(Atom(ii)%popu)) deallocate(Atom(ii)%popu)
        if (allocated(Atom(ii)%populte)) deallocate(Atom(ii)%populte)

        ! Eigenvalues/vectors
        if (allocated(Atom(ii)%eval)) deallocate(Atom(ii)%eval)
        if (allocated(Atom(ii)%evec)) deallocate(Atom(ii)%evec)

        ! Rho
        if (allocated(Atom(ii)%crho)) deallocate(Atom(ii)%crho)
        if (allocated(Atom(ii)%rhonull)) deallocate(Atom(ii)%rhonull)

        ! Damping
        if (allocated(Atom(ii)%damp)) deallocate(Atom(ii)%damp)
        if (allocated(Atom(ii)%ldamp)) deallocate(Atom(ii)%ldamp)

        ! Collisions
        if (allocated(Atom(ii)%Ccoeff)) deallocate(Atom(ii)%Ccoeff)
        if (allocated(Atom(ii)%CcoeffJ)) deallocate(Atom(ii)%CcoeffJ)
        if (allocated(Atom(ii)%gk)) deallocate(Atom(ii)%gk)
        if (allocated(Atom(ii)%icol)) deallocate(Atom(ii)%icol)

        ! Photoionization
        do jj=1,Atom(ii)%nphot
          if (allocated(Atom(ii)%phot(jj)%TEI)) &
            deallocate(Atom(ii)%phot(jj)%TEI)
        end do

        ! NCHLT
        if (allocated(Atom(ii)%NCHLT)) deallocate(Atom(ii)%NCHLT)

      end do

      ! For each passive atom
      do ii=1,nAb

        ! Population
        if (allocated(Atomb(ii)%n)) deallocate(Atomb(ii)%n)
        if (allocated(Atomb(ii)%popu)) deallocate(Atomb(ii)%popu)
        if (allocated(Atomb(ii)%populte)) &
          deallocate(Atomb(ii)%populte)

        ! Damping
        if (allocated(Atomb(ii)%damp)) deallocate(Atomb(ii)%damp)
        if (allocated(Atomb(ii)%ldamp)) deallocate(Atomb(ii)%ldamp)

      end do ! Passive atoms

      !
      ! Molecules
      !

      ! For each molecule
      do ii=1,nM

        ! Population, partition, and eq. const.
        if (allocated(Mol(ii)%n)) deallocate(Mol(ii)%n)
        if (allocated(Mol(ii)%pf)) deallocate(Mol(ii)%pf)
        if (allocated(Mol(ii)%eq)) deallocate(Mol(ii)%eq)

      end do ! Molecules

      !
      ! Atmosphere
      !

      ! Pointers
      nullify(Atmo%T,Atmo%vmi,Atmo%vx,Atmo%vy,Atmo%vz)
      if (associated(Atmo%Bx)) nullify(Atmo%Bx,Atmo%By,Atmo%Bz)

      !
      ! Magnetic field
      !
      if (allocated(Bfield%Bstrength)) &
        deallocate(Bfield%Bstrength,Bfield%Btheta,Bfield%Bphi)

      !
      ! Geometrical tensors
      !
      if (allocated(Geom%TS)) deallocate(Geom%TS,Geom%TB)

      return

      end subroutine free_cle_node

!#####################################################################
!#####################################################################
!#####################################################################

      !> Deallocate solutions in inversion pixel\n
      !!  Sol(Solution_F_class): Class with the data of the RT
      subroutine free_inv_solution(Sol)

      ! I/O
      type(Solution_F_class), intent(inout):: Sol

      ! Slaves, leave
      if (pid.gt.0) return

      ! Free all allocated
      if (allocated(Sol%i_J00)) deallocate(Sol%i_J00)
      if (allocated(Sol%i_J00C)) deallocate(Sol%i_J00C)
      if (allocated(Sol%e_tau1)) deallocate(Sol%e_tau1)
      if (allocated(Sol%i_J00P)) deallocate(Sol%i_J00P)
      if (allocated(Sol%e_Stk)) deallocate(Sol%e_Stk)
      if (allocated(Sol%e_Ctr)) deallocate(Sol%e_Ctr)
      if (allocated(Sol%i_StkI)) deallocate(Sol%i_StkI)
      if (allocated(Sol%i_Stk)) deallocate(Sol%i_Stk)
      if (allocated(Sol%i_JKQ)) deallocate(Sol%i_JKQ)
      if (allocated(Sol%i_JKQS)) deallocate(Sol%i_JKQS)
      if (allocated(Sol%i_JKQC)) deallocate(Sol%i_JKQC)
      if (allocated(Sol%i_rhoes)) deallocate(Sol%i_rhoes)
      if (allocated(Sol%i_J00_b)) deallocate(Sol%i_J00_b)
      if (allocated(Sol%i_J00C_b)) deallocate(Sol%i_J00C_b)
      if (allocated(Sol%e_tau1_b)) deallocate(Sol%e_tau1_b)
      if (allocated(Sol%i_J00P_b)) deallocate(Sol%i_J00P_b)
      if (allocated(Sol%e_Stk_b)) deallocate(Sol%e_Stk_b)
      if (allocated(Sol%e_Ctr_b)) deallocate(Sol%e_Ctr_b)
      if (allocated(Sol%i_StkI_b)) deallocate(Sol%i_StkI_b)
      if (allocated(Sol%i_Stk_b)) deallocate(Sol%i_Stk_b)
      if (allocated(Sol%i_JKQ_b)) deallocate(Sol%i_JKQ_b)
      if (allocated(Sol%i_JKQS_b)) deallocate(Sol%i_JKQS_b)
      if (allocated(Sol%i_JKQC_b)) deallocate(Sol%i_JKQC_b)
      if (allocated(Sol%i_rhoes_b)) deallocate(Sol%i_rhoes_b)
      if (allocated(Sol%i_J00_t)) deallocate(Sol%i_J00_t)
      if (allocated(Sol%i_J00C_t)) deallocate(Sol%i_J00C_t)
      if (allocated(Sol%e_tau1_t)) deallocate(Sol%e_tau1_t)
      if (allocated(Sol%i_J00P_t)) deallocate(Sol%i_J00P_t)
      if (allocated(Sol%e_Stk_t)) deallocate(Sol%e_Stk_t)
      if (allocated(Sol%e_Ctr_t)) deallocate(Sol%e_Ctr_t)
      if (allocated(Sol%i_StkI_t)) deallocate(Sol%i_StkI_t)
      if (allocated(Sol%i_Stk_t)) deallocate(Sol%i_Stk_t)
      if (allocated(Sol%i_JKQ_t)) deallocate(Sol%i_JKQ_t)
      if (allocated(Sol%i_JKQC_t)) deallocate(Sol%i_JKQC_t)
      if (allocated(Sol%i_JKQS_t)) deallocate(Sol%i_JKQS_t)
      if (allocated(Sol%i_rhoes_t)) deallocate(Sol%i_rhoes_t)

      end subroutine free_inv_solution

!#####################################################################
!#####################################################################
!#####################################################################

      end module free_mod
