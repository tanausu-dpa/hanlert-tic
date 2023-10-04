      !> Get CLE RT coefficients
      module getrtcle_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     10/xx/2022
!  Last version:
!     09/29/2023 V3.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.3 - Added arguments to Initcols (TdPA)
!
!     08/07/2023:    V3.0.2 - Added dummy variable for compatibility
!                             with modules (TdPA)
!
!     07/03/2023:    V3.0.1 - Adjusted for changes elsewhere in the
!                             code, it is untested (TdPA)
!
!     11/24/2022:    V3.0.0 - First version (TdPA)
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
!    Return RT coefficients (and Stokes parameters for a boundary)
!  for a given CLE point
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use background_mod
      use boundary_mod
      use broad_mod
      use btens_mod
      use chemic_mod
      use diagon_mod
      use free_mod
      use geometry_mod
      use initcols_mod
      use initphotoion_mod
      use initpopu_mod
      use jcalccle_mod
      use parameters_mod , only: B2L , TINYB
      use ratmo_mod
      use ratom_mod
      use rmol_mod
      use rtcoeff_mod
      use see_mod
      use stens_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return the RT coefficients for a given point along the LOS
      !! for a CLE calculation
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!           Atmo(Atmo_class): Structure with atmospheric
      !!                             data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with settings data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!       Geom(Geometry_class): Structure with quadrature data\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!                ix(integer): Index location along the LOS\n
      !!               if0(integer): Lower limit index for frequencies
      !!                             in this CPU\n
      !!               if1(integer): Upper limit index for frequencies
      !!                             in this CPU\n
      !!           batmo(double(:)): Atmosphere data\n
      !!            bion(double(:)): Read ionization fraction data\n
      !! ion_column_ind(integer(:)): Index of column in buffer for
      !!                             ionization data\n
      !!  ion_value_ind(integer(:)): Index of value in value array for
      !!                             ionization data\n
      !!       ion_value(double(:)): Numeric constant ionization
      !!                             fraction values\n
      !!         spect(spect_class): Structure with the input spectra
      !!                             data\n
      !!     chianti(chianti_class): Structure with the CHIANTI data\n
      !!         data1(dfloat(:,:)): Radiation transfer coefficients\n
      !!            indisk(logical): If the ray crossed the disk\n
      !!          isbottom(logical): If this is the bottom boundary
      subroutine getRTCLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec, &
                          Geom,Flgsg,fudge,kurucz,ix,if0,if1,batmo, &
                          bion,ion_column_ind,ion_value_ind, &
                          ion_value,spect,chianti,data1,indisk, &
                          isbottom)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Atmo_class):: Atmo
      type(Fctsg_class):: Flgsg
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(Frequency_class):: Frec
      type(Input_class):: Input
      type(MPI_class):: MPID
      type(spect_class):: spect
      type(Geometry_class):: Geom
      type(chianti_class):: chianti
      logical, intent(in):: indisk,isbottom
      integer, intent(in):: ix,if0,if1
      integer, dimension(:), intent(in):: ion_value_ind,ion_column_ind
      double precision, dimension(0:3,if0:if1,0:5):: data1
      double precision, dimension(:), intent(in):: ion_value,bion
      double precision, dimension(:), intent(in):: batmo

      ! Local
      type(Continuum_class):: Cont
      type(Coronapoint_class):: GeomP
      type(Bfield_class):: Bfield
      type(Geometry_class):: GeomS
      type(LTEline_class), dimension(:), allocatable:: dummy

      integer:: ia,it0,it1,ip0,ip1
      integer, dimension(:), allocatable:: nlte,depar

      double precision:: larmor
      double precision, dimension(nxphot,2):: Jphot

      complex(kind=8), dimension(-2:2,0:2,nxtran):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2), target:: TS
      complex(kind=8), dimension(:,:,:), pointer:: TB

      !
      ! Initialize
      !
      data1 = 0d0

      !
      ! Get atmosphere
      !
      call rAtmo_cle(batmo,Input,Atmo,ix)


      !
      ! Prepare atomic and molecular quantities depending on
      ! the atmospheric model
      !

      ! Prepare active atoms
      call prepareatom(Atom,nA)
      ! Prepare passive atoms
      if (nAb.gt.0) call prepareatom(Atomb,nAb)

      ! Prepare molecules
      if (Input%nM.gt.0) call preparemol(Mol,Input%nM)

      ! Control
      if (laborted) goto 1000

      !
      ! Get the local frame quantities
      !
      call getlocalframe(Atmo,Bfield,Atmo%z(ix),Atmo%mode,GeomP)

      !
      ! Compute TKQ geometrical tensors
      !
      call Stens(GeomP%geom(1),GeomP%geom(2),GeomP%geom(3),Flgsg,TS)

      !
      ! If magnetic field, rotate TKQ
      !
      if (Bfield%Bstrength(1).gt.TINYB) then

        ! Allocate
        allocate(TB(0:3,-2:2,0:2))
        call Btens(TS,TB,Flgsg,Bfield%Btheta(1),Bfield%Bphi(1))

        ! For each atom, diagonalize Hamiltonian
        do ia=1,nA
          call diagon(Atom(ia),Bfield,Input%zeeman_mode,Flgsg)
        end do

      ! No field
      else

        TB => TS

      end if

      ! Control
      if (laborted) return

      !
      ! Compute hydrogen/electron number density if not in input
      !
      if (Atmo%typo.gt.0.or.Input%keep_atmo) then

        ! Allocate nlte, depar as no-info
        allocate(nlte(Atmo%nele),depar(Atmo%nele))
        nlte = 0
        depar = 0

        ! Call equation of state
        call eqstate(Atmo,Atom,Atomb,nlte,depar)

        ! Deallocate
        deallocate(nlte,depar)

        ! Control
        if (laborted) goto 1000

        ! If the hydrogen densities were computed here, protect
        ! hydrogen from chemical equilibrium
        if (Atmo%typo.gt.0) then

          ! Check active atoms
          do ia=1,nA
            if (Atom(ia)%element.eq.' H') then
              Atom(ia)%mol_protect = Input%protect_H
              exit
            end if
          end do

          ! Check passive atoms
          do ia=1,nAb
            if (Atomb(ia)%element.eq.' H') then
              Atomb(ia)%mol_protect = Input%protect_H
              exit
            end if
          end do

        end if ! Compute eq. state
      end if ! Go trough eq. state at least for writing

      !
      ! Photoionization quantites, thermal part
      !
      do ia=1,nA
        if (Atom(ia)%nphot.lt.1) cycle
        call setphotoTEI(Atom(ia),Frec,Atmo%T,Atmo%ne,MPID)
      end do


      !
      ! Get atomic populations
      !

      ! Active atoms
      do ia=1,nA
        call Initpopu_CLE(Atom(ia),Atmo,bion,ion_column_ind(ia), &
                          ion_value,ion_value_ind(ia),chianti, &
                          ix,.True.)
        call Initcols(Atom(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.True.)
      end do

      ! Passive atoms
      do ia=1,nAb
        call Initpopu_CLE(Atomb(ia),Atmo,bion,ion_column_ind(ia), &
                          ion_value,ion_value_ind(ia),chianti, &
                          ix,.False.)
        call Initcols(Atomb(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.False.)
      end do

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
      call chemeq(Atom,Atomb,dummy,Mol,Atmo)

      ! Control
      if (laborted) goto 1000


      !
      ! Calculate line broadenings
      !

      ! For active atoms
      do ia=1,nA
        call broad(Atom(ia),Atmo,Input%folder,Input%keep_aparam)
      end do

      ! Control
      if (laborted) goto 1000

      ! And for background atoms
      if (Input%nAb.gt.0) then
        do ia=1,Input%nAb
          call broad(Atomb(ia),Atmo,Input%folder,.False.)
        end do
      end if

      ! Control
      if (laborted) goto 1000

      !
      ! Update LOS in Geom
      !
      Geom%L_mu(1) = cos(GeomP%geom(1))
      Geom%L_phi(1) = GeomP%geom(2)

      !
      ! Prepare Geometry for background
      !
      GeomS%nth = 1
      GeomS%nph = 1
      GeomS%nthlos = 0
      GeomS%nphlos = 0

      !
      ! Calculate background continuum quantities
      !
      call background(Atom,Atomb,Mol,Atmo,fudge,kurucz, &
                      Input,Frec%omega,Cont,GeomS, &
                      MPID,Flgsg)

      ! Control
      if (laborted) goto 1000

      !
      ! You can remove pressures now, if not using them later
      !
      if (allocated(Atmo%Pg)) then
        deallocate(Atmo%Pg)
        deallocate(Atmo%Pe)
        deallocate(Atmo%rho)
        deallocate(nlte)
        deallocate(depar)
      end if

      !
      ! Get radiation field
      !
      call getSEEJ(Atom,Atmo,Input%T_rad,Bfield,Flgsg,Frec,spect, &
                   Geom,GeomP,MPID,JKQC,JKQ,Jphot)

      !
      ! Solve SEE
      !

      ! For each atom
      do ia=1,nA

        ! Limiting indexes
        it0 = Atom(ia)%tshift + 1
        it1 = it0 + Atom(ia)%ntran - 1
        ip0 = Atom(ia)%pshift + 1
        ip1 = ip0 + Atom(ia)%nphot - 1

        ! Set magnetic data
        larmor = B2L*Bfield%Bstrength(1)

        ! Solve the SEE
        call SEE(Atom(ia),JKQ(:,:,it0:it1),JKQ(:,:,it0:it1), &
                 Jphot(ip0:ip1,:),larmor,Flgsg,1,-1)
      end do

      !
      ! RT coefficients
      !
      call RTcoeff_CLE(Frec,Atom,Atmo,MPID,Flgsg,Geom,GeomP,Bfield, &
                       TS,TB,if0,if1,JKQ,JKQC(:,:,if0:if1),spect, &
                       Cont%c(:,:,1,1),data1(:,:,0:4))

      !
      ! Boundary radiation
      !

      ! If bottom point
      if (isbottom) then

        ! If in the disk
        if (indisk) then

          ! If input spectra
          if (Input%lspect_input.and.spect%valid) then

            ! Get background Stokes from input spectra
            call get_bottom_spect(Atmo,Frec%omega_ou(if0:if1), &
                                  Atmo%vx(1),Atmo%vy(1),Atmo%vz(1), &
                                  GeomP,spect,data1(:,:,5))

          ! No input spectra, then Planckian
          else

            ! Initialize Stokes to Disk radiation
            call bottom(Frec%omega_ou,Input%T_rad, &
                        Atmo%vx(1),Atmo%vy(1), &
                        Atmo%vz(1),cos(GeomP%geom(1)), &
                        1d0,1d0,MPID,data1(:,:,5))

          end if ! If input spectra
        end if ! In disk
      end if ! If bottom

      ! Free TB
1000  if (Bfield%Bstrength(1).gt.TINYB) deallocate(TB)
      nullify(TB)

      ! Free memory allocated for this node
      call free_cle_node(Atom,Atomb,Mol,Atmo,Bfield,Geom)

      return

      end subroutine getRTCLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module getrtcle_mod
