      !> Get CLE RT coefficients
      module getrtcle_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     01/10/2022
!  Last version:
!     06/06/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/06/2025:    V4.0.3 - Added argument to chemeq to ensure
!                             compatibility with new input (TdPA)
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
!  getRTCLE
!    Calculate the RT coefficients for a given node in a CLE model. If
!  the node is in a boundary, compute the boundary condition as well
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
      use commons_mod
      use diagon_mod
      use free_mod
      use geometry_mod
      use initcols_mod
      use initphotoion_mod
      use initpopu_mod
      use jcalccle_mod
      use normalizer_mod
      use parameters_mod , only: B2L , TINYB
      use ratmo_mod
      use ratom_mod
      use rmol_mod
      use rtcoeff_mod
      use see_mod
      use stens_mod
      use strength_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return the RT coefficients for a given point along the LOS
      !! for a CLE calculation
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atomb(Atom_class(:)): Structures with atomic data for
      !!                             background atoms\n
      !!          Mol(Mol_class(:)): Structures with molecular data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!             Red(Red_class): Structure with redistribution
      !!                             input frequency data,
      !!                             redistribution function data, and
      !!                             profile or normalization data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!                ix(integer): Index of the location along the
      !!                             LOS\n
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
      !!         data1(double(:,:)): Radiation transfer coefficients\n
      !!            indisk(logical): If the ray crosses the disk\n
      !!          isbottom(logical): If this is the bottom boundary
      subroutine getRTCLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red, &
                          Geom,Flgsg,fudge,kurucz,ix,if0,if1,batmo, &
                          bion,ion_column_ind,ion_value_ind, &
                          ion_value,spect,chianti,data1,indisk, &
                          isbottom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(inout):: Frec
      type(Red_class), intent(in):: Red
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(spect_class), intent(inout):: spect
      type(Geometry_class), intent(inout):: Geom
      type(chianti_class), intent(in):: chianti
      logical, intent(in):: indisk,isbottom
      integer, intent(in):: ix,if0,if1
      integer, dimension(:), intent(in):: ion_value_ind,ion_column_ind
      double precision, dimension(0:3,if0:if1,0:5), &
                        intent(out):: data1
      double precision, dimension(:), intent(in):: ion_value,bion
      double precision, dimension(:), intent(inout):: batmo

      ! Local

      type(Continuum_class):: Cont
      type(Coronapoint_class):: GeomP
      type(Bfield_class):: Bfield
      type(Geometry_class):: GeomS
      type(LTEline_class), dimension(:), allocatable:: dummy

      integer:: ia,it0,it1,ip0,ip1
      integer, dimension(:), allocatable:: nlte,depar

      double precision:: larmor,vfac,ct,st,cc,sc
      double precision, dimension(nxphot,2):: Jphot

      complex(kind=8), dimension(-2:2,0:2,nxtran):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2), target:: TS
      complex(kind=8), dimension(:,:,:), pointer:: TB


      !
      ! Initialize
      data1 = 0d0


      !
      ! Get atmosphere
      call rAtmo_cle(batmo,Input,Atmo,ix)


      !
      ! Prepare atomic and molecular quantities depending on
      ! the atmospheric model
      !

      ! Prepare active atoms
      call prepareatom(Atom,nA)

      ! Prepare background atoms
      if (nAb.gt.0) call prepareatom(Atomb,nAb)

      ! Prepare molecules
      if (Input%nM.gt.0) call preparemol(Mol,Input%nM)

      ! Control
      if (laborted) goto 1001

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
        ! Not counting RAM because this is anecdotal
        allocate(TB(0:3,-2:2,0:2))
        call Btens(TS,TB,Flgsg,Bfield%Btheta(1),Bfield%Bphi(1))

        ! For each atom, diagonalize Hamiltonian and get transition
        ! strengths
        do ia=1,nA
          call diagon(Atom(ia),Bfield,Input%zeeman_mode,Flgsg)
          call strength_ev(Atom(ia),Bfield)
        end do

      ! No field
      else

        ! Just point to vertical reference frame tensors
        TB => TS

      end if

      ! Control
      if (laborted) goto 1000

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

      ! For each atom
      do ia=1,nA

        ! Skip if no photoionizations
        if (Atom(ia)%nphot.lt.1) cycle

        ! Set up thermal rate in SEE
        call setphotoTEI(Atom(ia),Frec,Atmo%T,Atmo%ne,.False.)

      end do ! Atoms


      !
      ! Get atomic populations
      !

      ! Active atoms
      do ia=1,nA

        ! Initialize populations
        call Initpopu_CLE(Atom(ia),Atmo,bion,ion_column_ind(ia), &
                          ion_value,ion_value_ind(ia),chianti, &
                          ix,.True.)

        ! Initialize collisions
        call Initcols(Atom(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.True.)

      end do ! Active atoms

      ! Passive atoms
      do ia=1,nAb

        ! Initialize populations
        call Initpopu_CLE(Atomb(ia),Atmo,bion,ion_column_ind(ia), &
                          ion_value,ion_value_ind(ia),chianti, &
                          ix,.False.)

        ! Initialize collisions
        call Initcols(Atomb(ia),Atmo,Input%folder,Flgsg, &
                      Input%keep_coll,.False.)

      end do ! Passive atoms

      !
      ! Protect atoms from chemical equilibrium
      !

      ! If protecting all atoms
      if (Input%chem_protect_all) then

        ! Protect active atoms
        do ia=1,nA
          Atom(ia)%mol_protect = .True.
        end do

        ! Protect background atoms
        do ia=1,nAb
          Atomb(ia)%mol_protect = .True.
        end do

      end if ! Protecting all atoms

      !
      ! Calculate chemical equilibrium
      !
      call chemeq(Atom,Atomb,dummy,Mol,Atmo,.False.)

      ! Control
      if (laborted) goto 1000


      !
      ! Calculate line broadenings
      !

      ! For active atoms
      do ia=1,nA

        ! Get line damping parameter
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
      ! Update LOS in Geom structure
      !
      Geom%L_mu(1) = cos(GeomP%geom(1))
      Geom%L_phi(1) = GeomP%geom(2)

      !
      ! Prepare Geometry for background quantities calculation
      !
      GeomS%nth = 1
      GeomS%nph = 1
      GeomS%nthlos = 0
      GeomS%nphlos = 0

      !
      ! Calculate background continuum quantities
      !

      ! If considering continuum sources
      if (Input%add_cont_cle) then

        ! Get background continuum
        call background(Atom,Atomb,Mol,Atmo,fudge,kurucz, &
                        Input,Frec%omega,Cont,GeomS, &
                        MPID,Flgsg)

      ! Neglecting continuum
      else

        ! Allocate a dummy array
        allocate(Cont%c(1,1,1,1))
        BRAMc = BRAMc + 1d-6*sizeof(Cont%c)
        Cont%c = 0d0

      end if ! Consider continuum sources

      ! Control
      if (laborted) goto 1000

      !
      ! You can remove pressures now, if not using them later
      !

      ! If gas pressure is allocated
      if (allocated(Atmo%Pg)) then

        ! Clean unused atmospheric variables
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pg)
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pe)
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%rho)
        deallocate(Atmo%Pg)
        deallocate(Atmo%Pe)
        deallocate(Atmo%rho)

        ! Clean auxiliar nlte and depar variables
        deallocate(nlte)
        deallocate(depar)

      end if

      !
      ! Get radiation field
      !
      call getSEEJ(Atom,Atmo,Input%T_rad,Input%use_allen, &
                   Input%flat_cle_in,Bfield,Flgsg,Frec,spect, &
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

        ! Set magnetic field constant
        larmor = B2L*Bfield%Bstrength(1)

        ! Solve the SEE
        call SEE(Atom(ia),JKQ(:,:,it0:it1),JKQ(:,:,it0:it1), &
                 Jphot(ip0:ip1,:),larmor,Flgsg,1)

      end do ! Active atoms

      ! Compute Doppler shift
      vfac = 1d0

      ! If dynamic
      if (dyn) then

        ! Trigonometry
        ct = cos(GeomP%geom(1))
        st = sqrt(1d0 - ct*ct)
        cc = cos(GeomP%geom(2))
        sc = sin(GeomP%geom(2))

        ! Doppler shift factor
        vfac = 1d0 - Atmo%vx(1)*st*cc - atmo%vy(1)*st*sc - &
                     Atmo%vz(1)*ct

      end if

      !
      ! Calculate RT coefficients
      !
      call RTcoeff_CLE(Frec,Red,Atom,Atmo,Flgsg,Geom,GeomP, &
                       vfac,Bfield,TS,TB,if0,if1, &
                       JKQ,JKQC(:,:,if0:if1), &
                       spect,Cont%c(:,:,1,1),Input%add_cont_cle, &
                       data1(:,:,0:4))

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

          ! If using Allen quantities
          else if (Input%use_allen) then

            ! Get background from Allen tabulation
            call get_bottom_allen(Atmo,Frec%omega_ou(if0:if1), &
                                  Atmo%vx(1),Atmo%vy(1),Atmo%vz(1), &
                                  GeomP,data1(:,:,5))

          ! No input spectra, then Planckian
          else

            ! Initialize Stokes to Disk radiation
            ! Static disk for the observer
            call bottom(Frec%omega_ou,Input%T_rad,1d0, &
                        if0,if1,data1(:,:,5))

          end if ! If input spectra
        end if ! In disk
      end if ! If bottom

      ! Free TB
1000  if (Bfield%Bstrength(1).gt.TINYB) deallocate(TB)
1001  nullify(TB)

      ! Free memory allocated for this node
      call free_cle_node(Atom,Atomb,Mol,Atmo,Bfield,Geom)

      return

      end subroutine getRTCLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module getrtcle_mod
