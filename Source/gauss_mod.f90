      !> Angular quadrature and geometrical tensors
      module gauss_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     18/04/2017
!  Last version:
!     09/02/2026 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/02/2026:    V4.0.3 - The magnetic field is checked for lack
!                             of axial symmetry only if multipoles
!                             other than K=0 are allowed (TdPA)
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
!  gauss
!    Generate the necessary directional quadratures, organize the
!  lines of sight, and initialize the geometrical tensors in the
!  vertical reference frame is necessary
!
!  setTS
!    Calculate the geometrical irreducible spherical tensors in
!  the vertical reference frame
!
!  setTB
!    Calculate the geometrical irreducible spherical tensors in the
!  magnetic field reference frame from the ones in the vertical
!  reference frame. It also checks the axiallity of the RTE
!
!  setTKLOS
!    Calculate the geometrical irreducible spherical tensors in the
!  vertical and magnetic field reference frames for a given line of
!  sight
!
!  check_axial
!    Check that the axiality, if existent, is consistent with the
!  model atmosphere
!
!  geom_index
!    Index the quadrature directions or the LOS directions to in a
!  contintiguous array
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use btens_mod
      use gaussaux_mod
      use stens_mod
      use parameters_mod , only : PI , TINYB

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate the necessary directional quadratures, organize the
      !! lines of sight, and initialize the geometrical tensors in the
      !! vertical reference frame is necessary\n
      !!     Input(Input_class): Structure with configuration data\n
      !!  GeomI(Geometry_class): Structure with geometric data for the
      !!                         intensity problem\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!          mode(integer): Type of run informing about the type
      !!                         of quadrature that is necessary:\n
      !!                             1: RT\n
      !!                             2: CLE\n
      !!            lp(logical): If solving the polarized problem\n
      !!            le(logical): If generating emergent profiles\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                         J-symbols\n
      subroutine gauss(Input,GeomI,Geom,mode,lp,le,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(inout):: Geom,GeomI
      logical, intent(in):: lp,le
      integer, intent(in):: mode

      ! Local

      integer:: ii,jj,kk

      complex(kind=8):: TS(0:3,-2:2,0:2)


      ! Nullify pointers
      nullify(GeomI%TS,GeomI%TSL,GeomI%TB,GeomI%TBL)
      nullify(Geom%TS,Geom%TSL,Geom%TB,Geom%TBL)

      !
      ! CLE
      !
      if (mode.eq.2) then

        !
        ! We only define the gaussian quadrature
        !

        ! Translate into Geom indexes
        Geom%nTh = Input%nTh
        Geom%nPh = Input%nPh
        GeomI%nTh = Geom%nTh
        GeomI%nPh = Geom%nPh

        ! Allocate Needed quantities
        allocate(Geom%V_gauss(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_gauss)
        allocate(Geom%W_gauss(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_gauss)
        allocate(Geom%V_mu(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_mu)
        allocate(Geom%V_mu_disk(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_mu_disk)
        allocate(Geom%V_theta(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_theta)
        allocate(Geom%W_mu(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_mu)

        ! Nodes and weights for simple gaussian quadrature
        call gaussaux(-1d0,1d0,Geom%V_gauss,Geom%W_gauss,Geom%nTh)

        ! Include the 1/2 factor (i.e., integrate to 1)
        Geom%W_gauss = Geom%W_gauss/sum(Geom%W_gauss)

        ! Dummy LOS angles
        allocate(Geom%L_theta(1))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%L_theta)
        allocate(Geom%L_mu(1))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%L_mu)
        allocate(Geom%L_phi(1))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%L_phi)
        Geom%L_mu(1) = 0d0
        Geom%L_theta(1) = 0d0
        Geom%L_phi(1) = 0d0

      !
      ! Non-CLE modes
      !
      else if (mode.eq.1) then

        ! Get from Input
        GeomI%nTh = Input%nThI
        GeomI%nPh = Input%nPhI
        Geom%nTh = Input%nTh
        Geom%nPh = Input%nPh
        GeomI%nThLOS = Input%nThLOS
        GeomI%nPhLOS = Input%nPhLOS
        Geom%nThLOS = Input%nThLOS
        Geom%nPhLOS = Input%nPhLOS

        !
        !
        ! Polar quadrature
        ! Define a global quadrature from a quadrature in each
        ! hemisphere
        !
        !

        !
        ! Intensity
        !

        ! The actual number of nodes is twice the input
        GeomI%nTh = GeomI%nTh*2
        ! Vector with the cos of the nodes
        allocate(GeomI%V_mu(GeomI%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_mu)
        ! Vector with the weights of the integral
        allocate(GeomI%W_mu(GeomI%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%W_mu)
        ! Vector with the angles of the nodes
        allocate(GeomI%V_theta(GeomI%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_theta)

        ! Get gaussian polar quadrature
        call fullgauss(GeomI%nTh,GeomI%V_mu,GeomI%W_mu)

        ! Store the actual angle in each node
        do ii=1,GeomI%nTh
          GeomI%V_theta(ii) = acos(GeomI%V_mu(ii))
        end do

        !
        ! Polarization
        !

        ! Define a global quadrature from a quadrature in each
        ! hemisphere
        ! The actual number of nodes is twice the input
        Geom%nTh = Geom%nTh*2
        ! Vector with the cos of the nodes
        allocate(Geom%V_mu(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_mu)
        ! Vector with the weights of the integral
        allocate(Geom%W_mu(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_mu)
        ! Vector with the angles of the nodes
        allocate(Geom%V_theta(Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_theta)

        ! Get gaussian polar quadrature
        call fullgauss(Geom%nTh,Geom%V_mu,Geom%W_mu)

        ! Store the actual angle in each node
        do ii=1,Geom%nTh
          Geom%V_theta(ii) = acos(Geom%V_mu(ii))
        end do

      end if ! Running mode

      !
      ! Axial quadrature
      !

      ! Intensity is non-axially symmetric
      if(GeomI%nPh.ge.1)then

        ! Flag
        GeomI%axial = .False.

        ! Nodes
        GeomI%nPh = GeomI%nPh*4
        GeomI%nPh2 = GeomI%nPh

        ! Vector with the azimuthal angle
        allocate(GeomI%V_phi(GeomI%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_phi)
        ! Vector with the cos of the azimuthal angle
        allocate(GeomI%V_mux(GeomI%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_mux)
        ! Sign of the sin of the azimuthal angle
        allocate(GeomI%V_muy(GeomI%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_muy)
        ! Weight of the azimuth integral
        allocate(GeomI%W_mux(GeomI%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%W_mux)
        ! Weight of the azimuth integral in emiss2
        allocate(GeomI%W_mux2(GeomI%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%W_mux2)

        ! Get trapezoidal quadrature for azimuth
        call fullazimuth(GeomI%nPh,GeomI%V_phi,GeomI%V_mux, &
                         GeomI%V_muy,GeomI%W_mux)

        ! The internal and external quadratures are the same
        GeomI%W_mux2 = GeomI%W_mux

      ! Intensity is axially symmetric
      else

        ! Flag
        GeomI%axial = .True.

        ! Nodes
        GeomI%nPh = 1
        ! This quantity is only for PRD AD
        if (AVI) then
          GeomI%nPh2 = 1
        else
          GeomI%nPh2 = 8
        end if

        ! Vector with the azimuthal angle
        allocate(GeomI%V_phi(GeomI%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_phi)
        ! Vector with the cos of the azimuthal angle
        allocate(GeomI%V_mux(GeomI%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_mux)
        ! Sign of the sin of the azimuthal angle
        allocate(GeomI%V_muy(GeomI%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_muy)
        ! Weight of the azimuth integral
        allocate(GeomI%W_mux(GeomI%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%W_mux)
        ! Weight of the azimuth integral in emiss2
        allocate(GeomI%W_mux2(GeomI%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(GeomI%W_mux2)

        ! There is no integral for the formal solution
        GeomI%W_mux = 0d0
        GeomI%W_mux(1) = 1d0

        ! Get trapezoidal quadrature for azimuth
        call fullazimuth(GeomI%nPh2,GeomI%V_phi,GeomI%V_mux, &
                         GeomI%V_muy,GeomI%W_mux2)

      end if ! axial symmetry

      ! Polarization is non-axially symmetric
      if(Geom%nPh.ge.1)then

        ! Flag
        Geom%axial = .False.

        ! Nodes
        Geom%nPh = Geom%nPh*4
        Geom%nPh2 = Geom%nPh

        ! Vector with the azimuthal angle
        allocate(Geom%V_phi(Geom%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_phi)
        ! Vector with the cos of the azimuthal angle
        allocate(Geom%V_mux(Geom%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_mux)
        ! Sign of the sin of the azimuthal angle
        allocate(Geom%V_muy(Geom%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_muy)
        ! Weight of the azimuth integral
        allocate(Geom%W_mux(Geom%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_mux)
        ! Weight of the azimuth integral in emiss2
        allocate(Geom%W_mux2(Geom%nPh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_mux2)

        ! Get trapezoidal quadrature for azimuth
        call fullazimuth(Geom%nPh,Geom%V_phi,Geom%V_mux, &
                         Geom%V_muy,Geom%W_mux)

        ! The internal and external quadratures are the same
        Geom%W_mux2 = Geom%W_mux

      ! Polarization is axially symmetric
      else

        ! Flag
        Geom%axial = .True.

        ! Nodes
        Geom%nPh = 1
        ! This quantity is only for PRD AD
        if (AV) then
          Geom%nPh2 = 1
        else
          Geom%nPh2 = 8
        end if

        ! Vector with the azimuthal angle
        allocate(Geom%V_phi(Geom%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_phi)
        ! Vector with the cos of the azimuthal angle
        allocate(Geom%V_mux(Geom%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_mux)
        ! Sign of the sin of the azimuthal angle
        allocate(Geom%V_muy(Geom%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%V_muy)
        ! Weight of the azimuth integral
        allocate(Geom%W_mux(Geom%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_mux)
        ! Weight of the azimuth integral in emiss2
        allocate(Geom%W_mux2(Geom%nPh2))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%W_mux2)

        ! There is no integral for the formal solution
        Geom%W_mux = 0d0
        Geom%W_mux(1) = 1d0

        ! Get trapezoidal quadrature for azimuth
        call fullazimuth(Geom%nPh2,Geom%V_phi,Geom%V_mux, &
                         Geom%V_muy,Geom%W_mux2)

      end if ! axial symmetry

      !
      ! The CLE case if done here
      !
      if (mode.eq.2) then

        ! Nullify pointers
        nullify(Geom%TS,Geom%TSo,Geom%TSL)
        nullify(Geom%TB,Geom%TBo,Geom%TBL)

        ! Control and leave
        call control
        return

      end if ! CLE case

      ! Store the flag in the common variable
      axiali = GeomI%axial
      axial = Geom%axial


      !
      ! Define and store TKQ tensors
      !

      ! If doing polarization
      if (lp) then

        ! TKQ in the vertical reference frame
        allocate(Geom%TS(0:3,-2:2,0:2,Geom%nPh2*Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%TS)

        ! The gamma angle is taken as 0 (Q>0 radial)
        Geom%gam = 0d0

        ! Initialize directional index
        kk = 0

        ! For each polar direction
        do ii=1,Geom%nTh

          ! For each azimuthal direction
          do jj=1,Geom%nPh2

            ! Advance direction index
            kk = kk + 1

            ! Get geometrical tensor
            call Stens(Geom%V_theta(ii),Geom%V_phi(jj),Geom%gam, &
                       Flgsg,TS)
            Geom%TS(:,:,:,kk) = TS

          end do ! Azimuthal direction
        end do ! Polar direction

        ! Different for out and in
        if (PRD.and..not.AV.and.axial) then

          ! TKQ in the vertical reference frame
          allocate(Geom%TSo(0:3,-2:2,0:2,Geom%nTh))
          MRAMc = MRAMc + 1d-6*sizeof(Geom%TSo)

          ! For each polar direction
          do ii=1,Geom%nTh

            ! Get geometrical tensor
            call Stens(Geom%V_theta(ii),Geom%V_phi(1),Geom%gam, &
                       Flgsg,TS)
            Geom%TSo(:,:,:,ii) = TS

          end do ! Polar direction

        ! Same
        else

          ! Just point
          Geom%TSo => Geom%TS

        end if ! Same or different angular grids
      end if ! Polarization

      !
      ! Transform LOS angles
      !

      ! If calculating emergent profiles
      if (le) then

        ! If doing polarization
        if (lp) then

          ! If there are LOS angles
          if (Geom%nThLOS.gt.0) then

            ! Allocate angles for LOS directions
            allocate(Geom%L_theta(Geom%nThLOS))
            MRAMc = MRAMc + 1d-6*sizeof(Geom%L_theta)
            allocate(Geom%L_mu(Geom%nThLOS))
            MRAMc = MRAMc + 1d-6*sizeof(Geom%L_mu)

          end if ! There are LOS angles

          ! If there are LOS azimuthal angles
          if (Geom%nPhLOS.gt.0) then

            ! Allocate angles for azimuth directions
            allocate(Geom%L_phi(Geom%nPhLOS))
            MRAMc = MRAMc + 1d-6*sizeof(Geom%L_phi)

          end if ! There are LOS angles

          ! For each polar angle LOS
          do ii=1,Geom%nThLOS

            ! Get data from Input
            Geom%L_mu(ii) = Input%L_mu(ii)
            Geom%L_theta(ii) = acos(Geom%L_mu(ii))

          end do ! Polar LOS

          ! For each azimuthal angle LOS
          do jj=1,Geom%nPhLOS

            ! Get data from Input
            Geom%L_phi(jj) =  Input%L_phi(jj)*pi/180d0

          end do ! Azimuthal LOS

        end if ! Doing polarization

        ! No polarization or inversion (always because of
        ! background)

        ! If there are LOS angles
        if (GeomI%nThLOS.gt.0) then

          ! Angles for LOS directions
          allocate(GeomI%L_theta(GeomI%nThLOS))
          MRAMc = MRAMc + 1d-6*sizeof(GeomI%L_theta)
          allocate(GeomI%L_mu(GeomI%nThLOS))
          MRAMc = MRAMc + 1d-6*sizeof(GeomI%L_mu)

        end if ! LOS polar angles

        ! If there are azimuthal angles
        if (GeomI%nPhLOS.gt.0) then

          ! Angles for azimuth directions
          allocate(GeomI%L_phi(GeomI%nPhLOS))
          MRAMc = MRAMc + 1d-6*sizeof(GeomI%L_phi)

        end if ! LOS azimuthal angles

        ! For each polar angle
        do ii=1,GeomI%nThLOS

          ! Get data from input
          GeomI%L_mu(ii) = Input%L_mu(ii)
          GeomI%L_theta(ii) = acos(GeomI%L_mu(ii))

        end do ! Polar LOS

        ! For each azimuthal angle
        do jj=1,GeomI%nPhLOS

          ! Get data from input
          GeomI%L_phi(jj) =  Input%L_phi(jj)*pi/180d0

        end do ! Azimuthal LOS

      end if ! Calculing emergent profiles

      !
      ! Polar quadrature for AA redistribution function
      !

      ! Define a global quadrature from a quadrature in each
      ! hemisphere

      ! The actual number of nodes is twice the input
      GeomI%nThAA = Input%nThAAI*2

      ! Vector with the weights of the integral
      allocate(GeomI%W_muAA(GeomI%nThAA))
      MRAMc = MRAMc + 1d-6*sizeof(GeomI%W_muAA)
      ! Vector with the cos and sin of the nodes
      allocate(GeomI%V_muAA(GeomI%nThAA))
      MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_muAA)
      allocate(GeomI%V_siAA(GeomI%nThAA))
      MRAMc = MRAMc + 1d-6*sizeof(GeomI%V_siAA)

      ! Get gaussian quadrature for the scattering angle in the
      ! redistribution functio
      call fullgauss(GeomI%nThAA,GeomI%V_muAA,GeomI%W_muAA)

      ! For every Gaussian node
      do ii=1,GeomI%nThAA

        ! If zero angle
        if (GeomI%V_muAA(ii).ge.1d0) then

          ! Exact sin
          GeomI%V_siAA(ii) = 0d0

        ! If pi angle
        else if (GeomI%V_muAA(ii).le.-1d0) then

          ! Exact sin
          GeomI%V_siAA(ii) = 0d0

        ! General
        else

          ! Compute sin
          GeomI%V_siAA(ii) = sqrt(1d0 - &
                                  GeomI%V_muAA(ii)*GeomI%V_muAA(ii))

        end if ! Exact 0 or pi angle

      end do ! Gaussian nodes


      ! Define a global quadrature from a quadrature in each
      ! hemisphere

      ! The actual number of nodes is twice the input
      Geom%nThAA = Input%nThAA*2

      ! Vector with the weights of the integral
      allocate(Geom%W_muAA(Geom%nThAA))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%W_muAA)
      ! Vector with the cos and sin of the nodes
      allocate(Geom%V_muAA(Geom%nThAA))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%V_muAA)
      allocate(Geom%V_siAA(Geom%nThAA))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%V_siAA)

      ! Get gaussian quadrature for the scattering angle in the
      ! redistribution functio
      call fullgauss(Geom%nThAA,Geom%V_muAA,Geom%W_muAA)

      ! For every Gaussian node
      do ii=1,Geom%nThAA

        ! If zero angle
        if (Geom%V_muAA(ii).ge.1d0) then

          ! Exact sin
          Geom%V_siAA(ii) = 0d0

        ! If pi angle
        else if (Geom%V_muAA(ii).le.-1d0) then

          ! Exact sin
          Geom%V_siAA(ii) = 0d0

        ! General
        else

          ! Compute sin
          Geom%V_siAA(ii) = sqrt(1d0 - &
                                 Geom%V_muAA(ii)*Geom%V_muAA(ii))

        end if ! Exact 0 or pi angle

      end do ! Gaussian nodes

      !
      ! Sanity check inputs
      !

      ! When we need to keep Stokes
      if (KSTK) then

        ! Intensity is not axial, but polarization is ??
        if (axial.and..not.axiali) then

          ! Error message
          umsg = 'Axial polarization and non-axial intensity '// &
                 'is not allowed and does not make much sense'
          urou = 'gauss'
          call aborted

        end if ! Axial polarization, but not intensity

        ! If both or them are not axial, nodes must coincide
        if (.not.axial.and..not.axiali) then

          ! If number of nodes ir different
          if (GeomI%nPh.ne.Geom%nPh) then

            ! Error message
            umsg = 'You specified two different non-axially '// &
                   'simmetric azimuthal quadratures, they '// &
                   'must coincide if Stokes are to be kept.'
            urou = 'gauss'
            call aborted

          end if ! Different quadratures for non-axial problems
        end if ! Both quadratures are non-axial

        ! Different polar quadratures
        if (GeomI%nTh.ne.Geom%nTh) then

          ! Error message
          umsg = 'Polar quadratures must coincide when '// &
                 'Stokes parameters must be kept.'
          urou = 'gauss'
          call aborted

        end if ! Different polar quadratures
      end if ! Keeping Stokes

      ! Check if everything is fine
      call control
      return

      end subroutine gauss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the geometrical irreducible spherical tensors in
      !! the vertical reference frame\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols
      subroutine setTS(Geom,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(inout):: Geom

      ! Local

      integer:: ii,jj,kk

      complex(kind=8):: TS(0:3,-2:2,0:2)


      ! TKQ in the vertical reference frame
      allocate(Geom%TS(0:3,-2:2,0:2,Geom%nPh2*Geom%nTh))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%TS)

      ! The gamma angle is taken as 0 (Q>0 radial)
      Geom%gam = 0d0

      ! Initialize directional index
      kk = 0

      ! For each polar direction
      do ii=1,Geom%nTh

        ! For each azimuthal direction
        do jj=1,Geom%nPh2

          ! Advance direction index
          kk = kk + 1

          ! Get geometrical tensor
          call Stens(Geom%V_theta(ii),Geom%V_phi(jj),Geom%gam, &
                     Flgsg,TS)
          Geom%TS(:,:,:,kk) = TS

        end do ! Azimuthal direction
      end do ! Polar direction

      end subroutine setTS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the geometrical irreducible spherical tensors in
      !! the magnetic field reference frame from the ones in the
      !! vertical reference frame. It also checks the axiallity of
      !! the RTE\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine setTB(Geom,Flgsg,Bfield)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(Geometry_class), intent(inout):: Geom

      ! Local

      integer:: ii,jj,kk,iz

      complex(kind=8):: TB(0:3,-2:2,0:2)


      ! Allocate TKQ in the magnetic reference frame
      allocate(Geom%TB(0:3,-2:2,0:2,Geom%nPh2*Geom%nTh,Rz0:Rz1))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%TB)

      ! Initialize directional index
      kk = 0

      ! For each polar direction
      do ii=1,Geom%nTh

        ! For each azimuthal direction
        do jj=1,Geom%nPh2

          ! Advance direction index
          kk = kk + 1

          ! For each height node
          do iz=Rz0,Rz1

            ! If there is magnetic field
            if (Bfield%Bstrength(iz).gt.TINYB) then

              ! Rotate to the magnetic field reference frame
              call Btens(Geom%TS(:,:,:,kk),TB,Flgsg, &
                         Bfield%Btheta(iz),Bfield%Bphi(iz))

            ! No magnetic field
            else

              ! Copy vertical reference frame tensor
              TB = Geom%TS(:,:,:,kk)

            end if ! Magnetic field

            ! Save result
            Geom%TB(:,:,:,kk,iz) = TB

          end do ! Heights
        end do ! Azimuths
      end do ! Polar angles

      ! If PRD-AD and axial
      if (PRD.and..not.AV.and.axial) then

        ! Allocate TKQ in the magnetic reference frame
        allocate(Geom%TBo(0:3,-2:2,0:2,Geom%nTh,Rz0:Rz1))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%TBo)

        ! Initialize
        kk = 0

        ! For each polar direction
        do ii=1,Geom%nTh

          ! Advance
          kk = kk + 1

          ! Copy
          Geom%TBo(:,:,:,ii,Rz0:Rz1) = Geom%TB(:,:,:,kk,Rz0:Rz1)

          ! Advance
          kk = kk - 1 + Geom%nPh2

        end do ! Polar direction

      ! Same grid
      else

        ! Just point
        Geom%TBo => Geom%TB

      end if ! Same or different angular grid

      end subroutine setTB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the geometrical irreducible spherical tensors in
      !! the vertical and magnetic field reference frames for a given
      !! line of sight\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine setTKQLOS(Geom,Flgsg,Bfield,ii,jj)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(Geometry_class), intent(inout):: Geom
      integer, intent(in):: ii,jj

      ! Local

      integer:: iz


      ! Clean the TSL variable
      if (associated(Geom%TSL)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TSL)
        deallocate(Geom%TSL)
      end if

      ! Allocate vertical reference frame TKQ
      allocate(Geom%TSL(0:3,-2:2,0:2,1))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%TSL)

      ! Clean the TBL variable
      if (associated(Geom%TBL)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%TBL)
        deallocate(Geom%TBL)
      end if

      ! Allocate magnetic field reference frame TKQ
      allocate(Geom%TBL(0:3,-2:2,0:2,1,nZ))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%TBL)

      ! Calculate TKQ in the vertical reference frame for the LOS
      ! direction
      call Stens(Geom%L_theta(ii),Geom%L_phi(jj), &
                 Geom%gam,Flgsg,Geom%TSL(:,:,:,1))

      ! For every height
      do iz=1,nZ

        ! If there is magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          ! Rotate the vertical reference frame TKQ into the magnetic
          ! field reference frame
          call Btens(Geom%TSL(:,:,:,1),Geom%TBL(:,:,:,1,iz),Flgsg, &
                     Bfield%Btheta(iz),Bfield%Bphi(iz))

        ! No magnetic field
        else

          ! Copy vertical reference frame TKQ
          Geom%TBL(:,:,:,:,iz) = Geom%TSL

        end if ! Magnetic field presence

      end do ! Heights

      end subroutine setTKQLOS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check that the axiality, if existent, is consistent with the
      !!  model atmosphere\n
      !!  vx(double(:)): Velocity field vector along X\n
      !!  vy(double(:)): Velocity field vector along Y\n
      !!   t(double(:)): Magnetic field polar angle
      subroutine check_axial(vx,vy,t)

      ! I/O

      double precision, dimension(:), intent(in):: t,vx,vy

      ! Local

      integer:: iz


      !
      ! If there is horizontal velocity
      if (maxval(vx).gt.0d0.or.maxval(vy).gt.0d0) then

        ! Error message
        umsg = 'You specified axial symmetry in '// &
               'input, but there are horizontal '// &
               'velocities'
        urou = 'hanlert'
        call aborted

      end if ! Horizontal velocity

      ! If multipoles
      if (Kcut.gt.0.or.Krad.gt.0) then

        ! For each height
        do iz=1,nz

          ! If magnetic field is not vertical
          if (t(iz).gt.0d0.and.t(iz).lt.PI) then

            ! Error message
            umsg = 'You specified axial symmetry in '// &
                   'input, but there are non-vertical '// &
                   'magnetic fields'
            urou = 'hanlert'
            call aborted

          end if ! Non-vertical magnetic field

        end do ! Heights

      end if ! Actually including multipoles

      end subroutine check_axial

!#####################################################################
!#####################################################################
!#####################################################################

      !> Index the quadrature directions or the LOS directions to in a
      !! contintiguous array\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!          los(logical): If performing last formal solution
      subroutine geom_index(Geom,los)

      ! I/O

      type(Geometry_class), intent(inout):: Geom
      logical, intent(in):: los

      ! Local

      integer:: jdir,ith,iph,njdir


      ! Clean indexing
      if (allocated(Geom%i_geom)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%i_geom)
        deallocate(Geom%i_geom)
      end if

      ! LOS
      if (los) then

        ! Number of directions
        njdir = Geom%nPhLOS*Geom%nThLOS
        Geom%njdir = njdir

        ! Allocate indexing
        allocate(Geom%i_geom(Geom%nPhLOS,Geom%nThLOS))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%i_geom)

      ! Quadrature
      else

        ! Number of directions
        njdir = Geom%nPh*Geom%nTh
        Geom%njdir = njdir

        ! Allocate indexing
        allocate(Geom%i_geom(Geom%nPh,Geom%nTh))
        MRAMc = MRAMc + 1d-6*sizeof(Geom%i_geom)

      end if


      ! Initialize indexing
      Geom%i_geom = 0

      ! Clean index of polar direction
      if (allocated(Geom%ithv)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%ithv)
        deallocate(Geom%ithv)
      end if

      ! Allocate index of polar direction
      allocate(Geom%ithv(njdir))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%ithv)

      ! Clen index of azimuthal direction
      if (allocated(Geom%iphv)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%iphv)
        deallocate(Geom%iphv)
      end if

      ! Allocate index of azimuthal direction
      allocate(Geom%iphv(njdir))
      MRAMc = MRAMc + 1d-6*sizeof(Geom%iphv)

      ! LOS
      if (los) then

        ! Initialize continguous index
        jdir = 0

        ! For each polar direction
        do ith=1,Geom%nThLOS

          ! For each azimuth
          do iph=1,Geom%nPhLOS

            ! Advance index
            jdir = jdir + 1

            ! Save continuous
            Geom%i_geom(iph,ith) = jdir

            ! Save inverse mappings
            Geom%ithv(jdir) = ith
            Geom%iphv(jdir) = iph

          end do ! Azimuths (LOS)
        end do ! Polar angles (LOS)

      ! Quadrature
      else

        ! Initialize continguous index
        jdir = 0

        ! For each polar direction
        do ith=1,Geom%nTh

          ! For each azimuth
          do iph=1,Geom%nPh

            ! Advance index
            jdir = jdir + 1

            ! Save continuous
            Geom%i_geom(iph,ith) = jdir

            ! Save inverse mappings
            Geom%ithv(jdir) = ith
            Geom%iphv(jdir) = iph

          end do ! Azimuths (quadrature)
        end do ! Polar angles (quadrature)

      end if ! LOS/quadrature

      end subroutine geom_index

!#####################################################################
!#####################################################################
!#####################################################################

      end module gauss_mod
