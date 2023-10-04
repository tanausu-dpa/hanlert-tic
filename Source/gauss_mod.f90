      !> Angular quadrature and geometrical tensors
      module gauss_mod
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
!     04/18/2017
!  Last version:
!     04/12/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     04/12/2023:    V3.0.5 - Bugfix: GeomI is always used in
!                             subroutine background. So it should be
!                             defined in the synthesis of
!                             polarization (HL)
!
!     03/08/2023:    V3.0.4 - The emergence variables in GeomI are
!                             always defined in inversion mode (TdPA)
!
!     02/14/2023:    V3.0.3 - Split the definition of the geometry
!                             between the intensity and the
!                             polarization problems (TdPA)
!                           - Moved into a module the algorithmic
!                             part of quadrature definitions (TdPA)
!
!     11/24/2022:    V3.0.2 - Added a branch for CLE mode where a
!                             normal gaussian quadrature is defined
!                             for 0.5*Int_{-1}^{1} (TdPA)
!                           - Added an early return for the CLE run
!                             mode (TdPA)
!
!     10/25/2022:    V3.0.1 - Changed the range in height loops (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case added the
!                             initialization of the TB geometrical
!                             tensor is in the new setTB routine
!                             instead of in gauss. Moreover,
!                             Atmo%v has changed to Atmo%vx,%vy, and
!                             %vz (TdPA)
!                           - The LOS geometrical tensors are now
!                             initialized just before they are
!                             needed by calling setTKQLOS, instead
!                             of doing it from the beginning in the
!                             gauss routine (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added routine check_axial to sanity
!                             check the axial input and its
!                             consistency (TdPA)
!
!     12/10/2019:    V1.1.2 - Now admits no LOS angles (TdPA)
!
!     11/19/2019:    V1.1.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!                           - Specific TINY variable (TdPA)
!
!     01/30/2019:    V1.0.6 - RTaxial is defined here (TdPA)
!                           - Added dependence on commons, which was
!                             missing (TdPA)
!
!     11/28/2018:    V1.0.5 - Removed B2L from used parameters and
!                             added TINY100 (TdPA)
!
!     11/26/2018:    V1.0.4 - Only rotate TS to TB if there is
!                             magnetic field (TdPA)
!
!     05/16/2018:    V1.0.3 - Bugfix: Azimuthal weights from AD
!                             integral in emiss2ord was summing over
!                             azimuthal nodes in the atmosphere, nPh,
!                             instead of the specific for the
!                             redistribution, nPh2 (TdPA)
!
!     07/20/2017:    V1.0.2 - If axial and AA, nPh2=1 too (TdPA)
!
!     06/08/2017:    V1.0.1 - The number of quadrature directions for
!                             the angle average is now an input (TdPA)
!
!     04/18/2017:    V1.0.0 - First version (TdPA)
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
!  gauss
!    This routines initializes the angular quadratures, angles and
!  geometrical tensors in the vertical reference frame
!
!  setTB
!    Calculate the geometrical tensors in the magnetic field reference
!  frame from the ones in the vertical reference frame. It also checks
!  the axiallity of the RTE
!
!  setTKLOS
!    Calculate the geometrical tensors in the vertical and magnetic
!  field reference frames for a given line of sight
!
!  check_axial
!    If the user has chosen axial symmetry, this routine checks if it
!  is consistent with the velocity and magnetic field vectors
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

      !> Generates directional quadratures, with nodes and weights,
      !! and computes geometrical tensors\n
      !!       Input(Input_class): Structure with settings data\n
      !!    GeomI(Geometry_class): Structure with geometry data for
      !!     Geom(Geometry_class): Structure with geometry data\n
      !!                           the intensity problem\n
      !!            mode(integer): Identify the use of the quadrature
      !!                           to define:\n
      !!                             1: RT
      !!                             2: CLE
      !!              lp(logical): Doing the polarization problem\n
      !!              le(logical): There will be emergence in this
      !!                           geometry\n
      !!       Flgsg(Fctsg_class): Structure with factorials and signs
      subroutine gauss(Input,GeomI,Geom,mode,lp,le,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(inout):: Geom,GeomI
      logical, intent(in):: lp,le
      integer, intent(in):: mode

      ! Local

      integer:: ii, jj

      complex(kind=8):: TS(0:3,-2:2,0:2)

      !
      ! In CLE, we only define the gaussian quadrature
      !
      if (mode.eq.2) then

        ! Translate into Geom indexes
        Geom%nTh = Input%nTh
        Geom%nPh = Input%nPh

        ! Allocate Needed quantities
        allocate(Geom%V_gauss(Geom%nTh))
        allocate(Geom%W_gauss(Geom%nTh))
        allocate(Geom%V_mu(Geom%nTh))
        allocate(Geom%V_mu_disk(Geom%nTh))
        allocate(Geom%V_theta(Geom%nTh))
        allocate(Geom%W_mu(Geom%nTh))

        ! Nodes and weights for simple gaussian quadrature
        call gaussaux(-1d0,1d0,Geom%V_gauss,Geom%W_gauss,Geom%nTh)

        ! Include the 1/2 factor (i.e., integrate to 1)
        Geom%W_gauss = Geom%W_gauss/sum(Geom%W_gauss)

        ! Dummy LOS angles
        allocate(Geom%L_theta(Geom%nTh))
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
        ! Polar quadrature
        !

        ! Define a global quadrature from a quadrature in each
        ! hemisphere
        ! The actual number of nodes is twice the input
        GeomI%nTh = GeomI%nTh*2
        ! Vector with the cos of the nodes
        allocate(GeomI%V_mu(GeomI%nTh))
        ! Vector with the weights of the integral
        allocate(GeomI%W_mu(GeomI%nTh))
        ! Vector with the angles of the nodes
        allocate(GeomI%V_theta(GeomI%nTh))

        ! Get quadrature
        call fullgauss(GeomI%nTh,GeomI%V_mu,GeomI%W_mu)

        ! Store the actual angle
        do ii=1,GeomI%nTh

          GeomI%V_theta(ii) = acos(GeomI%V_mu(ii))

        end do

        ! Define a global quadrature from a quadrature in each
        ! hemisphere
        ! The actual number of nodes is twice the input
        Geom%nTh = Geom%nTh*2
        ! Vector with the cos of the nodes
        allocate(Geom%V_mu(Geom%nTh))
        ! Vector with the weights of the integral
        allocate(Geom%W_mu(Geom%nTh))
        ! Vector with the angles of the nodes
        allocate(Geom%V_theta(Geom%nTh))

        ! Get quadrature
        call fullgauss(Geom%nTh,Geom%V_mu,Geom%W_mu)

        ! Store the actual angle
        do ii=1,Geom%nTh

          Geom%V_theta(ii) = acos(Geom%V_mu(ii))

        end do

      end if ! Running mode

      !
      ! Axial quadrature
      !

      ! If not axial symmetry
      if(GeomI%nPh.ge.1)then

        GeomI%axial = .False.
        GeomI%nPh = GeomI%nPh*4
        GeomI%nPh2 = GeomI%nPh
        ! Vector with the azimuthal angle
        allocate(GeomI%V_phi(GeomI%nPh))
        ! Vector with the cos of the azimuthal angle
        allocate(GeomI%V_mux(GeomI%nPh))
        ! Sign of the sin of the azimuthal angle
        allocate(GeomI%V_muy(GeomI%nPh))
        ! Weight of the azimuth integral
        allocate(GeomI%W_mux(GeomI%nPh))
        ! Weight of the azimuth integral in emiss2
        allocate(GeomI%W_mux2(GeomI%nPh))

        ! Get quadrature
        call fullazimuth(GeomI%nPh,GeomI%V_phi,GeomI%V_mux, &
                         GeomI%V_muy,GeomI%W_mux)

        ! They are the same
        GeomI%W_mux2 = GeomI%W_mux

      ! If axial symmetry
      else

        GeomI%axial = .True.
        GeomI%nPh = 1
        ! This quantity is only for PRD AD
        if (AVI) then
          GeomI%nPh2 = 1
        else
          GeomI%nPh2 = 8
        end if
        ! Vector with the azimuthal angle
        allocate(GeomI%V_phi(GeomI%nPh2))
        ! Vector with the cos of the azimuthal angle
        allocate(GeomI%V_mux(GeomI%nPh2))
        ! Sign of the sin of the azimuthal angle
        allocate(GeomI%V_muy(GeomI%nPh2))
        ! Weight of the azimuth integral
        allocate(GeomI%W_mux(GeomI%nPh2))
        ! Weight of the azimuth integral in emiss2
        allocate(GeomI%W_mux2(GeomI%nPh2))

        ! There is no integral for the formal solution
        GeomI%W_mux = 0d0
        GeomI%W_mux(1) = 1d0

        ! Get quadrature
        call fullazimuth(GeomI%nPh2,GeomI%V_phi,GeomI%V_mux, &
                         GeomI%V_muy,GeomI%W_mux2)

      end if ! axial symmetry

      ! If not axial symmetry
      if(Geom%nPh.ge.1)then

        Geom%axial = .False.
        Geom%nPh = Geom%nPh*4
        Geom%nPh2 = Geom%nPh
        ! Vector with the azimuthal angle
        allocate(Geom%V_phi(Geom%nPh))
        ! Vector with the cos of the azimuthal angle
        allocate(Geom%V_mux(Geom%nPh))
        ! Sign of the sin of the azimuthal angle
        allocate(Geom%V_muy(Geom%nPh))
        ! Weight of the azimuth integral
        allocate(Geom%W_mux(Geom%nPh))
        ! Weight of the azimuth integral in emiss2
        allocate(Geom%W_mux2(Geom%nPh))

        ! Get quadrature
        call fullazimuth(Geom%nPh,Geom%V_phi,Geom%V_mux, &
                         Geom%V_muy,Geom%W_mux)

        ! They are the same
        Geom%W_mux2 = Geom%W_mux

      ! If axial symmetry
      else

        Geom%axial = .True.
        Geom%nPh = 1
        ! This quantity is only for PRD AD
        if (AV) then
          Geom%nPh2 = 1
        else
          Geom%nPh2 = 8
        end if
        ! Vector with the azimuthal angle
        allocate(Geom%V_phi(Geom%nPh2))
        ! Vector with the cos of the azimuthal angle
        allocate(Geom%V_mux(Geom%nPh2))
        ! Sign of the sin of the azimuthal angle
        allocate(Geom%V_muy(Geom%nPh2))
        ! Weight of the azimuth integral
        allocate(Geom%W_mux(Geom%nPh2))
        ! Weight of the azimuth integral in emiss2
        allocate(Geom%W_mux2(Geom%nPh2))

        ! There is no integral for the formal solution
        Geom%W_mux = 0d0
        Geom%W_mux(1) = 1d0

        ! Get quadrature
        call fullazimuth(Geom%nPh2,Geom%V_phi,Geom%V_mux, &
                         Geom%V_muy,Geom%W_mux2)

      end if ! axial symmetry

      !
      ! The CLE case if done here
      !
      if (mode.eq.2) then

        ! Control
        call control

        ! And leave
        return

      end if ! CLE case

      ! Store the flag in the common variable
      axiali = GeomI%axial
      axial = Geom%axial


      !
      ! Define and store TKQ tensors
      !

      !
      ! Polarization
      !

      if (lp) then

        ! TKQ in the vertical reference frame
        allocate(Geom%TS(0:3,-2:2,0:2,Geom%nPh2,Geom%nTh))

        ! The gamma angle is taken as 0
        Geom%gam = 0d0

        ! Calculate TKQ in the vertical reference frame
        do ii=1,Geom%nTh
          do jj=1,Geom%nPh2

            call Stens(Geom%V_theta(ii),Geom%V_phi(jj),Geom%gam, &
                       Flgsg,TS)
            Geom%TS(:,:,:,jj,ii) = TS

          end do
        end do

      end if

      !
      ! Transform LOS angles

      ! If emergence
      if (le) then

        ! If polarization
        if (lp) then

          ! If there are LOS angles
          if (Geom%nThLOS.gt.0) then
            ! Angles for LOS directions
            allocate(Geom%L_theta(Geom%nThLOS))
            allocate(Geom%L_mu(Geom%nThLOS))
          end if

          ! If there are azimuthal angles
          if (Geom%nPhLOS.gt.0) then
            ! Angles for azimuth directions
            allocate(Geom%L_phi(Geom%nPhLOS))
          end if

          ! For each polar angle
          do ii=1,Geom%nThLOS
            Geom%L_mu(ii) = Input%L_mu(ii)
            Geom%L_theta(ii) = acos(Geom%L_mu(ii))
          end do
          ! For each azimuthal angle
          do jj=1,Geom%nPhLOS
            Geom%L_phi(jj) =  Input%L_phi(jj)*pi/180d0
          end do

        end if

        ! No polarization or inversion (always because of
        ! background)

        ! If there are LOS angles
        if (GeomI%nThLOS.gt.0) then
          ! Angles for LOS directions
          allocate(GeomI%L_theta(GeomI%nThLOS))
          allocate(GeomI%L_mu(GeomI%nThLOS))
        end if

        ! If there are azimuthal angles
        if (GeomI%nPhLOS.gt.0) then
          ! Angles for azimuth directions
          allocate(GeomI%L_phi(GeomI%nPhLOS))
        end if

        ! For each polar angle
        do ii=1,GeomI%nThLOS
          GeomI%L_mu(ii) = Input%L_mu(ii)
          GeomI%L_theta(ii) = acos(GeomI%L_mu(ii))
        end do
        ! For each azimuthal angle
        do jj=1,GeomI%nPhLOS
          GeomI%L_phi(jj) =  Input%L_phi(jj)*pi/180d0
        end do

      end if ! Emergence

      !
      ! Polar quadrature for AA redistribution function
      !

      ! Define a global quadrature from a quadrature in each
      ! hemisphere
      ! The actual number of nodes is twice the input
      GeomI%nThAA = Input%nThAAI*2
      ! Vector with the cos of the nodes
      allocate(GeomI%W_muAA(GeomI%nThAA))
      ! Vector with the weights of the integral
      allocate(GeomI%V_thetaAA(GeomI%nThAA))

      ! Get quadrature
      call fullgauss(GeomI%nThAA,GeomI%V_thetaAA,GeomI%W_muAA)
      GeomI%V_thetaAA = acos(GeomI%V_thetaAA)


      ! Define a global quadrature from a quadrature in each
      ! hemisphere
      ! The actual number of nodes is twice the input
      Geom%nThAA = Input%nThAA*2
      ! Vector with the cos of the nodes
      allocate(Geom%W_muAA(Geom%nThAA))
      ! Vector with the weights of the integral
      allocate(Geom%V_thetaAA(Geom%nThAA))


      ! Get quadrature
      call fullgauss(Geom%nThAA,Geom%V_thetaAA,Geom%W_muAA)
      Geom%V_thetaAA = acos(Geom%V_thetaAA)

      !
      ! Sanity check inputs
      !

      ! When we need to keep Stokes
      if (KSTK) then

        ! Intensity is not axial, but polarization is ??
        if (axial.and..not.axiali) then

          umsg = 'Axial polarization and non-axial intensity '// &
                 'is not allowed and does not make much sense'
          urou = 'gauss'
          call aborted

        end if ! Axial polarization, but not intensity

        ! If both or them are not axial, nodes must coincide
        if (.not.axial.and..not.axiali) then

          if (GeomI%nPh.ne.Geom%nPh) then
            umsg = 'You specified two different non-axially '// &
                   'simmetric azimuthal quadratures, they '// &
                   'must coincide if Stokes are to be kept.'
            urou = 'gauss'
            call aborted
          end if ! Different quadratures

        end if ! Both quadratures are non-axial

        ! Different polar quadratures
        if (GeomI%nTh.ne.Geom%nTh) then

          umsg = 'Polar quadratures must coincide when '// &
                 'Stokes parameters must be kept.'
          urou = 'gauss'
          call aborted

        end if ! Different polar quadratures

      end if ! Keeping Stokes

      ! Check if everything is fine
      call control

      return

      end subroutine Gauss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set the magnetic field geometrical tensors\n
      !!  Geom(Geometry_class): Structure with geometry data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine setTB(Geom,Flgsg,Bfield)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(inout):: Bfield
      type(Geometry_class), intent(inout):: Geom

      ! Local
      integer:: ii,jj,iz
      complex(kind=8):: TB(0:3,-2:2,0:2)

      ! TKQ in the magnetic reference frame
      allocate(Geom%TB(0:3,-2:2,0:2,Geom%nPh2,Geom%nTh,Rz0:Rz1))

      ! Calculate TKQ in the magnetic field reference frame
      do ii=1,Geom%nTh
        do jj=1,Geom%nPh2

          ! Rotate them to obtain the TKQ in the magnetic reference
          ! frame
          do iz=Rz0,Rz1

            ! If there is magnetic field
            if (Bfield%Bstrength(iz).gt.TINYB) then
              call Btens(Geom%TS(:,:,:,jj,ii),TB,Flgsg, &
                         Bfield%Btheta(iz),Bfield%Bphi(iz))
            ! No magnetic field
            else
              TB = Geom%TS(:,:,:,jj,ii)
            end if

            Geom%TB(:,:,:,jj,ii,iz) = TB

          end do

        end do
      end do

      !
      ! Check RT axiality
      !
      RTaxial = axial.and. &
                maxval(Bfield%Bstrength).le.TINYB

      end subroutine setTB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set TKQ tensors for an specific LOS\n
      !!  Geom(Geometry_class): Structure with geometry data\n
      !!    Flgsg(Fctsg_class): Structure with factorials and
      !!                        signs\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine setTKQLOS(Geom,Flgsg,Bfield,ii,jj)

      ! I/O
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(inout):: Bfield
      type(Geometry_class), intent(inout):: Geom
      integer, intent(in):: ii,jj

      ! Local
      integer:: iz

      ! TKQ in the vertical reference frame for LOS directions
      if (.not.allocated(Geom%TSL)) &
        allocate(Geom%TSL(0:3,-2:2,0:2))

      if (.not.allocated(Geom%TBL)) &
        allocate(Geom%TBL(0:3,-2:2,0:2,nZ))

      ! Calculate TKQ in the vertical reference frame for the LOS
      ! directions
      call Stens(Geom%L_theta(ii),Geom%L_phi(jj), &
                 Geom%gam,Flgsg,Geom%TSL)

      ! Rotate them to obtain the TKQ in the magnetic reference
      ! frame for the LOS directions
      do iz=1,nZ

        ! If there is magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then
          call Btens(Geom%TSL,Geom%TBL(:,:,:,iz),Flgsg, &
                     Bfield%Btheta(iz),Bfield%Bphi(iz))
        ! No magnetic field
        else
          Geom%TBL(:,:,:,iz) = Geom%TSL
        end if

      end do

      end subroutine setTKQLOS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks that axial conditions are respected.\n
      !!  vx(dfloat(:)): Velocity field vector along X\n
      !!  vy(dfloat(:)): Velocity field vector along Y\n
      !!   t(dfloat(:)): Magnetic field polar angle
      subroutine check_axial(vx,vy,t)

      ! I/O
      double precision, dimension(:), intent(in):: t,vx,vy

      ! Local
      integer:: iz

      ! Check velocity vector
      if (maxval(vx).gt.0d0.or.maxval(vy).gt.0d0) then
        umsg = 'You specified axial symmetry in '// &
               'input, but there are horizontal '// &
               'velocities'
        urou = 'hanlert'
        call aborted
      end if

      ! For each height
      do iz=1,nz

        if (t(iz).gt.0d0.and.t(iz).lt.PI) then
          umsg = 'You specified axial symmetry in '// &
                 'input, but there are non-vertical '// &
                 'magnetic fields'
          urou = 'hanlert'
          call aborted
        end if

      end do ! Heights

      end subroutine check_axial

!#####################################################################
!#####################################################################
!#####################################################################

      end module gauss_mod
