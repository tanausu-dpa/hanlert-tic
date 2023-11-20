      !> Manage the model atmosphere in the inversion
      module model_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/17/2023
!  Last version:
!     14/11/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     14/11/2023:    V3.0.5 - Redid the reference frame rotations
!                             because HL found an issue (TdPA)
!
!     10/04/2023:    V3.0.4 - The azimuth in the LOS is set to zero
!                             if there is no transversal component
!                             when transforming (TdPA)
!
!     10/03/2023:    V3.0.3 - Bugfix: wrong index for the ad-hoc
!                             asymmetry (HL)
!
!     09/28/2023:    V3.0.2 - Just formatting, not advancing version
!                             for this (TdPA)
!
!     07/03/2023:    V3.0.2 - Added Intpol_Var to interpolate the
!                             nodes into the model atmosphere for
!                             a given variable (TdPA)
!                           - Added Intpol_glob to manage the
!                             interpolation of global variables into
!                             the model atmosphere (TdPA)
!                           - Added code for no hydrostatic
!                             equilibrium (TdPA)
!                           - Generalized the transformation of the
!                             magnetic field between vertical and LOS
!                             frames to any LOS azimuth (TdPA)
!                           - Added vconversion and v2vlos to
!                             transform the velocity vector from
!                             cartesian to LOS and viceversa (TdPA)
!
!     03/15/2023:    V3.0.1 - Intpol_Atmo_* do not need the Flgsg
!                             argument (TdPA)
!                           - The Blos variables are in the same
!                             structure than the polar ones (TdPA)
!                           - Bugfix: The nodes by correction were
!                             not reseted for the LOS/POS magnetic
!                             field nodes (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/17/2023:    V0.0.0 - Started from 05/12/2020
!                             TIC@model_mod.f90 revision from
!                             Hao (TdPA)
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
!    Intpol_Atmo_all:
!      Manages calls to interpolate from nodes into atmosphere
!
!    Intpol_Var:
!      Interpolate a variable from the nodes into the model
!    atmosphere
!
!    Intpol_glob:
!      Manages the interpolation of global quantities from nodes
!    into atmosphere
!
!    Intpol_Atmo:
!      Manages the interpolation of thermal quantities from nodes
!    into atmosphere
!
!    Intpol_Bfield:
!      Manages the interpolation of magnetic quantities from nodes
!    into atmosphere
!
!    Bconversion:
!      Transform magnetic field from LOS to vertical
!
!    B2Blos:
!      Transform magnetic field from vertical to LOS
!
!    vconversion:
!      Transform velocity from LOS to cartesian
!
!    v2vlos:
!      Transform velocity from vertical to LOS
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use commons_mod
      use hydrostatic_mod
      use inter_mod
      use parameters_mod, only: PI , TINYA, TINYDP
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate all quantities from nodes into the atmosphere
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!           Atmo(Atmo_class): Structure with the model\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!         Input(Input_class): Structure with settings data\n
      !!         fudge(fudge_class): Structure with fudge data
      subroutine Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield,Atom,Atomb, &
                                 Mol,Input,fudge)

      !IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(fudge_class):: fudge

      ! If thermodynamic inversion
      if (Inf_Nodes%Nodes_type.eq.0) then

        ! Call interpolation of thermal quantities
        call Intpol_Atmo(Inf_Nodes,Atmo,Atom,Atomb,Mol,Input,fudge)

      ! If magnetic inversion
      else if (Inf_Nodes%Nodes_type.eq.1) then

        ! Call interpolation of magnetic quantities
        call Intpol_Bfield(Inf_Nodes,Atmo,Bfield)

      ! Invert all
      else

        ! Call interpolation of thermal quantities
        call Intpol_Atmo(Inf_Nodes,Atmo,Atom,Atomb,Mol,Input,fudge)

        ! Call interpolation of magnetic quantities
        call Intpol_Bfield(Inf_Nodes,Atmo,Bfield)

      end if ! Type of interpolation

      return

      end subroutine Intpol_Atmo_all

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate a variable from nodes into the model atmosphere\n
      !!  Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!        o_var(double(:)): Output stratification\n
      !!            z(double(:)): Optical depth axis\n
      !!             nn(integer): Dimension of optical depth axis\n
      !!           indx(integer): Parameter index
      subroutine Intpol_Var(Inf_Nodes,o_var,z,nn,indx)

      ! I/O
      type(Nodes_class), intent(inout):: Inf_Nodes
      integer, intent(in):: nn, indx
      double precision, dimension(:), intent(in):: z
      double precision, dimension(:), intent(inout):: o_var

      ! Local
      double precision, dimension(nn):: i_var


      ! If inverting
      if (Inf_Nodes%Nodes_Flags(indx)) then

        ! If contain values
        if (Inf_Nodes%Node_Type(indx).le.3) then

          ! Interpolate
          call Intpol(Inf_Nodes%Node(indx)%H, &
                      Inf_Nodes%Node(indx)%Var, &
                      Inf_Nodes%Num_Nodes(indx), z, o_var, nn, &
                      Inf_Nodes%Interpolation, 3)

        ! If contain corrections
        else

          ! Inverpolate
          call Intpol(Inf_Nodes%Node(indx)%H, &
                      Inf_Nodes%Node(indx)%Var, &
                      Inf_Nodes%Num_Nodes(indx), z, i_var, nn, &
                      Inf_Nodes%Interpolation, 3)

          ! Add correction
          o_var = o_var + i_var

          ! Reset correction
          Inf_Nodes%Node(indx)%Var = 0d0

        end if ! Type of node content
      end if ! Inverting

      end subroutine Intpol_Var

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate global quantities from nodes into the model
      !! atmosphere\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!           Atmo(Atmo_class): Structure with the model
      subroutine Intpol_glob(Inf_Nodes,Atmo)

      ! IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Atmo_class), intent(inout):: Atmo

      ! If inverting diffuse light, update
      if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_f)) &
        Atmo%f_diff = Inf_Nodes%Node(Inf_Nodes%index_f)%Var(1)

      end subroutine Intpol_glob

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate thermal quantities from nodes into the model
      !! atmosphere\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!           Atmo(Atmo_class): Structure with the model\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!         Input(Input_class): Structure with settings data\n
      !!         fudge(fudge_class): Structure with fudge data
      subroutine Intpol_Atmo(Inf_Nodes,Atmo,Atom,Atomb, &
                             Mol,Input,fudge)

      ! IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(fudge_class):: fudge

      ! Local
      double precision:: Pg
      double precision, dimension(:), allocatable:: Z


      ! Allocate auxiliar quantities
      allocate(Z(Atmo%nZ))

      ! Get log(tau)
      Z = log10(Atmo%z)

      ! Temperature
      call Intpol_Var(Inf_Nodes,Atmo%T,Z,Atmo%nz,Inf_Nodes%index_T)

      ! Type of velocity
      select case (Inf_Nodes%vtype)

        ! Vertical
        case(0)

          ! X velocity
          call Intpol_Var(Inf_Nodes,Atmo%vx,Z,Atmo%nz, &
                          Inf_Nodes%index_vx)

          ! Y velocity
          call Intpol_Var(Inf_Nodes,Atmo%vy,Z,Atmo%nz, &
                          Inf_Nodes%index_vy)

          ! Vertical velocity
          call Intpol_Var(Inf_Nodes,Atmo%vz,Z,Atmo%nz, &
                          Inf_Nodes%index_vz)

        ! LOS
        case(1)

          ! POS velocity
          call Intpol_Var(Inf_Nodes,Atmo%vpos,Z,Atmo%nz, &
                          Inf_Nodes%index_vx)

          ! Azimuth velocity
          call Intpol_Var(Inf_Nodes,Atmo%vphi,Z,Atmo%nz, &
                          Inf_Nodes%index_vy)

          ! LOS velocity
          call Intpol_Var(Inf_Nodes,Atmo%vlos,Z,Atmo%nz, &
                          Inf_Nodes%index_vz)


          ! Convert to vertical
          call vconversion(Atmo%nZ, Inf_Nodes%mu, Inf_Nodes%azimuth, &
                           Atmo%vlos, Atmo%vpos, Atmo%vphi, &
                           Atmo%vx, Atmo%vy, Atmo%vz)

        ! Error
        case default

          ! Abort
          umsg = 'The index of vtype is not correct.'// &
                 'Only 0 and 1 are accepted.'
          urou = 'Intpol_Atom'
          call aborted
          return

      end select

      ! Microturbulent velocity
      call Intpol_Var(Inf_Nodes,Atmo%vmi,Z, &
                      Atmo%nz,Inf_Nodes%index_vm)

      ! J21R
      call Intpol_Var(Inf_Nodes,Atmo%JKQin(4*nz+1:5*nz),Z, &
                      Atmo%nz,Inf_Nodes%index_J21R)
      ! J21I
      call Intpol_Var(Inf_Nodes,Atmo%JKQin(5*nz+1:6*nz),Z, &
                      Atmo%nz,Inf_Nodes%index_J21I)
      ! J22R
      call Intpol_Var(Inf_Nodes,Atmo%JKQin(6*nz+1:7*nz),Z, &
                      Atmo%nz,Inf_Nodes%index_J22R)
      ! J22I
      call Intpol_Var(Inf_Nodes,Atmo%JKQin(7*nz+1:8*nz),Z, &
                      Atmo%nz,Inf_Nodes%index_J22I)


      ! If hydrostatic
      if (Inf_Nodes%hydroeq) then

        ! If inverting the gas pressure
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg)) then

          ! Take from node
          Pg = Inf_Nodes%Node(Inf_Nodes%index_Pg)%Var(1)

        ! Not inverting gas pressure
        else

          ! Take from input
          Pg = Inf_Nodes%Pg_Bound

        end if

        ! If hydrostatic equilibrium
        if (Inf_Nodes%hydros) then

          ! Compute pressure stratification
          call Compute_Pressure_all(Atmo,Atom,Atomb,Mol,Input, &
                                    fudge,Pg)

        end if ! Hydrostatic equilibrium required due to changes

      ! No hydrostatic equilibrium
      else

        ! Pg
        call Intpol_Var(Inf_Nodes,Atmo%Pg,Z, &
                        Atmo%nz,Inf_Nodes%index_Pg)

      end if ! Doing hydrostatic equilibrium

      ! Deallocate auxiliar variables
      deallocate(Z)

      ! Diffuse light
      call Intpol_glob(Inf_Nodes,Atmo)

      return

      end subroutine Intpol_Atmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate thermal quantities from nodes into the model
      !! atmosphere\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!           Atmo(Atmo_class): Structure with the model\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data
      subroutine Intpol_Bfield(Inf_Nodes, Atmo, Bfield)

      ! IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield

      ! Local
      double precision, dimension(:), allocatable:: Z


      ! Allocate auxiliar variables
      allocate(Z(Atmo%nZ))

      ! Get logtau axis
      Z = log10(Atmo%z)

      ! Select type of magnetic field
      select case(Inf_Nodes%Btype)

        ! Vertical
        case(0)

          ! Bstrength
          call Intpol_Var(Inf_Nodes,Bfield%Bstrength,Z, &
                          Atmo%nz,Inf_Nodes%Index_B)

          ! Btheta
          call Intpol_Var(Inf_Nodes,Bfield%Btheta,Z, &
                          Atmo%nz,Inf_Nodes%Index_Bt)

          ! Bphi
          call Intpol_Var(Inf_Nodes,Bfield%Bphi,Z, &
                          Atmo%nz,Inf_Nodes%Index_Bp)

        ! Magnetic field in LOS
        case(1)

          ! Blos
          call Intpol_Var(Inf_Nodes,Bfield%Blos,Z, &
                          Atmo%nz,Inf_Nodes%Index_B)

          ! Bpos
          call Intpol_Var(Inf_Nodes,Bfield%Bpos,Z, &
                          Atmo%nz,Inf_Nodes%Index_Bt)

          ! Azimuth
          call Intpol_Var(Inf_Nodes,Bfield%Azimuth,Z, &
                          Atmo%nz,Inf_Nodes%Index_Bp)

          ! Convert to vertical
          call Bconversion(Atmo%nZ, Inf_Nodes%mu, Inf_Nodes%azimuth, &
                           Bfield%Blos, Bfield%Bpos, Bfield%Azimuth, &
                           Bfield%Bstrength, Bfield%Btheta, &
                           Bfield%Bphi)

        ! Error
        case default

          ! Abort
          umsg = 'The index of btype is not correct.'// &
                 'Only 0 and 1 are accepted.'
          urou = 'Intpol_Bfield'
          call aborted
          return

      end select ! Type of magnetic field representation

      ! Deallocate auxiliars
      deallocate(Z)

      ! Diffuse light
      call Intpol_glob(Inf_Nodes,Atmo)

      return

      end subroutine Intpol_Bfield

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform from LOS to vertical\n
      !!         Num(integer): Number of nodes\n
      !!           Mu(double): Cosine of the heliocentric angle\n
      !!          phi(double): Azimuth\n
      !!      blos(double(:)): Magnetic field along the LOS\n
      !!      bpos(double(:)): Magnetic field in the POS\n
      !!   azimuth(double(:)): Magnetic azimuth in the POS\n
      !! bstrength(double(:)): Magnetic field strength\n
      !!    btheta(double(:)): Magnetic field inclination\n
      !!      bphi(double(:)): Magnetic field azimuth
      subroutine Bconversion(Num,Mu,phi,blos,bpos,azimuth, &
                             bstrength,btheta,bphi)

      ! IO
      integer, intent(in):: Num
      double precision, intent(in):: Mu,phi
      double precision, dimension(:),intent(in):: blos
      double precision, dimension(:),intent(in):: bpos
      double precision, dimension(:),intent(in):: azimuth
      double precision, dimension(:),intent(inout):: bstrength
      double precision, dimension(:),intent(inout):: btheta
      double precision, dimension(:),intent(inout):: bphi

      ! Local
      integer:: i
      double precision:: Mup, cp, sp
      double precision, dimension(Num):: bx, by, bz, bxp, byp, bzp


      ! Get sin
      if (abs(Mu-1d0).lt.TINYA) then
        Mup = 0d0
      else
        Mup = sqrt(1d0 - Mu*Mu)
      end if

      ! Get cos-sin phi
      if (abs(phi).lt.TINYA) then
        cp =  1d0
        sp =  0d0
      else if (abs(phi-0.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp =  1d0
      else if (abs(phi-PI).lt.TINYA) then
        cp = -1d0
        sp =  0d0
      else if (abs(phi-1.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp = -1d0
      else
        cp = cos(phi)
        sp = sin(phi)
      end if

      ! Get module
      bstrength = sqrt(blos*blos + bpos*bpos)

      ! Get coordinates in POS
      bxp = bpos*cos(azimuth)
      byp = bpos*sin(azimuth)
      bzp = blos

      ! Auxiliar
      by = (bxp*Mu + bzp*Mup)

      ! Get coordinates in vertical
      bx =  cp*by - byp*sp
      by =  sp*by + byp*cp
      bz = -bxp*Mup + bzp*Mu

      where (abs(bx).lt.TINYDP)
        bx =  0d0
      end where
      where (abs(by).lt.TINYDP)
        by =  0d0
      end where
      where (abs(bz).lt.TINYDP)
        bz =  0d0
      end where

      ! For each position
      do i=1,Num

        ! If there is a field
        if (bstrength(i).gt.0) then

          ! Compute angles
          btheta(i) = acos(bz(i)/bstrength(i))
          bphi(i) = atan2(by(i),bx(i))

          ! Round
          if (btheta(i).lt.TINYA) btheta(i) = 0d0
          if (btheta(i).gt.PI-TINYA) btheta(i) = PI

          ! If negative phi, put in 0,2pi
          if (bphi(i).lt.0d0) bphi(i) = bphi(i) + 2d0*PI
          if (bphi(i).lt.TINYA) bphi(i) = 0d0

        ! If no field
        else

          ! No angles
          btheta(i) = 0d0
          bphi(i) = 0d0

        end if ! Magnetic field

      end do ! Positions

      return

      end subroutine Bconversion

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform from vertical to LOS\n
      !!         Num(integer): Number of nodes\n
      !!           Mu(double): Cosine of the heliocentric angle\n
      !!          phi(double): Azimuth\n
      !! bstrength(double(:)): Magnetic field strength\n
      !!    btheta(double(:)): Magnetic field inclination\n
      !!      bphi(double(:)): Magnetic field azimuth\n
      !!      blos(double(:)): Magnetic field along the LOS\n
      !!      bpos(double(:)): Magnetic field in the POS\n
      !!   azimuth(double(:)): Magnetic azimuth in the POS
      subroutine B2Blos(Num, Mu, phi, bstrength, btheta, bphi, blos, &
                        bpos, azimuth)

      ! I/O
      integer, intent(in):: Num
      double precision, intent(in):: Mu, phi
      double precision, dimension(:), intent(in):: bstrength
      double precision, dimension(:), intent(in):: btheta
      double precision, dimension(:), intent(in):: bphi
      double precision, dimension(:), intent(inout):: blos
      double precision, dimension(:), intent(inout):: bpos
      double precision, dimension(:), intent(inout):: azimuth

      ! Local
      integer:: ii

      double precision:: Mup, cp, sp
      double precision, dimension(Num):: bx, by, bz, bxp, byp, st


      ! Get sin
      if (abs(Mu-1d0).lt.TINYA) then
        Mup = 0d0
      else
        Mup = sqrt(1d0 - Mu*Mu)
      end if

      ! Get cos-sin phi
      if (abs(phi).lt.TINYA) then
        cp =  1d0
        sp =  0d0
      else if (abs(phi-0.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp =  1d0
      else if (abs(phi-PI).lt.TINYA) then
        cp = -1d0
        sp =  0d0
      else if (abs(phi-1.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp = -1d0
      else
        cp = cos(phi)
        sp = sin(phi)
      end if

      ! Get cartesian
      st = sin(btheta)
      bx = bstrength*st*cos(bphi)
      by = bstrength*st*sin(bphi)
      bz = bstrength*cos(btheta)

      ! Rotate to LOS
      blos = (bx*cp + by*sp)
      bxp  =  Mu*blos - bz*Mup
      byp  = -bx*sp + by*cp
      blos = Mup*blos + bz*Mu

      where (abs(blos).lt.TINYDP)
        blos = 0d0
      end where
      where (abs(bxp).lt.TINYDP)
        bxp = 0d0
      end where
      where (abs(byp).lt.TINYDP)
        byp = 0d0
      end where

      ! Get POS
      bpos = sqrt(bxp*bxp + byp*byp)

      ! For each point
      do ii=1,Num

        if (abs(bxp(ii)).gt.0d0) then
          azimuth(ii) = atan2(byp(ii),bxp(ii))
        else if (abs(byp(ii)).gt.0d0) then
          if (byp(ii).gt.0d0) then
            azimuth(ii) = PI*0.5d0
          else
            azimuth(ii) = PI*1.5d0
          end if
        else
          azimuth(ii) = 0d0
        end if

      end do

      return

      end subroutine B2Blos

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform from LOS to vertical\n
      !!         Num(integer): Number of nodes\n
      !!           Mu(double): Cosine of the heliocentric angle\n
      !!          phi(double): Azimuth\n
      !!    vlos(double(:)): Velocity along the LOS\n
      !!    vpos(double(:)): Velocity in the POS\n
      !! azimuth(double(:)): Velocity azimuth in the POS\n
      !!      vx(double(:)): Velocity x component\n
      !!      vy(double(:)): Velocity y component\n
      !!      vz(double(:)): Velocity z component
      subroutine vconversion(Num,Mu,phi,vlos,vpos,azimuth,vx,vy,vz)

      ! IO
      integer, intent(in):: Num
      double precision, intent(in):: Mu,phi
      double precision, dimension(:),intent(in):: vlos
      double precision, dimension(:),intent(in):: vpos
      double precision, dimension(:),intent(in):: azimuth
      double precision, dimension(:),intent(inout):: vx,vy,vz

      ! Local
      double precision:: Mup, cp, sp
      double precision, dimension(Num):: vxp, vyp, vzp


      ! Get sin
      if (abs(Mu-1d0).lt.TINYA) then
        Mup = 0d0
      else
        Mup = sqrt(1d0 - Mu*Mu)
      end if

      ! Get cos-sin phi
      if (abs(phi).lt.TINYA) then
        cp =  1d0
        sp =  0d0
      else if (abs(phi-0.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp =  1d0
      else if (abs(phi-PI).lt.TINYA) then
        cp = -1d0
        sp =  0d0
      else if (abs(phi-1.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp = -1d0
      else
        cp = cos(phi)
        sp = sin(phi)
      end if

      ! Get coordinates in POS
      vxp = vpos*cos(azimuth)
      vyp = vpos*sin(azimuth)
      vzp = vlos

      ! Auxiliar
      vy = (vxp*Mu + vzp*Mup)

      ! Get coordinates in vertical
      vx =  cp*vy - vyp*sp
      vy =  sp*vy + vyp*cp
      vz = -vxp*Mup + vzp*Mu

      where (abs(vx).lt.TINYDP)
        vx =  0d0
      end where
      where (abs(vy).lt.TINYDP)
        vy =  0d0
      end where
      where (abs(vz).lt.TINYDP)
        vz =  0d0
      end where

      return

      end subroutine vconversion

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform from cartesian to LOS\n
      !!       Num(integer): Number of nodes\n
      !!         Mu(double): Cosine of the heliocentric angle\n
      !!        phi(double): Azimuth\n
      !!      vx(double(:)): Velocity x component\n
      !!      vy(double(:)): Velocity y component\n
      !!      vz(double(:)): Velocity z component\n
      !!    vlos(double(:)): Velocity along the LOS\n
      !!    vpos(double(:)): Velocity in the POS\n
      !! azimuth(double(:)): Velocity azimuth in the POS
      subroutine v2vlos(Num, Mu, phi, vx, vy, vz, vlos, &
                        vpos, azimuth)

      ! I/O
      integer, intent(in):: Num
      double precision, intent(in):: Mu, phi
      double precision, dimension(:), intent(in):: vx,vy,vz
      double precision, dimension(:), intent(inout):: vlos
      double precision, dimension(:), intent(inout):: vpos
      double precision, dimension(:), intent(inout):: azimuth

      ! Local
      integer:: ii

      double precision:: Mup, cp, sp
      double precision, dimension(Num):: vxp, vyp


      ! Get sin
      if (abs(Mu-1d0).lt.TINYA) then
        Mup = 0d0
      else
        Mup = sqrt(1d0 - Mu*Mu)
      end if

      ! Get cos-sin phi
      if (abs(phi).lt.TINYA) then
        cp =  1d0
        sp =  0d0
      else if (abs(phi-0.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp =  1d0
      else if (abs(phi-PI).lt.TINYA) then
        cp = -1d0
        sp =  0d0
      else if (abs(phi-1.5d0*PI).lt.TINYA) then
        cp =  0d0
        sp = -1d0
      else
        cp = cos(phi)
        sp = sin(phi)
      end if

      ! Rotate to LOS
      vlos = (vx*cp + vy*sp)
      vxp  =  Mu*vlos - vz*Mup
      vyp  = -vx*sp + vy*cp
      vlos = Mup*vlos + vz*Mu

      where (abs(vlos).lt.TINYDP)
        vlos = 0d0
      end where
      where (abs(vxp).lt.TINYDP)
        vxp = 0d0
      end where
      where (abs(vyp).lt.TINYDP)
        vyp = 0d0
      end where

      ! Get POS
      vpos = sqrt(vxp*vxp + vyp*vyp)

      ! For each point
      do ii=1,Num

        if (abs(vxp(ii)).gt.0d0) then
          azimuth(ii) = atan2(vyp(ii),vxp(ii))
        else if (abs(vyp(ii)).gt.0d0) then
          if (vyp(ii).gt.0d0) then
            azimuth(ii) = PI*0.5d0
          else
            azimuth(ii) = PI*1.5d0
          end if
        else
          azimuth(ii) = 0d0
        end if

      end do

      return

      end subroutine v2vlos

!#####################################################################
!#####################################################################
!#####################################################################

      end module model_mod
