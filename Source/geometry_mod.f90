      !> Get envelope geometry
      module geometry_mod
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
!     10/04/2024 V3.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/04/2024:    V3.0.1 - The height in the slab models is now
!                             expected to not include the stellar
!                             radius (TdPA)
!                           - The slab model now expects the velocity
!                             in polar coordinates in the local
!                             reference frame (TdPA)
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
!    Get geometry for the local frame for a CLE point
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only: pi , TINYANG , TINYB
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return the RT coefficients for a given point along the LOS
      !! for a CLE calculation
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Bfield(Bfield_class): Structure with magnetic field
      !!                           data\n
      !!                x(dfloat): Coordinate along the LOS\n
      !!            mode(integer): Type of input atmospheric model\n
      !! GeomP(Coronapoint_class): Geometric data for a point in
      !!                           CLE
      subroutine getlocalframe(Atmo,Bfield,x,mode,GeomP)

      ! I/O
      type(Atmo_class):: Atmo
      type(Bfield_class), intent(out):: Bfield
      type(Coronapoint_class), intent(out):: GeomP
      integer, intent(in):: mode
      double precision, intent(in):: x


      ! Local
      double precision:: y,z
      double precision:: r,rk,h,alpha,beta,SY2,CY2,CY3,bxl,byl,bzl
      double precision, dimension(3):: vl,bl


      !
      ! Allocate magnetic field
      !
      allocate(Bfield%Bstrength(1))
      allocate(Bfield%Btheta(1))
      allocate(Bfield%Bphi(1))


      !
      ! 3D model
      !
      if (mode.eq.0.or.mode.eq.2) then

        ! Position
        y = Atmo%ypos
        z = Atmo%zpos

        ! Radius
        r = sqrt(x*x + y*y + z*z)

        ! Sky radius
        rk = sqrt(y*y + z*z)

        ! Height over the surface
        h = r - 1d0

        ! Theta
        GeomP%theta = acos(z/r)

        ! Phi
        if (abs(x).lt.TINYANG) then
          if (y.gt.0) then
            GeomP%phi = .5d0*pi
          else
            GeomP%phi = 1.5d0*pi
          end if
        else
          if (abs(y).lt.TINYANG) then
            if (x.gt.0) then
              GeomP%phi = 0d0
            else
              GeomP%phi = pi
            end if
          else
            GeomP%phi = atan(y/x)
            if (x.lt.0d0) then
              GeomP%phi = GeomP%phi + pi
            else if (y.lt.0d0) then
              GeomP%phi = GeomP%phi + 2d0*pi
            end if
          end if
        end if

        !
        ! Get coordinates of the LOS and the magnetic field in the
        ! local frame
        !
        alpha = - atan2(x,rk)           ! Copying plasma.f from Roberto
        beta = atan2(y,z)               ! Copying plasma.f from Roberto
        GeomP%geom(1) = .5d0*pi + alpha ! theta
        GeomP%geom(2) = 0d0             ! phi
        GeomP%geom(3) = beta + .5d0*pi  ! gamma

        ! Cosines and sines
        GeomP%CA = cos(alpha)
        GeomP%SA = sin(alpha)
        GeomP%CB = cos(beta)
        GeomP%SB = sin(beta)

        !
        ! Get velocity field in local frame
        !
        vl(1) =   GeomP%CA*Atmo%vx(1) + &
                  GeomP%SA*GeomP%SB*Atmo%vy(1) + &
                  GeomP%SA*GeomP%CB*Atmo%vz(1)
        vl(2) =   GeomP%CB*Atmo%vy(1) - &
                  GeomP%SB*Atmo%vz(1)
        vl(3) = - GeomP%SA*Atmo%vx(1) + &
                  GeomP%CA*GeomP%SB*Atmo%vy(1) + &
                  GeomP%CA*GeomP%CB*Atmo%vz(1)

        !
        ! Magnetic field
        !

        ! Strength
        bl(1) = sqrt(Atmo%bx(1)*Atmo%bx(1) + &
                     Atmo%by(1)*Atmo%by(1) + &
                     Atmo%bz(1)*Atmo%bz(1))

        ! If there is a magnetic field
        if (bl(1).gt.0d0) then

          !
          ! Get magnetic field in local frame
          !
          bxl =   GeomP%CA*Atmo%bx(1) + &
                  GeomP%SA*GeomP%SB*Atmo%by(1) + &
                  GeomP%SA*GeomP%CB*Atmo%bz(1)
          byl =   GeomP%CB*Atmo%by(1) - &
                  GeomP%SB*Atmo%bz(1)
          bzl = - GeomP%SA*Atmo%bx(1) + &
                  GeomP%CA*GeomP%SB*Atmo%by(1) + &
                  GeomP%CA*GeomP%CB*Atmo%bz(1)

          ! ThetaB
          if (bl(1).lt.abs(bzl)) then
            bl(2) = 0d0
          else
            bl(2) = acos(bzl/bl(1))
          end if

          ! PhiB
          if (abs(bxl).lt.TINYANG) then
            if (byl.gt.0) then
              bl(3) = .5d0*pi
            else
              bl(3) = 1.5d0*pi
            end if
          else
            if (abs(byl).lt.TINYANG) then
              if (bxl.gt.0) then
                bl(3) = 0d0
              else
                bl(3) = pi
              end if
            else
              bl(3) = atan(byl/bxl)
              if (bxl.lt.0d0) then
                bl(3) = bl(3) + pi
              else if (byl.lt.0d0) then
                bl(3) = bl(3) + 2d0*pi
              end if
            end if
          end if

        ! There is no magnetic field
        else

          bl(2) = 0d0
          bl(3) = 0d0

        end if

        Bfield%Bstrength = bl(1)
        Bfield%Btheta = bl(2)
        Bfield%Bphi = bl(3)
        Atmo%vx = vl(1)
        Atmo%vy = vl(2)
        Atmo%vz = vl(3)

      ! Slab model
      else if (mode.eq.1) then

        h = x
        r = h + 1d0
        GeomP%theta = Atmo%ypos
        GeomP%phi = 0d0

        !
        ! Get coordinates of the LOS and the magnetic field in the
        ! local frame
        !
        GeomP%geom(1) = GeomP%theta    ! theta
        GeomP%geom(2) = GeomP%phi      ! phi
        GeomP%geom(3) = .5d0*pi        ! gamma
        alpha = GeomP%geom(1) - .5d0*pi
        beta = GeomP%geom(3) - .5d0*pi

        GeomP%CA = cos(alpha)
        GeomP%SA = sin(alpha)
        GeomP%CB = cos(beta)
        GeomP%SB = sin(beta)

        Bfield%Bstrength = Atmo%bx(1)
        Bfield%Btheta = Atmo%by(1)
        Bfield%Bphi = Atmo%bz(1)

        vl(1) = Atmo%vx(1)
        vl(2) = Atmo%vy(1)
        vl(3) = Atmo%vz(1)

        Atmo%vx(1) = vl(1)*sin(vl(2))*cos(vl(3))
        Atmo%vy(1) = vl(1)*sin(vl(2))*sin(vl(3))
        Atmo%vz(1) = vl(1)*cos(vl(2))

      end if

      !
      ! Define if dynamic or axial
      !
      if ((abs(Atmo%vx(1)) + abs(Atmo%vy(1)) + &
           abs(Atmo%vz(1))).le.0d0) then
        axial = .True.
        dyn = .False.
      else if ((abs(Atmo%vy(1)) + abs(Atmo%vz(1))).le.0d0) then
        axial = .True.
        dyn = .True.
      else
        axial = .False.
        dyn = .True.
      end if


      !
      ! Get the CLV constants given the heights
      !

      ! Gamma angle (maximum angle from the point that intersects the
      ! solar surface
      GeomP%SY = 1d0/r
      SY2 = GeomP%SY*GeomP%SY
      GeomP%CY = sqrt(1d0 - SY2)
      CY2 = GeomP%CY*GeomP%CY
      CY3 = CY2*GeomP%CY

      ! CLV constants
      GeomP%CLV(1) = 1d0 - GeomP%CY
      GeomP%CLV(2) = GeomP%CY - 0.5d0 - 0.5d0*CY2* &
                     log((1d0+GeomP%SY)/GeomP%CY)/GeomP%SY
      GeomP%CLV(3) = (GeomP%CY + 2d0)*(GeomP%CY - 1d0)/ &
                     3d0/(GeomP%CY + 1d0)
      GeomP%CLV(4) = GeomP%CY*SY2
      GeomP%CLV(5) = 0.125d0*(8d0*CY3 - 3d0*CY2 - 8d0*GeomP%CY + &
                     2d0 + (4d0 - 3d0*CY2)*CY2* &
                     log((1d0+GeomP%SY)/GeomP%CY)/GeomP%SY)
      GeomP%CLV(6) = (GeomP%CY - 1d0)* &
                     (9d0*CY3 + 18d0*CY2 + 7d0*GeomP%CY - 4d0)/ &
                     15d0/(GeomP%CY + 1d0)

      end subroutine getlocalframe

!#####################################################################
!#####################################################################
!#####################################################################

      end module geometry_mod
