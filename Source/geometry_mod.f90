      !> Get envelope geometry
      module geometry_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     01/10/2022
!  Last version:
!     03/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  getlocalframe
!    Get the geometrical quantities for the local reference frame in
!  a node of the CLE calculation
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only: pi , TINYANG , TINYB , TINYVEL
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get the geometrical quantities for the local reference frame
      !! in a node of the CLE calculation\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Bfield(Bfield_class): Structure with magnetic field
      !!                           data\n
      !!                x(double): Coordinate along the LOS measured
      !!                           from the solar center\n
      !!            mode(integer): Type of input atmospheric model\n
      !! GeomP(Coronapoint_class): Structure with geometric data for a
      !!                           CLE node
      subroutine getlocalframe(Atmo,Bfield,x,mode,GeomP)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
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

      ! Not counting RAM because it is anecdotal
      allocate(Bfield%Bstrength(1))
      allocate(Bfield%Btheta(1))
      allocate(Bfield%Bphi(1))


      !
      ! 3D model (cartesian or not)
      !
      if (mode.eq.0.or.mode.eq.2) then

        ! Position in plane of the sky
        y = Atmo%ypos
        z = Atmo%zpos

        ! Radius from solar center
        r = sqrt(x*x + y*y + z*z)

        ! Sky radius
        rk = sqrt(y*y + z*z)

        ! Height over the surface
        h = r - 1d0

        ! Theta (heliocentric angle)
        GeomP%theta = acos(z/r)

        !
        ! Phi
        !

        ! If x is basically zero
        if (abs(x).lt.TINYANG) then

          ! Positive y
          if (y.gt.0) then

            ! Phi must be 90º
            GeomP%phi = .5d0*pi

          ! Negative y
          else

            ! Phi must be 270º
            GeomP%phi = 1.5d0*pi

          end if ! y sign

        ! x is not zero
        else

          ! If y is basically zero
          if (abs(y).lt.TINYANG) then

            ! Positive x
            if (x.gt.0) then

              ! Phi must be 0º
              GeomP%phi = 0d0

            ! Negative x
            else

              ! Phi must be 180º
              GeomP%phi = pi

            end if ! x sign

          ! Both x and y are different from 0
          else

            ! Get phi from arctan
            GeomP%phi = atan(y/x)

            ! If negative x
            if (x.lt.0d0) then

              ! Put phi in correct quadrant
              GeomP%phi = GeomP%phi + pi

            ! If negative y
            else if (y.lt.0d0) then

              ! Put phi in correct quadrant
              GeomP%phi = GeomP%phi + 2d0*pi

            end if ! y sign
          end if ! y magnitude
        end if ! x magnitude

        !
        ! Get coordinates of the LOS and the magnetic field in the
        ! local frame
        !

        ! alpha y beta are from plasma.f from Roberto
        alpha = - atan2(x,rk)
        beta = atan2(y,z)

        ! Get theta
        GeomP%geom(1) = .5d0*pi + alpha

        ! Get phi
        GeomP%geom(2) = 0d0

        ! Get gamma
        GeomP%geom(3) = beta + .5d0*pi

        ! Cosine and sine of alpha
        GeomP%CA = cos(alpha)
        GeomP%SA = sin(alpha)

        ! Cosine and sine of beta
        GeomP%CB = cos(beta)
        GeomP%SB = sin(beta)

        !
        ! Get velocity field in local frame
        !

        ! Apply rotation matrix
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

          ! Apply rotation matrix
          bxl =   GeomP%CA*Atmo%bx(1) + &
                  GeomP%SA*GeomP%SB*Atmo%by(1) + &
                  GeomP%SA*GeomP%CB*Atmo%bz(1)
          byl =   GeomP%CB*Atmo%by(1) - &
                  GeomP%SB*Atmo%bz(1)
          bzl = - GeomP%SA*Atmo%bx(1) + &
                  GeomP%CA*GeomP%SB*Atmo%by(1) + &
                  GeomP%CA*GeomP%CB*Atmo%bz(1)

          !
          ! ThetaB
          !

          ! If module is smaller than z component
          if (bl(1).le.abs(bzl)) then

            ! Longitudinal
            bl(2) = 0d0

          ! Module is larger than z component
          else

            ! We can compute arccos
            bl(2) = acos(bzl/bl(1))

          end if ! Relation between module and z component

          !
          ! PhiB
          !

          ! If the x component is too small
          if (abs(bxl).lt.TINYANG) then

            ! Positive y component
            if (byl.gt.0) then

              ! PhiB must be 90º
              bl(3) = .5d0*pi

            ! Negative y component
            else

              ! PhiB must be 270º
              bl(3) = 1.5d0*pi

            end if ! y component sign

          ! If the x component is NOT small
          else

            ! If the y component is small
            if (abs(byl).lt.TINYANG) then

              ! Positive x component
              if (bxl.gt.0) then

                ! PhiB must be 0º
                bl(3) = 0d0

              ! Negative x component
              else

                ! PhiB must be 180º
                bl(3) = pi

              end if ! x component sign

            ! If the y component is NOT small
            else

              ! Get PhiB from arctan
              bl(3) = atan(byl/bxl)

              ! If negative x component
              if (bxl.lt.0d0) then

                ! Put in correct quadrant
                bl(3) = bl(3) + pi

              ! If negative y component
              else if (byl.lt.0d0) then

                ! Put in correct quadrant
                bl(3) = bl(3) + 2d0*pi

              end if ! x and y component signs
            end if ! x component magnitude
          end if ! y component magnitude

        ! There is no magnetic field
        else

          ! Angles are zero
          bl(2) = 0d0
          bl(3) = 0d0

        end if

        ! Save magnetic field
        Bfield%Bstrength = bl(1)
        Bfield%Btheta = bl(2)
        Bfield%Bphi = bl(3)

        ! Save velocity
        Atmo%vx = vl(1)
        Atmo%vy = vl(2)
        Atmo%vz = vl(3)

      ! Slab model
      else if (mode.eq.1) then

        ! Height
        h = x

        ! Radius, measured from center
        r = h + 1d0

        ! Heliocentric angle and azimuth
        GeomP%theta = Atmo%ypos
        GeomP%phi = 0d0

        !
        ! Get coordinates of the LOS and the magnetic field in the
        ! local frame
        !

        ! Get theta
        GeomP%geom(1) = GeomP%theta

        ! Get phi
        GeomP%geom(2) = GeomP%phi

        ! Get gamma
        GeomP%geom(3) = .5d0*pi

        ! Get alpha and beta angles
        alpha = GeomP%geom(1) - .5d0*pi
        beta = GeomP%geom(3) - .5d0*pi

        ! Get cosine and sine of alpha
        GeomP%CA = cos(alpha)
        GeomP%SA = sin(alpha)

        ! Get cosine and sine of beta
        GeomP%CB = cos(beta)
        GeomP%SB = sin(beta)

        ! Save magnetic field
        Bfield%Bstrength = Atmo%bx(1)
        Bfield%Btheta = Atmo%by(1)
        Bfield%Bphi = Atmo%bz(1)

        ! Get velocity in spherical
        vl(1) = Atmo%vx(1)
        vl(2) = Atmo%vy(1)
        vl(3) = Atmo%vz(1)

        ! Transform to cartesian
        Atmo%vx(1) = vl(1)*sin(vl(2))*cos(vl(3))
        Atmo%vy(1) = vl(1)*sin(vl(2))*sin(vl(3))
        Atmo%vz(1) = vl(1)*cos(vl(2))

      end if ! Type of atmospheric model

      !
      ! Define if dynamic or axial
      !

      ! If the velocity is small
      if (sqrt(Atmo%vx(1)*Atmo%vx(1) + &
               Atmo%vy(1)*Atmo%vy(1) + &
               Atmo%vz(1)*Atmo%vz(1)).le.TINYVEL) then

        ! The problem is axial and static
        axial = .True.
        dyn = .False.

      ! If the transversal velocity is small
      else if ((abs(Atmo%vy(1)) + abs(Atmo%vz(1))).le.0d0) then

        ! The problem is axial but dynamic
        axial = .True.
        dyn = .True.

      ! If all components
      else

        ! The problem is not axial and is dynamic
        axial = .False.
        dyn = .True.

      end if ! Velocity amplitude and direction


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
