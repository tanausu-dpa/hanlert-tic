      !> Interpolation routines
      module inter_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Hao Li (IAC)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     09/08/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/08/2023:    V3.0.2 - Change the control value in subroutine
!                             Intpol_Quadratic_Bezier (HL)
!                           - Do Cubic Bezier interpolation if node
!                             number is larger than 2 (HL)
!
!     03/08/2023:    V3.0.1 - Added Intpol, Intpol_Cubic_Hermite,
!                             Intpol_Quadratic_Bezier,
!                             Intpol_Cubic_Bezier_C,
!                             Intpol_Quadratic_Bezier_C,
!                             Intpol_Lin, Intpol_Lin_mD, parabolic,
!                             and DERIV from the TIC@inter_mod.f90
!                             module (TdPA)
!                           - Added Intpol_Lin_stk used to interpolate
!                             Stokes profiles, avoiding extra
!                             calculations or generating an additional
!                             copy of the input array in the intensity
!                             case (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     07/23/2019:    V1.1.3 - In colinter, removed the scaling with
!                             the x variable and added option to
!                             force linear (TdPA)
!
!     05/22/2019:    V1.1.2 - Bugfix: In bilinear, one of the checks
!                             for the boundaries was wrong, it said
!                             z1.ge.x2(1), and it must be z2.ge.x2(1)
!                             instead (TdPA)
!
!     03/13/2019:    V1.1.1 - Now the warning about interpolation is
!                             given outside (TdPA)
!                           - Do not need aborted_mod (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     09/08/2017:    V1.0.1 - If the splines give negative collisional
!                             rates, it does linear instead (TdPA)
!                           - Had to add access to commons_mod (TdPA)
!                           - Removed dE from colinter call (TdPA)
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
!    Contains interpolation routines
!
! colinter:
!   Interpolates inelastic collisions
!
! bilinear:
!   Does a bilinear interpolation for a single output. Out of
!   boundaries is made constant
!
! linear:
!   Does a linear interpolation for a single output. Out of
!   boundaries is made constant
!
! spline:
!   Calculate the coefficients for cubic spline interpolation
!
! ispline:
!   Evaluates cubic spline interpolation
!   Alex G: January 2010
!
! spline_2d:
!   Bicubic spline interpolation
!
! Intpol:
!   Interpolation manager for the inversion
!
! Intpol_Cubic_Hermite:
!   Cubic hermite interpolation for the inversion
!
! Intpol_Quadratic_Bezier:
!   Quadratic bezier interpolation for the inversion
!
! Intpol_Cubic_Bezier_C:
!   Cubic bezier interpolation for the inversion, with
!   coefficients
!
! Intpol_Quadratic_Bezier_C:
!   Get Quadratic bezier interpolation for the inversion,
!   with coefficients
!
! Intpol_Lin:
!   Linear interpolation for the inversion
!
! Intpol_Lin_stk:
!   Linear interpolation for the inversion in Stokes parameters
!
! Intpol_Lin_mD:
!   Linear interpolation for the inversion for several variables
!
! Parabolic:
!   Gives X value for the minimum of Y in parabolic interpolation
!
! DERIV:
!   Calculate derivatives for interpolation for the inversion
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolator for collisional rates\n
      !!    xin(dfloat(:)): Input temperature axis\n
      !!    yin(dfloat(:)): Input collisional rates\n
      !!      nin(integer): Size of xin and yin\n
      !!   xout(dfloat(:)): Output temperature axis\n
      !!   yout(dfloat(:)): Interpolated rates\n
      !!     nout(integer): Size of xout and yout\n
      !!  flinear(logical): Bool to indicate if the interpolation
      !!                    must be linear\n
      !!  llinear(logical): Bool the indicates if the interpolation
      !!                    had to be linear
      subroutine colinter(xin,yin,nin,xout,yout,nout,flinear,llinear)

      ! I/O

      logical,intent(in):: flinear
      logical,intent(inout):: llinear
      integer, intent(in):: nin, nout
      double precision, dimension(:):: xin, yin, xout, yout

      ! Local

      integer:: iout

      double precision, dimension(nin):: b, c, d

      ! Checks if we were forced to use linear when splines
      ! were requested
      llinear = .False.

      ! If only one input value
      if (nin.eq.1) then

        yout = yin(1)

      ! If there are only two, it has to be linear
      else if (nin.eq.2) then

        do iout=1,nout
          call linear(xin,yin,xout(iout),yout(iout))
        end do

      ! If forcing linear interpolation
      else if (flinear) then

        do iout=1,nout
          call linear(xin,yin,xout(iout),yout(iout))
        end do

      ! If not, use splines
      else

        call spline(xin,yin,b,c,d,nin)
        do iout=1,nout

          yout(iout) = ispline(xout(iout),xin,yin,b,c,d,nin)

          if (yout(iout).lt.0) then
            llinear = .True.
            exit
          end if

        end do

        if (llinear) then

          do iout=1,nout
            call linear(xin,yin,xout(iout),yout(iout))
          end do

        end if

      end if

      end subroutine colinter

!#####################################################################
!#####################################################################
!#####################################################################

      !> Bilinear interpolator\n
      !!     x1(dfloat(:)): Input coordinates axis 1\n
      !!     x2(dfloat(:)): Input coordinates axis 2\n
      !!      y(dfloat(:)): Data to interpolate\n
      !!     z1(dfloat(:)): Output coordinates axis 1\n
      !!     z2(dfloat(:)): Output coordinates axis 2\n
      !!   yout(dfloat(:)): Interpolated data
      subroutine bilinear(x1,x2,y,z1,z2,yout)

      ! I/O

      double precision, intent(in):: z1,z2
      double precision, intent(out):: yout
      double precision, dimension(:), intent(in):: x1, x2
      double precision, dimension(:,:), intent(in):: y

      ! Local

      integer:: n1,n2,ii,i10,i11,i20,i21

      double precision:: x10,x11,x20,x21,yA,yB,dx1,dx2

      n1 = size(x1)
      n2 = size(x2)

      ! Within bounds
      if (z1.ge.x1(1).and.z1.le.x1(n1).and. &
          z2.ge.x2(1).and.z2.le.x2(n2)) then

        do ii=1,n1-1
          if (z1.ge.x1(ii).and.z1.le.x1(ii+1)) then
            i10 = ii
            i11 = ii+1
            x10 = x1(i10)
            x11 = x1(i11)
            dx1 = 1d0/(x11 - x10)
            exit
          end if
        end do
        do ii=1,n2-1
          if (z2.ge.x2(ii).and.z2.le.x2(ii+1)) then
            i20 = ii
            i21 = ii+1
            x20 = x2(i20)
            x21 = x2(i21)
            dx2 = 1d0/(x21 - x20)
            exit
          end if
        end do

        yA = ((y(i10,i21) - y(i10,i20))*z2 + &
              y(i10,i20)*x21 - y(i10,i21)*x20)*dx2
        yB = ((y(i11,i21) - y(i11,i20))*z2 + &
              y(i11,i20)*x21 - y(i11,i21)*x20)*dx2
        yout = ((yB - yA)*z1 + yA*x11 - yB*x10)*dx1

      ! Out of x2 bounds
      else if (z1.ge.x1(1).and.z1.le.x1(n1)) then

        do ii=1,n1-1
          if (z1.ge.x1(ii).and.z1.le.x1(ii+1)) then
            i10 = ii
            i11 = ii+1
            x10 = x1(i10)
            x11 = x1(i11)
            dx1 = 1d0/(x11 - x10)
            exit
          end if
        end do

        if (z2.lt.x2(1)) then

          i20 = 1

        else if (z2.gt.x2(n2)) then

          i20 = n2

        end if

        yout = ((y(i11,i20) - y(i10,i20))*z1 + &
                y(i10,i20)*x11 - y(i11,i20)*x10)*dx1

      ! Out of x1 bounds
      else if (z2.ge.x2(1).and.z2.le.x2(n2)) then

        do ii=1,n2-1
          if (z2.ge.x2(ii).and.z2.le.x2(ii+1)) then
            i20 = ii
            i21 = ii+1
            x20 = x2(i20)
            x21 = x2(i21)
            dx2 = 1d0/(x21 - x20)
            exit
          end if
        end do

        if (z1.lt.x1(1)) then

          i10 = 1

        else if (z1.gt.x1(n1)) then

          i10 = n1

        end if

        yout = ((y(i10,i21) - y(i10,i20))*z2 + &
                y(i10,i20)*x21 - y(i10,i21)*x20)*dx2

      ! Completely out of bounds
      else

        if (z1.lt.x1(1)) then

          i10 = 1

        else if (z1.gt.x1(n1)) then

          i10 = n1

        end if

        if (z2.lt.x2(1)) then

          i20 = 1

        else if (z2.gt.x2(n2)) then

          i20 = n2

        end if

        yout = y(i10,i20)

      end if

      end subroutine bilinear

!#####################################################################
!#####################################################################
!#####################################################################

      !> Linear interpolator\n
      !!    x(dfloat(:)): Input x axis\n
      !!    y(dfloat(:)): Input y axis\n
      !!       z(dfloat): Output x value\n
      !!    yout(dfloat): Interpolated value
      subroutine linear(x,y,z,yout)

      ! I/O

      double precision, intent(in):: z
      double precision, intent(out):: yout
      double precision, dimension(:), intent(in):: x,y

      ! Local

      integer:: n,ii,i0,i1

      double precision:: y0,y1,x0,x1,dx

      n = size(x)

      ! Out of lower bound
      if (z.le.x(1)) then

        yout = y(1)

      ! Out of upper bound
      else if (z.ge.x(n)) then

        yout = y(n)

      ! Between bounds
      else

        do ii=1,n-1
          if (z.ge.x(ii).and.z.le.x(ii+1)) then
            i0 = ii
            i1 = ii+1
            x0 = x(i0)
            x1 = x(i1)
            y0 = y(i0)
            y1 = y(i1)
            dx = 1d0/(x1 - x0)
            exit
          end if
        end do

        yout = ((y1 - y0)*z + y0*x1 - y1*x0)*dx

      end if

      end subroutine linear

!#####################################################################
!#####################################################################
!#####################################################################

      !> Cubic spline interpolator. Gets the interpolation
      !! coefficients\n
      !!    x(dfloat(:)): Input x axis\n
      !!    y(dfloat(:)): Input y axis\n
      !!    b(dfloat(:)): Coefficient of spline interpolation\n
      !!    c(dfloat(:)): Coefficient of spline interpolation\n
      !!    d(dfloat(:)): Coefficient of spline interpolation\n
      !!      n(integer): Size of x,y,b,c,d
      subroutine spline(x,y,b,c,d,n)

      ! I/O

      integer:: n
      double precision, dimension(:):: x, y, b, c, d

      ! Local

      integer:: i,j,k

      double precision:: a

      k = n - 1

      ! Need at least three points
      if(n.lt.2)return

      ! If just three points, linear interpolation
      if(n.lt.3)then
        b(1) = (y(2)-y(1))/(x(2)-x(1))
        c(1) = 0d0
        d(1) = 0d0
        b(2) = b(1)
        c(2) = 0d0
        d(2) = 0d0
        return
      end if

      ! Prepare
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i=2,k
        d(i) = x(i+1) - x(i)
        b(i) = 2d0*(d(i-1) + d(i))
        c(i+1) = (y(i+1) - y(i))/d(i)
        c(i) = c(i+1) - c(i)
      end do

      ! End conditions
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0d0
      c(n) = 0d0
      if(n.ne.3)then
        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)*d(1)/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      end if

      ! Forward elimination
      do i=2,n
        a = d(i-1)/b(i-1)
        b(i) = b(i) - a*d(i-1)
        c(i) = c(i) - a*c(i-1)
      end do

      ! Back substitution
      c(n) = c(n)/b(n)
      do j=1,k
        i = n-j
        c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do

      ! Coefficients
      b(n) = (y(n) - y(k))/d(k) + d(k)*(c(k) + 2d0*c(n))
      do i=1,k
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3d0*c(i)
      end do
      c(n) = 3d0*c(n)
      d(n) = d(n-1)

      end subroutine spline

!####################################################################
!####################################################################
!####################################################################

      !> Evaluates a cubic spline interpolation\n
      !!    u(dfloat(:)): Output x value\n
      !!    x(dfloat(:)): Input x axis\n
      !!    y(dfloat(:)): Input y axis\n
      !!    b(dfloat(:)): Coefficient of spline interpolation\n
      !!    c(dfloat(:)): Coefficient of spline interpolation\n
      !!    d(dfloat(:)): Coefficient of spline interpolation\n
      !!      n(integer): Size of x,y,b,c,d
      double precision function ispline(u,x,y,b,c,d,n)

      ! I/O

      integer:: n
      double precision:: u
      double precision, dimension(:):: x, y, b, c, d

      ! Local

      integer:: i,j,k

      double precision:: dx

      ! Constant if out of boundaries
      if(u.le.x(1))then
        ispline = y(1)
        return
      end if
      if(u.ge.x(n))then
        ispline = y(n)
        return
      end if

      ! Binary search for for i, such that x(i) <= u <= x(i+1)
      i = 1
      j = n + 1
      do while (j.gt.i+1)
        k = (i+j)/2
        if(u.lt.x(k))then
          j=k
        else
          i=k
        end if
      end do

      ! Evaluate spline interpolation
      dx = u - x(i)
      ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

      end function ispline

!####################################################################
!####################################################################
!####################################################################

      !> 2D cubic splines interpolator\n
      !!     x1(dfloat(:)): Input coordinates axis 1\n
      !!     x2(dfloat(:)): Input coordinates axis 2\n
      !!      y(dfloat(:)): Data to interpolate\n
      !!       n1(integer): Size of x1 and second column y\n
      !!       n2(integet): Size of x2 and first column y\n
      !!     z1(dfloat(:)): Output coordinates axis 1\n
      !!     z2(dfloat(:)): Output coordinates axis 2\n
      !!   yout(dfloat(:)): Interpolated data
      subroutine spline_2d(x1,x2,y,n1,n2,z1,z2,yout)

      ! I/O

      integer, intent(in):: n1, n2
      double precision, intent(in):: z1, z2
      double precision, intent(out):: yout
      double precision, dimension(:), intent(in):: x1
      double precision, dimension(:), intent(in):: x2
      double precision, dimension(:,:), intent(in):: y ! (n2,n1)

      ! Local

      integer:: i1,i2,iini,ifin,imed

      double precision:: d1,d2,dx1,dx2
      double precision, dimension(n1):: yslice1,y2slice1,V1
      double precision, dimension(n2):: yslice2,y2slice2,V2


      do i1=1,n1

        ! Take a slice of y in the first dimension
        yslice2 = y(:,i1)

        ! Calculate derivative of the slice
        y2slice2(1) = 0d0
        y2slice2(n2) = 0d0
        V2(1) = 0d0
        do i2=2,n2-1
          d1 = (x2(i2) - x2(i2-1))/(x2(i2+1) - x2(i2-1))
          d2 = d1*y2slice2(i2-1) + 2d0
          y2slice2(i2) = (d1 - 1d0)/d2
          V2(i2) = (6d0*((yslice2(i2+1) - yslice2(i2))/ &
                         (x2(i2+1) - x2(i2)) - &
                         (yslice2(i2) - yslice2(i2-1))/ &
                         (x2(i2) - x2(i2-1)))/ &
                        (x2(i2+1) - x2(i2-1)) - d1*V2(i2-1))/d2
        end do
        do i2=n2-1,1,-1
          y2slice2(i2) = y2slice2(i2)*y2slice2(i2+1) + V2(i2)
        end do

        ! Get the part of the interpolation in axis 2
        ! Find the neighbours by split search
        iini = 1
        ifin = n2
        do while ((ifin-iini).gt.1)
          imed = (iini + ifin)/2
          if (x2(imed).gt.z2) then
            ifin = imed
          else
            iini = imed
          end if
        end do
        ! Interpolate
        dx2 = x2(ifin) - x2(iini)
        d1 = (x2(ifin) - z2)/dx2
        d2 = (z2 - x2(iini))/dx2
        yslice1(i1) = d1*yslice2(iini) + d2*yslice2(ifin) + &
                      ((d1*(d1*d1 - 1d0))*y2slice2(iini) + &
                       (d2*(d2*d2 - 1d0))*y2slice2(ifin))* &
                      dx2*dx2/6d0
      end do

      ! Calculate derivative on the other direction
      y2slice1(1) = 0d0
      y2slice1(n1) = 0d0
      V1(1) = 0d0
      do i1=2,n1-1
        d1 = (x1(i1) - x1(i1-1))/(x1(i1+1) - x1(i1-1))
        d2 = d1*y2slice1(i1-1) + 2d0
        y2slice1(i1) = (d1 - 1d0)/d2
        V1(i1) = (6d0*((yslice1(i1+1) - yslice1(i1))/ &
                       (x1(i1+1) - x1(i1)) - &
                       (yslice1(i1) - yslice1(i1-1))/ &
                       (x1(i1) - x1(i1-1)))/ &
                      (x1(i1+1) - x1(i1-1)) -d1*V1(i1-1))/d2
      end do
      do i1=n1-1,1,-1
        y2slice1(i1) = y2slice1(i1)*y2slice1(i1+1) + V1(i1)
      end do

      ! Interpolate in the other dimension
      ! Find the neighbours by split search
      iini = 1
      ifin = n1
      do while ((ifin-iini).gt.1)
        imed = (iini + ifin)/2
        if (x1(imed).gt.z1) then
          ifin = imed
        else
          iini = imed
        end if
      end do
      ! Interpolate
      dx1 = x1(ifin) - x1(iini)
      d1 = (x1(ifin) - z1)/dx1
      d2 = (z1 - x1(iini))/dx1
      yout = d1*yslice1(iini) + d2*yslice1(ifin) + &
             ((d1*(d1*d1 - 1d0))*y2slice1(iini) + &
              (d2*(d2*d2 - 1d0))*y2slice1(ifin))* &
             (dx1*dx1/6d0)

      end subroutine spline_2d

!####################################################################
!####################################################################
!####################################################################

      !> Manages the interpolation of arrays in the inversion\n
      !!           x(double(:)): Input x axis\n
      !!           y(double(:)): Input y axis\n
      !!             n(integer): Size of x and y arrays\n
      !!          xx(double(:)): Output x axis\n
      !!          yy(double(:)): Output y axis\n
      !!            nn(integer): Size of xx and yy arrays\n
      !!   Indx_Intpol(integer): Type of interpolation:\n
      !!                          0 linear\n
      !!                          1 quadratic bezier\n
      !!                          2 cubic bezier\n
      !! Indx_Extrapol(integer): Type of extrapolation
      subroutine Intpol(x,y,n,xx,yy,nn,Indx_Intpol,Indx_Extrapol)

      ! IO
      integer, intent(in):: n, nn, Indx_Intpol, Indx_Extrapol
      double precision, dimension(n), intent(in):: x, y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local
      logical:: Increased

      integer:: i, mini, maxi

      double precision:: a0, b0, a1, b1, xmin, xmax, ymin, ymax


      ! Choose the interpolation
      select case(Indx_Intpol)

        ! Linear
        case(0)

          call Intpol_Lin(x, y, n, xx, yy, nn)

        ! Quadratic bezier
        case(1)

          ! Only if enough elements
          if(n.gt.2) then
            call Intpol_Quadratic_Bezier(x, y, n, xx, yy, nn)
          else
            call Intpol_Lin(x, y, n, xx, yy, nn)
          endif

        ! Cubic bezier
        case(2)

          ! Only if enough elements
          if(n.gt.2) then
            call Intpol_Cubic_Hermite(x, y, n, xx, yy, nn)
          else
            call Intpol_Lin(x, y, n, xx, yy, nn)
          endif

        case default

          yy = 0d0

      end select ! Type of interpolation

      ! If dimension was 1, return already
      if(n.eq.1) return

      ! If not extrapolating, return
      if (Indx_Extrapol.eq.3) return

      ! Initialize
      mini = 0d0
      maxi = 0d0

      ! Determine ordering
      if (x(n).gt.x(1)) then
        increased = .True.
        xmin = x(1)
        ymin = y(1)
        xmax = x(n)
        ymax = y(n)
      else
        increased = .False.
        xmin = x(n)
        ymin = y(n)
        xmax = x(1)
        ymax = y(1)
      end if

      ! Type of extrapolation
      select case(indx_extrapol)

        ! Zero out
        case(0)

          ! Check out of bounds
          do i=1,nn
            if (xx(i).lt.xmin) yy(i) = 0d0
            if (xx(i).gt.xmax) yy(i) = 0d0
          end do

        ! Extend extremes
        case(1)

          ! Check out bounds
          do i=1,nn
            if (xx(i).lt.xmin) yy(i) = ymin
            if (xx(i).gt.xmax) yy(i) = ymax
          end do

        ! Linear extrapolation
        case(2)

          ! If increasing
          if (increased) then

            ! Look for the minimum point in range
            do i=1,nn

              ! If found a point in range, select and exit
              if (xx(i).ge.x(1)) then
                mini = i
                exit
              end if

            end do ! All points

            ! Look for the maximum point in range
            do i=nn,1,-1

              ! If found a point in range, select and exit
              if (xx(i).le.x(n)) then
                maxi = i
                exit
              endif

            end do ! All points

            ! Linear interpolation between extremes
            a0 = (yy(mini+1) - yy(mini))/(xx(mini+1) - xx(mini))
            b0 = yy(mini) - a0*xx(mini)
            a1 = (yy(maxi) - yy(maxi-1))/(xx(maxi) - xx(maxi-1))
            b1 = yy(maxi) - a1*xx(maxi)

          ! If decreasing
          else

            ! Look for the minimum point in range
            do i=1,nn

              ! If found a point in range, select and exit
              if (xx(i).le.x(1)) then
                mini = i
                exit
              end if

            end do ! All points

            ! Look for the maximum point in range
            do i=nn,1,-1

              ! If found a point in range, select and exit
              if (xx(i).ge.x(n)) then
                maxi = i
                exit
              end if

            end do ! All points

            ! Get linear interpolation coefficients at both extremes
            a0 = (yy(mini) - yy(mini-1))/(xx(mini) - xx(mini-1))
            b0 = yy(mini) - a0*xx(mini)
            a1 = (yy(maxi+1) - yy(maxi))/(xx(maxi+1) - xx(maxi))
            b1 = yy(maxi) - a1*xx(maxi)

          end if ! Increasing or decreasing

          ! Exrapolate
          do i=1,nn
            if (xx(i).lt.xmin) yy(i) = a0*xx(i) + b0
            if (xx(i).gt.xmax) yy(i) = a1*xx(i) + b1
          end do

        ! No extrapolation
        case(3)

          return

        ! Error
        case default

          yy = 0d0

      end select ! Type of extrapolation

      return

      end subroutine Intpol

!####################################################################
!####################################################################
!####################################################################

      !> Perform cubic hermine interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !! Reference: Auer (2003)
      subroutine Intpol_Cubic_Hermite(x,y,n,xx,yy,nn)

      ! IO
      integer, intent(in):: n, nn
      double precision, dimension(n), intent(in):: x, y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local
      integer:: indx, i, j

      double precision:: D0, D1, H0, H1, YP0, YP1, U, UU, UUU


      H0 = x(2) - x(1)
      D0 = (y(2) - y(1))/H0
      YP0 = D0

      indx = 1
      do i=2,n

        if (i.lt.n) then

          H1 = x(i+1) - x(i)
          D1 = (y(i+1) - y(i))/H1
          YP1 = DERIV(H0, H1, D0, D1)

        else

          YP1 = (y(i) - y(i-1))/(x(i) - x(i-1))

        end if

        do j=indx,nn

          if ((xx(j) - x(i-1))*(xx(j) - x(i)).le.0d0) then

            U = (xx(j) - x(i-1))/H0
            UU = U*U
            UUU = UU*U
            yy(j) = (1d0 - 3d0*UU + 2d0*UUU)*y(i-1) + &
                    (3d0*UU - 2d0*UUU)*y(i) + &
                    (UUU - 2d0*UU + U)*H0*YP0 + (UUU - UU)*H0*YP1
            indx = indx + 1

          end if
        end do

        D0 = D1
        H0 = H1
        YP0 = YP1

      end do

      return

      end subroutine Intpol_Cubic_Hermite

!####################################################################
!####################################################################
!####################################################################

      !> Perform quadratic bezier interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !! Reference: Auer (2003)
      subroutine Intpol_Quadratic_Bezier(x, y, n, xx, yy, nn)

      ! IO
      integer, intent(in):: n, nn
      double precision, dimension(n), intent(in):: x, y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local
      integer:: indx, i, j

      double precision:: D0, D1, H0, H1, YP0, YP1, U, UU, C


      H0 = x(2) - x(1)
      D0 = (y(2) - y(1))/H0
      YP0 = D0

      indx = 1
      do i=2,n

        if (i.lt.n) then
          H1 = x(i+1) - x(i)
          D1 = (y(i+1) - y(i))/H1
          YP1 = DERIV(H0, H1, D0, D1)

        else

          YP1 = (y(i) - y(i-1))/(x(i) - x(i-1))

        end if

        !C = 0.5d0*(y(i-1) + y(i)) + H0*(YP0 - YP1)/6d0
        C = y(i-1)+H0*YP0/3d0

        do j=indx,nn

          if ((xx(j) - x(i-1))*(xx(j) - x(i)).le.0)then

            U = (xx(j) - x(i-1))/H0
            UU = U*U
            yy(j) = (1d0 - U)*(1d0 - U)*y(i-1) + &
                    UU*y(i) + C*2d0*U*(1d0 - U)
            indx = indx + 1

          end if
        end do

        D0 = D1
        H0 = H1
        YP0 = YP1

      end do

      return

      end subroutine Intpol_Quadratic_Bezier

!####################################################################
!####################################################################
!####################################################################

      !> Perform quadratic bezier interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !!  c0(double(:)): Interpolation coefficients\n
      !!  c1(double(:)): Interpolation coefficients\n
      !! Reference: Auer (2003)
      subroutine Intpol_Cubic_Bezier_C(x, y, n, xx, yy, nn, c0, c1)

      ! IO
      integer, intent(in):: n, nn
      double precision, dimension(n), intent(in):: x, y, c0, c1
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local
      integer:: indx, i, j

      double precision:: U, UU, UUU


      yy = 0
      indx = 1
      do i=1,n-1

        do j=indx,nn

          if ((xx(j) - X(i))*(xx(j) - X(i+1)).le.0d0)then

            U = (xx(j) - X(i))/(x(i) - x(i+1))
            UU = U*U
            UUU = UU*U
            yy(j) = (1d0 - U)*(1d0 - U)*(1d0 - U)*y(i) + &
                    UUU*y(i+1) + 3d0*U*(1d0 - U)*(1d0-U)*c0(i) + &
                    3d0*UU*(1d0 - U)*c1(i+1)
            indx = indx + 1

          end if
        end do
      end do

      return

      end subroutine Intpol_Cubic_Bezier_C

!####################################################################
!####################################################################
!####################################################################

      !> Perform quadratic bezier interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !!  c0(double(:)): Interpolation coefficients\n
      !!  c1(double(:)): Interpolation coefficients\n
      !! Reference: Auer (2003)
      subroutine Intpol_Quadratic_Bezier_C(x,y,n,xx,yy,nn,c0)

      ! IO
      integer, intent(in):: n, nn
      double precision, dimension(n), intent(in):: x, y, c0
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local
      integer:: indx, i, j

      double precision:: U, UU


      yy = 0
      indx = 1
      do i=1,n-1

        do j=indx,nn

          if ((xx(j) - X(i))*(xx(j) - X(i+1)).le.0d0) then

            U = (xx(j) - X(i))/(x(i) - x(i+1))
            UU = U*U
            yy(j) = (1d0 - U)*(1d0 - U)*y(i) + &
                    UU*y(i+1) + c0(i)*2d0*U*(1d0 - U)
            indx = indx + 1

          end if
        end do
      end do

      return

      end subroutine Intpol_Quadratic_Bezier_C

!####################################################################
!####################################################################
!####################################################################

      !> Perform linear interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis
      subroutine Intpol_Lin(x,y,n,xx,yy,nn)

      ! IO
      integer, intent(in):: n, nn
      double precision, dimension(n), intent(in):: x, y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local
      integer:: i, j, indx

      double precision:: a, b

      ! Single point
      if (n.eq.1) then

        yy = y(1)

      ! More than one point
      else if(n.gt.1) then

        indx = 1
        do i=1,n-1

          a = (y(i+1) - y(i))/(x(i+1) - x(i))
          b = y(i) - a*x(i)

          do j=indx,nn

            if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0) then

              yy(j) = a*xx(j) + b
              indx = indx + 1

            end if
          end do
        end do

      ! Error
      else

        yy = 0d0

      end if

      return

      end subroutine Intpol_Lin

!####################################################################
!####################################################################
!####################################################################

      !> Perform linear interpolation in Stokes parameters\n
      !!   x(double(:)): Input x axis\n
      !! y(double(:,:)): Input y axis\n
      !!     m(integer): Dimensionality of yy arrays\n
      !!     x(integer): Dimensionality of y arrays\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !! yy(double(:,:)): Output y axis
      subroutine Intpol_Lin_stk(x,y,m,k,n,xx,yy,nn)

      ! IO
      integer, intent(in):: m, n, nn, k
      double precision, dimension(n), intent(in):: x
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(k,n), intent(in):: y
      double precision, dimension(m,nn), intent(out):: yy

      ! Local
      integer:: i, j, indx

      double precision, dimension(k):: a, b


      ! Single point
      if (n.eq.1) then

        do i=1,n
          yy(1:k,i) = y(:,1)
          yy(k+1:m,i) = 0d0
        end do

      ! More than one point
      else if(n.gt.1) then

        indx = 1
        do i=1,n-1

          a = (y(:,i+1) - y(:,i))/(x(i+1) - x(i))
          b = y(:,i) - a*x(i)

          do j=indx,nn

            if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0) then

              yy(1:k,j) = a*xx(j) + b
              yy(k+1:m,j) = 0d0
              indx = indx + 1

            end if
          end do
        end do

      ! Error
      else

        yy = 0d0

      end if

      return

      end subroutine Intpol_Lin_stk

!####################################################################
!####################################################################
!####################################################################

      !> Perform linear interpolation in several variables\n
      !!   x(double(:)): Input x axis\n
      !! y(double(:,:)): Input y axis\n
      !!     m(integer): Dimensionality of y arrays\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !! yy(double(:,:)): Output y axis
      subroutine Intpol_Lin_mD(x,y,m,n,xx,yy,nn)

      ! IO
      integer, intent(in):: m, n, nn
      double precision, dimension(n), intent(in):: x
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(m,n), intent(in):: y
      double precision, dimension(m,nn), intent(out):: yy

      ! Local
      integer:: i, j, indx

      double precision, dimension(m):: a, b


      ! Single point
      if (n.eq.1) then

        do i=1,n
          yy(:,i) = y(:,1)
        end do

      ! More than one point
      else if(n.gt.1) then

        indx = 1
        do i=1,n-1

          a = (y(:,i+1) - y(:,i))/(x(i+1) - x(i))
          b = y(:,i) - a*x(i)

          do j=indx,nn

            if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0) then

              yy(:,j) = a*xx(j) + b
              indx = indx + 1

            end if
          end do
        end do

      ! Error
      else

        yy = 0d0

      end if

      return

      end subroutine Intpol_Lin_mD

!####################################################################
!####################################################################
!####################################################################

      !> Get value of X at the minimum of Y in parabolic
      !! parabolic interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!  Indx(integer): Index minimum Y\n
      !!  X_min(double): x value for minimum Y
      Subroutine Parabolic(X,Y,Indx,X_min)

      ! IO
      integer, intent(in):: Indx
      double precision,intent(out):: X_min
      double precision, dimension(:), intent(in):: X, Y

      ! Local
      integer:: Length

      double precision:: X_tmp1, X_tmp2, Y_tmp1, Y_tmp2
      double precision, dimension(:), allocatable:: XP


      Length = size(X)
      allocate(XP(Length))

      XP(1:Indx+1) = log(X(1:Indx+1))

      X_tmp1 = XP(Indx) - XP(Indx-1)
      X_tmp2 = XP(Indx) - XP(Indx+1)
      Y_tmp1 = Y(Indx) - Y(Indx-1)
      Y_tmp2 = Y(Indx) - Y(Indx+1)

     !a = (Y_tmp2/X_tmp2-Y_tmp1/X_tmp1)/(XP(Indx+1)-XP(Indx-1))
     !b = Y_tmp1/X_tmp1-a*(XP(Indx)+XP(Indx-1))
     !c = Y(Indx)-a*XP(Indx)**2-b*XP(Indx)

      X_min = XP(Indx) - 0.5d0*(X_tmp1*X_tmp1*Y_tmp2 - &
                                X_tmp2*X_tmp2*Y_tmp1)/ &
                         (X_tmp1*Y_tmp2 - X_tmp2*Y_tmp1)

      X_min = exp(X_min)

      deallocate(XP)

      return

      end subroutine Parabolic

!####################################################################
!####################################################################
!####################################################################

      !> Get derivative for interpolations in inversion
      !!  H0(double): Back step\n
      !!  H1(double): Forward step\n
      !!  D0(double): Back slope\n
      !!  D1(double): Forward slope\n
      !! YP1(double): Output value\n
      !! Reference: Auer (2003)
      function DERIV(H0,H1,D0,D1) result(YP1)

      ! IO
      double precision:: H0, H1, D0, D1
      double precision:: YP1

      ! Local
      double precision:: ALPHA


      if (D0*D1.gt.0) then

        ALPHA = (1d0 + H1/(H0+H1))/3d0
        YP1 = (D0*D1)/(ALPHA*D1 + (1d0 - ALPHA)*D0)
      else

        YP1 = 0d0

      end if

      end function DERIV

!####################################################################
!####################################################################
!####################################################################

      end module inter_mod
