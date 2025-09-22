      !> Interpolation routines
      module inter_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!  Start:
!     18/04/2017
!  Last version:
!     29/08/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     29/08/2025:    V4.0.1 - The way of extrapolating in Intpol can
!                             be different for top and bottom (TdPA)
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
!  colinter
!    Interpolate inelastic collisions using cubic splines and, if
!  resulting in any negative rate, changing to linear interpolation
!
!  bilinear
!    Two-dimensional bilinear interpolation. Out of boundary values
!  are taken as extended from the tabulation data as constant. Admits
!  a single set of output coordinates
!
!  linear
!    One-dimensional linear interpolation. Out of boundary values are
!  taken as extended from the tabulation data as constant. Admits a
!  single output coordinate
!
!  spline
!    Calculate the coefficients for one-dimensional cubic spline
!  interpolation
!
!  ispline
!    Evaluate the cubic spline interpolation given the coefficients.
!  Out of boundary values are taken as extended from the tabulation
!  data as constant. Admits a single output coordinate
!  Alex G: January 2010
!
!  spline_2d
!    Two-dimensional cubic spline interpolation. Admits a single set
!  of output coordinates 
!
!  Intpol
!    General interpolation routine for the model atmospheres in the
!  inversion mode
!
!  Intpol_Cubic_Hermite
!    One-dimensional cubic hermite interpolation. Does not have
!  control over out of bounds coordinates
!
!  Intpol_Quadratic_Bezier
!    One-dimensional quadratic bezier interpolation. Does not have
!  control over out of bounds coordinates
!
!  Intpol_Cubic_Bezier_C
!    One-dimensional cubic bezier interpolation with given
!  coefficients. Does not have control over out of bounds coordinates
!
!  Intpol_Quadratic_Bezier_C
!    One-dimensional quadratic bezier interpolation with given
!  coefficients. Does not have control over out of bounds coordinates
!
!  Intpol_Lin
!    One-dimensional linear interpolation. Does not have control over
!  out of bounds coordinates
!
!  Intpol_Lin_stk
!    Linear interpolation for Stokes parameters in wavelength. Does
!  not have control over out of bounds coordinates
!
!  Intpol_Lin_mD
!    Linear interpolation for multiple data arrays. Does not have
!  control over out of bounds coordinates
!
!  Parabolic
!    Calculate the coordinate value which gives the minimum value
!  for a parabolic interpolation
!
!  DERIV
!    Calculate derivatives for several interpolation functions in the
!  module
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

      !> Interpolate inelastic collisions using cubic splines and, if
      !! resulting in any negative rate, changing to linear
      !! interpolation\n
      !!    xin(double(:)): Input temperature axis\n
      !!    yin(double(:)): Input collisional rate tabulation\n
      !!      nin(integer): Size of xin and yin\n
      !!   xout(double(:)): Output temperature axis\n
      !!   yout(double(:)): Interpolated rates\n
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
      double precision, dimension(:), intent(in):: xin, yin
      double precision, dimension(:), intent(out):: xout, yout

      ! Local

      integer:: iout

      double precision, dimension(nin):: b, c, d


      ! Checks if we were forced to use linear when splines
      ! were requested
      llinear = .False.

      ! If only one input value
      if (nin.eq.1) then

        ! Output is the same constant
        yout = yin(1)

      ! If there are only two inputs, it has to be linear
      else if (nin.eq.2) then

        ! For each output temperature
        do iout=1,nout

          ! Do linear interpolation
          call linear(xin,yin,xout(iout),yout(iout))

        end do ! Output temperatures

      ! If forcing linear interpolation
      else if (flinear) then

        ! For each temperature
        do iout=1,nout

          ! Do linear interpolation
          call linear(xin,yin,xout(iout),yout(iout))

        end do ! Output temperatures

      ! Usual path
      else

        ! Prepare spline interpolation
        call spline(xin,yin,b,c,d,nin)

        ! For each output temperature
        do iout=1,nout

          ! Spline interpolation
          yout(iout) = ispline(xout(iout),xin,yin,b,c,d,nin)

          ! If negative result
          if (yout(iout).lt.0) then

            ! Flag linear interpolation and abort splines
            llinear = .True.
            exit

          end if ! Negative result

        end do ! Output temperatures

        ! The interpolation has to be linear afterall
        if (llinear) then

          ! For each output temperature
          do iout=1,nout

            ! Perform linear interpolation
            call linear(xin,yin,xout(iout),yout(iout))

          end do ! Output temperatures

        end if ! We need linear interpolation

      end if ! Dimensionality and forcing linear interpolation

      end subroutine colinter

!#####################################################################
!#####################################################################
!#####################################################################

      !> Two-dimensional bilinear interpolation. Out of boundary
      !! values are taken as extended from the tabulation data as
      !! constant. Admits a single set of output coordinates\n
      !!     x1(double(:)): Input coordinates axis 1\n
      !!     x2(double(:)): Input coordinates axis 2\n
      !!    y(double(:,:)): Data to interpolate\n
      !!        z1(double): Output coordinates axis 1\n
      !!        z2(double): Output coordinates axis 2\n
      !!      yout(double): Interpolated data
      subroutine bilinear(x1,x2,y,z1,z2,yout)

      ! I/O

      double precision, intent(in):: z1,z2
      double precision, intent(out):: yout
      double precision, dimension(:), intent(in):: x1, x2
      double precision, dimension(:,:), intent(in):: y

      ! Local

      integer:: n1,n2,ii,i10,i11,i20,i21

      double precision:: x10,x11,x20,x21,yA,yB,dx1,dx2


      ! Get sizes pf tabiñatopm
      n1 = size(x1)
      n2 = size(x2)

      ! If within bounds
      if (z1.ge.x1(1).and.z1.le.x1(n1).and. &
          z2.ge.x2(1).and.z2.le.x2(n2)) then

        ! For each index in axis 1
        do ii=1,n1-1

          ! If output in current range
          if (z1.ge.x1(ii).and.z1.le.x1(ii+1)) then

            ! Identify position and factors for linear interpolation
            i10 = ii
            i11 = ii+1
            x10 = x1(i10)
            x11 = x1(i11)
            dx1 = 1d0/(x11 - x10)

            ! And finish search
            exit

          end if ! Output in current range

        end do ! Index in axis 1

        ! For each index in axis 2
        do ii=1,n2-1

          ! If output in current range
          if (z2.ge.x2(ii).and.z2.le.x2(ii+1)) then

            ! Identify position and factors for linear interpolation
            i20 = ii
            i21 = ii+1
            x20 = x2(i20)
            x21 = x2(i21)
            dx2 = 1d0/(x21 - x20)

            ! And finish search
            exit

          end if ! Output in current range

        end do ! Index in axis 2

        ! Interpolate in axis 2 at left axis 1 position
        yA = ((y(i10,i21) - y(i10,i20))*z2 + &
              y(i10,i20)*x21 - y(i10,i21)*x20)*dx2

        ! Interpolate in axis 2 at right axis 1 position
        yB = ((y(i11,i21) - y(i11,i20))*z2 + &
              y(i11,i20)*x21 - y(i11,i21)*x20)*dx2

        ! Interpolate in axis 1
        yout = ((yB - yA)*z1 + yA*x11 - yB*x10)*dx1

      ! Out of x2 bounds but inside x1 bounds
      else if (z1.ge.x1(1).and.z1.le.x1(n1)) then

        ! For each index in axis 1
        do ii=1,n1-1

          ! If output in current range
          if (z1.ge.x1(ii).and.z1.le.x1(ii+1)) then

            ! Identify position and factors for linear interpolation
            i10 = ii
            i11 = ii+1
            x10 = x1(i10)
            x11 = x1(i11)
            dx1 = 1d0/(x11 - x10)

            ! And finish search
            exit

          end if ! Output in current range

        end do ! Index in axis 1

        ! If out of axis 2 at the bottom
        if (z2.lt.x2(1)) then

          ! Get first number
          i20 = 1

        ! If out of axis 2 at the top
        else if (z2.gt.x2(n2)) then

          ! Get last number
          i20 = n2

        end if ! Closest axis 2 boundary

        ! Interpolate in axis 1 along one axis 2 boundary
        yout = ((y(i11,i20) - y(i10,i20))*z1 + &
                y(i10,i20)*x11 - y(i11,i20)*x10)*dx1

      ! Out of x1 bounds but inside x2 bounds
      else if (z2.ge.x2(1).and.z2.le.x2(n2)) then

        ! For each index in axis 2
        do ii=1,n2-1

          ! If output in current range
          if (z2.ge.x2(ii).and.z2.le.x2(ii+1)) then

            ! Identify position and factors for linear interpolation
            i20 = ii
            i21 = ii+1
            x20 = x2(i20)
            x21 = x2(i21)
            dx2 = 1d0/(x21 - x20)

            ! And finish search
            exit

          end if ! Output in current range

        end do ! Index in axis 2

        ! If out of axis 1 at the bottom
        if (z1.lt.x1(1)) then

          ! Get first number
          i10 = 1

        ! If out of axis 2 at the top
        else if (z1.gt.x1(n1)) then

          ! Get last number
          i10 = n1

        end if ! Closest axis 1 boundary

        ! Interpolate in axis 2 along one axis 1 boundary
        yout = ((y(i10,i21) - y(i10,i20))*z2 + &
                y(i10,i20)*x21 - y(i10,i21)*x20)*dx2

      ! Completely out of bounds
      else

        ! If out of axis 1 at the bottom
        if (z1.lt.x1(1)) then

          ! Get first number
          i10 = 1

        ! If out of axis 1 at the top
        else if (z1.gt.x1(n1)) then

          ! Get last number
          i10 = n1

        end if ! Closest axis 1 boundary

        ! If out of axis 2 at the bottom
        if (z2.lt.x2(1)) then

          ! Get first number
          i20 = 1

        ! If out of axis 2 at the top
        else if (z2.gt.x2(n2)) then

          ! Get last number
          i20 = n2

        end if ! Closest axis 2 boundary

        ! Get closest corner
        yout = y(i10,i20)

      end if ! Inside or outside which boundaries in the tabulation

      end subroutine bilinear

!#####################################################################
!#####################################################################
!#####################################################################

      !> Linear interpolator\n
      !> One-dimensional linear interpolation. Out of boundary values
      !! are taken as extended from the tabulation data as constant.
      !! Admits a single output coordinate\n
      !!    x(double(:)): Input x axis\n
      !!    y(double(:)): Input y axis\n
      !!       z(double): Output x value\n
      !!    yout(double): Interpolated value
      subroutine linear(x,y,z,yout)

      ! I/O

      double precision, intent(in):: z
      double precision, intent(out):: yout
      double precision, dimension(:), intent(in):: x,y

      ! Local

      integer:: n,ii,i0,i1

      double precision:: y0,y1,x0,x1,dx


      ! Size of tabulation
      n = size(x)

      ! Out of lower bound
      if (z.le.x(1)) then

        ! Take leftmost value
        yout = y(1)

      ! Out of upper bound
      else if (z.ge.x(n)) then

        ! Take rightmost value
        yout = y(n)

      ! Between bounds
      else

        ! For every input coordinate
        do ii=1,n-1

          ! Check if output in current range
          if (z.ge.x(ii).and.z.le.x(ii+1)) then

            ! Take indexes and data for interpolation
            i0 = ii
            i1 = ii+1
            x0 = x(i0)
            x1 = x(i1)
            y0 = y(i0)
            y1 = y(i1)
            dx = 1d0/(x1 - x0)

            ! Finish search
            exit

          end if ! If within current range

        end do ! Input coordinates

        ! Linear interpolation
        yout = ((y1 - y0)*z + y0*x1 - y1*x0)*dx

      end if ! If within or out of which boundary

      end subroutine linear

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the coefficients for one-dimensional cubic spline
      !! interpolation\n
      !!    x(double(:)): Input x axis\n
      !!    y(double(:)): Input y axis\n
      !!    b(double(:)): Coefficient of spline interpolation\n
      !!    c(double(:)): Coefficient of spline interpolation\n
      !!    d(double(:)): Coefficient of spline interpolation\n
      !!      n(integer): Size of x,y,b,c,d
      subroutine spline(x,y,b,c,d,n)

      ! I/O

      integer, intent(in):: n
      double precision, dimension(:), intent(in):: x,y
      double precision, dimension(:), intent(out):: b,c,d

      ! Local

      integer:: i,j,k

      double precision:: a


      ! Get dimension minus 1
      k = n - 1

      ! If we do not have at least two points, abort
      if (n.lt.2) return

      ! If exactly two points
      if(n.lt.3)then

        ! Prepare linear interpolation
        b(1) = (y(2)-y(1))/(x(2)-x(1))
        c(1) = 0d0
        d(1) = 0d0
        b(2) = b(1)
        c(2) = 0d0
        d(2) = 0d0

        ! And done
        return

      end if ! Exactly two points

      !
      ! Prepare spline coefficients
      !

      ! Initial step and slope
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)

      ! For each intermediate point
      do i=2,k

        ! Get step
        d(i) = x(i+1) - x(i)

        ! Combine current and last step
        b(i) = 2d0*(d(i-1) + d(i))

        ! Next slope
        c(i+1) = (y(i+1) - y(i))/d(i)

        ! Slope difference
        c(i) = c(i+1) - c(i)

      end do ! Intermediate points

      ! Reverse first and last b
      b(1) = -d(1)
      b(n) = -d(n-1)

      ! Reset extremal slopes
      c(1) = 0d0
      c(n) = 0d0

      ! If more than three points
      if(n.ne.3)then

        ! Build derivatives at the boundaries
        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)*d(1)/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))

      end if ! More than three points

      ! From the second to last point
      do i=2,n

        ! Forward elimination
        a = d(i-1)/b(i-1)
        b(i) = b(i) - a*d(i-1)
        c(i) = c(i) - a*c(i-1)

      end do ! From second to last

      ! Scale last point
      c(n) = c(n)/b(n)

      ! From first to second-to-last point
      do j=1,k

        ! Back substitution
        i = n-j
        c(i) = (c(i) - d(i)*c(i+1))/b(i)

      end do ! From first to second-to-last

      ! Finalize b coefficient at last boundary
      b(n) = (y(n) - y(k))/d(k) + d(k)*(c(k) + 2d0*c(n))

      ! From first to the second-to-last point
      do i=1,k

        ! Finaline coefficients
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3d0*c(i)

      end do ! From first to the second-to-last point

      ! Finalize c and d at last boundary
      c(n) = 3d0*c(n)
      d(n) = d(n-1)

      end subroutine spline

!#####################################################################
!#####################################################################
!#####################################################################

      !> Evaluate the cubic spline interpolation given the
      !! coefficients. Out of boundary values are taken as extended
      !! from the tabulation data as constant. Admits a single output
      !! coordinate\n
      !!       u(double): Output x value\n
      !!    x(double(:)): Input x axis\n
      !!    y(double(:)): Input y axis\n
      !!    b(double(:)): Coefficient of spline interpolation\n
      !!    c(double(:)): Coefficient of spline interpolation\n
      !!    d(double(:)): Coefficient of spline interpolation\n
      !!      n(integer): Size of x,y,b,c,d
      double precision function ispline(u,x,y,b,c,d,n)

      ! I/O

      integer, intent(in):: n
      double precision, intent(in):: u
      double precision, dimension(:), intent(in):: x, y, b, c, d

      ! Local

      integer:: i,j,k

      double precision:: dx


      ! If below bottom boundary
      if(u.le.x(1))then

        ! Extend constant
        ispline = y(1)

        ! And return
        return

      end if ! Below bottom boundary

      ! If above upper boundary
      if(u.ge.x(n))then

        ! Extend constant
        ispline = y(n)

        ! And return
        return

      end if ! Above upper boundary

      !
      ! Binary search for for i, such that x(i) <= u <= x(i+1)
      !

      ! Initialize
      i = 1
      j = n + 1

      ! While not found
      do while (j.gt.i+1)

        ! Get middle index
        k = (i+j)/2

        ! If output coordinate below this one
        if (u.lt.x(k)) then

          ! New j value
          j=k

        ! If output coordinate above this one
        else

          ! New i value
          i=k

        end if ! Middle index coordinate value above or below output

      end do ! While not found

      ! Distance to previous coordinate
      dx = u - x(i)

      ! Evaluate spline interpolation
      ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

      return

      end function ispline

!#####################################################################
!#####################################################################
!#####################################################################

      !> Two-dimensional cubic spline interpolation. Admits a single
      !! set of output coordinates 
      !!     x1(double(:)): Input coordinates axis 1\n
      !!     x2(double(:)): Input coordinates axis 2\n
      !!    y(double(:,:)): Data to interpolate\n
      !!       n1(integer): Size of x1 and second column y\n
      !!       n2(integet): Size of x2 and first column y\n
      !!        z1(double): Output coordinates axis 1\n
      !!        z2(double): Output coordinates axis 2\n
      !!      yout(double): Interpolated data
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


      ! For each position in x1
      do i1=1,n1

        ! Take a slice of y in the first dimension
        yslice2 = y(:,i1)

        !
        ! Calculate coefficients in dimension 2 of the slice
        !

        ! Initialize
        y2slice2(1) = 0d0
        y2slice2(n2) = 0d0
        V2(1) = 0d0

        ! For intermediate positions in second axis
        do i2=2,n2-1

          ! Get derivative
          d1 = (x2(i2) - x2(i2-1))/(x2(i2+1) - x2(i2-1))
          d2 = d1*y2slice2(i2-1) + 2d0
          y2slice2(i2) = (d1 - 1d0)/d2
          V2(i2) = (6d0*((yslice2(i2+1) - yslice2(i2))/ &
                         (x2(i2+1) - x2(i2)) - &
                         (yslice2(i2) - yslice2(i2-1))/ &
                         (x2(i2) - x2(i2-1)))/ &
                        (x2(i2+1) - x2(i2-1)) - d1*V2(i2-1))/d2

        end do ! Intermediate positions

        ! From last to first point
        do i2=n2-1,1,-1

          ! Complete calculation of coefficient
          y2slice2(i2) = y2slice2(i2)*y2slice2(i2+1) + V2(i2)

        end do ! From last to first point

        ! Initialize search for range of output
        iini = 1
        ifin = n2

        ! While not found
        do while ((ifin-iini).gt.1)

          ! Get current middle point
          imed = (iini + ifin)/2

          ! If output is below
          if (x2(imed).gt.z2) then

            ! Update right index
            ifin = imed

          ! If output is above
          else

            ! Update left index
            iini = imed

          end if ! Above or below output

        end do ! While searching

        ! Get coefficients
        dx2 = x2(ifin) - x2(iini)
        d1 = (x2(ifin) - z2)/dx2
        d2 = (z2 - x2(iini))/dx2

        ! Interpolate
        yslice1(i1) = d1*yslice2(iini) + d2*yslice2(ifin) + &
                      ((d1*(d1*d1 - 1d0))*y2slice2(iini) + &
                       (d2*(d2*d2 - 1d0))*y2slice2(ifin))* &
                      dx2*dx2/6d0

      end do ! For each position in axis 1

      !
      ! Calculate coefficients in direction 1
      !

      ! Initialize
      y2slice1(1) = 0d0
      y2slice1(n1) = 0d0
      V1(1) = 0d0

      ! For intermediate positions in first axis
      do i1=2,n1-1

        ! Get derivative
        d1 = (x1(i1) - x1(i1-1))/(x1(i1+1) - x1(i1-1))
        d2 = d1*y2slice1(i1-1) + 2d0
        y2slice1(i1) = (d1 - 1d0)/d2
        V1(i1) = (6d0*((yslice1(i1+1) - yslice1(i1))/ &
                       (x1(i1+1) - x1(i1)) - &
                       (yslice1(i1) - yslice1(i1-1))/ &
                       (x1(i1) - x1(i1-1)))/ &
                      (x1(i1+1) - x1(i1-1)) -d1*V1(i1-1))/d2

      end do ! Intermediate positions

      ! From last to first point
      do i1=n1-1,1,-1

        ! Complete calculation of coefficient
        y2slice1(i1) = y2slice1(i1)*y2slice1(i1+1) + V1(i1)

      end do ! From last to first point

      ! Initialize search for range of output
      iini = 1
      ifin = n1

      ! While not found
      do while ((ifin-iini).gt.1)

        ! Get current middle point
        imed = (iini + ifin)/2

        ! If output is below
        if (x1(imed).gt.z1) then

          ! Update right index
          ifin = imed

        ! If output is above
        else

          ! Update left index
          iini = imed

        end if ! Above or below output

      end do ! While searching

      ! Get coefficients
      dx1 = x1(ifin) - x1(iini)
      d1 = (x1(ifin) - z1)/dx1
      d2 = (z1 - x1(iini))/dx1

      ! Interpolate
      yout = d1*yslice1(iini) + d2*yslice1(ifin) + &
             ((d1*(d1*d1 - 1d0))*y2slice1(iini) + &
              (d2*(d2*d2 - 1d0))*y2slice1(ifin))* &
             (dx1*dx1/6d0)

      end subroutine spline_2d

!#####################################################################
!#####################################################################
!#####################################################################

      !> General interpolation routine for the model atmospheres in
      !! the inversion mode\n
      !!            x(double(:)): Input x axis\n
      !!            y(double(:)): Input y axis\n
      !!              n(integer): Size of x and y arrays\n
      !!           xx(double(:)): Output x axis\n
      !!           yy(double(:)): Output y axis\n
      !!             nn(integer): Size of xx and yy arrays\n
      !!    Indx_Intpol(integer): Type of interpolation:\n
      !!                           0 linear\n
      !!                           1 quadratic bezier\n
      !!                           2 cubic bezier\n
      !!  Indx_Extrapol(integer): Type of extrapolation
      subroutine Intpol(x,y,n,xx,yy,nn,Indx_Intpol,Indx_Extrapol)

      ! I/O

      integer, intent(in):: n,nn,Indx_Intpol,Indx_Extrapol
      double precision, dimension(n), intent(in):: x,y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local

      logical:: Increased

      integer:: i,mini,maxi,extrap_top,extrap_bot

      double precision:: aa,bb,xmin,xmax,ymin,ymax


      ! Choose the interpolation type
      select case(Indx_Intpol)

        ! Linear
        case(0)

          ! Perform linear interpolation
          call Intpol_Lin(x,y,n,xx,yy,nn)

        ! Quadratic bezier
        case(1)

          ! Only if enough elements
          if(n.gt.2) then

            ! Perform quadratic Bezier interpolation
            call Intpol_Quadratic_Bezier(x,y,n,xx,yy,nn)

          ! Otherwise
          else

            ! Perform linear interpolation
            call Intpol_Lin(x,y,n,xx,yy,nn)

          endif ! Enough elements for quadratic bezier

        ! Cubic hermite
        case(2)

          ! Only if enough elements
          if(n.gt.2) then

            ! Perform cubic hermite interpolation
            call Intpol_Cubic_Hermite(x,y,n,xx,yy,nn)

          ! Otherwise
          else

            ! Perform linear interpolation
            call Intpol_Lin(x,y,n,xx,yy,nn)

          endif ! Enough elements for cubic hermite

        ! Unknown interpolation type
        case default

          ! Return 0
          yy = 0d0

      end select ! Type of interpolation

      ! If input dimension was 1, return already
      if (n.eq.1) return

      ! If not extrapolating, return already
      if (Indx_Extrapol.eq.0) return

      !
      ! Extrapolation
      !

      ! Get type for each boundary
      extrap_top = indx_extrapol/4
      extrap_bot = mod(indx_extrapol,4)

      ! Initialize
      mini = 0d0
      maxi = 0d0

      ! Determine ordering of input axis
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

      !
      ! Top boundary
      !

      ! If extrapolation
      if (extrap_top.gt.0.and.extrap_top.lt.4) then

        ! Type of extrapolation
        select case(extrap_top)

          ! Zero out
          case(1)

            ! Check out of bounds
            do i=1,nn
              if (xx(i).lt.xmin) yy(i) = 0d0
            end do

          ! Extend extremes as constant
          case(2)

            ! Check out bounds
            do i=1,nn
              if (xx(i).lt.xmin) yy(i) = ymin
            end do

          ! Linear extrapolation
          case(3)

            ! If increasing
            if (increased) then

              ! Look for the minimum point in range
              do i=1,nn

                ! If found a point in range
                if (xx(i).ge.x(1)) then

                  ! Select and exit
                  mini = i
                  exit

                end if

              end do ! All points

              ! Coefficients for linear interpolation between extremes
              aa = (yy(mini+1) - yy(mini))/(xx(mini+1) - xx(mini))
              bb = yy(mini) - aa*xx(mini)

            ! If decreasing
            else

              ! Look for the minimum point in range
              do i=1,nn

                ! If found a point in range
                if (xx(i).le.x(1)) then

                  ! Select and exit
                  mini = i
                  exit

                end if

              end do ! All points

              ! Get linear interpolation coefficients at both extremes
              aa = (yy(mini) - yy(mini-1))/(xx(mini) - xx(mini-1))
              bb = yy(mini) - aa*xx(mini)

            end if ! Increasing or decreasing

            ! Extrapolate
            do i=1,nn
              if (xx(i).lt.xmin) yy(i) = aa*xx(i) + bb
            end do

        end select ! Type of extrapolation

      end if ! Extrapolation top

      !
      ! Bottom boundary
      !

      ! If extrapolation
      if (extrap_bot.gt.0.and.extrap_bot.lt.4) then

      ! Type of extrapolation
      select case(extrap_bot)

        ! Zero out
        case(1)

          ! Check out of bounds
          do i=1,nn
            if (xx(i).gt.xmax) yy(i) = 0d0
          end do

        ! Extend extremes as constant
        case(2)

          ! Check out bounds
          do i=1,nn
            if (xx(i).gt.xmax) yy(i) = ymax
          end do

        ! Linear extrapolation
        case(3)

          ! If increasing
          if (increased) then

            ! Look for the maximum point in range
            do i=nn,1,-1

              ! If found a point in range
              if (xx(i).le.x(n)) then

                ! Select and exit
                maxi = i
                exit

              endif

            end do ! All points

            ! Coefficients for linear interpolation between extremes
            aa = (yy(maxi) - yy(maxi-1))/(xx(maxi) - xx(maxi-1))
            bb = yy(maxi) - aa*xx(maxi)

          ! If decreasing
          else

            ! Look for the maximum point in range
            do i=nn,1,-1

              ! If found a point in range
              if (xx(i).ge.x(n)) then

                ! Select and exit
                maxi = i
                exit

              end if

            end do ! All points

            ! Get linear interpolation coefficients at both extremes
            aa = (yy(maxi+1) - yy(maxi))/(xx(maxi+1) - xx(maxi))
            bb = yy(maxi) - aa*xx(maxi)

          end if ! Increasing or decreasing

          ! Extrapolate
          do i=1,nn
            if (xx(i).gt.xmax) yy(i) = aa*xx(i) + bb
          end do

        end select ! Type of extrapolation

      end if ! Extrapolating bottom

      return

      end subroutine Intpol

!#####################################################################
!#####################################################################
!#####################################################################

      !> One-dimensional cubic hermite interpolation. Does not have
      !! control over out of bounds coordinates\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !! Reference: Auer (2003)
      subroutine Intpol_Cubic_Hermite(x,y,n,xx,yy,nn)

      ! I/O

      integer, intent(in):: n,nn
      double precision, dimension(n), intent(in):: x,y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local

      integer:: indx,i,j

      double precision:: D0,D1,H0,H1,YP0,YP1,U,UU,UUU


      ! Step and slope
      H0 = x(2) - x(1)
      D0 = (y(2) - y(1))/H0
      YP0 = D0

      ! Initialize search index
      indx = 1

      ! From the second to last point
      do i=2,n

        ! If not the last
        if (i.lt.n) then

          ! Step and slope
          H1 = x(i+1) - x(i)
          D1 = (y(i+1) - y(i))/H1

          ! Get derivative
          YP1 = DERIV(H0, H1, D0, D1)

        ! Last point
        else

          ! First order derivative
          YP1 = (y(i) - y(i-1))/(x(i) - x(i-1))

        end if ! Last point or not

        ! From indx to the last output
        do j=indx,nn

          ! If output within current range
          if ((xx(j) - x(i-1))*(xx(j) - x(i)).le.0d0) then

            ! Get scaled distance and its powers
            U = (xx(j) - x(i-1))/H0
            UU = U*U
            UUU = UU*U

            ! Cubic hermite interpolation
            yy(j) = (1d0 - 3d0*UU + 2d0*UUU)*y(i-1) + &
                    (3d0*UU - 2d0*UUU)*y(i) + &
                    (UUU - 2d0*UU + U)*H0*YP0 + (UUU - UU)*H0*YP1

            ! Advance search index
            indx = j + 1

          end if ! Output within current range

        end do ! From indx to the last output

        ! Update left values
        D0 = D1
        H0 = H1
        YP0 = YP1

      end do ! From second to last point

      return

      end subroutine Intpol_Cubic_Hermite

!#####################################################################
!#####################################################################
!#####################################################################

      !> One-dimensional quadratic bezier interpolation. Does not have
      !! control over out of bounds coordinates\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !! Reference: Auer (2003)
      subroutine Intpol_Quadratic_Bezier(x,y,n,xx,yy,nn)

      ! I/O

      integer, intent(in):: n,nn
      double precision, dimension(n), intent(in):: x,y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local

      integer:: indx,i,j

      double precision:: D0,D1,H0,H1,YP0,YP1,U,UU,C


      ! Step and slope
      H0 = x(2) - x(1)
      D0 = (y(2) - y(1))/H0
      YP0 = D0

      ! Initialize search index
      indx = 1

      ! From second to last input points
      do i=2,n

        ! If not the last
        if (i.lt.n) then

          ! Get step and derivative
          H1 = x(i+1) - x(i)
          D1 = (y(i+1) - y(i))/H1
          YP1 = DERIV(H0, H1, D0, D1)

        ! Last point
        else

          ! First order derivative
          YP1 = (y(i) - y(i-1))/(x(i) - x(i-1))

        end if

       !C = 0.5d0*(y(i-1) + y(i)) + H0*(YP0 - YP1)/6d0
        C = y(i-1)+H0*YP0/3d0

        ! From indx to the last output
        do j=indx,nn

          ! If output within current range
          if ((xx(j) - x(i-1))*(xx(j) - x(i)).le.0)then

            ! Get scaled distance and its square
            U = (xx(j) - x(i-1))/H0
            UU = U*U

            ! Quadratic bezier interpolation
            yy(j) = (1d0 - U)*(1d0 - U)*y(i-1) + &
                    UU*y(i) + C*2d0*U*(1d0 - U)

            ! Advance search index
            indx = j + 1

          end if ! Output within current range

        end do ! From indx to the last output

        ! Update left values
        D0 = D1
        H0 = H1
        YP0 = YP1

      end do ! From second to last point

      return

      end subroutine Intpol_Quadratic_Bezier

!#####################################################################
!#####################################################################
!#####################################################################

      !> One-dimensional cubic bezier interpolation with given
      !! coefficients. Does not have control over out of bounds
      !! coordinates\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !!  c0(double(:)): Interpolation coefficient\n
      !!  c1(double(:)): Interpolation coefficient\n
      !! Reference: Auer (2003)
      subroutine Intpol_Cubic_Bezier_C(x,y,n,xx,yy,nn,c0,c1)

      ! I/O

      integer, intent(in):: n, nn
      double precision, dimension(n), intent(in):: x,y,c0,c1
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local

      integer:: indx,i,j

      double precision:: U,UU,UUU


      ! Initialize output
      yy = 0d0

      ! Initialize search index
      indx = 1

      ! From first to second-to-last input points
      do i=1,n-1

        ! From indx to the last output
        do j=indx,nn

          ! If output within current range
          if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0d0)then

            ! Get scaled distance and its powers
            U = (xx(j) - x(i))/(x(i) - x(i+1))
            UU = U*U
            UUU = UU*U

            ! Cubic bezier interpolation
            yy(j) = (1d0 - U)*(1d0 - U)*(1d0 - U)*y(i) + &
                    UUU*y(i+1) + 3d0*U*(1d0 - U)*(1d0-U)*c0(i) + &
                    3d0*UU*(1d0 - U)*c1(i+1)

            ! Advance search index
            indx = j + 1

          end if ! Output within current range

        end do ! From indx to the last output
      end do ! From first to second-to-last input position

      return

      end subroutine Intpol_Cubic_Bezier_C

!#####################################################################
!#####################################################################
!#####################################################################

      !> One-dimensional quadratic bezier interpolation with given
      !! coefficients. Does not have control over out of bounds
      !! coordinates\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays\n
      !!  c0(double(:)): Interpolation coefficient\n
      !!  c1(double(:)): Interpolation coefficient\n
      !! Reference: Auer (2003)
      subroutine Intpol_Quadratic_Bezier_C(x,y,n,xx,yy,nn,c0)

      ! I/O

      integer, intent(in):: n,nn
      double precision, dimension(n), intent(in):: x,y,c0
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local

      integer:: indx,i,j

      double precision:: U,UU


      ! Initialize output
      yy = 0

      ! Initialize search index
      indx = 1

      ! From first to second-to-last input points
      do i=1,n-1

        ! From indx to the last output
        do j=indx,nn

          ! If output within current range
          if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0d0) then

            ! Get scaled distance and its power
            U = (xx(j) - x(i))/(x(i) - x(i+1))
            UU = U*U

            ! Quadratic bezier interpolation
            yy(j) = (1d0 - U)*(1d0 - U)*y(i) + &
                    UU*y(i+1) + c0(i)*2d0*U*(1d0 - U)

            ! Advance search index
            indx = j + 1

          end if ! Output within current range

        end do ! From indx to the last output
      end do ! From first to second-to-last input position

      return

      end subroutine Intpol_Quadratic_Bezier_C

!#####################################################################
!#####################################################################
!#####################################################################

      !> One-dimensional linear interpolation. Does not have control
      !! over out of bounds coordinates\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!     n(integer): Size of x and y arrays\n
      !!  xx(double(:)): Output x axis\n
      !!  yy(double(:)): Output y axis\n
      !!    nn(integer): Size of xx and yy arrays
      subroutine Intpol_Lin(x,y,n,xx,yy,nn)

      ! I/O

      integer, intent(in):: n,nn
      double precision, dimension(n), intent(in):: x,y
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(nn), intent(out):: yy

      ! Local

      integer:: i, j, indx

      double precision:: a, b


      ! Single input point
      if (n.eq.1) then

        ! Take constant value
        yy = y(1)

      ! More than one point
      else if(n.gt.1) then

        ! Initialize search index
        indx = 1

        ! For every input position except the last
        do i=1,n-1

          ! Get linear interpolation coefficients
          a = (y(i+1) - y(i))/(x(i+1) - x(i))
          b = y(i) - a*x(i)

          ! From current index to the last in the output
          do j=indx,nn

            ! If output in current range
            if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0) then

              ! Linear interpolation
              yy(j) = a*xx(j) + b

              ! Advance search index
              indx = j + 1

            end if ! Output in current range

          end do ! From search index to the last in the output
        end do ! Every input but the last

      ! Error
      else

        ! Default zero
        yy = 0d0

      end if ! Input dimension

      return

      end subroutine Intpol_Lin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Linear interpolation for Stokes parameters in wavelength.
      !! Does not have control over out of bounds coordinates\n
      !!     x(double(:)): Input x axis\n
      !!   y(double(:,:)): Input y axis\n
      !!       m(integer): First dimension of yy arrays\n
      !!       k(integer): First dimension of y\n
      !!       n(integer): Size of x and second dimension of y\n
      !!    xx(double(:)): Output x axis\n
      !!  yy(double(:,:)): Output y
      subroutine Intpol_Lin_stk(x,y,m,k,n,xx,yy,nn)

      ! I/O

      integer, intent(in):: m,n,nn,k
      double precision, dimension(n), intent(in):: x
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(k,n), intent(in):: y
      double precision, dimension(m,nn), intent(out):: yy

      ! Local

      integer:: i,j,indx

      double precision, dimension(k):: a,b


      ! Single point
      if (n.eq.1) then

        ! For every output point
        do i=1,nn

          ! Make equal to the single input data
          yy(1:k,i) = y(:,1)

          ! Remaining first dimension set to zero
          yy(k+1:m,i) = 0d0

        end do ! Output points

      ! More than one point
      else if (n.gt.1) then

        ! Initialize search index
        indx = 1

        ! For all points but the last
        do i=1,n-1

          ! Linear interpolation coefficients
          a = (y(:,i+1) - y(:,i))/(x(i+1) - x(i))
          b = y(:,i) - a*x(i)

          ! From the current search index to the last output
          do j=indx,nn

            ! If in current range
            if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0) then

              ! Linear interpolation
              yy(1:k,j) = a*xx(j) + b

              ! Remaining first dimension set to zero
              yy(k+1:m,j) = 0d0

              ! Advance search index
              indx = j + 1

            end if ! In current range

          end do ! From current search index to the last output
        end do ! All input points but the last

      ! Error
      else

        ! Default to zero
        yy = 0d0

      end if ! Input dimension

      return

      end subroutine Intpol_Lin_stk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Linear interpolation for multiple data arrays. Does not have
      !! control over out of bounds coordinates\n
      !!     x(double(:)): Input x axis\n
      !!   y(double(:,:)): Input y axis\n
      !!       m(integer): First dimension of y\n
      !!       n(integer): Size of x and second dimension of y\n
      !!    xx(double(:)): Output x axis\n
      !!  yy(double(:,:)): Output y
      subroutine Intpol_Lin_mD(x,y,m,n,xx,yy,nn)

      ! I/O

      integer, intent(in):: m,n,nn
      double precision, dimension(n), intent(in):: x
      double precision, dimension(nn), intent(in):: xx
      double precision, dimension(m,n), intent(in):: y
      double precision, dimension(m,nn), intent(out):: yy

      ! Local

      integer:: i,j,indx

      double precision, dimension(m):: a,b


      ! Single point
      if (n.eq.1) then

        ! For every output
        do i=1,nn

          ! Set equal to input
          yy(:,i) = y(:,1)

        end do ! Output points

      ! More than one point
      else if(n.gt.1) then

        ! Initialize search index
        indx = 1

        ! For every input but the last
        do i=1,n-1

          ! Linear interpolation coefficients
          a = (y(:,i+1) - y(:,i))/(x(i+1) - x(i))
          b = y(:,i) - a*x(i)

          ! From the current search index to the last output
          do j=indx,nn

            ! If in current range
            if ((xx(j) - x(i))*(xx(j) - x(i+1)).le.0) then

              ! Linear interpolation
              yy(:,j) = a*xx(j) + b

              ! Advance search index
              indx = j + 1

            end if ! In current range

          end do ! From current search index to the last output
        end do ! All input points but the last

      ! Error
      else

        ! Default to zero
        yy = 0d0

      end if ! Input dimension

      return

      end subroutine Intpol_Lin_mD

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the coordinate value which gives the minimum value
      !! for a parabolic interpolation\n
      !!   x(double(:)): Input x axis\n
      !!   y(double(:)): Input y axis\n
      !!  indx(integer): Index of minimum y\n
      !!  x_min(double): x value for minimum y
      subroutine parabolic(x,y,indx,x_min)

      ! I/O

      integer, intent(in):: indx
      double precision, intent(out):: x_min
      double precision, dimension(:), intent(in):: x,y

      ! Local

      integer:: length

      double precision:: x_tmp1,x_tmp2,y_tmp1,y_tmp2
      double precision, dimension(:), allocatable:: xp


      ! Get size of x axis
      length = size(x)

      ! Allocate space of same length
      allocate(xp(length))

      ! Calculate logarithm from beginning to the position with
      ! current minimum value
      xp(1:indx+1) = log(x(1:indx+1))

      ! Calculate steps in x and y around the current minimum
      x_tmp1 = xp(indx) - xp(indx-1)
      x_tmp2 = xp(indx) - xp(indx+1)
      y_tmp1 = y(indx) - y(indx-1)
      y_tmp2 = y(indx) - y(indx+1)

     !a = (Y_tmp2/X_tmp2-Y_tmp1/X_tmp1)/(XP(Indx+1)-XP(Indx-1))
     !b = Y_tmp1/X_tmp1-a*(XP(Indx)+XP(Indx-1))
     !c = Y(Indx)-a*XP(Indx)**2-b*XP(Indx)

      ! Get coordinate for the minimum value of the interpolated
      ! parabola in logarithm
      x_min = xp(indx) - 0.5d0*(x_tmp1*x_tmp1*y_tmp2 - &
                                x_tmp2*x_tmp2*y_tmp1)/ &
                         (x_tmp1*y_tmp2 - x_tmp2*y_tmp1)

      ! Transform to linear scale
      x_min = exp(x_min)

      ! Free space
      deallocate(xp)

      return

      end subroutine Parabolic

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate derivatives for several interpolation functions in
      !! the module\n
      !!   H0(double): Back step\n
      !!   H1(double): Forward step\n
      !!   D0(double): Back slope\n
      !!   D1(double): Forward slope\n
      !!  YP1(double): Output value\n
      !! Reference: Auer (2003)
      function DERIV(H0,H1,D0,D1) result(YP1)

      ! I/O

      double precision:: H0,H1,D0,D1
      double precision:: YP1

      ! Local

      double precision:: ALPHA


      ! If same sign slope
      if (D0*D1.gt.0) then

        ! Get derivative
        ALPHA = (1d0 + H1/(H0+H1))/3d0
        YP1 = (D0*D1)/(ALPHA*D1 + (1d0 - ALPHA)*D0)

      ! If different sign
      else

        ! Extrema
        YP1 = 0d0

      end if ! Signs of slopes

      end function DERIV

!#####################################################################
!#####################################################################
!#####################################################################

      end module inter_mod
