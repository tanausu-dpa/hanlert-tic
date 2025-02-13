      !> Auxiliar routines for short characteristics
      module rtstepaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     27/04/2017
!  Last version:
!     19/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     19/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  RTomega
!    Calculate the coefficients for the integral of the source
!  function in the short-characteristics with BESSER
!
!  psi_lin
!    Calculate the coefficients for the integral of the source
!  function in the short-characteristics with linear interpolation
!
!  psi_par
!    Calculate the coefficients for the integral of the source
!  function in the short-characteristics with parabolic interpolation
!
!  ftau
!   Calculate the optical distance between two points
!
!  QBezierC0
!    Calculate the correction parameter for the integral of the
!  source function in the short-characteristics with BESSER
!
!  ybetwab
!    Determine if inputs (a,y,b) are in monotonic order
!
!  correctyab
!    Determine the y which makes (a,y,b) monotonic
!
!  MatVec_2d
!    Calculate c = A*b with "A" a 2x2 matrix and "b" a 2 vector
!
!  MatVec
!    Calculate c = A*b with "A" a 4x4 matrix and "b" a 4 vector
!
!  Matinv_2d
!    Calculate the inverse of a 2x2 matrix with the symmetry
!  properties of the propagation matrix for only I and Q
!
!  Matinv
!    Calculate the inverse of a 4x4 matrix with the symmetry
!  properties of the propagation matrix
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only : vacuum

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the coefficients for the integral of the source
      !! function in the short-characteristics with BESSER\n
      !!       ex(double): Exponencial of the optical distance\n
      !!        t(double): Optical distance\n
      !!  omega_m(double): First omega parameter\n
      !!  omega_o(double): Second omega parameter\n
      !!  omega_c(double): Third omega parameter
      subroutine RTomega(ex,t,omega_m,omega_o,omega_c)

      ! I/O

      double precision, intent(in):: ex,t
      double precision, intent(out):: omega_m,omega_o,omega_c

      ! If large enough optical distance
      if (t.gt.0.14d0) then

        ! Compute
        omega_m = ((2d0 - ex*(t*t + 2d0*t + 2d0))/(t*t))

      ! Small optical distance
      else

        ! Compute
        omega_m = ((t*(t*(t*(t*(t*(t*((140d0 - 18d0*t)*t - &
                  945d0) + 5400d0) - 25200d0) + 90720d0) - &
                  226800d0) + 302400d0))/907200d0)

      end if ! Optical distance


      ! If large enough optical distance
      if (t.gt.0.18d0) then

        ! Compute
        omega_o = (1d0 - 2d0*(t + ex - 1d0)/(t*t))
        omega_c = (2d0*(t - 2d0 + ex*(t + 2d0))/(t*t))

      ! Small optical distance
      else

        ! Compute
        omega_o = ((t*(t*(t*(t*(t*(t*((10d0 - t)*t - 90d0) + &
                  720d0) - 5040d0) + 30240d0) - 151200d0) + &
                  604800d0))/1814400d0)
        omega_c = ((t*(t*(t*(t*(t*(t*((35d0 - 4d0*t)*t - &
                  270d0) + 1800d0) - 10080d0) + 45360d0) - &
                  151200d0) + 302400d0))/907200d0)

      end if ! Optical distance

      end subroutine RTomega

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the coefficients for the integral of the source
      !! function in the short-characteristics with linear
      !! interpolation\n
      !!    ex(double): Exponencial of the optical distance\n
      !!     t(double): Optical distance\n
      !!  psim(double): First coefficient\n
      !!  psio(double): Second coefficient
      subroutine psi_lin(ex,t,psim,psio)

      ! I/O

      double precision, intent(in):: ex,t
      double precision, intent(out):: psim,psio

      ! If optical distance is large enough
      if (t.gt.0.11d0) then

        ! Compute
        psim = ((1d0-ex*(1d0+t))/t)
        psio = ((ex+t-1d0)/t)

      ! If optical distance is small enough
      else

        ! Compute
        psim = &
          ((t*(t*(t*(t*(t*(t*((63d0 - 8d0*t)*t - 432d0) + 2520d0) - &
            12096d0) + 45360d0) - 120960d0) + 181440d0))/362880d0)
        psio = ((t*(t*(t*(t*(t*(t*((9d0 - t)*t - 72d0) + 504d0) - &
               3024d0) + 15120d0) - 60480d0) + 181440d0))/362880d0)

      end if ! Optical distance

      end subroutine psi_lin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the coefficients for the integral of the source
      !! function in the short-characteristics with parabolic
      !! interpolation\n
      !!    ex(double): Exponencial of the optical distance\n
      !!    tm(double): Optical distance M-O\n
      !!    tp(double): Optical distance O-P\n
      !!  psim(double): First coefficient\n
      !!  psio(double): Second coefficient\n
      !!  psip(double): Third coefficient
      subroutine psi_par(ex,tm,tp,psim,psio,psip)

      ! I/O

      double precision, intent(in):: ex,tm,tp
      double precision, intent(out):: psim,psio,psip

      ! Local

      double precision:: w0, w1, w2


      ! Initialize
      w0 = 1d0 - ex
      w1 = w0 - tm*ex
      w2 = 2d0*w1 - tm*tm*ex

      ! Large enough optical distance
      if (tm.gt.0.11d0) then

        ! Compute
        psim = (w2 + w1*tp)/(tm*(tm + tp))
        psio = w0 + (w1*(tm - tp) - w2)/(tm*tp)
        psip = (w2 - w1*tm)/(tp*(tm + tp))

      ! Small optical distance
      else

        ! Compute
        psim = (tm*tm*(tm*(tm*(tm*(tm*(tm*(7d0*tm - 48d0) + 280.0) - &
                       1344d0) + 5040d0) - 13440d0) + 20160d0)*tp + &
                tm*tm*tm*(tm*(tm*(tm*((240d0 - 42d0*tm)*tm - &
                                      1120d0) + 4032d0) - 10080d0) + &
                          13440d0)) / 40320d0/(tm*(tm + tp))
        psio = (tm*(tm*(tm*(tm*(tm*(tm*((8d0 - tm)*tm - 56d0) + &
                                    336d0) - 1680d0) + 6720d0) - &
                        20160d0) + 40320d0)) / 40320d0 + &
               ((tm*tm*(tm*(tm*(tm*(tm*((48d0 - 7d0*tm)*tm - &
                                        280d0) + 1344d0) - &
                                5040d0) + 13440d0) - 20160d0)*tp + &
                 tm*tm*tm*(tm*(tm*(tm*((40d0 - 6d0*tm)*tm - &
                                       224d0) + 1008d0) - 3360d0) + &
                           6720d0)) / 40320d0)/(tm*tp)
        psip = ((tm*tm*tm*(tm*(tm*(tm*(tm*(3d0*tm - 20d0) + &
                                       112d0) - 504d0) + 1680d0) - &
                           3360d0)) / 20160d0) / (tp*(tm + tp))

      end if ! Optical distance

      end subroutine psi_par

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the optical distance between two points\n
      !!     opx(double): Opacity point 1\n
      !!     opu(double): Opacity point 2\n
      !!      ds(double): Geometrical distance\n
      !!  output(double): Optical distance
      subroutine ftau(opx,opu,ds,output)

      ! I/O

      double precision, intent(in):: opx,opu,ds
      double precision, intent(out):: output

      ! Linear interpolation
      output = 0.5d0*ds*(opx + opu + 2d0*vacuum)

      end subroutine ftau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the correction parameter for the integral of the
      !! source function in the short-characteristics with BESSER\n
      !!      h0(double): Optical distance backward\n
      !!      h1(double): Optical distance forward\n
      !!      ym(double): Source function previous point\n
      !!      yo(double): Source function current point\n
      !!      yp(double): Source function next point\n
      !!  output(double): C0 parameter
      subroutine QBezierC0(h0,h1,ym,yo,yp,output)

      ! I/O

      double precision, intent(in):: h0,h1,ym,yo,yp
      double precision, intent(out):: output

      ! Local

      logical:: cond0,cond1,cond2

      double precision:: dm,dp,yder,c0,c1,c2


      ! Both optical distances are positive and non-zero
      if (h0.gt.0d0.and.h1.gt.0d0) then

        ! Get derivatives
        dm = (yo - ym)/h0
        dp = (yp - yo)/h1

      ! One of them is zero
      else

        ! Keep Source at current position
        output = yo
        return

      endif ! Positive optical distances

      ! If different slopes
      if (dm*dp.le.0d0) then

        ! Keep Source at current position
        output = yo
        return

      endif ! Different slopes

      ! Get derivative
      yder = 0.5d0*(h0*dp + h1*dm)/(h0 + h1)

      ! Get correction to both sides
      c0 = yo - h0*yder
      c1 = yo + h1*yder

      ! Check if monotonic
      call ybetwab(c0,ym,yo,cond0)
      call ybetwab(c1,yo,yp,cond1)

      ! If both monotonic
      if (cond0.and.cond1) then

        ! Get borrection
        output = c0
        return

      ! If c0 non-monotonic
      elseif (.not.cond0) then

        ! Correct C0
        call correctyab(c0,ym,yo,output)

      ! If c1 non-monotonic
      elseif(.not.cond1)then

        ! Correct c1
        call correctyab(c1,yo,yp,c2)

        ! Correct c0
        yder = (c2 - yo)/h1
        c0 = yo - h0*yder

        ! Check if not it is monotonic
        call ybetwab(c0,ym,yo,cond2)

        ! If c0 non-monotonic
        if(.not.cond2)then

          ! Correct c0
          call correctyab(c0,ym,yo,output)

        ! If c0 monotonic
        else

          ! Get this c0
          output = c0

        endif ! Corrected c0 monotonic
      endif ! Monotonic

      return

      end subroutine QBezierC0

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine if inputs (a,y,b) are in monotonic order\n
      !!         y(double): Input central number\n
      !!         a(double): Input extreme number\n
      !!         b(double): Input extreme number\n
      !!   output(logical): If (a,y,b) is monotonic
      subroutine ybetwab(y,a,b,output)

      ! I/O
      logical, intent(out):: output
      double precision, intent(in):: y,a,b


      ! Check monotonic
      output = (a.le.b.and.y.ge.a.and.y.le.b).or. &
               (a.ge.b.and.y.le.a.and.y.ge.b)

      return

      end subroutine ybetwab

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine the y which makes (a,y,b) monotonic\n
      !!       y(double): Input central number\n
      !!       a(double): Input extreme number\n
      !!       b(double): Input extreme number\n
      !!  output(double): y that makes (a,y,b) monotonic
      subroutine correctyab(y,a,b,output)

      ! I/O

      double precision, intent(in):: y,a,b
      double precision, intent(out):: output

      ! Local

      double precision:: mini,maxi


      ! Get extremes
      mini = minval((/ a, b /))
      maxi = maxval((/ a, b /))

      ! If y between minimum and maximum
      if (y.ge.mini.and.y.le.maxi) then

        ! Was correct
        output = y

      ! If below minimum
      else if (y.lt.mini) then

        ! Set to minimum
        output = mini

      ! If above maximum
      else

        ! Set to maximum
        output = maxi

      end if ! If in range our outside which range

      return

      end subroutine correctyab

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate c = A*b with "A" a 2x2 matrix and "b" a 2 vector\n
      !!  A(double(:,:)): Matrix\n
      !!    b(double(:)): Vector\n
      !!    c(double(:)): Vector A*b
      subroutine MatVec_2d(A,b,c)

      ! I/O

      double precision, dimension(:,:), intent(in):: A
      double precision, dimension(:), intent(in):: b
      double precision, dimension(:), intent(out):: c

      ! Product
      c(1) = A(1,1)*b(1) + A(1,2)*b(2)
      c(2) = A(2,1)*b(1) + A(2,2)*b(2)

      end subroutine MatVec_2d

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate c = A*b with "A" a 4x4 matrix and "b" a 4 vector\n
      !!  A(double(:,:)): Matrix\n
      !!    b(double(:)): Vector\n
      !!    c(double(:)): Vector A*b
      subroutine MatVec(A,b,c)

      ! I/O

      double precision, dimension(:,:), intent(in):: A
      double precision, dimension(:), intent(in):: b
      double precision, dimension(:), intent(out):: c

      ! Product
      c(1) = A(1,1)*b(1) + A(1,2)*b(2) + A(1,3)*b(3) + A(1,4)*b(4)
      c(2) = A(2,1)*b(1) + A(2,2)*b(2) + A(2,3)*b(3) + A(2,4)*b(4)
      c(3) = A(3,1)*b(1) + A(3,2)*b(2) + A(3,3)*b(3) + A(3,4)*b(4)
      c(4) = A(4,1)*b(1) + A(4,2)*b(2) + A(4,3)*b(3) + A(4,4)*b(4)

      end subroutine MatVec

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the inverse of a 2x2 matrix with the symmetry
      !! properties of the propagation matrix for only I and Q\n
      !!  M(double(:,:)): Matrix to invert and inverted matrix
      subroutine Matinv_2d(M)

      ! I/O

      double precision, dimension(2,2), intent(inout):: M

      ! Local

      double precision:: idet

      ! Inverse determinant
      idet = 1d0/(1d0 - M(1,2)*M(1,2))

      ! Get inverse
      M(1,2) = -M(1,2)
      M(2,1) =  M(1,2)

      ! Scale
      M = M*idet

      return

      end subroutine Matinv_2d

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the inverse of a 4x4 matrix with the symmetry
      !! properties of the propagation matrix\n
      !!  M(double(:,:)): Matrix to invert and inverted matrix
      subroutine Matinv(M)

      ! I/O

      double precision, dimension(4,4), intent(inout):: M

      ! Local

      double precision:: a,b,c,d,e,f,o
      double precision:: aa,ab,ac,ad,ae,af
      double precision:: bb,bc,bd,be,bf
      double precision:: cc,cd,ce,cf
      double precision:: dd,de,df
      double precision:: ee,ef
      double precision:: ff
      double precision:: T0
      double precision:: P1,P2,P3,P4,P5,P6
      double precision:: S1,S2,S3,S4,S5,S6
      double precision:: idet

      ! Get minimal products
      a = M(1,2)
      b = M(1,3)
      c = M(1,4)
      d = M(2,3)
      e = M(2,4)
      f = M(3,4)
      o = 1d0
      aa = a*a
      ab = a*b
      ac = a*c
      ad = a*d
      ae = a*e
      af = a*f
      bb = b*b
      bc = b*c
      bd = b*d
      be = b*e
      bf = b*f
      cc = c*c
      cd = c*d
      ce = c*e
      cf = c*f
      dd = d*d
      de = d*e
      df = d*f
      ee = e*e
      ef = e*f
      ff = f*f
      T0 = (af - be + cd)
      P1 = T0*a
      P2 = T0*b
      P3 = T0*c
      P4 = T0*d
      P5 = T0*e
      P6 = T0*f
      S1 = bc - de
      S2 = ac + df
      S3 = ab - ef
      S4 = ae + bf
      S5 = ad - cf
      S6 = bd + ce

      ! Get inverse
      M(1,1) =   o + dd + ee + ff
      M(2,1) = - a + S6 - P6
      M(3,1) = - b - S5 + P5
      M(4,1) = - c - S4 - P4
      M(1,2) = - a - S6 - P6
      M(2,2) =   o - bb - cc + ff
      M(3,2) =   d + S3 - P3
      M(4,2) =   e + S2 + P2
      M(1,3) = - b + S5 + P5
      M(2,3) = - d + S3 + P3
      M(3,3) =   o - aa - cc + ee
      M(4,3) =   f + S1 - P1
      M(1,4) = - c + S4 - P4
      M(2,4) = - e + S2 - P2
      M(3,4) = - f + S1 + P1
      M(4,4) =   o - aa - bb + dd

      ! Get determinant
      idet = 1d0/(M(1,1) + a*M(2,1) + b*M(3,1) + c*M(4,1))

      ! Scale inverse
      M = M*idet

      return

      end subroutine Matinv

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtstepaux_mod
