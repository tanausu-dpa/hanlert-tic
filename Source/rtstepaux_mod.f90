      !> Auxiliar routines for short characteristics
      module rtstepaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/27/2017
!  Last version:
!     09/21/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/21/2023:    V3.0.2 - Added psi_par, MatVec_2d, and Matinv_2d
!                             routines (TdPA)
!
!     10/27/2022:    V3.0.1 - Introduced a factor 1/2 in yder
!                             definition in QBezierC0 to save between
!                             one and three products (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     02/20/2019:    V1.1.0 - Does not need aborted dependence, nor
!                             TINY variable (TdPA)
!
!     04/27/2017:    V1.0.0 - First version (TdPA)
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
!  RTomega:
!    Coefficients for short characteristics with besser interpolation
!
!  psi_lin:
!    Coefficients of short characteristics with linear interpolation
!
!  psi_par:
!    Coefficients of short characteristics with parabolic
!  interpolation
!
!  ftau:
!   Calculate optical distance between two points
!
!  QBezierC0
!    Coefficient for the P part of short characteristics in Besser
!  interpolation
!
!  ybetwab
!    Returns True if monotonic a-y-b
!
!  correctyab:
!    Takes the y that keeps the a,y,b monotonic
!
!  MatVec_2d:
!    Computes c = a*b with a matrix and b vector of size 2
!
!  MatVec:
!    Computes c = a*b with a matrix and b vector of size 4
!
!  Matinv_2d:
!    Explicit inverse of a 2x2 matrix with certain symmetry properties
!
!  Matinv:
!    Explicit inverse of a 4x4 matrix with certain symmetry properties
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

      !> Computes omega coefficients for BESSER short
      !! characteristics\n
      !!       ex(dfloat): Exponencial of the optical distance\n
      !!        t(dfloat): Optical distance\n
      !!  omega_m(dfloat): First omega parameter\n
      !!  omega_o(dfloat): Second omega parameter\n
      !!  omega_c(dfloat): Third omega parameter
      subroutine RTomega(ex,t,omega_m,omega_o,omega_c)

      ! I/O

      double precision, intent(in):: ex,t
      double precision, intent(out):: omega_m,omega_o,omega_c

      if (t.gt.0.14d0) then

        omega_m = ((2d0 - ex*(t*t + 2d0*t + 2d0))/(t*t))

      else

        omega_m = ((t*(t*(t*(t*(t*(t*((140d0 - 18d0*t)*t - &
                  945d0) + 5400d0) - 25200d0) + 90720d0) - &
                  226800d0) + 302400d0))/907200d0)

      end if


      if (t.gt.0.18d0) then

        omega_o = (1d0 - 2d0*(t + ex - 1d0)/(t*t))
        omega_c = (2d0*(t - 2d0 + ex*(t + 2d0))/(t*t))

      else

        omega_o = ((t*(t*(t*(t*(t*(t*((10d0 - t)*t - 90d0) + &
                  720d0) - 5040d0) + 30240d0) - 151200d0) + &
                  604800d0))/1814400d0)
        omega_c = ((t*(t*(t*(t*(t*(t*((35d0 - 4d0*t)*t - &
                  270d0) + 1800d0) - 10080d0) + 45360d0) - &
                  151200d0) + 302400d0))/907200d0)

      end if

      end subroutine RTomega

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes coefficients for linear short characteristics\n
      !!       ex(dfloat): Exponencial of the optical distance\n
      !!        t(dfloat): Optical distance\n
      !!     psim(dfloat): First coefficient\n
      !!     psio(dfloat): Second coefficient
      subroutine psi_lin(ex,t,psim,psio)

      ! I/O

      double precision, intent(in):: ex,t
      double precision, intent(out):: psim,psio

      if (t.gt.0.11d0) then

        psim = ((1d0-ex*(1d0+t))/t)
        psio = ((ex+t-1d0)/t)

      else

        psim = &
          ((t*(t*(t*(t*(t*(t*((63d0 - 8d0*t)*t - 432d0) + 2520d0) - &
            12096d0) + 45360d0) - 120960d0) + 181440d0))/362880d0)
        psio = ((t*(t*(t*(t*(t*(t*((9d0 - t)*t - 72d0) + 504d0) - &
               3024d0) + 15120d0) - 60480d0) + 181440d0))/362880d0)

      end if

      end subroutine psi_lin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes coefficients for parabolic short characteristics\n
      !!       ex(dfloat): Exponencial of the optical distance\n
      !!       tm(dfloat): Optical distance M-O\n
      !!       tp(dfloat): Optical distance O-P\n
      !!     psim(dfloat): First coefficient\n
      !!     psio(dfloat): Second coefficient
      !!     psip(dfloat): Second coefficient
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

      if (tm.gt.0.11d0) then

        psim = (w2 + w1*tp)/(tm*(tm + tp))
        psio = w0 + (w1*(tm - tp) - w2)/(tm*tp)
        psip = (w2 - w1*tm)/(tp*(tm + tp))

      else

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

      end if

      end subroutine psi_par

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes optical distances\n
      !!     opx(dfloat): Opacity point 1\n
      !!     opu(dfloat): Opacity point 2\n
      !!      ds(dfloat): Geometrical distance\n
      !!  output(dfloat): Optical distance
      subroutine ftau(opx,opu,ds,output)

      ! I/O

      double precision, intent(in):: opx,opu,ds
      double precision, intent(out):: output

      output = 0.5d0*ds*(opx + opu + 2d0*vacuum)

      end subroutine ftau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes C0 parameter of BESSER algorithm\n
      !!      h0(dfloat): Optical distance backward\n
      !!      h1(dfloat): Optical distance forward\n
      !!      ym(dfloat): Source function previous point\n
      !!      yo(dfloat): Source function current point\n
      !!      yp(dfloat): Source function next point\n
      !!  output(dfloat): C0 parameter
      subroutine QBezierC0(h0,h1,ym,yo,yp,output)

      ! I/O

      double precision, intent(in):: h0,h1,ym,yo,yp
      double precision, intent(out):: output

      ! Local

      logical:: cond0,cond1,cond2

      double precision:: dm,dp,yder,c0,c1,c2


      if (h0.gt.0d0.and.h1.gt.0d0) then

        dm = (yo - ym)/h0
        dp = (yp - yo)/h1

      else

        output = yo
        return

      endif

      if (dm*dp.le.0d0) then

        output = yo
        return

      endif

      yder = 0.5d0*(h0*dp + h1*dm)/(h0 + h1)
      c0 = yo - h0*yder
      c1 = yo + h1*yder

      call ybetwab(c0,ym,yo,cond0)
      call ybetwab(c1,yo,yp,cond1)

      if(cond0.and.cond1)then

        output = c0
        return

      elseif(.not.cond0)then

        call correctyab(c0,ym,yo,output)

      elseif(.not.cond1)then

        call correctyab(c1,yo,yp,c2)
        yder = (c2 - yo)/h1
        c0 = yo - h0*yder
        call ybetwab(c0,ym,yo,cond2)

        if(.not.cond2)then

          call correctyab(c0,ym,yo,output)

        else

          output = c0

        endif

      endif

      return

      end subroutine QBezierC0

!#####################################################################
!#####################################################################
!#####################################################################

      !> Returns true if monotonic (a,y,b)\n
      !!     y,a,b(dfloat): Input numbers\n
      !!   output(logical): True if (a,y,b) monotonic
      subroutine ybetwab(y,a,b,output)

      ! I/O

      logical, intent(out):: output
      double precision, intent(in):: y,a,b

      output = (a.le.b.and.y.ge.a.and.y.le.b).or. &
               (a.ge.b.and.y.le.a.and.y.ge.b)

      return

      end subroutine ybetwab

!#####################################################################
!#####################################################################
!#####################################################################

      !> Takes the y that keeps (a,y,b) monotonic\n
      !!     y,a,b(dfloat): Input numbers\n
      !!   output(logical): y that makes (a,y,b) monotonic
      subroutine correctyab(y,a,b,output)

      ! I/O

      double precision, intent(in):: y,a,b
      double precision, intent(out):: output

      ! Local

      double precision:: mini,maxi

      mini = minval((/ a, b /))
      maxi = maxval((/ a, b /))

      if (y.ge.mini.and.y.le.maxi) then

        output = y

      else if (y.lt.mini) then

        output = mini

      else

        output = maxi

      end if

      return

      end subroutine correctyab

!#####################################################################
!#####################################################################
!#####################################################################

      !> Product matrix and vector 2x2\n
      !!    a(dfloat(:,:)): Matrix\n
      !!      b(dfloat(:)): Vector\n
      !!      c(dfloat(:)): Vector A*b
      subroutine MatVec_2d(a,b,c)

      ! I/O

      double precision, dimension(:,:), intent(in):: a
      double precision, dimension(:), intent(in):: b
      double precision, dimension(:), intent(out):: c

      c(1) = a(1,1)*b(1) + a(1,2)*b(2)
      c(2) = a(2,1)*b(1) + a(2,2)*b(2)

      end subroutine MatVec_2d

!#####################################################################
!#####################################################################
!#####################################################################

      !> Product matrix and vector 4x4\n
      !!    a(dfloat(:,:)): Matrix\n
      !!      b(dfloat(:)): Vector\n
      !!      c(dfloat(:)): Vector A*b
      subroutine MatVec(a,b,c)

      ! I/O

      double precision, dimension(:,:), intent(in):: a
      double precision, dimension(:), intent(in):: b
      double precision, dimension(:), intent(out):: c

      c(1) = a(1,1)*b(1) + a(1,2)*b(2) + a(1,3)*b(3) + a(1,4)*b(4)
      c(2) = a(2,1)*b(1) + a(2,2)*b(2) + a(2,3)*b(3) + a(2,4)*b(4)
      c(3) = a(3,1)*b(1) + a(3,2)*b(2) + a(3,3)*b(3) + a(3,4)*b(4)
      c(4) = a(4,1)*b(1) + a(4,2)*b(2) + a(4,3)*b(3) + a(4,4)*b(4)

      end subroutine MatVec

!#####################################################################
!#####################################################################
!#####################################################################

      !> Inverts a matrix 2x2 with the absorption matrix symmetry
      !! properties\n
      !!     M(dfloat(:,:)): Matrix to invert and inverted matrix
      subroutine Matinv_2d(M)

      ! I/O

      double precision, dimension(2,2), intent(inout):: M

      ! Local

      double precision:: idet

      ! Inverse determinant
      idet = 1d0/(1d0 - M(1,2)*M(1,2))

      M(1,2) = -M(1,2)
      M(2,1) =  M(1,2)

      M = M*idet

      return

      end subroutine Matinv_2d

!#####################################################################
!#####################################################################
!#####################################################################

      !> Inverts a matrix 4x4 with the absorption matrix symmetry
      !! properties\n
      !!     M(dfloat(:,:)): Matrix to invert and inverted matrix
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

      idet = 1d0/(M(1,1) + a*M(2,1) + b*M(3,1) + c*M(4,1))

      M = M*idet

      return

      end subroutine Matinv

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtstepaux_mod
