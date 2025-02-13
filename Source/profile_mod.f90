      !> Module for Voigt and redistribution profiles
      module profile_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     20/04/2017
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  voigt
!    Calculate Voigt and Faraday-Voigt profiles. Real algebra version
!
!  Wfunc
!    Calculate the redistribution function in the laboratory frame
!    with non-coherent lower term and infinitely sharp lower levels
!
!  voigtI
!    Calculate the real Voigt profile
!
!  voigtI_A
!    Calculate the real Voigt profile with Armstrong (1967) algorithm
!
!  voigtI_AK1
!     Branch 1 for Armstrong (1967) algorithm
!
!  voigtI_AK2
!     Branch 2 for Armstrong (1967) algorithm
!
!  voigtI_AK3
!     Branch 3 for Armstrong (1967) algorithm
!
!  voigtI_H
!    Calculate the real Voigt profile with Hui et al. (1978) algorithm
!
!  voigtI_HC
!     Calculate the real Voigt profiles with Humlicek (1982)
!  algorithm. Real algebra version
!
!  gaussianI
!     Calculate a Gaussian profile
!
!  WfuncI
!    Calculate the real part of the redistribution function in the
!    laboratory frame with non-coherent lower term and infinitely
!    sharp lower levels
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use parameters_mod , only : PI , IPI, cImag , Wbiggauss , &
                                  WbiggaussI , smallexp , bigexp

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate Voigt and Faraday-Voigt profiles. Real algebra
      !! version\n
      !!    v(double): Normalized frequency\n
      !!    a(double): Normalized damping parameter\n
      !!   t(dcomplx): Complex Voigt and Faraday-Voigt profiles
      subroutine voigt(v,a,t)

      ! I/O

      double precision, intent(in):: v, a
      complex(kind=8), intent(out):: t

      ! Local

      real:: sa, sv, s, d

      double precision:: r1,i1,r2,i2,ru,iu,rx,ix,de

      ! Get single precision argumens
      sa = sngl(a)
      sv = sngl(v)

      ! Get parameters to decide region
      s = abs(sv)+sa
      d = .195e0*abs(sv)-.176e0

      ! First branch
      if (s.ge..15e2) then

        ! Output
        ! t = .5641896d0*z/(.5d0+z*z)

        ! (0.5 + z^2)* <- ans
        r1 = 0.5d0 + a*a - v*v
        i1 = 2d0*a*v

        ! Denominator
        de = .5641896d0/(r1*r1 + i1*i1)

        ! Output ans*z*de
        t = dcmplx(a*r1 + v*i1,a*i1 - v*r1)*de

      ! Other branches
      else

        ! Second branch
        if (s.ge..55e1) then

          ! Output
          ! u = z*z
          ! t = z*(.1410474d1 + .5641896d0*u)/(.75d0 + u*(.3d1 + u))

          ! u = z*z
          ru = a*a - v*v
          iu = -2d0*a*v

          ! 0.75 + u(3 + u) <- ans
          r1 =  0.75d0 + ru*(ru + 3d0) - iu*iu
          i1 =  -iu*(3d0 + 2d0*ru)

          ! Denominator
          de = 1d0/(r1*r1 + i1*i1)

          ! z*conj(ans) <- ans
          r2 = a*r1 + v*i1
          i2 = a*i1 - v*r1

          ! .1410474d1 + .5641896d0*u <- bns
          ru = ru*.5641896d0 + .1410474d1
          iu = iu*.5641896d0

          ! Output ans*bns*de
          t = dcmplx(r2*ru - i2*iu,r2*iu + ru*i2)*de

        ! Other branches
        else

          ! Thid branch
          if (sa.ge.d) then

            ! Output
            ! nt = .164955d2 + z*(.2020933d2 + z*(.1196482d2 + &
            !      z*(.3778987d1 + .5642236d0*z)))
            ! dt = .164955d2 + z*(.3882363d2 + z*(.3927121d2 + &
            !      z*(.2169274d2 + z*(.6699398d1 + z))))
            ! t = nt/dt

            ! .3778987d1 + .5642236d0*z <- ans
            r1 = .3778987d1 + .5642236d0*a
            i1 = -.5642236d0*v

            ! .1196482d2 + z*ans <- ans
            r2 = a*r1 + v*i1 + .1196482d2
            i2 = a*i1 - v*r1

            ! .2020933d2 + z*ans <- ans
            r1 = a*r2 + v*i2 + .2020933d2
            i1 = a*i2 - v*r2

            ! u = .164955d2 + z*ans
            ru = a*r1 + v*i1 + .164955d2
            iu = a*i1 - v*r1

            !
            ! .6699398d1 + z <- ans
            r1 = .6699398d1 + a
            i1 = -v

            ! .2169274d2 + z*ans <- ans

            r2 = a*r1 + v*i1 + .2169274d2
            i2 = a*i1 - v*r1

            ! .3927121d2 + z*ans <- ans
            r1 = a*r2 + v*i2 + .3927121d2
            i1 = a*i2 - v*r2

            ! .3882363d2 + z*ans <- ans
            r2 = a*r1 + v*i1 + .3882363d2
            i2 = a*i1 - v*r1

            ! conj(.164955d2 + z*ans) <- ans
            r1 =   a*r2 + v*i2 + .164955d2
            i1 = -(a*i2 - v*r2)

            ! Denominator
            de = 1d0/(r1*r1 + i1*i1)

            ! Output ans*u*de
            t = dcmplx(ru*r1 - iu*i1,ru*i1 + r1*iu)*de

          ! Fourth branch
          else

            ! Output
            ! u = z*z
            ! x = z*(.3618331d5 - u*(.33219905d4 - u*(.1540787d4 - &
            !     u*(.2190313d3 - u*(.3576683d2 - u*(.1320522d1 - &
            !     .56419d0*u))))))
            ! y = .320666d5 - u*(.2432284d5 - u*(.9022228d4 - &
            !     u*(.2186181d4 - u*(.3642191d3 - u*(.6157037d2 - &
            !     u*(.1841439d1 - u))))))
            ! t = exp(u) - x/y

            ! u = z*z
            ru =  a*a - v*v
            iu = -2d0*a*v

            !
            ! .1320522d1 - .56419d0*u <- ans
            r1 =  .1320522d1 - .56419d0*ru
            i1 = -.56419d0*iu

            ! .3576683d2 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .3576683d2
            i2 = -(ru*i1 + iu*r1)

            ! .2190313d3 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .2190313d3
            i1 = -(ru*i2 + iu*r2)

            ! .1540787d4 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .1540787d4
            i2 = -(ru*i1 + iu*r1)

            ! .33219905d4 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .33219905d4
            i1 = -(ru*i2 + iu*r2)

            ! .3618331d5 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .3618331d5
            i2 = -(ru*i1 + iu*r1)

            ! x = z*ans
            rx = a*r2 + v*i2
            ix = a*i2 - v*r2

            !
            ! .1841439d1 - u <- ans
            r1 = .1841439d1 - ru
            i1 = -iu

            ! .6157037d2 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .6157037d2
            i2 = -(ru*i1 + iu*r1)

            ! .3642191d3 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .3642191d3
            i1 = -(ru*i2 + iu*r2)

            ! .2186181d4 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .2186181d4
            i2 = -(ru*i1 + iu*r1)

            ! .9022228d4 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .9022228d4
            i1 = -(ru*i2 + iu*r2)

            ! .2432284d5 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .2432284d5
            i2 = -(ru*i1 + iu*r1)

            ! conj(.320666d5 - u*ans) <- ans
            r1 = -(ru*r2 - iu*i2) + .320666d5
            i1 =  (ru*i2 + iu*r2)

            ! denominator
            de = 1d0/(r1*r1 + i1*i1)

            !
            ! x*ans <- ans
            r2 = (r1*rx - i1*ix)
            i2 = (r1*ix + i1*rx)

            ! Output exp(u) - x/y
            t = exp(dcmplx(ru,iu)) - dcmplx(r2,i2)*de

          end if ! Third or Fourth
        end if ! Second or Third/Fourth
      end if ! First or Second/Third/Fourth

      end subroutine voigt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the redistribution function in the laboratory frame
      !! with non-coherent lower term and infinitely sharp lower
      !! levels\n
      !!  omega1(double): Input frequency\n
      !!   omega(double): Output frequency\n
      !!      Dw(double): Output transition Doppler width\n
      !!     Dw1(double): Input transition Doppler width\n
      !!      el(double): Level l energy\n
      !!     el1(double): Level l' energy\n
      !!      eu(double): Level u energy\n
      !!     eu1(double): Level u' energy\n
      !!      ef(double): Level f energy\n
      !!      al(double): Term l inverse lifetime\n
      !!      au(double): Term u inverse lifetime\n
      !!      af(double): Term f inverse lifetime\n
      !!     aul(double): Input transition elastic width\n
      !!     auf(double): Output transition elastic width\n
      !!      C1(double): Cosine of scattering angle\n
      !!      S1(double): Sine of scattering angle\n
      !!  stype(integer): Type of scattering (geometry wise)
      function Wfunc(omega1,omega,Dw,Dw1,el,el1,eu,eu1,ef,al,au,af, &
                     aul,auf,C1,S1,stype)

      ! I/O

      integer, intent(in):: stype
      double precision, intent(in):: omega, omega1, Dw, Dw1
      double precision, intent(in):: el, el1, eu, eu1, ef, C1, S1
      double precision, intent(in):: al, au, af, aul, auf

      ! Local

      double precision:: iDw0,xil,xif,kkp,kkm,norm,inorm
      double precision:: vul,vu1l,vul1,vu1l1,wuf,wu1f
      double precision:: atl,atf,earg,earg1,efact,efact1

      complex(kind=8):: prof1,prof2,prof3,prof4,Wfunc


      ! Inverse composed Doppler width
      iDw0 = 1d0/(sqrt(Dw*Dw + Dw1*Dw1 - 2d0*Dw*Dw1*C1))

      ! Doppler widths factors
      xil = Dw1*iDw0
      xif = Dw*iDw0

      ! Branching ratios for resonances
      kkp = .5d0*(1d0 + xif*xif - xil*xil)
      kkm = 1d0 - kkp

      ! Resonances in input frequency
      vul = (omega1 - (eu-el))*iDw0
      vu1l = (omega1 - (eu1-el))*iDw0
      vul1 = (omega1 - (eu-el1))*iDw0
      vu1l1 = (omega1 - (eu1-el1))*iDw0

      ! Resonances in output frequency
      wuf = (omega - (eu-ef))*iDw0
      wu1f = (omega - (eu1-ef))*iDw0

      ! Exponential arguments
      earg = abs(wuf - vul)
      earg1 = abs(wuf - vul1)

      ! If argument too big
      if (earg.gt.Wbiggauss) then

        ! Make zero
        efact = 0d0

      ! Normal argument
      else

        ! Square
        earg = earg*earg

        ! If not too small
        if (earg.gt.smallexp) then

          ! Calculate exponential
          efact = exp(-earg)

        ! Too small
        else

          ! Taylor series
          efact = 1d0 + earg*(0.5d0*earg - 1d0)

        end if ! Normal or too small argument
      end if ! Too big argument

      ! If argument too big
      if (earg1.gt.Wbiggauss) then

        ! Make zero
        efact1 = 0d0

      ! Normal argument
      else

        ! Square
        earg1 = earg1*earg1

        ! If not too small
        if (earg1.gt.smallexp) then

          ! Calculate exponential
          efact1 = exp(-earg1)

        ! Too small
        else

          ! Taylor series
          efact1 = 1d0 + earg1*(0.5d0*earg1 - 1d0)

        end if ! Normal or too small argument
      end if ! Too big argument

      ! Get damping parameters
      atf = (au+af+auf)*iDw0
      atl = (au+al+aul)*iDw0

      ! General scattering
      if (stype.eq.0) then

        ! Normalization
        inorm = 1d0/(S1*xil*xif)

        ! Complete Damping parameter factors
        atf = atf*inorm
        atl = atl*inorm

        ! Compute Voigt profiles
        call voigt((kkp*vul+kkm*wuf)*inorm,atf,prof1)
        call voigt((kkp*vu1l+kkm*wu1f)*inorm,atl,prof2)
        call voigt((kkp*vul1+kkm*wuf)*inorm,atl,prof3)
        call voigt((kkp*vu1l1+kkm*wu1f)*inorm,atf,prof4)

        ! Complete normalization factor
        norm = PI*iDw0*iDw0*inorm

      ! Backward scattering
      else

        ! Get profiles
        prof1 = cImag/(kkp*vul   + kkm*wuf  + cImag*atf)
        prof2 = cImag/(kkp*vu1l  + kkm*wu1f + cImag*atl)
        prof3 = cImag/(kkp*vul1  + kkm*wuf  + cImag*atl)
        prof4 = cImag/(kkp*vu1l1 + kkm*wu1f + cImag*atf)

        ! Complete normalization factor
        norm =  sqrt(PI)*xil*xif

      end if ! Type of scattering

      ! Get redistribution function
      Wfunc = ( efact*(prof1 + conjg(prof2)) + &
               efact1*(prof3 + conjg(prof4)))*norm

      end function Wfunc

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the real Voigt profile\n
      !!   v(double): Normalized frequency\n
      !!   a(double): Normalized damping parameter\n
      !!  ou(double): Voigt profile
      subroutine voigtI(v,a,ou)

      ! I/O

      double precision, intent(in):: v, a
      double precision, intent(out):: ou

      ! Type of Voigt function
      select case (VOITY)

        ! Humlicek 1982 (accuracy to 10^-4 but with Faraday)
        case (0)

          ! Call relevant routine
          call voigtI_HC(v,a,ou)

        ! Armstrong 1967 (accurate to 6 figs and slow for a > 1.5)
        case (1)

          ! Call relevant routine
          call voigtI_A(v,a,ou)

        ! Hui et al. 1978 (1% accurate larger v, faster large a)
        case (2)

          ! Call relevant routine
          call voigtI_H(v,a,ou)

        ! Gaussian
        case (3)

          ! Call relevant routine
          call gaussianI(v,ou)

        ! Error
        case default

          ! Unknown
          ou = 0d0

      end select

      return

      end subroutine voigtI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the real Voigt profile with Armstrong (1967)
      !! algorithm\n
      !!   v(double): Normalized frequency\n
      !!   a(double): Normalized damping parameter\n
      !!  ou(double): Voigt profile
      subroutine voigtI_A(v,a,ou)

      ! I/O

      double precision, intent(in):: v, a
      double precision, intent(out):: ou

      ! Local

      double precision:: vi


      ! Get absolute value of normalized frequency
      if (v.lt.0d0) then
        vi = -v
      else
        vi = v
      end if

      ! First branch
      if ((a.lt.1d0.and.vi.lt.4d0).or.(a*(1d0+vi).lt.1.8d0)) then

        ! First branch
        ou = voigtI_AK1(vi,a)

      ! Second branch
      else if (a.lt.2.5d0.and.vi.lt.4d0) then

        ! Second branch
        ou = voigtI_AK2(vi,a)

      ! Thid branch
      else

        ! Third branch
        ou = voigtI_AK3(vi,a)

      end if ! Calculation branch

      return

      end subroutine voigtI_A

!#####################################################################
!#####################################################################
!#####################################################################

      !> Branch 1 of Armstrong (1967) algorithm\n
      !!   v(double): Normalized frequency\n
      !!   a(double): Normalized damping parameter
      double precision function voigtI_AK1(v,a)

      ! I/O

      double precision, intent(in):: v,a

      ! Local

      integer, parameter:: nn = 34

      double precision, dimension(nn), parameter:: cc = &
        (/ 0.1999999999972224d0, -0.1840000000029998d0, &
           0.1558399999965025d0, -0.1216640000043988d0, &
           0.0877081599940391d0, -0.0585141248086907d0, &
           0.0362157301623914d0, -0.0208497654398036d0, &
           0.0111960116346270d0, -0.56231896167109d-2, &
           0.26487634172265d-2, -0.11732670757704d-2, &
           0.4899519978088d-3, -0.1933630801528d-3, &
           0.722877446788d-4, -0.256555124979d-4, &
           0.86620736841d-5, -0.27876379719d-5, &
           0.8566873627d-6, -0.2518433784d-6, 0.709360221d-7, &
           -0.191732257d-7, 0.49801256d-8, -0.12447734d-8, &
           0.2997777d-9, -0.696450d-10, 0.156262d-10, -0.33897d-11, &
           0.7116d-12, -0.1447d-12, 0.285d-13, -0.55d-14, &
           0.10d-14, -0.2d-15 /)

      integer:: n

      double precision:: a2,v2,u1,d1,d2,dn,iv2,ff,an
      double precision:: q,g,coeff,b1,b2,bn,v1


      ! Square
      a2 = a*a
      v2 = v*v

      ! Too big argument
      if ((v2-a2).gt.70d0) then

        ! Make zero
        u1 = 0d0

      ! Not too big
      else

        ! Get exponential
        u1 = exp(a2 - v2)*cos(2d0*v*a)

      end if

      ! Wing
      if (v.gt.5d0) then

        ! Compute coefficients
        iv2 = 1d0/v2
        d1 = -iv2*(0.5d0 + iv2*(0.75d0 + iv2*(1.875d0 +  &
             iv2*(6.5625d0 + iv2*(29.53125d0 + &
             iv2*(1162.4218d0 + iv2*1055.7421d0))))))
        d2 = 0.5d0*(1d0 - d1)/v

      ! Core
      else

        ! Compute coefficients
        b1 = 0d0
        b2 = 0d0
        v1 = 0.2d0*v
        coeff = 4d0*v1*v1 - 2d0
        do n=nn,1,-1
          bn = coeff*b1 - b2 + cc(n)
          b2 = b1
          b1 = bn
        end do
        d2 = v1*(bn - b2)
        d1 = 1d0 - 2d0*v*d2

      end if ! Wing or core

      ! Initial ff
      ff = a*d1

      ! If damping large enough
      if (a.gt.1d-8) then

        ! Start fit
        q = 1d0
        an = a

        ! Series
        do n=2,50

          ! Update coefficients
          dn = (v*d1 + d2)*(-2d0/dble(n))
          d2 = d1
          d1 = dn

          ! Even elements
          if (mod(n,2).gt.0) then

            ! Update ff
            q = -q
            an = an*a2
            g = dn*an
            ff = ff + q*g

            ! If ratio small
            if (abs(g/ff).le.1d-8) then

              ! Finish calculation
              voigtI_AK1 = u1 - 1.12837917d0*ff
              return

            end if ! Ratio g/ff small
          end if ! Even elements

        end do ! Series

      end if ! Large enough damping

      ! Complete calculation
      voigtI_AK1 = u1 - 1.12837917d0*ff

      return

      end function voigtI_AK1

!#####################################################################
!#####################################################################
!#####################################################################

      !> Branch 2 of Armstrong (1967) algorithm\n
      !!  v(double): Normalized frequency\n
      !!  a(double): Normalized damping parameter
      double precision function voigtI_AK2(v,a)

      ! I/O

      double precision, intent(in):: v, a

      ! Local

      integer:: n
      double precision:: g,r,s,a2,ia

      integer, parameter:: nn = 10

      double precision, dimension(nn), parameter:: t = &
        (/ 0.2453407083d0, 0.7374737285d0, 1.2340762153d0, &
           1.7385377121d0, 2.2549740020d0, 2.7888060584d0, &
           3.3478545673d0, 3.9447640401d0, 4.6036824495d0, &
           5.3874808900d0 /)

      double precision, dimension(nn), parameter:: w = &
        (/ 4.6224366960d-01, 2.8667550536d-01, 1.0901720602d-01, &
           2.4810520887d-02, 3.2437733422d-03, 2.2833863601d-04, &
           7.8025564785d-06, 1.0860693707d-07, 4.3993409922d-10, &
           2.2293936455d-13 /)

      ! Square and inverse
      a2 = a*a
      ia = 1d0/a

      ! Initialize
      g = 0d0

      ! Series
      do n=1,nn

        ! Add contribution
        r = t(n) - v
        s = t(n) + v
        g = g + (4d0*t(n)*t(n) - 2d0)*(r*atan(r*ia) + &
            s*atan(s*ia) - 0.5d0*a*(log(a2 + r*r) + &
            log(a2 + s*s)))*w(n)

      end do ! Series

      ! Finish calculation
      voigtI_AK2 = g*IPI

      return

      end function voigtI_AK2

!#####################################################################
!#####################################################################
!#####################################################################

      !> Branch 3 of Armstrong (1967) algorithm\n
      !!  v(double): Normalized frequency\n
      !!  a(double): Normalized damping parameter
      double precision function voigtI_AK3(v,a)

      ! I/O

      double precision, intent(in):: v,a

      ! Local

      integer:: n

      double precision:: g,a2

      integer, parameter:: nn = 10

      double precision, dimension(nn), parameter:: t = &
        (/ 0.2453407083d0, 0.7374737285d0, 1.2340762153d0, &
           1.7385377121d0, 2.2549740020d0, 2.7888060584d0, &
           3.3478545673d0, 3.9447640401d0, 4.6036824495d0, &
           5.3874808900d0 /)

      double precision, dimension(nn), parameter:: w = &
        (/ 4.6224366960d-01, 2.8667550536d-01, 1.0901720602d-01, &
           2.4810520887d-02, 3.2437733422d-03, 2.2833863601d-04, &
           7.8025564785d-06, 1.0860693707d-07, 4.3993409922d-10, &
           2.2293936455d-13 /)

      ! Square
      a2 = a*a

      ! Initialize
      g = 0d0

      ! Series
      do n=1,nn

        ! Add contribution
        g = g + (1d0/((v-t(n))*(v-t(n)) + a2) + &
                 1d0/((v+t(n))*(v+t(n)) + a2))*w(n)

      end do ! Series

      ! Complete calculation
      voigtI_AK3 = a*g*IPI

      return

      end function voigtI_AK3

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the real Voigt profile with Hui et al. (1978)
      !! algorithm\n
      !!   v(double): Normalized frequency\n
      !!   a(double): Normalized damping parameter\n
      !!  ou(double): Voigt profile
      subroutine voigtI_H(v,a,ou)

      ! I/O

      double precision, intent(in):: v,a
      double precision, intent(out):: ou

      ! Local

      complex(kind=8):: z,W


      ! Make complex
      z = dcmplx(a, -v)

      ! Calculate
      W = (122.607931777104326d0 + z*(214.382388694706425d0 + &
           z*(181.928533092181549d0 + z*(93.155580458138441d0 + &
           z*(30.180142196210589d0 + z*(5.912626209773153d0 + &
           z*0.564189583562615d0))))))/ &
           (122.607931773875350d0 + z*(352.730625110963558d0 + &
           z*(457.334478783897737d0 + z*(348.703917719495792d0 + &
           z*(170.354001821091472d0 + z*(53.992906912940207d0 + &
           z*(10.479857114260399d0 + z)))))))

      ! Return real
      ou = dble(W)

      return

      end subroutine voigtI_H

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the real Voigt profiles with Humlicek (1982)
      !! algorithm. Real algebra version\n
      !!   v(double): Normalized frequency\n
      !!   a(double): Normalized damping parameter\n
      !!  ou(double): Voigt profile
      subroutine voigtI_HC(v,a,ou)

      ! I/O

      double precision, intent(in):: v,a
      double precision, intent(out):: ou

      ! Local

      real:: sa,sv
      real:: s,d

      double precision:: r1,i1,r2,i2,ru,iu,rx,ix,de

      ! Get single precision arguments
      sa = sngl(a)
      sv = sngl(v)

      ! Compute quantities to decide branching
      s = abs(sv)+sa
      d = .195e0*abs(sv)-.176e0

      ! First branch
      if (s.ge..15e2) then

        ! Output
        ! t = .5641896d0*z/(.5d0+z*z)

        ! (0.5 + z^2)* <- ans
        r1 = 0.5d0 + a*a - v*v
        i1 = 2d0*a*v

        ! Denominator
        de = .5641896d0/(r1*r1 + i1*i1)

        ! Output ans*z*de
        ou = (a*r1 + v*i1)*de

      ! Rest of branches
      else

        ! Second branch
        if (s.ge..55e1) then

          ! Output
          ! u = z*z
          ! t = z*(.1410474d1 + .5641896d0*u)/(.75d0 + u*(.3d1 + u))

          ! u = z*z
          ru = a*a - v*v
          iu = -2d0*a*v

          ! 0.75 + u(3 + u) <- ans
          r1 =  0.75d0 + ru*(ru + 3d0) - iu*iu
          i1 =  -iu*(3d0 + 2d0*ru)

          ! Denominator
          de = 1d0/(r1*r1 + i1*i1)

          ! z*conj(ans) <- ans
          r2 = a*r1 + v*i1
          i2 = a*i1 - v*r1

          ! .1410474d1 + .5641896d0*u <- bns
          ru = ru*.5641896d0 + .1410474d1
          iu = iu*.5641896d0

          ! Output ans*bns*de
          ou = (r2*ru - i2*iu)*de

        ! Rest of branches
        else

          ! Third branch
          if (sa.ge.d) then

            ! Output
            ! nt = .164955d2 + z*(.2020933d2 + z*(.1196482d2 + &
            !      z*(.3778987d1 + .5642236d0*z)))
            ! dt = .164955d2 + z*(.3882363d2 + z*(.3927121d2 + &
            !      z*(.2169274d2 + z*(.6699398d1 + z))))
            ! t = nt/dt

            ! .3778987d1 + .5642236d0*z <- ans
            r1 = .3778987d1 + .5642236d0*a
            i1 = -.5642236d0*v

            ! .1196482d2 + z*ans <- ans
            r2 = a*r1 + v*i1 + .1196482d2
            i2 = a*i1 - v*r1

            ! .2020933d2 + z*ans <- ans
            r1 = a*r2 + v*i2 + .2020933d2
            i1 = a*i2 - v*r2

            ! u = .164955d2 + z*ans
            ru = a*r1 + v*i1 + .164955d2
            iu = a*i1 - v*r1

            !
            ! .6699398d1 + z <- ans
            r1 = .6699398d1 + a
            i1 = -v

            ! .2169274d2 + z*ans <- ans

            r2 = a*r1 + v*i1 + .2169274d2
            i2 = a*i1 - v*r1

            ! .3927121d2 + z*ans <- ans
            r1 = a*r2 + v*i2 + .3927121d2
            i1 = a*i2 - v*r2

            ! .3882363d2 + z*ans <- ans
            r2 = a*r1 + v*i1 + .3882363d2
            i2 = a*i1 - v*r1

            ! conj(.164955d2 + z*ans) <- ans
            r1 =   a*r2 + v*i2 + .164955d2
            i1 = -(a*i2 - v*r2)

            ! Denominator
            de = 1d0/(r1*r1 + i1*i1)

            ! Output ans*u*de
            ou = (ru*r1 - iu*i1)*de

          ! Fourth brach
          else

            ! Output
            ! u = z*z
            ! x = z*(.3618331d5 - u*(.33219905d4 - u*(.1540787d4 - &
            !     u*(.2190313d3 - u*(.3576683d2 - u*(.1320522d1 - &
            !     .56419d0*u))))))
            ! y = .320666d5 - u*(.2432284d5 - u*(.9022228d4 - &
            !     u*(.2186181d4 - u*(.3642191d3 - u*(.6157037d2 - &
            !     u*(.1841439d1 - u))))))
            ! t = exp(u) - x/y

            ! u = z*z
            ru =  a*a - v*v
            iu = -2d0*a*v

            !
            ! .1320522d1 - .56419d0*u <- ans
            r1 =  .1320522d1 - .56419d0*ru
            i1 = -.56419d0*iu

            ! .3576683d2 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .3576683d2
            i2 = -(ru*i1 + iu*r1)

            ! .2190313d3 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .2190313d3
            i1 = -(ru*i2 + iu*r2)

            ! .1540787d4 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .1540787d4
            i2 = -(ru*i1 + iu*r1)

            ! .33219905d4 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .33219905d4
            i1 = -(ru*i2 + iu*r2)

            ! .3618331d5 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .3618331d5
            i2 = -(ru*i1 + iu*r1)

            ! x = z*ans
            rx = a*r2 + v*i2
            ix = a*i2 - v*r2

            !
            ! .1841439d1 - u <- ans
            r1 = .1841439d1 - ru
            i1 = -iu

            ! .6157037d2 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .6157037d2
            i2 = -(ru*i1 + iu*r1)

            ! .3642191d3 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .3642191d3
            i1 = -(ru*i2 + iu*r2)

            ! .2186181d4 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .2186181d4
            i2 = -(ru*i1 + iu*r1)

            ! .9022228d4 - u*ans <- ans
            r1 = -(ru*r2 - iu*i2) + .9022228d4
            i1 = -(ru*i2 + iu*r2)

            ! .2432284d5 - u*ans <- ans
            r2 = -(ru*r1 - iu*i1) + .2432284d5
            i2 = -(ru*i1 + iu*r1)

            ! conj(.320666d5 - u*ans) <- ans
            r1 = -(ru*r2 - iu*i2) + .320666d5
            i1 =  (ru*i2 + iu*r2)

            ! denominator
            de = 1d0/(r1*r1 + i1*i1)

            !
            ! x*ans <- ans
            r2 = (r1*rx - i1*ix)

            ! Output exp(u) - x/y
            ou = dble(exp(dcmplx(ru,iu))) - r2*de

          end if ! Third or fourth branch
        end if ! Second or Third/fourth branches
      end if ! First or Second/third/fourth branches

      end subroutine voigtI_HC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate a Gaussian profile\n
      !!   v(double): Normalized frequency\n
      !!  ou(double): Voigt profile
      subroutine gaussianI(v,ou)

      ! I/O

      double precision, intent(in):: v
      double precision, intent(out):: ou

      ! Local

      double precision:: v2


      ! True argument
      v2 = v*v

      ! Control overflow
      if (v2.gt.bigexp) then

        ! Zero due to big argument
        ou = 0d0

      ! Normal
      else if (v2.gt.smallexp) then

        ! Get exponential
        ou = exp(-v2)

      ! Control underflow
      else

        ! Second order Taylor
        ou = 1d0 - v2 + 0.5d0*v2*v2

      end if

      end subroutine gaussianI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the real part of the redistribution function in the
      !! laboratory frame with non-coherent lower term and infinitely
      !! sharp lower levels\n
      !!  omega1(double): Input frequency\n
      !!   omega(double): Output frequency\n
      !!      Dw(double): Output transition Doppler width\n
      !!     Dw1(double): Input transition Doppler width\n
      !!     atl(double): Input transition width\n
      !!     atf(double): Output transition width\n
      !!      C1(double): Cosine of scattering angle\n
      !!      S1(double): Sine of scattering angle\n
      !!  stype(integer): Type of scattering (geometry wise)
      function WfuncI(omega1,omega,Dw,Dw1,atl,atf,C1,S1,stype)

      ! I/O

      integer, intent(in):: stype
      double precision, intent(in):: omega,omega1,Dw,Dw1
      double precision, intent(in):: atf,atl,C1,S1

      ! Local

      double precision:: iDw0,xil,xif,kkp,kkm,norm,inorm
      double precision:: x,a1,a2,vul,wuf,earg,efact
      double precision:: WfuncI,prof1,prof2


      ! Inverse composed Doppler width
      iDw0 = 1d0/(sqrt(Dw*Dw + Dw1*Dw1 - 2d0*Dw*Dw1*C1))

      ! Doppler widths factors
      xil = Dw1*iDw0
      xif = Dw*iDw0

      ! Branching ratios for resonances
      kkp = .5d0*(1d0 + xif*xif - xil*xil)
      kkm = 1d0 - kkp

      ! Resonances in input frequency
      vul = omega1*iDw0
      wuf = omega*iDw0

      ! Exponential argument
      earg = abs(wuf - vul)

      ! If argument too big
      if (earg.gt.Wbiggauss) then

        ! Make zero
        efact = 0d0

      ! Normal argument
      else

        ! Square
        earg = earg*earg

        ! If not too small
        if (earg.gt.smallexp) then

          ! Calculate exponential
          efact = exp(-earg)

        ! Too small
        else

          ! Taylor series
          efact = 1d0 + earg*(0.5d0*earg - 1d0)

        end if ! Normal or too small argument
      end if ! Too big argument

      ! Get resonance and damping parameters
      x = kkp*vul + kkm*wuf
      a1 = atf*iDw0
      a2 = atl*iDw0

      ! General scattering
      if (stype.eq.0) then

        ! Normalization factor
        inorm = 1d0/(S1*xil*xif)

        ! Normalization
        x = x*inorm
        a1 = a1*inorm
        a2 = a2*inorm

        ! Get Voigt profiles
        call voigtI(x,a1,prof1)
        call voigtI(x,a2,prof2)

        ! Complete normalization factor
        norm = PI*iDw0*iDw0*inorm

      ! 3-term or 2-term backward scattering
      else

        ! Get profiles
        prof1 = dble(cImag/(x + cImag*a1))
        prof2 = dble(cImag/(x + cImag*a2))

        ! Complete normalization factor
        norm =  sqrt(PI)*xil*xif

      endif

      ! Get redistribution function
      WfuncI = 2d0*efact*(prof1 + prof2)*norm

      end function WfuncI

!#####################################################################
!#####################################################################
!#####################################################################

      end module profile_mod
