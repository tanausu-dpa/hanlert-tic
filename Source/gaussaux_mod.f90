      !> Gauss quadrature
      module gaussaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     18/04/2017
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
!  fullgauss
!    Set-up a gaussian quadrature as the combination of two gaussian
!  quadratures between [-1,0] and [0,1]
!
!  fullaxial
!    Set-up a trapezoidal quadrature in [0,2pi)
!
!  gaussaux
!    Set-up a gaussian quadrature in an arbitrary range
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use parameters_mod , only : TINYM , PI

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up a gaussian quadrature as the combination of two
      !! gaussian quadratures between [-1,0] and [0,1]\n
      !!     nth(integer): Number of nodes\n
      !!  V_mu(double(:)): Node coordinate values in mu\n
      !!  W_mu(double(:)): Node weight values 
      subroutine fullgauss(nTh,V_mu,W_mu)

      ! I/O

      integer, intent(in):: nTh
      double precision, dimension(:), intent(out):: V_mu,W_mu

      ! Local

      integer:: ii,nang,da,n

      double precision:: a,b,ren
      double precision, dimension(nTh/2):: v_mu_l, w_mu_l, xx, ww


      ! Initialize dimensions, indexes, and auxiliar quantities
      n = nTh/2
      nang = (n + 1)/2
      da = mod(n,2) - 1
      a = -1d0
      b =  1d0

      ! Calculate gaussian quadrature nodes and weights in [-1,1]
      call gaussaux(a, b, xx, ww, n)

      ! Second half
      do ii=n,nang-da, -1

        ! Reorder the nodes and weights
        v_mu_l(n - ii + 1) = xx(ii)
        w_mu_l(n - ii + 1) = ww(ii)

      end do

      ! First half
      do ii=1,nang-da-1

        ! Reorder the nodes and weights
        v_mu_l(nang + ii) = xx(nang - da - ii)
        w_mu_l(nang + ii) = ww(nang - da - ii)

      end do

      ! Initialize normalization variable
      ren = 0d0

      ! Integrate
      do ii=1,n

        ! Add contribution
        ren = ren + w_mu_l(ii)

      end do

      ! For each node
      do ii=1,n

        ! Normalize weights
        w_mu_l(ii) = w_mu_l(ii)/ren

      end do

      ! Redefine nodes and weights with a variable change to split
      ! the integral in [-1,1] into two integrals in [-1,0] and
      ! [0,1]
      do ii=1,n

        ! Node and weight in [-1,0]
        V_mu(ii) = .5d0*(v_mu_l(n + 1 - ii) - 1d0)
        W_mu(ii) = .5d0*w_mu_l(n + 1 - ii)

        ! Node and weight in [0,1]
        V_mu(2*n + 1 - ii) = .5d0*(v_mu_l(ii) + 1d0)
        W_mu(2*n + 1 - ii) = .5d0*w_mu_l(ii)

      end do

      end subroutine fullgauss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up a trapezoidal quadrature in [0,2pi)\n
      !!      nph(integer): Number of nodes\n
      !!  V_phi(double(:)): Node coordinate values in angle\n
      !!  V_mux(double(:)): Node coordinate values in cosine\n
      !!  V_muy(double(:)): Node sign of the sines\n
      !!  W_mux(double(:)): Node weight values
      subroutine fullazimuth(nPh,V_phi,V_mux,V_muy,W_mux)

      ! I/O

      integer, intent(in):: nPh
      double precision, dimension(:), intent(out):: V_phi,V_mux
      double precision, dimension(:), intent(out):: V_muy,W_mux

      ! Local

      integer:: ii

      double precision:: ren


      ! Fpr each node
      do ii=1,nPh

        ! Get angle
        V_phi(ii) = 2d0*PI*(dble(ii) - .5D0)/dble(nPh)

        ! Cosine
        V_mux(ii) = cos(V_phi(ii))

        ! And weight
        W_mux(ii) = 1d0/dble(nPh)

      end do ! nodes

      !
      ! Define the sign of sin(phi)
      !

      ! First half positive
      do ii=1,nPh/2

        V_muy(ii) = 1d0

      end do

      ! Second half positive
      do ii=nPh/2+1,nPh

        V_muy(ii) = -1d0

      end do

      !
      ! Normalize
      !

      ! Reset variable
      ren = 0d0

      ! Integrate in azimuth
      do ii=1,nPh

        ! Add contribution
        ren = ren + W_mux(ii)

      end do

      ! For each node
      do ii=1,nPh

        ! Normalize the weights
        W_mux(ii) = W_mux(ii)/ren

      enddo

      end subroutine fullazimuth

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set-up a gaussian quadrature in an arbitrary range\n
      !!     x0(double): Initial coordinate\n
      !!     x1(double): Final coordinate\n
      !!  xx(double(:)): Node coordinate values\n
      !!  ww(double(:)): Node weight values\n
      !!    n1(integer): Number of nodes
      subroutine gaussaux(x0,x1,xx,ww,n1)

      ! I/O

      integer, intent(in):: n1
      double precision, intent(in):: x0, x1
      double precision, dimension(:), intent(out):: xx, ww

      ! Local

      integer:: i1, i2, n2

      double precision:: d1, d2, d3, d4, xm, xp, mu0, mu1


      ! Get upper limit of summation
      n2 = (n1+1)/2

      ! Get mid point and half size
      xp = .5d0*(x1+x0)
      xm = .5d0*(x1-x0)

      ! For every node (half)
      do i1=1,n2

        ! Initial mu value
        mu0 = cos(PI*(dble(i1) - .25d0)/(dble(n1)+.5d0))

        ! Iterate until convergence
        do while (.True.)

          ! Initialize
          d1 = 1d0
          d2 = 0d0

          ! Inner summation
          do i2=1,n1

            ! Get factor numerator
            d3 = d2
            d2 = d1
            d1 = ((2d0*dble(i2)-1d0)*mu0*d2 - (dble(i2)-1d0)*d3)/i2

          end do ! Inner summation

          ! Get factor denomination
          d4 = dble(n1)*(mu0*d1-d2)/(mu0*mu0-1d0)

          ! Save old value and get new one
          mu1 = mu0
          mu0 = mu1 - d1/d4

          ! If converged, exit
          if (abs(mu1 - mu0).le.TINYM) exit

        end do ! Iteration

        ! Get node value and weight
        xx(i1) = xp - xm*mu0
        ww(i1) = 2d0*xm/((1d0-mu0*mu0)*d4*d4)

        ! Get reflected node value and weight
        xx(n1+1-i1) = xp + xm*mu0
        ww(n1+1-i1) = ww(i1)

      end do ! Outer nodes (half)

      end subroutine gaussaux

!#####################################################################
!#####################################################################
!#####################################################################

      end module gaussaux_mod
