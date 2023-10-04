      !> Gauss quadrature
      module gaussaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     02/13/2023 V3.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     02/13/2023:    V3.0.1 - Added fullgauss and fullazimuth
!                             subroutines (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     02/20/2019:    V1.1.0 - Use specific TINY (TdPA)
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
!  fullgauss:
!    Set-up the polar quadrature with split hemispheres
!
!  fullaxial:
!    Set-up axial quadrature
!
!  gaussaux:
!    Auxiliar function to define gauss quadrature
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

      !> Returns gaussian quadrature\n
      !!        nth(integer): Nodes\n
      !!     V_mu(double(:)): Output node values in mu\n
      !!     W_mu(double(:)): Output weights
      subroutine fullgauss(nTh,V_mu,W_mu)

      ! I/O
      integer, intent(in):: nTh
      double precision, dimension(:), intent(out):: V_mu,W_mu

      ! Local
      integer:: ii,nang,da,n
      double precision:: a,b,ren
      double precision, dimension(nTh/2):: v_mu_l, w_mu_l, xx, ww

      ! Initialize indexes
      n = nTh/2
      nang = (n + 1)/2
      da = mod(n,2) - 1
      a = -1d0
      b =  1d0

      ! Calculate nodes and weights
      call gaussaux(a, b, xx, ww, n)

      ! Reorder the nodes and weights
      ! Second half
      do ii=n,nang-da, -1

        v_mu_l(n - ii + 1) = xx(ii)
        w_mu_l(n - ii + 1) = ww(ii)

      end do

      ! First half
      do ii=1,nang-da-1

        v_mu_l(nang + ii) = xx(nang - da - ii)
        w_mu_l(nang + ii) = ww(nang - da - ii)

      end do

      ! Initialize variable
      ren = 0d0

      ! Integrate
      do ii=1,n

        ren = ren + w_mu_l(ii)

      end do

      ! Normalize
      do ii=1,n

        w_mu_l(ii) = w_mu_l(ii)/ren

      end do

      ! Redefine nodes and weights with a variable change to split
      ! the integral from -1 to 1 into two integrals
      do ii=1,n

        V_mu(ii) = .5d0*(v_mu_l(n + 1 - ii) - 1d0)
        W_mu(ii) = .5d0*w_mu_l(n + 1 - ii)

        V_mu(2*n + 1 - ii) = .5d0*(v_mu_l(ii) + 1d0)
        W_mu(2*n + 1 - ii) = .5d0*w_mu_l(ii)

      end do

      end subroutine fullgauss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Returns azimuthal quadrature for RT problem
      !!        nph(integer): Nodes\n
      !!     V_phi(double(:)): Output node values in angle\n
      !!     V_mux(double(:)): Output node value in cosine\n
      !!     V_muy(double(:)): Output node sign of sine\n
      !!     W_mux(double(:)): Output weights
      subroutine fullazimuth(nPh,V_phi,V_mux,V_muy,W_mux)

      ! I/O
      integer, intent(in):: nPh
      double precision, dimension(:), intent(out):: V_phi,V_mux
      double precision, dimension(:), intent(out):: V_muy,W_mux

      ! Local
      integer:: ii
      double precision:: ren

      ! Define phi, cos(phi) and weight with trapezoidal rule
      do ii=1,nPh

        V_phi(ii) = 2d0*PI*(dble(ii) - .5D0)/dble(nPh)
        V_mux(ii) = cos(V_phi(ii))
        W_mux(ii) = 1d0/dble(nPh)

      end do

      ! Define the sign of sin(phi)
      ! First half positive
      do ii=1,nPh/2

        V_muy(ii) = 1d0

      end do

      ! Second half positive
      do ii=nPh/2+1,nPh

        V_muy(ii) = -1d0

      end do

      ! Reset variable
      ren = 0d0

      ! Integrate in azimuth
      do ii=1,nPh

        ren = ren + W_mux(ii)

      end do

      ! Normalize the weights
      do ii=1,nPh

        W_mux(ii) = W_mux(ii)/ren

      enddo

      end subroutine fullazimuth

!#####################################################################
!#####################################################################
!#####################################################################

      !> Returns gaussian nodes and weights within a given interval\n
      !!       x0(dfloat): Initial point\n
      !!       x1(dfloat): Final point\n
      !!    xx(dfloat(:)): Node coordinates\n
      !!    ww(dfloat(:)): Node weights\n
      !!      n1(integer): Number of nodes
      subroutine gaussaux(x0,x1,xx,ww,n1)

      ! I/O

      integer, intent(in):: n1
      double precision, intent(in):: x0, x1
      double precision, dimension(:), intent(out):: xx, ww

      ! Local

      integer:: i1, i2, n2

      double precision:: d1, d2, d3, d4, xm, xp, mu0, mu1

      n2 = (n1+1)/2

      xp = .5d0*(x1+x0)
      xm = .5d0*(x1-x0)

      do i1=1,n2

        mu0 = cos(PI*(dble(i1) - .25d0)/(dble(n1)+.5d0))

        do while (.True.)

          d1 = 1d0
          d2 = 0d0

          do i2=1,n1

            d3 = d2
            d2 = d1
            d1 = ((2d0*dble(i2)-1d0)*mu0*d2 - (dble(i2)-1d0)*d3)/i2

          end do

          d4 = dble(n1)*(mu0*d1-d2)/(mu0*mu0-1d0)
          mu1 = mu0
          mu0 = mu1 - d1/d4

          if(abs(mu1 - mu0).le.TINYM)exit

        end do

        xx(i1) = xp - xm*mu0
        ww(i1) = 2d0*xm/((1d0-mu0*mu0)*d4*d4)

        xx(n1+1-i1) = xp + xm*mu0
        ww(n1+1-i1) = ww(i1)

      end do

      end subroutine gaussaux

!#####################################################################
!#####################################################################
!#####################################################################

      end module gaussaux_mod
