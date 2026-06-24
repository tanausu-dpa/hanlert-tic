      !> Background opacity contributions
      module backgroundaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     23/02/2017
!  Last version:
!     15/06/2026 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/06/2026:    V4.0.1 - Bugfix: When calling rayleigh for an
!                             active atom from chi_freq, the input
!                             index was undefined (TdPA)
!                           - Improvements in rayleigh function when
!                             the call is for an active atom (TdPA)
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
!    This module calculates the continuum quantities
!
!  Thomson
!    Computes Thomson scattering coefficient
!
!  Hminus_bf
!    Computes H- bound-free absorptivity and emissivity
!
!  Hminus_ff
!    Computes H- free-free absorptivity
!
!  HI_bf
!    Computes H bound-free absorptivity and emissivity
!
!  gHI_bf
!    Computes Hydrogen bound-free Gaunt factor
!
!  HI_ff
!    Computes H free-free absorptivity
!
!  gHI_ff
!    Computes Hydrogen free-free Gaunt factor
!
!  back_bf
!    Computes atomic bound-free absorptivity and emissivity
!
!  back_bb
!    Computes atomic bound-bound absorptivity and emissivity
!
!  backH_bb
!    Computes hydrogen bound-bound absorptivity and emissivity with
!  the hard-coded model
!
!  rayleigh
!    Computes Rayleigh scattering coefficient for hydrogen
!
!  rayleigh_H
!    Computes Rayleigh scattering coefficient for H2
!
!  OH_bf
!    Computes OH bound-free absorptivity and emissivity
!
!  CH_bf
!    Computes CH bound-free absorptivity and emissivity
!
!  H2m_ff
!    Computes H2- free-free absorptivity
!
!  HHp_ff
!    Computes H + p+ free-free absorptivity
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use funnj_mod
      use inter_mod
      use math_mod
      use parameters_mod , ONLY : c2 , convF , ergtoev , c , eps0 , &
                                  me , PI , pi4eps0 , qel , kb , &
                                  ryd , fktoJ , hplanck , ktoev
      use profile_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Thomson scattering coefficient\n
      !!  ne(double(:)): Electron number density\n
      !!   iz0(integer): First height index to consider\n
      !!   iz1(integer): Last height index to consider\n
      !! sig(double(:)): Thomson scattering coefficient
      subroutine Thomson(ne,iz0,iz1,sig)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, dimension(iz0:iz1), intent(in):: ne
      double precision, dimension(iz0:iz1), intent(out):: sig

      ! Local

      integer:: iz

      double precision:: tsig


      ! Thomson scattering cross section in cm^2
      tsig = 6.6524587158d-25

      ! For each height, multiply by electron density
      do iz=iz0,iz1
        sig(iz) = ne(iz)*tsig
      end do

      end subroutine Thomson

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes H- bound-free absorptivity and emissivity\n
      !!    freq(double): Frequency\n
      !!  nhm(double(:)): H- number density\n
      !!    T(double(:)): Temperature\n
      !!    iz0(integer): First height index to consider\n
      !!    iz1(integer): Last height index to consider\n
      !!  eta(double(:)): Absorptivity\n
      !!  eps(double(:)): Emissivity
      subroutine Hminus_bf(freq,nhm,T,iz0,iz1,eta,eps)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: nhm, T
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer, parameter:: nin = 34
      integer:: iz

      double precision:: lambda, alpha, exu, pE, arg
      double precision, dimension(nin):: ba,ca,da,x,y


      !
      ! Geltman (1962) table 3
      !

      ! Wavelength [nm]
      x = (/ 0d0, 5d1, 1d2, 1.5d2, 2d2, 2.5d2, 3d2, 3.5d2, 4d2, &
             4.5d2, 5d2, 5.5d2, 6d2, 6.5d2, 7d2, 7.5d2, 8d2, 8.5d2, &
             9d2, 9.5d2, 1d3, 1.050d3, 1.1d3, 1.15d3, 1.2d3, 1.25d3, &
             1.3d3, 1.35d3, 1.4d3, 1.45d3, 1.5d3, 1.55d3, 1.6d3, &
             1.6419d3 /)

      ! Cross section [1d19 m^2]
      y = (/ 0d0, 1.5d-1, 3.3d-1, 5.7d-1, 8.5d-1, 1.17d0, 1.52d0, &
             1.89d0, 2.23d0, 2.55d0, 2.84d0, 3.11d0, 3.35d0, &
             3.56d0, 3.71d0, 3.83d0, 3.92d0, 3.95d0, 3.93d0, 3.85d0, &
             3.73d0, 3.58d0, 3.38d0, 3.14d0, 2.85d0, 2.54d0, 2.20d0, &
             1.83d0, 1.46d0, 1.06d0, 7.1d-1, 4d-1, 1.7d-1, 0d0 /)

      ! Get the corresponding wavelength in nm
      lambda = 1d2/freq

      ! If out of limits
      if (lambda.le.x(1).or.lambda.ge.x(nin)) then

        ! No contribution
        eta = 0d0
        eps = 0d0

      ! If within limits
      else

        ! Frequency quantities
        arg = c2*freq*1d4
        pE = convF*2d21*c*freq*freq*freq

        ! Interpolate cross section with cubic splines
        call spline(x,y,ba,ca,da,nin)
        alpha = ispline(lambda,x,y,ba,ca,da,nin)*1d-17

        ! For each height to consider
        do iz=iz0,iz1

          ! Get inverse exponential
          exu = arg/T(iz)
          exu = diexp(exu)

          ! Compute absorptivity and emissivity
          eta(iz) = nhm(iz)*(1d0 - exu)*alpha
          eps(iz) = nhm(iz)*pE*exu*alpha

        end do ! Heights

      end if ! Within wavelength limits

      end subroutine Hminus_bf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes H- free-free absorptivity\n
      !!    freq(double): Frequency\n
      !!   nh(double(:)): H number density\n
      !!   ne(double(:)): Electron density\n
      !!    T(double(:)): Temperature\n
      !!    iz0(integer): First height index to consider\n
      !!    iz1(integer): Last height index to consider\n
      !!  eta(double(:)): Absorptivity
      subroutine Hminus_ff(freq,nh,ne,T,iz0,iz1,eta)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: nh, ne, T
      double precision, dimension(iz0:iz1), intent(out):: eta

      ! Local

      integer, parameter:: n1 = 16, n2 = 17, n3 = 6
      integer:: ii,iz

      double precision:: lambda, ilambda, tev, cc, C0
      double precision, dimension(n1):: x1
      double precision, dimension(n2):: x2
      double precision, dimension(n1,n2):: y
      double precision, dimension(n3) :: ya, yb, yc, yd, ye, yf, yg


      !
      ! Stilley & Callaway (1970) table 1
      !

      ! Theta [eV]
      x1 = (/ 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, 1d0, 1.1d0, 1.2d0, &
              1.3d0, 1.4d0, 1.5d0, 1.6d0, 1.7d0, 1.8d0, 1.9d0, &
              2.0d0 /)

      ! wavelength [nm]
      x2 = (/ 0d0, 303.8d0, 455.6d0, 506.3d0, 569.5d0, 650.9d0, &
              759.4d0, 911.3d0, 1013d0, 1139d0, 1302d0, 1519d0, &
              1823d0, 2278d0, 3038d0, 4556d0, 9113d0 /)

      ! cross section [1d29 m^2]
      y(:,1) = 0d0
      y(:,2) = (/ 3.44d-2, 4.18d-2, 4.91d-2, 5.65d-2, 6.39d-2, &
                  7.13d-2, 7.87d-2, 8.62d-2, 9.36d-2, 1.01d-1, &
                  1.08d-1, 1.16d-1, 1.23d-1, 1.30d-1, 1.38d-1, &
                  1.45d-1 /)
      y(:,3) = (/ 7.80d-2, 9.41d-2, 1.10d-1, 1.25d-1, 1.40d-1, &
                  1.56d-1, 1.71d-1, 1.86d-1, 2.01d-1, 2.16d-1, &
                  2.31d-1, 2.45d-1, 2.60d-1, 2.75d-1, 2.89d-1, &
                  3.03d-1 /)
      y(:,4) = (/ 9.59d-2, 1.16d-1, 1.35d-1, 1.53d-1, 1.72d-1, &
                  1.90d-1, 2.08d-1, 2.25d-1, 2.43d-1, 2.61d-1, &
                  2.78d-1, 2.96d-1, 3.13d-1, 3.30d-1, 3.47d-1, &
                  3.64d-1 /)
      y(:,5) = (/ 1.21d-1, 1.45d-1, 1.69d-1, 1.92d-1, 2.14d-1, &
                  2.36d-1, 2.58d-1, 2.80d-1, 3.01d-1, 3.22d-1, &
                  3.43d-1, 3.64d-1, 3.85d-1, 4.06d-1, 4.26d-1, &
                  4.46d-1 /)
      y(:,6) = (/ 1.56d-1, 1.88d-1, 2.18d-1, 2.47d-1, 2.76d-1, &
                  3.03d-1, 3.31d-1, 3.57d-1, 3.84d-1, 4.10d-1, &
                  4.36d-1, 4.62d-1, 4.87d-1, 5.12d-1, 5.37d-1, &
                  5.62d-1 /)
      y(:,7) = (/ 2.10d-1, 2.53d-1, 2.93d-1, 3.32d-1, 3.69d-1, &
                  4.06d-1, 4.41d-1, 4.75d-1, 5.09d-1, 5.43d-1, &
                  5.76d-1, 6.08d-1, 6.40d-1, 6.72d-1, 7.03d-1, &
                  7.34d-1 /)
      y(:,8) = (/ 2.98d-1, 3.59d-1, 4.16d-1, 4.70d-1, 5.22d-1, &
                  5.73d-1, 6.21d-1, 6.68d-1, 7.15d-1, 7.60d-1, &
                  8.04d-1, 8.47d-1, 8.90d-1, 9.32d-1, 9.73d-1, &
                  1.01d0 /)
      y(:,9) = (/ 3.65d-1, 4.39d-1, 5.09d-1, 5.75d-1, 6.39d-1, &
                  7.00d-1, 7.58d-1, 8.15d-1, 8.71d-1, 9.25d-1, &
                  9.77d-1, 1.03d0, 1.08d0, 1.13d0, 1.18d0, 1.23d0 /)
      y(:,10) = (/ 4.58d-1, 5.50d-1, 6.37d-1, 7.21d-1, 8.00d-1, &
                   8.76d-1, 9.49d-1, 1.02d0, 1.09d0, 1.15d0, &
                   1.22d0, 1.28d0, 1.34d0, 1.40d0, 1.46d0, 1.52d0 /)
      y(:,11) = (/ 5.92d-1, 7.11d-1, 8.24d-1, 9.31d-1, 1.03d0, &
                   1.13d0, 1.23d0, 1.32d0, 1.40d0, 1.49d0, 1.57d0, &
                   1.65d0, 1.73d0, 1.80d0, 1.88d0, 1.95d0 /)
      y(:,12) = (/ 7.98d-1, 9.58d-1, 1.11d0, 1.25d0, 1.39d0, &
                   1.52d0, 1.65d0, 1.77d0, 1.89d0, 2.00d0, 2.11d0, &
                   2.21d0, 2.32d0, 2.42d0, 2.51d0, 2.61d0 /)
      y(:,13) = (/ 1.14d0, 1.36d0, 1.58d0, 1.78d0, 1.98d0, 2.17d0, &
                   2.34d0, 2.52d0, 2.68d0, 2.84d0, 3.00d0, 3.15d0, &
                   3.29d0, 3.43d0, 3.57d0, 3.70d0 /)
      y(:,14) = (/ 1.77d0, 2.11d0, 2.44d0, 2.75d0, 3.05d0, 3.34d0, &
                   3.62d0, 3.89d0, 4.14d0, 4.39d0, 4.63d0, 4.86d0, &
                   5.08d0, 5.30d0, 5.51d0, 5.71d0 /)
      y(:,15) = (/ 3.10d0, 3.71d0, 4.29d0, 4.84d0, 5.37d0, 5.87d0, &
                   6.36d0, 6.83d0, 7.28d0, 7.72d0, 8.14d0, 8.55d0, &
                   8.95d0, 9.33d0, 9.71d0, 1.01d1 /)
      y(:,16) = (/ 6.92d0, 8.27d0, 9.56d0, 1.08d1, 1.19d1, 1.31d1, &
                   1.42d1, 1.52d1, 1.62d1, 1.72d1, 1.82d1, 1.91d1, &
                   2.00d1, 2.09d1, 2.17d1, 2.25d1 /)
      y(:,17) = (/ 2.75d1, 3.29d1, 3.80d1, 4.28d1, 4.75d1, 5.19d1, &
                   5.62d1, 6.04d1, 6.45d1, 6.84d1, 7.23d1, 7.60d1, &
                   7.97d1, 8.32d1, 8.67d1, 9.01d1 /)

      !
      ! John (1988) table 1
      !

      ! Coefficient for long wavelength cross section calculation
      ya = (/ 0d0, 2.483346d3, -3.449889d3, 2.20004d3, -6.96271d2, &
              8.8283d1 /)
      yb = (/ 0d0, 2.85827d2, -1.158382d3, 2.427719d3, -1.8414d3, &
              4.44517d2 /)
      yc = (/ 0d0, -2.054291d3, 8.746523d3, -1.3651105d4, &
              8.624970d3, -1.863864d3 /)
      yd = (/ 0d0, 2.827776d3, -1.1485632d4, 1.6755524d4, &
              -1.0051530d4, 2.095288d3 /)
      ye = (/ 0d0, -1.341537d3, 5.303609d3, -7.510494d3, &
              4.400067d3, -9.01788d2 /)
      yf = (/ 0d0, 2.08952d2, -8.12939d2, 1.132738d3, -6.55020d2, &
              1.32985d2 /)

      ! Get the correspondent wavelength in nm
      lambda = 1d2/freq

      ! Initialize absorptivity
      eta = 0d0

      ! If below limits, no contribution
      if (lambda.lt.x2(1)) then

        ! Just return with the zero
        return

      ! If within limits, use table 1 from Stilley & Callaway (1970)
      else if (lambda.le.x2(n2)) then

        ! For each height
        do iz=iz0,iz1

          ! Theta [eV]
          tev = ktoev/T(iz)

          ! Energy is between the table limits
          if (tev.ge.x1(1).and.tev.le.x1(n1)) then

            ! Bilinear interpolation in tabulation
            call bilinear(x1,x2,y,tev,lambda,cc)

          ! Energy is above table limits
          else if (tev.lt.x1(1)) then

            ! Linear interpolation along the border of the tabulation
            call linear(x2,y(1,:),lambda,cc)

          ! Energy is below table limits
          else if (tev.gt.x1(n1)) then

            ! Linear interpolation along the border of the tabulation
            call linear(x2,y(n1,:),lambda,cc)

          end if ! Within the limits in the first dimension

          ! Calculate absorptivity
          ! 1d-19 = 1d-29 (m^2) * 1d12 (cm^-3**2 -> m^-3**2) *
          !         1d-2 (m^-1 -> cm^-1)
          eta(iz) = nh(iz)*1d-19*cc*ne(iz)*kb*T(iz)

        end do ! Heights to consider

      ! If above the wavelength limit, use table 1 from John (1988)
      else

        ! Convert lambda to microns
        lambda = lambda*1d-3
        ! And store its inverse
        ilambda = 1d0/lambda

        ! Initialize
        yg = 0d0

        ! Compute the cross section Theta coefficients using the table
        do ii=1,n3
          yg(ii) = lambda*lambda*ya(ii) + yb(ii) + ilambda* &
                   (yc(ii) + ilambda*(yd(ii) + ilambda* &
                    (ye(ii) + ilambda*yf(ii))))
        end do

        ! Multiplicative constant
        ! 1d-22 = 1d-34 cm^-1 * 1d12 cm^-3 -> m^-3
        C0 = kb*ktoev*1d-22

        ! At each height
        do iz=iz0,iz1

          ! Get tempeature in eV
          tev = sqrt(ktoev/T(iz))

          ! Add contributions to absorptivity
          cc = 1d0
          do ii=2,n3
            cc = cc*tev
            eta(iz) = eta(iz) + cc*yg(ii)
          end do

          ! Complete with number densities and units constant
          eta(iz) = eta(iz)*nh(iz)*ne(iz)*C0

        end do ! Heights

      end if ! Invalid, short, or long wavelength

      end subroutine Hminus_ff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes H bound-free absorptivity and emissivity\n
      !!      freq(double): Frequency\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!      T(double(:)): Temperature\n
      !!      iz0(integer): First height index to consider\n
      !!      iz1(integer): Last height index to consider\n
      !!    eta(double(:)): Absorptivity\n
      !!    eps(double(:)): Emissivity\n
      subroutine HI_bf(freq,Atom,T,iz0,iz1,eta,eps)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: T
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer:: iz,iterm,iphot,ilevel,iJ,nfr

      double precision:: neff,gbf,c0,sig,arg,gij,clambda,exu,pE
      double precision, dimension(iz0:iz1):: nps,np,nts,nt


      ! Frequency related constant quantities
      arg = c2*freq*1d4
      pE = convF*2d21*c*freq*freq*freq

      ! Cross section constant part
      ! 1d-5 = 1d-9 (real_c -> c) * 1d4 (m^2 -> cm^2)
      c0 = 16d-5*qel*qel*hplanck/sqrt(27d0)/pi4eps0/me/c/ryd/fktoJ

      ! LTE population of HII (it must be last level for Hydrogen
      ! atom)
      nps = Atom%populte(Atom%nlevel,iz0:iz1)

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! For each photoionization
      do iphot=1,Atom%nphot

        ! Get the index with the maximum frequency depending on the
        ! type of input
        nfr = (1 - Atom%phot(iphot)%mode)*Atom%phot(iphot)%nfreq + &
              Atom%phot(iphot)%mode

        ! If the frequency is out of the range of this transition,
        ! skip it
        if (freq.lt.Atom%phot(iphot)%edge.or. &
            freq.gt.Atom%phot(iphot)%infreq(nfr)) cycle

        ! Get lower level of this photoionization
        ilevel = Atom%phot(iphot)%ilevell

        ! LTE population of the lower level
        nts = Atom%populte(ilevel,iz0:iz1)

        ! Current population of the lower level
        nt = Atom%popu(ilevel,iz0:iz1)

        ! Get the term and sublevel indexes of the lower level
        iterm = Atom%term(ilevel)
        iJ = Atom%sublevel(ilevel)

        ! Calculate the effective principal quantum number
        neff = sqrt(ryd/ &
                    (Atom%FSfreq(Atom%nJ(Atom%nMulti),Atom%nMulti) - &
                     Atom%FSfreq(iJ,iterm)))

        ! Gaunt factor
        gbf = gHI_bf(freq,neff,1d0)

        ! Quotient between edge and current frequencies
        clambda = Atom%phot(iphot)%edge/freq

        ! Compute cross section
        sig = c0*neff*gbf*clambda*clambda*clambda

        ! Curent Population of HII
        np = Atom%popu(Atom%nlevel,iz0:iz1)

        ! For each height
        do iz=iz0,iz1

          ! Get exponential
          exu = arg/T(iz)
          exu = diexp(exu)

          ! Compute absortivity and emissivity
          gij = nts(iz)*exu/nps(iz)
          eta(iz) = sig*(1d0 - exu)*nt(iz) + eta(iz)
          eps(iz) = pE*gij*sig*np(iz) + eps(iz)

        end do ! Heights
      end do ! Photionizations

      end subroutine HI_bf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Hydrogen bound-free Gaunt factor\n
      !!  freq(double): Frequency\n
      !!     n(double): Principal quantum number\n
      !!     Z(double): Ion charge
      double precision function gHI_bf(freq,n,Z)

      ! I/O

      double precision, intent(in):: freq,n,Z

      ! Local

      double precision:: x, y


      ! Conver to Eq. units
      x = freq/ryd/Z/Z
      y = 1d0/x/n/n
      x = x**(1d0/3d0)

      ! Calculate Gaunt factor
      gHI_bf = 1d0 + .1728d0*x*(1d0 - 2d0*y) - .0496d0*x*x* &
               (1d0 - (1d0 - y)*2d0*y/3d0)

      end function gHI_bf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes H free-free absorptivity\n
      !!      freq(dfloat): Frequency\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!      T(double(:)): Temperature\n
      !!     ne(double(:)): Electron number density\n
      !!      iz0(integer): First height index to consider\n
      !!      iz1(integer): Last height index to consider\n
      !!    eta(double(:)): Absorptivity
      subroutine HI_ff(freq,Atom,T,ne,iz0,iz1,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: T,ne
      double precision, dimension(iz0:iz1), intent(out):: eta

      ! Local

      integer:: iz

      double precision:: freq3,arg,C0,exu,Tin,gff
      double precision, dimension(iz0:iz1):: np


      ! Frequency related quantities
      arg = c2*freq*1d4
      freq3 = freq*freq*freq*c*c*c*1d48
      freq3 = 1d0/freq3

      ! Constant part of cross section
      c0 = qel*qel/pi4eps0/sqrt(me)
      ! 1d1 = 1d-9 (real_c -> c) * 1d12 (m^6 -> cm^6) *
      !       1d-2 (m^-1 -> cm^-1)
      c0 = sqrt(32d0*pi/27d0/kb)*c0*c0*c0*1d1/hplanck/c

      ! Population of HII (for Hydrogen model it must be last level)
      np = Atom%popu(Atom%nlevel,iz0:iz1)

      ! For each height to consider
      do iz=iz0,iz1

        ! Calculate exponential
        Tin = 1d0/T(iz)
        exu = arg*Tin
        exu = diexp(exu)

        ! Get Gaunt factor
        gff = gHI_ff(freq,1d0,T(iz))

        ! Calculate absorptivity
        eta(iz) = c0*gff*(1d0 - exu)*freq3*ne(iz)*np(iz)*sqrt(Tin)

      end do ! Heights

      end subroutine HI_ff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Hydrogen free-free Gaunt factor\n
      !!   freq(double): Frequency\n
      !!      Z(double): Ion charge\n
      !!      T(double): Temperature
      double precision function gHI_ff(freq,Z,T)

      ! I/O

      double precision, intent(in):: freq, Z, T

      ! Local

      double precision:: x, y


      ! Conver to Eq. units
      x = freq/ryd/Z/Z
      x = x**(1d0/3d0)
      y = 2d0*kb*T*1d7/convF/freq

      ! Calculate Gaunt factor
      gHI_ff = 1d0 + .1728d0*x*(1d0 + y) - .0496d0*x*x* &
               (1d0 + (1d0 + y)*y/3d0)

      ! Cannot be smaller than one
      if (gHI_ff.lt.1d0) gHI_ff = 1d0

      end function gHI_ff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes atomic bound-free absorptivity and emissivity\n
      !!      freq(double): Frequency\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!      T(double(:)): Temperature\n
      !!      iz0(integer): First height index to consider\n
      !!      iz1(integer): Last height index to consider\n
      !!    eta(double(:)): Absorptivity\n
      !!    eps(double(:)): Emissivity
      subroutine back_bf(freq,Atom,T,iz0,iz1,eta,eps)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: T
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer:: iz,iterm,iterm1,iphot,ilevel,ilevel1
      integer:: iJ,iJ1,nfr

      double precision:: neff,gbf,gbfe,sig,arg,gij,clambda, pE, Z
      double precision, dimension(iz0:iz1):: nps,np,nts,nt,exu
      double precision, dimension(:), allocatable:: sb,sc,sd,sx,sy

      ! Frequency quantities
      arg = c2*freq*1d4
      pE = convF*2d21*c*freq*freq*freq

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! Calculate exponential for every considered height
      do iz=iz0,iz1
        Z = arg/T(iz)
        exu(iz) = diexp(Z)
      end do

      ! For each photoionization
      do iphot=1,Atom%nphot

        ! Identify levels
        ilevel1 = Atom%phot(iphot)%ilevelu
        ilevel = Atom%phot(iphot)%ilevell

        ! Get the index of maximum frequency depending on the input
        ! mode
        nfr = (1 - Atom%phot(iphot)%mode)*Atom%phot(iphot)%nfreq + &
              Atom%phot(iphot)%mode

        ! If the frequency is out of the ranges of this transition,
        ! skip
        if (freq.lt.Atom%phot(iphot)%edge.or. &
            freq.gt.Atom%phot(iphot)%infreq(nfr)) cycle

        ! Determine population lower level
        nts = Atom%populte(ilevel,iz0:iz1)
        nt = Atom%popu(ilevel,iz0:iz1)

        ! Determine population upper level
        nps = Atom%populte(ilevel1,iz0:iz1)
        np = Atom%popu(ilevel1,iz0:iz1)

        ! If the input was explicit
        if (Atom%phot(iphot)%mode.eq.0) then

          ! Ensure we have space for the spline coefficients
          if (allocated(sx)) then
            if (size(sx).lt.nfr) then
              deallocate(sx)
              deallocate(sy)
              deallocate(sb)
              deallocate(sc)
              deallocate(sd)
              allocate(sx(nfr))
              allocate(sy(nfr))
              allocate(sb(nfr))
              allocate(sc(nfr))
              allocate(sd(nfr))
            end if
          else
            allocate(sx(nfr))
            allocate(sy(nfr))
            allocate(sb(nfr))
            allocate(sc(nfr))
            allocate(sd(nfr))
          end if

          ! Store data in sx and sy variables
          sx(nfr:1:-1) = 1d2/Atom%phot(iphot)%infreq
          sy(nfr:1:-1) = Atom%phot(iphot)%inalpha

          ! Interpolate cross section with cubic splines in cm^2
          call spline(sx(1:nfr),sy(1:nfr), &
                      sb(1:nfr),sc(1:nfr),sd(1:nfr),nfr)
          sig = ispline(1d2/freq,sx(1:nfr),sy(1:nfr), &
                        sb(1:nfr),sc(1:nfr),sd(1:nfr),nfr)*1d4

          ! If negative cross-section
          if (sig.lt.0) then

            ! Use linear interpolation
            call linear(sx(1:nfr),sy(1:nfr), &
                        1d2/freq,sig)
            sig = sig*1d4

          end if ! Failed splines

        ! If the input was hidrogenic
        else

          ! Find the term and J indexes for the levels involved
          iterm1 = Atom%term(ilevel1)
          iJ1 = Atom%sublevel(ilevel1)
          iterm = Atom%term(ilevel)
          iJ = Atom%sublevel(ilevel)

          ! Get charge of the upper level ion
          Z = dble(Atom%stage(iterm1) - 1)

          ! Compute effective principal quantum number
          neff = Z*sqrt(ryd/ &
                 (Atom%FSfreq(iJ1,iterm1) - Atom%FSfreq(iJ,iterm)))

          ! Compute Gaunt factors at this frequency and at edge
          gbf  = gHI_bf(freq,neff,Z)
          gbfe = gHI_bf(Atom%phot(iphot)%edge,neff,Z)

          ! Quotient between edge and current frequencies
          clambda = Atom%phot(iphot)%edge/freq

          ! Cross section at this frequency
          sig = Atom%phot(iphot)%inalpha(1)*1d4* &
                clambda*clambda*clambda*gbf/gbfe

        end if ! Explicit or hydrogenic

        ! For each height to consider
        do iz=iz0,iz1

          ! Compute absorptivity
          eta(iz) = sig*(1d0 - exu(iz))*nt(iz) + eta(iz)

          ! If there is population in the upper level
          if (nps(iz).gt.0d0) then

            ! Compute emissivity
            gij = nts(iz)*exu(iz)/nps(iz)
            eps(iz) = pE*gij*sig*np(iz) + eps(iz)

          end if ! Population in upper level

        end do ! Heights
      end do ! Photoionizations

      ! Free
      if (allocated(sx)) deallocate(sx,sy,sb,sc,sd)

      return

      end subroutine back_bf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes atomic bound-bound absorptivity and emissivity\n
      !!        freq(double): Frequency\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!        T(double(:)): Temperature\n
      !!      vmi(double(:)): Microturbulent velocity\n
      !!         DwT(double): Thermal Doppler width to estimate
      !!                      ranges\n
      !!     vfac(double(:)): Doppler shift factor\n
      !!        iz0(integer): First height index to consider\n
      !!        iz1(integer): Last height index to consider\n
      !!      fline(logical): Bool that tells if a line was found in
      !!                      this frequency\n
      !!      eta(double(:)): Absorptivity\n
      !!      eps(double(:)): Emissivity
      subroutine back_bb(freq,Atom,Flgsg,T,vmi,DwT,vfac,iz0,iz1, &
                         fline,eta,eps)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(out):: fline
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq, DwT
      double precision, dimension(iz0:iz1), intent(in):: T, vmi, vfac
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer:: itran,fitran,iterml,itermu,iJl,iJu,ill,ilu,iz

      double precision:: freq0,dfreq,Aul,gul
      double precision:: rJl,rJu,rLu,rLl,rS
      double precision:: c0,c1,pEl,sqrtpi,Dw,adamp,sig,prof

      ! Initialize the line found variable
      fline = .False.

      ! hc/4pi and square root of pi
      c0 = hplanck*c*1d12*.25d0/PI
      sqrtpi = sqrt(PI)

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! For every transition
      do itran=1,Atom%ntran

        ! Get terms
        itermu = Atom%fst(itran)%itermu
        iterml = Atom%fst(itran)%iterml

        ! Get quantum numbers of the terms
        rLu = Atom%rLval(itermu)
        rLl = Atom%rLval(iterml)
        rS = Atom%Sval(itermu)

        ! If ML atom
        if (Atom%ML) then

          ! Get frequency of the transition
          freq0 = Atom%FSfreq(1,itermu) - &
                  Atom%FSfreq(1,iterml)

          ! Get the range this transition fills
          dfreq = freq0*Atom%Dwvl(itran)*DwT

          ! If the frequency is out of range, skip
          if (abs(freq - freq0).gt.dfreq) cycle

          ! If not, we have found a line
          fline = .True.

          ! Einstein coefficient, quotient of statistical
          ! weights, and energy
          Aul = Atom%Ecoeff(itermu,iterml)
          gul = (2d0*rLu + 1d0)/(2d0*rLl + 1d0)
          pEl = convF*2d21*c*freq0*freq0*freq0

          ! Constant part of the RT coefficients
          ! 1d14 = 1d9 (m^3 -> cm^3) * 1d5 (10^5 cm^-1 -> cm^-1)
          c1 = c0*freq0*1d14*Aul/sqrtpi

          ! For each height
          do iz=iz0,iz1

            ! Doppler width
            Dw = Atom%cDopp*sqrt(T(iz))
            Dw = freq0*sqrt(Dw*Dw + vmi(iz)**2d0)

            ! Damping parameter
            adamp = (Atom%damp(itermu,iz) + &
                     Atom%damp(iterml,iz) + &
                     Atom%ldamp(itran,iz))/Dw

            ! Voigt profile
            call voigtI((freq0 - freq*vfac(iz))/Dw,adamp,prof)

            ! 'cross-section'
            sig = c1*prof/Dw

            ! Absorptivity
            eta(iz) = eta(iz) + sig* &
                      (Atom%popu(iterml,iz)*gul - &
                       Atom%popu(itermu,iz))/pEl

            ! Emissivity
            eps(iz) = eps(iz) + sig*Atom%popu(itermu,iz)

          end do !heights

        ! If MT atom
        else

          ! Check it is allowed (spin conserved)
          if (abs(rS - atom%Sval(iterml)).gt..1d0) cycle

          ! For each FS transition
          do fitran=1,Atom%fst(itran)%nt

            ! Get sublevel indexes
            iJu = Atom%fst(itran)%ilevelu(fitran)
            iJl = Atom%fst(itran)%ilevell(fitran)

            ! Get level index
            ilu = Atom%irho(itermu)%irho_ij(iJu)
            ill = Atom%irho(iterml)%irho_ij(iJl)

            ! Get total angular momentum
            rJu = Atom%rJval(iJu,itermu)
            rJl = Atom%rJval(iJl,iterml)

            ! Get frequency of the transition
            freq0 = Atom%FSfreq(iju,itermu) - &
                    Atom%FSfreq(ijl,iterml)

            ! Get the range this transition holds
            dfreq = freq0*Atom%Dwvl(itran)*DwT

            ! If the frequency is out of range, skip
            if (abs(freq - freq0).gt.dfreq) cycle

            ! If not, we have found a line
            fline = .True.

            ! Einstein coefficient
            Aul = fun6j(rLu,rLl,1d0,rJl,rJu,rS,Flgsg)
            Aul = (2d0*rLu+1d0)*(2d0*rJl+1d0)*Aul*Aul* &
                  Atom%Ecoeff(itermu,iterml)
           !Aul = Atom%fst(itran)%Aul(iJu,iJl)

            ! Quotient statistical weights and energy
            gul = (2d0*rJu + 1d0)/(2d0*rJl + 1d0)
            pEl = convF*2d21*c*freq0*freq0*freq0

            ! Constant part of the RT coefficients
            ! 1d14 = 1d9 (m^3 -> cm^3) * 1d5 (10^5 cm^-1 -> cm^-1)
            c1 = c0*freq0*1d14*Aul/sqrtpi

            ! For each height
            do iz=iz0,iz1

              ! Doppler width
              Dw = Atom%cDopp*sqrt(T(iz))
              Dw = freq0*sqrt(Dw*Dw + vmi(iz)**2d0)

              ! Damping parameter
              adamp = (Atom%damp(itermu,iz) + &
                       Atom%damp(iterml,iz) + &
                       Atom%ldamp(itran,iz))/Dw

              ! Voigt profile
              call voigtI((freq0 - freq*vfac(iz))/Dw,adamp,prof)

              ! 'cross-section'
              sig = c1*prof/Dw

              ! Absorptivity
              eta(iz) = eta(iz) + sig* &
                        (Atom%popu(ill,iz)*gul - &
                         Atom%popu(ilu,iz))/pEl

              ! Emissivity
              eps(iz) = eps(iz) + sig*Atom%popu(ilu,iz)

            end do !heights
          end do ! FS transition

        end if ! ML or MT

      end do ! Transitions

      return

      end subroutine back_bb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes hydrogen bound-bound absorptivity and emissivity
      !! with the hard-coded model\n
      !!      freq(double): Frequency\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!      T(double(:)): Temperature\n
      !!    vmi(double(:)): Microturbulent velocity\n
      !!       DwT(double): Thermal Doppler width to estimate ranges\n
      !!   vfac(double(:)): Doppler shift factor\n
      !!      iz0(integer): First height index to consider\n
      !!      iz1(integer): Last height index to consider\n
      !!    fline(logical): Bool that tells if a line was found in
      !!                    this frequency\n
      !!    eta(double(:)): Absorptivity\n
      !!    eps(double(:)): Emissivity
      subroutine backH_bb(freq,Atom,T,vmi,DwT,vfac,iz0,iz1,fline, &
                          eta,eps)

      ! I/O

      type(Atom_class), intent(in):: Atom
      logical, intent(out):: fline
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq, DwT
      double precision, dimension(iz0:iz1), intent(in):: T, vmi, vfac
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer:: itran,iterml,itermu,iz

      double precision:: freq0,dfreq,Aul,gul
      double precision:: c0,c1,pEl,sqrtpi,Dw,adamp,sig,prof


      ! Initialize the line found variable
      fline = .False.

      ! hc/4pi and square root of pi
      c0 = hplanck*c*1d12*.25d0/PI
      sqrtpi = sqrt(PI)

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! Run over transitions
      do itran=1,Atom%ntran

        ! Get terms
        itermu = Atom%fst(itran)%itermu
        iterml = Atom%fst(itran)%iterml

        ! Get frequency of the transition
        freq0 = Atom%Dfreq(itran)

        ! Get the range this transition fills
        dfreq = freq0*Atom%Dwvl(itran)*DwT

        ! If the frequency is out of range, skip
        if (abs(freq - freq0).gt.dfreq) cycle

        ! If not, we have found a line
        fline = .True.

        ! Einstein coefficient, quotient of statistical weights, and
        ! energy
        Aul = Atom%Ecoeff(itermu,iterml)
        gul = (2d0*Atom%rJval(1,itermu) + 1d0)/ &
              (2d0*Atom%rJval(1,iterml) + 1d0)
        pEl = convF*2d21*c*freq0*freq0*freq0

        ! Constant part of the RT coefficients
        ! 1d14 = 1d9 (m^3 -> cm^3) * 1d5 (10^5 cm^-1 -> cm^-1)
        c1 = c0*freq0*1d14*Aul/sqrtpi

        ! For each height
        do iz=iz0,iz1

          ! Doppler width
          Dw = Atom%cDopp*sqrt(T(iz))
          Dw = freq0*sqrt(Dw*Dw + vmi(iz)**2d0)

          ! Damping parameter
          adamp = (Atom%damp(itermu,iz) + &
                   Atom%damp(iterml,iz) + &
                   Atom%ldamp(itran,iz))/Dw

          ! Voigt profile
          call voigtI((freq0 - freq*vfac(iz))/Dw,adamp,prof)

          ! 'cross-section'
          sig = c1*prof/Dw

          ! Absorptivity
          eta(iz) = eta(iz) + sig* &
                    (Atom%popu(iterml,iz)*gul - &
                     Atom%popu(itermu,iz))/pEl

          ! Emissivity
          eps(iz) = eps(iz) + sig*Atom%popu(itermu,iz)

        end do ! Heights
      end do ! Transitions

      return

      end subroutine backH_bb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Rayleigh scattering coefficient for hydrogen\n
      !!      freq(double): Frequency\n
      !!    ifreq(integer): Frequency index\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!       DwT(double): Thermal Doppler width to estimate ranges\n
      !!      iz0(integer): First height index to consider\n
      !!      iz1(integer): Last height index to consider\n
      !!      act(logical): Bool that tells if hydrogen is active\n
      !!    sig(double(:)): Scattering coefficient
      subroutine rayleigh(freq,ifreq,Atom,DwT,iz0,iz1,act,sig)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz0,iz1,ifreq,act
      double precision, intent(in):: freq, DwT
      double precision, dimension(iz0:iz1), intent(out):: sig

      ! Local

      integer:: iterm1, iz, itran, ij

      double precision:: c0,CfA,f,fij,deg0,deg1
      double precision:: lambda0,freq0,flimit,icfreq
      double precision, dimension(iz0:iz1):: n0
      double precision, dimension(Atom%ntran):: freqR


      ! Initialize RT coefficient
      sig = 0d0

      ! If there are no transitions, no contribution either
      if (Atom%ntran.eq.0) return

      ! Upper limit of frequency
      flimit = 1d-4

      ! If the Hydrogen is not active, look for the red limits
      ! of the lines of the Lyman series in the hard-coded atom
      if (act.eq.0) then

        ! For every upper term
        do iterm1=2,Atom%nMulti

          ! If not connected to ground term
          if (Atom%irad(1,iterm1).lt.1) cycle

          ! Get transition index
          itran = Atom%irad(1,iterm1)

          ! Frequency of the transition
          freq0 = Atom%Dfreq(itran)

          ! Red side frequency
          freqR(itran) = freq0/(1d0 + DwT*Atom%Dwvl(itran))

          ! Update the limit
          flimit = max(flimit,freqR(itran))

        end do ! Upper terms

      ! If the atom is active
      else

        ! For every upper term
        do iterm1=2,Atom%nMulti

          ! If not connected to ground term
          if (Atom%irad(1,iterm1).lt.1) cycle

          ! Get transition index
          itran = Atom%irad(1,iterm1)

          ! Update the limit
          flimit = max(flimit,Atom%freqR(itran))

        end do ! Upper terms

      end if ! Hard-coded or active

      ! If the frequency is larger than the larger red wing, no
      ! contribution
      if (freq.ge.flimit) return

      ! Initialize the population of the ground term
      n0 = 0d0

      ! Calculate the population of the ground term
      do ij=1,Atom%nJ(1)
        n0 = n0 + Atom%popu(ij,iz0:iz1)
      end do

      ! Constants
      CfA = 2d0*pi*qel*qel*1d-9/eps0/me/c
      c0 = qel*qel*1d-18/me/c/c/pi4eps0
      c0 = 8d0*pi*c0*c0/3d0

      ! Degeneration of the ground term
      deg0 = (2d0*Atom%rLval(1) + 1d0)*(2d0*Atom%Sval(1) + 1d0)

      ! Initialize Rayleigh coefficient
      fij = 0d0

      ! For each upper term in the Lyman series
      do iterm1=2,Atom%nMulti

        ! Get transition
        itran = Atom%irad(1,iterm1)

        ! If no transition, skip
        if (itran.lt.1) cycle

        ! Get Transition frequency
        freq0 = Atom%Dfreq(itran)

        ! If not active
        if (act.eq.0) then

          ! Check that the frequency is lower than the red wing of
          ! the line previously determined
          if (freq.ge.freqR(itran)) cycle

        ! If it is active
        else

          ! Check that the frequency is lower than the transition
          ! frequency and that the line is absent in this frequency
          if (Atom%fflag(itran)%absent.and.freq.gt.freq0) cycle
          if (ifreq.gt.0) then
            if (ifreq.ge.Atom%if0(itran)) cycle
          else
            if (freq.ge.Atom%freqR(itran)) cycle
          end if

        end if ! Active or passive atom

        ! Degeneration of upper term
        deg1 = (2d0*Atom%rLval(iterm1) + 1d0)* &
               (2d0*Atom%Sval(iterm1) + 1d0)

        ! Frequency quantities entering the calculation
        icfreq = freq0/freq
        icfreq = 1d0/(icfreq*icfreq - 1d0)
        lambda0 = 1d-7/freq0

        ! Oscillator strength
        f = Atom%Ecoeff(iterm1,1)*1d8* &
            lambda0*lambda0*deg1/deg0/CfA

        ! Add to Rayleigh coefficient
        fij = fij + f*icfreq*icfreq

      end do ! Upper terms

      ! Compute Rayleigh cross section
      ! 1d4 (m^2 -> cm^2)
      sig = fij*c0*1d4

      ! For each height, compute scattering contribution
      do iz=iz0,iz1
        sig(iz) = sig(iz)*n0(iz)
      end do

      return

      end subroutine rayleigh

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes Rayleigh scattering coefficient for H2\n
      !!     freq(double): Frequency\n
      !!   nH2(double(:)): H2 density\n
      !!     iz0(integer): First height index to consider\n
      !!     iz1(integer): Last height index to consider\n
      !!   sig(double(:)): Scattering coefficient
      subroutine rayleigh_H2(freq,nH2,iz0,iz1,sig)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: nH2
      double precision, dimension(iz0:iz1), intent(out):: sig

      ! Local

      integer, parameter:: n=21
      integer:: iz

      double precision:: C1,C2,C3
      double precision:: lambda, CC
      double precision, dimension(n):: x
      double precision, dimension(n):: y


      ! Initialize coefficients
      C1 = 8.779d1
      C2 = 1.323d6
      C3 = 2.245d10

      ! Lambda [nm]
      x = (/ 121.57d0, 130d0, 140d0, 150d0, 160d0, 170d0, &
             185.46d0, 186.27d0, 193.58d0, 199.05d0, 230.29d0, &
             237.91d0, 253.56d0, 275.36d0, 296.81d0, 334.24d0, &
             404.77d0, 407.9d0, 435.96d0, 546.23d0, 632.8d0 /)

      ! Cross section [1d18 cm^2]
      y = (/ 2.35d-6, 1.22d-6, 6.8d-7, 4.24d-7, 2.84d-7, &
             2.00d-7, 1.25d-7, 1.22d-7, 1d-7, 8.7d-8, 4.29d-8, &
             3.68d-8, 2.75d-8, 1.89d-8, 1.36d-8, 8.11d-9, &
             3.6d-9, 3.48d-9, 2.64d-9, 1.04d-9, 5.69d-10 /)

      ! Get the wavelength in nm
      lambda = 1d2/freq

      ! Initialize RT coefficient
      sig = 0d0

      ! If below limits, no contribution
      if (lambda.lt.x(1)) return

      ! If within limits
      if (lambda.le.x(n)) then

        ! Linear interpolation of cross section
        call linear(x,y,lambda,CC)

      ! If above limits, use formula
      else

        lambda = lambda*lambda
        lambda = 1d0/lambda

        CC = lambda*lambda*(C1 + lambda*(C2 + lambda*C3))

      end if ! Within or above the limit

      ! Convert cross section to cm^2
      CC = CC*1d-18

      ! For each height compute scattering contribution
      do iz=iz0,iz1
        sig(iz) = CC*nH2(iz)
      end do

      return

      end subroutine rayleigh_H2

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes OH bound-free absorptivity and emissivity\n
      !!    freq(double): Frequency\n
      !!  nOH(double(:)): OH density\n
      !!    T(double(:)): Temperature\n
      !!    iz0(integer): First height index to consider\n
      !!    iz1(integer): Last height index to consider\n
      !!  eta(double(:)): Absorptivity\n
      !!  eps(double(:)): Emissivity
      subroutine OH_bf(freq,nOH,T,iz0,iz1,eta,eps)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: nOH, T
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer, parameter:: n2=15, n1=130
      integer:: iz

      double precision:: arg, pE, pEeV, cc, exu
      double precision, dimension(n1):: x1
      double precision, dimension(n2):: x2
      double precision, dimension(n2,n1):: y


      ! Kurucz, van Dishoeck, and Tarafdar (1987) table 1
      x1 = (/ 2.1d0,  2.2d0,  2.3d0,  2.4d0,  2.5d0,  2.6d0,  2.7d0, &
              2.8d0,  2.9d0,  3.0d0, 3.1d0,  3.2d0,  3.3d0,  3.4d0, &
              3.5d0,  3.6d0,  3.7d0,  3.8d0,  3.9d0,  4.0d0, 4.1d0, &
              4.2d0,  4.3d0,  4.4d0,  4.5d0,  4.6d0,  4.7d0,  4.8d0, &
              4.9d0,  5.0d0, 5.1d0,  5.2d0,  5.3d0,  5.4d0,  5.5d0, &
              5.6d0,  5.7d0,  5.8d0,  5.9d0,  6.0d0, 6.1d0,  6.2d0, &
              6.3d0,  6.4d0,  6.5d0,  6.6d0,  6.7d0,  6.8d0,  6.9d0, &
              7.0d0, 7.1d0,  7.2d0,  7.3d0,  7.4d0,  7.5d0,  7.6d0, &
              7.7d0,  7.8d0,  7.9d0,  8.0d0, 8.1d0,  8.2d0,  8.3d0, &
              8.4d0,  8.5d0,  8.6d0,  8.7d0,  8.8d0,  8.9d0,  9.0d0, &
              9.1d0,  9.2d0,  9.3d0,  9.4d0,  9.5d0,  9.6d0,  9.7d0, &
              9.8d0,  9.9d0, 10.0d0, 10.1d0, 10.2d0, 10.3d0, 10.4d0, &
              10.5d0, 10.6d0, 10.7d0, 10.8d0, 10.9d0, 11.0d0, &
              11.1d0, 11.2d0, 11.3d0, 11.4d0, 11.5d0, 11.6d0, &
              11.7d0, 11.8d0, 11.9d0, 12.0d0, 12.1d0, 12.2d0, &
              12.3d0, 12.4d0, 12.5d0, 12.6d0, 12.7d0, 12.8d0, &
              12.9d0, 13.0d0, 13.1d0, 13.2d0, 13.3d0, 13.4d0, &
              13.5d0, 13.6d0, 13.7d0, 13.8d0, 13.9d0, 14.0d0, &
              14.1d0, 14.2d0, 14.3d0, 14.4d0, 14.5d0, 14.6d0, &
              14.7d0, 14.8d0, 14.9d0, 15.0d0 /)
      x2 = (/ 2d3, 2.5d3, 3d3, 3.5d3, 4d3, 4.5d3, 5d3, 5.5d3, 6d3, &
              6.5d3, 7d3, 7.5d3, 8d3, 8.5d3, 9d3 /)
      y(:,1) = (/ -30.855d0,-29.121d0,-27.976d0,-27.166d0,-26.566d0, &
                  -26.106d0,-25.742d0,-25.448d0,-25.207d0,-25.006d0, &
                  -24.836d0,-24.691d0,-24.566d0,-24.457d0,-24.363d0 /)
      y(:,2) = (/ -30.494d0,-28.760d0,-27.615d0,-26.806d0,-26.206d0, &
                  -25.745d0,-25.381d0,-25.088d0,-24.846d0,-24.645d0, &
                  -24.475d0,-24.330d0,-24.205d0,-24.097d0,-24.002d0 /)
      y(:,3) = (/ -30.157d0,-28.425d0,-27.280d0,-26.472d0,-25.872d0, &
                  -25.411d0,-25.048d0,-24.754d0,-24.513d0,-24.312d0, &
                  -24.142d0,-23.997d0,-23.872d0,-23.764d0,-23.669d0 /)
      y(:,4) = (/ -29.848d0,-28.117d0,-26.974d0,-26.165d0,-25.566d0, &
                  -25.105d0,-24.742d0,-24.448d0,-24.207d0,-24.006d0, &
                  -23.836d0,-23.692d0,-23.567d0,-23.458d0,-23.364d0 /)
      y(:,5) = (/ -29.567d0,-27.837d0,-26.693d0,-25.885d0,-25.286d0, &
                  -24.826d0,-24.462d0,-24.169d0,-23.928d0,-23.727d0, &
                  -23.557d0,-23.412d0,-23.287d0,-23.179d0,-23.084d0 /)
      y(:,6) = (/ -29.307d0,-27.578d0,-26.436d0,-25.628d0,-25.029d0, &
                  -24.569d0,-24.205d0,-23.912d0,-23.671d0,-23.470d0, &
                  -23.300d0,-23.155d0,-23.031d0,-22.922d0,-22.828d0 /)
      y(:,7) = (/ -29.068d0,-27.341d0,-26.199d0,-25.391d0,-24.792d0, &
                  -24.332d0,-23.969d0,-23.676d0,-23.435d0,-23.234d0, &
                  -23.064d0,-22.920d0,-22.795d0,-22.687d0,-22.592d0 /)
      y(:,8) = (/ -28.820d0,-27.115d0,-25.978d0,-25.172d0,-24.574d0, &
                  -24.115d0,-23.752d0,-23.459d0,-23.218d0,-23.017d0, &
                  -22.848d0,-22.703d0,-22.579d0,-22.470d0,-22.376d0 /)
      y(:,9) = (/ -28.540d0,-26.891d0,-25.768d0,-24.968d0,-24.372d0, &
                  -23.914d0,-23.552d0,-23.259d0,-23.019d0,-22.818d0, &
                  -22.649d0,-22.504d0,-22.380d0,-22.272d0,-22.177d0 /)
      y(:,10) = (/-28.275d0,-26.681d0,-25.574d0,-24.779d0,-24.186d0, &
                  -23.729d0,-23.368d0,-23.076d0,-22.836d0,-22.636d0, &
                  -22.467d0,-22.322d0,-22.198d0,-22.090d0,-21.996d0 /)
      y(:,11) = (/-27.993d0,-26.470d0,-25.388d0,-24.602d0,-24.014d0, &
                  -23.560d0,-23.200d0,-22.909d0,-22.669d0,-22.470d0, &
                  -22.301d0,-22.157d0,-22.033d0,-21.925d0,-21.831d0 /)
      y(:,12) = (/-27.698d0,-26.252d0,-25.204d0,-24.433d0,-23.851d0, &
                  -23.401d0,-23.043d0,-22.754d0,-22.515d0,-22.316d0, &
                  -22.148d0,-22.005d0,-21.881d0,-21.773d0,-21.679d0 /)
      y(:,13) = (/-27.398d0,-26.026d0,-25.019d0,-24.267d0,-23.696d0, &
                  -23.251d0,-22.896d0,-22.609d0,-22.372d0,-22.174d0, &
                  -22.007d0,-21.864d0,-21.741d0,-21.634d0,-21.540d0 /)
      y(:,14) = (/-27.100d0,-25.791d0,-24.828d0,-24.102d0,-23.543d0, &
                  -23.106d0,-22.756d0,-22.472d0,-22.238d0,-22.041d0, &
                  -21.875d0,-21.733d0,-21.611d0,-21.504d0,-21.411d0 /)
      y(:,15) = (/-26.807d0,-25.549d0,-24.631d0,-23.933d0,-23.391d0, &
                  -22.964d0,-22.621d0,-22.341d0,-22.109d0,-21.915d0, &
                  -21.751d0,-21.610d0,-21.488d0,-21.383d0,-21.290d0 /)
      y(:,16) = (/-26.531d0,-25.310d0,-24.431d0,-23.761d0,-23.238d0, &
                  -22.823d0,-22.488d0,-22.214d0,-21.986d0,-21.795d0, &
                  -21.633d0,-21.494d0,-21.374d0,-21.269d0,-21.178d0 /)
      y(:,17) = (/-26.239d0,-25.066d0,-24.225d0,-23.585d0,-23.082d0, &
                  -22.681d0,-22.356d0,-22.089d0,-21.866d0,-21.679d0, &
                  -21.520d0,-21.383d0,-21.265d0,-21.162d0,-21.072d0 /)
      y(:,18) = (/-25.945d0,-24.824d0,-24.017d0,-23.405d0,-22.923d0, &
                  -22.538d0,-22.223d0,-21.964d0,-21.748d0,-21.565d0, &
                  -21.410d0,-21.276d0,-21.160d0,-21.059d0,-20.970d0 /)
      y(:,19) = (/-25.663d0,-24.587d0,-23.810d0,-23.222d0,-22.761d0, &
                  -22.391d0,-22.088d0,-21.838d0,-21.629d0,-21.452d0, &
                  -21.300d0,-21.170d0,-21.057d0,-20.958d0,-20.872d0 /)
      y(:,20) = (/-25.372d0,-24.350d0,-23.603d0,-23.038d0,-22.596d0, &
                  -22.241d0,-21.950d0,-21.710d0,-21.508d0,-21.337d0, &
                  -21.190d0,-21.064d0,-20.954d0,-20.858d0,-20.774d0 /)
      y(:,21) = (/-25.076d0,-24.111d0,-23.396d0,-22.853d0,-22.429d0, &
                  -22.088d0,-21.809d0,-21.578d0,-21.384d0,-21.220d0, &
                  -21.078d0,-20.957d0,-20.851d0,-20.758d0,-20.676d0 /)
      y(:,22) = (/-24.779d0,-23.870d0,-23.189d0,-22.669d0,-22.261d0, &
                  -21.934d0,-21.667d0,-21.445d0,-21.259d0,-21.101d0, &
                  -20.965d0,-20.848d0,-20.746d0,-20.656d0,-20.578d0 /)
      y(:,23) = (/-24.486d0,-23.629d0,-22.983d0,-22.486d0,-22.095d0, &
                  -21.781d0,-21.524d0,-21.311d0,-21.132d0,-20.980d0, &
                  -20.850d0,-20.737d0,-20.639d0,-20.553d0,-20.478d0 /)
      y(:,24) = (/-24.183d0,-23.382d0,-22.774d0,-22.302d0,-21.928d0, &
                  -21.627d0,-21.381d0,-21.177d0,-21.005d0,-20.859d0, &
                  -20.734d0,-20.625d0,-20.531d0,-20.449d0,-20.376d0 /)
      y(:,25) = (/-23.867d0,-23.127d0,-22.561d0,-22.116d0,-21.761d0, &
                  -21.474d0,-21.238d0,-21.043d0,-20.878d0,-20.738d0, &
                  -20.617d0,-20.513d0,-20.423d0,-20.344d0,-20.274d0 /)
      y(:,26) = (/-23.538d0,-22.862d0,-22.340d0,-21.926d0,-21.592d0, &
                  -21.320d0,-21.096d0,-20.909d0,-20.751d0,-20.617d0, &
                  -20.502d0,-20.402d0,-20.315d0,-20.239d0,-20.172d0 /)
      y(:,27) = (/-23.234d0,-22.604d0,-22.120d0,-21.734d0,-21.422d0, &
                  -21.166d0,-20.953d0,-20.776d0,-20.625d0,-20.497d0, &
                  -20.387d0,-20.291d0,-20.208d0,-20.135d0,-20.071d0 /)
      y(:,28) = (/-22.934d0,-22.347d0,-21.898d0,-21.541d0,-21.250d0, &
                  -21.010d0,-20.811d0,-20.643d0,-20.500d0,-20.378d0, &
                  -20.273d0,-20.182d0,-20.102d0,-20.033d0,-19.971d0 /)
      y(:,29) = (/-22.637d0,-22.092d0,-21.676d0,-21.345d0,-21.075d0, &
                  -20.853d0,-20.666d0,-20.508d0,-20.374d0,-20.259d0, &
                  -20.159d0,-20.073d0,-19.997d0,-19.931d0,-19.872d0 /)
      y(:,30) = (/-22.337d0,-21.835d0,-21.452d0,-21.147d0,-20.899d0, &
                  -20.693d0,-20.520d0,-20.373d0,-20.247d0,-20.139d0, &
                  -20.046d0,-19.964d0,-19.892d0,-19.830d0,-19.774d0 /)
      y(:,31) = (/-22.049d0,-21.584d0,-21.230d0,-20.950d0,-20.721d0, &
                  -20.531d0,-20.372d0,-20.236d0,-20.119d0,-20.019d0, &
                  -19.931d0,-19.855d0,-19.788d0,-19.729d0,-19.676d0 /)
      y(:,32) = (/-21.768d0,-21.337d0,-21.011d0,-20.754d0,-20.544d0, &
                  -20.370d0,-20.223d0,-20.098d0,-19.991d0,-19.898d0, &
                  -19.817d0,-19.746d0,-19.683d0,-19.628d0,-19.579d0 /)
      y(:,33) = (/-21.494d0,-21.096d0,-20.796d0,-20.559d0,-20.367d0, &
                  -20.208d0,-20.074d0,-19.960d0,-19.861d0,-19.776d0, &
                  -19.701d0,-19.636d0,-19.578d0,-19.527d0,-19.482d0 /)
      y(:,34) = (/-21.233d0,-20.861d0,-20.585d0,-20.368d0,-20.193d0, &
                  -20.048d0,-19.926d0,-19.821d0,-19.732d0,-19.654d0, &
                  -19.586d0,-19.526d0,-19.473d0,-19.426d0,-19.384d0 /)
      y(:,35) = (/-20.983d0,-20.635d0,-20.380d0,-20.181d0,-20.021d0, &
                  -19.889d0,-19.778d0,-19.683d0,-19.602d0,-19.531d0, &
                  -19.469d0,-19.415d0,-19.367d0,-19.324d0,-19.286d0 /)
      y(:,36) = (/-20.743d0,-20.418d0,-20.182d0,-19.999d0,-19.853d0, &
                  -19.733d0,-19.633d0,-19.547d0,-19.474d0,-19.410d0, &
                  -19.354d0,-19.305d0,-19.261d0,-19.223d0,-19.189d0 /)
      y(:,37) = (/-20.515d0,-20.210d0,-19.991d0,-19.824d0,-19.690d0, &
                  -19.581d0,-19.490d0,-19.413d0,-19.347d0,-19.290d0, &
                  -19.240d0,-19.196d0,-19.157d0,-19.122d0,-19.092d0 /)
      y(:,38) = (/-20.297d0,-20.011d0,-19.808d0,-19.654d0,-19.532d0, &
                  -19.434d0,-19.352d0,-19.282d0,-19.223d0,-19.172d0, &
                  -19.127d0,-19.088d0,-19.054d0,-19.023d0,-18.996d0 /)
      y(:,39) = (/-20.090d0,-19.822d0,-19.633d0,-19.491d0,-19.381d0, &
                  -19.291d0,-19.218d0,-19.156d0,-19.103d0,-19.057d0, &
                  -19.018d0,-18.983d0,-18.952d0,-18.925d0,-18.901d0 /)
      y(:,40) = (/-19.893d0,-19.642d0,-19.467d0,-19.337d0,-19.236d0, &
                  -19.155d0,-19.089d0,-19.034d0,-18.987d0,-18.946d0, &
                  -18.912d0,-18.881d0,-18.854d0,-18.831d0,-18.810d0 /)
      y(:,41) = (/-19.705d0,-19.472d0,-19.309d0,-19.190d0,-19.098d0, &
                  -19.025d0,-18.966d0,-18.917d0,-18.876d0,-18.840d0, &
                  -18.810d0,-18.783d0,-18.760d0,-18.739d0,-18.721d0 /)
      y(:,42) = (/-19.527d0,-19.310d0,-19.161d0,-19.051d0,-18.968d0, &
                  -18.903d0,-18.851d0,-18.807d0,-18.771d0,-18.740d0, &
                  -18.713d0,-18.690d0,-18.670d0,-18.653d0,-18.637d0 /)
      y(:,43) = (/-19.357d0,-19.159d0,-19.022d0,-18.922d0,-18.847d0, &
                  -18.789d0,-18.743d0,-18.704d0,-18.673d0,-18.646d0, &
                  -18.623d0,-18.603d0,-18.586d0,-18.571d0,-18.558d0 /)
      y(:,44) = (/-19.195d0,-19.016d0,-18.892d0,-18.803d0,-18.736d0, &
                  -18.684d0,-18.643d0,-18.610d0,-18.583d0,-18.560d0, &
                  -18.540d0,-18.523d0,-18.509d0,-18.496d0,-18.485d0 /)
      y(:,45) = (/-19.042d0,-18.883d0,-18.772d0,-18.693d0,-18.634d0, &
                  -18.589d0,-18.553d0,-18.525d0,-18.501d0,-18.481d0, &
                  -18.465d0,-18.451d0,-18.438d0,-18.428d0,-18.419d0 /)
      y(:,46) = (/-18.894d0,-18.758d0,-18.662d0,-18.593d0,-18.542d0, &
                  -18.503d0,-18.473d0,-18.448d0,-18.428d0,-18.412d0, &
                  -18.398d0,-18.386d0,-18.376d0,-18.367d0,-18.359d0 /)
      y(:,47) = (/-18.752d0,-18.639d0,-18.559d0,-18.501d0,-18.458d0, &
                  -18.426d0,-18.400d0,-18.380d0,-18.363d0,-18.350d0, &
                  -18.338d0,-18.328d0,-18.320d0,-18.313d0,-18.306d0 /)
      y(:,48) = (/-18.611d0,-18.523d0,-18.460d0,-18.415d0,-18.381d0, &
                  -18.355d0,-18.334d0,-18.318d0,-18.304d0,-18.293d0, &
                  -18.284d0,-18.276d0,-18.269d0,-18.263d0,-18.258d0 /)
      y(:,49) = (/-18.471d0,-18.408d0,-18.362d0,-18.329d0,-18.304d0, &
                  -18.285d0,-18.269d0,-18.257d0,-18.247d0,-18.238d0, &
                  -18.231d0,-18.224d0,-18.219d0,-18.214d0,-18.210d0 /)
      y(:,50) = (/-18.330d0,-18.290d0,-18.261d0,-18.239d0,-18.223d0, &
                  -18.211d0,-18.201d0,-18.192d0,-18.185d0,-18.179d0, &
                  -18.174d0,-18.169d0,-18.165d0,-18.162d0,-18.159d0 /)
      y(:,51) = (/-18.190d0,-18.168d0,-18.154d0,-18.143d0,-18.135d0, &
                  -18.129d0,-18.124d0,-18.120d0,-18.116d0,-18.112d0, &
                  -18.109d0,-18.106d0,-18.104d0,-18.102d0,-18.100d0 /)
      y(:,52) = (/-18.055d0,-18.047d0,-18.043d0,-18.042d0,-18.040d0, &
                  -18.039d0,-18.039d0,-18.038d0,-18.037d0,-18.036d0, &
                  -18.035d0,-18.034d0,-18.033d0,-18.033d0,-18.032d0 /)
      y(:,53) = (/-17.929d0,-17.931d0,-17.935d0,-17.939d0,-17.943d0, &
                  -17.946d0,-17.948d0,-17.950d0,-17.952d0,-17.953d0, &
                  -17.955d0,-17.956d0,-17.957d0,-17.958d0,-17.959d0 /)
      y(:,54) = (/-17.818d0,-17.826d0,-17.834d0,-17.842d0,-17.849d0, &
                  -17.855d0,-17.860d0,-17.865d0,-17.869d0,-17.872d0, &
                  -17.875d0,-17.878d0,-17.881d0,-17.883d0,-17.886d0 /)
      y(:,55) = (/-17.724d0,-17.736d0,-17.747d0,-17.758d0,-17.767d0, &
                  -17.775d0,-17.782d0,-17.788d0,-17.793d0,-17.798d0, &
                  -17.803d0,-17.807d0,-17.811d0,-17.815d0,-17.819d0 /)
      y(:,56) = (/-17.651d0,-17.665d0,-17.678d0,-17.690d0,-17.701d0, &
                  -17.710d0,-17.718d0,-17.725d0,-17.732d0,-17.738d0, &
                  -17.744d0,-17.749d0,-17.755d0,-17.760d0,-17.765d0 /)
      y(:,57) = (/-17.601d0,-17.615d0,-17.629d0,-17.642d0,-17.653d0, &
                  -17.663d0,-17.672d0,-17.680d0,-17.688d0,-17.695d0, &
                  -17.701d0,-17.708d0,-17.714d0,-17.720d0,-17.726d0 /)
      y(:,58) = (/-17.572d0,-17.587d0,-17.602d0,-17.614d0,-17.626d0, &
                  -17.636d0,-17.645d0,-17.654d0,-17.662d0,-17.670d0, &
                  -17.677d0,-17.684d0,-17.691d0,-17.698d0,-17.704d0 /)
      y(:,59) = (/-17.565d0,-17.581d0,-17.595d0,-17.607d0,-17.619d0, &
                  -17.629d0,-17.638d0,-17.647d0,-17.656d0,-17.664d0, &
                  -17.671d0,-17.679d0,-17.686d0,-17.693d0,-17.700d0 /)
      y(:,60) = (/-17.580d0,-17.594d0,-17.608d0,-17.620d0,-17.630d0, &
                  -17.640d0,-17.650d0,-17.658d0,-17.667d0,-17.675d0, &
                  -17.682d0,-17.690d0,-17.697d0,-17.704d0,-17.711d0 /)
      y(:,61) = (/-17.613d0,-17.626d0,-17.639d0,-17.649d0,-17.659d0, &
                  -17.669d0,-17.677d0,-17.686d0,-17.694d0,-17.701d0, &
                  -17.709d0,-17.716d0,-17.723d0,-17.730d0,-17.737d0 /)
      y(:,62) = (/-17.663d0,-17.675d0,-17.685d0,-17.695d0,-17.703d0, &
                  -17.711d0,-17.719d0,-17.727d0,-17.734d0,-17.741d0, &
                  -17.748d0,-17.755d0,-17.761d0,-17.768d0,-17.774d0 /)
      y(:,63) = (/-17.728d0,-17.737d0,-17.745d0,-17.752d0,-17.759d0, &
                  -17.766d0,-17.772d0,-17.778d0,-17.785d0,-17.791d0, &
                  -17.797d0,-17.803d0,-17.808d0,-17.814d0,-17.820d0 /)
      y(:,64) = (/-17.803d0,-17.809d0,-17.814d0,-17.818d0,-17.823d0, &
                  -17.828d0,-17.832d0,-17.837d0,-17.842d0,-17.847d0, &
                  -17.852d0,-17.856d0,-17.861d0,-17.866d0,-17.871d0 /)
      y(:,65) = (/-17.884d0,-17.886d0,-17.888d0,-17.889d0,-17.891d0, &
                  -17.893d0,-17.896d0,-17.899d0,-17.902d0,-17.905d0, &
                  -17.908d0,-17.912d0,-17.915d0,-17.919d0,-17.922d0 /)
      y(:,66) = (/-17.966d0,-17.964d0,-17.961d0,-17.959d0,-17.958d0, &
                  -17.958d0,-17.958d0,-17.959d0,-17.960d0,-17.961d0, &
                  -17.963d0,-17.964d0,-17.966d0,-17.968d0,-17.970d0 /)
      y(:,67) = (/-18.040d0,-18.034d0,-18.028d0,-18.023d0,-18.019d0, &
                  -18.016d0,-18.013d0,-18.012d0,-18.010d0,-18.010d0, &
                  -18.009d0,-18.009d0,-18.009d0,-18.009d0,-18.010d0 /)
      y(:,68) = (/-18.096d0,-18.087d0,-18.078d0,-18.071d0,-18.065d0, &
                  -18.059d0,-18.055d0,-18.051d0,-18.047d0,-18.045d0, &
                  -18.042d0,-18.040d0,-18.039d0,-18.037d0,-18.036d0 /)
      y(:,69) = (/-18.125d0,-18.115d0,-18.105d0,-18.097d0,-18.089d0, &
                  -18.082d0,-18.076d0,-18.070d0,-18.065d0,-18.061d0, &
                  -18.057d0,-18.053d0,-18.051d0,-18.048d0,-18.046d0 /)
      y(:,70) = (/-18.120d0,-18.112d0,-18.103d0,-18.095d0,-18.087d0, &
                  -18.079d0,-18.072d0,-18.066d0,-18.060d0,-18.055d0, &
                  -18.050d0,-18.046d0,-18.042d0,-18.039d0,-18.036d0 /)
      y(:,71) = (/-18.083d0,-18.078d0,-18.071d0,-18.064d0,-18.057d0, &
                  -18.050d0,-18.044d0,-18.037d0,-18.032d0,-18.026d0, &
                  -18.022d0,-18.017d0,-18.014d0,-18.010d0,-18.007d0 /)
      y(:,72) = (/-18.025d0,-18.022d0,-18.017d0,-18.012d0,-18.006d0, &
                  -18.000d0,-17.994d0,-17.989d0,-17.984d0,-17.979d0, &
                  -17.975d0,-17.971d0,-17.968d0,-17.965d0,-17.963d0 /)
      y(:,73) = (/-17.957d0,-17.955d0,-17.952d0,-17.948d0,-17.943d0, &
                  -17.938d0,-17.934d0,-17.929d0,-17.925d0,-17.922d0, &
                  -17.918d0,-17.916d0,-17.913d0,-17.911d0,-17.910d0 /)
      y(:,74) = (/-17.890d0,-17.889d0,-17.886d0,-17.882d0,-17.879d0, &
                  -17.875d0,-17.871d0,-17.867d0,-17.864d0,-17.862d0, &
                  -17.860d0,-17.858d0,-17.857d0,-17.856d0,-17.855d0 /)
      y(:,75) = (/-17.831d0,-17.829d0,-17.826d0,-17.822d0,-17.819d0, &
                  -17.815d0,-17.812d0,-17.810d0,-17.807d0,-17.806d0, &
                  -17.804d0,-17.803d0,-17.803d0,-17.803d0,-17.803d0 /)
      y(:,76) = (/-17.786d0,-17.782d0,-17.777d0,-17.773d0,-17.769d0, &
                  -17.766d0,-17.763d0,-17.761d0,-17.759d0,-17.758d0, &
                  -17.757d0,-17.757d0,-17.757d0,-17.758d0,-17.759d0 /)
      y(:,77) = (/-17.753d0,-17.747d0,-17.741d0,-17.735d0,-17.731d0, &
                  -17.727d0,-17.724d0,-17.722d0,-17.721d0,-17.720d0, &
                  -17.720d0,-17.720d0,-17.721d0,-17.722d0,-17.724d0 /)
      y(:,78) = (/-17.733d0,-17.724d0,-17.716d0,-17.709d0,-17.703d0, &
                  -17.699d0,-17.696d0,-17.694d0,-17.693d0,-17.692d0, &
                  -17.692d0,-17.693d0,-17.694d0,-17.695d0,-17.697d0 /)
      y(:,79) = (/-17.723d0,-17.711d0,-17.700d0,-17.691d0,-17.685d0, &
                  -17.680d0,-17.676d0,-17.674d0,-17.673d0,-17.672d0, &
                  -17.673d0,-17.673d0,-17.675d0,-17.676d0,-17.678d0 /)
      y(:,80) = (/-17.718d0,-17.702d0,-17.689d0,-17.679d0,-17.672d0, &
                  -17.667d0,-17.663d0,-17.660d0,-17.659d0,-17.659d0, &
                  -17.659d0,-17.660d0,-17.661d0,-17.663d0,-17.665d0 /)
      y(:,81) = (/-17.713d0,-17.695d0,-17.681d0,-17.670d0,-17.662d0, &
                  -17.656d0,-17.653d0,-17.650d0,-17.649d0,-17.649d0, &
                  -17.649d0,-17.650d0,-17.651d0,-17.653d0,-17.655d0 /)
      y(:,82) = (/-17.705d0,-17.686d0,-17.671d0,-17.660d0,-17.652d0, &
                  -17.647d0,-17.643d0,-17.641d0,-17.640d0,-17.640d0, &
                  -17.640d0,-17.641d0,-17.643d0,-17.645d0,-17.647d0 /)
      y(:,83) = (/-17.690d0,-17.671d0,-17.657d0,-17.647d0,-17.640d0, &
                  -17.635d0,-17.632d0,-17.630d0,-17.630d0,-17.630d0, &
                  -17.631d0,-17.632d0,-17.634d0,-17.636d0,-17.639d0 /)
      y(:,84) = (/-17.667d0,-17.649d0,-17.637d0,-17.629d0,-17.623d0, &
                  -17.619d0,-17.618d0,-17.617d0,-17.617d0,-17.618d0, &
                  -17.619d0,-17.621d0,-17.623d0,-17.626d0,-17.628d0 /)
      y(:,85) = (/-17.635d0,-17.621d0,-17.611d0,-17.605d0,-17.601d0, &
                  -17.600d0,-17.599d0,-17.599d0,-17.601d0,-17.602d0, &
                  -17.604d0,-17.607d0,-17.609d0,-17.612d0,-17.615d0 /)
      y(:,86) = (/-17.596d0,-17.585d0,-17.579d0,-17.576d0,-17.575d0, &
                  -17.575d0,-17.576d0,-17.578d0,-17.580d0,-17.582d0, &
                  -17.585d0,-17.588d0,-17.591d0,-17.595d0,-17.598d0 /)
      y(:,87) = (/-17.550d0,-17.544d0,-17.542d0,-17.542d0,-17.544d0, &
                  -17.546d0,-17.548d0,-17.552d0,-17.555d0,-17.558d0, &
                  -17.562d0,-17.566d0,-17.570d0,-17.573d0,-17.577d0 /)
      y(:,88) = (/-17.501d0,-17.500d0,-17.501d0,-17.504d0,-17.508d0, &
                  -17.513d0,-17.517d0,-17.521d0,-17.526d0,-17.530d0, &
                  -17.535d0,-17.539d0,-17.544d0,-17.548d0,-17.553d0 /)
      y(:,89) = (/-17.449d0,-17.452d0,-17.457d0,-17.463d0,-17.470d0, &
                  -17.476d0,-17.482d0,-17.488d0,-17.493d0,-17.499d0, &
                  -17.504d0,-17.509d0,-17.514d0,-17.519d0,-17.524d0 /)
      y(:,90) = (/-17.396d0,-17.403d0,-17.412d0,-17.420d0,-17.429d0, &
                  -17.437d0,-17.444d0,-17.451d0,-17.458d0,-17.464d0, &
                  -17.470d0,-17.476d0,-17.481d0,-17.487d0,-17.492d0 /)
      y(:,91) = (/-17.344d0,-17.355d0,-17.366d0,-17.377d0,-17.387d0, &
                  -17.396d0,-17.405d0,-17.413d0,-17.420d0,-17.427d0, &
                  -17.434d0,-17.440d0,-17.446d0,-17.452d0,-17.458d0 /)
      y(:,92) = (/-17.295d0,-17.307d0,-17.321d0,-17.333d0,-17.345d0, &
                  -17.355d0,-17.365d0,-17.373d0,-17.382d0,-17.389d0, &
                  -17.397d0,-17.404d0,-17.410d0,-17.417d0,-17.423d0 /)
      y(:,93) = (/-17.249d0,-17.264d0,-17.278d0,-17.292d0,-17.304d0, &
                  -17.316d0,-17.326d0,-17.335d0,-17.344d0,-17.352d0, &
                  -17.360d0,-17.368d0,-17.375d0,-17.382d0,-17.389d0 /)
      y(:,94) = (/-17.209d0,-17.225d0,-17.241d0,-17.255d0,-17.268d0, &
                  -17.280d0,-17.291d0,-17.301d0,-17.310d0,-17.319d0, &
                  -17.327d0,-17.335d0,-17.343d0,-17.350d0,-17.357d0 /)
      y(:,95) = (/-17.177d0,-17.194d0,-17.210d0,-17.225d0,-17.239d0, &
                  -17.251d0,-17.262d0,-17.272d0,-17.282d0,-17.291d0, &
                  -17.300d0,-17.308d0,-17.316d0,-17.324d0,-17.331d0 /)
      y(:,96) = (/-17.154d0,-17.172d0,-17.189d0,-17.204d0,-17.218d0, &
                  -17.230d0,-17.242d0,-17.252d0,-17.262d0,-17.272d0, &
                  -17.280d0,-17.289d0,-17.298d0,-17.306d0,-17.314d0 /)
      y(:,97) = (/-17.144d0,-17.162d0,-17.179d0,-17.194d0,-17.208d0, &
                  -17.220d0,-17.232d0,-17.242d0,-17.253d0,-17.262d0, &
                  -17.271d0,-17.280d0,-17.289d0,-17.297d0,-17.306d0 /)
      y(:,98) = (/-17.146d0,-17.164d0,-17.181d0,-17.196d0,-17.210d0, &
                  -17.222d0,-17.234d0,-17.245d0,-17.255d0,-17.265d0, &
                  -17.274d0,-17.283d0,-17.292d0,-17.301d0,-17.309d0 /)
      y(:,99) = (/-17.163d0,-17.180d0,-17.197d0,-17.212d0,-17.225d0, &
                  -17.237d0,-17.249d0,-17.260d0,-17.270d0,-17.280d0, &
                  -17.289d0,-17.298d0,-17.307d0,-17.316d0,-17.325d0 /)
      y(:,100) =(/-17.193d0,-17.211d0,-17.227d0,-17.241d0,-17.254d0, &
                  -17.266d0,-17.277d0,-17.288d0,-17.298d0,-17.308d0, &
                  -17.317d0,-17.327d0,-17.336d0,-17.345d0,-17.353d0 /)
      y(:,101) =(/-17.239d0,-17.256d0,-17.271d0,-17.284d0,-17.297d0, &
                  -17.309d0,-17.320d0,-17.330d0,-17.340d0,-17.350d0, &
                  -17.359d0,-17.369d0,-17.378d0,-17.387d0,-17.395d0 /)
      y(:,102) =(/-17.299d0,-17.315d0,-17.329d0,-17.342d0,-17.354d0, &
                  -17.365d0,-17.376d0,-17.386d0,-17.396d0,-17.405d0, &
                  -17.415d0,-17.424d0,-17.433d0,-17.442d0,-17.451d0 /)
      y(:,103) =(/-17.373d0,-17.388d0,-17.402d0,-17.414d0,-17.425d0, &
                  -17.436d0,-17.446d0,-17.456d0,-17.466d0,-17.475d0, &
                  -17.484d0,-17.493d0,-17.502d0,-17.511d0,-17.520d0 /)
      y(:,104) =(/-17.462d0,-17.476d0,-17.489d0,-17.500d0,-17.511d0, &
                  -17.521d0,-17.531d0,-17.541d0,-17.550d0,-17.559d0, &
                  -17.569d0,-17.578d0,-17.587d0,-17.595d0,-17.604d0 /)
      y(:,105) =(/-17.567d0,-17.581d0,-17.592d0,-17.603d0,-17.613d0, &
                  -17.623d0,-17.632d0,-17.641d0,-17.651d0,-17.660d0, &
                  -17.669d0,-17.678d0,-17.686d0,-17.695d0,-17.704d0 /)
      y(:,106) =(/-17.689d0,-17.701d0,-17.712d0,-17.722d0,-17.732d0, &
                  -17.741d0,-17.750d0,-17.759d0,-17.768d0,-17.777d0, &
                  -17.786d0,-17.795d0,-17.803d0,-17.812d0,-17.821d0 /)
      y(:,107) =(/-17.829d0,-17.840d0,-17.851d0,-17.860d0,-17.869d0, &
                  -17.878d0,-17.887d0,-17.896d0,-17.904d0,-17.913d0, &
                  -17.922d0,-17.930d0,-17.939d0,-17.948d0,-17.956d0 /)
      y(:,108) =(/-17.988d0,-18.000d0,-18.010d0,-18.019d0,-18.028d0, &
                  -18.036d0,-18.045d0,-18.053d0,-18.062d0,-18.070d0, &
                  -18.079d0,-18.087d0,-18.096d0,-18.104d0,-18.112d0 /)
      y(:,109) =(/-18.171d0,-18.183d0,-18.192d0,-18.201d0,-18.210d0, &
                  -18.218d0,-18.227d0,-18.235d0,-18.243d0,-18.252d0, &
                  -18.260d0,-18.268d0,-18.277d0,-18.285d0,-18.293d0 /)
      y(:,110) =(/-18.381d0,-18.393d0,-18.403d0,-18.413d0,-18.422d0, &
                  -18.430d0,-18.438d0,-18.447d0,-18.455d0,-18.463d0, &
                  -18.471d0,-18.479d0,-18.487d0,-18.495d0,-18.503d0 /)
      y(:,111) =(/-18.625d0,-18.638d0,-18.650d0,-18.660d0,-18.669d0, &
                  -18.678d0,-18.687d0,-18.695d0,-18.703d0,-18.711d0, &
                  -18.719d0,-18.726d0,-18.734d0,-18.742d0,-18.750d0 /)
      y(:,112) =(/-18.912d0,-18.929d0,-18.943d0,-18.955d0,-18.966d0, &
                  -18.975d0,-18.984d0,-18.993d0,-19.001d0,-19.008d0, &
                  -19.016d0,-19.023d0,-19.031d0,-19.038d0,-19.045d0 /)
      y(:,113) =(/-19.260d0,-19.283d0,-19.303d0,-19.320d0,-19.333d0, &
                  -19.345d0,-19.355d0,-19.364d0,-19.372d0,-19.380d0, &
                  -19.387d0,-19.394d0,-19.400d0,-19.407d0,-19.413d0 /)
      y(:,114) =(/-19.704d0,-19.740d0,-19.771d0,-19.796d0,-19.816d0, &
                  -19.832d0,-19.845d0,-19.855d0,-19.863d0,-19.870d0, &
                  -19.876d0,-19.882d0,-19.887d0,-19.892d0,-19.897d0 /)
      y(:,115) =(/-20.339d0,-20.386d0,-20.424d0,-20.454d0,-20.476d0, &
                  -20.492d0,-20.502d0,-20.509d0,-20.513d0,-20.516d0, &
                  -20.518d0,-20.520d0,-20.521d0,-20.523d0,-20.524d0 /)
      y(:,116) =(/-21.052d0,-21.075d0,-21.093d0,-21.105d0,-21.114d0, &
                  -21.120d0,-21.123d0,-21.125d0,-21.126d0,-21.127d0, &
                  -21.128d0,-21.130d0,-21.131d0,-21.133d0,-21.135d0 /)
      y(:,117) =(/-21.174d0,-21.203d0,-21.230d0,-21.255d0,-21.278d0, &
                  -21.299d0,-21.320d0,-21.339d0,-21.357d0,-21.375d0, &
                  -21.392d0,-21.408d0,-21.424d0,-21.439d0,-21.454d0 /)
      y(:,118) =(/-21.285d0,-21.317d0,-21.346d0,-21.372d0,-21.395d0, &
                  -21.416d0,-21.435d0,-21.452d0,-21.468d0,-21.483d0, &
                  -21.497d0,-21.511d0,-21.524d0,-21.536d0,-21.548d0 /)
      y(:,119) =(/-21.396d0,-21.429d0,-21.459d0,-21.486d0,-21.511d0, &
                  -21.532d0,-21.551d0,-21.569d0,-21.585d0,-21.600d0, &
                  -21.614d0,-21.627d0,-21.640d0,-21.652d0,-21.663d0 /)
      y(:,120) =(/-21.516d0,-21.549d0,-21.580d0,-21.609d0,-21.635d0, &
                  -21.658d0,-21.678d0,-21.696d0,-21.713d0,-21.728d0, &
                  -21.742d0,-21.755d0,-21.767d0,-21.779d0,-21.790d0 /)
      y(:,121) =(/-21.651d0,-21.681d0,-21.711d0,-21.738d0,-21.763d0, &
                  -21.785d0,-21.804d0,-21.821d0,-21.837d0,-21.851d0, &
                  -21.864d0,-21.876d0,-21.887d0,-21.898d0,-21.908d0 /)
      y(:,122) =(/-21.810d0,-21.831d0,-21.853d0,-21.874d0,-21.893d0, &
                  -21.910d0,-21.925d0,-21.938d0,-21.950d0,-21.961d0, &
                  -21.971d0,-21.980d0,-21.989d0,-21.998d0,-22.006d0 /)
      y(:,123) =(/-22.009d0,-22.016d0,-22.026d0,-22.037d0,-22.048d0, &
                  -22.058d0,-22.066d0,-22.074d0,-22.081d0,-22.088d0, &
                  -22.094d0,-22.099d0,-22.105d0,-22.111d0,-22.117d0 /)
      y(:,124) =(/-22.353d0,-22.317d0,-22.296d0,-22.284d0,-22.276d0, &
                  -22.270d0,-22.266d0,-22.262d0,-22.260d0,-22.258d0, &
                  -22.257d0,-22.257d0,-22.257d0,-22.258d0,-22.259d0 /)
      y(:,125) =(/-22.705d0,-22.609d0,-22.552d0,-22.515d0,-22.488d0, &
                  -22.468d0,-22.451d0,-22.438d0,-22.427d0,-22.418d0, &
                  -22.410d0,-22.405d0,-22.400d0,-22.397d0,-22.395d0 /)
      y(:,126) =(/-22.889d0,-22.791d0,-22.731d0,-22.690d0,-22.659d0, &
                  -22.634d0,-22.612d0,-22.594d0,-22.579d0,-22.566d0, &
                  -22.555d0,-22.546d0,-22.539d0,-22.533d0,-22.528d0 /)
      y(:,127) =(/-23.211d0,-23.109d0,-23.041d0,-22.989d0,-22.945d0, &
                  -22.906d0,-22.872d0,-22.842d0,-22.816d0,-22.793d0, &
                  -22.774d0,-22.757d0,-22.743d0,-22.732d0,-22.722d0 /)
      y(:,128) =(/-25.312d0,-24.669d0,-24.250d0,-23.959d0,-23.746d0, &
                  -23.587d0,-23.463d0,-23.366d0,-23.288d0,-23.225d0, &
                  -23.173d0,-23.131d0,-23.095d0,-23.066d0,-23.041d0 /)
      y(:,129) =(/-25.394d0,-24.752d0,-24.333d0,-24.041d0,-23.829d0, &
                  -23.669d0,-23.546d0,-23.449d0,-23.371d0,-23.308d0, &
                  -23.256d0,-23.214d0,-23.178d0,-23.149d0,-23.124d0 /)
      y(:,130) =(/-25.430d0,-24.787d0,-24.369d0,-24.077d0,-23.865d0, &
                  -23.705d0,-23.582d0,-23.484d0,-23.407d0,-23.344d0, &
                  -23.292d0,-23.249d0,-23.214d0,-23.185d0,-23.160d0 /)

      ! Frequency quantities
      pEeV = convF*freq
      pE = pEeV*2d21*c*freq*freq
      pEeV = pEeV*ergtoev
      arg = c2*freq*1d4

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! If out of limits, no contribution
      if (pEeV.lt.x1(1).or.pEeV.gt.x1(n1)) then

        ! Return with zeros
        return

      ! Within limits in axis 1
      else

        ! For each height
        do iz=iz0,iz1

          ! Within temperature limits
          if (T(iz).ge.x2(1).and.T(iz).le.x2(n2)) then

            ! Bilinear interpolation of cross-section
            call bilinear(x2,x1,y,T(iz),pEeV,cc)

            ! It was logarithmic
            cc = 1d1**cc

            ! Compute exponential
            exu = arg/T(iz)
            exu = diexp(exu)

            ! Compute absorptivity and emissivity
            eta(iz) = nOH(iz)*(1d0 - exu)*cc
            eps(iz) = nOH(iz)*pE*exu*cc

          end if ! Within temperature limits

        end do ! Heights

      end if ! Within limits in axis 1

      end subroutine OH_bf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes CH bound-free absorptivity and emissivity\n
      !!    freq(double): Frequency\n
      !!  nCH(double(:)): CH density\n
      !!    T(double(:)): Temperature\n
      !!    iz0(integer): First height index to consider\n
      !!    iz1(integer): Last height index to consider\n
      !!  eta(double(:)): Absorptivity\n
      !!  eps(double(:)): Emissivity
      subroutine CH_bf(freq,nCH,T,iz0,iz1,eta,eps)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: nCH, T
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      integer, parameter:: n2=15, n1=105
      integer:: iz

      double precision:: arg, pE, pEeV, cc, exu
      double precision, dimension(n1):: x1
      double precision, dimension(n2):: x2
      double precision, dimension(n2,n1):: y


      ! Kurucz, van Dishoeck, and Tarafdar (1987) table 2
      x1 = (/ 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, &
              0.8d0, 0.9d0, 1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0, &
              1.5d0, 1.6d0, 1.7d0, 1.8d0, 1.9d0, 2.0d0, 2.1d0, &
              2.2d0, 2.3d0, 2.4d0, 2.5d0, 2.6d0, 2.7d0, 2.8d0, &
              2.9d0, 3.0d0, 3.1d0, 3.2d0, 3.3d0, 3.4d0, 3.5d0, &
              3.6d0, 3.7d0, 3.8d0, 3.9d0, 4.0d0, 4.1d0, 4.2d0, &
              4.3d0, 4.4d0, 4.5d0, 4.6d0, 4.7d0, 4.8d0, 4.9d0, &
              5.0d0, 5.1d0, 5.2d0, 5.3d0, 5.4d0, 5.5d0, 5.6d0, &
              5.7d0, 5.8d0, 5.9d0, 6.0d0, 6.1d0, 6.2d0, 6.3d0, &
              6.4d0, 6.5d0, 6.6d0, 6.7d0, 6.8d0, 6.9d0, 7.0d0, &
              7.1d0, 7.2d0, 7.3d0, 7.4d0, 7.5d0, 7.6d0, 7.7d0, &
              7.8d0, 7.9d0, 8.0d0, 8.1d0, 8.2d0, 8.3d0, 8.4d0, &
              8.5d0, 8.6d0, 8.7d0, 8.8d0, 8.9d0, 9.0d0, 9.1d0, &
              9.2d0, 9.3d0, 9.4d0, 9.5d0, 9.6d0, 9.7d0, 9.8d0, &
              9.9d0, 10.0d0, 10.1d0, 10.2d0, 10.3d0, 10.4d0, 10.5d0 /)
      x2 = (/ 2d3, 2.5d3, 3d3, 3.5d3, 4d3, 4.5d3, 5d3, 5.5d3, 6d3, &
              6.5d3, 7d3, 7.5d3, 8d3, 8.5d3, 9d3 /)
      y(:,1) = (/ -38d0,-38d0,-38d0,-38d0,-38d0,-38d0,-38d0,-38d0, &
                  -38d0,-38d0,-38d0,-38d0,-38d0,-38d0,-38d0 /)
      y(:,2) = (/ -32.727d0, -31.151d0, -30.133d0, -29.432d0, &
                  -28.925d0, -28.547d0, -28.257d0, -28.030d0, &
                  -27.848d0, -27.701d0, -27.580d0, -27.479d0, &
                  -27.395d0, -27.322d0, -27.261d0 /)
      y(:,3) = (/ -31.588d0, -30.011d0, -28.993d0, -28.290d0, &
                  -27.784d0, -27.405d0, -27.115d0, -26.887d0, &
                  -26.705d0, -26.558d0, -26.437d0, -26.336d0, &
                  -26.251d0, -26.179d0, -26.117d0 /)
      y(:,4) = (/ -30.407d0, -28.830d0, -27.811d0, -27.108d0, &
                  -26.601d0, -26.223d0, -25.932d0, -25.705d0, &
                  -25.523d0, -25.376d0, -25.255d0, -25.154d0, &
                  -25.069d0, -24.997d0, -24.935d0 /)
      y(:,5) = (/ -29.513d0, -27.937d0, -26.920d0, -26.218d0, &
                  -25.712d0, -25.334d0, -25.043d0, -24.816d0, &
                  -24.635d0, -24.487d0, -24.366d0, -24.266d0, &
                  -24.181d0, -24.109d0, -24.047d0 /)
      y(:,6) = (/ -28.910d0, -27.341d0, -26.327d0, -25.628d0, &
                  -25.123d0, -24.746d0, -24.457d0, -24.230d0, &
                  -24.049d0, -23.902d0, -23.782d0, -23.681d0, &
                  -23.597d0, -23.525d0, -23.464d0 /)
      y(:,7) = (/ -28.517d0, -26.961d0, -25.955d0, -25.261d0, &
                  -24.760d0, -24.385d0, -24.098d0, -23.873d0, &
                  -23.694d0, -23.548d0, -23.429d0, -23.329d0, &
                  -23.245d0, -23.174d0, -23.113d0 /)
      y(:,8) = (/ -28.213d0, -26.675d0, -25.680d0, -24.993d0, &
                  -24.497d0, -24.127d0, -23.843d0, -23.620d0, &
                  -23.443d0, -23.299d0, -23.181d0, -23.082d0, &
                  -22.999d0, -22.929d0, -22.869d0 /)
      y(:,9) = (/ -27.942d0, -26.427d0, -25.446d0, -24.769d0, &
                  -24.280d0, -23.915d0, -23.635d0, -23.416d0, &
                  -23.241d0, -23.100d0, -22.983d0, -22.887d0, &
                  -22.805d0, -22.736d0, -22.677d0 /)
      y(:,10) = (/ -27.706d0, -26.210d0, -25.241d0, -24.572d0, &
                   -24.088d0, -23.728d0, -23.451d0, -23.235d0, &
                   -23.063d0, -22.923d0, -22.808d0, -22.713d0, &
                   -22.633d0, -22.565d0, -22.507d0 /)
      y(:,11) = (/ -27.475d0, -26.000d0, -25.043d0, -24.382d0, &
                   -23.905d0, -23.548d0, -23.275d0, -23.062d0, &
                   -22.891d0, -22.753d0, -22.640d0, -22.546d0, &
                   -22.467d0, -22.400d0, -22.343d0 /)
      y(:,12) = (/ -27.221d0, -25.783d0, -24.844d0, -24.193d0, &
                   -23.723d0, -23.372d0, -23.102d0, -22.892d0, &
                   -22.724d0, -22.588d0, -22.476d0, -22.384d0, &
                   -22.306d0, -22.240d0, -22.184d0 /)
      y(:,13) = (/ -26.863d0, -25.506d0, -24.607d0, -23.979d0, &
                   -23.523d0, -23.182d0, -22.919d0, -22.714d0, &
                   -22.550d0, -22.417d0, -22.309d0, -22.218d0, &
                   -22.142d0, -22.078d0, -22.023d0 /)
      y(:,14) = (/ -26.685d0, -25.347d0, -24.457d0, -23.835d0, &
                   -23.382d0, -23.044d0, -22.784d0, -22.580d0, &
                   -22.418d0, -22.286d0, -22.178d0, -22.089d0, &
                   -22.014d0, -21.950d0, -21.896d0 /)
      y(:,15) = (/ -26.085d0, -24.903d0, -24.105d0, -23.538d0, &
                   -23.120d0, -22.805d0, -22.561d0, -22.370d0, &
                   -22.217d0, -22.093d0, -21.991d0, -21.906d0, &
                   -21.835d0, -21.775d0, -21.723d0 /)
      y(:,16) = (/ -25.902d0, -24.727d0, -23.936d0, -23.376d0, &
                   -22.964d0, -22.654d0, -22.415d0, -22.227d0, &
                   -22.076d0, -21.955d0, -21.855d0, -21.772d0, &
                   -21.702d0, -21.644d0, -21.593d0 /)
      y(:,17) = (/ -25.215d0, -24.196d0, -23.510d0, -23.019d0, &
                   -22.655d0, -22.378d0, -22.163d0, -21.992d0, &
                   -21.855d0, -21.744d0, -21.653d0, -21.577d0, &
                   -21.513d0, -21.459d0, -21.412d0 /)
      y(:,18) = (/ -24.914d0, -23.937d0, -23.284d0, -22.820d0, &
                   -22.475d0, -22.212d0, -22.007d0, -21.845d0, &
                   -21.715d0, -21.609d0, -21.522d0, -21.449d0, &
                   -21.388d0, -21.336d0, -21.292d0 /)
      y(:,19) = (/ -24.519d0, -23.637d0, -23.039d0, -22.606d0, &
                   -22.281d0, -22.030d0, -21.834d0, -21.678d0, &
                   -21.552d0, -21.450d0, -21.365d0, -21.295d0, &
                   -21.236d0, -21.185d0, -21.142d0 /)
      y(:,20) = (/ -24.086d0, -23.222d0, -22.650d0, -22.246d0, &
                   -21.948d0, -21.722d0, -21.546d0, -21.407d0, &
                   -21.296d0, -21.205d0, -21.131d0, -21.070d0, &
                   -21.018d0, -20.974d0, -20.937d0 /)
      y(:,21) = (/ -23.850d0, -23.018d0, -22.472d0, -22.088d0, &
                   -21.805d0, -21.590d0, -21.422d0, -21.289d0, &
                   -21.182d0, -21.095d0, -21.024d0, -20.964d0, &
                   -20.914d0, -20.872d0, -20.835d0 /)
      y(:,22) = (/ -23.136d0, -22.445d0, -21.994d0, -21.676d0, &
                   -21.440d0, -21.259d0, -21.117d0, -21.004d0, &
                   -20.912d0, -20.837d0, -20.775d0, -20.723d0, &
                   -20.679d0, -20.642d0, -20.611d0 /)
      y(:,23) = (/ -23.199d0, -22.433d0, -21.927d0, -21.573d0, &
                   -21.314d0, -21.119d0, -20.969d0, -20.851d0, &
                   -20.758d0, -20.682d0, -20.621d0, -20.571d0, &
                   -20.529d0, -20.493d0, -20.463d0 /)
      y(:,24) = (/ -22.696d0, -22.020d0, -21.585d0, -21.286d0, &
                   -21.071d0, -20.912d0, -20.791d0, -20.697d0, &
                   -20.622d0, -20.563d0, -20.514d0, -20.475d0, &
                   -20.442d0, -20.414d0, -20.391d0 /)
      y(:,25) = (/ -22.119d0, -21.557d0, -21.194d0, -20.943d0, &
                   -20.761d0, -20.624d0, -20.518d0, -20.434d0, &
                   -20.367d0, -20.313d0, -20.268d0, -20.231d0, &
                   -20.201d0, -20.175d0, -20.153d0 /)
      y(:,26) = (/ -21.855d0, -21.300d0, -20.931d0, -20.673d0, &
                   -20.485d0, -20.344d0, -20.235d0, -20.151d0, &
                   -20.084d0, -20.031d0, -19.988d0, -19.953d0, &
                   -19.924d0, -19.900d0, -19.880d0 /)
      y(:,27) = (/ -21.126d0, -20.673d0, -20.382d0, -20.184d0, &
                   -20.044d0, -19.943d0, -19.868d0, -19.811d0, &
                   -19.769d0, -19.736d0, -19.710d0, -19.690d0, &
                   -19.674d0, -19.662d0, -19.652d0 /)
      y(:,28) = (/ -20.502d0, -20.150d0, -19.922d0, -19.766d0, &
                   -19.657d0, -19.578d0, -19.520d0, -19.478d0, &
                   -19.446d0, -19.422d0, -19.404d0, -19.390d0, &
                   -19.379d0, -19.371d0, -19.365d0 /)
      y(:,29) = (/ -20.030d0, -19.724d0, -19.530d0, -19.399d0, &
                   -19.309d0, -19.245d0, -19.199d0, -19.166d0, &
                   -19.142d0, -19.125d0, -19.112d0, -19.103d0, &
                   -19.096d0, -19.091d0, -19.088d0 /)
      y(:,30) = (/ -19.640d0, -19.364d0, -19.189d0, -19.074d0, &
                   -18.996d0, -18.943d0, -18.906d0, -18.881d0, &
                   -18.863d0, -18.852d0, -18.844d0, -18.839d0, &
                   -18.837d0, -18.836d0, -18.836d0 /)
      y(:,31) = (/ -19.333d0, -19.092d0, -18.939d0, -18.838d0, &
                   -18.770d0, -18.725d0, -18.695d0, -18.675d0, &
                   -18.662d0, -18.655d0, -18.651d0, -18.649d0, &
                   -18.649d0, -18.651d0, -18.653d0 /)
      y(:,32) = (/ -19.070d0, -18.880d0, -18.756d0, -18.674d0, &
                   -18.621d0, -18.585d0, -18.562d0, -18.548d0, &
                   -18.540d0, -18.536d0, -18.536d0, -18.537d0, &
                   -18.539d0, -18.542d0, -18.546d0 /)
      y(:,33) = (/ -18.851d0, -18.708d0, -18.617d0, -18.558d0, &
                   -18.521d0, -18.498d0, -18.484d0, -18.477d0, &
                   -18.475d0, -18.476d0, -18.478d0, -18.482d0, &
                   -18.487d0, -18.493d0, -18.498d0 /)
      y(:,34) = (/ -18.709d0, -18.599d0, -18.533d0, -18.494d0, &
                   -18.471d0, -18.459d0, -18.454d0, -18.454d0, &
                   -18.457d0, -18.462d0, -18.469d0, -18.476d0, &
                   -18.483d0, -18.490d0, -18.498d0 /)
      y(:,35) = (/ -18.656d0, -18.572d0, -18.524d0, -18.497d0, &
                   -18.485d0, -18.480d0, -18.482d0, -18.486d0, &
                   -18.493d0, -18.501d0, -18.510d0, -18.519d0, &
                   -18.527d0, -18.536d0, -18.544d0 /)
      y(:,36) = (/ -18.670d0, -18.613d0, -18.582d0, -18.566d0, &
                   -18.561d0, -18.562d0, -18.568d0, -18.575d0, &
                   -18.583d0, -18.592d0, -18.601d0, -18.610d0, &
                   -18.619d0, -18.627d0, -18.635d0 /)
      y(:,37) = (/ -18.728d0, -18.700d0, -18.687d0, -18.683d0, &
                   -18.685d0, -18.691d0, -18.698d0, -18.706d0, &
                   -18.715d0, -18.723d0, -18.731d0, -18.739d0, &
                   -18.745d0, -18.752d0, -18.758d0 /)
      y(:,38) = (/ -18.839d0, -18.835d0, -18.836d0, -18.842d0, &
                   -18.849d0, -18.857d0, -18.865d0, -18.872d0, &
                   -18.878d0, -18.883d0, -18.888d0, -18.892d0, &
                   -18.895d0, -18.898d0, -18.900d0 /)
      y(:,39) = (/ -19.034d0, -19.041d0, -19.049d0, -19.057d0, &
                   -19.064d0, -19.069d0, -19.071d0, -19.071d0, &
                   -19.070d0, -19.068d0, -19.065d0, -19.061d0, &
                   -19.058d0, -19.054d0, -19.051d0 /)
      y(:,40) = (/ -19.372d0, -19.378d0, -19.382d0, -19.380d0, &
                   -19.372d0, -19.359d0, -19.341d0, -19.321d0, &
                   -19.300d0, -19.280d0, -19.261d0, -19.243d0, &
                   -19.227d0, -19.212d0, -19.199d0 /)
      y(:,41) = (/ -19.780d0, -19.777d0, -19.763d0, -19.732d0, &
                   -19.686d0, -19.631d0, -19.573d0, -19.517d0, &
                   -19.465d0, -19.419d0, -19.379d0, -19.344d0, &
                   -19.314d0, -19.288d0, -19.265d0 /)
      y(:,42) = (/ -20.151d0, -20.133d0, -20.087d0, -20.009d0, &
                   -19.911d0, -19.810d0, -19.715d0, -19.631d0, &
                   -19.559d0, -19.497d0, -19.446d0, -19.402d0, &
                   -19.365d0, -19.333d0, -19.306d0 /)
      y(:,43) = (/ -20.525d0, -20.454d0, -20.312d0, -20.138d0, &
                   -19.970d0, -19.825d0, -19.705d0, -19.607d0, &
                   -19.528d0, -19.464d0, -19.411d0, -19.367d0, &
                   -19.330d0, -19.300d0, -19.274d0 /)
      y(:,44) = (/ -20.869d0, -20.655d0, -20.366d0, -20.104d0, &
                   -19.894d0, -19.731d0, -19.604d0, -19.505d0, &
                   -19.426d0, -19.363d0, -19.312d0, -19.271d0, &
                   -19.236d0, -19.208d0, -19.184d0 /)
      y(:,45) = (/ -21.179d0, -20.768d0, -20.380d0, -20.081d0, &
                   -19.856d0, -19.686d0, -19.556d0, -19.454d0, &
                   -19.375d0, -19.311d0, -19.260d0, -19.218d0, &
                   -19.184d0, -19.155d0, -19.131d0 /)
      y(:,46) = (/ -21.167d0, -20.601d0, -20.206d0, -19.925d0, &
                   -19.719d0, -19.565d0, -19.447d0, -19.355d0, &
                   -19.283d0, -19.226d0, -19.180d0, -19.143d0, &
                   -19.112d0, -19.087d0, -19.066d0 /)
      y(:,47) = (/ -20.918d0, -20.348d0, -19.976d0, -19.720d0, &
                   -19.536d0, -19.401d0, -19.299d0, -19.220d0, &
                   -19.159d0, -19.112d0, -19.073d0, -19.043d0, &
                   -19.018d0, -18.998d0, -18.981d0 /)
      y(:,48) = (/ -20.753d0, -20.204d0, -19.847d0, -19.602d0, &
                   -19.427d0, -19.299d0, -19.203d0, -19.129d0, &
                   -19.072d0, -19.028d0, -18.993d0, -18.965d0, &
                   -18.942d0, -18.924d0, -18.909d0 /)
      y(:,49) = (/ -20.456d0, -19.987d0, -19.677d0, -19.460d0, &
                   -19.302d0, -19.186d0, -19.098d0, -19.030d0, &
                   -18.978d0, -18.937d0, -18.904d0, -18.878d0, &
                   -18.857d0, -18.841d0, -18.827d0 /)
      y(:,50) = (/ -20.154d0, -19.734d0, -19.461d0, -19.272d0, &
                   -19.136d0, -19.035d0, -18.960d0, -18.902d0, &
                   -18.858d0, -18.824d0, -18.797d0, -18.775d0, &
                   -18.759d0, -18.745d0, -18.735d0 /)
      y(:,51) = (/ -19.941d0, -19.544d0, -19.288d0, -19.114d0, &
                   -18.992d0, -18.903d0, -18.837d0, -18.788d0, &
                   -18.751d0, -18.723d0, -18.701d0, -18.684d0, &
                   -18.671d0, -18.661d0, -18.654d0 /)
      y(:,52) = (/ -19.657d0, -19.321d0, -19.104d0, -18.956d0, &
                   -18.853d0, -18.779d0, -18.724d0, -18.684d0, &
                   -18.655d0, -18.632d0, -18.615d0, -18.602d0, &
                   -18.592d0, -18.585d0, -18.579d0 /)
      y(:,53) = (/ -19.388d0, -19.109d0, -18.930d0, -18.810d0, &
                   -18.725d0, -18.664d0, -18.620d0, -18.586d0, &
                   -18.562d0, -18.543d0, -18.529d0, -18.518d0, &
                   -18.510d0, -18.503d0, -18.498d0 /)
      y(:,54) = (/ -19.201d0, -18.953d0, -18.794d0, -18.686d0, &
                   -18.611d0, -18.556d0, -18.515d0, -18.485d0, &
                   -18.462d0, -18.446d0, -18.433d0, -18.423d0, &
                   -18.416d0, -18.410d0, -18.406d0 /)
      y(:,55) = (/ -18.923d0, -18.719d0, -18.588d0, -18.500d0, &
                   -18.439d0, -18.396d0, -18.365d0, -18.344d0, &
                   -18.328d0, -18.318d0, -18.311d0, -18.307d0, &
                   -18.304d0, -18.303d0, -18.302d0 /)
      y(:,56) = (/ -18.614d0, -18.458d0, -18.361d0, -18.298d0, &
                   -18.258d0, -18.232d0, -18.216d0, -18.206d0, &
                   -18.202d0, -18.201d0, -18.202d0, -18.205d0, &
                   -18.208d0, -18.213d0, -18.218d0 /)
      y(:,57) = (/ -18.419d0, -18.295d0, -18.222d0, -18.178d0, &
                   -18.153d0, -18.139d0, -18.132d0, -18.131d0, &
                   -18.133d0, -18.138d0, -18.143d0, -18.150d0, &
                   -18.157d0, -18.164d0, -18.172d0 /)
      y(:,58) = (/ -18.296d0, -18.201d0, -18.148d0, -18.118d0, &
                   -18.101d0, -18.094d0, -18.091d0, -18.093d0, &
                   -18.096d0, -18.101d0, -18.107d0, -18.113d0, &
                   -18.120d0, -18.126d0, -18.132d0 /)
      y(:,59) = (/ -18.021d0, -17.992d0, -17.977d0, -17.970d0, &
                   -17.967d0, -17.968d0, -17.970d0, -17.974d0, &
                   -17.978d0, -17.983d0, -17.989d0, -17.994d0, &
                   -18.000d0, -18.005d0, -18.011d0 /)
      y(:,60) = (/ -17.694d0, -17.686d0, -17.686d0, -17.691d0, &
                   -17.698d0, -17.708d0, -17.718d0, -17.729d0, &
                   -17.740d0, -17.750d0, -17.761d0, -17.771d0, &
                   -17.781d0, -17.790d0, -17.798d0 /)
      y(:,61) = (/ -17.374d0, -17.384d0, -17.400d0, -17.420d0, &
                   -17.440d0, -17.462d0, -17.483d0, -17.503d0, &
                   -17.523d0, -17.541d0, -17.558d0, -17.575d0, &
                   -17.590d0, -17.603d0, -17.616d0 /)
      y(:,62) = (/ -17.169d0, -17.199d0, -17.230d0, -17.262d0, &
                   -17.293d0, -17.323d0, -17.351d0, -17.378d0, &
                   -17.404d0, -17.427d0, -17.449d0, -17.469d0, &
                   -17.488d0, -17.505d0, -17.520d0 /)
      y(:,63) = (/ -17.151d0, -17.184d0, -17.217d0, -17.250d0, &
                   -17.282d0, -17.313d0, -17.342d0, -17.369d0, &
                   -17.395d0, -17.418d0, -17.440d0, -17.461d0, &
                   -17.480d0, -17.497d0, -17.513d0 /)
      y(:,64) = (/ -17.230d0, -17.260d0, -17.290d0, -17.320d0, &
                   -17.348d0, -17.375d0, -17.401d0, -17.425d0, &
                   -17.448d0, -17.469d0, -17.489d0, -17.508d0, &
                   -17.525d0, -17.541d0, -17.556d0 /)
      y(:,65) = (/ -17.379d0, -17.403d0, -17.425d0, -17.446d0, &
                   -17.467d0, -17.486d0, -17.505d0, -17.524d0, &
                   -17.541d0, -17.558d0, -17.574d0, -17.588d0, &
                   -17.602d0, -17.615d0, -17.627d0 /)
      y(:,66) = (/ -17.596d0, -17.604d0, -17.609d0, -17.612d0, &
                   -17.616d0, -17.622d0, -17.628d0, -17.636d0, &
                   -17.644d0, -17.652d0, -17.661d0, -17.670d0, &
                   -17.679d0, -17.687d0, -17.695d0 /)
      y(:,67) = (/ -17.846d0, -17.823d0, -17.795d0, -17.770d0, &
                   -17.750d0, -17.735d0, -17.725d0, -17.719d0, &
                   -17.716d0, -17.715d0, -17.716d0, -17.719d0, &
                   -17.722d0, -17.726d0, -17.730d0 /)
      y(:,68) = (/ -18.089d0, -18.015d0, -17.942d0, -17.882d0, &
                   -17.836d0, -17.802d0, -17.777d0, -17.760d0, &
                   -17.748d0, -17.740d0, -17.736d0, -17.734d0, &
                   -17.733d0, -17.734d0, -17.736d0 /)
      y(:,69) = (/ -18.299d0, -18.156d0, -18.038d0, -17.947d0, &
                   -17.881d0, -17.833d0, -17.798d0, -17.774d0, &
                   -17.757d0, -17.745d0, -17.738d0, -17.733d0, &
                   -17.730d0, -17.729d0, -17.729d0 /)
      y(:,70) = (/ -18.441d0, -18.243d0, -18.096d0, -17.991d0, &
                   -17.915d0, -17.860d0, -17.821d0, -17.792d0, &
                   -17.772d0, -17.757d0, -17.746d0, -17.738d0, &
                   -17.733d0, -17.730d0, -17.728d0 /)
      y(:,71) = (/ -18.474d0, -18.262d0, -18.111d0, -18.004d0, &
                   -17.926d0, -17.869d0, -17.826d0, -17.795d0, &
                   -17.771d0, -17.753d0, -17.740d0, -17.730d0, &
                   -17.722d0, -17.717d0, -17.713d0 /)
      y(:,72) = (/ -18.387d0, -18.191d0, -18.053d0, -17.952d0, &
                   -17.878d0, -17.823d0, -17.782d0, -17.752d0, &
                   -17.729d0, -17.711d0, -17.698d0, -17.689d0, &
                   -17.681d0, -17.676d0, -17.672d0 /)
      y(:,73) = (/ -18.161d0, -17.990d0, -17.874d0, -17.793d0, &
                   -17.736d0, -17.696d0, -17.668d0, -17.648d0, &
                   -17.634d0, -17.625d0, -17.619d0, -17.616d0, &
                   -17.614d0, -17.614d0, -17.615d0 /)
      y(:,74) = (/ -17.908d0, -17.774d0, -17.690d0, -17.637d0, &
                   -17.604d0, -17.583d0, -17.572d0, -17.567d0, &
                   -17.566d0, -17.568d0, -17.571d0, -17.576d0, &
                   -17.581d0, -17.587d0, -17.593d0 /)
      y(:,75) = (/ -17.681d0, -17.589d0, -17.540d0, -17.515d0, &
                   -17.506d0, -17.505d0, -17.511d0, -17.520d0, &
                   -17.530d0, -17.542d0, -17.554d0, -17.566d0, &
                   -17.578d0, -17.589d0, -17.600d0 /)
      y(:,76) = (/ -17.647d0, -17.606d0, -17.584d0, -17.575d0, &
                   -17.573d0, -17.576d0, -17.582d0, -17.589d0, &
                   -17.597d0, -17.605d0, -17.614d0, -17.623d0, &
                   -17.631d0, -17.639d0, -17.646d0 /)
      y(:,77) = (/ -17.300d0, -17.291d0, -17.291d0, -17.297d0, &
                   -17.307d0, -17.319d0, -17.333d0, -17.347d0, &
                   -17.361d0, -17.375d0, -17.389d0, -17.402d0, &
                   -17.415d0, -17.427d0, -17.438d0 /)
      y(:,78) = (/ -16.786d0, -16.802d0, -16.825d0, -16.853d0, &
                   -16.883d0, -16.914d0, -16.944d0, -16.974d0, &
                   -17.003d0, -17.030d0, -17.055d0, -17.079d0, &
                   -17.101d0, -17.122d0, -17.141d0 /)
      y(:,79) = (/ -16.489d0, -16.533d0, -16.579d0, -16.625d0, &
                   -16.670d0, -16.713d0, -16.754d0, -16.793d0, &
                   -16.830d0, -16.864d0, -16.896d0, -16.925d0, &
                   -16.952d0, -16.977d0, -17.000d0 /)
      y(:,80) = (/ -16.694d0, -16.724d0, -16.756d0, -16.789d0, &
                   -16.823d0, -16.856d0, -16.888d0, -16.919d0, &
                   -16.949d0, -16.976d0, -17.002d0, -17.026d0, &
                   -17.048d0, -17.069d0, -17.088d0 /)
      y(:,81) = (/ -16.935d0, -16.951d0, -16.971d0, -16.993d0, &
                   -17.016d0, -17.040d0, -17.064d0, -17.088d0, &
                   -17.111d0, -17.132d0, -17.153d0, -17.172d0, &
                   -17.190d0, -17.206d0, -17.222d0 /)
      y(:,82) = (/ -17.200d0, -17.208d0, -17.220d0, -17.235d0, &
                   -17.251d0, -17.269d0, -17.286d0, -17.304d0, &
                   -17.322d0, -17.338d0, -17.354d0, -17.369d0, &
                   -17.384d0, -17.397d0, -17.409d0 /)
      y(:,83) = (/ -17.597d0, -17.591d0, -17.589d0, -17.590d0, &
                   -17.594d0, -17.600d0, -17.608d0, -17.617d0, &
                   -17.626d0, -17.635d0, -17.645d0, -17.654d0, &
                   -17.662d0, -17.671d0, -17.679d0 /)
      y(:,84) = (/ -18.166d0, -18.134d0, -18.107d0, -18.085d0, &
                   -18.068d0, -18.056d0, -18.047d0, -18.041d0, &
                   -18.038d0, -18.036d0, -18.035d0, -18.035d0, &
                   -18.036d0, -18.038d0, -18.039d0 /)
      y(:,85) = (/ -19.000d0, -18.917d0, -18.838d0, -18.770d0, &
                   -18.714d0, -18.669d0, -18.632d0, -18.603d0, &
                   -18.579d0, -18.560d0, -18.545d0, -18.532d0, &
                   -18.522d0, -18.514d0, -18.507d0 /)
      y(:,86) = (/ -20.313d0, -19.982d0, -19.754d0, -19.592d0, &
                   -19.472d0, -19.380d0, -19.309d0, -19.253d0, &
                   -19.208d0, -19.172d0, -19.143d0, -19.119d0, &
                   -19.099d0, -19.083d0, -19.069d0 /)
      y(:,87) = (/ -19.751d0, -19.611d0, -19.520d0, -19.461d0, &
                   -19.423d0, -19.398d0, -19.382d0, -19.372d0, &
                   -19.366d0, -19.364d0, -19.363d0, -19.364d0, &
                   -19.366d0, -19.368d0, -19.371d0 /)
      y(:,88) = (/ -19.581d0, -19.431d0, -19.337d0, -19.277d0, &
                   -19.240d0, -19.218d0, -19.207d0, -19.202d0, &
                   -19.203d0, -19.207d0, -19.212d0, -19.220d0, &
                   -19.228d0, -19.236d0, -19.245d0 /)
      y(:,89) = (/ -19.685d0, -19.506d0, -19.389d0, -19.311d0, &
                   -19.258d0, -19.222d0, -19.199d0, -19.184d0, &
                   -19.175d0, -19.170d0, -19.168d0, -19.169d0, &
                   -19.171d0, -19.174d0, -19.177d0 /)
      y(:,90) = (/ -19.977d0, -19.756d0, -19.606d0, -19.501d0, &
                   -19.425d0, -19.370d0, -19.330d0, -19.300d0, &
                   -19.278d0, -19.262d0, -19.250d0, -19.241d0, &
                   -19.235d0, -19.230d0, -19.227d0 /)
      y(:,91) = (/ -20.445d0, -20.158d0, -19.958d0, -19.815d0, &
                   -19.711d0, -19.633d0, -19.574d0, -19.528d0, &
                   -19.493d0, -19.465d0, -19.442d0, -19.425d0, &
                   -19.410d0, -19.398d0, -19.389d0 /)
      y(:,92) = (/ -20.980d0, -20.625d0, -20.391d0, -20.229d0, &
                   -20.110d0, -20.020d0, -19.949d0, -19.892d0, &
                   -19.846d0, -19.807d0, -19.775d0, -19.748d0, &
                   -19.724d0, -19.704d0, -19.687d0 /)
      y(:,93) = (/ -21.404d0, -21.023d0, -20.771d0, -20.594d0, &
                   -20.461d0, -20.358d0, -20.274d0, -20.205d0, &
                   -20.148d0, -20.099d0, -20.058d0, -20.022d0, &
                   -19.991d0, -19.965d0, -19.942d0 /)
      y(:,94) = (/ -21.309d0, -20.970d0, -20.753d0, -20.603d0, &
                   -20.495d0, -20.412d0, -20.348d0, -20.295d0, &
                   -20.252d0, -20.215d0, -20.185d0, -20.158d0, &
                   -20.135d0, -20.115d0, -20.098d0 /)
      y(:,95) = (/ -21.221d0, -20.906d0, -20.707d0, -20.574d0, &
                   -20.480d0, -20.412d0, -20.361d0, -20.322d0, &
                   -20.292d0, -20.268d0, -20.249d0, -20.233d0, &
                   -20.221d0, -20.210d0, -20.201d0 /)
      y(:,96) = (/ -21.441d0, -21.097d0, -20.878d0, -20.728d0, &
                   -20.623d0, -20.546d0, -20.489d0, -20.446d0, &
                   -20.413d0, -20.387d0, -20.368d0, -20.352d0, &
                   -20.340d0, -20.330d0, -20.322d0 /)
      y(:,97) = (/ -21.668d0, -21.305d0, -21.071d0, -20.911d0, &
                   -20.797d0, -20.713d0, -20.650d0, -20.602d0, &
                   -20.565d0, -20.536d0, -20.514d0, -20.496d0, &
                   -20.481d0, -20.470d0, -20.460d0 /)
      y(:,98) = (/ -21.926d0, -21.556d0, -21.316d0, -21.150d0, &
                   -21.031d0, -20.942d0, -20.874d0, -20.822d0, &
                   -20.782d0, -20.750d0, -20.724d0, -20.704d0, &
                   -20.687d0, -20.674d0, -20.663d0 /)
      y(:,99) = (/ -22.319d0, -21.937d0, -21.686d0, -21.510d0, &
                   -21.380d0, -21.282d0, -21.206d0, -21.147d0, &
                   -21.099d0, -21.061d0, -21.031d0, -21.006d0, &
                   -20.985d0, -20.968d0, -20.954d0 /)
      y(:,100) = (/ -22.969d0, -22.561d0, -22.288d0, -22.092d0, &
                    -21.945d0, -21.832d0, -21.743d0, -21.672d0, &
                    -21.616d0, -21.570d0, -21.533d0, -21.503d0, &
                    -21.477d0, -21.457d0, -21.439d0 /)
      y(:,101) = (/ -24.001d0, -23.527d0, -23.199d0, -22.957d0, &
                    -22.772d0, -22.629d0, -22.516d0, -22.427d0, &
                    -22.355d0, -22.297d0, -22.250d0, -22.212d0, &
                    -22.180d0, -22.153d0, -22.131d0 /)
      y(:,102) = (/ -24.233d0, -23.774d0, -23.477d0, -23.273d0, &
                    -23.128d0, -23.022d0, -22.943d0, -22.883d0, &
                    -22.837d0, -22.802d0, -22.774d0, -22.752d0, &
                    -22.735d0, -22.721d0, -22.710d0 /)
      y(:,103) = (/ -24.550d0, -23.913d0, -23.521d0, -23.266d0, &
                    -23.094d0, -22.976d0, -22.893d0, -22.836d0, &
                    -22.796d0, -22.768d0, -22.750d0, -22.737d0, &
                    -22.730d0, -22.726d0, -22.725d0 /)
      y(:,104) = (/ -24.301d0, -23.665d0, -23.274d0, -23.019d0, &
                    -22.848d0, -22.730d0, -22.648d0, -22.591d0, &
                    -22.552d0, -22.525d0, -22.507d0, -22.495d0, &
                    -22.489d0, -22.485d0, -22.485d0 /)
      y(:,105) = (/ -24.519d0, -23.883d0, -23.491d0, -23.237d0, &
                    -23.065d0, -22.948d0, -22.866d0, -22.809d0, &
                    -22.770d0, -22.743d0, -22.724d0, -22.713d0, &
                    -22.706d0, -22.703d0, -22.702d0 /)

      ! Frequency quantities
      pEeV = convF*freq
      pE = pEeV*2d21*c*freq*freq
      pEeV = pEeV*ergtoev
      arg = c2*freq*1d4

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! If out of frequency limits, no contribution
      if (pEeV.lt.x1(1).or.pEeV.gt.x1(n1)) then

        ! Return with zeros
        return

      ! Within limits in axis 1
      else

        ! For every height
        do iz=iz0,iz1

          ! If within temperature limits
          if (T(iz).ge.x2(1).and.T(iz).le.x2(n2)) then

            ! Bilinear interpolation of cross-section
            call bilinear(x2,x1,y,T(iz),pEeV,cc)

            ! It was logarithmic
            cc = 1d1**cc

            ! Compute exponential
            exu = arg/T(iz)
            exu = diexp(exu)

            ! Compute absorptivity and emissivity
            eta(iz) = nCH(iz)*(1d0 - exu)*cc
            eps(iz) = nCH(iz)*pE*exu*cc

          end if ! Within temperature limits

        end do ! Heights

      end if ! Within limits in axis 1

      end subroutine CH_bf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes H2- free-free absorptivity\n
      !!    freq(double): Frequency\n
      !!  nH2(double(:)): H2 density\n
      !!    T(double(:)): Temperature\n
      !!   ne(double(:)): Electron density\n
      !!    iz0(integer): First height index to consider\n
      !!    iz1(integer): Last height index to consider\n
      !!  eta(double(:)): Absorptivity
      subroutine H2m_ff(freq,nH2,T,ne,iz0,iz1,eta)

      ! I/O

      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: nH2, T, ne
      double precision, dimension(iz0:iz1), intent(out):: eta

      ! Local

      integer, parameter:: n2=8, n1=19
      integer:: iz

      double precision:: lambda, pres, sig, tev
      double precision, dimension(n1):: x1
      double precision, dimension(n2):: x2
      double precision, dimension(n2,n1):: y


      ! Wavelength [nm]
      x1 = (/ 0.0d0, 350.5d0, 414.2d0, 506.3d0, 569.6d0, 650.9d0, &
            759.4d0, 911.3d0, 1139.1d0, 1518.8d0, 1822.6d0, &
            2278.3d0, 3037.7d0, 3645.2d0, 4556.5d0, 6075.3d0, &
            9113.0d0, 11391.3d0, 15188.3d0 /)

      ! Tempeature [eV]
      x2 = (/ 0.5d0, 0.8d0, 1.0d0, 1.2d0, 1.6d0, 2.0d0, &
              2.8d0, 3.6d0 /)

      ! Cross-section
      y(:,1)  = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
      y(:,2)  = (/ 4.17d-2, 6.10d-2, 7.34d-2, 8.59d-2, 1.11d-1, &
                   1.37d-1, 1.87d-1, 2.40d-1 /)
      y(:,3)  = (/ 5.84d-2, 8.43d-2, 1.01d-1, 1.17d-1, 1.49d-1, &
                   1.82d-1, 2.49d-1, 3.16d-1 /)
      y(:,4)  = (/ 8.70d-2, 1.24d-1, 1.46d-1, 1.67d-1, 2.10d-1, &
                   2.53d-1, 3.39d-1, 4.27d-1 /)
      y(:,5)  = (/ 1.10d-1, 1.54d-1, 1.80d-1, 2.06d-1, 2.55d-1, &
                   3.05d-1, 4.06d-1, 5.07d-1 /)
      y(:,6)  = (/ 1.43d-1, 1.98d-1, 2.30d-1, 2.59d-1, 3.17d-1, &
                   3.75d-1, 4.92d-1, 6.09d-1 /)
      y(:,7)  = (/ 1.92d-1, 2.64d-1, 3.03d-1, 3.39d-1, 4.08d-1, &
                   4.76d-1, 6.13d-1, 7.51d-1 /)
      y(:,8)  = (/ 2.73d-1, 3.71d-1, 4.22d-1, 4.67d-1, 5.52d-1, &
                   6.33d-1, 7.97d-1, 9.63d-1 /)
      y(:,9)  = (/ 4.20d-1, 5.64d-1, 6.35d-1, 6.97d-1, 8.06d-1, &
                   9.09d-1, 1.11d0, 1.32d0 /)
      y(:,10) = (/ 7.36d-1, 9.75d-1, 1.09d0, 1.18d0, 1.34d0,  &
                   1.48d0, 1.74d0, 2.01d0 /)
      y(:,11) = (/ 1.05d0, 1.39d0, 1.54d0, 1.66d0, 1.87d0, &
                   2.04d0, 2.36d0, 2.68d0 /)
      y(:,12) = (/ 1.63d0, 2.14d0, 2.36d0, 2.55d0, 2.84d0, &
                   3.07d0, 3.49d0, 3.90d0 /)
      y(:,13) = (/ 2.89d0, 3.76d0, 4.14d0, 4.44d0, 4.91d0, &
                   5.28d0, 5.90d0, 6.44d0 /)
      y(:,14) = (/ 4.15d0, 5.38d0, 5.92d0, 6.35d0, 6.99d0, &
                   7.50d0, 8.32d0, 9.02d0 /)
      y(:,15) = (/ 6.47d0, 8.37d0, 9.20d0, 9.84d0, 1.08d1, &
                   1.16d1, 1.28d1, 1.38d1 /)
      y(:,16) = (/ 1.15d1,1.48d1, 1.63d1, 1.74d1, 1.91d1,  &
                   2.04d1, 2.24d1, 2.40d1 /)
      y(:,17) = (/ 2.58d1, 3.33d1, 3.65d1, 3.90d1, 4.27d1,  &
                   4.54d1, 4.98d1, 5.33d1 /)
      y(:,18) = (/ 4.03d1, 5.20d1, 5.70d1, 6.08d1, 6.65d1,  &
                   7.08d1, 7.76d1, 8.30d1 /)
      y(:,19) = (/ 7.16d1, 9.23d1, 1.01d2, 1.08d2, 1.18d2,  &
                   1.26d2, 1.38d2, 1.47d2 /)

      ! Get wavelength in nm
      lambda = 1d2/freq

      ! Initialize RT coefficient
      eta = 0d0

      ! If above limits, no contribution
      if (lambda.gt.x1(n1)) return

      ! For each height
      do iz=iz0,iz1

        ! Check that there is some H2
        if (nH2(iz).lt.1d-30) cycle

        ! Transform T to eV
        tev = ktoev/T(iz)

        ! Bilinear interpolation of cross section
        call bilinear(x2,x1,y,tev,lambda,sig)

        ! Electron pressure
        pres = ne(iz)*kb*T(iz)

        ! Compute absorptivity
        ! 1d-19 = 1d-23 (1d6 * 1d-29) * 1d6 (cm^-3 -> m^-3)
        !         1d-2 (m^-1 -> cm^-1)
        eta(iz) = nH2(iz)*1d-19*pres*sig

      end do ! Heights

      return

      end subroutine H2m_ff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes H + p+ free-free absorptivity\n
      !!      freq(double): Frequency\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!      T(double(:)): Temperature\n
      !!      iz0(integer): First height index to consider\n
      !!      iz1(integer): Last height index to consider\n
      !!    eta(double(:)): Absorptivity
      subroutine HHp_ff(freq,Atom,T,iz0,iz1,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: T
      double precision, dimension(iz0:iz1), intent(out):: eta

      ! Local

      integer, parameter:: n2=10, n1=15
      integer:: iz, ij

      double precision:: lambda, sig
      double precision, dimension(n1):: x1
      double precision, dimension(n2):: x2
      double precision, dimension(n2,n1):: y
      double precision, dimension(iz0:iz1):: np, n0


      ! Wavelength [nm]
      x1 = (/ 0d0, 384.6d0, 555.6d0, 833.3d0, 1111.1d0, &
              1428.6d0, 1666.7d0, 2d3, 2.5d3, 2857.1d0, &
              3333.3d0, 4d3, 5d3, 6666.7d0, 1d4 /)

      ! Temperature [K]
      x2 = (/ 2.5d3, 3d3, 3.5d3, 4d3, 5d3, 6d3, 7d3, 8d3, &
              1d4, 1.2d4 /)

      ! Cross-section
      y(:,1)  = (/ 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, &
                   0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0 /)
      y(:,2)  = (/ 0.46d0, 0.46d0, 0.42d0, 0.39d0, 0.36d0, &
                   0.33d0, 0.32d0, 0.30d0, 0.27d0, 0.25d0 /)
      y(:,3)  = (/ 0.70d0, 0.62d0, 0.59d0, 0.56d0, 0.51d0, &
                   0.43d0, 0.41d0, 0.39d0, 0.35d0, 0.34d0 /)
      y(:,4)  = (/ 0.92d0, 0.86d0, 0.80d0, 0.76d0, 0.70d0, &
                   0.64d0, 0.59d0, 0.55d0, 0.48d0, 0.43d0 /)
      y(:,5)  = (/ 1.11d0, 1.04d0, 0.96d0, 0.91d0, 0.82d0, &
                   0.74d0, 0.68d0, 0.62d0, 0.53d0, 0.46d0 /)
      y(:,6)  = (/ 1.26d0, 1.19d0, 1.09d0, 1.02d0, 0.90d0, &
                   0.80d0, 0.72d0, 0.66d0, 0.55d0, 0.48d0 /)
      y(:,7)  = (/ 1.37d0, 1.25d0, 1.15d0, 1.07d0, 0.93d0, &
                   0.83d0, 0.74d0, 0.67d0, 0.56d0, 0.49d0 /)
      y(:,8)  = (/ 1.44d0, 1.32d0, 1.21d0, 1.12d0, 0.97d0, &
                   0.84d0, 0.75d0, 0.67d0, 0.56d0, 0.48d0 /)
      y(:,9)  = (/ 1.54d0, 1.39d0, 1.26d0, 1.15d0, 0.98d0, &
                   0.85d0, 0.75d0, 0.67d0, 0.55d0, 0.46d0 /)
      y(:,10) = (/ 1.58d0, 1.42d0, 1.27d0, 1.16d0, 0.98d0, &
                   0.84d0, 0.74d0, 0.66d0, 0.54d0, 0.45d0 /)
      y(:,11) = (/ 1.62d0, 1.43d0, 1.28d0, 1.15d0, 0.97d0, &
                   0.83d0, 0.72d0, 0.64d0, 0.52d0, 0.44d0 /)
      y(:,12) = (/ 1.63d0, 1.43d0, 1.27d0, 1.14d0, 0.95d0, &
                   0.80d0, 0.70d0, 0.62d0, 0.50d0, 0.42d0 /)
      y(:,13) = (/ 1.62d0, 1.40d0, 1.23d0, 1.10d0, 0.90d0, &
                   0.77d0, 0.66d0, 0.59d0, 0.48d0, 0.39d0 /)
      y(:,14) = (/ 1.55d0, 1.33d0, 1.16d0, 1.03d0, 0.84d0, &
                   0.71d0, 0.60d0, 0.53d0, 0.43d0, 0.36d0 /)
      y(:,15) = (/ 1.39d0, 1.18d0, 1.02d0, 0.90d0, 0.73d0, &
                   0.60d0, 0.52d0, 0.46d0, 0.37d0, 0.31d0 /)

      ! Get wavelength in nm
      lambda = 1d2/freq

      ! Initialize RT coefficient
      eta = 0d0

      ! If above limits, no contribution
      if (lambda.gt.x1(n1)) return

      ! Get population of HII (must be last level)
      np = Atom%popu(Atom%nlevel,iz0:iz1)

      ! Compute population of ground term
      n0 = 0d0
      do ij=1,Atom%nJ(1)
        n0 = n0 + Atom%popu(ij,iz0:iz1)
      end do

      ! For each height
      do iz=iz0,iz1

        ! Bilinear interpolation of cross section
        call bilinear(x2,x1,y,T(iz),lambda,sig)

        ! Compute absorptivity
        ! 1d-39 = 1d-23 (1d6 * 1d-29) * 1d-14 (1d6 * 1d-20)
        !      1d-2 (m^-1 -> cm^-1)
        eta(iz) = n0(iz)*np(iz)*sig*1d-39

      end do ! Heights

      return

      end subroutine HHp_ff

!#####################################################################
!#####################################################################
!#####################################################################

      end module backgroundaux_mod
