      !> Radiation transfer coefficients
      module rtcoeffaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Contributors:
!     John Dennis (NCAR)
!  Start:
!     27/04/2017
!  Last version:
!     03/10/2025 V4.0.6
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/10/2025:    V4.0.6 - Bugfix: Wrong index in the integral
!                             of JKQC with ad-hoc JKQ (TdPA)
!                           - Bugfix: The Q index started from 0 when
!                             applying the JKQ symmetries (TdPA)
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
!  Info:
!
!   - Units for absorption are cm^-1 (true absorption coefficient).
!     Needs to be multiplied by the actual atomic density.
!
!   - Units for emissivity are given in number of photons per unit
!     interval of time (s) and normalized frequency (in units of
!     Doppler width), emitted by a unit volume of gas (cm^-3) of unit
!     atomic density, within one steradian. In order to compute
!     photometric values of the intensity, the output needs be
!     multiplied by the actual atomic density.
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  rt1ordNB
!    Calculate the absorptivity and emissivity of a given atomic
!  line in the absence of magnetic field
!
!  absorbNB
!    Calculate the absorptivity of a given atomic line in the absence
!  of magnetic field
!
!  emissNB
!    Calculate the emissivity of a given atomic line in the absence
!  of magnetic field. This one is used to get the emissivity of a PRD
!  line and thus calculates the emissivity in the comoving frame for
!  all directions
!
!  emiss2ordNB
!    Calculate the second order emissivity of a given atomic line in
!  the absence of magnetic field. This subroutine only computes the
!  positive contribution of the coherent scattering and the negative
!  contribution for the flat-spectrum, i.e., the result needs to be
!  added to the product of emissNB to get the actual emissivity
!
!  rt1ord
!    Calculate the absorptivity and emissivity of a given atomic
!  line in the presence of a magnetic field
!
!  absorb
!    Calculate the absorptivity of a given atomic line in the presence
!  of a magnetic field
!
!  emiss
!    Calculate the emissivity of a given atomic line in the presence
!  of a magnetic field. This one is used to get the emissivity of a
!  PRD line and thus calculates the emissivity in the comoving frame
!  for all directions
!
!  emiss2ord
!    Calculate the second order emissivity of a given atomic line in
!  the presence of a magnetic field. This subroutine only computes the
!  positive contribution of the coherent scattering and the negative
!  contribution for the flat-spectrum, i.e., the result needs to be
!  added to the product of emiss to get the actual emissivity
!
!  absorb1
!    Calculate the two most internal loops in the calculation of the
!  absorption coefficients in the presence of a magnetic field
!
!  emiss1
!    Calculate the two most internal loops in the calculation of the
!  emission coefficients in the presence of a magnetic field
!
!  get_Warr
!    Calculate the redistribution function
!
!  abs_inner_atomic_loop
!    Calculate the inner atomic loops in the calculation of the
!  absorptivity in the presence of a magnetic field
!
!  emi_inner_atomic_loop
!    Calculate the inner atomic loops in the calculation of the
!  emissivity in the presence of a magnetic field
!
!  emi2_inner_atomic_loop
!    Calculate the inner atomic loop in the calculation of the
!  second order emissivity in the presence of a magnetic field
!
!  rt1ordLTE
!    Calculate the absorptivity and emissivity of a given LTE line in
!  the presence of a magnetic field
!
!  absorbLTE
!    Calculate the absorptivity of a given LTE line in the presence of
!  a magnetic field
!
!  photoeps
!    Calculate the emissivity of a given recombination transition
!
!  photoepsS
!    Calculate the emissivity of a given recombination transition with
!  frequency quantities stored in RAM
!
!  getJKQstar
!    Calculate the frequency dependent JKQ for the angle-averaged
!  second order emissivity in the presence of velocities in the
!  comoving frame
!
!  getJKQADasym
!    Calculate the frequency dependent JKQ for the angle-dependent
!  second order emissivity
!
!  getStkinnu
!    Interpolate the Stokes parameters into the requested frequency
!
!  getStkin
!    Interpolate the Stokes parameters into the input frequency axis
!
!  getJKQinnu
!    Interpolate the frequency dependent JKQ into the requested
!  frequency
!
!  getJKQin
!    Interpolate the frequency dependent JKQ into the input frequency
!  axis
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use funnj_mod
      use math_mod
      use profile_mod
      use parameters_mod , only : cZero, TINYJS , TINYEV , TINYCO , &
                                  TINYER , TINYB , TINYO , TINYWAR, &
                                  sqrt3, IPI, c , c2 , convF , &
                                  cSaha , fktoJ , kb, IPI41, IPI42, &
                                  vrfrac , sqrt5 , B2LK
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity and emissivity of a given atomic
      !! line in the absence of magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  TS(dcomplx(:,:,:)): Geometrical tensors in the vertical
      !!                      reference frame\n
      !!    omega(double(:)): Frequency array\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!      itran(integer): Index of transition to compute\n
      !!     itermu(integer): Upper term of the transition\n
      !!     iterml(integer): Lower term of the transition\n
      !!         iz(integer): Height index\n
      !!        if0(integer): First frequency index for this
      !!                      transition\n
      !!        if1(integer): Last frequency index for this
      !!                      transition\n
      !!   Norma(Prof_class): Normalization factors for Voigt
      !!                      profiles or Voigt profiles\n
      !!          Dw(double): Doppler width of the transition\n
      !!        vfac(double): Doppler shift factor\n
      !!        absK(double): Unit transformation factor\n
      !!     eta0(double(:)): Intensity absorptivity\n
      !!     eta1(double(:)): Q absorptivity\n
      !!     eta2(double(:)): U absorptivity\n
      !!     eta3(double(:)): V absorptivity\n
      !!     rha1(double(:)): Q dichroic absorptivity\n
      !!     rha2(double(:)): U dichroic absorptivity\n
      !!     rha3(double(:)): V dichroic absorptivity\n
      !!     eps0(double(:)): Intensity emissivity\n
      !!     eps1(double(:)): Q emissivity\n
      !!     eps2(double(:)): U emissivity\n
      !!     eps3(double(:)): V emissivity\n
      !!     rhs1(double(:)): Q 'dichroic' emissivity\n
      !!     rhs2(double(:)): U 'dichroic' emissivity\n
      !!     rhs3(double(:)): V 'dichroic' emissivity
      subroutine rt1ordNB(Atom,TS,omega,Flgsg,itran,itermu,iterml, &
                          iz,if0,if1,Norma,Dw,vfac,absK, &
                          eta0,eta1,eta2,eta3,rha1,rha2,rha3, &
                          eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iz,if0,if1
      double precision, intent(in):: Dw,absK,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TS

      ! Local

      integer:: ifreq,K,iQ,iU,iU1,iL,iL1,iR,indU,indL,indK,Kmin,Kmax

      double precision:: rLu,rLl,S,rJu,rJu1,rJl,rJl1
      double precision:: iDw,at,Dfreq,vfacw,tempRe,tempRa
      double precision:: f61,f62a,f62e,f63,f64,eu,el,rK,au,al,aul

      complex(kind=8), dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

      eta0 = 0d0
      eta1 = 0d0
      eta2 = 0d0
      eta3 = 0d0
      rha1 = 0d0
      rha2 = 0d0
      rha3 = 0d0
      eps0 = 0d0
      eps1 = 0d0
      eps2 = 0d0
      eps3 = 0d0
      rhs1 = 0d0
      rhs2 = 0d0
      rhs3 = 0d0

      !
      ! Get terms and transition quantities
      !

      ! Inverse Doppler width
      iDw = 1d0/Dw

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)*iDw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Doppler shift in doppler units
      vfacw = vfac*iDw

      !
      ! Common part
      !

      ! For each Ju
      do iU=1,Atom%nJ(itermu)

        ! Get Ju
        rJu = Atom%rJval(iU,itermu)

        ! Get eigenvalue upper level
        eu = Atom%FSfreq(iU,itermu)

        ! Level index
        indU = Atom%irho(itermu)%irho_ij(iU)

        ! Factor Ju
        f61 = 2d0*rJu + 1d0

        ! For each Jl
        do iL=1,Atom%nJ(iterml)

          ! Get Jl
          rJl = Atom%rJval(iL,iterml)

          ! 6-J
          f62e = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

          ! Small J-symbol
          if (abs(f62e).lt.TINYJS) cycle

          ! Absorption factor
          f62a = f62e*f61*Flgsg%sg(nint(1d0+rJu+rJl))* &
                 sqrt(2d0*rJl+1d0)

          ! Emission factor
          f62e = f62e*sqrt(f61)*(2d0*rJl+1d0)

          ! Get eigenvalue lower level
          el = Atom%FSfreq(iL,iterml)

          ! Level and transition indexes
          indL = Atom%irho(iterml)%irho_ij(iL)
          indK = Atom%trano(itran)%indNB(indL,indU)

          !
          ! Compute profile
          !

          ! If stored
          if (Norma%VRAM) then

            ! Copy stored profile
            prof = Norma%cp(:,indK)

          ! Not stored
          else

            ! Shift term
            Dfreq = (eu - el)*iDw

            ! For each frequency
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreq - omega(ifreq)*vfacw,at,prof(ifreq))

            end do ! frequencies

            ! Normalize profile
            prof = dcmplx(dble(prof)*Norma%Norm(indK), dimag(prof))

          end if ! Storing

          !
          ! Absorption
          !

          ! For each Jl'
          do iL1=1,Atom%nJ(iterml)

            ! Get Jl'
            rJl1 = Atom%rJval(iL1,iterml)

            ! 6-J
            f63 = fun6j(rLu,rLl,1d0,rJl1,rJu,S,Flgsg)

            ! Small J-symbol
            if (abs(f63).lt.TINYJS) cycle

            ! Combine factors
            f63 = f62a*f63*sqrt(2d0*rJl1+1d0)

            ! Determine the limits in K
            Kmin = nint(abs(rJl-rJl1))
            Kmax = min(nint(rJl+rJl1),Atom%Kcut(iterml),2)

            ! For each K
            do K=Kmin,Kmax

              ! Get the real number
              rK = dble(K)

              ! For each Q
              do iQ=-K,K

                ! Get the SEE index
                iR = Atom%irho(iterml)%Jrho(iL1,iL)%kq(iQ,K)

                ! If flagged as small, skip
                if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

                ! Racah algebra
                f64 = fun6j(1d0,1d0,rK,rJl1,rJl,rJu,Flgsg)

                ! Combine everything but geometrical tensor
                profk = f64*f63*Flgsg%sg(K)*prof*Atom%crho(iR,iz)

                ! Absorptivity
                eta0 = eta0 + dble(TS(0,iQ,K)*profk)
                eta1 = eta1 + dble(TS(1,iQ,K)*profk)
                eta2 = eta2 + dble(TS(2,iQ,K)*profk)
                eta3 = eta3 + dble(TS(3,iQ,K)*profk)

                ! Dispersion
                rha1 = rha1 + dimag(TS(1,iQ,K)*profk)
                rha2 = rha2 + dimag(TS(2,iQ,K)*profk)
                rha3 = rha3 + dimag(TS(3,iQ,K)*profk)

              end do ! Q
            end do ! K
          end do ! iL1

          !
          ! Emission
          !

          ! For each Ju'
          do iU1=1,Atom%nJ(itermu)

            ! Ju'
            rJu1 = Atom%rJval(iU1,itermu)

            ! 6-J
            f63 = fun6j(rLu,rLl,1d0,rJl,rJu1,S,Flgsg)

            ! Small J-symbol
            if (abs(f63).lt.TINYJS) cycle

            ! Combine factors
            f63 = f63*f62e*sqrt(2d0*rJu1+1d0)* &
                  Flgsg%sg(nint(1d0+rJl+rJu1))

            ! Determine the limits in K
            Kmin = nint(abs(rJu-rJu1))
            Kmax = min(nint(rJu+rJu1),Atom%Kcut(itermu),2)

            ! For each K
            do K=Kmin,Kmax

              ! Get the real number
              rK = dble(K)

              ! For each Q
              do iQ=-K,K

                ! Get the SEE index
                iR = Atom%irho(itermu)%Jrho(iU,iU1)%kq(iQ,K)

                ! If flagged as small, skip
                if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

                ! 6-J
                f64 = fun6j(1d0,1d0,rK,rJu,rJu1,rJl,Flgsg)

                ! Combine everything but geometrical tensor
                profk = f64*f63*prof*Atom%crho(iR,iz)

                ! Emissivity
                eps0 = eps0 + dble(TS(0,iQ,K)*profK)
                eps1 = eps1 + dble(TS(1,iQ,K)*profK)
                eps2 = eps2 + dble(TS(2,iQ,K)*profK)
                eps3 = eps3 + dble(TS(3,iQ,K)*profK)

                ! Dispersion
                rhs1 = rhs1 + dimag(TS(1,iQ,K)*profK)
                rhs2 = rhs2 + dimag(TS(2,iQ,K)*profK)
                rhs3 = rhs3 + dimag(TS(3,iQ,K)*profK)

              end do ! Q
            end do ! K
          end do ! iU1

          ! Return to common loop

        end do ! iU
      end do ! iL

      ! Common parts for coefficients
      tempRe = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
               Atom%Ecoeff(itermu,iterml)*iDw

      ! Coefficient for absorption
      tempRa = tempRe/absK

      ! Final values
      eta0 = tempRa*eta0
      eta1 = tempRa*eta1
      eta2 = tempRa*eta2
      eta3 = tempRa*eta3
      rha1 = tempRa*rha1
      rha2 = tempRa*rha2
      rha3 = tempRa*rha3
      eps0 = tempRe*eps0
      eps1 = tempRe*eps1
      eps2 = tempRe*eps2
      eps3 = tempRe*eps3
      rhs1 = tempRe*rhs1
      rhs2 = tempRe*rhs2
      rhs3 = tempRe*rhs3

      end subroutine rt1ordNB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity of a given atomic line in the
      !! absence of magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  TS(dcomplx(:,:,:)): Geometrical tensors in the vertical
      !!                      reference frame\n
      !!    omega(double(:)): Frequency array\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!      itran(integer): Index of transition to compute\n
      !!     itermu(integer): Upper term of the transition\n
      !!     iterml(integer): Lower term of the transition\n
      !!         iz(integer): Height index\n
      !!        if0(integer): First frequency index for this
      !!                      transition\n
      !!        if1(integer): Last frequency index for this
      !!                      transition\n
      !!   Norma(Prof_class): Normalization factors for Voigt
      !!                      profiles or Voigt profiles\n
      !!          Dw(double): Doppler width of the transition\n
      !!        vfac(double): Doppler shift factor\n
      !!        absK(double): Unit transformation factor\n
      !!     eta0(double(:)): Intensity absorptivity\n
      !!     eta1(double(:)): Q absorptivity\n
      !!     eta2(double(:)): U absorptivity\n
      !!     eta3(double(:)): V absorptivity\n
      !!     rha1(double(:)): Q dichroic absorptivity\n
      !!     rha2(double(:)): U dichroic absorptivity\n
      !!     rha3(double(:)): V dichroic absorptivity
      subroutine absorbNB(Atom,TS,omega,Flgsg,itran,itermu,iterml, &
                          iz,if0,if1,Norma,Dw,vfac,absK, &
                          eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iz,if0,if1
      double precision, intent(in):: Dw,absK,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TS

      ! Local

      integer:: ifreq,K,iQ,iU,iL,iL1,iR,Kmin,Kmax,indU,indL,indK

      double precision:: rLu,rLl,S,rJu,rJl,rJl1,f61,f62,f63,f64
      double precision:: eu,el,rK,au,al,aul,iDw,at,Dfreq,vfacw,tempR

      complex(kind=8), dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

      eta0 = 0d0
      eta1 = 0d0
      eta2 = 0d0
      eta3 = 0d0
      rha1 = 0d0
      rha2 = 0d0
      rha3 = 0d0

      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)*iDw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Doppler shift in doppler units
      vfacw = vfac*iDw


      !
      ! Compute absorptivity
      !

      ! For each Ju
      do iU=1,Atom%nJ(itermu)

        ! Get Ju
        rJu = Atom%rJval(iU,itermu)

        ! Get eigenvalue upper level
        eu = Atom%FSfreq(iU,itermu)

        ! Factor Ju
        f61 = 2d0*rJu + 1d0

        ! Level index
        indU = Atom%irho(itermu)%irho_ij(iU)

        ! For each Jl
        do iL=1,Atom%nJ(iterml)

          ! Get Jl
          rJl = Atom%rJval(iL,iterml)

          ! 6-J
          f62 = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

          ! Small J-symbol
          if (abs(f62).lt.TINYJS) cycle

          ! Complete factor
          f62 = f62*f61*Flgsg%sg(nint(1d0+rJu+rJl))*sqrt(2d0*rJl+1d0)

          ! Get eigenvalue lower level
          el = Atom%FSfreq(iL,iterml)

          ! Level and transition index
          indL = Atom%irho(iterml)%irho_ij(iL)
          indK = Atom%trano(itran)%indNB(indL,indU)

          !
          ! Compute profile
          !

          ! If stored
          if (Norma%VRAM) then

            ! Copy stored profile
            prof = Norma%cp(:,indK)

          ! Not stored
          else

            ! Shift term
            Dfreq = (eu - el)*iDw

            ! For each frequency
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreq - omega(ifreq)*vfacw,at,prof(ifreq))

            end do ! frequencies

            ! Normalize profile
            prof = dcmplx(dble(prof)*Norma%Norm(indK), &
                          dimag(prof))

          end if ! Storing

          ! For each Jl'
          do iL1=1,Atom%nJ(iterml)

            ! Get Jl
            rJl1 = Atom%rJval(iL1,iterml)

            ! 6-j
            f63 = fun6j(rLu,rLl,1d0,rJl1,rJu,S,Flgsg)

            ! Small J-symbol
            if (abs(f63).lt.TINYJS) cycle

            ! Complete factor
            f63 = f62*f63*sqrt(2d0*rJl1+1d0)

            ! Determine the limits in K
            Kmin = nint(abs(rJl-rJl1))
            Kmax = min(nint(rJl+rJl1),Atom%Kcut(iterml),2)

            ! For each K
            do K=Kmin,Kmax

              ! Get the real number
              rK = dble(K)

              ! For each Q
              do iQ=-K,K

                ! Get the SEE index
                iR = Atom%irho(iterml)%Jrho(iL1,iL)%kq(iQ,K)

                ! If flagged as small, skip
                if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

                ! Racah algebra
                f64 = fun6j(1d0,1d0,rK,rJl1,rJl,rJu,Flgsg)

                ! Combine everything but geometrical tensor
                profk = f64*f63*Flgsg%sg(K)*prof*Atom%crho(iR,iz)

                ! Absorptivity
                eta0 = eta0 + dble(TS(0,iQ,K)*profk)
                eta1 = eta1 + dble(TS(1,iQ,K)*profk)
                eta2 = eta2 + dble(TS(2,iQ,K)*profk)
                eta3 = eta3 + dble(TS(3,iQ,K)*profk)

                ! Dispersion
                rha1 = rha1 + dimag(TS(1,iQ,K)*profk)
                rha2 = rha2 + dimag(TS(2,iQ,K)*profk)
                rha3 = rha3 + dimag(TS(3,iQ,K)*profk)

              end do ! Q
            end do ! K
          end do ! iL1
        end do ! iL
      end do ! iU

      ! Common part for the two coefficients
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)*iDw/absK

      ! Final values
      eta0 = tempR*eta0
      eta1 = tempR*eta1
      eta2 = tempR*eta2
      eta3 = tempR*eta3
      rha1 = tempR*rha1
      rha2 = tempR*rha2
      rha3 = tempR*rha3

      return

      end subroutine absorbNB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity of a given atomic line in the
      !! absence of magnetic field. This one is used to get the
      !! emissivity of a PRD line and thus calculates the emissivity
      !! in the comoving frame for all directions\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!  TS(dcomplx(:,:,:,:)): Geometrical tensors in the vertical
      !!                        reference frame\n
      !!      omega(double(:)): Frequency array\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!        itran(integer): Index of transition to compute\n
      !!       itermu(integer): Upper term of the transition\n
      !!       iterml(integer): Lower term of the transition\n
      !!           iz(integer): Height index\n
      !!          if0(integer): First frequency index for this
      !!                        transition\n
      !!          if1(integer): Last frequency index for this
      !!                        transition\n
      !!        njdir(integer): Number of directions\n
      !!     Norma(Prof_class): Normalization factors for Voigt
      !!                        profiles or Voigt profiles\n
      !!            Dw(double): Doppler width of the transition\n
      !!       eps0(double(:)): Intensity emissivity\n
      !!       eps1(double(:)): Q emissivity\n
      !!       eps2(double(:)): U emissivity\n
      !!       eps3(double(:)): V emissivity
      subroutine emissNB(Atom,TS,omega,Flgsg,itran,itermu, &
                         iterml,iz,if0,if1,njdir,Norma,Dw, &
                         eps0,eps1,eps2,eps3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iz,if0,if1,njdir
      double precision, intent(in):: Dw
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(njdir,if0:if1), intent(out):: eps0
      double precision, dimension(njdir,if0:if1), intent(out):: eps1
      double precision, dimension(njdir,if0:if1), intent(out):: eps2
      double precision, dimension(njdir,if0:if1), intent(out):: eps3
      complex(kind=8), dimension(0:3,-2:2,0:2,njdir), intent(in):: TS

      ! Local

      integer:: ifreq,K,iQ,iL,iU,iU1,iR,Kmin,Kmax,indU,indL,indK,jdir

      double precision:: rLl,rLu,S,rJl,rJu,rJu1,rK,f61,f62,f63,f64
      double precision:: el,eu,al,au,aul,tempR,iDw,at,Dfreq

      complex(kind=8), dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

      eps0 = 0d0
      eps1 = 0d0
      eps2 = 0d0
      eps3 = 0d0

      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)*iDw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)


      !
      ! Compute emissivity
      !

      ! For each Jl
      do iL=1,Atom%nJ(iterml)

        ! Get eigenvalue lower level
        el = Atom%FSfreq(iL,iterml)

        ! Jl
        rJl = Atom%rJval(iL,iterml)

        ! Jl factor
        f61 = 2d0*rJl+1d0

        ! Get index
        indL = Atom%irho(iterml)%irho_ij(iL)

        ! For each Ju
        do iU=1,Atom%nJ(itermu)

          ! Ju
          rJu = Atom%rJval(iU,itermu)

          ! 6-J
          f62 = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

          ! Small J-symbol
          if (abs(f62).lt.TINYJS) cycle

          ! Complete factor
          f62 = f62*f61*sqrt(2d0*rJu+1d0)

          ! Get eigenvalue upper level
          eu = Atom%FSfreq(iU,itermu)

          ! Get indexes
          indU = Atom%irho(itermu)%irho_ij(iU)
          indK = Atom%trano(itran)%indNB(indL,indU)

          !
          ! Compute profile
          !

          ! If stored
          if (Norma%VRAM) then

            ! Copy stored profile
            prof = Norma%cp(:,indK)

          ! Not stored
          else

            ! Shift term
            Dfreq = eu - el

            ! For each frequency
            do ifreq=if0,if1

              ! Calculate profile
              call voigt((Dfreq - omega(ifreq))*iDw,at,prof(ifreq))

            end do ! frequencies

            ! Normalize
            prof = dcmplx(dble(prof)*Norma%Norm(indK), &
                          dimag(prof))

          end if ! Storing

          ! For each Ju'
          do iU1=1,Atom%nJ(itermu)

            ! Ju'
            rJu1 = Atom%rJval(iU1,itermu)

            ! 6-J
            f63 = fun6j(rLu,rLl,1d0,rJl,rJu1,S,Flgsg)

            ! Small J-symbol
            if (abs(f63).lt.TINYJS) cycle

            ! Complete factor
            f63 = f63*f62*sqrt(2d0*rJu1+1d0)* &
                  Flgsg%sg(nint(1d0+rJl+rJu1))

            ! Determine the limits in K
            Kmin = nint(abs(rJu-rJu1))
            Kmax = min(nint(rJu+rJu1),Atom%Kcut(itermu),2)

            ! For each K
            do K=Kmin,Kmax

              ! Get the real number
              rK = dble(K)

              ! For each Q
              do iQ=-K,K

                ! Get the SEE index
                iR = Atom%irho(itermu)%Jrho(iU,iU1)%kq(iQ,K)

                ! If flagged as small, skip
                if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

                ! 6-J
                f64 = fun6j(1d0,1d0,rK,rJu,rJu1,rJl,Flgsg)

                ! Everything but geometric tensor
                profk = f64*f63*prof*Atom%crho(iR,iz)

                ! Directions
                do jdir=1,njdir

                  ! Emissivity
                  eps0(jdir,:) = eps0(jdir,:) + &
                                 dble(TS(0,iQ,K,jdir)*profK)
                  eps1(jdir,:) = eps1(jdir,:) + &
                                 dble(TS(1,iQ,K,jdir)*profK)
                  eps2(jdir,:) = eps2(jdir,:) + &
                                 dble(TS(2,iQ,K,jdir)*profK)
                  eps3(jdir,:) = eps3(jdir,:) + &
                                 dble(TS(3,iQ,K,jdir)*profK)

                end do ! Directions
              end do ! Q
            end do ! K
          end do ! iU1
        end do ! iU
      end do ! iL

      ! Factor for emissivity
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)*iDw

      ! Final values
      eps0 = tempR*eps0
      eps1 = tempR*eps1
      eps2 = tempR*eps2
      eps3 = tempR*eps3

      return

      end subroutine emissNB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the second order emissivity of a given atomic line
      !! in the absence of magnetic field. This subroutine only
      !! computes the positive contribution of the coherent scattering
      !! and the negative contribution for the flat-spectrum, i.e.,
      !! the result needs to be added to the product of emissNB to get
      !! the actual emissivity\n
      !!         Atom(Atom_class): Structure with atomic data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!               vx(double): Velocity vector along X\n
      !!               vy(double): Velocity vector along Y\n
      !!               vz(double): Velocity vector along Z\n
      !!            lvel(logical): If dynamic node\n
      !!         omega(double(:)): Frequency array\n
      !!          Fed(Reda_class): Structure with redistribution
      !!                           output frequency data\n
      !!          Red(Redb_class): Structure with redistribution input
      !!                           frequency data\n
      !!       RWarr(Redb2_class): Structure with redistribution
      !!                           function data\n
      !!        Norma(Prof_class): Normalization factors for Voigt
      !!                           profiles or Voigt profiles\n
      !!       Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                           and J-symbols\n
      !!           jtran(integer): Output transition index\n
      !!          itermu(integer): Upper term of output transition\n
      !!          itermf(integer): Lower term of output transition\n
      !!              iz(integer): Height index\n
      !!             if0(integer): First frequency index for this
      !!                           transition\n
      !!             if1(integer): Last frequency index for this
      !!                           transition\n
      !!            Mif0(integer): First frequency for this CPU\n
      !!            Mif1(integer): Last frequency for this CPU\n
      !!              DwT(double): Thermal part of Doppler width\n
      !!               Dw(double): Doppler width of the output
      !!                           transition\n
      !!              vmi(double): Microturbulent velocity\n
      !!  TSout(docmplx(:,:,:,:)): Geometrical tensor in the vertical
      !!                           reference frame\n
      !!  Stokes(double(:,:,:,:)): Stokes parameters\n
      !!    JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!     JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                           over the absorption profile\n
      !!    JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                           frequency dependence\n
      !!         eps20(double(:)): Intensity emissivity\n
      !!         eps21(double(:)): Q emissivity\n
      !!         eps22(double(:)): U emissivity\n
      !!         eps23(double(:)): V emissivity
      subroutine emiss2ordNB(Atom,Geom,vx,vy,vz,lvel,omega,Fed,Red, &
                             RWarr,Norma,Flgsg,jtran,itermu,itermf, &
                             iz,if0,if1,Mif0,Mif1,DwT,Dw,vmi,TSout, &
                             Stokes,JKQa,JKQ,JKQC, &
                             eps20,eps21,eps22,eps23)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      type(Redb2_class), intent(inout):: RWarr
      type(Prof_class), intent(in):: Norma
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: lvel
      integer, intent(in):: jtran,itermu,itermf,iz,if0,if1,Mif0,Mif1
      double precision, intent(in):: DwT,Dw,vmi,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(Geom%njdir,if0:if1), &
                        intent(out):: eps20,eps21,eps22,eps23
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), &
                       target, intent(in):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2,Geom%njdir), &
                       intent(in):: TSout

      ! Local

      logical:: PRDc,LPRAM,cohIn,conj,asym

      integer:: ith1,iph1,ish1,iterml,itran,iti,icom,jdir,idir
      integer:: mF,iU,iU1,iL,iL1,iifreq,iran,ifreq,iR,ishift
      integer:: K,iQ,K1,iQ1,iQl,Kmin,Kmax,Kl,kwfreq0,llfreq0
      integer:: jjfreq,jjfreq0,kkfreq0,kkfreq0b,nmfreq
      integer:: indF,indU,indU1,indL,indL1,indK,nfs,i0,i1

      double precision:: rLu,rLl,rLf,S,rJu,rJu1,rJl,rJl1,rJf
      double precision:: rK,Q,rK1,Q1,Ql,rKl,wlf
      double precision:: eu,eu1,el,el1,ef,au,af,al,auf,aul,at
      double precision:: f61,f62,f63,f64,f0tmp,ftmp,f1tmp,fKfj
      double precision:: Dw1,Norme0,Norme1
      double precision:: omegai,daux,sig,Dfreq,hau
      double precision:: iDw,cost,sint,cosc,sinc,vfac1
      double precision, dimension(0:3):: StokesM
      double precision, dimension(:,:), allocatable:: Stokesin

      complex(kind=8):: hanleden,prof,rhoc,intgr,Norme2,y0,PRDin
      complex(kind=8), dimension(if0:if1):: CRD0,CRD,PRDaux
      complex(kind=8), dimension(if0:if1,Geom%njdir):: tmp0
      complex(kind=8), dimension(if0:if1,Geom%njdir):: tmp1
      complex(kind=8), dimension(if0:if1,Geom%njdir):: tmp2
      complex(kind=8), dimension(if0:if1,Geom%njdir):: tmp3
      complex(kind=8), dimension(-2:2,0:2):: Jrad
      complex(kind=8), dimension(:), allocatable, target:: Warr2
      complex(kind=8), dimension(:), allocatable:: Warr2xW
      complex(kind=8), dimension(:), allocatable:: intergrin
      complex(kind=8), dimension(:,:), allocatable:: PRD
      complex(kind=8), dimension(:,:,:), allocatable, target:: JKQinMV
      complex(kind=8), dimension(:,:,:), allocatable, target:: JradC

      ! Pointers

      type(Redc_class), pointer:: p_red
      type(Redc2_class), pointer:: p_rwarr
      type(Redc2_class), target:: p_dummy
      integer, pointer:: p_mfreq
      complex(kind=8), dimension(:), pointer:: p_warr2
      complex(kind=8), dimension(:), pointer:: p_JKQ
      complex(kind=8), dimension(:), pointer:: p_JKQC
      complex(kind=8), dimension(:,:,:,:), pointer:: TKQo


      ! Routine name
      urou = 'emiss2ordNB'

      ! Initialize pointers
      nullify(p_red,p_rwarr,p_mfreq,p_warr2,p_JKQ,p_JKQC,TKQo)

      ! If dynamic and angle-average
      if (AV.and.lvel) then

        ! Point to correct tensors
        TKQo => Geom%TS

        ! Get JKQ in comoving frame
        call getJKQstar(Red,Geom,iz,Atom%ntran,Atom%tif0, &
                        Atom%tif1,DwT,vx,vy,vz,omega,Flgsg, &
                        TKQo,Stokes,JKQa,JKQinMV)

        ! Free pointer
        nullify(TKQo)

      ! If angle-dependent
      else if (.not.AV) then

        ! Check if ad-hoc asymmetries
        asym = size(JKQa).gt.10

        ! If there are ad-hoc asymmetries
        if (asym) then

          ! Point to correct tensors
          TKQo => Geom%TS

          ! Get JKQ correction for asymmetry
          call getJKQADasym(Red,Geom,iz,Atom%ntran,Atom%tif0, &
                            Atom%tif1,DwT,vx,vy,vz, &
                            .False.,0d0,0d0,omega,Flgsg, &
                            TKQo,Stokes,JKQC,JKQa,JKQinMV)

        end if ! Ad-hoc asymmetries
      end if ! Dynamic and AA

      !
      ! Allocate PRD
      !

      ! If angle-average
      if (AV) then

        ! Only one direction
        allocate(PRD(if0:if1,1))

      ! Angle-dependent
      else

        ! Allocate all output directions
        allocate(PRD(if0:if1,Geom%njdir))

      end if ! AA or AD

      !
      ! Get terms and transition quantities
      !

      ! Inverse Doppler width
      iDw = 1d0/Dw

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      at = (au + af + auf)*iDw

      ! Lifetime Hanle factor
      hau = 2d0*(au+auf)*iDw

      ! Units normalization factor for CRD profile
      Norme0 = (1d5*iDw)*.5d0*sqrt(IPI)

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLf = Atom%rLval(itermf)


      !
      ! Initialize the emission coefficient
      !
      eps20 = 0d0
      eps21 = 0d0
      eps22 = 0d0
      eps23 = 0d0

      ! If there are frequencies and angle-dependent
      if (.not.AV.and.Red%mxfreq.gt.0) then

        ! Allocate auxialiar arrays
        allocate(Warr2xW(Red%mxfreq))
        allocate(intergrin(Red%mxfreq))

      end if ! Frequencies to integrate and angle-dependent


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible low->up transitions
      do iti=1,Atom%trano(jtran)%nt

        ! Get transition index and lower term
        itran = Atom%trano(jtran)%indT(iti)
        iterml = Atom%fst(itran)%iterml

        ! Point to input transition
        p_red => Red%trani(iti)

        ! If PRAM
        if (PRAM) then

          ! Point to input transition and get if storing
          p_rwarr => RWarr%trani(iti)
          LPRAM = PRAM.and.p_rwarr%RAM

        ! If not, nothing stored
        else

          ! Not storing locally either
          p_rwarr => p_dummy
          LPRAM = .False.

        end if

        ! Predict size of interpolation block
        nmfreq = sum(p_red%mfreq)

        ! If angle-averaged
        if (AV) then

          ! If dynamic
          if (lvel) then

            ! Get input radiation field
            call getJKQin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                          JradC,JKQinMV)

          ! If static
          else

            ! Get input radiation field
            call getJKQin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                          JradC,JKQC(:,:,Red%ggf0:Red%ggf1))

          end if ! Dynamic

        ! If angle-dependent
        else

          ! If Rayleigh scattering and there is coherent
          if (jtran.eq.itran.and.Geom%V_CScatt(1).ge.1d0) then

            ! Scale dimension
            nmfreq = nmfreq*(Geom%nScatt-1)
            nfs = 1

          ! Raman scattering
          else

            ! Scale dimension
            nmfreq = nmfreq*Geom%nScatt
            nfs = 0

          end if ! Rayleigh/Raman scattering

          ! Get interpolated intensity
          call getStkin(Geom,p_red,Fed,Red,Mif0,Mif1,omega, &
                        vx,vy,vz,lvel,Stokesin,Stokes)

          ! If asymmetries, get input radiation field
          if (asym) &
            call getJKQin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                          JradC,JKQinMV)

        end if ! AA/AD

        ! Get the 'flat' JKQ for this input transition
        JRad = JKQ(:,:,itran)

        ! Doppler width for the input transition
        Dw1 = Atom%Dfreq(itran)*sqrt(DwT*DwT + vmi**2d0)

        ! Damping parameter input lower level and input transition
        al = Atom%damp(iterml,iz)
        aul = Atom%ldamp(itran,iz)

        ! Angular momentum input lower level
        rLl = Atom%rLval(iterml)

        ! Initialize kkfreq index
        kwfreq0 = 0

        ! For each Jf
        do mF=1,Atom%nJ(itermf)

          ! Get eigenvalue final lower level
          ef = Atom%FSfreq(mF,itermf) - Atom%TRfreq(itermf)

          ! Get index
          indF = Atom%irho(itermf)%irho_ij(mF)

          ! Get Jf
          rJf = Atom%rJval(mF,itermf)

          ! For each Ju
          do iU=1,Atom%nJ(itermu)

            ! Get Ju
            rJu = Atom%rJval(iU,itermu)

            ! Compute J-symbol
            f61 = fun6j(rJu,rJf,1d0,rLf,rLu,S,Flgsg)

            ! Small J symbol
            if (abs(f61).lt.TINYJS) cycle

            ! Get eigenvalue upper level
            eu = Atom%FSfreq(iU,itermu) - Atom%TRfreq(itermu)

            ! Get indexes
            indU = Atom%irho(itermu)%irho_ij(iU)
            indK = Atom%trano(jtran)%indNB(indF,indU)

            ! Complete factor
            f61 = f61*(2d0*rJu+1d0)*(2d0*rJf+1d0)


            !
            ! Flat contribution. Implicit branching
            !

            ! If stored
            if (Norma%VRAM) then

              ! Get first order profile
              CRD0 = Norme0*conjg(Norma%cp(:,indK))

            ! If not stored
            else

              ! Shift term
              Dfreq  = eu  - ef + Atom%Dfreq(jtran)

              ! Normalization
              Norme1 = Norme0*Norma%Norm(indK)

              ! For each frequency
              do ifreq=if0,if1

                ! Calculate profile u-f
                call voigt((Dfreq - omega(ifreq))*iDw,at,prof)

                ! Flat spectrum contribution, normalizing
                CRD0(ifreq) = dcmplx(Norme1*dble(prof), &
                                    -Norme0*dimag(prof))

              end do ! frequencies

            end if ! Storing Voigt

            ! For each Ju'
            do iU1=1,Atom%nJ(itermu)

              ! Get Ju'
              rJu1 = Atom%rJval(iU1,itermu)

              ! 6-J symbol
              f62 = fun6j(rJu1,rJf,1d0,rLf,rLu,S,Flgsg)

              ! Small J symbol
              if (abs(f62).lt.TINYJS) cycle

              ! Get eigenvalue upper' level
              eu1 = Atom%FSfreq(iU1,itermu) - Atom%TRfreq(itermu)

              ! Get indexes
              indU1 = Atom%irho(itermu)%irho_ij(iU1)
              indK = Atom%trano(jtran)%indNB(indF,indU1)

              ! Complete factor
              f62 = f62*(2d0*rJu1+1d0)*Flgsg%sg(nint(rJu1+rJf+1d0))

              !
              ! Flat contribution. Implicit branching
              !

              ! If stored
              if (Norma%VRAM) then

                ! Flat spectrum contribution
                CRD = CRD0 + Norme0*Norma%cp(:,indK)

              ! If not stored
              else

                ! Shift term
                Dfreq = eu1 - ef + Atom%Dfreq(jtran)

                ! Normalization factor
                Norme1 = Norme0*Norma%Norm(indK)

                ! For each frequency
                do ifreq=if0,if1

                  ! Calculate profile u'-f
                  call voigt((Dfreq - omega(ifreq))*iDw,at,prof)

                  ! Flat spectrum contribution
                  CRD(ifreq) = CRD0(ifreq) + &
                               dcmplx(Norme1*dble(prof), &
                                      Norme0*dimag(prof))

                end do ! frequencies

              end if ! Storing Voigt

              !
              ! Continue with the 2nd order emissivity
              !

              ! For each Jl
              do iL=1,Atom%nJ(iterml)

                ! Get Jl
                rJl = Atom%rJval(iL,iterml)

                ! 6-J symbol
                f63 = fun6j(rJu,rJl,1d0,rLl,rLu,S,Flgsg)

                ! Small J symbol
                if (abs(f63).lt.TINYJS) cycle

                ! Complete factor
                f63 = f63*sqrt(2d0*rJl+1d0)

                ! Get eigenvalue of input lower level
                el = Atom%FSfreq(iL,iterml) - Atom%TRfreq(iterml)

                ! Get indexes
                indL = Atom%irho(iterml)%irho_ij(iL)

                ! For each Jl'
                do iL1=1,Atom%nJ(iterml)

                  ! Get Jl1
                  rJl1 = Atom%rJval(iL1,iterml)

                  ! 6-J symbol
                  f64 = fun6j(rJu1,rJl1,1d0,rLl,rLu,S,Flgsg)

                  ! Small J-symbol
                  if (abs(f64).lt.TINYJS) cycle

                  ! Get eigenvalue of input lower' level
                  el1 = Atom%FSfreq(iL1,iterml) - Atom%TRfreq(iterml)

                  ! Get index
                  indL1 = Atom%irho(iterml)%irho_ij(iL1)

                  ! Complete factor
                  f64 = f64*sqrt(2d0*rJl1+1d0)

                  ! Hanle factor
                  ! TODO ATTENTION TO THIS
                  hanleden = dcmplx(hau,(eu-eu1)*iDw)

                  ! Initialize temporal variable
                  tmp0 = cZero
                  tmp1 = cZero
                  tmp2 = cZero
                  tmp3 = cZero

        !
        ! Reset indent
        !


        ! If storing redistribution
        if (LPRAM) then

          ! Get the component index
          icom = Atom%trano(jtran)%trani(iti)% &
                      indNB(indL1,indL,indF,indU1,indU)

          ! Get if need to compute redistribution
          PRDc = p_rwarr%iPPRD(icom)

        ! If not storing
        else

          ! Always need to calculate
          PRDc = .True.

        end if ! Storing redistribution function


        !
        ! Create array of Wfunc
        !

        ! If need to compute it and there are frequencies
        if (PRDc.and.nmfreq.gt.0) then

          ! Calculate redistribution function Warr2
          call get_Warr(Atom,Geom,Fed,p_red,p_rwarr,LPRAM,Mif0,Mif1, &
                        itran,jtran,icom,kwfreq0,nmfreq,omega, &
                        Dw,Dw1,el,el1,eu,eu1,ef, &
                        al,au,af,aul,auf,Warr2)

        end if ! Initialized

        ! Initialize factor
        f0tmp = f61*f62*f63*f64

        ! For each K
        do K=0,Krad

          ! Get real value
          rK = dble(K)

          ! Racah algebra
          ftmp = f0tmp*sqrt(2d0*rK+1d0)

          ! For each Q value
          do iQ=-K,K

            ! Get real value
            Q = dble(iQ)

      !
      ! Reset identation
      !

      !
      ! Integral over input frequencies
      !

      ! If storing Warr
      if (LPRAM) then

        ! If there are frequencies
        if (nmfreq.gt.0) then

          ! Copy from RAM
          allocate(p_warr2(nmfreq))
          p_warr2 = dcmplx(p_rwarr%Pwarr2(kwfreq0+1:kwfreq0+nmfreq))

        end if ! There are frequencies

      ! If not storing
      else

        ! Just point to the one calculated here
        if (allocated(Warr2)) p_warr2 => Warr2

      end if ! Storing Warr

      ! Initialize integral
      PRD = cZero

      ! Difference between l and f energies
      wlf = el - ef

      !
      ! Angle-average (Integral)
      !
      if (AV) then

        ! Check if coherent
        cohIn = abs(wlf).gt.0d0

        ! Check if conjugated
        conj = iQ.lt.0
        sig = Flgsg%sg(iQ)

        ! If coherent
        if (minval(p_red%mfreq).lt.1) then

          ! If dynamic
          if (lvel) then

            ! Just point
            if (conj) then
              p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(-iQ,K,:)
            else
              p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(iQ,K,:)
            end if

          ! If static
          else

            ! Just point
            if (conj) then
              p_JKQC(Red%ggf0:Red%ggf1) => &
                                        JKQC(-iQ,K,Red%ggf0:Red%ggf1)
            else
              p_JKQC(Red%ggf0:Red%ggf1) => &
                                         JKQC(iQ,K,Red%ggf0:Red%ggf1)
            end if

          end if ! Dynamic
        end if ! Coherent

        ! Initialize frequency indexes
        jjfreq = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Fed%nran
          do ifreq=Fed%if0(iran),Fed%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Manage MPI
            if (iifreq.lt.Mif0) cycle
            if (iifreq.gt.Mif1) exit

            ! Point to dimension
            p_mfreq => p_red%mfreq(iifreq)

            ! Initialize
            PRDin = cZero

            ! If coherent wing
            if (p_mfreq.lt.1) then

              ! Interpolate
              if (cohIn) then

                ! Input frequency
                omegai = omega(ifreq) - wlf

                ! Shift array beginning
                ishift = 1 - Red%ggf0

                ! Get interpolated JKQ
                y0 = getJKQinnu(omega(Red%ggf0:Red%ggf1), &
                                p_JKQC, &
                                ifreq+ishift, &
                                Atom%tif0(itran)+ishift, &
                                Atom%tif1(itran)+ishift, &
                                omegai)

                ! Fully coherent contribution
                if (conj) then
                  PRDin = sig*conjg(y0)
                else
                  PRDin = y0
                end if

              ! Just same frequency
              else

                ! Fully coherent contribution
                if (conj) then
                  PRDin = sig*conjg(p_JKQC(ifreq))
                else
                  PRDin = p_JKQC(ifreq)
                end if

              end if ! If need to interpolate

            ! Non-coherent
            else

              ! Compute norm
              Norme2 = sum(p_warr2(jjfreq+1:jjfreq+p_mfreq)* &
                           p_red%W_freq(jjfreq+1:jjfreq+p_mfreq))

              ! If valid norm
              if (abs(Norme2).gt.0d0) then

                ! If conjugate
                if (conj) then

                  ! Point to positive Q
                  p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,-iQ,K)

                  ! Integrate
                  PRDin = sum(conjg(p_JKQ)* &
                              p_red%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                              p_warr2(jjfreq+1:jjfreq+p_mfreq))*sig

                ! Not conjugate
                else

                  ! Point
                  p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,iQ,K)

                  ! Integrate
                  PRDin = sum(p_JKQ* &
                              p_red%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                              p_warr2(jjfreq+1:jjfreq+p_mfreq))

                end if


                ! Normalize to the first order profile (the
                ! normalization is completed later when multiplying
                ! by CRD)
                PRDin = PRDin/Norme2

              end if ! Valid norm

              ! Update index
              jjfreq = jjfreq + p_mfreq

            end if ! Coherent or non-coherent

            ! Substract the flat spectrum part due to just
            ! radiative excitation
            PRD(ifreq,1) = PRDin - Jrad(iQ,K)

          end do ! Output frequencies
        end do ! Output frequencies ranges

        ! Clean JKQ
        nullify(p_JKQ)
        nullify(p_JKQC)

      !
      ! Angle-dependent (Integral)
      !
      else

        ! Check if dynamic or with shift
        cohIn = lvel.or.abs(wlf).gt.0d0

        ! If asymmetries
        if (asym) then

          ! Check if conjugated
          conj = iQ.lt.0
          sig = Flgsg%sg(iQ)

          ! Just point
          if (conj) then
            p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(-iQ,K,:)
          else
            p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(iQ,K,:)
          end if
        end if ! Asymmetries

        ! For each output directions
        do jdir=1,Geom%njdir

          ! Initialize frequency indexes
          jjfreq = 0 ! Only for asymmetries
          jjfreq0 = 0
          llfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! Manage MPI
              if (iifreq.lt.Mif0) cycle
              if (iifreq.gt.Mif1) exit

              ! Initialize
              PRDin = cZero

              ! Point to dimension
              p_mfreq => p_red%mfreq(iifreq)

              ! Input direction index
              idir = 0

              ! For each polar direction
              do ith1=1,Geom%nTh

                ! For each azimuthal direction
                do iph1=1,Geom%nPh2

                  ! Advance index
                  idir = idir + 1

                  ! Scattering index
                  ish1 = Geom%i_scatt(iph1,ith1,jdir)

                  ! Special treatment if forward for two terms
                  if ((jtran.eq.itran.and. &
                       Geom%V_CScatt(ish1).ge.1d0).or. &
                      (p_mfreq.lt.1)) then

                    ! Interpolate
                    if (cohIn) then

                      ! Input frequency
                      omegai = omega(ifreq) - wlf

                      ! If axial
                      if (axial) then

                        ! If there are dynamics
                        if (lvel) then

                          ! Get director cosine
                          cost = Geom%V_mu(ith1)

                          ! Calculate Doppler shift factor
                          vfac1 = 1d0 - vz*cost

                          ! We will be using the inverse
                          vfac1 = 1d0/vfac1

                          ! Target frequency
                          omegai = omegai*vfac1

                        end if

                        ! Interpolate
                        StokesM = getStkinnu(omega, &
                                             Stokes(:,:,1,ith1), &
                                             ifreq, &
                                             Atom%tif0(itran), &
                                             Atom%tif1(itran), &
                                             omegai)

                      ! If non-axial
                      else

                        ! If there are dynamics
                        if (lvel) then

                          ! Get directional trigonometry
                          cost = Geom%V_mu(ith1)
                          sint = sqrt(1d0 - cost*cost)
                          cosc = Geom%v_mux(iph1)
                          sinc = Geom%v_muy(iph1)* &
                                 sqrt(1d0 - cosc*cosc)

                          ! Calculate Doppler shift factor
                          vfac1 = 1d0 - vx*sint*cosc - &
                                        vy*sint*sinc - &
                                        vz*cost

                          ! We will be using the inverse
                          vfac1 = 1d0/vfac1

                          ! Target frequency
                          omegai = omegai*vfac1

                        end if

                        ! Interpolate
                        StokesM = getStkinnu(omega, &
                                             Stokes(:,:,iph1,ith1), &
                                             ifreq, &
                                             Atom%tif0(itran), &
                                             Atom%tif1(itran), &
                                             omegai)

                      end if ! Axiality

                      ! Asymmetry
                      if (asym) then

                        ! Get interpolated JKQ
                        ishift = 1 - Red%ggf0
                        y0 = getJKQinnu(omega(Red%ggf0:Red%ggf1), &
                                        p_JKQC, &
                                        ifreq+ishift, &
                                        Atom%tif0(itran)+ishift, &
                                        Atom%tif1(itran)+ishift, &
                                        omegai)
                        ! Conjugate
                        if (conj) y0 = sig*conjg(y0)

                      end if ! Asymmetry

                    ! Fully coherent
                    else

                      ! Get value at frequency
                      if (axial) then
                        StokesM = Stokes(:,ifreq,1,ith1)
                      else
                        StokesM = Stokes(:,ifreq,iph1,ith1)
                      end if

                      ! Asymmetry
                      if (asym) then

                        ! Get value
                        y0 = p_JKQC(ifreq)

                        ! Conjugate
                        if (conj) y0 = sig*conjg(y0)

                      end if ! Asymmetry

                    end if ! Full coherence

                    ! Sum over Stokes parameters
                    intgr = sum(Geom%TS(:,iQ,K,idir)* &
                                StokesM)

                    ! Asymmetry
                    if (asym) intgr = intgr + y0

                    ! Add to integral
                    PRDin = PRDin + intgr*Geom%W_mu(ith1)* &
                                          Geom%W_mux2(iph1)

                  ! Non-forward 2-term scattering
                  else

                    ! Shift in indexes
                    kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                    ! Multiply Warr2 and weights
                    Warr2xW(1:p_mfreq) = &
                          p_warr2(kkfreq0b+1:kkfreq0b+p_mfreq)* &
                          p_red%W_freq(llfreq0+1:llfreq0+p_mfreq)

                    ! Compute norm denominator
                    Norme2 = sum(Warr2xW(1:p_mfreq))

                    ! If valid norm
                    if (abs(Norme2).gt.0d0) then

                      ! Sum Stokes
                      intergrin(1:p_mfreq) = &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,0)* &
                              Geom%TS(0,iQ,K,idir) + &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,1)* &
                              Geom%TS(1,iQ,K,idir) + &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,2)* &
                              Geom%TS(2,iQ,K,idir) + &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,3)* &
                              Geom%TS(3,iQ,K,idir)

                      ! If asymmetry
                      if (asym) then

                        ! Conjugate
                        if (conj) then

                          ! Point to positive Q
                          p_JKQ => &
                                  JradC(jjfreq+1:jjfreq+p_mfreq,-iQ,K)

                          ! Add
                          intergrin(1:p_mfreq) = &
                                      intergrin(1:p_mfreq) + &
                                      sig*conjg(p_JKQ)

                        ! Not conjugate
                        else

                          ! Point to positive Q
                          p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,iQ,K)

                          ! Add
                          intergrin(1:p_mfreq) = &
                                          intergrin(1:p_mfreq) + p_JKQ

                        end if ! Conjugate
                      end if ! Asymmetry

                      ! Normalize to the first order profile
                      ! and add the directional weights (the
                      ! normalization will be completed later
                      ! multiplying by CRD)
                      PRDin = PRDin + &
                              sum(Warr2xW(1:p_mfreq)* &
                                  intergrin(1:p_mfreq))* &
                              Geom%W_mu(ith1)* &
                              Geom%W_mux2(iph1)/Norme2

                    end if ! Valid norm
                  end if ! Type of scattering

                  ! Advance index if not axial
                  if (.not.axial) jjfreq0 = jjfreq0 + p_mfreq

                end do ! azimuthal nodes

                ! Update index if axial
                if (axial) jjfreq0 = jjfreq0 + p_mfreq

              end do ! polar nodes

              ! Update llfreq and kkfreq
              llfreq0 = llfreq0 + p_mfreq
              kkfreq0 = kkfreq0 + p_mfreq*(Geom%nScatt-nfs)
              jjfreq = jjfreq + p_mfreq ! Asymmetries

              ! Subtract the flat spectrum part due to just
              ! radiative excitation
              PRD(ifreq,jdir) = PRDin - Jrad(iQ,K)

            end do ! output frequencies
          end do ! output frequencies ranges
        end do ! Output directions

      end if ! AA/AD (Integral)

      ! Clean p_warr2
      if (LPRAM.and.nmfreq.gt.0) deallocate(p_warr2)
      nullify(p_warr2)

      ! For each K'
      do K1=0,2

        ! Get real value
        rK1 = dble(K1)

        ! Racah algebra
        f1tmp = fun6j(rK1,rJu,rJu1,rJf,1d0,1d0,Flgsg)

        ! If forbidden (6j-sym=0), skip
        if(abs(f1tmp).lt.TINYJS) cycle

        f1tmp = sqrt(2d0*rK1+1d0)*f1tmp

        ! For each Q'
        do iQ1=-K1,K1

          ! Get real value
          Q1 = dble(iQ1)

          ! Get Ql
          iQl = iQ + iQ1
          Ql = dble(iQl)

          ! Determine the limits in K
          Kmin = max(abs(iQl),nint(abs(rK-rK1)), &
                     nint(abs(rJl-rJl1)))
          Kmax = min(nint(rJl+rJl1),nint(rK+rK1),Kcut)

          ! Initialize sum
          rhoc = cZero

          ! For each Kl
          do Kl=Kmin,Kmax

            ! Check sum K + Kl
            if ((Kl+K).gt.Kcut) cycle

            ! Get the real number
            rKl = dble(Kl)

            ! Get the SEE index
            iR = Atom%irho(iterml)%Jrho(iL1,iL)%kq(iQl,Kl)

            ! If flagged as small, skip
            if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

            ! Racah algebra
            fKfj = sqrt(2d0*rKl+1d0)*Flgsg%sg(Kl-iQl)* &
                   fun3j(rK,rK1,rKl,Q,Q1,-Ql,Flgsg)* &
                   fun9j(rK,rK1,rKl,1d0,rJu1,rJl1, &
                         1d0,rJu,rJl,Flgsg)

            ! Small J symbols
            if (abs(fKfj).lt.TINYJS) cycle

            ! Accumulate in the sum
            rhoc = fKfj*Atom%crho(iR,iz) + rhoc

          end do ! Kl

          ! If no rhoKQ, skip
          if (abs(rhoc).lt.TINYER) cycle

          ! Complete sum
          rhoc = rhoc*ftmp*f1tmp

          !
          ! Angle-average (tmp)
          !
          if (AV) then

            ! For each output frequency
            iifreq = 0
            do iran=1,Fed%nran

              ! Limits in range
              i0 = iifreq + 1
              i1 = i0 + Fed%if1(iran) - Fed%if0(iran)

              ! Update iifreq
              iifreq = i1

              ! Check MPI
              if (i1.lt.Mif0) cycle
              if (i0.gt.Mif1) exit

              ! Get limits
              if (Mif0.gt.i0) then
                i0 = Fed%if0(iran) + Mif0 - i0
              else
                i0 = Fed%if0(iran)
              end if
              if (Mif1.lt.i1) then
                i1 = Fed%if1(iran) + Mif1 - i1
              else
                i1 = Fed%if1(iran)
              end if

              ! Auxiliar (complete normalization)
              PRDaux(i0:i1) = PRD(i0:i1,1)*CRD(i0:i1)*rhoc

              ! Output directions
              do jdir=1,Geom%njdir

                ! Add to tmp
                tmp0(i0:i1,jdir) = tmp0(i0:i1,jdir) + &
                                   TSout(0,iQ1,K1,jdir)*PRDaux(i0:i1)
                tmp1(i0:i1,jdir) = tmp1(i0:i1,jdir) + &
                                   TSout(1,iQ1,K1,jdir)*PRDaux(i0:i1)
                tmp2(i0:i1,jdir) = tmp2(i0:i1,jdir) + &
                                   TSout(2,iQ1,K1,jdir)*PRDaux(i0:i1)
                tmp3(i0:i1,jdir) = tmp3(i0:i1,jdir) + &
                                   TSout(3,iQ1,K1,jdir)*PRDaux(i0:i1)

              end do ! Output directions
            end do ! Output frequency ranges

          !
          ! Angle-dependent (tmp)
          !
          else

            ! For each output frequency
            iifreq = 0
            do iran=1,Fed%nran

              ! Limits in range
              i0 = iifreq + 1
              i1 = i0 + Fed%if1(iran) - Fed%if0(iran)

              ! Update iifreq
              iifreq = i1

              ! Check MPI
              if (i1.lt.Mif0) cycle
              if (i0.gt.Mif1) exit

              ! Get limits
              if (Mif0.gt.i0) then
                i0 = Fed%if0(iran) + Mif0 - i0
              else
                i0 = Fed%if0(iran)
              end if
              if (Mif1.lt.i1) then
                i1 = Fed%if1(iran) + Mif1 - i1
              else
                i1 = Fed%if1(iran)
              end if

              ! Output directions
              do jdir=1,Geom%njdir

                ! Auxiliar (complete normalization)
                PRDaux(i0:i1) = PRD(i0:i1,jdir)*CRD(i0:i1)*rhoc

                ! Add to tmp
                tmp0(i0:i1,jdir) = tmp0(i0:i1,jdir) + &
                                   TSout(0,iQ1,K1,jdir)*PRDaux(i0:i1)
                tmp1(i0:i1,jdir) = tmp1(i0:i1,jdir) + &
                                   TSout(1,iQ1,K1,jdir)*PRDaux(i0:i1)
                tmp2(i0:i1,jdir) = tmp2(i0:i1,jdir) + &
                                   TSout(2,iQ1,K1,jdir)*PRDaux(i0:i1)
                tmp3(i0:i1,jdir) = tmp3(i0:i1,jdir) + &
                                   TSout(3,iQ1,K1,jdir)*PRDaux(i0:i1)

              end do ! Output directions
            end do ! Output frequency ranges

          end if ! AA/AD (tmp)

        end do !Q'
      end do !K'

            !
            ! Recover indentation
            !

          end do ! Q
        end do ! K

        ! Update redistribution initial index
        kwfreq0 = kwfreq0 + nmfreq

        ! Apply hanle factor and Einstein coefficient
        daux = (2d0*rLl+1d0)*Atom%Ecoeff(iterml,itermu)
        eps20 = eps20 + dble(transpose(tmp0)/hanleden)*daux
        eps21 = eps21 + dble(transpose(tmp1)/hanleden)*daux
        eps22 = eps22 + dble(transpose(tmp2)/hanleden)*daux
        eps23 = eps23 + dble(transpose(tmp3)/hanleden)*daux

                  !
                  ! Recover indentation
                  !

                end do ! iL1
              end do ! iL
            end do ! iU1
          end do ! iU
        end do ! mF
      end do ! Terms

      ! Apply common factor
      daux = 3d0*.5d0*IPI42*(2d0*rLu+1d0)* &
             Atom%Ecoeff(itermu,itermf)*1d-10/(c*Dw)
      eps20 = eps20*daux
      eps21 = eps21*daux
      eps22 = eps22*daux
      eps23 = eps23*daux

      ! Clean pointers
      if (associated(p_red)) nullify(p_red)
      if (associated(p_rwarr)) nullify(p_rwarr)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_warr2)) nullify(p_warr2)
      if (associated(p_JKQ)) nullify(p_JKQ)
      if (associated(p_JKQC)) nullify(p_JKQC)

      return

      end subroutine emiss2ordNB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity and emissivity of a given atomic
      !! line in the presence of a magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  TB(dcomplx(:,:,:)): Geometrical tensors in the magnetic
      !!                      field reference frame\n
      !!    omega(double(:)): Frequency array\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!      itran(integer): Index of transition to compute\n
      !!     itermu(integer): Upper term of the transition\n
      !!     iterml(integer): Lower term of the transition\n
      !!         iz(integer): Height index\n
      !!        if0(integer): First frequency index for this
      !!                      transition\n
      !!        if1(integer): Last frequency index for this
      !!                      transition\n
      !!   Norma(Prof_class): Normalization factors for Voigt
      !!                      profiles or Voigt profiles\n
      !!          Dw(double): Doppler width of the transition\n
      !!        vfac(double): Doppler shift factor\n
      !!        absK(double): Unit transformation factor\n
      !!     eta0(double(:)): Intensity absorptivity\n
      !!     eta1(double(:)): Q absorptivity\n
      !!     eta2(double(:)): U absorptivity\n
      !!     eta3(double(:)): V absorptivity\n
      !!     rha1(double(:)): Q dichroic absorptivity\n
      !!     rha2(double(:)): U dichroic absorptivity\n
      !!     rha3(double(:)): V dichroic absorptivity\n
      !!     eps0(double(:)): Intensity emissivity\n
      !!     eps1(double(:)): Q emissivity\n
      !!     eps2(double(:)): U emissivity\n
      !!     eps3(double(:)): V emissivity\n
      !!     rhs1(double(:)): Q 'dichroic' emissivity\n
      !!     rhs2(double(:)): U 'dichroic' emissivity\n
      !!     rhs3(double(:)): V 'dichroic' emissivity
      subroutine rt1ord(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                        iz,if0,if1,Norma,Dw,vfac,absK, &
                        eta0,eta1,eta2,eta3,rha1,rha2,rha3, &
                        eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iz,if0,if1
      double precision, intent(in):: Dw,absK,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      logical:: lNCHLT

      integer:: ifreq,K,iq,iq1,iQQ,indL,indU,indK
      integer:: nMu,nMl,iMu,iMl,iMu1,iMl1,iJl1,iJlb,iJu1,iJub
      integer:: iU,iU1,kU1,kUb,iL,iL1,kL1,kLb

      double precision:: rLu,rLl,S,rJu1,rJub,rJl1,rJlb
      double precision:: rJumax,rJlmax,rMu,rMu1,rMl,rMl1
      double precision:: eu,el,rK,QQ,q,q1,au,al,aul,ftmp
      double precision:: Cu1,Cub,Cl1,Clb,EVul,EVu1l,EVul1
      double precision:: iDw,at,Dfreq,vfacw,tempRe,tempRa

      complex(kind=8):: tK,rhoc
      complex(kind=8), dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

      eta0 = 0d0
      eta1 = 0d0
      eta2 = 0d0
      eta3 = 0d0
      rha1 = 0d0
      rha2 = 0d0
      rha3 = 0d0
      eps0 = 0d0
      eps1 = 0d0
      eps2 = 0d0
      eps3 = 0d0
      rhs1 = 0d0
      rhs2 = 0d0
      rhs3 = 0d0

      ! If to consider non-coherent lower term and this atom has the
      ! data allocated
      if (NCHLT.and.allocated(Atom%NCHLT)) then

        ! Check if this transition has non-coherent lower term for
        ! this height
        lNCHLT = Atom%NCHLT(iz,itran)

      ! Coherent lower term
      else

        ! Flag non-coherent as false
        lNCHLT = .False.

      end if ! non-coherent lower term

      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)*iDw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Determine the maximum angular momentum and the number
      ! of magnetic sublevels for that maximum momentum
      rJumax = rLu+S
      nMu = nint(2d0*rJumax+1d0)
      rJlmax = rLl + S
      nMl = nint(2d0*rJlmax+1d0)

      ! Doppler shift in doppler units
      vfacw = vfac*iDw


      !
      ! Common part
      !

      ! For each Ml
      do iMl=1,nMl

        ! Value of Ml
        rMl = -rJlmax + dble(iMl-1)

        ! For each mu_l
        do iL=1,Atom%nblk(iMl,iterml)

          ! Get eigenvalue lower level
          el = Atom%eval(iL,iMl,iterml,iz)

          ! Level index
          indL = Atom%irho(iterml)%jM(iL,iMl)

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! If not pi nor sigma, skip
            if (abs(rMu-rMl).gt.1.1d0) cycle

            ! Get difference between M momenta in integer
            q = rMu-rMl
            iq = nint(q)

            ! For each mu_u
            do iU=1,Atom%nblk(iMu,itermu) ! sum over mu_u

              ! Dipole strength
              EVul = Atom%rdipev(iz)%rdipev(itran)% &
                          rdip(iq,iU,iL,iMu,iMl)

              ! Check if small
              if (abs(EVul).lt.TINYEV) cycle

              ! Get eigenvalue upper level
              eu = Atom%eval(iU,iMu,itermu,iz)

              ! Get indexes
              indU = Atom%irho(itermu)%jM(iU,iMu)
              indK = Atom%trano(itran)%indB(indL,indU)

              !
              ! Compute profile
              !

              ! If stored
              if (Norma%VRAM) then

                ! Copy stored profile
                prof = Norma%cp(:,indK)

              ! Not stored
              else

                ! Shift term
                Dfreq = (eu - el + Atom%Dfreq(itran))*iDw

                ! For each frequency
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies

                ! Normalize
                prof = dcmplx(dble(prof)*Norma%Norm(indK), &
                              dimag(prof))

              end if ! Storing

              ! For each possible K
              do K=0,2

                ! Get the real number
                rK = dble(K)

                ! For each Q
                do iQQ=-K,K

                  ! Get q' index
                  iq1 = iQQ + iq

                  ! Check selection rules
                  if (abs(iq1).gt.1) cycle

                  ! Get q' and Q' values
                  QQ = dble(iQQ)
                  q1 = dble(iq1)

                  ! Racah algebra
                  ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

                  ! If not allowed (3J-sym=0) skip
                  if (abs(ftmp).lt.TINYJS) cycle

                  ! Complete the factor
                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)*EVul

                  !
                  ! Reset indentation
                  !

      !
      ! Emission
      !

      ! Initialize tK
      tK = cZero

      ! For each Mu'
      do iMu1=1,nMu

        ! Value of Mu'
        rMu1 = -rJumax + dble(iMu1-1)

        ! Check valid J symbols
        if (nint(rMu1-rMl).ne.iq1) cycle

        ! For each mu_u'
        do iU1=1,Atom%nblk(iMu1,itermu)

          ! Get dipole strength
          EVu1l = Atom%rdipev(iz)%rdipev(itran)% &
                          rdip(iq1,iU1,iL,iMu1,iMl)

          ! Check if small
          if (abs(EVu1l).lt.TINYEV) cycle

          ! For each Ju'
          do kU1=1,Atom%nblk(iMu1,itermu)

            ! Get eigenvector
            Cu1 = Atom%evec(kU1,iU1,iMu1,itermu,iz)

            ! If coefficient too small, skip
            if (abs(Cu1).lt.TINYEV) cycle

            ! Get J level index
            iJu1 = Atom%iJval(kU1,iMu1,itermu)

            ! Get angular momentum
            rJu1 = Atom%rJval(iJu1,itermu)

            ! Add sign
            Cu1 = Cu1*Flgsg%sg(nint(rJu1 - rMu1))

            ! For each Jub
            do kUb=1,Atom%nblk(iMu,itermu) ! sum Jub

              ! Get eigenvector for upper level b
              Cub = Atom%evec(kUb,iU,iMu,itermu,iz)

              ! If coefficient too small, skip
              if (abs(Cub).lt.TINYEV) cycle

              ! Get J level index
              iJub = Atom%iJval(kUb,iMu,itermu)

              ! Get angular momentum
              rJub = Atom%rJval(iJub,itermu)

              ! Sum over (Kl,Ql)
              call emiss1(Atom,Flgsg,iz,itermu,iJu1, &
                          iJub,rJu1,rJub,rMu1,rMu,rhoc)

              ! Uncomment the following line for Zeeman
             !if (iJub.ne.iJu1) rhoc=cZero

              ! If no 'population', skip
              if (abs(rhoc).lt.TINYER) cycle

              ! Accumulate into tK
              tK = ftmp*EVu1l*Cu1*Cub*rhoc + tK

            end do ! kUb
          end do ! kU1
        end do ! iU1
      end do ! iMu1

      ! Add the profile
      profK = prof*tK

      ! Emissivity
      eps0 = eps0 + dble(TB(0,iQQ,K)*profK)
      eps1 = eps1 + dble(TB(1,iQQ,K)*profK)
      eps2 = eps2 + dble(TB(2,iQQ,K)*profK)
      eps3 = eps3 + dble(TB(3,iQQ,K)*profK)

      ! Dispersion
      rhs1 = rhs1 + dimag(TB(1,iQQ,K)*profK)
      rhs2 = rhs2 + dimag(TB(2,iQQ,K)*profK)
      rhs3 = rhs3 + dimag(TB(3,iQQ,K)*profK)

      !
      ! Absorption
      !

      ! Initialize tK
      tK = cZero

      ! For each Ml'
      do iMl1=1,nMl

        ! If non-coherent lower term, skip if not diagonal
        if (lNCHLT) then
          if (iMl1.ne.iMl) cycle
        end if

        ! Value of Ml'
        rMl1 = -rJlmax + dble(iMl1-1)

        ! Check valid J symbols
        if (nint(rMu-rMl1).ne.iq1) cycle

        ! For each mu_l'
        do iL1=1,Atom%nblk(iMl1,iterml)

          ! Get dipole strength
          EVul1 = Atom%rdipev(iz)%rdipev(itran)% &
                          rdip(iq1,iU,iL1,iMu,iMl1)

          ! Check if small
          if (abs(EVul1).lt.TINYEV) cycle

          ! For each Jl'
          do kL1=1,Atom%nblk(iMl1,iterml)

            ! Eigenvector
            Cl1 = Atom%evec(kL1,iL1,iMl1,iterml,iz)

            ! If coefficient too small, skip
            if (abs(Cl1).lt.TINYEV) cycle

            ! Get J level index
            iJl1 = Atom%iJval(kL1,iMl1,iterml)

            ! Get angular momentum
            rJl1 = Atom%rJval(iJl1,iterml)

            ! For each Jlb
            do kLb=1,Atom%nblk(iMl,iterml)

              ! Get eigenvector for lower level b
              Clb = Atom%evec(kLb,iL,iMl,iterml,iz)

              ! If coefficient too small, skip
              if (abs(Clb).lt.TINYEV) cycle

              ! Get J level index
              iJlb = Atom%iJval(kLb,iMl,iterml)

              ! Get angular momentum
              rJlb = Atom%rJval(iJlb,iterml)

              ! Sum over (Kl,Ql)
              call absorb1(Atom,Flgsg,iz,iterml,iJlb, &
                           iJl1,rJlb,rJl1,rMl,rMl1,0, &
                           rhoc)

              ! If no population, skip
              if (abs(rhoc).lt.TINYER) cycle

              ! Accumulate into tK
              tK = tK + &
                   ftmp*EVul1*Cl1*Clb*rhoc*Flgsg%sg(nint(rJlb-rMl))

            end do ! kLb
          end do ! kL1
        end do ! iL1
      end do ! iMl1

      ! Add the profile
      profK = prof*tk

      ! Absorptivity
      eta0 = eta0 + dble(TB(0,iQQ,K)*profk)
      eta1 = eta1 + dble(TB(1,iQQ,K)*profk)
      eta2 = eta2 + dble(TB(2,iQQ,K)*profk)
      eta3 = eta3 + dble(TB(3,iQQ,K)*profk)

      ! Dispersion
      rha1 = rha1 + dimag(TB(1,iQQ,K)*profk)
      rha2 = rha2 + dimag(TB(2,iQQ,K)*profk)
      rha3 = rha3 + dimag(TB(3,iQQ,K)*profk)

                  ! Restore identation

                end do ! iQQ
              end do ! K
            end do ! iU
          end do ! iMu
        end do ! iL
      end do ! iMl

      ! Common parts for coefficients
      tempRe = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
               Atom%Ecoeff(itermu,iterml)/Dw
      tempRa = tempRe/absK

      ! Final values
      eps0 = tempRe*eps0
      eps1 = tempRe*eps1
      eps2 = tempRe*eps2
      eps3 = tempRe*eps3
      rhs1 = tempRe*rhs1
      rhs2 = tempRe*rhs2
      rhs3 = tempRe*rhs3
      eta0 = tempRa*eta0
      eta1 = tempRa*eta1
      eta2 = tempRa*eta2
      eta3 = tempRa*eta3
      rha1 = tempRa*rha1
      rha2 = tempRa*rha2
      rha3 = tempRa*rha3

      end subroutine rt1ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient.\n
      !> Calculate the absorptivity of a given atomic line in the
      !! presence of a magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  TB(dcomplx(:,:,:)): Geometrical tensors in the magnetic
      !!                      field reference frame\n
      !!    omega(double(:)): Frequency array\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!      itran(integer): Index of transition to compute\n
      !!     itermu(integer): Upper term of the transition\n
      !!     iterml(integer): Lower term of the transition\n
      !!         iz(integer): Height index\n
      !!        if0(integer): First frequency index for this
      !!                      transition\n
      !!        if1(integer): Last frequency index for this
      !!                      transition\n
      !!   Norma(Prof_class): Normalization factors for Voigt
      !!                      profiles or Voigt profiles\n
      !!          Dw(double): Doppler width of the transition\n
      !!        vfac(double): Doppler shift factor\n
      !!        absK(double): Unit transformation factor\n
      !!     eta0(double(:)): Intensity absorptivity\n
      !!     eta1(double(:)): Q absorptivity\n
      !!     eta2(double(:)): U absorptivity\n
      !!     eta3(double(:)): V absorptivity\n
      !!     rha1(double(:)): Q dichroic absorptivity\n
      !!     rha2(double(:)): U dichroic absorptivity\n
      !!     rha3(double(:)): V dichroic absorptivity
      subroutine absorb(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                        iz,if0,if1,Norma,Dw,vfac,absK, &
                        eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iz,if0,if1
      double precision, intent(in):: Dw,absK,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      logical:: lNCHLT

      integer:: ifreq,K,iq,iq1,iQQ,nMu,nMl,iMu,iMl,iMl1,iU,iL
      integer:: indU,indL,indK

      double precision:: rLu,rLl,S,rJumax,rJlmax,rMu,rMl,rMl1
      double precision:: eu,el,rK,QQ,q,q1,au,al,aul,ftmp,tempR
      double precision:: iDw,at,Dfreq,vfacw

      complex(kind=8):: tK
      complex(kind=8), dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

      eta0 = 0d0
      eta1 = 0d0
      eta2 = 0d0
      eta3 = 0d0
      rha1 = 0d0
      rha2 = 0d0
      rha3 = 0d0

      ! If to consider non-coherent lower term and this atom has the
      ! data allocated
      if (NCHLT.and.allocated(Atom%NCHLT)) then

        ! Check if this transition has non-coherent lower term for
        ! this height
        lNCHLT = Atom%NCHLT(iz,itran)

      ! Coherent lower term
      else

        ! Flag non-coherent as false
        lNCHLT = .False.

      end if ! non-coherent lower term

      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)*iDw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Determine the maximum angular momentum and the number
      ! of magnetic sublevels for that maximum momentum
      rJumax = rLu+S
      nMu = nint(2d0*rJumax+1d0)
      rJlmax = rLl+S
      nMl = nint(2d0*rJlmax+1d0)

      ! Doppler shift in doppler units
      vfacw = vfac*iDw


      !
      ! Compute absorptivity
      !

      ! For each Mu
      do iMu=1,nMu

        ! Value of Mu
        rMu = -rJumax + dble(iMu-1)

        ! For each mu_u
        do iU=1,Atom%nblk(iMu,itermu)

          ! Get eigenvalue upper level
          eu = Atom%eval(iU,iMu,itermu,iz)

          ! Level index
          indU = Atom%irho(itermu)%jM(iU,iMu)

          ! For each Ml
          do iMl=1,nMl

            ! Value of Ml
            rMl = -rJlmax + dble(iMl-1)

            ! If not pi nor sigma, skip
            if (nint(abs(rMu-rMl)).gt.1) cycle

            ! Get difference between M momentums in integer
            q = rMu-rMl
            iq = nint(q)

            ! For each mu_l
            do iL=1,Atom%nblk(iMl,iterml)

              ! Get eigenvalue lower level
              el = Atom%eval(iL,iMl,iterml,iz)

              ! Get indexes
              indL = Atom%irho(iterml)%jM(iL,iMl)
              indK = Atom%trano(itran)%indB(indL,indU)

              !
              ! Compute profile
              !

              ! If stored
              if (Norma%VRAM) then

                ! Copy stored profile
                prof = Norma%cp(:,indK)

              ! Not stored
              else

                ! Shift term
                Dfreq = (eu - el + Atom%Dfreq(itran))*iDw

                ! For each frequency
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies

                ! Normalize profile
                prof = dcmplx(dble(prof)*Norma%Norm(indK), &
                              dimag(prof))

              end if ! Storing

              ! For each Ml'
              do iMl1=1,nMl

                ! If non-coherent lower term, skip if not diagonal
                if (lNCHLT) then
                  if (iMl1.ne.iMl) cycle
                end if

                ! Value of Ml'
                rMl1 = -rJlmax + dble(iMl1-1)

                ! If not pi nor sigma, skip
                if (nint(abs(rMu-rMl1)).gt.1) cycle

                ! Get the difference between M momentums
                q1 = rMu - rMl1
                QQ = q1-q

                ! Make the difference integers
                iq1 = nint(q1)
                iQQ = nint(QQ)

                ! For each K
                do K=abs(iQQ),2

                  ! Get the real number
                  rK = dble(K)

                  ! Racah algebra
                  ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

                  ! If not allowed (3J-sym=0) skip
                  if (abs(ftmp).lt.TINYJS) cycle

                  ! Complete factor
                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)

                  ! Run inner atomic loop
                  call abs_inner_atomic_loop(Atom,Flgsg, &
                                             itermu,iMu,iU, &
                                             iterml,iMl,iL,iMl1, &
                                             iz,itran,iq,iq1, &
                                             rMl,rMl1,ftmp,tk)

                  ! Add the profile
                  profK = prof*tk

                  ! Absorptivity
                  eta0 = eta0 + dble(TB(0,iQQ,K)*profk)
                  eta1 = eta1 + dble(TB(1,iQQ,K)*profk)
                  eta2 = eta2 + dble(TB(2,iQQ,K)*profk)
                  eta3 = eta3 + dble(TB(3,iQQ,K)*profk)

                  ! Dispersion
                  rha1 = rha1 + dimag(TB(1,iQQ,K)*profk)
                  rha2 = rha2 + dimag(TB(2,iQQ,K)*profk)
                  rha3 = rha3 + dimag(TB(3,iQQ,K)*profk)

                end do ! K
              end do ! Ml1
            end do ! iL
          end do ! Ml
        end do ! iU
      end do ! Mu

      ! Multiplicative constant
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)*iDw/absK

      ! Final values
      eta0 = tempR*eta0
      eta1 = tempR*eta1
      eta2 = tempR*eta2
      eta3 = tempR*eta3
      rha1 = tempR*rha1
      rha2 = tempR*rha2
      rha3 = tempR*rha3

      return

      end subroutine absorb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity of a given atomic line in the
      !! presence of a magnetic field. This one is used to get the
      !! emissivity of a PRD line and thus calculates the emissivity
      !! in the comoving frame for all directions\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!  TB(dcomplx(:,:,:,:)): Geometrical tensors in the magnetic
      !!                        field reference frame\n
      !!      omega(double(:)): Frequency array\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!        itran(integer): Index of transition to compute\n
      !!       itermu(integer): Upper term of the transition\n
      !!       iterml(integer): Lower term of the transition\n
      !!           iz(integer): Height index\n
      !!          if0(integer): First frequency index for this
      !!                        transition\n
      !!          if1(integer): Last frequency index for this
      !!                        transition\n
      !!        njdir(integer): Number of directions\n
      !!     Norma(Prof_class): Normalization factors for Voigt
      !!                        profiles or Voigt profiles\n
      !!            Dw(double): Doppler width of the transition\n
      !!       eps0(double(:)): Intensity emissivity\n
      !!       eps1(double(:)): Q emissivity\n
      !!       eps2(double(:)): U emissivity\n
      !!       eps3(double(:)): V emissivity
      subroutine emiss(Atom,TB,omega,Flgsg,itran,itermu, &
                       iterml,iz,if0,if1,njdir,Norma,Dw, &
                       eps0,eps1,eps2,eps3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iz,if0,if1,njdir
      double precision, intent(in):: Dw
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(njdir,if0:if1), intent(out):: eps0
      double precision, dimension(njdir,if0:if1), intent(out):: eps1
      double precision, dimension(njdir,if0:if1), intent(out):: eps2
      double precision, dimension(njdir,if0:if1), intent(out):: eps3
      complex(kind=8), dimension(0:3,-2:2,0:2,njdir), intent(in):: TB

      ! Local

      integer:: ifreq,K,iq,iq1,iQQ,jdir,indU,indL,indK
      integer:: nMl,nMu,iMl,iMu,iMu1,iL,iU

      double precision:: rLl,rLu,S,rJlmax,rJumax,rMl,rMu,rMu1,el,eu
      double precision:: rK,QQ,q,q1,al,au,aul,ftmp,tempR,iDw,at,Dfreq

      complex(kind=8):: tK
      complex(kind=8), dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

      eps0 = 0d0
      eps1 = 0d0
      eps2 = 0d0
      eps3 = 0d0

      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)*iDw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Determine the maximum angular momentum and the number
      ! of magnetic sublevels for that maximum momentum
      rJumax = rLu+S
      nMu = nint(2d0*rJumax+1d0)
      rJlmax = rLl + S
      nMl = nint(2d0*rJlmax+1d0)

      ! Doppler shift in doppler units
      iDw = 1d0/Dw


      !
      ! Compute emissivity
      !

      ! For each Ml
      do iMl=1,nMl

        ! Value of Ml
        rMl = -rJlmax + dble(iMl-1)

        ! For each mu_l
        do iL=1,Atom%nblk(iMl,iterml)

          ! Get eigenvalue lower level
          el = Atom%eval(iL,iMl,iterml,iz)/Dw

          ! Level index
          indL = Atom%irho(iterml)%jM(iL,iMl)

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! If not pi nor sigma, skip
            if (nint(abs(rMu-rMl)).gt.1) cycle

            ! Get difference between M momentums in integer
            q = rMu-rMl
            iq = nint(q)

            ! For each mu_u
            do iU=1,Atom%nblk(iMu,itermu) ! sum over mu_u

              ! Get eigenvalue upper level
              eu = Atom%eval(iU,iMu,itermu,iz)/Dw

              ! Get indexes
              indU = Atom%irho(itermu)%jM(iU,iMu)
              indK = Atom%trano(itran)%indB(indL,indU)

              !
              ! Compute profile
              !

              ! If stored
              if (Norma%VRAM) then

                ! Copy stored profile
                prof = Norma%cp(:,indK)

              ! Not stored
              else

                ! Shift term
                Dfreq = eu - el + Atom%Dfreq(itran)

                ! For each frequency
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt((Dfreq - omega(ifreq))*iDw,at, &
                             prof(ifreq))

                end do ! frequencies

                ! Normalize
                prof = dcmplx(dble(prof)*Norma%Norm(indK), &
                              dimag(prof))

              end if ! Storing

              ! For each Mu'
              do iMu1=1,nMu

                ! Value of Mu'
                rMu1 = -rJumax + dble(iMu1-1)

                ! If not pi nor sigma, skip
                if (nint(abs(rMu1-rMl)).gt.1) cycle

                ! Get the difference between M momentums
                q1 = rMu1-rMl
                QQ = q1-q

                ! Make the difference integers
                iq1 = nint(q1)
                iQQ = nint(QQ)

                ! For each K
                do K=abs(iQQ),2

                  ! Get the real number
                  rK = dble(K)

                  ! Racah algebra
                  ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

                  ! If not allowed (3j-sym=0) skip
                  if (abs(ftmp).lt.TINYJS) cycle

                  ! Complete factor
                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)

                  ! Call inner atomic loop
                  call emi_inner_atomic_loop(Atom,Flgsg, &
                                             itermu,iMu,iU,iMu1, &
                                             iterml,iMl,iL, &
                                             iz,itran, &
                                             iq,iq1, &
                                             rMu,rMu1, &
                                             ftmp,tk)

                  ! Add the profile
                  profK = prof*tK

                  ! For each direction
                  do jdir=1,njdir

                    ! Emissivity
                    eps0(jdir,:) = eps0(jdir,:) + &
                                   dble(TB(0,iQQ,K,jdir)*profK)
                    eps1(jdir,:) = eps1(jdir,:) + &
                                   dble(TB(1,iQQ,K,jdir)*profK)
                    eps2(jdir,:) = eps2(jdir,:) + &
                                   dble(TB(2,iQQ,K,jdir)*profK)
                    eps3(jdir,:) = eps3(jdir,:) + &
                                   dble(TB(3,iQQ,K,jdir)*profK)

                  end do ! Directions
                end do ! K
              end do ! Mu1
            end do ! iU
          end do ! Mu
        end do ! iL
      end do ! Ml

      ! Constant factor
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)*iDw

      ! Final values
      eps0 = tempR*eps0
      eps1 = tempR*eps1
      eps2 = tempR*eps2
      eps3 = tempR*eps3

      return

      end subroutine emiss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the second order emissivity of a given atomic line
      !! in the presence of a magnetic field. This subroutine only
      !! computes the positive contribution of the coherent scattering
      !! and the negative contribution for the flat-spectrum, i.e.,
      !! the result needs to be added to the product of emiss to get
      !! the actual emissivity\n
      !!         Atom(Atom_class): Structure with atomic data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!               vx(double): Velocity vector along X\n
      !!               vy(double): Velocity vector along Y\n
      !!               vz(double): Velocity vector along Z\n
      !!            lvel(logical): If dynamic node\n
      !!         omega(double(:)): Frequency array\n
      !!          Fed(Reda_class): Structure with redistribution
      !!                           output frequency data\n
      !!          Red(Redb_class): Structure with redistribution input
      !!                           frequency data\n
      !!       RWarr(Redb2_class): Structure with redistribution
      !!                           function data\n
      !!        Norma(Prof_class): Normalization factors for Voigt
      !!                           profiles or Voigt profiles\n
      !!       Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                           and J-symbols\n
      !!           jtran(integer): Output transition index\n
      !!          itermu(integer): Upper term of output transition\n
      !!          itermf(integer): Lower term of output transition\n
      !!              iz(integer): Height index\n
      !!             if0(integer): First frequency index for this
      !!                           transition\n
      !!             if1(integer): Last frequency index for this
      !!                           transition\n
      !!            Mif0(integer): First frequency for this CPU\n
      !!            Mif1(integer): Last frequency for this CPU\n
      !!              DwT(double): Thermal part of Doppler width\n
      !!               Dw(double): Doppler width of the output
      !!                           transition\n
      !!     Bfield(Bfield_class): Structure with magnetic field
      !!                           data\n
      !!              vmi(double): Microturbulent velocity\n
      !!  TBout(docmplx(:,:,:,:)): Geometrical tensor in the magnetic
      !!                           field reference frame\n
      !!  Stokes(double(:,:,:,:)): Stokes parameters\n
      !!    JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!     JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                           over the absorption profile\n
      !!    JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                           frequency dependence\n
      !!         eps20(double(:)): Intensity emissivity\n
      !!         eps21(double(:)): Q emissivity\n
      !!         eps22(double(:)): U emissivity\n
      !!         eps23(double(:)): V emissivity
      subroutine emiss2ord(Atom,Geom,vx,vy,vz,lvel,omega,Fed,Red, &
                           RWarr,Norma,Flgsg,jtran,itermu,itermf,iz, &
                           if0,if1,Mif0,Mif1,DwT,Dw,Bfield,vmi, &
                           TBout,Stokes,JKQa,JKQ,JKQC, &
                           eps20,eps21,eps22,eps23)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      type(Redb2_class), intent(inout):: RWarr
      type(Prof_class), intent(in):: Norma
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: lvel
      integer, intent(in):: jtran,itermu,itermf,iz,if0,if1,Mif0,Mif1
      double precision, intent(in):: DwT,Dw,vmi,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(Geom%njdir,if0:if1), &
                        intent(out):: eps20,eps21,eps22,eps23
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(in):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2,Geom%njdir), &
                       intent(in):: TBout

      ! Local

      logical:: PRDc,LPRAM,lNCHLT,cohIn,conj,asym

      integer:: itran,iterml,ith1,iph1,ish1,jdir,idir
      integer:: ifreq,iifreq,iti,iran,ishift,nfs
      integer:: K,iQQ,K1,iPP,iq,iq1,ip,ip1
      integer:: nMl,nMu,nMf,iMl,iMl1,iMu,iMu1,iMf
      integer:: indF,indU,indU1,indL,indL1,indK,icom
      integer:: jjfreq0,jjfreq,kkfreq0,kkfreq0b,kkfreq,kwfreq0
      integer:: llfreq0,nmfreq,iL,iL1,iU,iU1,mF

      double precision:: omegai,rLl,rLu,rLf,S,rJlmax,rJumax,rJfmax
      double precision:: rMl,rMl1,rMu,rMu1,rMf
      double precision:: el,el1,eu,eu1,ef,wlf
      double precision:: QQ,rK1,PP,q,q1,p,p1
      double precision:: al,au,af,auf,aul,Dw1,hau,at,Dfreq,sig
      double precision:: iDw,cost,sint,cosc,sinc,vfac1
      double precision:: Norme0,Norme1,daux
      double precision, dimension(0:3):: StokesM
      double precision, dimension(:,:), allocatable:: Stokesin

      complex(kind=8):: hanleden,prof,y0,intgr,PRD,Norme2
      complex(kind=8), dimension(if0:if1):: CRD,CRD0
      complex(kind=8), dimension(Geom%njdir,if0:if1):: tmp0
      complex(kind=8), dimension(Geom%njdir,if0:if1):: tmp1
      complex(kind=8), dimension(Geom%njdir,if0:if1):: tmp2
      complex(kind=8), dimension(Geom%njdir,if0:if1):: tmp3
      complex(kind=8), dimension(-2:2,0:2):: Jrad
      complex(kind=8), dimension(:), allocatable:: tmpK
      complex(kind=8), dimension(:), allocatable, target:: Warr2
      complex(kind=8), dimension(:), allocatable:: Warr2xW
      complex(kind=8), dimension(:), allocatable:: intergrin
      complex(kind=8), dimension(:,:,:), allocatable, target:: JKQinMV
      complex(kind=8), dimension(:,:,:), allocatable, target:: JradC

      ! Pointers

      type(Redc_class), pointer:: p_red
      type(Redc2_class), pointer:: p_rwarr
      integer, pointer:: p_mfreq
      complex(kind=8), dimension(:), pointer:: p_warr2
      complex(kind=8), dimension(:), pointer:: p_JKQ
      complex(kind=8), dimension(:), pointer:: p_JKQC
      complex(kind=8), dimension(:,:,:,:), pointer:: TKQo


      ! Routine name
      urou = 'emiss2ord'

      ! Initialize pointers
      nullify(p_red,p_rwarr,p_mfreq,p_warr2,p_JKQ,p_JKQC,TKQo)


      !
      ! Construct mean intensity if needed
      !

      ! If angle average and dynamics
      if (lvel.and.AV) then

        ! Point to correct tensors
        TKQo => Geom%TB(:,:,:,:,iz)

        ! Get JKQ in comoving frame
        call getJKQstar(Red,Geom,iz,Atom%ntran,Atom%tif0, &
                        Atom%tif1,DwT,vx,vy,vz,omega,Flgsg, &
                        TKQo,Stokes,JKQa,JKQinMV)

        ! If ad-hoc asymmetries, JradC needs to be rotated
        if (size(JKQa).ge.10) &
          call fieldB_alt(JKQinMV,Red%ggf1-Red%ggf0+1, &
                          Flgsg,Bfield%Btheta(iz), &
                          Bfield%Bphi(iz),1)

        ! Free pointer
        nullify(TKQo)

      ! If angle average but static
      else if (AV) then

        ! Allocate and copy
        allocate(JKQinMV(-2:2,0:2,Red%ggf0:Red%ggf1))
        JKQinMV = JKQC(:,:,Red%ggf0:Red%ggf1)

        ! Rotate only if there is some magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) &
          call fieldB_alt(JKQinMV,Red%ggf1-Red%ggf0+1, &
                          Flgsg,Bfield%Btheta(iz), &
                          Bfield%Bphi(iz),1)

      ! If angle-dependent
      else

        ! Check if ad-hoc asymmetries
        asym = size(JKQa).gt.10

        ! If there are ad-hoc asymmetries
        if (asym) then

          ! Point to correct tensors
          TKQo => Geom%TB(:,:,:,:,iz)

          ! Get JKQ correction for asymmetry
          call getJKQADasym(Red,Geom,iz,Atom%ntran,Atom%tif0, &
                            Atom%tif1,DwT,vx,vy,vz, &
                            .True.,Bfield%Btheta(iz), &
                            Bfield%Bphi(iz),omega,Flgsg, &
                            TKQo,Stokes,JKQC,JKQa,JKQinMV)

        end if ! Ad-hoc asymmetries
      end if ! Dynamics


      !
      ! Get terms and transition quantities
      !

      ! Inverse Doppler width
      iDw = 1d0/Dw

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      at = (au + af + auf)*iDw
      hau = 2d0*(au+auf)*iDw

      ! Units normalization factor for CRD
      Norme0 = (1d5*iDw)*.5d0*sqrt(IPI)

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLf = Atom%rLval(itermf)

      ! Determine the maximum angular momentum and the number
      ! of magnetic sublevels for that maximum momentum
      rJumax = rLu + S
      nMu = nint(2d0*rJumax+1d0)
      rJfmax = rLf + S
      nMf = nint(2d0*rJfmax+1d0)


      !
      ! Initialize the emission coefficient
      !
      eps20 = 0d0
      eps21 = 0d0
      eps22 = 0d0
      eps23 = 0d0

      ! If there are frequencies and angle-dependent
      if (.not.AV.and.Red%mxfreq.gt.0) then

        ! Allocate auxialiar arrays
        allocate(Warr2xW(Red%mxfreq))
        allocate(intergrin(Red%mxfreq))

      end if ! Frequencies to integrate and angle-dependent


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the input transitions
      do iti=1,Atom%trano(jtran)%nt

        ! Get transition index and lower term
        itran = Atom%trano(jtran)%indT(iti)
        iterml = Atom%fst(itran)%iterml

        ! Point to input transition
        p_red => Red%trani(iti)

        ! If PRAM, point to the redistribution subblock
        if (PRAM) then

          ! Point to input transition and get if storing
          p_rwarr => RWarr%trani(iti)
          LPRAM = PRAM.and.p_rwarr%RAM

        ! If not, nothing stored
        else

          ! Not storing locally either
          LPRAM = .False.
          allocate(RWarr%trani(1))
          p_rwarr => Rwarr%trani(1)

        end if ! If storing

        ! If to consider non-coherent lower term and this atom has the
        ! data allocated
        if (NCHLT.and.allocated(Atom%NCHLT)) then

          ! Check if this transition has non-coherent lower term for
          ! this height
          lNCHLT = Atom%NCHLT(iz,itran)

        ! Coherent lower term
        else

          ! Flag non-coherent as false
          lNCHLT = .False.

        end if ! non-coherent lower term

        ! Predict size of interpolation block
        nmfreq = sum(p_red%mfreq)

        ! If angle-averaged
        if (AV) then

          ! Get input radiation field
          call getJKQin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                        JradC,JKQinMV)

        ! If angle-dendent
        else

          ! If Rayleigh scattering and there is coherent
          if (jtran.eq.itran.and.Geom%V_CScatt(1).ge.1d0) then

            ! Scale dimension
            nmfreq = nmfreq*(Geom%nScatt-1)
            nfs = 1

          ! Raman scattering
          else

            ! Scale dimension
            nmfreq = nmfreq*Geom%nScatt
            nfs = 0

          end if ! Rayleigh/Raman scattering

          ! Get interpolated intensity
          call getStkin(Geom,p_red,Fed,Red,Mif0,Mif1,omega, &
                        vx,vy,vz,lvel,Stokesin,Stokes)

          ! If asymmetries, get input radiation field
          if (asym) &
            call getJKQin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                          JradC,JKQinMV)

        end if ! AA or AD

        ! Get the 'flat' JKQ for this input transition
        JRad = JKQ(:,:,itran)

        ! Doppler width for the input transition
        Dw1 = Atom%Dfreq(itran)*sqrt(DwT*DwT + vmi**2d0)

        ! Damping parameter input lower level and input transition
        al = Atom%damp(iterml,iz)
        aul = Atom%ldamp(itran,iz)

        ! Angular momentum input lower level
        rLl = Atom%rLval(iterml)

        ! Determine maximum value of J and number of magnetic
        ! sublevels for this maximum J
        rJlmax = rLl + S
        nMl = nint(2d0*rJlmax+1d0)

        ! Initialize kkfreq index
        kwfreq0 = 0

        ! For each Mf
        do iMf=1,nMf

          ! Value of Mf
          rMf = -rJfmax + dble(iMf-1)

          ! For each mu_f
          do mF=1,Atom%nblk(iMf,itermf)

            ! Get eigenvalue final lower level
            ef = Atom%eval(mF,iMf,itermf,iz)

            ! Get indexes
            indF = Atom%irho(itermf)%jM(mF,iMf)

            ! For each Mu
            do iMu=1,nMu

              ! Value of Mu
              rMu = -rJumax + dble(iMu-1)

              ! Difference between M momentums, done integer
              q = rMu - rMf
              iq = nint(q)

              ! If not pi nor sigma, skip
              if(abs(iq).gt.1) cycle

              ! For each mu_u
              do iU=1,Atom%nblk(iMu,itermu)

                ! Get eigenvalue upper level
                eu = Atom%eval(iU,iMu,itermu,iz)

                ! Get indexes
                indU = Atom%irho(itermu)%jM(iU,iMu)
                indK = Atom%trano(jtran)%indB(indF,indU)

                ! If stored
                if (Norma%VRAM) then

                  ! Get first order profile
                  CRD0 = Norme0*conjg(Norma%cp(:,indK))

                ! If not stored
                else

                  ! Shift term
                  Dfreq = eu - ef + Atom%Dfreq(jtran)

                  ! Normalization
                  Norme1 = Norme0*Norma%Norm(indK)

                  !
                  ! Flat contribution. Implicit branching
                  !

                  ! For each frequency
                  do ifreq=if0,if1

                    ! Calculate profile u-f
                    call voigt((Dfreq - omega(ifreq))*iDw,at,prof)

                    ! Get first order profile
                    CRD0(ifreq) = dcmplx(Norme1*dble(prof), &
                                        -Norme0*dimag(prof))

                  end do ! frequencies

                end if ! Storing Voigt

                ! For each Mu'
                do iMu1=1,nMu

                  ! Value of Mu'
                  rMu1 = -rJumax + dble(iMu1-1)

                  ! Difference between M momentums
                  q1 = rMu1-rMf
                  QQ = q1-q

                  ! Convert to integers
                  iq1 = nint(q1)
                  iQQ = nint(QQ)

                  ! If not pi or sigma, skip
                  if(abs(iq1).gt.1) cycle

                  ! For each mu_u'
                  do iU1=1,Atom%nblk(iMu1,itermu)

                    ! Get eigenvalue upper' level
                    eu1 = Atom%eval(iU1,iMu1,itermu,iz)

                    ! Get indexes
                    indU1 = Atom%irho(itermu)%jM(iU1,iMu1)
                    indK = Atom%trano(jtran)%indB(indF,indU1)

                    ! If stored
                    if (Norma%VRAM) then

                      ! Flat spectrum contribution
                      CRD = CRD0 + Norme0*Norma%cp(:,indK)

                    ! If not stored
                    else

                      ! Shift term
                      Dfreq = eu1 - ef + Atom%Dfreq(jtran)

                      ! Normalization
                      Norme1 = Norme0*Norma%Norm(indK)

                      !
                      ! Flat contribution. Implicit branching
                      !

                      ! For each frequency
                      do ifreq=if0,if1

                        ! Calculate profile u'-f
                        call voigt((Dfreq - omega(ifreq))*iDw, &
                                   at,prof)

                        ! Flag spectrum contribution
                        CRD(ifreq) = CRD0(ifreq) + &
                                     dcmplx(Norme1*dble(prof), &
                                            Norme0*dimag(prof))

                      end do ! frequencies

                    end if ! Storing Voigt

                    ! Hanle factor
                    ! TODO ATTENTION TO THIS
                    hanleden = dcmplx(hau,(eu-eu1)*iDw)

                    !
                    ! Continue with the 2nd order emissivity
                    !

                    ! For each Ml
                    do iMl=1,nMl

                      ! Value of Ml
                      rMl = -rJlmax + dble(iMl-1)

                      ! Difference between M momentums, in integer
                      p = rMu-rMl
                      ip = nint(p)

                      ! If not pi nor sigma, skip
                      if (abs(ip).gt.1) cycle

                      ! For each mu_l
                      do iL=1,Atom%nblk(iMl,iterml)

                        ! Get eigenvalue of input lower level
                        el = Atom%eval(iL,iMl,iterml,iz)

                        ! Get indexes
                        indL = Atom%irho(iterml)%jM(iL,iMl)

                        ! For each Ml'
                        do iMl1=1,nMl

                          ! Value of Ml'
                          rMl1 = -rJlmax + dble(iMl1-1)

                          ! Difference between M momentums
                          p1 = rMu1-rMl1
                          PP = p1-p

                          ! Convert to integer
                          ip1 = nint(p1)
                          iPP = nint(PP)

                          ! If not pi nor sigma, skip
                          if (abs(ip1).gt.1) cycle

                          ! If non-coherent lower term, skip if not
                          ! diagonal
                          if (lNCHLT) then
                            if (iMl.ne.iMl1) cycle
                          end if

                          ! For each mu_l'
                          do iL1=1,Atom%nblk(iMl1,iterml)

                            ! If non-coherent lower term, skip if not
                            ! diagonal
                            if (lNCHLT) then
                              if (iL.ne.iL1) cycle
                            end if

                            ! Get eigenvalue of input lower' level
                            el1 = Atom%eval(iL1,iMl1,iterml,iz)

                            ! Get indexes
                            indL1 = Atom%irho(iterml)%jM(iL1,iMl1)

                            ! Initialize temporal variable
                            tmp0 = cZero
                            tmp1 = cZero
                            tmp2 = cZero
                            tmp3 = cZero

        !
        ! Reset indent
        !

        ! If storing redistribution
        if (LPRAM) then

          ! Get the component index
          icom = Atom%trano(jtran)%trani(iti)% &
                      indB(indL1,indL,indF,indU1,indU)

          ! Get if need to compute redistribution
          PRDc = p_rwarr%iPPRD(icom)

        ! If not storing
        else

          ! Always need to calculate
          PRDc = .True.

        end if ! Storing redistribution function


        !
        ! Create array of Wfunc
        !

        ! If need to compute it and there are frequencies
        if (PRDc.and.nmfreq.gt.0) then

          ! Calculate redistribution function Warr2
          call get_Warr(Atom,Geom,Fed,p_red,p_rwarr,LPRAM,Mif0,Mif1, &
                        itran,jtran,icom,kwfreq0,nmfreq,omega, &
                        Dw,Dw1,el,el1,eu,eu1,ef, &
                        al,au,af,aul,auf,Warr2)

        end if ! Initialized

        ! Allocate tmpK
        allocate(tmpK(abs(iQQ):2))

        ! For each K'
        do K1=abs(iPP),2

          ! Get real value
          rK1 = dble(K1)

          ! Call inner atomic loop
          call emi2_inner_atomic_loop(Atom,Flgsg, &
                                      itermu,iMu,iU,iMu1,iU1, &
                                      itermf,iMf,mF, &
                                      iterml,iMl,iL,iMl1,iL1, &
                                      iz,jtran,itran,K1, &
                                      iQQ,ip,ip1,iq,iq1, &
                                      rMl,rMl1,rK1, &
                                      p,p1,PP,q,q1,QQ, &
                                      tmpK)


          !
          ! Integral over input frequencies
          !

          ! If storing Warr
          if (LPRAM) then

            ! If there are frequencies
            if (nmfreq.gt.0) then

              ! Allocate
              allocate(p_warr2(nmfreq))

              ! Take from stored
              p_warr2 = dcmplx(p_rwarr% &
                                 Pwarr2(kwfreq0+1:kwfreq0+nmfreq))

            end if ! There are frequencies

          ! If not storing
          else

            ! Point to was calculated before
            if (allocated(Warr2)) p_warr2 => Warr2

          end if ! Storing Warr

          ! Difference between l and f energies
          wlf = el - ef

          !
          ! Angle-average (integral)
          !
          if (AV) then

            ! Check coherent shift
            cohIn = abs(wlf).gt.0d0

            !
            ! Check if conjugated
            conj = iPP.lt.0
            sig = Flgsg%sg(iPP)

            ! If coherent
            if (minval(p_red%mfreq).lt.1) then

              ! Just point
              if (conj) then
                p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(-iPP,K1,:)
              else
                p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(iPP,K1,:)
              end if

            end if ! Coherent

            ! Initialize frequency indexes
            jjfreq = 0
            kkfreq = 0

            ! For each output frequency
            iifreq = 0
            do iran=1,Fed%nran
              do ifreq=Fed%if0(iran),Fed%if1(iran)

                ! Advance index
                iifreq = iifreq + 1

                ! Manage MPI
                if (iifreq.lt.Mif0) cycle
                if (iifreq.gt.Mif1) exit

                ! Point to dimension
                p_mfreq => p_red%mfreq(iifreq)

                ! Initialize
                PRD = cZero

                ! If coherent wing
                if (p_mfreq.lt.1) then

                  ! Shifted frequency
                  if (cohIn) then

                    ! Input frequency
                    omegai = omega(ifreq) - wlf

                    ! Shift array beginning
                    ishift = 1 - Red%ggf0

                    ! Get interpolated JKQ
                    y0 = getJKQinnu(omega(Red%ggf0:Red%ggf1), &
                                    p_JKQC, &
                                    ifreq+ishift, &
                                    Atom%tif0(itran)+ishift, &
                                    Atom%tif1(itran)+ishift, &
                                    omegai)

                    ! Fully coherent contribution
                    if (conj) then
                      PRD = sig*conjg(y0)
                    else
                      PRD = y0
                    end if

                  ! Just same frequency
                  else

                    ! Fully coherent contribution
                    if (conj) then
                      PRD = sig*conjg(p_JKQC(ifreq))
                    else
                      PRD = p_JKQC(ifreq)
                    end if

                  end if ! Coherent shift

                ! Non-coherent
                else

                  ! Norm denominator
                  Norme2 = sum(p_warr2(kkfreq+1:kkfreq+p_mfreq)* &
                               p_red%W_freq(jjfreq+1:jjfreq+p_mfreq))

                  ! If valid norm
                  if (abs(Norme2).gt.0d0) then

                    ! If conjugate
                    if (conj) then

                      ! Point to positive Q
                      p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,-iPP,K1)

                      ! Integrate
                      PRD = sig*sum(conjg(p_JKQ)* &
                              p_red%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                              p_warr2(kkfreq+1:kkfreq+p_mfreq))

                    ! Not conjugate
                    else

                      ! Point to positive Q
                      p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,iPP,K1)

                      ! Integrate
                      PRD = sum(p_JKQ* &
                              p_red%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                              p_warr2(kkfreq+1:kkfreq+p_mfreq))

                    end if

                    ! Normalize to the first order profile (the
                    ! product with CRD happens later)
                    PRD = PRD/Norme2

                  end if ! Valid norm

                  ! Update indexes
                  jjfreq = jjfreq + p_mfreq
                  kkfreq = kkfreq + p_mfreq

                end if ! Coherent wing

                ! Substract the flat spectrum part due to just
                ! radiative excitation
                PRD = PRD - Jrad(iPP,K1)

                !
                ! Compute the main part of emiss2ord
                !

                ! Scale with CRD to finish normalization
                PRD = PRD*CRD(ifreq)

                ! For each K
                do K=abs(iQQ),2

                  ! Recycle prof as scaled PRD
                  prof = PRD*tmpK(K)

                  ! Initialize index
                  do jdir=1,Geom%njdir

                    ! Add TKQ to the PRD contribution and accumulate
                    tmp0(jdir,ifreq) = prof*TBout(0,-iQQ,K,jdir) + &
                                       tmp0(jdir,ifreq)
                    tmp1(jdir,ifreq) = prof*TBout(1,-iQQ,K,jdir) + &
                                       tmp1(jdir,ifreq)
                    tmp2(jdir,ifreq) = prof*TBout(2,-iQQ,K,jdir) + &
                                       tmp2(jdir,ifreq)
                    tmp3(jdir,ifreq) = prof*TBout(3,-iQQ,K,jdir) + &
                                       tmp3(jdir,ifreq)

                  end do ! Output directions
                end do ! K
              end do ! Output frequencies
            end do ! Output frequencies ranges

            ! Clean JKQ
            nullify(p_JKQ)
            nullify(p_JKQC)

          !
          ! Angle-dependent (Integral)
          !
          else

            ! Check if dynamic or with shift
            cohIn = lvel.or.abs(wlf).gt.0d0

            ! If asymmetries
            if (asym) then

              !
              ! Check if conjugated
              conj = iPP.lt.0
              sig = Flgsg%sg(iPP)

              ! Just point
              if (conj) then
                p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(-iPP,K1,:)
              else
                p_JKQC(Red%ggf0:Red%ggf1) => JKQinMV(iPP,K1,:)
              end if
            end if ! Asymmetries

            ! For each output direction
            do jdir=1,Geom%njdir

              ! Initialize frequency indexes
              jjfreq = 0 ! Only for asymmetries
              jjfreq0 = 0
              llfreq0 = 0
              kkfreq0 = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fed%nran
                do ifreq=Fed%if0(iran),Fed%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Manage MPI
                  if (iifreq.lt.Mif0) cycle
                  if (iifreq.gt.Mif1) exit

                  ! Initialize
                  PRD = cZero

                  ! Point to dimension
                  p_mfreq => p_red%mfreq(iifreq)

                  ! Input direction index
                  idir = 0

                  ! For each polar direction in quadrature
                  do ith1=1,Geom%nTh

                    ! For each azimuthal direction in quadrature
                    do iph1=1,Geom%nPh2

                      ! Advance index
                      idir = idir + 1

                      ! Scattering index
                      ish1 = Geom%i_scatt(iph1,ith1,jdir)

                      ! Special treatment if forward for two terms
                      if ((jtran.eq.itran.and. &
                           Geom%V_CScatt(ish1).ge.1d0).or. &
                          (p_mfreq.lt.1)) then

                        ! Need to interpolate
                        if (cohIn) then

                          ! Input frequency
                          omegai = omega(ifreq) - wlf

                          ! If axial
                          if (axial) then

                            ! If there are dynamics
                            if (lvel) then

                              ! Get director cosine
                              cost = Geom%V_mu(ith1)

                              ! Calculate Doppler shift factor
                              vfac1 = 1d0 - vz*cost

                              ! We will be using the inverse
                              vfac1 = 1d0/vfac1

                              ! Target frequency
                              omegai = omegai*vfac1

                            end if

                            ! Interpolate
                            StokesM = getStkinnu(omega, &
                                                 Stokes(:,:,1,ith1), &
                                                 ifreq, &
                                                 Atom%tif0(itran), &
                                                 Atom%tif1(itran), &
                                                 omegai)
                          ! If not axial
                          else

                            ! If there are dynamics
                            if (lvel) then

                              ! Get directional trigonimetric f.
                              cost = Geom%V_mu(ith1)
                              sint = sqrt(1d0 - cost*cost)
                              cosc = Geom%v_mux(iph1)
                              sinc = Geom%v_muy(iph1)* &
                                     sqrt(1d0 - cosc*cosc)

                              ! Calculate Doppler shift factor
                              vfac1 = 1d0 - vx*sint*cosc - &
                                            vy*sint*sinc - &
                                            vz*cost

                              ! We will be using the inverse
                              vfac1 = 1d0/vfac1

                              ! Shift
                              omegai = omegai*vfac1

                            end if

                            ! Interpol
                            StokesM = getStkinnu(omega, &
                                            Stokes(:,:,iph1,ith1), &
                                            ifreq, &
                                            Atom%tif0(itran), &
                                            Atom%tif1(itran), &
                                            omegai)

                          end if ! Axiality

                          ! Asymmetry
                          if (asym) then

                            ! Get interpolated JKQ
                            ishift = 1 - Red%ggf0
                            y0 = getJKQinnu(omega(Red%ggf0:Red%ggf1),&
                                            p_JKQC, &
                                            ifreq+ishift, &
                                            Atom%tif0(itran)+ishift, &
                                            Atom%tif1(itran)+ishift, &
                                            omegai)
                            ! Conjugate
                            if (conj) y0 = sig*conjg(y0)

                          end if ! Asymmetry

                        ! Fully coherent
                        else

                          ! Get value at frequency
                          if (axial) then
                            StokesM = Stokes(:,ifreq,1,ith1)
                          else
                            StokesM = Stokes(:,ifreq,iph1,ith1)
                          end if

                          ! Asymmetry
                          if (asym) then

                            ! Get value
                            y0 = p_JKQC(ifreq)

                            ! Conjugate
                            if (conj) y0 = sig*conjg(y0)

                          end if ! Asymmetry

                        end if ! Full coherence

                        ! Sum over Stokes parameters
                        intgr = sum(StokesM* &
                                    Geom%TB(:,iPP,K1,idir,iz))

                        ! Asymmetry
                        if (asym) intgr = intgr + y0

                        ! Add to PRD contribution
                        PRD = PRD + &
                              intgr*Geom%W_mu(ith1)*Geom%W_mux2(iph1)

                      ! Non-forward 2-term scattering
                      else

                        ! Shift in indexes
                        kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                        ! Multiply Warr2 and weights
                        Warr2xW(1:p_mfreq) = &
                              p_warr2(kkfreq0b+1:kkfreq0b+p_mfreq)* &
                              p_red%W_freq(llfreq0+1:llfreq0+p_mfreq)

                        ! Compute norm
                        Norme2 = sum(Warr2xW(1:p_mfreq))

                        ! If valid norm
                        if (abs(Norme2).gt.0d0) then

                          ! Sum Stokes
                          intergrin(1:p_mfreq) = &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,0)* &
                              Geom%TB(0,iPP,K1,idir,iz) + &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,1)* &
                              Geom%TB(1,iPP,K1,idir,iz) + &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,2)* &
                              Geom%TB(2,iPP,K1,idir,iz) + &
                              Stokesin(jjfreq0+1:jjfreq0+p_mfreq,3)* &
                              Geom%TB(3,iPP,K1,idir,iz)

                          ! If asymmetry
                          if (asym) then

                            ! Conjugate
                            if (conj) then

                              ! Point to positive Q
                              p_JKQ => &
                                JradC(jjfreq+1:jjfreq+p_mfreq,-iPP,K1)

                              ! Add
                              intergrin(1:p_mfreq) = &
                                          intergrin(1:p_mfreq) + &
                                          sig*conjg(p_JKQ)

                            ! Not conjugate
                            else

                              ! Point to positive Q
                              p_JKQ => &
                                 JradC(jjfreq+1:jjfreq+p_mfreq,iPP,K1)

                              ! Add
                              intergrin(1:p_mfreq) = &
                                          intergrin(1:p_mfreq) + p_JKQ

                            end if ! Conjugate
                          end if ! asymmetry

                          ! Normalize to the first order profile
                          ! and add the directional weights (the
                          ! product with CRD happens later)
                          PRD = PRD + &
                                sum(Warr2xW(1:p_mfreq)* &
                                    intergrin(1:p_mfreq))* &
                                Geom%W_mu(ith1)* &
                                Geom%W_mux2(iph1)/Norme2

                        end if ! Valid norm
                      end if ! Type of scattering

                      ! Advance index if not axial
                      if (.not.axial) jjfreq0 = jjfreq0 + p_mfreq

                    end do ! azimuthal nodes

                    ! Update jjfreq if axial
                    if (axial) jjfreq0 = jjfreq0 + p_mfreq

                  end do ! polar nodes

                  ! Update jjfreq
                  llfreq0 = llfreq0 + p_mfreq
                  kkfreq0 = kkfreq0 + p_mfreq*(Geom%nScatt-nfs)
                  jjfreq = jjfreq + p_mfreq ! Asymmetries

                  ! Subtract the flat spectrum part due to just
                  ! radiative excitation
                  PRD = PRD - Jrad(iPP,K1)

                  ! Scale with CRD profile, completing the
                  ! normalization
                  PRD = PRD*CRD(ifreq)

                  ! For each K
                  do K=abs(iQQ),2

                    ! Add TKQ to the PRD contribution and accumulate
                    tmp0(jdir,ifreq) = tmpK(K)*TBout(0,-iQQ,K,jdir)* &
                                       PRD + tmp0(jdir,ifreq)
                    tmp1(jdir,ifreq) = tmpK(K)*TBout(1,-iQQ,K,jdir)* &
                                       PRD + tmp1(jdir,ifreq)
                    tmp2(jdir,ifreq) = tmpK(K)*TBout(2,-iQQ,K,jdir)* &
                                       PRD + tmp2(jdir,ifreq)
                    tmp3(jdir,ifreq) = tmpK(K)*TBout(3,-iQQ,K,jdir)* &
                                       PRD + tmp3(jdir,ifreq)

                  end do ! K
                end do ! output frequencies
              end do ! output frequencies ranges
            end do ! Output directions

            ! Clean JKQ
            nullify(p_JKQ)
            nullify(p_JKQC)

          end if ! AA/AD (Integral)

          ! Clean p_warr2
          if (LPRAM.and.nmfreq.gt.0) deallocate(p_warr2)
          nullify(p_warr2)

        end do ! K1

        ! Free tmpK
        deallocate(tmpK)

        ! Update redistribution initial index
        kwfreq0 = kwfreq0 + nmfreq

        ! Apply hanle factor and Einstein coefficient
        daux = (2d0*rLl+1d0)*Atom%Ecoeff(iterml,itermu)
        eps20 = eps20 + dble(tmp0/hanleden)*daux
        eps21 = eps21 + dble(tmp1/hanleden)*daux
        eps22 = eps22 + dble(tmp2/hanleden)*daux
        eps23 = eps23 + dble(tmp3/hanleden)*daux

                            !
                            ! Recover indentation
                            !

                          end do ! iL1
                        end do ! Ml1
                      end do ! iL
                    end do ! Ml
                  end do ! iU1
                end do ! Mu1
              end do ! iU
            end do ! Mu
          end do ! mF
        end do ! Mf

        ! Not storing
        if (.not.PRAM) deallocate(Rwarr%trani)

      end do ! Input transitions

      ! Apply common factor
      daux = 3d0*.5d0*IPI42*(2d0*rLu+1d0)* &
             Atom%Ecoeff(itermu,itermf)*1d-10/(c*Dw)
      eps20 = eps20*daux
      eps21 = eps21*daux
      eps22 = eps22*daux
      eps23 = eps23*daux

      ! Clean pointers
      if (associated(p_red)) nullify(p_red)
      if (associated(p_rwarr)) nullify(p_rwarr)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_warr2)) nullify(p_warr2)
      if (associated(p_JKQ)) nullify(p_JKQ)
      if (associated(p_JKQC)) nullify(p_JKQC)

      return

      end subroutine emiss2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the two most internal loops in the calculation of
      !! the absorption coefficients in the presence of a magnetic
      !! field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                      and J-symbols\n
      !!         iz(integer): Height index\n
      !!      iterm(integer): Term index\n
      !!        iJ1(integer): J index of the level\n
      !!        iJ2(integer): J' index of the level\n
      !!        rJ1(integer): J angular momentum\n
      !!        rJ2(integer): J' angular momentum\n
      !!        rM1(integer): M magnetic quantum number\n
      !!        rM2(integer): M' magnetic quantum number\n
      !!        Kin(integer): K in the call from emiss2ord\n
      !!       summ(dcomplx): Summation result
      subroutine absorb1(Atom,Flgsg,iz,iterm,iJ1,iJ2,rJ1,rJ2,rM1, &
                         rM2,Kin,summ)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,iterm,iJ1,iJ2,Kin
      double precision, intent(in):: rJ1,rJ2,rM1,rM2
      complex(kind=8), intent(out):: summ

      ! Local

      integer:: K,iQ,Kmin,Kmax,iR

      double precision:: rK,Q,fKf3j


      ! Difference between magnetic momenta
      Q = rM1 - rM2

      ! Convert to integer
      iQ = nint(Q)

      ! Determine the limits in K
      Kmin = max(abs(iQ),nint(abs(rJ1-rJ2)))
      Kmax = min(nint(rJ1+rJ2),Atom%Kcut(iterm))

      ! Initialize the output
      summ = cZero

      ! For each K
      do K=Kmin,Kmax
     !do K=Kmin,0!Kmax

        ! Check for Kcut
        if (KcutAB.and.abs(K-Kin).gt.Atom%Kcut(iterm)) cycle

        ! Get the real number
        rK = dble(K)

        ! Get the SEE index
        iR = Atom%irho(iterm)%Jrho(iJ2,iJ1)%kq(iQ,K)

        ! If flagged as small, skip
        if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

        ! Racah algebra
        fKf3j = sqrt(2d0*rK+1d0)*fun3j(rJ1,rJ2,rK,rM1,-rM2,-Q,Flgsg)

        ! Accumulate in the sum
        summ = fKf3j*Atom%crho(iR,iz) + summ

      end do ! K

      end subroutine absorb1

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the two most internal loops in the calculation of
      !! the emission coefficients in the presence of a magnetic
      !! field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                      and J-symbols\n
      !!         iz(integer): Height index\n
      !!      iterm(integer): Term index\n
      !!        iJ1(integer): J index of the level\n
      !!        iJ2(integer): J' index of the level\n
      !!        rJ1(integer): J angular momentum\n
      !!        rJ2(integer): J' angular momentum\n
      !!        rM1(integer): M magnetic quantum number\n
      !!        rM2(integer): M' magnetic quantum number\n
      !!       summ(dcomplx): Summation result
      subroutine emiss1(Atom,Flgsg,iz,iterm,iJ1,iJ2,rJ1,rJ2, &
                        rM1,rM2,summ)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,iterm,iJ1,iJ2
      double precision, intent(in):: rJ1,rJ2,rM1,rM2
      complex(kind=8), intent(out):: summ

      ! Local

      integer:: K,iQ,Kmin,Kmax,iR

      double precision:: rK,Q,fKf3j


      ! Difference between magnetic momenta
      Q = rM1-rM2

      ! Convert to integer
      iQ = nint(Q)

      ! Determine the limits in K
      Kmin = max(abs(iQ),nint(abs(rJ1-rJ2)))
      Kmax = min(nint(rJ1+rJ2),Atom%Kcut(iterm))

      ! Initialize the output
      summ = cZero

      ! For each K
      do K=Kmin,Kmax
     !do K=Kmin,0

        ! Get the real number
        rK = dble(K)

        ! Get the SEE index
        iR = Atom%irho(iterm)%Jrho(iJ2,iJ1)%kq(iQ,K)

        ! If flagged as small, skip
        if (iR.le.0.or.Atom%rhonull(iR,iz)) cycle

        ! Racah algebra
        fKf3j = sqrt(2d0*rK+1d0)*fun3j(rJ1,rJ2,rK,rM1,-rM2,-Q,Flgsg)

        !!!!
        ! Uncomment the following statement for natural population
        !!!!
       !if (K.eq.0) summ = dcmplx(fKf3j*sqrt(2d0*rJ1+1d0),.0d0)

        ! Accumulate in the sum
        summ = fKf3j*Atom%crho(iR,iz) + summ

      end do ! K

      end subroutine emiss1

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the redistribution function\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!       Fed(Reda_class): Structure with redistribution output
      !!                        frequency data\n
      !!     p_red(Redc_class): Structure with redistribution input
      !!                        frequency data for a given input
      !!                        transition\n
      !!  p_rwarr(Redc2_class): Structure to store redistribution
      !!                        functions\n
      !!        LPRAM(logical): If storing the redistribution
      !!                        function in RAM\n
      !!         Mif0(integer): First frequency for this CPU\n
      !!         Mif1(integer): Last frequency for this CPU\n
      !!        itran(integer): Input transition index\n
      !!        jtran(integer): Output transition index\n
      !!         icom(integer): Index of the magnetic component of the
      !!                        combination of input and output
      !!                        transitions\n
      !!      kwfreq0(integer): Offset for the redistribution
      !!                        function of the icom component\n
      !!       nmfreq(integer): Total amount of input frequencies
      !!                        for this combination of output and
      !!                        input transitions in this CPU\n
      !!         omega(double): Frequency array\n
      !!            Dw(double): Doppler width of the output
      !!                        transition\n
      !!           Dw1(double): Doppler width of the input
      !!                        transition\n
      !!            el(double): Energy of the component l\n
      !!           el1(double): Energy of the component l'\n
      !!            eu(double): Energy of the component u\n
      !!           eu1(double): Energy of the component u'\n
      !!            ef(double): Energy of the component f\n
      !!            al(double): Inverse lifetime of the lower term\n
      !!            au(double): Inverse lifetime of the upper term\n
      !!            af(double): Inverse lifetime of the final term\n
      !!           aul(double): Damping parameter of the input
      !!                        transition\n
      !!           auf(double): Damping parameter of the output
      !!                        transition\n
      !!    Warr2(dcomplex(:)): Redistribution function
      subroutine get_Warr(Atom,Geom,Fed,p_red,p_rwarr,LPRAM, &
                          Mif0,Mif1,itran,jtran,icom,kwfreq0, &
                          nmfreq,omega,Dw,Dw1,el,el1,eu,eu1,ef, &
                          al,au,af,aul,auf,Warr2)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Reda_class), intent(in):: Fed
      type(Redc_class), intent(in), pointer:: p_red
      type(Redc2_class), intent(inout):: p_rwarr
      logical, intent(in):: LPRAM
      integer, intent(in):: icom,nmfreq,jtran,itran,kwfreq0,Mif0,Mif1
      double precision, intent(in):: el,el1,eu,eu1,ef,al,au,af,aul,auf
      double precision, intent(in):: Dw,Dw1
      double precision, dimension(:), intent(in):: omega
      complex(kind=8), dimension(:), &
                       allocatable, intent(inout):: Warr2

      ! Local

      integer:: kkfreq0,kkfreq,jjfreq0,jjfreq,iifreq,iran,ifreq,jfreq
      integer:: ith1,ish1,stype

      double precision:: omegao,omegai,rep,imp

      ! Pointers

      integer, pointer:: p_mfreq


      ! If Warr2 not allocated
      if (.not.allocated(Warr2)) then

        ! Allocate
        allocate(Warr2(nmfreq))

      ! If already allocated
      else

        ! If wrong size
        if (size(Warr2).ne.nmfreq) then

          ! Re-allocate
          deallocate(Warr2)
          allocate(Warr2(nmfreq))

        end if ! Wrong-size
      end if ! Allocated or not


      !
      ! Angle-average (Warr2)
      !
      if (AV) then

        ! Initialize Warr2
        Warr2 = cZero

        ! Initialize frequency index
        kkfreq = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Fed%nran
          do ifreq=Fed%if0(iran),Fed%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Manage MPI
            if (iifreq.lt.Mif0) cycle
            if (iifreq.gt.Mif1) exit

            ! Point to dimension
            p_mfreq => p_red%mfreq(iifreq)

            ! Skip coherent
            if (p_mfreq.lt.1) cycle

            ! Get output frequency
            omegao = omega(ifreq) - Atom%Dfreq(jtran)

            ! For each input frequency
            do jfreq=1,p_mfreq

              ! Advance indexes
              kkfreq = kkfreq + 1

              ! Get input frequency
              omegai = p_red%omega(kkfreq) - Atom%Dfreq(itran)

              ! For each direction in the integral AA quadrature
              do ith1=1,Geom%nThAA

                ! Add the contribution to the angular integral
                ! of the redistribution function
                Warr2(kkfreq) = Warr2(kkfreq) + &
                                Geom%W_muAA(ith1)* &
                                Wfunc(omegai,omegao, &
                                      Dw,Dw1,el,el1,eu,eu1,ef, &
                                             al,au,af,aul,auf, &
                                      Geom%V_muAA(ith1), &
                                      Geom%V_siAA(ith1),0)*IPI42

              end do  ! direction nodes
            end do ! input frequencies
          end do ! output frequencies
        end do ! output frequencies ranges


      !
      ! Angle-dependent (Warr2)
      !
      else

        ! Initialize frequency index
        jjfreq0 = 0
        kkfreq0 = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Fed%nran
          do ifreq=Fed%if0(iran),Fed%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Manage MPI
            if (iifreq.lt.Mif0) cycle
            if (iifreq.gt.Mif1) exit

            ! Point to dimension
            p_mfreq => p_red%mfreq(iifreq)

            ! Coherent wing
            if (p_mfreq.lt.1) cycle

            ! Get output frequency
            omegao = omega(ifreq) - Atom%Dfreq(jtran)

            ! For each scattering angle
            do ish1=1,Geom%nScatt

              ! Check forward scattering two-terms, skip
              if (itran.eq.jtran.and. &
                  Geom%V_CScatt(ish1).ge.1d0) cycle

              ! Check forward or backward
              if (Geom%V_SScatt(ish1).le.0d0) then
                stype = 1
              else
                stype = 0
              end if

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq
                kkfreq = kkfreq0 + jfreq

                ! Get input frequency
                omegai = p_red%omega(jjfreq) - &
                         Atom%Dfreq(itran)

                ! Calculate redistribution function
                Warr2(kkfreq) = Wfunc(omegai,omegao, &
                                      Dw,Dw1,el,el1,eu,eu1,ef, &
                                      al,au,af,aul,auf, &
                                      Geom%V_CScatt(ish1), &
                                      Geom%V_SScatt(ish1), &
                                      stype)*IPI42

              end do ! input frequencies

              ! Update kkfreq0
              kkfreq0 = kkfreq0 + p_mfreq

            end do  ! Scattering angles

            ! Update jjfreq0
            jjfreq0 = jjfreq0 + p_mfreq

          end do ! output frequencies
        end do ! output frequencies ranges

      end if ! AA/AD (Warr2)

      ! If storing
      if (LPRAM) then

        ! Signal that this one does not need to be calculated
        ! anymore
        p_rwarr%iPPRD(icom) = .False.

        ! For the just calculated frequencies
        do jfreq=kkfreq-nmfreq+1,kkfreq

          ! Split real and imaginary parts
          rep = dble(Warr2(jfreq))
          imp = dimag(Warr2(jfreq))

          ! Zero our small for single precision to avoid
          ! underflow
          if (rep.le.TINYWAR) rep = 0d0
          if (abs(imp).le.TINYWAR) imp = 0d0

          ! Store in single precision
          p_rwarr%PWarr2(kwfreq0+jfreq) = &
                                      cmplx(real(rep),real(imp))

        end do ! Calculated frequencies

      end if ! If storing

      end subroutine get_Warr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the inner atomic loops in the calculation of the
      !! absorptivity in the presence of a magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!     itermu(integer): Upper term of the transition\n
      !!        iMu(integer): Index of the upper level magnetic
      !!                      component\n
      !!         iU(integer): Index of the upper level component\n
      !!     iterml(integer): Lower term of the transition\n
      !!        iMl(integer): Index of the lower level magnetic
      !!                      component\n
      !!         iL(integer): Index of the lower level component\n
      !!       iMl1(integer): Index of the other lower level magnetic
      !!                      component\n
      !!         iz(integer): Height index\n
      !!      itran(integer): Index of transition to compute\n
      !!         iq(integer): q value\n
      !!        iq1(integer): q' value\n
      !!         rMl(double): Magnetic quantum number of the lower
      !!                      level magnetic component\n
      !!        rMl1(double): Magnetic quantum number of the other
      !!                      lower level magnetic component\n
      !!        ftmp(double): Multiplicative factor from the
      !!                      outer loops\n
      !!         tk(dcomplx): Result of the sum
      subroutine abs_inner_atomic_loop(Atom,Flgsg,itermu,iMu,iU, &
                                       iterml,iMl,iL,iMl1,iz,itran, &
                                       iq,iq1,rMl,rMl1,ftmp,tk)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: itermu,iMu,iU
      integer, intent(in):: iterml,iMl,iL,iMl1
      integer, intent(in):: iz,itran,iq,iq1
      double precision, intent(in):: rMl,rMl1,ftmp
      complex(kind=8), intent(out):: tK

      ! Local

      integer:: kU,kL,kU1,kL1,kLb,iJu,iJl,iJu1,iJl1,iJlb

      double precision:: Cu,Cl,CC,Cu1,CC1,Clb,rJu,rJl,rJu1,rJl1,rJlb

      complex(kind=8):: rhoc


      ! Initialize tK
      tK = cZero

      ! For each Ju
      do kU=1,Atom%nblk(iMu,itermu)

        ! Get eigenvector upper level
        Cu = Atom%evec(kU,iU,iMu,itermu,iz)

        ! If coefficient too small, skip
        if (abs(Cu).lt.TINYEV) cycle

        ! Get J level index
        iJu = Atom%iJval(kU,iMu,itermu)

        ! Get angular momentum
        rJu = Atom%rJval(iJu,itermu)

        ! For each Jl
        do kL=1,Atom%nblk(iMl,iterml)

          ! Get eigenvector lower level
          Cl = Atom%evec(kL,iL,iMl,iterml,iz)

          ! If coefficient too small, skip
          if (abs(Cl).lt.TINYEV) cycle

          ! Get J level index
          iJl = Atom%iJval(kL,iMl,iterml)

          ! Get angular momentum
          rJl = Atom%rJval(iJl,iterml)

          ! Coefficients times the dipolar matrix
          CC = Cl*Cu*Atom%rdip(itran)%rdip(iq,iMu,iMl,iJu,iJl)

          ! If coefficient small, skip
          if (abs(CC).lt.TINYEV) cycle

          ! For each Ju'
          do kU1=1,Atom%nblk(iMu,itermu)

            ! Get eigenvector upper level'
            Cu1 = Atom%evec(kU1,iU,iMu,itermu,iz)

            ! If coefficient too small, skip
            if (abs(Cu1).lt.TINYEV) cycle

            ! Get J level index
            iJu1 = Atom%iJval(kU1,iMu,itermu)

            ! Get angular momentum
            rJu1 = Atom%rJval(iJu1,itermu)

            ! For each Jl'
            do kL1=1,Atom%nblk(iMl1,iterml)

              ! Get J level index
              iJl1 = Atom%iJval(kL1,iMl1,iterml)

              ! Get angular momentum
              rJl1 = Atom%rJval(iJl1,iterml)

              ! Coefficient times dipolar matrix
              CC1 = Cu1*Atom%rdip(itran)%rdip(iq1,iMu,iMl1,iJu1,iJl1)

              ! If coefficient small, skip
              if (abs(CC1).lt.TINYCO) cycle

              ! For each Jlb
              do kLb=1,Atom%nblk(iMl,iterml)

                ! Get eigenvector for lower level b
                Clb = Atom%evec(kLb,iL,iMl,iterml,iz)

                ! If coefficient too small, skip
                if (abs(Clb).lt.TINYEV) cycle

                ! Get J level index
                iJlb = Atom%iJval(kLb,iMl,iterml)

                ! Get angular momentum
                rJlb = Atom%rJval(iJlb,iterml)

                ! Coefficient and sign
                Clb = Clb*Flgsg%sg(nint(rJlb-rMl))

                ! Sum over (Kl,Ql)
                call absorb1(Atom,Flgsg,iz,iterml,iJlb, &
                             iJl1,rJlb,rJl1,rMl,rMl1,0,rhoc)

                ! If no population, skip
                if (abs(rhoc).lt.TINYER) cycle

                ! Accumulate into tK
                tK = ftmp*CC*CC1*Clb*rhoc + tK

              end do ! kLb
            end do ! kL1
          end do ! kU1
        end do ! kL
      end do ! kU

      end subroutine abs_inner_atomic_loop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the inner atomic loops in the calculation of the
      !! emissivity in the presence of a magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!     itermu(integer): Upper term of the transition\n
      !!        iMu(integer): Index of the upper level magnetic
      !!                      component\n
      !!         iU(integer): Index of the upper level component\n
      !!       iMu1(integer): Index of the other upper level magnetic
      !!                      component\n
      !!     iterml(integer): Lower term of the transition\n
      !!        iMl(integer): Index of the lower level magnetic
      !!                      component\n
      !!         iL(integer): Index of the lower level component\n
      !!         iz(integer): Height index\n
      !!      itran(integer): Index of transition to compute\n
      !!         iq(integer): q value\n
      !!        iq1(integer): q' value\n
      !!         rMu(double): Magnetic quantum number of the upper
      !!                      level magnetic component\n
      !!        rMu1(double): Magnetic quantum number of the other
      !!                      upper level magnetic component\n
      !!        ftmp(double): Multiplicative factor from the
      !!                      outer loops\n
      !!         tk(dcomplx): Result of the sum
      subroutine emi_inner_atomic_loop(Atom,Flgsg,itermu,iMu,iU, &
                                       iMu1,iterml,iMl,iL,iz,itran, &
                                       iq,iq1,rMu,rMu1,ftmp,tk)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: itermu,iMu,iU,iMu1
      integer, intent(in):: iterml,iMl,iL
      integer, intent(in):: iz,itran,iq,iq1
      double precision, intent(in):: rMu,rMu1,ftmp
      complex(kind=8), intent(out):: tK

      ! Local

      integer:: kL,kU,kL1,kU1,kUb,iJl,iJu,iJl1,iJu1,iJub

      double precision:: Cl,Cu,CC,Cl1,CC1,Cub,rJl,rJu,rJl1,rJu1,rJub

      complex(kind=8):: rhoc


      ! Initialize tK
      tK = cZero

      ! For each Jl
      do kL=1,Atom%nblk(iMl,iterml)

        ! Get eigenvector lower level
        Cl = Atom%evec(kL,iL,iMl,iterml,iz)

        ! If coefficient too small, skip
        if (abs(Cl).lt.TINYEV) cycle

        ! Get J level index
        iJl = Atom%iJval(kL,iMl,iterml)

        ! Get angular momentum
        rJl = Atom%rJval(iJl,iterml)

        ! For each Ju
        do kU=1,Atom%nblk(iMu,itermu)

          ! Get eigenvector lower level
          Cu = Atom%evec(kU,iU,iMu,itermu,iz)

          ! If coefficient too small, skip
          if (abs(Cu).lt.TINYEV) cycle

          ! Get J level index
          iJu = Atom%iJval(kU,iMu,itermu)

          ! Get angular momentum
          rJu = Atom%rJval(iJu,itermu)

          ! Coefficients times the dipolar matrix
          CC = Cu*Cl*Atom%rdip(itran)%rdip(iq,iMu,iMl,iJu,iJl)

          ! If coefficient small, skip
          if (abs(CC).lt.TINYCO) cycle

          ! For each Jl'
          do kL1=1,Atom%nblk(iMl,iterml)

            ! Get eigenvector upper level'
            Cl1 = Atom%evec(kL1,iL,iMl,iterml,iz)

            ! If coefficient too small, skip
            if (abs(Cl1).lt.TINYEV) cycle

            ! Get J level index
            iJl1 = Atom%iJval(kL1,iMl,iterml)

            ! Get angular momentum
            rJl1 = Atom%rJval(iJl1,iterml)

            ! For each Ju'
            do kU1=1,Atom%nblk(iMu1,itermu)

              ! Get J level index
              iJu1 = Atom%iJval(kU1,iMu1,itermu)

              ! Get angular momentum
              rJu1 = Atom%rJval(iJu1,itermu)

              ! Coefficient times dipolar matrix
              CC1 = Cl1*Flgsg%sg(nint(rJu1-rMu1))* &
                    Atom%rdip(itran)%rdip(iq1,iMu1,iMl,iJu1,iJl1)

              ! If coefficient small, skip
              if (abs(CC1).lt.TINYCO) cycle

              ! For each Jub
              do kUb=1,Atom%nblk(iMu,itermu) ! sum Jub

                ! Get eigenvector for upper level b
                Cub = Atom%evec(kUb,iU,iMu,itermu,iz)

                ! If coefficient too small, skip
                if (abs(Cub).lt.TINYEV) cycle

                ! Get J level index
                iJub = Atom%iJval(kUb,iMu,itermu)

                ! Get angular momentum
                rJub = Atom%rJval(iJub,itermu)

                ! Sum over (Kl,Ql)
                call emiss1(Atom,Flgsg,iz,itermu,iJu1, &
                            iJub,rJu1,rJub,rMu1,rMu,rhoc)

                ! Uncomment the following line for Zeeman
               !if (iJub.ne.iJu1) rhoc=cZero

                ! If no population, skip
                if (abs(rhoc).lt.TINYER) cycle

                ! Accumulate into tK
                tK = ftmp*CC*CC1*Cub*rhoc + tK

              end do ! kUb
            end do ! kU1
          end do ! kL1
        end do ! kU
      end do ! kL

      end subroutine emi_inner_atomic_loop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the inner atomic loops in the calculation of the
      !! second order emissivity in the presence of a magnetic field\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!     itermu(integer): Upper term of the transition\n
      !!        iMu(integer): Index of the upper level magnetic
      !!                      component\n
      !!         iU(integer): Index of the upper level component\n
      !!       iMu1(integer): Index of the other upper level magnetic
      !!                      component\n
      !!        iU1(integer): Index of the other upper level
      !!                      component\n
      !!     itermf(integer): Final term of the transition\n
      !!        iMf(integer): Index of the final level magnetic
      !!                      component\n
      !!         mF(integer): Index of the final level component\n
      !!     iterml(integer): Lower term of the transition\n
      !!        iMl(integer): Index of the lower level magnetic
      !!                      component\n
      !!         iL(integer): Index of the lower level component\n
      !!       iMl1(integer): Index of the other lower level magnetic
      !!                      component\n
      !!        iL1(integer): Index of the other lower level
      !!                      component\n
      !!         iz(integer): Height index\n
      !!      jtran(integer): Index of the output transition\n
      !!      itran(integer): Index of the input transition\n
      !!         K1(integer): K' multipole integer value\n
      !!        iQQ(integer): Q integer value\n
      !!         ip(integer): p integer value\n
      !!        ip1(integer): p' integer value\n
      !!         iq(integer): q integer value\n
      !!        iq1(integer): q' integer value\n
      !!         rMl(double): Magnetic quantum number of the lower
      !!                      level magnetic component\n
      !!        rMl1(double): Magnetic quantum number of the other
      !!                      lower level magnetic component\n
      !!         rK1(double): K' multipole value\n
      !!           p(double): p value\n
      !!          p1(double): p' value\n
      !!          PP(double): P value\n
      !!           q(double): q value\n
      !!          q1(double): q' value\n
      !!          QQ(double): Q value\n
      !!    tmpk(dcomplx(:)): Result of the sum for each K
      subroutine emi2_inner_atomic_loop(Atom,Flgsg, &
                                        itermu,iMu,iU,iMu1,iU1, &
                                        itermf,iMf,mF, &
                                        iterml,iMl,iL,iMl1,iL1, &
                                        iz,jtran,itran,K1, &
                                        iQQ,ip,ip1,iq,iq1, &
                                        rMl,rMl1,rK1, &
                                        p,p1,PP,q,q1,QQ, &
                                        tmpK)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: itermu,iMu,iU,iMu1,iU1
      integer, intent(in):: itermf,iMf,mF
      integer, intent(in):: iterml,iMl,iL,iMl1,iL1
      integer, intent(in):: iz,jtran,itran,K1,ip,ip1,iq,iq1,iQQ
      double precision, intent(in):: rMl,rMl1
      double precision, intent(in):: rK1,p,p1,PP,q,q1,QQ
      complex(kind=8), dimension(abs(iQQ):2), intent(out):: tmpK

      ! Local

      integer:: kLb,iJlb,kLb1,iJlb1,kL,iJl,kL1,iJl1
      integer:: kU2,iJu2,kU3,iJu3,kU,iJu,kU1,iJu1
      integer:: kF,iJf,kF1,iJf1,K

      double precision:: rK,ftmp,f1tmp,rJlb,rJlb1,rJl,rJl1
      double precision:: rJu2,rJu3,rJu,rJu1,rJf,rJf1
      double precision:: Clb,Clb1,Cl,Cl1,Cu2,Cu3,CC2,CC3,Cu,Cu1
      double precision:: Cf,Cf1,CC,CC1

      complex(kind=8):: rhoc


      ! Initialize
      tmpK = cZero

      ! For each K
      do K=abs(iQQ),2

        ! Get real value
        rK = dble(K)

        ! Racah algebra
        ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

        ! If forbidden (3j-sym=0), skip
        if (abs(ftmp).lt.TINYJS) cycle

        ! Racah algebra
        f1tmp = fun3j(1d0,1d0,rK1,-p,p1,-PP,Flgsg)

        ! If forbidden (3j-sym=0), skip
        if (abs(f1tmp).lt.TINYJS) return

        ! Add the rest of the factor
        ftmp = ftmp*Flgsg%sg(iq)*sqrt(2d0*rK+1d0)
        f1tmp = f1tmp*Flgsg%sg(ip1)*sqrt(2d0*rK1+1d0)

        ! For each Jlb
        do kLb=1,Atom%nblk(iMl,iterml)

          ! Get eigenvector for lower b level
          Clb = Atom%evec(kLb,iL,iMl,iterml,iz)

          ! If coefficient too small, skip
          if (abs(Clb).lt.TINYEV) cycle

          ! Get J level index
          iJlb = Atom%iJval(kLb,iMl,iterml)

          ! Get angular momentum
          rJlb = Atom%rJval(iJlb,iterml)

          ! For each Jlb'
          do kLb1=1,Atom%nblk(iMl1,iterml)

            ! Get eigenvector for lower b level
            Clb1 = Atom%evec(kLb1,iL1,iMl1,iterml,iz)

            ! If coefficient too small, skip
            if (abs(Clb1).lt.TINYEV) cycle

            ! Get J level index
            iJlb1 = Atom%iJval(kLb1,iMl1,iterml)

            ! Get angular momentum
            rJlb1 = Atom%rJval(iJlb1,iterml)

            ! sum over (Kl,Ql)
            call absorb1(Atom,Flgsg,iz,iterml,iJlb,iJlb1, &
                         rJlb,rJlb1,rMl,rMl1,K1,rhoc)

            !!!!
            !  Uncomment the following line for Zeeman effect
           !if (iJlb.ne.iJlb1.of.iFlb.ne.iFlb1) rhoc=cZero
            !!!!

            !!!!
            ! Uncomment the following line for incoherent lower
            !  term
           !if (.not.(rMl1.eq.rMl.and.iL1.eq.iL)) rhoc=cZero
            !!!!

            ! If no rhoKQ, skip
            if (abs(rhoc).lt.TINYER) cycle

            ! Add coefficients to rhoKQ
            rhoc = Flgsg%sg(nint(rJlb-rMl))*Clb*Clb1*rhoc

            ! For each Jl
            do kL=1,Atom%nblk(iMl,iterml)

              ! Get eigenvector for lower level
              Cl = Atom%evec(kL,iL,iMl,iterml,iz)

              ! If coefficient too small, skip
              if (abs(Cl).lt.TINYEV) cycle

              ! Get J level index
              iJl = Atom%iJval(kL,iMl,iterml)

              ! Get angular momentum
              rJl = Atom%rJval(iJl,iterml)

              ! For each Ju''
              do kU2=1,Atom%nblk(iMu,itermu)

                ! Get eigenvector for upper'' level
                Cu2 = Atom%evec(kU2,iU,iMu,itermu,iz)

                ! If coefficient too small, skip
                if (abs(Cu2).lt.TINYEV) cycle

                ! Get J level index
                iJu2 = Atom%iJval(kU2,iMu,itermu)

                ! Get angular momentum
                rJu2 = Atom%rJval(iJu2,itermu)

                ! Add coefficient to dipole strength
                CC2 = Cl*Cu2*Atom%rdip(itran)% &
                                  rdip(ip,iMu,iMl,iJu2,iJl)

                ! If coefficient too small, skip
                if (abs(CC2).lt.TINYCO) cycle

                ! For each Jl'
                do kL1=1,Atom%nblk(iMl1,iterml)

                  ! Get eigenvector for lower' level
                  Cl1 = Atom%evec(kL1,iL1,iMl1,iterml,iz)

                  ! If coefficient too small, skip
                  if (abs(Cl1).lt.TINYEV) cycle

                  ! Get J level index
                  iJl1 = Atom%iJval(kL1,iMl1,iterml)

                  ! Get angular momentum
                  rJl1 = Atom%rJval(iJl1,iterml)

                  ! For each Ju'''
                  do kU3=1,Atom%nblk(iMu1,itermu)

                    ! Get eigenvector for upper''' level
                    Cu3 = Atom%evec(kU3,iU1,iMu1,itermu,iz)

                    ! If coefficient too small, skip
                    if (abs(Cu3).lt.TINYEV) cycle

                    ! Get J level index
                    iJu3 = Atom%iJval(kU3,iMu1,itermu)

                    ! Get angular momentum
                    rJu3 = Atom%rJval(iJu3,itermu)

                    ! KCUT
!                   if (abs(rJu3-rJu2).gt.rK1) cycle

                    ! Add coefficients to dipole strength
                    CC3 = Cl1*Cu3* &
                          Atom%rdip(itran)% &
                               rdip(ip1,iMu1,iMl1,iJu3,iJl1)

                    ! If coefficient too small, skip
                    if (abs(CC3).lt.TINYCO) cycle

        !
        ! Reset identation
        !

        ! For each Jf
        do kF=1,Atom%nblk(iMf,itermf)

          ! Get eigenvector for final level
          Cf = Atom%evec(kF,mF,iMf,itermf,iz)

          ! If coefficient too small, skip
          if (abs(Cf).lt.TINYEV) cycle

          ! Get J level index
          iJf = Atom%iJval(kF,iMf,itermf)

          ! Get angular momentum
          rJf = Atom%rJval(iJf,itermf)

          ! For each Ju
          do kU=1,Atom%nblk(iMu,itermu)

            ! Get eigenvector for upper level
            Cu = Atom%evec(kU,iU,iMu,itermu,iz)

            ! If coefficient too small, skip
            if (abs(Cu).lt.TINYEV) cycle

            ! Get J level index
            iJu = Atom%iJval(kU,iMu,itermu)

            ! Get angular momentum
            rJu = Atom%rJval(iJu,itermu)

            ! Add coefficients to dipole strength
            CC = Cf*Cu*Atom%rdip(jtran)% &
                            rdip(iq,iMu,iMf,iJu,iJf)

            ! If coefficient too small, skip
            if (abs(CC).lt.TINYCO) cycle

            ! For each Jf'
            do kF1=1,Atom%nblk(iMf,itermf)

              ! Get eigenvector for final' level
              Cf1 = Atom%evec(kF1,mF,iMf,itermf,iz)

              ! If coefficient too small, skip
              if (abs(Cf1).lt.TINYEV) cycle

              ! Get J level index
              iJf1 = Atom%iJval(kF1,iMf,itermf)

              ! Get angular momentum
              rJf1 = Atom%rJval(iJf1,itermf)

              ! For each Ju'
              do kU1=1,Atom%nblk(iMu1,itermu)

                ! Get eigenvector for upper' level
                Cu1 = Atom%evec(kU1,iU1,iMu1,itermu,iz)

                ! If coefficient too small, skip
                if (abs(Cu1).lt.TINYEV) cycle

                ! Get J level index
                iJu1 = Atom%iJval(kU1,iMu1,itermu)

                ! Get angular momentum
                rJu1 = Atom%rJval(iJu1,itermu)

                ! Add coefficients to dipole strength
                CC1 = Cf1*Cu1*Atom%rdip(jtran)% &
                                   rdip(iq1,iMu1,iMf,iJu1,iJf1)

                ! If coefficient big enough, add contribution to
                ! temporal variable
                if (abs(CC1).gt.TINYCO) &
                  tmpK(K) = f1tmp*ftmp*CC*CC1*CC2*CC3*rhoc + tmpK(K)

              end do ! kU1
            end do ! kF1
          end do ! kU
        end do ! kF
                  end do ! kU3
                end do ! kL1
              end do ! kU2
            end do ! kL
          end do ! kLb1
        end do ! kLb
      end do ! K

      end subroutine emi2_inner_atomic_loop

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity and emissivity of a given LTE
      !! line in the presence of a magnetic field\n
      !!  line(LTEline_class): Structure with LTE line data\n
      !!  TKQ(dcomplx(:,:,:)): Geometrical tensors in the suitable
      !!                       reference frame\n
      !!     omega(double(:)): Frequency array\n
      !!   Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                       J-symbols\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!    aprof(Prof_class): Voigt profiles\n
      !!           Dw(double): Doppler width of the transition\n
      !!         vfac(double): Doppler shift factor\n
      !!    Bstrength(double): Magnetic field strength\n
      !!           pE(double): Unit transformation factor\n
      !!      eta0(double(:)): Intensity absorptivity\n
      !!      eta1(double(:)): Q absorptivity\n
      !!      eta2(double(:)): U absorptivity\n
      !!      eta3(double(:)): V absorptivity\n
      !!      rha1(double(:)): Q dichroic absorptivity\n
      !!      rha2(double(:)): U dichroic absorptivity\n
      !!      rha3(double(:)): V dichroic absorptivity\n
      !!      eps0(double(:)): Intensity emissivity\n
      !!      eps1(double(:)): Q emissivity\n
      !!      eps2(double(:)): U emissivity\n
      !!      eps3(double(:)): V emissivity\n
      !!      rhs1(double(:)): Q 'dichroic' emissivity\n
      !!      rhs2(double(:)): U 'dichroic' emissivity\n
      !!      rhs3(double(:)): V 'dichroic' emissivity
      subroutine rt1ordLTE(line,TKQ,omega,Flgsg,iz,if0,if1,aprof,Dw, &
                           vfac,Bstrength,pE, &
                           eta0,eta1,eta2,eta3,rha1,rha2,rha3, &
                           eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Prof_class), intent(in):: aprof
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw,pE,vfac,Bstrength
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0,eta1
      double precision, dimension(if0:if1), intent(out):: eta2,eta3
      double precision, dimension(if0:if1), intent(out):: eps0,eps1
      double precision, dimension(if0:if1), intent(out):: eps2,eps3
      double precision, dimension(if0:if1), intent(out):: rha1,rha2
      double precision, dimension(if0:if1), intent(out):: rhs1,rhs2
      double precision, dimension(if0:if1), intent(out):: rha3,rhs3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQ

      ! Local

      integer:: ifreq,iMu,iMl,iq,icom

      double precision:: at,feta,feps,Dfreqw,dnubw,vfacw,iDw
      double precision:: rMu,rMl,ftmp,q,dtK0,dtK1,dtK2,dtK3

      complex(kind=8):: tK0,tK1,tK2,tK3,prof


      !
      ! Initialize variables
      !
      eta0 = 0d0
      eta1 = 0d0
      eta2 = 0d0
      eta3 = 0d0
      rha1 = 0d0
      rha2 = 0d0
      rha3 = 0d0

      ! Inverse Doppler width
      iDw = 1d0/Dw


      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE*iDw

      ! Emissivity factor
      feps = line%nu(iz)*1d3*IPI41*line%Aul*iDw

      ! Line quantities

      ! Damping parameter
      at = line%damp(iz)*iDw

      ! Energy
      Dfreqw = (line%Eu - line%El)*iDw

      ! Shift
      vfacw = vfac*iDw

      ! Initialize component index
      icom = 0

      ! For each Mu
      do iMu=1,line%nMu

        ! Mu
        rMu = -line%Ju + dble(iMu-1)

        ! For each Ml
        do iMl=1,line%nMl

          ! Ml
          rMl = -line%Jl + dble(iMl-1)

          ! q value
          q = rMl - rMu
          iq = nint(q)

          ! Selection rule
          if (abs(q).gt.1) cycle

          ! 3J
          ftmp = fun3j(line%Ju,line%Jl,1d0,-rMu,rMl,-q,Flgsg)
          ftmp = ftmp*ftmp
          ftmp = ftmp*Flgsg%sg(1+iq)*sqrt3

          ! K = 0
          tK0 = TKQ(0,0,0)*fun3j(1d0,1d0,0d0,q,-q,0d0,Flgsg)

          ! K = 1
          tK3 = sqrt3*TKQ(3,0,1)*fun3j(1d0,1d0,1d0,q,-q,0d0,Flgsg)

          ! K = 2
          tK2 = sqrt5*fun3j(1d0,1d0,2d0,q,-q,0d0,Flgsg)
          tK0 = tK0 + tK2*TKQ(0,0,2)
          tK1 = tK2*TKQ(1,0,2)
          tK2 = tK2*TKQ(2,0,2)

          ! Scale
          dtK0 = dble(tK0)*ftmp
          dtK1 = dble(tK1)*ftmp
          dtK2 = dble(tK2)*ftmp
          dtK3 = dble(tK3)*ftmp

          ! Advance component index
          icom = icom + 1

          ! If stored in RAM
          if (aprof%VRAM) then

            ! Copy profile
            eta0 = eta0 + dtK0*dble(aprof%cp(:,icom))
            eta1 = eta1 + dtK1*dble(aprof%cp(:,icom))
            eta2 = eta2 + dtK2*dble(aprof%cp(:,icom))
            eta3 = eta3 + dtK3*dble(aprof%cp(:,icom))
            rha1 = rha1 + dtK1*dimag(aprof%cp(:,icom))
            rha2 = rha2 + dtK2*dimag(aprof%cp(:,icom))
            rha3 = rha3 + dtK3*dimag(aprof%cp(:,icom))

          ! Not stored
          else

            ! Magnetic shift
            dnubw = B2LK*Bstrength* &
                    (line%gu*rMu - line%gl*rMl)*iDw

            ! For each frequency
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreqw + dnubw - omega(ifreq)*vfacw,at,prof)

              ! Compute
              eta0(ifreq) = eta0(ifreq) + dtK0*dble(prof)
              eta1(ifreq) = eta1(ifreq) + dtK1*dble(prof)
              eta2(ifreq) = eta2(ifreq) + dtK2*dble(prof)
              eta3(ifreq) = eta3(ifreq) + dtK3*dble(prof)
              rha1(ifreq) = rha1(ifreq) + dtK1*dimag(prof)
              rha2(ifreq) = rha2(ifreq) + dtK2*dimag(prof)
              rha3(ifreq) = rha3(ifreq) + dtK3*dimag(prof)

            end do ! frequencies

          end if ! Type of profile calculation

        end do ! Ml
      end do ! Mu

      ! Add units
      eps0 = eta0*feps
      eps1 = eta1*feps
      eps2 = eta2*feps
      eps3 = eta3*feps
      rhs1 = rha1*feps
      rhs2 = rha2*feps
      rhs3 = rha3*feps
      eta0 = eta0*feta
      eta1 = eta1*feta
      eta2 = eta2*feta
      eta3 = eta3*feta
      rha1 = rha1*feta
      rha2 = rha2*feta
      rha3 = rha3*feta

      return

      end subroutine rt1ordLTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity of a given LTE line in the
      !! presence of a magnetic field\n
      !!  line(LTEline_class): Structure with LTE line data\n
      !!  TKQ(dcomplx(:,:,:)): Geometrical tensors in the suitable
      !!                       reference frame\n
      !!     omega(double(:)): Frequency array\n
      !!   Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                       J-symbols\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!    aprof(Prof_class): Voigt profiles\n
      !!           Dw(double): Doppler width of the transition\n
      !!         vfac(double): Doppler shift factor\n
      !!    Bstrength(double): Magnetic field strength\n
      !!           pE(double): Unit transformation factor\n
      !!      eta0(double(:)): Intensity absorptivity\n
      !!      eta1(double(:)): Q absorptivity\n
      !!      eta2(double(:)): U absorptivity\n
      !!      eta3(double(:)): V absorptivity\n
      !!      rha1(double(:)): Q dichroic absorptivity\n
      !!      rha2(double(:)): U dichroic absorptivity\n
      !!      rha3(double(:)): V dichroic absorptivity
      subroutine absorbLTE(line,TKQ,omega,Flgsg,iz,if0,if1,aprof,Dw, &
                           vfac,Bstrength,pE, &
                           eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Prof_class), intent(in):: aprof
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw, pE, vfac, Bstrength
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0,eta1
      double precision, dimension(if0:if1), intent(out):: eta2,eta3
      double precision, dimension(if0:if1), intent(out):: rha1,rha2
      double precision, dimension(if0:if1), intent(out):: rha3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQ

      ! Local

      integer:: ifreq,iMu,iMl,iq,icom

      double precision:: at,feta,Dfreqw,dnubw,vfacw,iDw
      double precision:: rMu,rMl,ftmp,q,dtK0,dtK1,dtK2,dtK3

      complex(kind=8):: tK0,tK1,tK2,tK3,prof


      !
      ! Initialize variables
      !

      eta0 = 0d0
      eta1 = 0d0
      eta2 = 0d0
      eta3 = 0d0
      rha1 = 0d0
      rha2 = 0d0
      rha3 = 0d0

      ! Inverse Doppler width
      iDw = 1d0/Dw


      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE*iDw

      ! Line quantities

      ! Damping parameter
      at = line%damp(iz)*iDw

      ! Energy
      Dfreqw = (line%Eu - line%El)*iDw

      ! Shift
      vfacw = vfac*iDw

      ! Initialize component index
      icom = 0

      ! For each Mu
      do iMu=1,line%nMu

        ! Mu
        rMu = -line%Ju + dble(iMu-1)

        ! For each Ml
        do iMl=1,line%nMl

          ! Ml
          rMl = -line%Jl + dble(iMl-1)

          ! q value
          q = rMl - rMu
          iq = nint(q)

          ! Selection rule
          if (abs(q).gt.1) cycle

          ! 3J
          ftmp = fun3j(line%Ju,line%Jl,1d0,-rMu,rMl,-q,Flgsg)
          ftmp = ftmp*ftmp
          ftmp = ftmp*Flgsg%sg(1+iq)*sqrt3

          ! K = 0
          tK0 = TKQ(0,0,0)*fun3j(1d0,1d0,0d0,q,-q,0d0,Flgsg)

          ! K = 1
          tK3 = sqrt3*TKQ(3,0,1)*fun3j(1d0,1d0,1d0,q,-q,0d0,Flgsg)

          ! K = 2
          tK2 = sqrt5*fun3j(1d0,1d0,2d0,q,-q,0d0,Flgsg)
          tK0 = tK0 + tK2*TKQ(0,0,2)
          tK1 = tK2*TKQ(1,0,2)
          tK2 = tK2*TKQ(2,0,2)

          ! Scale
          dtK0 = dble(tK0)*ftmp
          dtK1 = dble(tK1)*ftmp
          dtK2 = dble(tK2)*ftmp
          dtK3 = dble(tK3)*ftmp

          ! Advance component index
          icom = icom + 1

          ! If stored in RAM
          if (aprof%VRAM) then

            ! Copy profiles
            eta0 = eta0 + dtK0*dble(aprof%cp(:,icom))
            eta1 = eta1 + dtK1*dble(aprof%cp(:,icom))
            eta2 = eta2 + dtK2*dble(aprof%cp(:,icom))
            eta3 = eta3 + dtK3*dble(aprof%cp(:,icom))
            rha1 = rha1 + dtK1*dimag(aprof%cp(:,icom))
            rha2 = rha2 + dtK2*dimag(aprof%cp(:,icom))
            rha3 = rha3 + dtK3*dimag(aprof%cp(:,icom))

          ! Not stored
          else

            ! Magnetic shift
            dnubw = B2LK*Bstrength* &
                    (line%gu*rMu - line%gl*rMl)*iDw

            ! For each frequency
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreqw + dnubw - omega(ifreq)*vfacw,at,prof)

              ! Add contribution
              eta0(ifreq) = eta0(ifreq) + dtK0*dble(prof)
              eta1(ifreq) = eta1(ifreq) + dtK1*dble(prof)
              eta2(ifreq) = eta2(ifreq) + dtK2*dble(prof)
              eta3(ifreq) = eta3(ifreq) + dtK3*dble(prof)
              rha1(ifreq) = rha1(ifreq) + dtK1*dimag(prof)
              rha2(ifreq) = rha2(ifreq) + dtK2*dimag(prof)
              rha3(ifreq) = rha3(ifreq) + dtK3*dimag(prof)

            end do ! frequencies

          end if ! Type of profile calculation

        end do ! Ml
      end do ! Mu

      ! Add units
      eta0 = eta0*feta
      eta1 = eta1*feta
      eta2 = eta2*feta
      eta3 = eta3*feta
      rha1 = rha1*feta
      rha2 = rha2*feta
      rha3 = rha3*feta

      return

      end subroutine absorbLTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity of a given recombination
      !! transition\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  omega(double(:)): Frequency array\n
      !!         T(double): Temperature\n
      !!        ne(double): Electron number density\n
      !!    itran(integer): Transition index\n
      !!  ilevelu(integer): Upper level index\n
      !!       iz(integer): Height index\n
      !!      if0(integer): First frequency index for this
      !!                    transition\n
      !!      if1(integer): Last frequency index for this
      !!                    transition\n
      !!    eps(double(:)): Emissivity\n
      !!    eta(double(:)): Stimulated emissivity
      subroutine photoeps(Atom,omega,T,ne,itran,ilevelu,iz,if0,if1, &
                          eps,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran,ilevelu,iz,if0,if1
      double precision, intent(in):: T,ne
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq,iJu,itermu,iR

      double precision:: c0,c1,exu,pE,Saha,rJu,rhou
      double precision, dimension(if0:if1):: omega3


      !
      ! Saha term constants
      !
      c0 = fktoJ/kb/T
      Saha = cSaha*ne*Atom%phot(itran)%glu* &
             exp(Atom%phot(itran)%edge*c0)/(T**(1.5d0))

      !
      ! Indexes
      !

      ! Get term index
      itermu = Atom%term(ilevelu)

      ! Get J level index
      iJu = Atom%sublevel(ilevelu)

      ! Get angular momentum
      rJu = Atom%rJval(iJu,itermu)

      ! Get SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! Get upper level population
      rhou = sqrt(2d0*rJu+1d0)*dble(Atom%crho(iR,iz))

      ! Apply Saha factor
      rhou = rhou*Saha

      ! Compute exponential argument constant
      c0 = c2*1d4/T

      ! Compute energy constant part
      c1 = 2d21*c*convF

      ! Cubic frequency
      omega3 = omega(if0:if1)
      omega3 = omega3*omega3*omega3

      ! For each frequency
      do ifreq=if0,if1

        ! Exponential
        exu = c0*omega(ifreq)
        exu = diexp(exu)

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu*rhou

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies

      return

      end subroutine photoeps

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity of a given recombination transition
      !! with frequency quantities stored in RAM\n
      !!   Atom(Atom_class): Structure with the atomic data\n
      !!  omega3(double(:)): Frequency array to the power of 3\n
      !!     exu(double(:)): Exponential in the stimulated part of the
      !!                     photoionization emissivity\n
      !!          T(double): Temperature\n
      !!         ne(double): Electron number density\n
      !!     itran(integer): Transition index\n
      !!   ilevelu(integer): Upper level index\n
      !!        iz(integer): Height index\n
      !!       if0(integer): First frequency index for this
      !!                     transition\n
      !!       if1(integer): Last frequency index for this
      !!                     transition\n
      !!     eps(double(:)): Emissivity\n
      !!     eta(double(:)): Stimulated emissivity
      subroutine photoepsS(Atom,omega3,exu,T,ne,itran, &
                           ilevelu,iz,if0,if1,eps,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran,ilevelu,iz,if0,if1
      double precision, intent(in):: T,ne
      double precision, dimension(if0:if1), intent(in):: omega3
      double precision, dimension(if0:if1), intent(in):: exu
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq,iJu,itermu,iR

      double precision:: c0,c1,pE,Saha,rJu,rhou


      !
      ! Saha term constants
      !
      c0 = fktoJ/kb/T
      Saha = cSaha*ne*Atom%phot(itran)%glu* &
             exp(Atom%phot(itran)%edge*c0)/(T**(1.5d0))

      !
      ! Indexes
      !

      ! Get term index
      itermu = Atom%term(ilevelu)

      ! Get J level index
      iJu = Atom%sublevel(ilevelu)

      ! Get angular momentum
      rJu = Atom%rJval(iJu,itermu)

      ! Get SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! Get upper level population
      rhou = sqrt(2d0*rJu+1d0)*dble(Atom%crho(iR,iz))

      ! Apply Saha factor
      rhou = rhou*Saha

      ! Compute exponential argument constant
      c0 = c2*1d4/T

      ! Compute energy constant part
      c1 = 2d21*c*convF

      ! For each frequency
      do ifreq=if0,if1

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu(ifreq)*rhou

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies

      return

      end subroutine photoepsS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the frequency dependent JKQ for the angle-averaged
      !! second order emissivity in the presence of velocities in the
      !! comoving frame\n
      !!          Red(Redb_class): Structure with redistribution input
      !!                           frequency data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!              iz(integer): Height index\n
      !!           ntran(integer): Number of transitions in the atom\n
      !!         tif0(integer(:)): Lower limit to search in Stokes
      !!                           interpolation\n
      !!         tif1(integer(:)): Upper limit to search in Stokes
      !!                           interpolation\n
      !!              DwT(double): Thermal part of the Doppler width\n
      !!               vx(double): Velocity vector along X\n
      !!               vy(double): Velocity vector along Y\n
      !!               vz(double): Velocity vector along Z\n
      !!         omega(double(:)): Frequency array\n
      !!       Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                           and J-symbols\n
      !!  Stokes(double(:,:,:,:)): Stokes parameters\n
      !!    JKQa(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                           field tensors\n
      !!   JRadC(dcomplex(:,:,:)): Frequency dependent JKQ in the
      !!                           comoving frame
      subroutine getJKQstar(Red,Geom,iz,ntran,tif0,tif1,DwT,vx,vy, &
                            vz,omega,Flgsg,TKQo,Stokes,JKQa,JradC)

      ! I/O

      type(Redb_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,ntran
      integer, dimension(:), intent(in):: tif0,tif1
      double precision, intent(in):: DwT,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(:,:,:), &
                       allocatable, intent(out):: JradC
      complex(kind=8), dimension(0:3,-2:2,0:2,1:Geom%nth*Geom%nph2), &
                       intent(in):: TKQo

      ! Local

      logical:: shift,asym

      integer:: ifreq,ith1,K,iph1,iQ,jz,jdir,jtran,if0,if1

      double precision:: vfac1,omegao,cost,sint,cosc,sinc
      double precision, dimension(0:3):: StokesM


      ! If velocity is above threshold
      shift = (vx*vx + vy*vy + vz*vz)*1d6*c.ge.vrfrac*DwT
      asym = size(JKQa).gt.10

      ! Allocate
      allocate(JradC(-2:2,0:2,Red%ggf0:Red%ggf1))

      !
      ! For each frequency, get mean intensity
      !

      ! Initialize
      vfac1 = 1d0
      JradC = cZero

      !
      ! If there are ad-hoc asymmetries
      !
      if (asym) then

        ! If forcing the ad-hoc asymmetry only anisotropy in vertical
        ! frame
        if (force_asym.or.axial) then

          ! Relevant frequencies
          do ifreq=Red%ggf0,Red%ggf1

            ! Initialize
            if0 = 1
            if1 = nfreq

            ! Find transition for frequency limits
            do jtran=1,ntran

              ! If out of limits, skip
              if (ifreq.lt.tif0(jtran)) cycle
              if (ifreq.gt.tif1(jtran)) cycle

              ! Found
              if0 = tif0(jtran)
              if1 = tif1(jtran)
              exit

            end do

            ! Initialize direction index
            jdir = 0

            ! For each polar direction
            do ith1=1,Geom%nTh

              ! Get director cosines
              if (shift) cost = Geom%V_mu(ith1)

              ! If axial
              if (axial) then

                ! Advance direction
                jdir = jdir + 1

                ! Calculate Doppler shift factor
                if (shift) then

                  ! Get shift
                  vfac1 = 1d0 + vz*cost

                  ! Get target frequency
                  omegao = omega(ifreq)*vfac1

                  ! Interpolate
                  StokesM = getStkinnu(omega, &
                                       Stokes(:,:,1,ith1), &
                                       ifreq,if0,if1,omegao)

                ! No shift
                else

                  ! Get value
                  StokesM = Stokes(:,ifreq,1,ith1)

                end if ! shift

                ! For each multipole
                do K=0,Krad

                  ! Add to integral
                  JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                     Geom%W_mu(ith1)* &
                                   sum(StokesM* &
                                       Geom%TS(:,0,K,jdir))
                end do ! Multipoles

                ! Advance non-used directiond
                jdir = jdir + Geom%nPh2 - 1

              ! Non-axial
              else

                ! Get other factor in inclination for shift
                if (shift) sint = sqrt(1d0 - cost*cost)

                ! For each azimuth
                do iph1=1,Geom%nPh

                  ! Advance direction
                  jdir = jdir + 1

                  ! If there is a significant enough shift
                  if (shift) then

                    ! Get total shift (inverse)
                    cosc = Geom%v_mux(iph1)
                    sinc = Geom%v_muy(iph1)* &
                           sqrt(1d0 - cosc*cosc)
                    vfac1 = 1d0 + vx*sint*cosc + &
                                  vy*sint*sinc + &
                                  vz*cost

                    ! Target frequency
                    omegao = omega(ifreq)*vfac1

                    ! Interpolate
                    StokesM = getStkinnu(omega, &
                                         Stokes(:,:,iph1,ith1), &
                                         ifreq,if0,if1,omegao)

                  ! No shift
                  else

                    ! Value
                    StokesM = Stokes(:,ifreq,iph1,ith1)

                  end if ! shift

                  ! For every multipole
                  do K=0,Krad

                    ! Add to integral
                    JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                       Geom%W_mu(ith1)* &
                                       Geom%W_mux(iph1)* &
                                       sum(StokesM* &
                                           Geom%TS(:,0,K,jdir))
                  end do ! Multipoles
                end do ! Azimuths

              end if ! Axial

            end do ! Polar directions
          end do ! Relevant frequencies

        ! Additive ad-hoc asymmetry
        else

          ! Relevant frequencies
          do ifreq=Red%ggf0,Red%ggf1

            ! Initialize
            if0 = 1
            if1 = nfreq

            ! Find transition for frequency limits
            do jtran=1,ntran

              ! If out of limits, skip
              if (ifreq.lt.tif0(jtran)) cycle
              if (ifreq.gt.tif1(jtran)) cycle

              ! Found
              if0 = tif0(jtran)
              if1 = tif1(jtran)
              exit

            end do

            ! Initialize direction
            jdir = 0

            ! For each polar direction
            do ith1=1,Geom%nTh

              ! Compute z shift
              if (shift) then
                cost = Geom%V_mu(ith1)
                sint = sqrt(1d0 - cost*cost)
              end if

              ! For each azimuth
              do iph1=1,Geom%nPh

                ! Advance direction
                jdir = jdir + 1

                ! If the shift is significant enough
                if (shift) then

                  ! Get shift (inverse)
                  cosc = Geom%v_mux(iph1)
                  sinc = Geom%v_muy(iph1)* &
                         sqrt(1d0 - cosc*cosc)
                  vfac1 = 1d0 + vx*sint*cosc + &
                                vy*sint*sinc + &
                                vz*cost

                  ! Target frequency
                  omegao = omega(ifreq)*vfac1

                  ! Interpolate
                  StokesM = getStkinnu(omega, &
                                       Stokes(:,:,iph1,ith1), &
                                       ifreq,if0,if1,omegao)

                ! No shift
                else

                  ! Value
                  StokesM = Stokes(:,ifreq,iph1,ith1)

                end if ! Shift

                ! For every independent multipole
                do K=0,Krad
                  do iQ=0,K

                    ! Add to integral
                    JradC(iQ,K,ifreq) = JradC(iQ,K,ifreq) + &
                                        Geom%W_mu(ith1)* &
                                        Geom%W_mux2(iph1)* &
                                        sum(StokesM* &
                                            Geom%TS(:,iQ,K,jdir))
                  end do ! Q
                end do ! K
              end do ! Azimuth
            end do ! Polar
          end do ! Frequencies

        end if ! Type of ad-hoc asymmetry

        ! Current height for JKQa
        jz = iz - Rz0 + 1

        !
        ! Now add ad-hoc assymetry
        !

        ! Relevant frequencies
        do ifreq=Red%ggf0,Red%ggf1

          ! Add ad-hoc asymmetry
          JradC(0:2,1:2,ifreq) = JradC(0:2,1:2,ifreq) + &
                                 JKQa(3:5,:,jz)*dble(JradC(0,0,ifreq))

        end do ! Relevant frequencies

      ! No Ad-hoc asymmetries
      else

        ! For each frequency, get JKQ
        do ifreq=Red%ggf0,Red%ggf1

          ! Initialize
          if0 = 1
          if1 = nfreq

          ! Find transition for frequency limits
          do jtran=1,ntran

            ! If out of limits, skip
            if (ifreq.lt.tif0(jtran)) cycle
            if (ifreq.gt.tif1(jtran)) cycle

            ! Found
            if0 = tif0(jtran)
            if1 = tif1(jtran)
            exit

          end do

          ! Initialize directions
          jdir = 0

          ! For each polar direction
          do ith1=1,Geom%nTh

            ! Get director cosin
            if (shift) cost = Geom%V_mu(ith1)

            ! If axial symmetric
            if (axial) then

              ! Advance direction
              jdir = jdir + 1

              ! If sifnigicant Doppler shift
              if (shift) then

                ! Get shift
                vfac1 = 1d0 + vz*cost

                ! Target frequency
                omegao = omega(ifreq)*vfac1

                ! Interpolate
                StokesM = getStkinnu(omega, &
                                     Stokes(:,:,1,ith1), &
                                     ifreq,if0,if1,omegao)

              ! No shift
              else

                ! Value
                StokesM = Stokes(:,ifreq,1,ith1)

              end if ! shift

              ! For every multipole
              do K=0,Krad

                ! Add to integral
                JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                   Geom%W_mu(ith1)* &
                                   sum(StokesM*TKQo(:,0,K,jdir))
              end do ! Multipoles

              ! Advance non-used directions
              jdir = jdir + Geom%nPh2 - 1

            ! Non axial symmetric
            else

              ! If there is shift, get sin
              if (shift) sint = sqrt(1d0 - cost*cost)

              ! For each azimuth
              do iph1=1,Geom%nPh

                ! Advance direction
                jdir = jdir + 1

                ! If significant enough shift
                if (shift) then

                  ! Get Doppler shift (inverse)
                  cosc = Geom%v_mux(iph1)
                  sinc = Geom%v_muy(iph1)* &
                         sqrt(1d0 - cosc*cosc)
                  vfac1 = 1d0 + vx*sint*cosc + &
                                vy*sint*sinc + &
                                vz*cost

                  ! Target frequency
                  omegao = omega(ifreq)*vfac1

                  ! Interpolate
                  StokesM = getStkinnu(omega, &
                                       Stokes(:,:,iph1,ith1), &
                                       ifreq,if0,if1,omegao)

                ! No shift
                else

                  ! Value
                  StokesM = Stokes(:,ifreq,iph1,ith1)

                end if ! Shift

                ! For every independent multipole
                do K=0,Krad
                  do iQ=0,K

                    ! Add to integral
                    JradC(iQ,K,ifreq) = JradC(iQ,K,ifreq) + &
                                          Geom%W_mu(ith1)* &
                                          Geom%W_mux2(iph1)* &
                                          sum(StokesM* &
                                              TKQo(:,iQ,K,jdir))
                  end do ! K
                end do ! Q
              end do ! Azimuth

            end if ! Axial symmetry

          end do ! Polar
        enddo ! Frequencies

      end if ! Are there ad-hoc asymmetries or not

      ! If not axial 
      if (.not.axial.or.force_asym) then

        ! For every dependent multipole
        do K=1,Krad
          do iQ=0,K

            ! Use relations
            JradC(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JradC(iQ,K,:))

          end do ! Q
        end do ! K

      end if ! Not axial

      end subroutine getJKQstar

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the frequency dependent JKQ for the angle-dependent
      !! second order emissivity\n
      !!          Red(Redb_class): Structure with redistribution input
      !!                           frequency data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!              iz(integer): Height index\n
      !!           ntran(integer): Number of transitions in the atom\n
      !!         tif0(integer(:)): Lower limit to search in Stokes
      !!                           interpolation\n
      !!         tif1(integer(:)): Upper limit to search in Stokes
      !!                           interpolation\n
      !!              DwT(double): Thermal part of the Doppler width\n
      !!               vx(double): Velocity vector along X\n
      !!               vy(double): Velocity vector along Y\n
      !!               vz(double): Velocity vector along Z\n
      !!          lfield(logical): If magnetic field\n
      !!           Btheta(double): Magnetic field polar angle\n
      !!             Bphi(double): Magnetic field azimuth angle\n
      !!         omega(double(:)): Frequency array\n
      !!       Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                           and J-symbols\n
      !!  Stokes(double(:,:,:,:)): Stokes parameters\n
      !!    JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                           frequency dependence\n
      !!    JKQa(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                           field tensors\n
      !!   JRadC(dcomplex(:,:,:)): Frequency dependent JKQ in the
      !!                           comoving frame
      subroutine getJKQADasym(Red,Geom,iz,ntran,tif0,tif1,DwT,vx,vy, &
                              vz,lfield,Btheta,Bphi,omega,Flgsg, &
                              TKQo,Stokes,JKQC,JKQa,JradC)

      ! I/O

      type(Redb_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: lfield
      integer, intent(in):: iz,ntran
      integer, dimension(:), intent(in):: tif0,tif1
      double precision, intent(in):: DwT,vx,vy,vz,Btheta,Bphi
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nfreq), &
                       target, intent(in):: JKQC
      complex(kind=8), dimension(:,:,:), &
                       allocatable, intent(out):: JradC
      complex(kind=8), dimension(0:3,-2:2,0:2,1:Geom%nth*Geom%nph2), &
                       intent(in):: TKQo

      ! Local

      logical:: shift

      integer:: ifreq,ith1,K,iph1,iQ,jz,jdir,jtran,if0,if1

      double precision:: vfac1,omegao,cost,sint,cosc,sinc
      double precision, dimension(0:3):: StokesM
      complex(kind=8), dimension(:,:,:), pointer:: JradCp


      ! If velocity is above threshold
      shift = (vx*vx + vy*vy + vz*vz)*1d6*c.ge.vrfrac*DwT

      ! Allocate
      allocate(JradC(-2:2,0:2,Red%ggf0:Red%ggf1))
      JradC = cZero
      if (lfield.or.shift) then
        allocate(JradCp(-2:2,0:2,Red%ggf0:Red%ggf1))
        JradCp = cZero
      end if

      ! Height for JKQa
      jz = iz - Rz0 + 1

      ! If dynamic
      if (shift) then

        ! Initialize
        vfac1 = 1d0

        ! For each frequency, get JKQ
        do ifreq=Red%ggf0,Red%ggf1

          ! Initialize
          if0 = 1
          if1 = nfreq

          ! Find transition for frequency limits
          do jtran=1,ntran

            ! If out of limits, skip
            if (ifreq.lt.tif0(jtran)) cycle
            if (ifreq.gt.tif1(jtran)) cycle

            ! Found
            if0 = tif0(jtran)
            if1 = tif1(jtran)
            exit

          end do

          ! Initialize directions
          jdir = 0

          ! For each polar direction
          do ith1=1,Geom%nTh

            ! Get director cosin
            cost = Geom%V_mu(ith1)

            ! If axial symmetric
            if (axial) then

              ! Advance direction
              jdir = jdir + 1

              ! Get shift
              vfac1 = 1d0 + vz*cost

              ! Target frequency
              omegao = omega(ifreq)*vfac1

              ! Interpolate
              StokesM = getStkinnu(omega, &
                                   Stokes(:,:,1,ith1), &
                                   ifreq,if0,if1,omegao)

              ! Value
              StokesM = Stokes(:,ifreq,1,ith1)

              ! For every multipole
              do K=0,Krad

                ! Add to integral
                JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                   Geom%W_mu(ith1)* &
                                   sum(StokesM*Geom%TS(:,0,K,jdir))
                JradCp(0,K,ifreq) = JradCp(0,K,ifreq) + &
                                    Geom%W_mu(ith1)* &
                                    sum(StokesM*TKQo(:,0,K,jdir))
              end do ! Multipoles

              ! Advance non-used directions
              jdir = jdir + Geom%nPh2 - 1

            ! Non axial symmetric
            else

              ! If there is shift, get sin
              sint = sqrt(1d0 - cost*cost)

              ! For each azimuth
              do iph1=1,Geom%nPh

                ! Advance direction
                jdir = jdir + 1

                ! Get Doppler shift (inverse)
                cosc = Geom%v_mux(iph1)
                sinc = Geom%v_muy(iph1)* &
                       sqrt(1d0 - cosc*cosc)
                vfac1 = 1d0 + vx*sint*cosc + &
                              vy*sint*sinc + &
                              vz*cost

                ! Target frequency
                omegao = omega(ifreq)*vfac1

                ! Interpolate
                StokesM = getStkinnu(omega, &
                                     Stokes(:,:,iph1,ith1), &
                                     ifreq,if0,if1,omegao)

                ! For every independent multipole
                do K=0,Krad
                  do iQ=0,K

                    ! Add to integral
                    JradC(iQ,K,ifreq) = JradC(iQ,K,ifreq) + &
                                          Geom%W_mu(ith1)* &
                                          Geom%W_mux2(iph1)* &
                                          sum(StokesM* &
                                              Geom%TS(:,iQ,K,jdir))
                    JradCp(iQ,K,ifreq) = JradCp(iQ,K,ifreq) + &
                                         Geom%W_mu(ith1)* &
                                         Geom%W_mux2(iph1)* &
                                         sum(StokesM* &
                                             TKQo(:,iQ,K,jdir))
                  end do ! K
                end do ! Q
              end do ! Azimuth

            end if ! Axial symmetry

          end do ! Polar
        end do ! Frequencies

        ! Force asymmetry
        if (force_asym) JradC(1:2,:,:) = cZero

        ! For each frequency
        do ifreq=Red%ggf0,Red%ggf1

          ! Add ad-hoc
          JradC(0:2,1:2,ifreq) = JradC(0:2,1:2,ifreq) + &
                                 JKQa(3:5,:,jz)*dble(JradC(0,0,ifreq))

        end do ! Relevant frequencies

        ! Ensure relations for every dependent multipole
        do K=1,Krad
          do iQ=1,K

            ! Use relations
            JradC(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JradC(iQ,K,:))
            JradCp(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JradCp(iQ,K,:))

          end do ! Q
        end do ! K

        ! If magnetic field, rotate
        if (lfield) &
          call fieldB_alt(JradC,Red%ggf1-Red%ggf0+1, &
                          Flgsg,Btheta,Bphi,1)

        ! Get difference
        JradC = JradC - JradCp

        ! Free
        deallocate(JradCp)
        nullify(JradCp)

      ! If static
      else

        ! Initialize
        JradC = JKQC(:,:,Red%ggf0:Red%ggf1)

        ! Force asymmetry
        if (force_asym) JradC(1:2,:,:) = cZero

        ! For each frequency
        do ifreq=Red%ggf0,Red%ggf1

          ! Add ad-hoc
          JradC(0:2,1:2,ifreq) = JradC(0:2,1:2,ifreq) + &
                                 JKQa(3:5,:,jz)*dble(JradC(0,0,ifreq))

        end do ! Relevant frequencies

        ! Ensure relations for every dependent multipole
        do K=1,Krad
          do iQ=1,K

            ! Use relations
            JradC(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JradC(iQ,K,:))

          end do ! Q
        end do ! K

        ! If magnetic field
        if (lfield) then

          ! Copy
          JradCp = JKQC(:,:,Red%ggf0:Red%ggf1)

          ! Rotate both
          call fieldB_alt(JradC,Red%ggf1-Red%ggf0+1, &
                          Flgsg,Btheta,Bphi,1)
          call fieldB_alt(JradCp,Red%ggf1-Red%ggf0+1, &
                          Flgsg,Btheta,Bphi,1)

        ! No magnetic field
        else

          ! Point
          JradCp(-2:2,0:2,Red%ggf0:Red%ggf1) => &
                                          JKQC(:,:,Red%ggf0:Red%ggf1)
        end if

        ! Get difference
        JradC = JradC - JradCp

        !
        ! Free
        !

        ! Magnetic
        if (lfield) then

          ! Deallocate and nullify
          deallocate(JradCp)
          nullify(JradCp)

        ! No magnetic field
        else

          ! Nullify
          nullify(JradCp)

        end if ! Magnetic or not
      end if ! Dynamic or static

      end subroutine getJKQADasym

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the Stokes parameters into the requested
      !! frequency\n
      !!     omega(double(:)): Frequency array\n
      !!  Stokes(double(:,:)): Stokes parameters\n
      !!       ifreq(integer): Frequency index of the output frequency
      !!                       associated to the requested input
      !!                       frequency\n
      !!         if0(integer): Lower limit to search in Stokes
      !!                       interpolation\n
      !!         if1(integer): Upper limit to search in Stokes
      !!                       interpolation\n
      !!            x(double): Frequency to interpolate into
      function getStkinnu(omega,Stokes,ifreq,if0,if1,x)

      ! I/O

      integer, intent(in):: ifreq,if0,if1
      double precision, intent(in):: x
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq), intent(in):: Stokes
      double precision, dimension(0:3):: getStkinnu

      ! Local

      integer:: jfreq

      double precision, dimension(0:3):: dxs, dys


      ! Initialize as equals
      getStkinnu = Stokes(:,ifreq)

      ! If omegai > omega(ifreq)
      if (x.gt.omega(ifreq)) then

        ! If out of right boundary
        if (x.ge.omega(if1)-TINYO) then

          ! Get boundary value
          getStkinnu = Stokes(:,if1)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,if1-1

            ! If this exact frequency is in output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              ! Copy
              getStkinnu = Stokes(:,jfreq)
              return

            ! If the input is between this output and the next
            else if (x.ge.omega(jfreq).and. &
                     x.lt.omega(jfreq+1)) then

              ! Linear interpolation
              dys = Stokes(:,jfreq+1) - Stokes(:,jfreq)
              dxs = x - omega(jfreq)
              getStkinnu = dxs*dys/(omega(jfreq+1) - omega(jfreq)) + &
                         Stokes(:,jfreq)
              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries

      ! If omegai < omegao
      else if (x.lt.omega(ifreq)) then

        ! If out of left boundary
        if (x.le.omega(if0)+TINYO) then

          ! Get boundary value
          getStkinnu = Stokes(:,if0)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,if0+1,-1

            ! If this exact frequency is in output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getStkinnu = Stokes(:,jfreq)

              return

            ! If the input is between this output and the next
            else if (x.ge.omega(jfreq-1).and. &
                     x.lt.omega(jfreq)) then

              ! Linear interpolation
              dys = Stokes(:,jfreq) - Stokes(:,jfreq-1)
              dxs = x - omega(jfreq-1)
              getStkinnu = dxs*dys/(omega(jfreq) - omega(jfreq-1)) + &
                         Stokes(:,jfreq-1)
              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries
      end if ! omegai > omegao

      end function getStkinnu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the Stokes parameters into the input frequency
      !! axis\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!     p_red(Redc_class): Structure with redistribution input
      !!                        frequency data for a given input
      !!                        transition\n
      !!       Fed(Reda_class): Structure with redistribution
      !!                        output frequency data\n
      !!       Red(Redb_class): Structure with redistribution input
      !!                        frequency data\n
      !!         Mif0(integer): First frequency for this CPU\n
      !!         Mif1(integer): Last frequency for this CPU\n
      !!      omega(double(:)): Frequency array\n
      !!            vx(double): Velocity vector along X\n
      !!            vy(double): Velocity vector along Y\n
      !!            vz(double): Velocity vector along Z\n
      !!         lvel(logical): If dynamic node\n
      !!    Stkin(double(:,:)): Interpolated Stokes parameters\n
      !!  Stk(double(:,:,:,:)): Original Stokes parameters
      subroutine getStkin(Geom,p_red,Fed,Red,Mif0,Mif1,omega, &
                          vx,vy,vz,lvel,Stkin,Stk)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(Redc_class), intent(in), pointer:: p_red
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      logical, intent(in):: lvel
      integer, intent(in):: Mif0,Mif1
      double precision, intent(in):: vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:,:,:,:), intent(in):: Stk
      double precision, dimension(:,:), &
                        allocatable, intent(out):: Stkin

      ! Local

      integer:: ith1,iph1,nmfreq,iran,ifreq,iifreq,lifreq,ibfreq
      integer:: jfreq,jjfreq0,jjfreq,kkfreq0,kkfreq,nblock

      double precision:: dx,vfac1,cost,sint,cosc,sinc
      double precision, dimension(4):: y0,dy

      ! Pointer

      integer, pointer:: p_mfreq


      ! Nullify
      nullify(p_mfreq)

      ! Get J size
      nmfreq = 0

      ! Run over all output frequencies
      iifreq = 0
      do iran=1,Fed%nran
        do ifreq=Fed%if0(iran),Fed%if1(iran)

          ! Advance index
          iifreq = iifreq + 1

          ! MPI
          if (iifreq.gt.Mif1) exit
          if (iifreq.lt.Mif0) cycle

          ! Input frequency number
          p_mfreq => p_red%mfreq(iifreq)

          ! Add size
          if (p_mfreq.gt.0) nmfreq = nmfreq + p_mfreq

        end do
      end do ! Output frequencies

      ! If no size, return
      if (nmfreq.lt.1) return

      ! Allocate Stokes parameters
      if (axial) then
        allocate(Stkin(nmfreq*Geom%nTh,0:3))
      else
        allocate(Stkin(nmfreq*(Geom%nPh2*Geom%nTh),0:3))
      end if

      ! If axial
      if (axial) then

        ! If dynamic
        if (lvel) then

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! For each input direction
              do ith1=1,Geom%nth

                ! Get director cosines
                cost = Geom%V_mu(ith1)

                ! Calculate Doppler shift factor
                vfac1 = 1d0 - vz*cost

                ! We will be using the inverse
                vfac1 = 1d0/vfac1

                ! Reset the search frequency
                lifreq = Red%ggf0

                ! For each input frequency
                do jfreq=1,p_mfreq

                  ! Advance indexes
                  jjfreq = jjfreq0 + jfreq
                  kkfreq = kkfreq0 + jfreq

                  ! If out of range, take the value at the
                  ! boundary
                  if (p_red%omega(jjfreq)*vfac1.le. &
                      omega(Red%ggf0)+TINYO) then

                    ! We are still looking in the first one
                    lifreq = Red%ggf0

                    ! The index to take is 1
                    Stkin(kkfreq,:) = Stk(:,Red%ggf0,1,ith1)

                  ! If out of range, take the value at the
                  ! boundary
                  else if (p_red%omega(jjfreq)*vfac1.ge. &
                           (omega(Red%ggf1) - TINYO)) then

                    ! We are in the last frequency
                    lifreq = Red%ggf1

                    ! The index to take is nfreq
                    Stkin(kkfreq,:) = Stk(:,Red%ggf1,1,ith1)

                  ! If within the boundaries
                  else

                    ! Search between the last found frequency and
                    ! all but the boundary
                    do ibfreq=lifreq,nfreq-1

                      ! If this exact frequency is in output
                      if (abs(p_red%omega(jjfreq)*vfac1 - &
                              omega(ibfreq)).lt.TINYO) then

                        ! We are in the found frequency
                        lifreq = ibfreq

                        ! This frequency gives us the value
                        Stkin(kkfreq,:) = Stk(:,lifreq,1,ith1)

                        exit

                      ! If the input is between this output and
                      ! the next
                      else if(p_red%omega(jjfreq)*vfac1.ge. &
                              omega(ibfreq).and. &
                              p_red%omega(jjfreq)*vfac1.lt. &
                              omega(ibfreq+1)) then

                        ! We found it in the index of the lower
                        lifreq = ibfreq

                        ! The first index is the lower
                        y0 = Stk(:,lifreq,1,ith1)

                        ! The second index is the upper
                        dy = Stk(:,lifreq+1,1,ith1) - y0

                        ! Inverse of the distance
                        ! between the two outputs
                        dx = (p_red%omega(jjfreq)*vfac1 - &
                              omega(lifreq))/ &
                             (omega(lifreq+1) - omega(lifreq))

                        ! Interpolate
                        Stkin(kkfreq,:) = dx*dy + y0

                        exit

                      end if ! Check output frequency

                    end do ! Run output frequencies

                  end if ! Check if out of limits

                end do ! Run input frequencies

                ! Update indexes
                kkfreq0 = kkfreq0 + p_mfreq

              end do ! Input polar

              ! Update indexes
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! Output frequencies
          end do ! Output frequency ranges

        ! If static
        else

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! Reset the search frequency
              lifreq = Red%ggf0

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq

                ! If out of range, take the value at the
                ! boundary
                if (p_red%omega(jjfreq).le. &
                    omega(Red%ggf0)+TINYO) then

                  ! We are still looking in the first one
                  lifreq = Red%ggf0

                  ! Inclinations
                  do ith1=1,Geom%nTh

                    ! Output index
                    kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                    ! The index to take is 1
                    Stkin(kkfreq,:) = Stk(:,Red%ggf0,1,ith1)

                  end do

                ! If out of range, take the value at the
                ! boundary
                else if (p_red%omega(jjfreq).ge. &
                         (omega(Red%ggf1) - TINYO)) then

                  ! We are in the last frequency
                  lifreq = Red%ggf1

                  ! Inclinations
                  do ith1=1,Geom%nTh

                    ! Output index
                    kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                    ! The index to take is nfreq
                    Stkin(kkfreq,:) = Stk(:,Red%ggf1,1,ith1)

                  end do

                ! If within the boundaries
                else

                  ! Search between the last found frequency and
                  ! all but the boundary
                  do ibfreq=lifreq,nfreq-1

                    ! If this exact frequency is in output
                    if (abs(p_red%omega(jjfreq)- &
                            omega(ibfreq)).lt.TINYO) then

                      ! We are in the found frequency
                      lifreq = ibfreq

                      ! Inclinations
                      do ith1=1,Geom%nTh

                        ! Output index
                        kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                        ! The index to take is nfreq
                        Stkin(kkfreq,:) = Stk(:,lifreq,1,ith1)

                      end do

                      exit

                    ! If the input is between this output and
                    ! the next
                    else if(p_red%omega(jjfreq).ge. &
                            omega(ibfreq).and. &
                            p_red%omega(jjfreq).lt. &
                            omega(ibfreq+1)) then

                      ! We found it in the index of the lower
                      lifreq = ibfreq

                      ! Inverse of the distance
                      ! between the two outputs
                      dx = (p_red%omega(jjfreq) - &
                            omega(lifreq))/ &
                           (omega(lifreq+1) - omega(lifreq))

                      ! Inclinations
                      do ith1=1,Geom%nTh

                        ! Output index
                        kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                        ! The first index is the lower
                        y0 = Stk(:,lifreq,1,ith1)

                        ! The second index is the upper
                        dy = Stk(:,lifreq+1,1,ith1) - y0

                        ! Interpolate
                        Stkin(kkfreq,:) = dx*dy + y0

                      end do

                      exit

                    end if ! Check output frequency

                  end do ! Run output frequencies

                end if ! Check if out of limits

              end do ! Run input frequencies

              ! Update frequency index
              jjfreq0 = jjfreq0 + p_mfreq

              ! Update Stokes index
              kkfreq0 = kkfreq0 + p_mfreq*Geom%nTh

            end do ! Output frequencies
          end do ! Output frequency ranges

        end if ! Not dynamic

      ! Non-axially symmetric
      else

        ! If dynamic
        if (lvel) then

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! For each input direction
              do ith1=1,Geom%nth
                do iph1=1,Geom%nph

                  ! Get director cosines
                  cost = Geom%V_mu(ith1)
                  sint = sqrt(1d0 - cost*cost)
                  cosc = Geom%v_mux(iph1)
                  sinc = Geom%v_muy(iph1)*sqrt(1d0 - cosc*cosc)

                  ! Calculate Doppler shift factor
                  vfac1 = 1d0 - vx*sint*cosc - &
                                vy*sint*sinc - &
                                vz*cost

                  ! We will be using the inverse
                  vfac1 = 1d0/vfac1

                  ! Reset the search frequency
                  lifreq = Red%ggf0

                  ! For each input frequency
                  do jfreq=1,p_mfreq

                    ! Advance indexes
                    jjfreq = jjfreq0 + jfreq
                    kkfreq = kkfreq0 + jfreq

                    ! If out of range, take the value at the
                    ! boundary
                    if (p_red%omega(jjfreq)*vfac1.le. &
                        omega(Red%ggf0)+TINYO) then

                      ! We are still looking in the first one
                      lifreq = Red%ggf0

                      ! The index to take is 1
                      Stkin(kkfreq,:) = Stk(:,Red%ggf0,iph1,ith1)

                    ! If out of range, take the value at the
                    ! boundary
                    else if (p_red%omega(jjfreq)*vfac1.ge. &
                             (omega(Red%ggf1) - TINYO)) then

                      ! We are in the last frequency
                      lifreq = Red%ggf1

                      ! The index to take is nfreq
                      Stkin(kkfreq,:) = Stk(:,Red%ggf1,iph1,ith1)

                    ! If within the boundaries
                    else

                      ! Search between the last found frequency and
                      ! all but the boundary
                      do ibfreq=lifreq,nfreq-1

                        ! If this exact frequency is in output
                        if (abs(p_red%omega(jjfreq)*vfac1 - &
                                omega(ibfreq)).lt.TINYO) then

                          ! We are in the found frequency
                          lifreq = ibfreq

                          ! This frequency gives us the value
                          Stkin(kkfreq,:) = Stk(:,lifreq,iph1,ith1)

                          exit

                        ! If the input is between this output and
                        ! the next
                        else if(p_red%omega(jjfreq)*vfac1.ge. &
                                omega(ibfreq).and. &
                                p_red%omega(jjfreq)*vfac1.lt. &
                                omega(ibfreq+1)) then

                          ! We found it in the index of the lower
                          lifreq = ibfreq

                          ! The first index is the lower
                          y0 = Stk(:,lifreq,iph1,ith1)

                          ! The second index is the upper
                          dy = Stk(:,lifreq+1,iph1,ith1) - y0

                          ! Inverse of the distance
                          ! between the two outputs
                          dx = (p_red%omega(jjfreq)*vfac1 - &
                                omega(lifreq))/ &
                               (omega(lifreq+1) - omega(lifreq))

                          ! Interpolate
                          Stkin(kkfreq,:) = dx*dy + y0

                          exit

                        end if ! Check output frequency

                      end do ! Run output frequencies

                    end if ! Check if out of limits

                  end do ! Run input frequencies

                  ! Update index
                  kkfreq0 = kkfreq0 + p_mfreq

                end do ! Input azimuth
              end do ! Input polar

              ! Update index
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! Output frequencies
          end do ! Output frequency ranges

        ! If static
        else

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! Reset the search frequency
              lifreq = Red%ggf0

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq

                ! If out of range, take the value at the
                ! boundary
                if (p_red%omega(jjfreq).le. &
                    omega(Red%ggf0)+TINYO) then

                  ! We are still looking in the first one
                  lifreq = Red%ggf0

                  ! Block counter
                  nblock = -1

                  ! Directions
                  do ith1=1,Geom%nTh
                    do iph1=1,Geom%nPh

                      ! Add block
                      nblock = nblock + 1

                      ! Output index
                      kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                      ! The index to take is 1
                      Stkin(kkfreq,:) = Stk(:,Red%ggf0,iph1,ith1)

                    end do
                  end do

                ! If out of range, take the value at the
                ! boundary
                else if (p_red%omega(jjfreq).ge. &
                         (omega(Red%ggf1) - TINYO)) then

                  ! We are in the last frequency
                  lifreq = Red%ggf1

                  ! Block counter
                  nblock = -1

                  ! Directions
                  do ith1=1,Geom%nTh
                    do iph1=1,Geom%nPh

                      ! Add block
                      nblock = nblock + 1

                      ! Output index
                      kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                      ! The index to take is nfreq
                      Stkin(kkfreq,:) = Stk(:,Red%ggf1,iph1,ith1)

                    end do
                  end do

                ! If within the boundaries
                else

                  ! Search between the last found frequency and
                  ! all but the boundary
                  do ibfreq=lifreq,nfreq-1

                    ! If this exact frequency is in output
                    if (abs(p_red%omega(jjfreq)- &
                            omega(ibfreq)).lt.TINYO) then

                      ! We are in the found frequency
                      lifreq = ibfreq

                      ! Block counter
                      nblock = -1

                      ! Directions
                      do ith1=1,Geom%nTh
                        do iph1=1,Geom%nPh

                          ! Add block
                          nblock = nblock + 1

                          ! Output index
                          kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                          ! This frequency gives us the value
                          Stkin(kkfreq,:) = Stk(:,lifreq,iph1,ith1)

                        end do
                      end do

                      exit

                    ! If the input is between this output and
                    ! the next
                    else if(p_red%omega(jjfreq).ge. &
                            omega(ibfreq).and. &
                            p_red%omega(jjfreq).lt. &
                            omega(ibfreq+1)) then

                      ! We found it in the index of the lower
                      lifreq = ibfreq

                      ! Inverse of the distance
                      ! between the two outputs
                      dx = (p_red%omega(jjfreq) - &
                            omega(lifreq))/ &
                           (omega(lifreq+1) - omega(lifreq))

                      ! Block counter
                      nblock = -1

                      ! Directions
                      do ith1=1,Geom%nTh
                        do iph1=1,Geom%nPh

                          ! Add block
                          nblock = nblock + 1

                          ! The first index is the lower
                          y0 = Stk(:,lifreq,iph1,ith1)

                          ! The second index is the upper
                          dy = Stk(:,lifreq+1,iph1,ith1) - y0

                          ! Output index
                          kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                          ! Interpolate
                          Stkin(kkfreq,:) = dx*dy + y0

                        end do
                      end do

                      exit

                    end if ! Check output frequency

                  end do ! Run output frequencies

                end if ! Check if out of limits

              end do ! Run input frequencies

              ! Update frequency index
              jjfreq0 = jjfreq0 + p_mfreq

              ! Update Stokes index
              kkfreq0 = kkfreq0 + p_mfreq*Geom%nTh*Geom%nPh

            end do ! Output frequencies
          end do ! Output frequency ranges

        end if ! Static
      end if ! Non-axial

      ! Nullify
      nullify(p_mfreq)

      end subroutine getStkin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the frequency dependent JKQ into the requested
      !! frequency\n
      !!  omega(double(:)): Frequency array\n
      !!  JKQ(dcomplex(:)): Frequency dependent radiation field
      !!                    tensors\n
      !!    kfreq(integer): Frequency index of the output frequency
      !!                    associated to the requested input
      !!                    frequency (shifted to limited vector)\n
      !!      jf0(integer): Lower limit to search in Stokes
      !!                    interpolation\n
      !!      jf1(integer): Upper limit to search in Stokes
      !!                    interpolation\n
      !!         x(double): Frequency to interpolate into
      function getJKQinnu(omega,JKQ,kfreq,jf0,jf1,x)

      ! I/O

      integer, intent(in):: kfreq,jf0,jf1
      double precision, intent(in):: x
      double precision, dimension(:), intent(in):: omega
      complex(kind=8), dimension(:), intent(in):: JKQ

      complex(kind=8):: getJKQinnu

      ! Local

      integer:: jfreq,ifreq,if0,if1

      double precision:: dxs

      complex(kind=8):: dys


      ! Copy line limits
      if0 = jf0
      if1 = jf1

      ! Fix line limits
      if (if0.lt.1) if0 = 1
      if (if1.gt.size(omega)) if1 = size(omega)

      !
      ! Check limits
      !

      ! If below lower limit
      if (kfreq.lt.if0) then

        ! Set to lower limit
        ifreq = if0

      ! If above upper limit
      else if (kfreq.gt.if1) then

        ! Set to upper limit
        ifreq = if1

      ! Otherwise
      else

        ! Copy shifted rest index
        ifreq = kfreq

      end if ! Below/above/in limits

      ! Initialize as "equal"
      getJKQinnu = JKQ(ifreq)

      ! If omegai > omega(ifreq)
      if (x.gt.omega(ifreq)) then

        ! If out of right boundary
        if (x.ge.omega(if1)-TINYO) then

          ! Get boundary value
          getJKQinnu = JKQ(if1)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,if1-1

            ! If this exact frequency is in output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              ! Get exact value
              getJKQinnu = JKQ(jfreq)
              return

            ! If the input is between this output and the next
            else if (x.ge.omega(jfreq).and. &
                     x.lt.omega(jfreq+1)) then

              ! Linear interpolation
              dys = JKQ(jfreq+1) - JKQ(jfreq)
              dxs = x - omega(jfreq)
              getJKQinnu = dxs*dys/(omega(jfreq+1) - omega(jfreq)) + &
                           JKQ(jfreq)
              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries

      ! If omegai < omegao
      else if (x.lt.omega(ifreq)) then

        ! If out of left boundary
        if (x.le.omega(if0)+TINYO) then

          ! Get boundary value
          getJKQinnu = JKQ(if0)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,if0+1,-1

            ! If this exact frequency is in output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              ! Get exact value
              getJKQinnu = JKQ(jfreq)
              return

            ! If the input is between this output and the next
            else if (x.ge.omega(jfreq-1).and. &
                     x.lt.omega(jfreq)) then

              ! Linear interpolation
              dys = JKQ(jfreq) - JKQ(jfreq-1)
              dxs = x - omega(jfreq-1)
              getJKQinnu = dxs*dys/(omega(jfreq) - omega(jfreq-1)) + &
                         JKQ(jfreq-1)
              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries
      end if ! omegai > omegao

      end function getJKQinnu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the frequency dependent JKQ into the input
      !! frequency axis\n
      !!    p_red(Redc_class): Structure with redistribution input
      !!                       frequency data for a given input
      !!                       transition\n
      !!      Fed(Reda_class): Structure with redistribution output
      !!                       frequency data\n
      !!      Red(Redb_class): Structure with redistribution input
      !!                       frequency data\n
      !!        Mif0(integer): First frequency for this CPU\n
      !!        Mif1(integer): Last frequency for this CPU\n
      !!      nmfreq(integer): Size of frequency space\n
      !!     omega(double(:)): Frequency array\n
      !!  Jin(dcomplx(:,:,:)): Interpolated frequency dependent JKQ
      !!                       tensors\n
      !!  JKQ(dcomplx(:,:,:)): Original frequency dependent JKQ
      !!                       tensors
      subroutine getJKQin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                          Jin,JKQ)

      ! I/O

      type(Redc_class), intent(in), pointer:: p_red
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      integer, intent(in):: nmfreq,Mif0,Mif1
      double precision, dimension(:), intent(in):: omega
      complex(kind=8), dimension(:,:,:), &
                       allocatable, intent(out):: Jin
      complex(kind=8), dimension(-2:2,0:2,Red%ggf0:Red%ggf1), &
                       intent(in):: JKQ

      ! Local

      integer:: lifreq,ibfreq,iran,ifreq,iifreq,jfreq,jjfreq0,jjfreq

      double precision:: dx

      complex(kind=8), dimension(0:2,0:2):: y0, dy

      ! Pointer

      integer, pointer:: p_mfreq


      ! Nullify
      nullify(p_mfreq)

      ! If no size, return
      if (nmfreq.lt.1) return

      ! Allocate J
      allocate(Jin(nmfreq,0:2,0:2))

      ! Initialize index
      jjfreq0 = 0

      ! For each output frequency
      iifreq = 0
      do iran=1,Fed%nran
        do ifreq=Fed%if0(iran),Fed%if1(iran)

          ! Advance index
          iifreq = iifreq + 1

          ! MPI
          if (iifreq.gt.Mif1) exit
          if (iifreq.lt.Mif0) cycle

          ! Input frequency number
          p_mfreq => p_red%mfreq(iifreq)

          ! Reset the search frequency
          lifreq = Red%ggf0

          ! For each input frequency
          do jfreq=1,p_mfreq

            ! Advance indexes
            jjfreq = jjfreq0 + jfreq

            ! If out of range, take the value at the boundary
            if (p_red%omega(jjfreq).le.omega(Red%ggf0)+TINYO) then

              ! We are still looking in the first one
              lifreq = Red%ggf0

              ! The index to take is 1
              Jin(jjfreq,0:2,0:2) = JKQ(0:2,0:2,Red%ggf0)

            ! If out of range, take the value at the boundary
            else if (p_red%omega(jjfreq).ge. &
                     (omega(Red%ggf1) - TINYO)) then

              ! We are in the last frequency
              lifreq = Red%ggf1

              ! The index to take is nfreq
              Jin(jjfreq,0:2,0:2) = JKQ(0:2,0:2,Red%ggf1)

            ! If within the boundaries
            else

              ! Search between the last found frequency and all but
              ! the boundary
              do ibfreq=lifreq,Red%ggf1-1

                ! If this exact frequency is in output
                if (abs(p_red%omega(jjfreq) - &
                        omega(ibfreq)).lt.TINYO) then

                  ! We are in the found frequency
                  lifreq = ibfreq

                  ! This frequency gives us the value
                  Jin(jjfreq,0:2,0:2) = JKQ(0:2,0:2,lifreq)
                  exit

                ! If the input is between this output and the next
                else if(p_red%omega(jjfreq).ge. &
                        omega(ibfreq).and. &
                        p_red%omega(jjfreq).lt. &
                        omega(ibfreq+1)) then

                  ! We found it in the index of the lower
                  lifreq = ibfreq

                  ! Inverse of the distance between the two outputs
                  dx = (p_red%omega(jjfreq) - omega(lifreq))/ &
                       (omega(lifreq+1) - omega(lifreq))

                  ! The first index is the lower
                  y0 = JKQ(0:2,0:2,ibfreq)

                  ! Difference with next
                  dy = JKQ(0:2,0:2,ibfreq+1) - y0

                  ! Interpolate
                  Jin(jjfreq,0:2,0:2) = dy*dx + y0
                  exit

                end if ! Check output frequency

              end do ! Run output frequencies

            end if ! Check if out of limits

          end do ! Run input frequencies

          ! Update index in general
          jjfreq0 = jjfreq0 + p_mfreq

        end do ! Output frequencies
      end do ! Output frequency ranges

      ! Free pointers
      nullify(p_mfreq)

      end subroutine getJKQin

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeffaux_mod
