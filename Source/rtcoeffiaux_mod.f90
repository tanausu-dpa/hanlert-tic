      !> Absorptivities and emissivities
      module rtcoeffiaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Contributors:
!     John Dennis (NCAR)
!  Start:
!     04/20/2017
!  Last version:
!     08/07/2023 V3.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.4 - Added rt1ordI, absorbILTE, emissILTE,
!                             and rt1ordILTE (TdPA)
!
!     03/23/2023:    V3.0.3 - Bugfix: Ensured the OpenMP version
!                             compiles after some years of changes,
!                             albeit did not test it works (TdPA)
!
!     02/14/2023:    V3.0.2 - AV y axial now need to be AVI and
!                             axiali instead (TdPA)
!
!     10/26/2022:    V3.0.1 - Changed the storage structure of the
!                             rdip variable (TdPA)
!                           - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case Atmo%v has
!                             changed to Atmo%vx,%vy, and %vz (TdPA)
!
!     06/21/2022:    V2.1.0 - Added coherent scattering in the
!                             observers frame (TdPA)
!                             NOTE: Limited testing (AA, static, and
!                             non-magnetic)
!                           - Added function getJin to interpolate J
!                             to a given frequency (TdPA)
!
!     03/22/2021:    V2.0.1 - One of the changes in the angle-averaged
!                             branch in emissI2ord led to a
!                             measurable decrease of performance. I
!                             assume it was related to cache misses,
!                             so I came back to the legacy structure
!                             of that branch (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Deep changes to include OpenMP (TdPA)
!
!     10/26/2020:   V1.13.1 - Now index1, index2, and dx are pointers
!                             to avoid copying (TdPA)
!
!     09/11/2020:   V1.13.0 - Completely changed the Frec and Red
!                             structures (see V1.15.0 in
!                             omegabuild_mod). Applied the changes to
!                             correctly use them (TdPA)
!                           - Compute interpolation data when it is
!                             not stored in RAM (TdPA)
!                           - Check for underflows before storing
!                             the redistribution function in single
!                             precision (TdPA)
!                           - Instead of having a block of code for
!                             storing or not storing, get a pointer
!                             (or a copy is storing) to the
!                             redistribution data (TdPA)
!
!     07/31/2020:   V1.12.2 - Bugfix: There was a pointer directed
!                             in the wrong place, just above its
!                             corresponding loop (TdPA)
!                           - Removed unused code (TdPA)
!
!     07/02/2020:   V1.12.1 - In emiss2ordI, it is necessary to
!                             specify the limits of the CRD profile
!                             when in RAM in case it does not coincide
!                             with the range of indexes for PRD
!                             calculations (TdPA)
!
!     06/08/2020:   V1.12.0 - Split the interpolation stored indexes
!                             into two vectors to speed-up the
!                             emiss2ord routine (JD)
!                           - Added the same extra pointers than in
!                             emiss2ord(NB) in emiss2ordI (TdPA)
!
!     12/10/2019:   V1.11.2 - emissI2ord only received the JKQC range
!                             that is needed (TdPA)
!                           - Experimental changes while trying to
!                             improve the interpolation computational
!                             cost in emissI2ord (TdPA)
!                           - Now we distinguish between
!                             angle-averaged and angle-dependent
!                             Warr2 variables (TdPA)
!
!     11/19/2019:   V1.11.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - The redistribution is now normalized
!                             when computed and not when it is
!                             integrated (TdPA)
!
!     11/13/2019:   V1.11.0 - Added support for voigt profiles read
!                             from a file (TdPA)
!                           - Now routines use the normalized
!                             populations instead of the density
!                             matrix elements (TdPA)
!
!     10/16/2019:   V1.10.0 - The cosine and sine of the scattering
!                             angle for the redistribution function
!                             are computed at the beginning, and
!                             passed directly to WfuncI. This saves a
!                             lot of calls to cos() and sin() (TdPA)
!
!     10/03/2019:    V1.9.1 - trano is now indexed in Red_class (TdPA)
!
!     09/13/2019:    V1.9.0 - Added branch in emissI2ord for angle-
!                             average redistribution for the dynamic
!                             case, following the strategy of angle-
!                             average in the comoving frame (TdPA)
!
!     06/04/2019:    V1.8.0 - Simplifications due to the change of the
!                             iterative scheme (TdPA)
!
!     05/31/2019:    V1.7.3 - Added the emerging variable to the
!                             second order emissivity to distinguish
!                             the angles that are being used to
!                             compute the scattering angle (TdPA)
!                           - Added an underflow check for the
!                             redistribution profiles before type
!                             casting them into single precision in
!                             its structure (TdPA)
!
!     05/08/2019:    V1.7.2 - Introduced a check before storing
!                             redistribution profiles because the type
!                             casting from double to single could
!                             end up in NaN values (TdPA)
!                           - The interpolator getStkinI initializes
!                             as equal instead of waiting to the last
!                             momment (TdPA)
!
!     03/18/2019:    V1.7.1 - Removed inverse exponential from
!                             photoepsIS (TdPA)
!
!     02/20/2019:    V1.7.0 - New verbosity (TdPA)
!                           - Using specific TINY variables (TdPA)
!                           - Using diexp function (TdPA)
!
!     02/11/2019:    V1.6.3 - Bugfix: The interpolation for the
!                             dynamic case was wrong, moved part of
!                             the interpolation to omegabuild (TdPA)
!
!     02/08/2019:    V1.6.2 - Bugfix: Introduced a check for when the
!                             redistribution function integrates to 0
!                             to avoid 0/0 (TdPA)
!
!     09/04/2018:    V1.6.1 - Passing iexu to photoepsIS (TdPA)
!                           - Avoiding to divide by 0 with the
!                             exponential in photoepsI (TdPA)
!
!     08/06/2018:    V1.6.0 - Added the possibility to get Voigt
!                             profiles from RAM (TdPA)
!                           - Changed the argument Norm in the b-b
!                             subroutines (TdPA)
!                           - Added photoepsIS subroutine, that takes
!                             pre-computed frequency quantities from
!                             the RAM (TdPA)
!                           - Made changes to emiss2ord that allows
!                             for partial storage of Wfunc2 (TdPA)
!
!     05/16/2018:    V1.5.0 - Added proper angle-dependent case (TdPA)
!                           - Split the angle-dependent part in axial
!                             and non-axial symmetryc (TdPA)
!
!     12/05/2017:    V1.4.0 - Can neglect Raman scattering (TdPA)
!
!     10/23/2017:    V1.3.0 - Warr2 is stored into a single precision
!                             variable (TdPA)
!
!     09/08/2017:    V1.2.0 - Introduced pointers in emiss2ord to
!                             reduce verbosity and to improve the
!                             performance when compiling without
!                             optimization flags, important for the
!                             debugging mode (TdPA)
!                           - Big bug in angle dependent version,
!                             there was a jdir where there should be
!                             a jtran (TdPA)
!
!     07/21/2017:    V1.1.0 - omega and Wfreq are one step higher in
!                             Frec (TdPA)
!
!     07/19/2017:    V1.0.6 - Using the range structure also in the
!                             variable Warr2 (TdPA)
!
!     06/16/2017:    V1.0.5 - Changed RAM to IRAM (TdPA)
!                           - Changes in emiss2ord towards
!                             optimization (TdPA)
!
!     06/13/2017:    V1.0.4 - Avoid copying Warr2 (TdPA)
!
!     06/12/2017:    V1.0.3 - The limits for the transition are
!                             passed as arguments (TdPA)
!                           - Warr2 can be stored (TdPA)
!
!     06/08/2017:    V1.0.2 - Applied last change from last version
!                             in first order too (TdPA)
!                           - Changed second order emissivity to
!                             calculate the ALI part (TdPA)
!
!     05/12/2017:    V1.0.1 - Bugfix: WfuncI does not need the
!                             energies, because they are already
!                             subtracted in the call to the routine,
!                             there were subtracted twice (TdPA)
!                           - The input frequencies now go with FS
!                             transitions (TdPA)
!                           - Some extra variables to avoid repeating
!                             products when calling voigt in emiss2
!                             (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!  absorbI:
!    This subroutine calculates the absorption coefficients for
!    b-b transitions.
!    Units are cm^-1 (true absorption coefficient)
!    Need be multiplied by the actual atomic density
!
!  emissI:
!    This subroutine calculates the emission Stokes vector
!    The output value is given in number of photons per unit interval
!    of time (s) and normalized frequency (in units of Doppler width),
!    emitted by a unit volume of gas (cm^-3) of unit atomic density,
!    within one steradian.
!    In order to compute photometric values of the intensity, the
!    output needs be multiplied by the actual atomic density.
!
!  rt1ordI:
!    This subroutine calculates the absorption and emission
!    coefficients for b-b transitions
!    Units (abs) are cm^-1 (true absorption coefficient)
!    The output value (emiss) is given in number of photons per unit
!    interval of time (s) and normalized frequency (in units of
!    Doppler width), emitted by a unit volume of gas (cm^-3) of unit
!    atomic density, within one steradian.
!    Need be multiplied by the actual atomic density
!
!  emissI2ord:
!    This subroutine calculates the coherent scattering coefficients
!
!  absorbILTE:
!    This subroutine calculates the absorption coefficients for
!    an LTE line.
!    Units are cm^-1 (true absorption coefficient)
!    Need be multiplied by the actual atomic density
!
!  emissILTE:
!    This subroutine calculates the emissivity of an LTE line
!    The output value is given in number of photons per unit interval
!    of time (s) and normalized frequency (in units of Doppler width),
!    emitted by a unit volume of gas (cm^-3) of unit atomic density,
!    within one steradian.
!    In order to compute photometric values of the intensity, the
!    output needs be multiplied by the actual atomic density.
!
!  rt1ordILTE:
!    This subroutine calculates the absorption and emission
!    coefficients for an LTE line.
!    Units (abs) are cm^-1 (true absorption coefficient)
!    The output (emis) is given in number of photons per unit interval
!    of time (s) and normalized frequency (in units of Doppler width),
!    emitted by a unit volume of gas (cm^-3) of unit atomic density,
!    within one steradian.
!    Need be multiplied by the actual atomic density
!
!  photoabsI:
!    This subroutine calculates the absorption coefficients for b-f
!    transitions
!
!  photoepsI:
!    This subroutine calculates the emission coefficients for b-f
!    transitions
!
!  photoepsIS:
!    This subroutine calculates the emission coefficients for b-f
!    transitions using some precomputed quantities.
!
!  getStkinI:
!    This function interpolates Stokes I for the forward 2-level
!    scattering case
!
!  getJin:
!    This function interpolates J for the coherent scattering in the
!    observers frame case
!
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use math_mod
      use omp_mod
      use parameters_mod , only : IPI, c, c2, convF, cSaha, kb, &
                                  fktoJ, cZero, bigexp, IPI41, &
                                  IPI42, IPI2, IPI4, TINYO, vrfrac
      use profile_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!              iJu(integer): Upper level of transition\n
      !!              iJl(integer): Lower level of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!                pE(dfloat): Unit transformation factor\n
      !!          aprof(dfloat(:)): Absorption profile\n
      !!            eta(dfloat(:)): Intensity absorptivity
      subroutine absorbI(Atom,omega,itran,fitran,itermu,iterml, &
                         iJu,iJl,iz,if0,if1,Norma,Dw,vfac,pE, &
                         aprof,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran,fitran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw, pE, vfac
      double precision, dimension(:), intent(in):: omega, aprof
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq, i

      double precision:: eu,el,au,al,aul,rho,at,feta
      double precision:: Dfreqw,vfacw,prof


      !
      ! Get population factor
      !

      ! Level index
      i = Atom%irho(iterml)%irho_ij(iJl)

      ! Population
      rho = dble(Atom%popu(i,iz))

      ! Absorptibity factor
      feta = rho*1d3*IPI41*Atom%fst(itran)%Blu(iJl,iJu)*pE/Dw

      ! If reading a file
      if (vifil) then

!$omp parallel workshare default(none) &
!$omp shared(feta,aprof,eta)
        eta = feta*aprof
!$omp end parallel workshare

      ! If stored in RAM
      else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(fitran,feta,Norma,eta)
        eta = feta*Norma%prof(fitran,1,1,1)%p
!$omp end parallel workshare

      ! Not stored
      else

        ! Level quantities

        ! Damping parameter
        al = Atom%damp(iterml,iz)/Dw

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Damping parameter
        au = Atom%damp(itermu,iz)/Dw

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)/Dw

        ! Intermediate quantities
        at = au+al+aul
        Dfreqw = (eu - el)/Dw
        vfacw = vfac/Dw
        feta = feta*Norma%Norm(fitran,1,1,1)

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feta,Dfreqw,omega,vfacw,at,eta)
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          eta(ifreq) = feta*prof

        end do ! frequencies
!$omp end parallel do

      end if ! Type of profile calculation

      return

      end subroutine absorbI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emission coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!              iJu(integer): Upper level of transition\n
      !!              iJl(integer): Lower level of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!          aprof(dfloat(:)): Emission profile\n
      !!            eps(dfloat(:)): Intensity emissivity\n
      !!              rhou(dfloat): Factor for Lambda operator
      subroutine emissI(Atom,omega,itran,fitran,itermu,iterml, &
                        iJu,iJl,iz,if0,if1,Norma,Dw,vfac, &
                        aprof,eps,rhou)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran,fitran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw, vfac
      double precision, intent(out):: rhou
      double precision, dimension(:), intent(in):: omega, aprof
      double precision, dimension(if0:if1), intent(out):: eps

      ! Local

      integer:: ifreq,i,iR

      double precision:: el,eu,al,au,aul,rho,at
      double precision:: Dfreqw,vfacw,prof,feps


      !
      ! Get population factor
      !

      ! Level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Population
      rho = Atom%popu(i,iz)

      ! Emissivity factor
      feps = rho*1d3*IPI41*Atom%fst(itran)%Aul(iJu,iJl)/Dw

      ! If reading a file
      if (vifil) then

!$omp parallel workshare default(none) &
!$omp shared(feps,aprof,eps)
        eps = feps*aprof
!$omp end parallel workshare

      ! If stored in RAM
      else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(fitran,feps,Norma,eps)
        eps = feps*Norma%prof(fitran,1,1,1)%p
!$omp end parallel workshare

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)/Dw

        ! Level quantities

        ! Damping parameter
        au = Atom%damp(itermu,iz)/Dw

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Damping parameter
        al = Atom%damp(iterml,iz)/Dw

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Intermediate quantities
        at = au+al+aul
        Dfreqw = (eu - el)/Dw
        vfacw = vfac/Dw
        feps = feps*Norma%Norm(fitran,1,1,1)

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feps,Dfreqw,omega,vfacw,at,eps)
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          eps(ifreq) = feps*prof

        end do ! frequencies
!$omp end parallel do

      end if ! Type of profile calculation

      return

      end subroutine emissI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!              iJu(integer): Upper level of transition\n
      !!              iJl(integer): Lower level of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!                pE(dfloat): Unit transformation factor\n
      !!          aprof(dfloat(:)): Absorption profile\n
      !!            eta(dfloat(:)): Intensity absorptivity
      !!            eps(dfloat(:)): Intensity emissivity\n
      !!              rhou(dfloat): Factor for Lambda operator
      subroutine rt1ordI(Atom,omega,itran,fitran,itermu,iterml, &
                         iJu,iJl,iz,if0,if1,Norma,Dw,vfac,pE, &
                         aprof,eta,eps,rhou)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran,fitran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw, pE, vfac
      double precision, intent(out):: rhou
      double precision, dimension(:), intent(in):: omega, aprof
      double precision, dimension(if0:if1), intent(out):: eta, eps

      ! Local

      integer:: ifreq, i, iR

      double precision:: eu,el,au,al,aul,rho,at,feta,feps
      double precision:: Dfreqw,vfacw,prof

      !
      ! Get population factor
      !

      !
      ! Lower

      ! Level index
      i = Atom%irho(iterml)%irho_ij(iJl)

      ! Population
      rho = dble(Atom%popu(i,iz))

      ! Absorptibity factor
      feta = rho*1d3*IPI41*Atom%fst(itran)%Blu(iJl,iJu)*pE/Dw

      !
      ! Upper

      ! Level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Population
      rho = Atom%popu(i,iz)

      ! Emissivity factor
      feps = rho*1d3*IPI41*Atom%fst(itran)%Aul(iJu,iJl)/Dw

      ! If reading a file
      if (vifil) then

!$omp parallel workshare default(none) &
!$omp shared(feps,feta,aprof,eps,eta)
        eps = feps*aprof
        eta = feta*aprof
!$omp end parallel workshare

      ! If stored in RAM
      else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(fitran,feps,feta,Norma,eps,eta)
        eps = feps*Norma%prof(fitran,1,1,1)%p
        eta = feta*Norma%prof(fitran,1,1,1)%p
!$omp end parallel workshare

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)/Dw

        ! Level quantities

        ! Damping parameter
        au = Atom%damp(itermu,iz)/Dw

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Damping parameter
        al = Atom%damp(iterml,iz)/Dw

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Intermediate quantities
        at = au+al+aul
        Dfreqw = (eu - el)/Dw
        vfacw = vfac/Dw
        feps = feps*Norma%Norm(fitran,1,1,1)

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feps,feta,Dfreqw,omega,vfacw,at,eps,eta)
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          eps(ifreq) = feps*prof
          eta(ifreq) = feta*prof

        end do ! frequencies
!$omp end parallel do

      end if ! Type of profile calculation

      return

      end subroutine rt1ordI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the second order emission coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!         emerging(logical): Indicates if emergence solution\n
      !!                vx(dfloat): Velocity vector along X\n
      !!                vy(dfloat): Velocity vector along Y\n
      !!                vz(dfloat): Velocity vector along Z\n
      !!          omega(dfloat(:)): Frequency array\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!            stype(integer): Type of scattering (geometry
      !!                            wise)\n
      !!    Fin(Frequencyc2_class): Structure with the input frequency
      !!                            information\n
      !!             nbth(integer): Maximum index in polar angles for
      !!                            the input directions\n
      !!             nbph(integer): Maximum index in azimuth angles
      !!                            for the input directions\n
      !!              nfs(integer): Number of forward scattering
      !!                            directions\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!            jtran(integer): Output transition index term
      !!                            wise\n
      !!           fjtran(integer): Output transition index level
      !!                            wise\n
      !!           itermu(integer): Upper term of output transition\n
      !!           itermf(integer): Lower term of output transition\n
      !!              iJu(integer): Upper level of output transition\n
      !!              iJf(integer): Lower level of output transition\n
      !!               iz(integer): Height index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              ith(integer): Output direction polar index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!               DwT(dfloat): Thermal part of Doppler width\n
      !!                Dw(dfloat): Doppler width output transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!               vmi(dfloat): Microturbulent velocity\n
      !!     Stokes(dfloat(:,:,:)): Stokes parameters\n
      !!            JKQ(dfloat(:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!           JKQC(dfloat(:)): Mean intensity with frequency
      !!                            dependence\n
      !!          aprof(dfloat(:)): Absorption profile\n
      !!           eps2(dfloat(:)): Intensity emissivity\n
      !!            rpf(dfloat(:)): Factor for lambda operator
      subroutine emissI2ord(Atom,Geom,emerging,vx,vy,vz,omega, &
                            Red,stype,Fin,nbth,nbph,nfs,Norma, &
                            jtran,fjtran,itermu,itermf,iJu,iJf, &
                            iz,iph,ith,if0,if1,DwT,Dw,vfac,vmi, &
                            Stokes,JKQ,JKQC,aprof,eps2,rpf)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Frequencyc2_class), intent(inout):: Fin
      type(Nindex_class), intent(in):: Norma
      type(Redc2_class), intent(inout):: Red
      logical, intent(in):: emerging
      integer, intent(in):: jtran,fjtran,itermu,itermf,iJu,iJf
      integer, intent(in):: iz,ith,iph,if0,if1
      integer, intent(in):: nfs,nbth,nbph
      integer, dimension(:,:), intent(in):: stype
      double precision, intent(in):: DwT,Dw,vfac,vmi,vx,vy,vz
      double precision, dimension(:), intent(in):: omega, aprof
      double precision, dimension(:), intent(in):: JKQ
      double precision, dimension(Fin%ggf0:Fin%ggf1), &
                        intent(in), target:: JKQC
      double precision, dimension(:,:,:), intent(in):: Stokes
      double precision, dimension(if0:if1), intent(out):: eps2,rpf

      ! Local

      logical:: PRDc,LIRAM,shift,tointerp

      integer:: i,j,iR,itran,fitran,iterml,ffjtran,ffitran,ffktran
      integer:: ith1,iph1,ifreq,iifreq,jfreq,iran
#ifdef _OPENMP
      integer:: llfreq,iidir,tid
#endif
      integer:: jjfreq,jjfreq0,kkfreq,kkfreq0,nmfreq,lifreq,ibfreq
      integer:: iJl,iti,ios

      double precision:: auxb,norme2,cost,sint,cosc,sinc
      double precision:: rLl,rLu,rLf,S,rJl,rJu,rJf,vfac1
      double precision:: el,eu,ef,al,au,af,auf,aul,Dw1
      double precision:: StokesM,Jrad,Dfreq,Dfreq1
      double precision:: PRDin,prof,hanleden,rhoc,rhou,dy,y0
      double precision:: vfacw,Dfreqw,atfw,atf,atl,omegao,omegai
      double precision, dimension(Geom%nPh2,Geom%nTh):: ThK
      double precision, dimension(if0:if1):: PRD,CRD
      double precision, dimension(:), allocatable:: StokesMV
      double precision, dimension(:), allocatable, target:: JKQinMV
      double precision, dimension(:), allocatable:: Wcos, Wsin
      double precision, dimension(:), allocatable, target:: Warr2
      double precision, dimension(:,:), allocatable:: TWcos, TWsin
#ifdef _OPENMP
      double precision, dimension(:,:), allocatable:: PRDdir
#endif

      ! Pointers
      type(Frequencyd_class), pointer:: p_frec
      type(Redd_class), pointer:: p_red
      integer, pointer:: p_mfreq
      integer, dimension(:), pointer:: p_index1, p_index2
      double precision, dimension(:), pointer:: p_dx
      double precision, dimension(:), pointer:: p_JKQ
      double precision, dimension(:), pointer:: p_warr2


      ! Routine name
      urou = 'emissI2ord'

      ! Initialize pointers
      nullify(p_frec)
      nullify(p_red)
      nullify(p_mfreq)
      nullify(p_JKQ)
      nullify(p_warr2)
      nullify(p_index1)
      nullify(p_index2)
      nullify(p_dx)


      !
      ! Scattering angles
      !

      ! If angle-averaged redistribution
      if (AVI) then

        !
        ! Compute cosines and sines
        !

        ! Allocate arrays
        allocate(Wcos(Geom%nThAA))
        allocate(Wsin(Geom%nThAA))

        ! For each direction in the integral AA quadrature
        do ith1=1,Geom%nThAA
          Wcos(ith1) = cos(Geom%V_thetaAA(ith1))
          Wsin(ith1) = sin(Geom%V_thetaAA(ith1))
        end do

        !
        ! Construct mean intensity if needed
        !

        ! If dynamics
        if (dyn) then

          ! Allocate mean intensity
          allocate(JKQinMV(Fin%ggf0:Fin%ggf1))

          ! If velocity is below threshold
          shift = (vx*vx + vy*vy + vz*vz)*1d6*c.ge. &
                  vrfrac*DwT
          vfac1 = 1d0

          ! Initialize
          JKQinMV = 0d0

          ! For each frequency, get mean intensity
!$omp parallel do default(none) &
!$omp private(ifreq,ith1,cost,vfac1,omegao,StokesM,sint,iph1) &
!$omp private(cosc,sinc) &
!$omp shared(Fin,Geom,shift,axiali,vx,vy,vz,omega,JKQinMV,Stokes)
          do ifreq=Fin%ggf0,Fin%ggf1

            ! For each polar direction
            do ith1=1,Geom%nTh

              ! Get director cosines
              if (shift) cost = Geom%V_mu(ith1)

              ! If axial symmetric
              if (axiali) then

                ! Calculate Doppler shift factor
                if (shift) then
                  vfac1 = 1d0 + vz*cost
                  omegao = omega(ifreq)*vfac1
                  StokesM = getStkinI(omega, &
                                      Stokes(:,1,ith1), &
                                      ifreq,omegao)
                  JKQinMV(ifreq) = JKQinMV(ifreq) + &
                                   StokesM*Geom%W_mu(ith1)
                else
                  JKQinMV(ifreq) = JKQinMV(ifreq) + &
                                   Stokes(ifreq,1,ith1)* &
                                   Geom%W_mu(ith1)
                end if

              ! Non axial symmetric
              else

                if (shift) sint = sqrt(1d0 - cost*cost)

                ! For each azimuth
                do iph1=1,Geom%nPh

                  ! Calculate Doppler shift factor
                  if (shift) then
                    cosc = Geom%v_mux(iph1)
                    sinc = Geom%v_muy(iph1)* &
                           sqrt(1d0 - cosc*cosc)
                    vfac1 = 1d0 + vx*sint*cosc + &
                                  vy*sint*sinc + &
                                  vz*cost
                    omegao = omega(ifreq)*vfac1
                    StokesM = getStkinI(omega, &
                                        Stokes(:,iph1,ith1), &
                                        ifreq,omegao)
                    JKQinMV(ifreq) = JKQinMV(ifreq) + &
                                     StokesM*Geom%W_mu(ith1)* &
                                     Geom%W_mux2(iph1)
                  else
                    JKQinMV(ifreq) = JKQinMV(ifreq) + &
                                     Stokes(ifreq,iph1,ith1)* &
                                     Geom%W_mu(ith1)* &
                                     Geom%W_mux2(iph1)
                  end if

                end do ! Azimuth

              end if ! Axial symmetry

            end do ! Polar
          enddo ! Frequencies
!$omp end parallel do

        end if ! Dynamics

      ! If angle-dependent redistribution
      else

        ! Allocate arrays for cosines and sines
        allocate(TWcos(Geom%nPh2,Geom%nTh))
        allocate(TWsin(Geom%nPh2,Geom%nTh))

        ! If doing an emerging direction
        if (emerging) then

          ! For each polar quadrature angle
          do ith1=1,Geom%nTh

            ! For each azimuthal quadrature angle
            do iph1=1,Geom%nPh2

              ! Calculate scattering angle between the quadrature
              ! direction and the LOS direction
              ThK(iph1,ith1) = atom2lab(Geom%L_theta(ith), &
                                        Geom%L_phi(iph), &
                                        Geom%V_theta(ith1), &
                                        Geom%V_phi(iph1))

              ! Compute trigonometric functions
              TWcos(iph1,ith1) = cos(ThK(iph1,ith1))
              TWsin(iph1,ith1) = sin(ThK(iph1,ith1))

            end do ! azimuthal nodes
          end do ! polar nodes

        ! If NOT doing an emerging direction
        else

          ! For each polar quadrature angle
          do ith1=1,Geom%nTh

            ! For each azimuthal quadrature angle
            do iph1=1,Geom%nPh2

              ! Calculate scattering angle between the quadrature
              ! directions and the current one
              ThK(iph1,ith1) = atom2lab(Geom%V_theta(ith), &
                                        Geom%V_phi(iph), &
                                        Geom%V_theta(ith1), &
                                        Geom%V_phi(iph1))

              ! Compute trigonometric functions
              TWcos(iph1,ith1) = cos(ThK(iph1,ith1))
              TWsin(iph1,ith1) = sin(ThK(iph1,ith1))

            end do ! azimuthal nodes
          end do ! polar nodes

        end if ! emergence direction

        ! If axial, prepare vector to store Stokes of input
        if (axiali) then

          allocate(StokesMV(Fin%mxfreq))

        end if ! Axial symmetry

#ifdef _OPENMP
        ! Allocate auxiliar variable to deal with OpenMP
        allocate(PRDdir(Geom%nTh*Geom%nPh2,if0:if1))
#endif
      end if ! angle-dependent


      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      atf = au + af + auf

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLf = Atom%rLval(itermf)

      ! Get FS global index of output transition
      ffjtran = Atom%ifst_ij(fjtran,jtran)

      ! Energies
      ef = Atom%FSfreq(iJf,itermf)
      eu = Atom%FSfreq(iJu,itermu)

      ! Angular momentums
      rJf = Atom%rJval(iJf,itermf)
      rJu = Atom%rJval(iJu,itermu)

      ! Transition frequency
      Dfreq = eu - ef

      ! Upper level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! Upper level density
      rhou = Atom%popu(i,iz)

      ! Trano index
      ffktran = Atom%itrano(ffjtran)

      !
      ! Initializations
      !
      eps2 = 0d0


      !
      ! Flat contribution. Implicit branching
      !

      ! If in file
      if (vifil) then

!$omp parallel workshare default(none) &
!$omp shared(CRD,aprof,Dw)
        CRD = aprof*1d5*sqrt(IPI)/Dw
!$omp end parallel workshare

      ! If stored
      else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(CRD,Norma,fjtran,Dw,if0,if1)
        CRD = Norma%prof(fjtran,1,1,1)%p(if0:if1)*1d5*sqrt(IPI)/Dw
!$omp end parallel workshare

      ! If not stored
      else

        Dfreqw = Dfreq/Dw
        atfw = atf/Dw
        vfacw = vfac/Dw
        norme2 = Norma%Norm(fjtran,1,1,1)*1d5*sqrt(IPI)/Dw

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(CRD,if0,if1,norme2,Dfreqw,vfacw,atfw,omega)
        do ifreq=if0,if1

          ! Calculate profile u-f
          call voigtI(Dfreqw - omega(ifreq)*vfacw,atfw,prof)

          ! Flat spectrum contribution
          CRD(ifreq) = prof*norme2

        end do ! frequencies
!$omp end parallel do

      end if ! Storing Voigt


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible lower terms
      do i=1,Atom%nMulti-1

        ! If there is no transition or this term is larger
        ! than the upper term of the output transition, skip
        if(i.ge.itermu.or.Atom%irad(i,itermu).lt.1) cycle

        ! Store the input lower term index
        iterml = i

        ! Get index of input transition
        itran = Atom%irad(iterml,itermu)

        ! Damping parameter lower level input transition
        al = Atom%damp(iterml,iz)
        aul = Atom%ldamp(itran,iz)
        atl = au + al + aul

        ! Angular momentum input lower level
        rLl = Atom%rLval(iterml)

        ! For each lower level
        do iJl=1,Atom%nJ(iterml)

          ! Check if there is a FS transition
          if(Atom%fst(itran)%irad(iJu,iJl).lt.1) cycle

          ! Energy input lower level
          el = Atom%FSfreq(iJl,iterml)

          ! Angular momentum
          rJl = Atom%rJval(iJl,iterml)

          ! Set input transition indexes
          fitran = Atom%fst(itran)%irad(iJu,iJl)
          ffitran = Atom%ifst_ij(fitran,itran)

          ! Get index of input transition in structure
          ios = -1
          do iti=1,Atom%trano(ffktran)%nt
            if (Atom%trano(ffktran)%ind(iti).eq.ffitran) then
              ios = 1
              exit
            end if
          end do
          if (ios.lt.0) cycle

          ! If IRAM, point to the redistribution subblock
          if (IRAM) then

            p_red => Red%trani(iti)

            LIRAM = IRAM.and.p_red%RAM

          ! If not, nothing stored
          else

            LIRAM = .False.

          end if

          ! Point to input transition
          p_frec => Fin%trani(iti)

          ! Predict size of interpolation block
          nmfreq = sum(p_frec%mfreq)

          ! If angle-dependent
          if (.not.AVI) then

            ! If dynamic extra dimensions, if static just
            ! frequencies
            if (dyn) then

              ! For axial problems
              if (axiali) then

                ! Size is just polar
                nmfreq = nmfreq*Geom%nTh

              ! For non-axial problems
              else

                ! Skip backward rayleigh
                nmfreq = nmfreq*(Geom%nTh*Geom%nPh2 - nfs)

              end if ! Axial
            end if ! Dynamic
          end if ! AD

          ! If interpolation data
          if (p_frec%RAM) then

            p_index1 => p_frec%index1
            p_index2 => p_frec%index2
            p_dx => p_frec%dx

          ! No interpolation data
          else if (nmfreq.gt.0) then

            ! Allocate
            allocate(p_index1(nmfreq))
            allocate(p_index2(nmfreq))
            allocate(p_dx(nmfreq))

!$omp parallel default(none) &
!$omp private(iifreq,lifreq,iran,ifreq,jfreq,p_mfreq) &
!$omp private(ith1,iph1,cost,sint,cosc,sinc,vfac1) &
!$omp private(jjfreq,kkfreq,jjfreq0,kkfreq0,ibfreq) &
!$omp shared(ffjtran,ffitran,p_index1,p_index2,p_dx,p_frec) &
!$omp shared(dyn,AVI,axiali,nbth,nbph,vx,vy,vz,stype,omega,Fin,Geom) &
!$omp shared(nfreq)

            ! Initialize index
            jjfreq0 = 0
            kkfreq0 = 0

            ! For each output frequency
            iifreq = 0
            do iran=1,Fin%nran
              do ifreq=Fin%if0(iran),Fin%if1(iran)

                ! Advance index
                iifreq = iifreq + 1

                ! Input frequency number
                p_mfreq => p_frec%mfreq(iifreq)

                ! For each input direction
                do ith1=1,nbth
                  do iph1=1,nbph

                    ! If dynamics and AD
                    if (dyn.and..not.AVI) then

                      ! For axial problems
                      if (axiali) then

                        ! Automatically skip extra azimuths
                        if (iph1.gt.1) cycle

                        ! Get director cosines
                        cost = Geom%V_mu(ith1)

                        ! Calculate Doppler shift factor
                        vfac1 = 1d0 - vz*cost

                        ! We will be using the inverse
                        vfac1 = 1d0/vfac1

                      ! For non-axial problems
                      else

                        ! If angle-dependent, check backward Rayleigh
                        ! scattering
                        if (ffjtran.eq.ffitran.and. &
                            stype(iph1,ith1).lt.0) cycle

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

                      end if ! Axial

                    ! Not dynamic or AV
                    else

                      ! No shift
                      vfac1 = 1d0

                      ! Only one direction
                      if (iph1.gt.1.or.ith1.gt.1) cycle

                    end if ! Dynamics

                    ! Reset the search frequency
                    lifreq = Fin%ggf0

                    ! For each input frequency
!$omp do
                    do jfreq=1,p_mfreq

                      ! Advance indexes
                      jjfreq = jjfreq0 + jfreq
                      kkfreq = kkfreq0 + jfreq

                      ! If out of range, take the value at the
                      ! boundary
                      if (p_frec%omega(jjfreq)*vfac1.le. &
                          omega(1)+TINYO) then

                        ! We are still looking in the first one
                        lifreq = 1

                        ! The index to take is 1
                        p_index1(kkfreq) = 1

                        ! The index to take is 1
                        p_index2(kkfreq) = 1

                        ! We do not need this number
                        p_dx(kkfreq) = 0d0

                      ! If out of range, take the value at the
                      ! boundary
                      else if (p_frec%omega(jjfreq)*vfac1.ge. &
                               (omega(nfreq) - TINYO)) then

                        ! We are in the last frequency
                        lifreq = nfreq

                        ! The index to take is nfreq
                        p_index1(kkfreq) = nfreq

                        ! The index to take is nfreq
                        p_index2(kkfreq) = nfreq

                        ! We do not need this number
                        p_dx(kkfreq) = 0d0

                      ! If within the boundaries
                      else

                        ! Search between the last found frequency and
                        ! all but the boundary
                        do ibfreq=lifreq,nfreq-1

                          ! If this exact frequency is in output
                          if (abs(p_frec%omega(jjfreq)*vfac1 - &
                                  omega(ibfreq)).lt.TINYO) then

                            ! We are in the found frequency
                            lifreq = ibfreq

                            ! This frequency gives us the value
                            p_index1(kkfreq) = lifreq

                            ! This frequency gives us the value
                            p_index2(kkfreq) = lifreq

                            ! We do not need this number
                            p_dx(kkfreq) = 0d0

                            exit

                          ! If the input is between this output and
                          ! the next
                          else if(p_frec%omega(jjfreq)*vfac1.ge. &
                                  omega(ibfreq).and. &
                                  p_frec%omega(jjfreq)*vfac1.lt. &
                                  omega(ibfreq+1)) then

                            ! We found it in the index of the lower
                            lifreq = ibfreq

                            ! The first index is the lower
                            p_index1(kkfreq) = lifreq

                            ! The second index is the upper
                            p_index2(kkfreq) = lifreq+1

                            ! Store the inverse of the distance
                            ! between the two outputs
                            p_dx(kkfreq) = &
                                (p_frec%omega(jjfreq)*vfac1 - &
                                 omega(lifreq))/ &
                                (omega(lifreq+1) - omega(lifreq))

                            exit

                          end if ! Check output frequency

                        end do ! Run output frequencies

                      end if ! Check if out of limits

                    end do ! Run input frequencies
!$omp end do
                    ! Update index in general
                    kkfreq = kkfreq0 + p_mfreq
                    jjfreq = jjfreq0 + p_mfreq

                    ! Actually update index
                    if (.not.AVI.and..not.axiali.and.dyn) &
                      kkfreq0 = kkfreq

                  end do ! Input azimuth

                  ! Update kkfreq
                  if (.not.AVI.and.dyn) kkfreq0 = kkfreq

                end do ! Input polar

                ! Update index
                jjfreq0 = jjfreq
                kkfreq0 = kkfreq

              end do ! Output frequencies
            end do ! Output frequency ranges
!$omp end parallel

          end if ! Interpolation data

          ! Flat spectrum J00
          JRad = JKQ(ffitran)

          ! Input transition frequency
          Dfreq1 = eu - el

          ! Doppler width for the input transition
          Dw1 = Dfreq1*sqrt(DwT*DwT + vmi**2d0)

          ! Hanle factor
          ! TODO ATTENTION TO THIS
          hanleden = 2d0*(au+auf)/Dw

          if (LIRAM) then
            PRDc = p_red%iIPRD
          else
            PRDc = .True.
          end if

          !
          ! Create array of Wfunc
          !

          if (PRDc) then

            ! If angle-averaged
            if(AVI)then

              ! Initialize array
              nmfreq = sum(p_frec%mfreq)
              if (nmfreq.gt.0) then
                if (.not.allocated(Warr2)) then
                  allocate(Warr2(nmfreq))
                else
                  if (size(Warr2).ne.nmfreq) then
                    deallocate(Warr2)
                    allocate(Warr2(nmfreq))
                  end if
                end if
              end if

!$omp parallel default(none) &
!$omp private(jjfreq,iifreq,iran,ifreq,jfreq,omegao,p_mfreq,Norme2) &
!$omp private(ith1,omegai,tid,kkfreq) &
!$omp shared(Fin,omega,vfac,Dfreq,Dfreq1,p_frec,Geom,Dw,Dw1,atl,atf) &
!$omp shared(Wcos,Wsin,Warr2,omp,LIRAM,p_red,nmfreq)

#ifdef _OPENMP
              tid = omp_get_thread_num() + 1
#endif

              ! Inititlize redistribution
!$omp workshare
              Warr2 = 0d0
!$omp end workshare

              ! Initialize frequency index
              jjfreq = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fin%nran
                do ifreq=Fin%if0(iran),Fin%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Point to dimension
                  p_mfreq => p_frec%mfreq(iifreq)

                  ! Skip coherent
                  if (p_mfreq.lt.1) cycle
#ifdef _OPENMP
                  ! If multi-threading
                  if (omp) then

                    ! If out of range below
                    if (iifreq.lt.Fin%oif0(tid)) then
                      jjfreq = jjfreq + p_mfreq
                      cycle
                    end if
                    ! If out of range above
                    if (iifreq.gt.Fin%oif1(tid)) exit
                  end if
#endif
                  ! Get output frequency
                  omegao = omega(ifreq)*vfac - Dfreq

                  ! Initialize norm
                  Norme2 = 0d0

                  ! For each input frequency
                  do kkfreq=jjfreq+1,jjfreq+p_mfreq

                    ! Get input frequency
                    omegai = p_frec%omega(kkfreq) - Dfreq1

                    ! For each direction in the integral AA quadrature
                    do ith1=1,Geom%nThAA

                      ! Add the contribution to the angular integral
                      ! of the redistribution function
                      Warr2(kkfreq) = Warr2(kkfreq) + &
                                                  Geom%W_muAA(ith1)* &
                                WfuncI(omegai,omegao,Dw,Dw1,atl,atf, &
                                       Wcos(ith1),Wsin(ith1),0)*IPI42

                    end do  ! direction nodes

                    ! Store with weight
                    Warr2(kkfreq) = Warr2(kkfreq)* &
                                    p_frec%w_freq(kkfreq)

                    ! Add to norm
                    Norme2 = Norme2 + Warr2(kkfreq)

                  end do ! input frequencies

                  ! Update jjfreq
                  jjfreq = jjfreq + p_mfreq

                  ! Apply normalization
                  if (Norme2.gt.0d0) then
                    Warr2(jjfreq-p_mfreq+1:jjfreq) = &
                                    Warr2(jjfreq-p_mfreq+1:jjfreq)/ &
                                    Norme2
                  end if

                end do ! output frequencies
              end do ! output frequencies ranges
!$omp barrier
!$omp flush (Warr2)

              ! If storing
              if (LIRAM) then
!$omp single
                p_red%iIPRD = .False.
!$omp end single
!$omp do
                do jfreq=1,nmfreq
                  if (Warr2(jfreq).gt.1e-30) then
                    p_red%IWarr2(jfreq) = real(Warr2(jfreq))
                  else
                    p_red%IWarr2(jfreq) = 0e0
                  end if
                end do
!$omp end do
              end if ! If storing
!$omp end parallel

            ! If angle-dependent
            else

              ! Initialize array
              if (ffitran.eq.ffjtran) then
                nmfreq = sum(p_frec%mfreq)*(Geom%nTh*Geom%nPh2-nfs)
              else
                nmfreq = sum(p_frec%mfreq)*Geom%nTh*Geom%nPh2
              end if
              if (nmfreq.gt.0) then
                if (.not.allocated(Warr2)) then
                  allocate(Warr2(nmfreq))
                else
                  if (size(Warr2).ne.nmfreq) then
                    deallocate(Warr2)
                    allocate(Warr2(nmfreq))
                  end if
                end if
              end if

!$omp parallel default(none) &
!$omp private(jjfreq0,jjfreq,kkfreq0,kkfreq,iifreq,iran,ifreq) &
!$omp private(ith1,iph1,llfreq,omegao,omegai,Norme2,p_mfreq,tid) &
!$omp shared(omega,vfac,Dfreq,ffitran,ffjtran,stype,Warr2,Fin,Geom) &
!$omp shared(Dw,Dw1,atl,atf,TWcos,TWsin,p_frec,LIRAM,p_red,nmfreq) &
!$omp shared(Dfreq1,omp)

#ifdef _OPENMP
              tid = omp_get_thread_num() + 1
              llfreq = 0
#endif

              ! Initialize frequency indexes
              jjfreq0 = 0
              kkfreq0 = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fin%nran
                do ifreq=Fin%if0(iran),Fin%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Point to dimension
                  p_mfreq => p_frec%mfreq(iifreq)

                  ! Coherent wing
                  if (p_mfreq.lt.1) cycle

                  ! Get output frequency
                  omegao = omega(ifreq)*vfac - Dfreq

                  ! For each polar direction
                  do ith1=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph1=1,Geom%nPh2

                      ! Check backward scattering two-terms
                      if (ffitran.eq.ffjtran.and. &
                          stype(iph1,ith1).lt.0) cycle
#ifdef _OPENMP
                      ! If multi-threading
                      if (omp) then

                        ! Advance llfreq index
                        llfreq = llfreq + 1

                        ! If out of range below
                        if (llfreq.lt.Fin%oif0(tid)) then
                          kkfreq0 = kkfreq0 + p_mfreq
                          cycle
                        end if
                        ! If out of range above
                        if (llfreq.gt.Fin%oif1(tid)) exit
                      end if
#endif
                      ! Inititlize Norm
                      Norme2 = 0d0

                      ! For each input frequency
                      do jfreq=1,p_mfreq

                        ! Advance indexes
                        jjfreq = jjfreq0 + jfreq
                        kkfreq = kkfreq0 + jfreq

                        ! Get input frequency
                        omegai = p_frec%omega(jjfreq) - Dfreq1

                        ! Calculate redistribution function
                        ! and apply weight
                        Warr2(kkfreq) = WfuncI(omegai,omegao, &
                                               Dw,Dw1,atl,atf, &
                                               TWcos(iph1,ith1), &
                                               TWsin(iph1,ith1), &
                                               stype(iph1,ith1))* &
                                        IPI42*p_frec%w_freq(jjfreq)

                        ! Add to norm
                        Norme2 = Norme2 + Warr2(kkfreq)

                      end do ! mfreq

                      ! Normalize
                      if (Norme2.gt.0d0) &
                        Warr2(kkfreq-p_mfreq+1:kkfreq) = &
                                 Warr2(kkfreq-p_mfreq+1:kkfreq)/Norme2

                      ! Update kkfreq0
                      kkfreq0 = kkfreq0 + p_mfreq

                    end do  ! Azimuthal directions
                  end do  ! Polar directions

                  ! Update jjfreq0
                  jjfreq0 = jjfreq0 + p_mfreq

                end do ! Output frequencies
              end do  ! Output frequency ranges
!$omp barrier
!$omp flush (Warr2)

              ! If storing
              if (LIRAM) then

!$omp single
                p_red%iIPRD = .False.
!$omp end single
!$omp do
                do jfreq=1,nmfreq
                  if (Warr2(jfreq).gt.1e-30) then
                    p_red%Iwarr2(jfreq) = real(Warr2(jfreq))
                  else
                    p_red%Iwarr2(jfreq) = 0e0
                  end if
                end do
!$omp end do
              end if ! If storing
!$omp end parallel
            end if ! AV or AD
          end if ! initialized

          ! Lower input level index
          j = Atom%irho(iterml)%irho_ij(iJl)

          ! Lower input level SEE index
          iR = Atom%irho(iterml)%Jrho(iJl,iJl)%kq(0,0)
          rhoc = sqrt(2d0*rJl+1d0)*dble(Atom%crho(iR,iz))


          !
          ! Integral over input frequencies
          !

          ! If angle-averaged
          if(AVI)then

            ! If dynamic
            if (dyn) then

              p_JKQ(Fin%ggf0:Fin%ggf1) => JKQinMV

            ! If static
            else

              p_JKQ(Fin%ggf0:Fin%ggf1) => JKQC

            end if

            ! Get total frequency size
            if (LIRAM) nmfreq = sum(p_frec%mfreq)

!$omp parallel default(none) &
!$omp private(jjfreq,iifreq,kkfreq,iran,ifreq,tid,p_mfreq,dy,y0) &
!$omp private(omegai) &
!$omp shared(Fin,p_frec,p_JKQ,p_index1,p_index2,p_dx,p_warr2) &
!$omp shared(nmfreq,PRD,p_red,LIRAM,warr2,omp,dyn,omega,vfac) &
!$omp shared(CRD)

            ! If storing Warr
            if (LIRAM) then
              if (nmfreq.gt.0) then
!$omp single
                allocate(p_warr2(nmfreq))
!$omp end single
!$omp workshare
                p_warr2 = dble(p_red%IWarr2)
!$omp end workshare
              end if
            ! If not storing
            else
!$omp single
              if (allocated(Warr2)) p_warr2 => Warr2
!$omp end single
            end if

#ifdef _OPENMP
            tid = omp_get_thread_num() + 1
#endif
            ! Initialize frequency index
            jjfreq = 0

            ! For each output frequency
            iifreq = 0
            do iran=1,Fin%nran
              do ifreq=Fin%if0(iran),Fin%if1(iran)

                ! Advance index
                iifreq = iifreq + 1

                ! Point to dimension
                p_mfreq => p_frec%mfreq(iifreq)

                ! If coherent wing
                if (p_mfreq.lt.1) then

                  ! If dynamic
                  if (dyn) then

                    ! Input frequency
                    omegai = omega(ifreq)*vfac

                    ! Get JKQ
                    y0 = getJin(omega(Fin%ggf0:Fin%ggf1), &
                                p_JKQ, &
                                ifreq-Fin%ggf0+1, &
                                Fin%ggf1-Fin%ggf0+1,omegai)

                    ! Fully coherent contribution
                    PRD(ifreq) = CRD(ifreq)*y0

                  ! Static
                  else

                    ! Fully coherent contribution
                    PRD(ifreq) = CRD(ifreq)*p_JKQ(ifreq)

                  end if ! Dynamic

                  ! Skip rest
                  cycle

                end if

#ifdef _OPENMP
                ! If multi-threading
                if (omp) then

                  ! If out of range below
                  if (iifreq.lt.Fin%oif0(tid)) then
                    jjfreq = jjfreq + p_mfreq
                    cycle
                  end if
                  ! If out of range above
                  if (iifreq.gt.Fin%oif1(tid)) exit
                end if
#endif
                ! Initialize
                PRD(ifreq) = 0d0

                ! Input frequencies
                do jfreq=jjfreq+1,jjfreq+p_mfreq

                  ! Interpolate
                  y0 = p_JKQ(p_index1(jfreq))

                  ! Linear interpolation
                  dy = p_JKQ(p_index2(jfreq)) - y0

                  ! Add J00 interpolated
                  PRD(ifreq) = PRD(ifreq) + &
                               (dy*p_dx(jfreq) + y0)*p_warr2(jfreq)
                end do

                ! Update jjfreq
                jjfreq = jjfreq + p_mfreq

              end do ! output frequencies
            end do ! output frequencies ranges

!$omp end parallel

          ! If angle-dependent
          else

            ! Compute the sizes
            if (LIRAM) then
              if (ffitran.eq.ffjtran) then
                nmfreq = sum(p_frec%mfreq)*(Geom%nTh*Geom%nPh2-nfs)
              else
                nmfreq = sum(p_frec%mfreq)*Geom%nTh*Geom%nPh2
              end if
            end if

!$omp parallel default(none) &
!$omp private(jjfreq0,kkfreq,iifreq,iran,ifreq,p_mfreq,ith1) &
!$omp private(iph1,jjfreq,jfreq,y0,dy,StokesMV,cost,vfac1,omegai) &
!$omp private(StokesM,PRDin,auxb,sint,cosc,sinc,tid,llfreq) &
!$omp private(kkfreq0,iidir,tointerp) &
!$omp shared(LIRAM,nmfreq,Warr2,p_red,axiali,Fin,p_frec,Geom,p_dx) &
!$omp shared(p_index1,p_index2,Stokes,ffjtran,ffitran,stype,dyn) &
!$omp shared(vx,vy,vz,omega,vfac,omp,p_warr2) &
!$omp shared(PRDdir,PRD)

#ifdef _OPENMP
            tid = omp_get_thread_num() + 1
            llfreq = 0
!$omp workshare
            ! Initialize 2nd order part
            PRDdir = 0d0
            PRD = 0d0
!$omp end workshare
#else
            ! Initialize 2nd order part
            PRD = 0d0
#endif
            ! If storing Warr
            if (LIRAM) then
              if (nmfreq.gt.0) then
!$omp single
                allocate(p_warr2(nmfreq))
!$omp end single
!$omp workshare
                p_warr2 = dble(p_red%IWarr2)
!$omp end workshare
              end if
            else
!$omp single
              if(allocated(Warr2)) p_warr2 => Warr2
!$omp end single
            end if

            ! If axial symmetry
            if (axiali) then

              ! Initialize frequency indexes
              jjfreq0 = 0
              kkfreq0 = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fin%nran
                do ifreq=Fin%if0(iran),Fin%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Point to dimension
                  p_mfreq => p_frec%mfreq(iifreq)
#ifdef _OPENMP
                  ! Initialize direction index
                  iidir = 0
#endif
                  ! For each polar direction
                  do ith1=1,Geom%nTh

                    ! Initialize interpolated
                    tointerp = .True.

                    ! For each azimuthal direction
                    do iph1=1,Geom%nPh2
#ifdef _OPENMP
                      ! Advance direction and thread index
                      iidir = iidir + 1
                      llfreq = llfreq + 1
#endif
                      ! Special treatment if forward for two terms
                      if ((ffjtran.eq.ffitran.and. &
                           stype(iph1,ith1).lt.0).or. &
                          (p_mfreq.lt.1)) then
#ifdef _OPENMP
                        ! If multi-threading
                        if (omp) then
                          ! If out of range below
                          if (llfreq.lt.Fin%oif0(tid)) cycle
                          ! If out of range above
                          if (llfreq.gt.Fin%oif1(tid)) exit
                        end if
#endif
                        ! If there are dynamics
                        if (dyn) then

                          ! Get cosine of direction
                          cost = Geom%V_mu(ith1)

                          ! Calculate Doppler shift factor
                          vfac1 = 1d0 - vz*cost

                          ! We will be using the inverse
                          vfac1 = 1d0/vfac1

                          ! Input frequency
                          omegai = omega(ifreq)*vfac*vfac1

                          ! Get Stokes
                          StokesM = getStkinI(omega, &
                                              Stokes(:,1,ith1), &
                                              ifreq,omegai)
                        else

                          StokesM = Stokes(ifreq,1,ith1)

                        end if ! Dynamics

                        ! Add the directional weights
#ifdef _OPENMP
                        PRDdir(iidir,ifreq) = StokesM* &
                                              Geom%W_mu(ith1)* &
                                              Geom%W_mux2(iph1)
#else
                        PRD(ifreq) = PRD(ifreq) + &
                                     StokesM*Geom%W_mu(ith1)* &
                                     Geom%W_mux2(iph1)
#endif
                      else
#ifdef _OPENMP
                        ! If multi-threading
                        if (omp) then
                          ! If out of range below
                          if (llfreq.lt.Fin%oif0(tid)) then
                            kkfreq0 = kkfreq0 + p_mfreq
                            cycle
                          end if
                          ! If out of range above
                          if (llfreq.gt.Fin%oif1(tid)) exit
                        end if
#endif
                        ! If need to interpolate
                        if (tointerp) then

                          ! Update flag
                          tointerp = .False.

                          ! For each input frequency
                          do jfreq=1,p_mfreq

                            ! Advance index
                            jjfreq = jjfreq0 + jfreq

                            ! Initial
                            y0 = Stokes(p_index1(jjfreq),1,ith1)

                            ! Slope
                            dy = Stokes(p_index2(jjfreq),1,ith1) - y0

                            ! Linear interpolation
                            StokesMV(jfreq) = dy*p_dx(jjfreq) + y0

                          end do ! Input frequencies

                        end if ! Need to interpolate

                        ! Reset partial integral
                        PRDin = 0d0

                        ! For each input frequency
                        do jfreq=1,p_mfreq

                          ! Advance index
                          kkfreq = kkfreq0 + jfreq

                          ! Add contribution to the partial integral
                          PRDin = PRDin + &
                                  StokesMV(jfreq)*p_warr2(kkfreq)

                        end do ! Input frequencies

                        ! Add the directional weights
#ifdef _OPENMP
                        PRDdir(iidir,ifreq) = PRDin* &
                                              Geom%W_mu(ith1)* &
                                              Geom%W_mux2(iph1)
#else
                        PRD(ifreq) = PRD(ifreq) + PRDin* &
                                     Geom%W_mu(ith1)* &
                                     Geom%W_mux2(iph1)
#endif
                        ! Update kkfreq0
                        kkfreq0 = kkfreq0 + p_mfreq

                      end if ! Forward scattering two terms

                    end do ! azimuthal nodes

                    ! Update jjfreq
                    if (dyn) jjfreq0 = jjfreq0 + p_mfreq

                  end do ! polar nodes

                  ! Update jjfreq
                  if (.not.dyn) jjfreq0 = jjfreq0 + p_mfreq

                end do ! output frequencies
              end do ! output frequencies ranges

            ! If non-axial symmetryc
            else

              ! Initialize indexes
              jjfreq0 = 0
              kkfreq0 = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fin%nran
                do ifreq=Fin%if0(iran),Fin%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Point to dimension
                  p_mfreq => p_frec%mfreq(iifreq)
#ifdef _OPENMP
                  ! Initialize direction index and
                  iidir = 0
#endif
                  ! For each polar direction
                  do ith1=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph1=1,Geom%nPh2
#ifdef _OPENMP
                      ! Advance direction and thread index
                      iidir = iidir + 1
                      llfreq = llfreq + 1
#endif
                      ! Special treatment if forward for two terms
                      if ((ffjtran.eq.ffitran.and. &
                           stype(iph1,ith1).lt.0).or. &
                          (p_mfreq.lt.1)) then
#ifdef _OPENMP
                        ! If multi-threading
                        if (omp) then
                          ! If out of range below
                          if (llfreq.lt.Fin%oif0(tid)) cycle
                          ! If out of range above
                          if (llfreq.gt.Fin%oif1(tid)) exit
                        end if
#endif
                        ! If there are dynamics
                        if (dyn) then

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

                          ! Input frequency
                          omegai = omega(ifreq)*vfac*vfac1

                          StokesM = getStkinI(omega, &
                                              Stokes(:,iph1,ith1), &
                                              ifreq,omegai)

                        else

                          StokesM = Stokes(ifreq,iph1,ith1)

                        end if ! Dynamics
#ifdef _OPENMP
                        PRDdir(iidir,ifreq) = StokesM* &
                                              Geom%W_mu(ith1)* &
                                              Geom%W_mux2(iph1)
#else
                        PRD(ifreq) = PRD(ifreq) + &
                                     StokesM*Geom%W_mu(ith1)* &
                                     Geom%W_mux2(iph1)
#endif
                      else
#ifdef _OPENMP
                        ! If multi-threading
                        if (omp) then
                          ! If out of range below
                          if (llfreq.lt.Fin%oif0(tid)) then
                            kkfreq0 = kkfreq0 + p_mfreq
                            if (dyn) jjfreq0 = jjfreq0 + p_mfreq
                            cycle
                          end if
                          ! If out of range above
                          if (llfreq.gt.Fin%oif1(tid)) exit
                        end if
#endif
                        ! Reset partial integral
                        PRDin = 0d0

                        ! For each input frequency
                        do jfreq=1,p_mfreq

                          ! Advance indexes
                          jjfreq = jjfreq0 + jfreq
                          kkfreq = kkfreq0 + jfreq

                          ! Linear interpolation
                          y0 = Stokes(p_index1(jjfreq),iph1,ith1)
                          dy = Stokes(p_index2(jjfreq),iph1,ith1) - y0
                          StokesM = dy*p_dx(jjfreq) + y0

                          ! Add contribution to the partial integral
                          PRDin = PRDin + StokesM*p_warr2(kkfreq)

                        end do ! Input frequencies

                        ! Add the directional weights
#ifdef _OPENMP
                        PRDdir(iidir,ifreq) = PRDin* &
                                              Geom%W_mu(ith1)* &
                                              Geom%W_mux2(iph1)
#else
                        PRD(ifreq) = PRD(ifreq) + PRDin* &
                                     Geom%W_mu(ith1)* &
                                     Geom%W_mux2(iph1)
#endif
                        ! Update indexes
                        if (dyn) jjfreq0 = jjfreq0 + p_mfreq
                        kkfreq0 = kkfreq0 + p_mfreq

                      end if ! Forward scattering two terms

                    end do ! azimuthal nodes
                  end do ! polar nodes

                  ! Advance
                  if (.not.dyn) jjfreq0 = jjfreq0 + p_mfreq

                end do ! output frequencies
              end do ! output frequencies ranges

            end if ! Axial symmetry
!$omp end parallel

#ifdef _OPENMP
            ! For each output frequency
            do iran=1,Fin%nran
              do ifreq=Fin%if0(iran),Fin%if1(iran)

                ! Get proper PRD quantity
                PRD(ifreq) = PRD(ifreq) + sum(PRDdir(:,ifreq))

              end do ! Output frequencies
            end do ! Output frequency range
#endif
          endif ! AV or AD

          ! Clean warr2
          if (LIRAM) deallocate(p_warr2)
          nullify(p_warr2)

          ! Clean interpolation data
          if (.not.p_frec%RAM) deallocate(p_index1,p_index2,p_dx)
          nullify(p_index1,p_index2,p_dx)

          !
          ! Flat spectrum contribution. Implicit branching ratio
          !

          ! For each output frequency
          do iran=1,Fin%nran
            do ifreq=Fin%if0(iran),Fin%if1(iran)

              ! Substract the flat spectrum part due to just
              ! radiative excitation
              PRD(ifreq) = PRD(ifreq) - Jrad

            end do ! output frequencies
          end do ! output frequencies ranges

          ! Apply hanle factor and Einstein coefficient
          auxb = rhoc*Atom%fst(itran)%Blu(iJl,iJu)/hanleden
          eps2 = eps2 + PRD*auxb

        end do ! Lower J levels
      end do ! Lower terms

      ! Calculate factor for ALI and eps2
      auxb = 1d-8*IPI2/c/Dw
      rpF = eps2*auxb
      auxb = 1d-2*IPI4*Atom%fst(jtran)%Aul(iJu,iJf)
      eps2 = rpF*CRD*auxb
      rpF = rpF/rhou + 1d0

      ! Clean pointers
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_JKQ)) nullify(p_JKQ)
      if (associated(p_warr2)) nullify(p_warr2)

      return

      end subroutine emissI2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient for a LTE line.\n
      !!       line(LTEline_class): Structure with the LTE line data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!      aprof(LTEprof_class): Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!                pE(dfloat): Unit transformation factor\n
      !!            eta(dfloat(:)): Intensity absorptivity
      subroutine absorbILTE(line,omega,iz,if0,if1,aprof,Dw,vfac,pE, &
                            eta)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(LTEprof_class), intent(in):: aprof
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw, pE, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq

      double precision:: at,feta,Dfreqw,vfacw,prof


      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE/Dw

      ! If stored in RAM
      if (aprof%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(feta,aprof,eta)
        eta = feta*aprof%p
!$omp end parallel workshare

      ! Not stored
      else

        ! Level quantities

        ! Damping parameter
        at = line%damp(iz)/Dw

        ! Energy
        Dfreqw = (line%Eu - line%El)/Dw

        ! Shift
        vfacw = vfac/Dw

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feta,Dfreqw,omega,vfacw,at,eta)
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          eta(ifreq) = feta*prof

        end do ! frequencies
!$omp end parallel do

      end if ! Type of profile calculation

      return

      end subroutine absorbILTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emission coefficient.\n
      !!       line(LTEline_class): Structure with the LTE line data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!      aprof(LTEprof_class): Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!            eps(dfloat(:)): Intensity emissivity\n
      subroutine emissILTE(line,omega,iz,if0,if1,aprof,Dw,vfac,eps)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(LTEprof_class), intent(in):: aprof
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps

      ! Local

      integer:: ifreq

      double precision:: at,Dfreqw,vfacw,prof,feps


      !
      ! Get population factor
      !

      ! Emissivity factor
      feps = line%nu(iz)*1d3*IPI41*line%Aul/Dw

      ! If stored in RAM
      if (aprof%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(feps,aprof,eps)
        eps = feps*aprof%p
!$omp end parallel workshare

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        at = line%damp(iz)/Dw

        ! Energy
        Dfreqw = (line%Eu - line%El)/Dw

        ! Shift
        vfacw = vfac/Dw

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feps,Dfreqw,omega,vfacw,at,eps)
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          eps(ifreq) = feps*prof

        end do ! frequencies
!$omp end parallel do

      end if ! Type of profile calculation

      return

      end subroutine emissILTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient for a LTE line.\n
      !!       line(LTEline_class): Structure with the LTE line data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!      aprof(LTEprof_class): Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!                pE(dfloat): Unit transformation factor\n
      !!            eta(dfloat(:)): Intensity absorptivity\n
      !!            eps(dfloat(:)): Intensity emissivity
      subroutine rt1ordILTE(line,omega,iz,if0,if1,aprof,Dw,vfac,pE, &
                            eta,eps)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(LTEprof_class), intent(in):: aprof
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw, pE, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq

      double precision:: at,feta,feps,Dfreqw,vfacw,prof


      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE/Dw

      ! Emissivity factor
      feps = line%nu(iz)*1d3*IPI41*line%Aul/Dw

      ! If stored in RAM
      if (aprof%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(feta,feps,aprof,eta,eps)
        eta = feta*aprof%p
        eps = feps*aprof%p
!$omp end parallel workshare

      ! Not stored
      else

        ! Level quantities

        ! Damping parameter
        at = line%damp(iz)/Dw

        ! Energy
        Dfreqw = (line%Eu - line%El)/Dw

        ! Shift
        vfacw = vfac/Dw

        ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feta,feps,Dfreqw,omega,vfacw,at,eta,eps)
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          eta(ifreq) = feta*prof
          eps(ifreq) = feps*prof

        end do ! frequencies
!$omp end parallel do

      end if ! Type of profile calculation

      return

      end subroutine rt1ordILTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorptivity due to photoionization.\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!       itran(integer): Transition index\n
      !!     ilevell(integer): Lower level index\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!       eta(dfloat(:)): Absorptivity
      subroutine photoabsI(Atom,itran,ilevell,iz,if0,if1,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran, ilevell, iz, if0, if1
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq, iterml, iJl, iR

      double precision:: rJl,rhol


      !
      ! Get indexes
      !

      ! Term index
      iterml = Atom%term(ilevell)

      ! J level index
      iJl = Atom%sublevel(ilevell)

      ! Angular momentum
      rJl = Atom%rJval(iJl,iterml)

      ! SEE index
      iR = Atom%irho(iterml)%Jrho(iJl,iJl)%kq(0,0)

      ! Population lower level
      rhol = sqrt(2d0*rJl+1d0)*dble(Atom%crho(iR,iz))


      ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq) shared(eta,if0,if1,Atom,itran,rhol)
      do ifreq=if0,if1

        ! Compute absorptivity
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*rhol

      end do ! frequencies
!$omp end parallel do

      return

      end subroutine photoabsI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emissivity due to recombination\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!     omega(dfloat(:)): Frequency array\n
      !!            T(dfloat): Temperature\n
      !!           ne(dfloat): Electron density\n
      !!       itran(integer): Transition index\n
      !!     ilevelu(integer): Upper level index\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!       eps(dfloat(:)): Emissivity\n
      !!       eta(dfloat(:)): Stimulated emissivity\n
      !!         rhou(dfloat): Factor for Lambda operator
      subroutine photoepsI(Atom,omega,T,ne,itran,ilevelu,iz,if0,if1, &
                           eps,eta,rhou)

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran, ilevelu, iz, if0, if1
      double precision, intent(in):: T, ne
      double precision, dimension(:), intent(in):: omega
      double precision, intent(out):: rhou
      double precision, dimension(if0:if1), intent(out):: eta, eps

      integer:: ifreq, iJu, itermu, iR

      double precision:: c0,c1,exu,iexu,pE,Saha
      double precision:: rJu,tmp
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

      ! Get upper level rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Get upper level population
      tmp = sqrt(2d0*rJu+1d0)*rhou

      ! Apply Saha factor
      tmp = tmp*Saha

      ! Compute exponential argument constant
      c0 = c2*1d4/T

      ! Compute energy constant part
      c1 = 2d21*c*convF

      ! Compute numerator of rpf
      Saha = 1d0/Saha/sqrt(2d0*rJu+1d0)

      omega3 = omega(if0:if1)
      omega3 = omega3*omega3*omega3

      ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,iexu,exu,pE) &
!$omp shared(eta,eps,if0,if1,c0,c1,omega,omega3,tmp,Atom,itran)
      do ifreq=if0,if1

        iexu = c0*omega(ifreq)
        exu = diexp(iexu)

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu*tmp

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies
!$omp end parallel do

      return

      end subroutine photoepsI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emissivity due to recombination with
      !! precomputed quantities\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    omega3(dfloat(:)): Frequency array to the cube\n
      !!       exu(dfloat(:)): Exponential for emissivity\n
      !!      iexu(dfloat(:)): Inverse of the exponential for
      !!                       emissivity\n
      !!            T(dfloat): Temperature\n
      !!           ne(dfloat): Electron density\n
      !!       itran(integer): Transition index\n
      !!     ilevelu(integer): Upper level index\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!       eps(dfloat(:)): Emissivity\n
      !!       eta(dfloat(:)): Stimulated emissivity\n
      !!         rhou(dfloat): Factor for Lambda operator
      subroutine photoepsIS(Atom,omega3,exu,T,ne,itran, &
                            ilevelu,iz,if0,if1,eps,eta,rhou)

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran, ilevelu, iz, if0, if1
      double precision, intent(in):: T, ne
      double precision, dimension(if0:if1), intent(in):: omega3
      double precision, dimension(if0:if1), intent(in):: exu
      double precision, intent(out):: rhou
      double precision, dimension(if0:if1), intent(out):: eta, eps

      integer:: ifreq, iJu, itermu, iR

      double precision:: c0,c1,pE,Saha
      double precision:: rJu,tmp


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

      ! Get upper level rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Get upper level population
      tmp = sqrt(2d0*rJu+1d0)*rhou

      ! Apply Saha factor
      tmp = tmp*Saha

      ! Compute exponential argument constant
      c0 = c2*1d4/T

      ! Compute energy constant part
      c1 = 2d21*c*convF

      ! Compute numerator of rpf
      Saha = 1d0/Saha/sqrt(2d0*rJu+1d0)

      ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,pE) &
!$omp shared(eta,eps,if0,if1,c1,omega3,exu,tmp,Atom,itran)
      do ifreq=if0,if1

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu(ifreq)*tmp

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies
!$omp end parallel do

      return

      end subroutine photoepsIS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolates intensity to the requested frequency\n
      !!     omega(dfloat(:)): Frequency array\n
      !!    Stokes(dfloat(:)): Intensity\n
      !!       ifreq(integer): Frequency index of the output frequency
      !!                       associated to the requested input
      !!                       frequency\n
      !!            x(dfloat): Input frequency to interpolate into
      double precision function getStkinI(omega,Stokes,ifreq,x)

      ! I/O
      integer, intent(in):: ifreq
      double precision, intent(in):: x
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: Stokes

      ! Local
      integer:: jfreq
      double precision:: dx, dy

      ! Initialize as equals
      getStkinI = Stokes(ifreq)

      ! If omegai > omega(ifreq)
      if (x.gt.omega(ifreq)) then

        ! If out of right boundary
        if (x.ge.omega(nfreq)-TINYO) then

          getStkinI = Stokes(nfreq)

          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,nfreq-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getStkinI = Stokes(jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq).and. &
                    x.lt.omega(jfreq+1)) then

              dy = Stokes(jfreq+1) - Stokes(jfreq)

              dx = x - omega(jfreq)

              getStkinI = dx*dy/(omega(jfreq+1) - omega(jfreq)) + &
                          Stokes(jfreq)

              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries

      ! If omegai < omegao
      else if (x.lt.omega(ifreq)) then

        ! If out of left boundary
        if (x.le.omega(1)+TINYO) then

          getStkinI = Stokes(1)

          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,2,-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getStkinI = Stokes(jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq-1).and. &
                    x.lt.omega(jfreq)) then

              dy = Stokes(jfreq) - Stokes(jfreq-1)

              dx = x - omega(jfreq-1)

              getStkinI = dx*dy/(omega(jfreq) - omega(jfreq-1)) + &
                          Stokes(jfreq-1)

              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries
      end if ! omegai > omegao

      return

      end function getStkinI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolates J00 to the requested frequency\n
      !!  omega(dfloat(:)): Frequency array\n
      !!    JKQ(dfloat(:)): Radiation field tensor\n
      !!    ifreq(integer): Frequency index of the output frequency
      !!                    associated to the requested input
      !!                    frequency (shifted to limited vector)\n
      !!    nfreq(integer): Size of input vectors\n
      !!                    associated to the requested input
      !!                    frequency (shifted to limited vector)\n
      !!         x(dfloat): Input frequency to interpolate into
      function getJin(omega,JKQ,ifreq,nfreq,x)

      ! I/O
      integer, intent(in):: ifreq,nfreq
      double precision, intent(in):: x
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: JKQ

      double precision:: getJin

      ! Local
      integer:: jfreq
      double precision:: dxs,dys

      ! Initialize as equals
      getJin = JKQ(ifreq)

      ! If omegai > omega(ifreq)
      if (x.gt.omega(ifreq)) then

        ! If out of right boundary
        if (x.ge.omega(nfreq)-TINYO) then

          getJin = JKQ(nfreq)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,nfreq-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getJin = JKQ(jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq).and. &
                    x.lt.omega(jfreq+1)) then

              dys = JKQ(jfreq+1) - JKQ(jfreq)

              dxs = x - omega(jfreq)

              getJin = dxs*dys/(omega(jfreq+1) - omega(jfreq)) + &
                       JKQ(jfreq)

              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries

      ! If omegai < omegao
      else if (x.lt.omega(ifreq)) then

        ! If out of left boundary
        if (x.le.omega(1)+TINYO) then

          getJin = JKQ(1)

          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,2,-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getJin = JKQ(jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq-1).and. &
                    x.lt.omega(jfreq)) then

              dys = JKQ(jfreq) - JKQ(jfreq-1)

              dxs = x - omega(jfreq-1)

              getJin = dxs*dys/(omega(jfreq) - omega(jfreq-1)) + &
                       JKQ(jfreq-1)

              return

            end if ! Check output frequency

          end do ! Searching frequency

        end if ! Within boundaries
      end if ! omegai > omegao

      end function getJin

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeffiaux_mod
