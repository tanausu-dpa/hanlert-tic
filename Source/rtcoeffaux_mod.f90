      !> Radiation transfer coefficients
      module rtcoeffaux_mod
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
!     04/27/2017
!  Last version:
!     02/11/2025 V3.1.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     02/11/2025:    V3.1.5 - Bugfix: failing to consider forward
!                             scattering for Raman scattering in
!                             the second order emissivity (TdPA)
!
!     04/02/2024:    V3.1.4 - Debugged rt1ordNB and rt1ord (TdPA)
!
!     04/01/2024:    V3.1.3 - Added rt1ordNB and rt1ord to compute
!                             together absorption and emission. The
!                             latter is optional and requires using
!                             the energy eigenbasis for the dipole
!                             strength (TdPA)
!                           - Fixed some comments (TdPA)
!                           - Carried out one of the operations in
!                             tempR in emiss and emissNB (TdPA)
!
!     11/29/2023:    V3.1.2 - With coherent scattering, wrong use of
!                             the JKQ tensors with Q < 0 (TdPA)
!
!     11/24/2023:    V3.1.1 - Enter calculation of Warr2 and
!                             deallocation of p_warr2 pointer only if
!                             there were frequencies (TdPA)
!                           - Bugfix: getJKQinnu could be called
!                             with out of range indexes (TdPA)
!
!     11/14/2023:    V3.1.0 - Reworked emiss2ord significantly into
!                             emiss2ord_AA and emiss2ord_AD (TdPA)
!                           - Reworked emiss2ordNB significantly into
!                             emiss2ordNB_AA and emiss2ordNB_AD (TdPA)
!                           - Renamed getStokin to getStokinnu and
!                             added a new getStokin (TdPA)
!                           - Renamed getJKQin to getJKQinnu and added
!                             a new getJKQin (TdPA)
!
!     10/31/2023:    V3.0.8 - Placeholder change to use the
!                             pre-computed cosines and sines of
!                             scattering angles in angle-average, but
!                             a significant update is expected for
!                             next version (TdPA)
!
!     09/29/2023:    V3.0.7 - Updated to term- and transition-wise
!                             K cut limits (TdPA)
!                           - Avoid computing tensor components
!                             with negative Q values in function
!                             getJKQstar (TdPA)
!
!     08/07/2023:    V3.0.6 - Added absorbLTE and rt1ordLTE (TdPA)
!
!     03/23/2023:    V3.0.5 - Bugfix: Ensured the OpenMP version
!                             compiles after some years of changes,
!                             albeit did not test it works (TdPA)
!
!     11/24/2022:    V3.0.4 - Bugfix: Missing definition of cost
!                             for the additive ad-hoc assymetry part
!                             in getJKQstar (TdPA)
!
!     11/10/2022:    V3.0.3 - JKQa are new inputs to emiss2ord,
!                             emiss2ordNB, and getJKQstar (TdPA)
!                           - Added a rotation of the JKQ(k) in
!                             emiss2ord when the dynamic AA tensors
!                             are calculated and there are ad-hoc
!                             asymmetries (TdPA)
!                           - getJKQstar now adds the ad-hoc
!                             asymmetries as well. The resulting
!                             tensors in this case are in the
!                             vertical reference frame, and not in
!                             the magnetic field frame (TdPA)
!
!     10/26/2022:    V3.0.2 - Changed the storage structure of the
!                             rdip variable (TdPA)
!                           - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.1 - Nullify pointers as when starting each
!                             routine (TdPA)
!                           - Clean pointers at exit (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case Atmo%v has
!                             changed to Atmo%vx,%vy, and %vz (TdPA)
!
!     06/21/2022:    V2.2.0 - Added coherent scattering in the
!                             observers frame (TdPA)
!                             NOTE: Limited testing (AA, static, and
!                             non-magnetic)
!                           - Added function getJKQin to interpolate
!                             JKQ to a given frequency (TdPA)
!
!     05/24/2022:    V2.1.0 - Bugfix in a commented part of the code,
!                             so no actual change (TdPA)
!
!     05/17/2022:    V2.1.0 - Bugfix?: Changed the way the
!                             redistribution function integral is
!                             normalized. The old way (normalize only
!                             real part) had issues because it broke
!                             a necessary balance between the real and
!                             imaginary parts when combined with the
!                             outgoing geometrical tensors (TdPA)
!                             NOTE: This did not affect previous
!                             calculations with angle-averaged
!                             redistribution or angle-dependent
!                             assuming axial symmetry, as long as
!                             4 or 8 azimuthal nodes were considered.
!
!     02/11/2022:    V2.0.3 - Bugfix: In the non-magnetic absorptivity
!                             the MO coefficients were computed as the
!                             dichroic coefficients, with the real
!                             part instead of imaginary. Luckily, this
!                             could be only a problem if there was no
!                             magnetic field nor axial symmetry (TdPA)
!
!     03/24/2021:    V2.0.2 - Not using OpenMP to initialize leps
!                             varibles in emiss2ord (TdPA)
!                           - Bugfix: Comparing with wrong array in
!                             the iU1 loop of emiss2ord with activated
!                             OpenMP (TdPA)
!
!     03/23/2021:    V2.0.1 - One of the changes in the angle-averaged
!                             branch in emiss2ordNB led to a
!                             measurable decrease of performance. I
!                             assume it was related to cache misses,
!                             so I came back to the legacy structure
!                             of that branch (TdPA)
!                           - Bugfix: The auxiliar leps variables
!                             in emiss2ord used with OpenMP were
!                             initialized one loop above what they
!                             should, completely breaking any lambda
!                             type transition (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Deep changes to include OpenMP (TdPA)
!                           - Added routines with long sections
!                             of code repeated in both emiss2ord and
!                             emiss2ordNB: getJKQstar, getscatter,
!                             and getinterpolation (TdPA)
!
!     01/13/2021:   V1.15.2 - Bugfix: ensure that Atom%NCHLT is
!                             allocated in order to check it (TdPA)
!
!     10/26/2020:   V1.15.1 - Now index1, index2, and dx are pointers
!                             to avoid copying (TdPA)
!
!     09/11/2020:   V1.15.0 - Completely changed the Frec and Red
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
!                           - Do not touch the imaginary part when
!                             integrating, ever (TdPA)
!                           - Put back the explicit integrals in
!                             the emiss2ord and emiss2ordNB routines
!                             to make use of the rolling indexes in
!                             an easier way, thus removing the
!                             routines calcJ2ordInt, calcSF2ordInt,
!                             and calcS2ordInt (TdPA)
!
!     06/26/2020:   V1.14.0 - Changed when the non-coherent lower
!                             term approximation is applied in
!                             emiss2ord and added it to absorb as
!                             well (TdPA)
!                           - Added an option to apply the extra
!                             Kcut or not in absorb1 (TdPA)
!                           - Changed the extra cut of K orders in
!                             absorb1, because we (RC and myself)
!                             think it makes more sense (TdPA)
!
!     06/05/2020:   V1.13.0 - Changed the order of the JradC indexes
!                             to speed-up emiss2ord (JD)
!                           - Added extra pointers in the
!                             interpolation section (JD)
!                           - Split the interpolation stored indexes
!                             into two vectors to speed-up the
!                             emiss2ord routine (JD)
!                           - Created a new routine to carry out the
!                             radiation field interpolation (JD)
!                           - Extended all previous changes to all
!                             possible paths of the 2nd order
!                             emissivity (TdPA)
!
!     06/01/2020:   V1.12.2 - Moved the Hanle denominator, hanleden,
!                             out some loops. Also pre-computed its
!                             real part, to avoid the repetition of a
!                             multiplication and a division (TdPA)
!                           - Introduced the non-coherent lower term
!                             approximation in the redistribution
!                             function (TdPA)
!                           - The Kcut affecting rJu2 and rJu3 is
!                             commented right now (TdPA)
!
!     11/19/2019:   V1.12.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     11/13/2019:   V1.12.0 - Added support for voigt profiles read
!                             from a file (TdPA)
!
!     10/16/2019:   V1.11.0 - The cosine and sine of the scattering
!                             angle for the redistribution function
!                             are computed at the beginning, and
!                             passed directly to Wfunc. This saves a
!                             lot of calls to cos() and sin() (TdPA)
!
!     10/03/2019:   V1.10.1 - trano is now indexed in Red_class (TdPA)
!
!     09/13/2019:   V1.10.0 - Added branch in emiss2ord and
!                             emiss2ordNB for angle-average
!                             redistribution for the dynamic case,
!                             following the strategy of angle-average
!                             in the comoving frame (TdPA)
!
!     05/31/2019:    V1.9.4 - Added the emerging variable to the
!                             second order emissivity to distinguish
!                             the angles that are being used to
!                             compute the scattering angle (TdPA)
!                           - Moved the unit normalization factor for
!                             CRD out of the loops (TdPA)
!                           - Added an underflow check for the
!                             redistribution profiles before type
!                             casting them into single precision in
!                             its structure (TdPA)
!
!     05/28/2019:    V1.9.3 - Bugfix: Copied and pasted the change
!                             to getStkin from rtcoeffiaux_mod and
!                             forgot to change it (TdPA)
!                           - Forgot to change the declarations of
!                             prof1 and prof2 to just prof (TdPA)
!
!     05/08/2019:    V1.9.2 - Split computation of the two components
!                             of the CRD profile in emiss2ord to avoid
!                             computing the same profile more than
!                             once (TdPA)
!                           - Introduced a check before storing
!                             redistribution profiles because the type
!                             casting from double to single could
!                             end up in NaN values (TdPA)
!                           - The interpolator getStkin initializes
!                             as equal instead of waiting to the last
!                             momment (TdPA)
!
!     04/09/2019:    V1.9.1 - Wrong routine name in emiss2ord (TdPA)
!
!     02/20/2019:    V1.9.0 - New verbosity (TdPA)
!                           - Using specific TINY variables (TdPA)
!                           - Using diexp function (TdPA)
!
!     02/11/2019:    V1.8.3 - Bugfix: The interpolation for the
!                             dynamic case was wrong, moved part of
!                             the interpolation to omegabuild (TdPA)
!
!     02/08/2019:    V1.8.2 - Bugfix: Introduced a check for when the
!                             redistribution function integrates to 0
!                             to avoid 0/0 (TdPA)
!
!     10/01/2018:    V1.8.1 - Bugfix: Missing normalization factor
!                             for units in the CRD profile of
!                             emiss2ord when the option to store Voigt
!                             profiles was active (TdPA)
!
!     08/06/2018:    V1.8.0 - Added the possibility to get Voigt
!                             profiles from RAM (TdPA)
!                           - Changed the argument Norm in the b-b
!                             subroutines (TdPA)
!                           - Added photoepsS subroutine, that takes
!                             pre-computed frequency quantities from
!                             the RAM (TdPA)
!                           - Made changes to emiss2ord and
!                             emiss2ordNB that allows for partial
!                             storage of Wfunc2 (TdPA)
!
!     08/03/2018:    V1.7.3 - Changed the implementation of the Kcut
!                             once again. It was wrong to cut on the
!                             J differences (TdPA)
!                           - Bugfix: The K argument for absorb1 in
!                             absorb must be 0, not K (TdPA)
!
!     07/27/2018:    V1.7.2 - Changed the implementation of the Kcut.
!                             It is only really consistent for cuts
!                             to K=0 (TdPA)
!                           - Changed the implementation of the Kcut.
!
!     05/17/2018:    V1.7.1 - Forgot the CRD term in the collinear
!                             scattering case, because it is not
!                             exactly the same than the only
!                             intensity case (TdPA)
!
!     05/16/2018:    V1.7.0 - Added proper angle-dependent case (TdPA)
!                           - Split the angle-dependent part in axial
!                             and non-axial symmetryc (TdPA)
!
!     12/19/2017:    V1.6.1 - Bugfix: The imaginary part of the CRD
!                             profiles in emiss2ord and emiss2ordNB
!                             was not being normalized. At the point
!                             of writing this changelog, I do not know
!                             how severe the effect is (TdPA)
!
!     12/05/2017:    V1.6.0 - Can neglect Raman scattering (TdPA)
!
!     10/30/2017:    V1.5.1 - The limit of K for a non-magnetic rate
!                             is 2, given by the TKQ tensor (TdPA)
!
!     10/23/2017:    V1.5.0 - Warr2 is stored into a single precision
!                             variable (TdPA)
!
!     10/13/2017:    V1.4.1 - Avoids redundant calculations of the
!                             interpolation + integral. Swapped the
!                             K and K1 loops and added a flag (TdPA)
!
!     10/03/2017:    V1.4.0 - Non magnetic versions of all the b-b
!                             radiation transfer coefficients (TdPA)
!
!     09/22/2017:    V1.3.0 - Possibility to limit K (TdPA)
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
!     07/19/2017:    V1.0.5 - Using the range structure also in the
!                             variable Warr2 (TdPA)
!
!     06/19/2017:    V1.0.4 - Added storage of Warr2 (TdPA)
!                           - Changed inputs to reflect the changes in
!                             rtcoeff (TdPA)
!
!     06/16/2017:    V1.0.3 - Changes in emiss2ord towards
!                             optimization (TdPA)
!
!     06/15/2017:    V1.0.2 - Bugfix: The RTcoefficients need to be
!                             initialized to 0 in emiss and absorb
!                             because they are cumulative (TdPA)
!                           - Bugfix: cyclical reference at = (at +...
!                             -> at = (au + ... (TdPA)
!
!     06/12/2017:    V1.0.1 - The limits for the transition are
!                             passed as arguments (TdPA)
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
!  rt1ordNB:
!    This subroutine calculates the absorption and emission
!    coefficients without magnetic field.
!    Units for absorption are cm^-1 (true absorption coefficient)
!    and needs to be multiplied by the actual atomic density
!    The output value for emissivity is given in number of photons per
!    unit interval of time (s) and normalized frequency (in units of
!    Doppler width), emitted by a unit volume of gas (cm^-3) of unit
!    atomic density, within one steradian. In order to compute
!    photometric values of the intensity, the output needs be
!    multiplied by the actual atomic density.
!    This routine combines absorbNB and emissNB together
!
!  rt1ord:
!    Like rt1ordNB, but for the magnetic field case.
!    This routine combines absorb and emiss together, but it makes
!    use of the energy representation for the dipole strength in order
!    to achieve an optimal combination
!
!  absorb:
!    This subroutine calculates the absorption coefficients.
!    Units are cm^-1 (true absorption coefficient)
!    Needs to be multiplied by the actual atomic density
!
!  absorbNB:
!    Like absorb but in the unmagnetized case
!
!  absorb1:
!    This subroutine performs the two most internal loops in
!    the calculation of the absorption coefficients (needed by
!    absorb)
!
!  emiss:
!    This subroutine calculates the emission Stokes vector
!    The output value is given in number of photons per unit interval
!    of time (s) and normalized frequency (in units of Doppler width),
!    emitted by a unit volume of gas (cm^-3) of unit atomic density,
!    within one steradian.
!    In order to compute photometric values of the intensity, the
!    output needs be multiplied by the actual atomic density.
!
!  emissNB:
!    Like emiss but in the unmagnetized case
!
!  emiss1:
!    This subroutine performs the two most internal loops in
!    the calculation of the emission coefficients (needed by
!    emiss)
!
!  emiss2ord_AA:
!    This subroutine calculates the coherent scattering coefficients
!  in the angle-averaged approximation
!
!  emiss2ord_AD:
!    This subroutine calculates the coherent scattering coefficients
!  with angle-dependent redistribution
!
!  emiss2ordNB_AA:
!    Like emiss2ord_AA but in the unmagnetized case
!
!  emiss2ordNB_AD:
!    Like emiss2ord_AD but in the unmagnetized case
!
!  absorbLTE:
!    Calculate the absorption coefficient of an LTE line
!    Units are cm^-1 (true absorption coefficient)
!    Needs to be multiplied by the actual atomic density
!
!  rt1ordLTE:
!    Calculate the absorption and emission coefficient of an LTE
!    line.
!    Units (absorb) are cm^-1 (true absorption coefficient)
!    The output value (eimss) is given in number of photons per unit
!    interval of time (s) and normalized frequency (in units of
!    Doppler width), emitted by a unit volume of gas (cm^-3) of unit
!    atomic density, within one steradian.
!    In order to compute photometric values of the intensity, the
!    output needs be multiplied by the actual atomic density.
!
!  photoabs:
!    This subroutine calculates the absorption coefficients for b-f
!    transitions
!
!  photoeps:
!    This subroutine calculates the emission coefficients for b-f
!    transitions
!
!  getJKQstar:
!    Compute JKQ for angle-averaged with velocities from the Stokes
!    parameters
!
!  getscatter:
!    Get scattering angles for angle-dependent redistribution
!
!  getinterpolation:
!    Gets the coefficients to linearly interpolate the input radiation
!    field when they cannot be stored
!
!  getStkinnu:
!    This function interpolates Stokes parameters for the forward
!    2-level scattering case
!
!  getStkin:
!    This function interpolates Stokes parameters into the input
!    frequency axis
!
!  getJKQinnu:
!    This function interpolates JKQ tensor components for the forward
!    coherent scattering in the observers frame case
!
!  getJKQin:
!    This function interpolates JKQ tensor components into the input
!    frequency axis
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
                                  TINYER , TINYB , TINYO , sqrt3, &
                                  IPI, c , c2 , convF , cSaha , &
                                  fktoJ , kb, IPI41, IPI42, vrfrac , &
                                  sqrt5 , B2LK
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption and emission coefficients in absence
      !! of magnetic fields.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity\n
      !!           eps0(dfloat(:)): Intensity emissivity\n
      !!           eps1(dfloat(:)): Q emissivity\n
      !!           eps2(dfloat(:)): U emissivity\n
      !!           eps3(dfloat(:)): V emissivity\n
      !!           rhs1(dfloat(:)): Q 'dichroic' emissivity\n
      !!           rhs2(dfloat(:)): U 'dichroic' emissivity\n
      !!           rhs3(dfloat(:)): V 'dichroic' emissivity
      subroutine rt1ordNB(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                         iz,if0,if1,Norma,Dw,vfac,absK,aprof, &
                         eta0,eta1,eta2,eta3,rha1,rha2,rha3, &
                         eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, absK, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB


      ! Local

      integer:: ifreq,K,iQ,iU,iU1,iL,iL1,iR
      integer:: Kmin,Kmax

      double precision:: rLu,rLl,S,rJu,rJu1,rJl,rJl1
      double precision:: f61,f62a,f62e,f63,f64
      double precision:: eu,el,rK,au,al,aul
      double precision:: at,Dfreq,vfacw,tempRe,tempRa

      complex(kind=8),dimension(if0:if1):: prof, profK


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

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

!$omp parallel default(none) &
!$omp private(iU,iU1,iL,Dfreq,ifreq,iL1,K,iQ,rJu,eu,f61,rJl,el) &
!$omp private(rJl1,rJu1,f63,Kmin,Kmax,rK,iR,f64,f62a,f62e) &
!$omp private(tempRe,tempRa) &
!$omp shared(Atom,Flgsg,Dw,at,S,rLu,rLl,vfacw,vpfil,prof,itran) &
!$omp shared(profk,eta0,eta1,eta2,eta3,TB,Norma,rha1,rha2,rha3) &
!$omp shared(itermu,iterml,aprof,if0,if1,omega,iz,absK)

      !
      ! Common part
      !

      ! For each Ju
      do iU=1,Atom%nJ(itermu)

        ! Get Ju
        rJu = Atom%rJval(iU,itermu)

        ! Get eigenvalue upper level
        eu = Atom%FSfreq(iU,itermu)/Dw

        ! Factor Ju
        f61 = 2d0*rJu + 1d0

        ! For each Jl
        do iL=1,Atom%nJ(iterml)

          ! Get Jl
          rJl = Atom%rJval(iL,iterml)

          ! 6-j
          f62e = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

          if (abs(f62e).lt.TINYJS) cycle

          ! Absorption factor
          f62a = f62e*f61*Flgsg%sg(nint(1d0+rJu+rJl))* &
                 sqrt(2d0*rJl+1d0)

          ! Emission factor
          f62e = f62e*sqrt(f61)*(2d0*rJl+1d0)

          ! Get eigenvalue lower level
          el = Atom%FSfreq(iL,iterml)/Dw

          !
          ! Compute profile
          !

          ! If in file
          if (vpfil) then

!$omp workshare
            prof = aprof(:,Atom%i_Vind(itran)%indNB(iL,iU))
!$omp end workshare

          ! If stored
          else if (Norma%VRAM) then

!$omp workshare
            prof = Norma%prof(iL,iU,1,1)%cp
!$omp end workshare

          ! Not stored
          else

            ! Shift term
            Dfreq = eu - el

            ! For each frequency
!$omp do
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreq - omega(ifreq)*vfacw,at,prof(ifreq))

            end do ! frequencies
!$omp end do

            ! Normalize profile
!$omp workshare
            prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,1,1), &
                          dimag(prof))
!$omp end workshare

          end if ! Storing

          !
          ! Absorption
          !

          ! For each Jl'
          do iL1=1,Atom%nJ(iterml)

            ! Get Jl
            rJl1 = Atom%rJval(iL1,iterml)

            ! 6-j
            f63 = fun6j(rLu,rLl,1d0,rJl1,rJu,S,Flgsg)

            if (abs(f63).lt.TINYJS) cycle

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

!$omp workshare
                profk = f64*f63*Flgsg%sg(K)*prof*Atom%crho(iR,iz)
!$omp end workshare

!$omp workshare
                ! Absorptivity
                eta0 = eta0 + dble(TB(0,iQ,K)*profk)
                eta1 = eta1 + dble(TB(1,iQ,K)*profk)
                eta2 = eta2 + dble(TB(2,iQ,K)*profk)
                eta3 = eta3 + dble(TB(3,iQ,K)*profk)

                ! Dispersion
                rha1 = rha1 + dimag(TB(1,iQ,K)*profk)
                rha2 = rha2 + dimag(TB(2,iQ,K)*profk)
                rha3 = rha3 + dimag(TB(3,iQ,K)*profk)
!$omp end workshare

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

            ! 6-j
            f63 = fun6j(rLu,rLl,1d0,rJl,rJu1,S,Flgsg)

            if (abs(f63).lt.TINYJS) cycle

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

                ! 6-j
                f64 = fun6j(1d0,1d0,rK,rJu,rJu1,rJl,Flgsg)

!$omp workshare
                profk = f64*f63*prof*Atom%crho(iR,iz)
!$omp end workshare

!$omp workshare
                ! Absorptivity
                eps0 = eps0 + dble(TB(0,iQ,K)*profK)
                eps1 = eps1 + dble(TB(1,iQ,K)*profK)
                eps2 = eps2 + dble(TB(2,iQ,K)*profK)
                eps3 = eps3 + dble(TB(3,iQ,K)*profK)

                ! Dispersion
                rhs1 = rhs1 + dimag(TB(1,iQ,K)*profK)
                rhs2 = rhs2 + dimag(TB(2,iQ,K)*profK)
                rhs3 = rhs3 + dimag(TB(3,iQ,K)*profK)
!$omp end workshare

              end do ! Q
            end do ! K
          end do ! iU1

          ! Return to common loop

        end do ! iU
      end do ! iL

      ! Common parts for coefficients
      tempRe = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
               Atom%Ecoeff(itermu,iterml)/Dw
      tempRa = tempRe/absK

      ! Final values
!$omp workshare
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
!$omp end workshare

!$omp end parallel

      end subroutine rt1ordNB

#ifdef RDIPEV
!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption and emission coefficients.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!                            profile\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity\n
      !!           eps0(dfloat(:)): Intensity emissivity\n
      !!           eps1(dfloat(:)): Q emissivity\n
      !!           eps2(dfloat(:)): U emissivity\n
      !!           eps3(dfloat(:)): V emissivity\n
      !!           rhs1(dfloat(:)): Q 'dichroic' emissivity\n
      !!           rhs2(dfloat(:)): U 'dichroic' emissivity\n
      !!           rhs3(dfloat(:)): V 'dichroic' emissivity
      subroutine rt1ord(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                        iz,if0,if1,Norma,Dw,vfac,absK,aprof, &
                        eta0,eta1,eta2,eta3,rha1,rha2,rha3, &
                        eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, absK, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      logical:: lNCHLT

      integer:: ifreq,K,iq,iq1,iQQ
      integer:: nMu,nMl,iMu,iMl,iMu1,iMl1
      integer:: iJl1,iJlb,iJu1,iJub
      integer:: iU,iU1,kU,kU1,kUb,iL,iL1,kL,kL1,kLb

      double precision:: rLu,rLl,S,rJu,rJu1,rJub,rJl,rJl1,rJlb
      double precision:: rJumax,rJlmax,rMu,rMu1,rMl,rMl1
      double precision:: eu,el,rK,QQ,q,q1,au,al,aul,ftmp
      double precision:: Cu1,Cub,Cl1,Clb
      double precision:: at,Dfreq,vfacw,tempRe,tempRa
      double precision:: EVul,EVu1l,EVul1

      complex(kind=8):: tK,rhoc
      complex(kind=8),dimension(if0:if1):: prof, profK


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

      if (NCHLT.and.allocated(Atom%NCHLT)) then
        lNCHLT = Atom%NCHLT(iz,itran)
      else
        lNCHLT = .False.
      end if

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

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
      vfacw = vfac/Dw

!$omp parallel default(none) &
!$omp private(iMu,rMu,iU,eu,iMl,rMl,q,iq,iL,el,Dfreq,ifreq,iMu1) &
!$omp private(rMu1,iMl1,rMl1,q1,QQ,iq1,iQQ,K,rK,ftmp,tK,kU) &
!$omp private(iJu,rJu,kL,iJl,rJl,kU1,Cu1,Cl1,iJu1,rJu1,kL1) &
!$omp private(iJl1,rJl1,kUb,kLb,Cub,Clb,iJub,iJlb,rJub,rJlb) &
!$omp private(rhoc,tempR,EVul,EVu1l,EVul1) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3,at,rLu,rLl,vfacw) &
!$omp shared(eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3,lNCHLT) &
!$omp shared(nMu,nMl,rJumax,rJlmax,Atom,vpfil,prof,aprof,Norma) &
!$omp shared(iterml,itermu,profk,iz,Dw,itran,if0,if1) &
!$omp shared(omega,Flgsg,TB)

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
          el = Atom%eval(iL,iMl,iterml,iz)/Dw

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! If not pi nor sigma, skip
            if (abs(rMu-rMl).gt.1.1d0) cycle

            ! Get difference between M momentums in integer
            q = rMu-rMl
            iq = nint(q)

            ! For each mu_u
            do iU=1,Atom%nblk(iMu,itermu) ! sum over mu_u

              ! Get eigenvalue upper level
              eu = Atom%eval(iU,iMu,itermu,iz)/Dw

              ! Dipole strength
              EVul = Atom%rdipev(iz)%rdipev(itran)% &
                          rdip(iq,iU,iL,iMu,iMl)

              ! Check if small
              if (abs(EVul).lt.TINYEV) cycle

              !
              ! Compute profile
              !

              ! If in file
              if (vpfil) then
!$omp workshare
                prof = aprof(:,Atom%i_Vind(itran)%ind(iL,iMl,iU,iMu))
!$omp end workshare
              ! If stored
              else if (Norma%VRAM) then
!$omp workshare
                prof = Norma%prof(iL,iU,iMl,iMu)%cp
!$omp end workshare
              ! Not stored
              else

                ! Shift term
                Dfreq = eu - el + Atom%Dfreq(itran)/Dw

                ! For each frequency
!$omp do
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies
!$omp end do
                ! Normalize
!$omp workshare
                prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,iMl,iMu), &
                              dimag(prof))
!$omp end workshare

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

                  ! If not allowed (3j-sym=0) skip
                  if (abs(ftmp).lt.TINYJS) cycle

                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)*EVul

      ! Reset identation
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

            ! Sign
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

              ! If no population, skip
              if (abs(rhoc).lt.TINYER) cycle

              ! Accumulate into tK
              tK = ftmp*EVu1l*Cu1*Cub*rhoc + tK

             end do ! kUb
           end do ! kU1
         end do ! iU1
       end do ! iMu1

       ! Add the profile
!$omp workshare
       profK = prof*tK
!$omp end workshare

!$omp workshare
       ! Emissivity
       eps0 = eps0 + dble(TB(0,iQQ,K)*profK)
       eps1 = eps1 + dble(TB(1,iQQ,K)*profK)
       eps2 = eps2 + dble(TB(2,iQQ,K)*profK)
       eps3 = eps3 + dble(TB(3,iQQ,K)*profK)

       ! Dispersion
       rhs1 = rhs1 + dimag(TB(1,iQQ,K)*profK)
       rhs2 = rhs2 + dimag(TB(2,iQQ,K)*profK)
       rhs3 = rhs3 + dimag(TB(3,iQQ,K)*profK)
!$omp end workshare

      !
      ! Absorption
      !

      ! Initialize tK
      tK = cZero

      ! For each Ml'
      do iMl1=1,nMl

        ! NCHLT
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
!$omp workshare
      profK = prof*tk
!$omp end workshare

!$omp workshare
      ! Absorptivity
      eta0 = eta0 + dble(TB(0,iQQ,K)*profk)
      eta1 = eta1 + dble(TB(1,iQQ,K)*profk)
      eta2 = eta2 + dble(TB(2,iQQ,K)*profk)
      eta3 = eta3 + dble(TB(3,iQQ,K)*profk)

      ! Dispersion
      rha1 = rha1 + dimag(TB(1,iQQ,K)*profk)
      rha2 = rha2 + dimag(TB(2,iQQ,K)*profk)
      rha3 = rha3 + dimag(TB(3,iQQ,K)*profk)
!$omp end workshare

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
!$omp workshare
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
!$omp end workshare
!$omp end parallel

      end subroutine rt1ord

#endif
!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!                            profile\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity
      subroutine absorb(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                        iz,if0,if1,Norma,Dw,vfac,absK,aprof, &
                        eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, absK, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      logical:: lNCHLT

      integer:: ifreq,K,iq,iq1,iQQ
      integer:: nMu,nMl,iMu,iMl,iMl1,iJu,iJu1,iJl,iJl1,iJlb
      integer:: iU,kU,kU1,iL,kL,kL1,kLb

      double precision:: rLu,rLl,S,rJu,rJu1,rJl,rJl1,rJlb
      double precision:: rJumax,rJlmax,rMu,rMl,rMl1
      double precision:: eu,el,rK,QQ,q,q1,au,al,aul,ftmp,tempR
      double precision:: Cu,Cu1,Cl,Clb,CC,CC1
      double precision:: at,Dfreq,vfacw

      complex(kind=8):: tK,rhoc
      complex(kind=8),dimension(if0:if1):: prof, profK


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

      if (NCHLT.and.allocated(Atom%NCHLT)) then
        lNCHLT = Atom%NCHLT(iz,itran)
      else
        lNCHLT = .False.
      end if

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

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
      vfacw = vfac/Dw

!$omp parallel default(none) &
!$omp private(iMu,rMu,iU,eu,iMl,rMl,q,iq,iL,el,Dfreq,ifreq,iMl1) &
!$omp private(rMl1,q1,QQ,iq1,iQQ,K,rK,ftmp,tK,kU,Cu,iJu,rJu) &
!$omp private(kL,Cl,iJl,rJl,CC,kU1,Cu1,iJu1,rJu1,kL1) &
!$omp private(iJl1,rJl1,CC1,kLb,Clb,iJlb,rJlb,rhoc,tempR) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3,at,rLu,rLl,vfacw) &
!$omp shared(nMu,nMl,rJumax,rJlmax,Atom,vpfil,prof,aprof,Norma) &
!$omp shared(lNCHLT,iterml,itermu,profk,iz,Dw,itran,if0,if1) &
!$omp shared(omega,Flgsg,TB,absK)

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
          eu = Atom%eval(iU,iMu,itermu,iz)/Dw

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
              el = Atom%eval(iL,iMl,iterml,iz)/Dw

              !
              ! Compute profile
              !

              ! If in file
              if (vpfil) then
!$omp workshare
                prof = aprof(:,Atom%i_Vind(itran)%ind(iL,iMl,iU,iMu))
!$omp end workshare
              ! If stored
              else if (Norma%VRAM) then
!$omp workshare
                prof = Norma%prof(iL,iU,iMl,iMu)%cp
!$omp end workshare
              ! Not stored
              else

                ! Shift term
                Dfreq = eu - el + Atom%Dfreq(itran)/Dw

                ! For each frequency
!$omp do
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies
!$omp end do
                ! Normalize profile
!$omp workshare
                prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,iMl,iMu), &
                              dimag(prof))
!$omp end workshare

              end if ! Storing

              ! For each Ml'
              do iMl1=1,nMl

                ! NCHLT
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

                  ! If not allowed (3j-sym=0) skip
                  if (abs(ftmp).lt.TINYJS) cycle

                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)

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
                      CC = Cl*Cu*Atom%rdip(itran)% &
                                      rdip(iq,iMu,iMl,iJu,iJl)

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
                          CC1 = Cu1*Atom%rdip(itran)% &
                                         rdip(iq1,iMu,iMl1,iJu1,iJl1)

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
                                         iJl1,rJlb,rJl1,rMl,rMl1,0, &
                                         rhoc)

                            ! If no population, skip
                            if (abs(rhoc).lt.TINYER) cycle

                            ! Accumulate into tK
                            tK = ftmp*CC*CC1*Clb*rhoc + tK

                          end do ! kLb
                        end do ! kL1
                      end do ! kU1
                    end do ! kL
                  end do ! kU

                  ! Add the profile
!$omp workshare
                  profK = prof*tk
!$omp end workshare

                  ! For each Stokes parameter

!$omp workshare
                  ! Absorptivity
                  eta0 = eta0 + dble(TB(0,iQQ,K)*profk)
                  eta1 = eta1 + dble(TB(1,iQQ,K)*profk)
                  eta2 = eta2 + dble(TB(2,iQQ,K)*profk)
                  eta3 = eta3 + dble(TB(3,iQQ,K)*profk)

                  ! Dispersion
                  rha1 = rha1 + dimag(TB(1,iQQ,K)*profk)
                  rha2 = rha2 + dimag(TB(2,iQQ,K)*profk)
                  rha3 = rha3 + dimag(TB(3,iQQ,K)*profk)
!$omp end workshare

                end do ! K
              end do ! Ml1
            end do ! iL
          end do ! Ml
        end do ! iU
      end do ! Mu

      ! Common part for the two coefficients
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)/(Dw*absK)

      ! Final values
!$omp workshare
      eta0 = tempR*eta0
      eta1 = tempR*eta1
      eta2 = tempR*eta2
      eta3 = tempR*eta3
      rha1 = tempR*rha1
      rha2 = tempR*rha2
      rha3 = tempR*rha3
!$omp end workshare
!$omp end parallel

      return

      end subroutine absorb

#ifdef _OPENMP
!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient, prepared to split in
      !! components with OpenMP.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!                            profile\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!      omp_1c(omp_1c_class): Indexes for OpenMP\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity
      subroutine absorb_c(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                          iz,if0,if1,Norma,Dw,vfac,absK,aprof,ompc1, &
                          eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      type(omp_1c_class), intent(in):: ompc1
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, absK, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      logical:: lNCHLT

      integer:: ifreq,K,iq,iq1,iQQ
      integer:: nMu,nMl,iMu,iMl,iMl1,iJu,iJu1,iJl,iJl1,iJlb
      integer:: iU,kU,kU1,iL,kL,kL1,kLb,tid,icom

      double precision:: rLu,rLl,S,rJu,rJu1,rJl,rJl1,rJlb
      double precision:: rJumax,rJlmax,rMu,rMl,rMl1
      double precision:: eu,el,rK,QQ,q,q1,au,al,aul,ftmp,tempR
      double precision:: Cu,Cu1,Cl,Clb,CC,CC1
      double precision:: at,Dfreq,vfacw

      complex(kind=8):: tK,rhoc
      complex(kind=8),dimension(if0:if1):: prof, profK


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

      if (NCHLT.and.allocated(Atom%NCHLT)) then
        lNCHLT = Atom%NCHLT(iz,itran)
      else
        lNCHLT = .False.
      end if

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

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
      vfacw = vfac/Dw

!$omp parallel default(None) &
!$omp private(icom,tid,iMu,rMu,iU,eu,iMl,rMl,q,iq,iL,el,prof,Dfreq) &
!$omp private(ifreq,iMl1,rMl1,q1,QQ,iq1,iQQ,K,rK,ftmp,tK,Ku,Cu,iJu) &
!$omp private(rJu,kL,Cl,iJl,rJl,CC,kU1,Cu1,iJu1,rJu1) &
!$omp private(kL1,iJl1,rJl1,CC1,kLb,Clb,iJlb,rJlb,rhoc) &
!$omp private(profK) &
!$omp shared(lNCHLT,at,rLu,rLl,nMu,nMl,rJumax,rJlmax,vfacw,Atom) &
!$omp shared(aprof,Norma,vpfil,if0,if1,Flgsg,iz,iterml,TB,ompc1) &
!$omp shared(itermu,Dw,itran,omega) &
!$omp shared(pid) &
!$omp reduction(+: eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! Initialize component count
      icom = 0

      ! Get thread index
      tid = omp_get_thread_num() + 1


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
          eu = Atom%eval(iU,iMu,itermu,iz)/Dw

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

              ! Advance index and check
              icom = icom + 1

              ! If smaller, skip
              if (icom.lt.ompc1%if0(tid)) cycle
              ! If larger, finish
              if (icom.gt.ompc1%if1(tid)) exit

              ! Get eigenvalue lower level
              el = Atom%eval(iL,iMl,iterml,iz)/Dw

              !
              ! Compute profile
              !

              ! If in file
              if (vpfil) then

                prof = aprof(:,Atom%i_Vind(itran)%ind(iL,iMl,iU,iMu))

              ! If stored
              else if (Norma%VRAM) then

                prof = Norma%prof(iL,iU,iMl,iMu)%cp

              ! Not stored
              else

                ! Shift term
                Dfreq = eu - el + Atom%Dfreq(itran)/Dw

                ! For each frequency
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies

                ! Normalize profile
                prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,iMl,iMu), &
                              dimag(prof))

              end if ! Storing

              ! For each Ml'
              do iMl1=1,nMl

                ! NCHLT
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

                  ! If not allowed (3j-sym=0) skip
                  if (abs(ftmp).lt.TINYJS) cycle

                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)

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
                      CC = Cl*Cu*Atom%rdip(itran)% &
                                      rdip(iq,iMu,iMl,iJu,iJl)

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
                          CC1 = Cu1*Atom%rdip(itran)% &
                                         rdip(iq1,iMu,iMl1,iJu1,iJl1)

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

                            ! Sum over (Ku,Qu)
                            call absorb1(Atom,Flgsg,iz,iterml,iJlb, &
                                         iJl1,rJlb,rJl1,rMl,rMl1,0, &
                                         rhoc)

                            ! If no population, skip
                            if (abs(rhoc).lt.TINYER) cycle

                            ! Accumulate into tK
                            tK = ftmp*CC*CC1*Clb*rhoc + tK

                          end do ! kLb
                        end do ! kL1
                      end do ! kU1
                    end do ! kL
                  end do ! kU

                  ! Add the profile
                  profK = prof*tk

                  ! For each Stokes parameter

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
!$omp end parallel

      ! Common part for the two coefficients
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)/(Dw*absK)

      ! Final values
!$omp parallel workshare default(none) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3,tempR)
      eta0 = tempR*eta0
      eta1 = tempR*eta1
      eta2 = tempR*eta2
      eta3 = tempR*eta3
      rha1 = tempR*rha1
      rha2 = tempR*rha2
      rha3 = tempR*rha3
!$omp end parallel workshare

      return

      end subroutine absorb_c
#endif

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient in absence of magnetic
      !! fields.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity
      subroutine absorbNB(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                          iz,if0,if1,Norma,Dw,vfac,absK,aprof, &
                          eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, absK, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0
      double precision, dimension(if0:if1), intent(out):: eta1,rha1
      double precision, dimension(if0:if1), intent(out):: eta2,rha2
      double precision, dimension(if0:if1), intent(out):: eta3,rha3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      integer:: ifreq,K,iQ,iU,iL,iL1,iR
      integer:: Kmin,Kmax

      double precision:: rLu,rLl,S,rJu,rJl,rJl1
      double precision:: f61,f62,f63,f64
      double precision:: eu,el,rK,au,al,aul
      double precision:: at,Dfreq,vfacw,tempR

      complex(kind=8),dimension(if0:if1):: prof, profK


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

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

!$omp parallel default(none) &
!$omp private(iU,iL,Dfreq,ifreq,iL1,K,iQ,rJu,eu,f61,rJl,f62,el) &
!$omp private(rJl1,f63,Kmin,Kmax,rK,iR,f64,tempR) &
!$omp shared(Atom,Flgsg,Dw,at,S,rLu,rLl,vfacw,vpfil,prof,itran) &
!$omp shared(profk,eta0,eta1,eta2,eta3,TB,Norma,rha1,rha2,rha3) &
!$omp shared(itermu,iterml,aprof,if0,if1,omega,iz,absK)


      !
      ! Compute absorptivity
      !

      ! For each Ju
      do iU=1,Atom%nJ(itermu)

        ! Get Ju
        rJu = Atom%rJval(iU,itermu)

        ! Get eigenvalue upper level
        eu = Atom%FSfreq(iU,itermu)/Dw

        ! Factor Ju
        f61 = 2d0*rJu + 1d0

        ! For each Jl
        do iL=1,Atom%nJ(iterml)

          ! Get Jl
          rJl = Atom%rJval(iL,iterml)

          ! 6-j
          f62 = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

          if (abs(f62).lt.TINYJS) cycle

          f62 = f62*f61*Flgsg%sg(nint(1d0+rJu+rJl))*sqrt(2d0*rJl+1d0)

          ! Get eigenvalue lower level
          el = Atom%FSfreq(iL,iterml)/Dw

          !
          ! Compute profile
          !

          ! If in file
          if (vpfil) then

!$omp workshare
            prof = aprof(:,Atom%i_Vind(itran)%indNB(iL,iU))
!$omp end workshare

          ! If stored
          else if (Norma%VRAM) then

!$omp workshare
            prof = Norma%prof(iL,iU,1,1)%cp
!$omp end workshare

          ! Not stored
          else

            ! Shift term
            Dfreq = eu - el

            ! For each frequency
!$omp do
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreq - omega(ifreq)*vfacw,at,prof(ifreq))

            end do ! frequencies
!$omp end do

            ! Normalize profile
!$omp workshare
            prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,1,1), &
                          dimag(prof))
!$omp end workshare

          end if ! Storing

          ! For each Jl'
          do iL1=1,Atom%nJ(iterml)

            ! Get Jl
            rJl1 = Atom%rJval(iL1,iterml)

            ! 6-j
            f63 = fun6j(rLu,rLl,1d0,rJl1,rJu,S,Flgsg)

            if (abs(f63).lt.TINYJS) cycle

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

!$omp workshare
                profk = f64*f63*Flgsg%sg(K)*prof*Atom%crho(iR,iz)
!$omp end workshare

!$omp workshare
                ! Absorptivity
                eta0 = eta0 + dble(TB(0,iQ,K)*profk)
                eta1 = eta1 + dble(TB(1,iQ,K)*profk)
                eta2 = eta2 + dble(TB(2,iQ,K)*profk)
                eta3 = eta3 + dble(TB(3,iQ,K)*profk)

                ! Dispersion
                rha1 = rha1 + dimag(TB(1,iQ,K)*profk)
                rha2 = rha2 + dimag(TB(2,iQ,K)*profk)
                rha3 = rha3 + dimag(TB(3,iQ,K)*profk)
!$omp end workshare

              end do ! Q
            end do ! K
          end do ! iL1
        end do ! iL
      end do ! iU

      ! Common part for the two coefficients
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)/(Dw*absK)

      ! Final values
!$omp workshare
      eta0 = tempR*eta0
      eta1 = tempR*eta1
      eta2 = tempR*eta2
      eta3 = tempR*eta3
      rha1 = tempR*rha1
      rha2 = tempR*rha2
      rha3 = tempR*rha3
!$omp end workshare

!$omp end parallel

      return

      end subroutine absorbNB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the inner summation on rhoKQ(l)\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!            iterm(integer): Term to sum over\n
      !!              iJ1(integer): J index of the level to sum over\n
      !!              iJ2(integer): J' index of the level to sum
      !!                            over\n
      !!              rJ1(integer): J of the level to sum over\n
      !!              rJ2(integer): J' of the level to sum over\n
      !!              rM1(integer): M of the level to sum over\n
      !!              rM2(integer): M' of the level to sum over\n
      !!              Kin(integer): K of the inner part of emiss2ord
      !!              summ(dcmplx): Summation result
      subroutine absorb1(Atom,Flgsg,iz,iterm,iJ1,iJ2,rJ1,rJ2,rM1, &
                         rM2,Kin,summ)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz, iterm, iJ1, iJ2, Kin
      double precision, intent(in):: rJ1, rJ2, rM1, rM2
      complex(kind=8), intent(out):: summ

      ! Local

      integer:: K,iQ,Kmin,Kmax,iR

      double precision:: rK,Q,fKf3j


      ! Difference between magnetic momentums
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

      !> Computes the emission coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!           eps0(dfloat(:)): Intensity emissivity\n
      !!           eps1(dfloat(:)): Q emissivity\n
      !!           eps2(dfloat(:)): U emissivity\n
      !!           eps3(dfloat(:)): V emissivity\n
      !!           rhs1(dfloat(:)): Q 'dichroic' emissivity\n
      !!           rhs2(dfloat(:)): U 'dichroic' emissivity\n
      !!           rhs3(dfloat(:)): V 'dichroic' emissivity
      subroutine emiss(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                       iz,if0,if1,Norma,Dw,vfac,aprof, &
                       eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      integer:: ifreq,K,iq,iq1,iQQ
      integer:: nMl,nMu,iMl,iMu,iMu1,iJl,iJl1,iJu,iJu1,iJub
      integer:: iL,kL,kL1,iU,kU,kU1,kUb

      double precision:: rLl,rLu,S,rJl,rJl1,rJu,rJu1,rJub
      double precision:: rJlmax,rJumax,rMl,rMu,rMu1,el,eu
      double precision:: rK,QQ,q,q1,al,au,aul,ftmp,tempR
      double precision:: Cl,Cl1,Cu,Cub,CC,CC1
      double precision:: at,Dfreq,vfacw

      complex(kind=8):: tK,rhoc
      complex(kind=8),dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

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

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

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
      vfacw = vfac/Dw

!$omp parallel default(none) &
!$omp private(iMu,rMu,iU,eu,iMl,rMl,q,iq,iL,el,Dfreq,ifreq,iMu1) &
!$omp private(rMu1,q1,QQ,iq1,iQQ,K,rK,ftmp,tK,kU,Cu,iJu,rJu) &
!$omp private(kL,Cl,iJl,rJl,CC,kU1,Cl1,iJu1,rJu1,kL1) &
!$omp private(iJl1,rJl1,CC1,kUb,Cub,Ijub,rJub,rhoc,tempR) &
!$omp shared(eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3,at,rLu,rLl,vfacw) &
!$omp shared(nMu,nMl,rJumax,rJlmax,Atom,vpfil,prof,aprof,Norma) &
!$omp shared(iterml,itermu,profk,iz,Dw,itran,if0,if1) &
!$omp shared(omega,Flgsg,TB)


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

              !
              ! Compute profile
              !

              ! If in file
              if (vpfil) then
!$omp workshare
                prof = aprof(:,Atom%i_Vind(itran)%ind(iL,iMl,iU,iMu))
!$omp end workshare
              ! If stored
              else if (Norma%VRAM) then
!$omp workshare
                prof = Norma%prof(iL,iU,iMl,iMu)%cp
!$omp end workshare
              ! Not stored
              else

                ! Shift term
                Dfreq = eu - el + Atom%Dfreq(itran)/Dw

                ! For each frequency
!$omp do
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies
!$omp end do
                ! Normalize
!$omp workshare
                prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,iMl,iMu), &
                              dimag(prof))
!$omp end workshare

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

                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)

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
                      CC = Cu*Cl*Atom%rdip(itran)% &
                                      rdip(iq,iMu,iMl,iJu,iJl)

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
                                Atom%rdip(itran)% &
                                     rdip(iq1,iMu1,iMl,iJu1,iJl1)

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

                  ! Add the profile
!$omp workshare
                  profK = prof*tK
!$omp end workshare

!$omp workshare
                  ! Emissivity
                  eps0 = eps0 + dble(TB(0,iQQ,K)*profK)
                  eps1 = eps1 + dble(TB(1,iQQ,K)*profK)
                  eps2 = eps2 + dble(TB(2,iQQ,K)*profK)
                  eps3 = eps3 + dble(TB(3,iQQ,K)*profK)

                  ! Dispersion
                  rhs1 = rhs1 + dimag(TB(1,iQQ,K)*profK)
                  rhs2 = rhs2 + dimag(TB(2,iQQ,K)*profK)
                  rhs3 = rhs3 + dimag(TB(3,iQQ,K)*profK)
!$omp end workshare

                end do ! K
              end do ! Mu1
            end do ! iU
          end do ! Mu
        end do ! iL
      end do ! Ml

      ! Common part for the two coefficients
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)/Dw

      ! Final values
!$omp workshare
      eps0 = tempR*eps0
      eps1 = tempR*eps1
      eps2 = tempR*eps2
      eps3 = tempR*eps3
      rhs1 = tempR*rhs1
      rhs2 = tempR*rhs2
      rhs3 = tempR*rhs3
!$omp end workshare
!$omp end parallel

      return

      end subroutine emiss

#ifdef _OPENMP
!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emission coefficient, prepared to split in
      !! components with OpenMP.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!      omp_1c(omp_1c_class): Indexes for OpenMP\n
      !!           eps0(dfloat(:)): Intensity emissivity\n
      !!           eps1(dfloat(:)): Q emissivity\n
      !!           eps2(dfloat(:)): U emissivity\n
      !!           eps3(dfloat(:)): V emissivity\n
      !!           rhs1(dfloat(:)): Q 'dichroic' emissivity\n
      !!           rhs2(dfloat(:)): U 'dichroic' emissivity\n
      !!           rhs3(dfloat(:)): V 'dichroic' emissivity
      subroutine emiss_c(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                         iz,if0,if1,Norma,Dw,vfac,aprof,ompc1, &
                         eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      type(omp_1c_class), intent(in):: ompc1
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      integer:: ifreq,K,iq,iq1,iQQ
      integer:: nMl,nMu,iMl,iMu,iMu1,iJl,iJl1,iJu,iJu1,iJub
      integer:: iL,kL,kL1,iU,kU,kU1,kUb,tid,icom

      double precision:: rLl,rLu,S,rJl,rJl1,rJu,rJu1,rJub
      double precision:: rJlmax,rJumax,rMl,rMu,rMu1,el,eu
      double precision:: rK,QQ,q,q1,al,au,aul,ftmp,tempR
      double precision:: Cl,Cl1,Cu,Cub,CC,CC1
      double precision:: at,Dfreq,vfacw

      complex(kind=8):: tK,rhoc
      complex(kind=8),dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

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

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

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
      vfacw = vfac/Dw

!$omp parallel default(None) &
!$omp private(icom,tid,iMl,rMl,iL,el,iMu,rMu,q,iq,iU,eu,prof,Dfreq) &
!$omp private(ifreq,iMu1,rMu1,q1,QQ,iq1,iQQ,K,rK,ftmp,tK,Kl,Cl,iJl) &
!$omp private(rJl,kU,Cu,iJu,rJu,CC,kL1,Cl1,iJl1,rJl1) &
!$omp private(kU1,iJu1,rJu1,CC1,kUb,Cub,iJub,rJub,rhoc) &
!$omp private(profK) &
!$omp shared(NCHLT,at,rLu,rLl,nMu,nMl,rJumax,rJlmax,vfacw,Atom) &
!$omp shared(aprof,Norma,vpfil,if0,if1,Flgsg,iz,itermu,TB,ompc1) &
!$omp shared(iterml,Dw,itran,omega) &
!$omp shared(pid) &
!$omp reduction(+: eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! Initialize component count
      icom = 0

      ! Get thread index
      tid = omp_get_thread_num() + 1

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

              ! Advance index and check
              icom = icom + 1

              ! If smaller, skip
              if (icom.lt.ompc1%if0(tid)) cycle
              ! If larger, finish
              if (icom.gt.ompc1%if1(tid)) exit

              ! Get eigenvalue upper level
              eu = Atom%eval(iU,iMu,itermu,iz)/Dw

              !
              ! Compute profile
              !

              ! If in file
              if (vpfil) then

                prof = aprof(:,Atom%i_Vind(itran)%ind(iL,iMl,iU,iMu))

              ! If stored
              else if (Norma%VRAM) then

                prof = Norma%prof(iL,iU,iMl,iMu)%cp

              ! Not stored
              else

                ! Shift term
                Dfreq = eu - el + Atom%Dfreq(itran)/Dw

                ! For each frequency
                do ifreq=if0,if1

                  ! Calculate profile
                  call voigt(Dfreq - omega(ifreq)*vfacw,at, &
                             prof(ifreq))

                end do ! frequencies

                ! Normalize
                prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,iMl,iMu), &
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

                  ftmp = ftmp*Flgsg%sg(iq1+1)*sqrt(2d0*rK+1d0)

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
                      CC = Cu*Cl*Atom%rdip(itran)% &
                                      rdip(iq,iMu,iMl,iJu,iJl)

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
                                Atom%rdip(itran)% &
                                     rdip(iq1,iMu1,iMl,iJu1,iJl1)

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

                  ! Add the profile
                  profK = prof*tK

                  ! Absorptivity
                  eps0 = eps0 + dble(TB(0,iQQ,K)*profK)
                  eps1 = eps1 + dble(TB(1,iQQ,K)*profK)
                  eps2 = eps2 + dble(TB(2,iQQ,K)*profK)
                  eps3 = eps3 + dble(TB(3,iQQ,K)*profK)

                  ! Dispersion
                  rhs1 = rhs1 + dimag(TB(1,iQQ,K)*profK)
                  rhs2 = rhs2 + dimag(TB(2,iQQ,K)*profK)
                  rhs3 = rhs3 + dimag(TB(3,iQQ,K)*profK)

                end do ! K
              end do ! Mu1
            end do ! iU
          end do ! Mu
        end do ! iL
      end do ! Ml
!$omp end parallel

      ! Common part for the two coefficients
      tempR = 1d8*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)*(1d-5/Dw)

      ! Final values
!$omp parallel workshare default(none) &
!$omp shared(eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3,tempR)
      eps0 = tempR*eps0
      eps1 = tempR*eps1
      eps2 = tempR*eps2
      eps3 = tempR*eps3
      rhs1 = tempR*rhs1
      rhs2 = tempR*rhs2
      rhs3 = tempR*rhs3
!$omp end parallel workshare

      return

      end subroutine emiss_c
#endif

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emission coefficient in absence of magnetic
      !! fields.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            itran(integer): Index of transition to compute\n
      !!           itermu(integer): Upper term of transition\n
      !!           iterml(integer): Lower term of transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!              absK(dfloat): Unit transformation factor\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!           eps0(dfloat(:)): Intensity emissivity\n
      !!           eps1(dfloat(:)): Q emissivity\n
      !!           eps2(dfloat(:)): U emissivity\n
      !!           eps3(dfloat(:)): V emissivity\n
      !!           rhs1(dfloat(:)): Q 'dichroic' emissivity\n
      !!           rhs2(dfloat(:)): U 'dichroic' emissivity\n
      !!           rhs3(dfloat(:)): V 'dichroic' emissivity
      subroutine emissNB(Atom,TB,omega,Flgsg,itran,itermu,iterml, &
                         iz,if0,if1,Norma,Dw,vfac,aprof, &
                         eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Nindex_class), intent(in):: Norma
      integer, intent(in):: itran, itermu, iterml, iz, if0, if1
      double precision, intent(in):: Dw, vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps0
      double precision, dimension(if0:if1), intent(out):: eps1,rhs1
      double precision, dimension(if0:if1), intent(out):: eps2,rhs2
      double precision, dimension(if0:if1), intent(out):: eps3,rhs3
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      integer:: ifreq,K,iQ,iL,iU,iU1,iR
      integer:: Kmin,Kmax

      double precision:: rLl,rLu,S,rJl,rJu,rJu1,rK
      double precision:: f61,f62,f63,f64
      double precision:: el,eu,al,au,aul,tempR
      double precision:: at,Dfreq,vfacw

      complex(kind=8),dimension(if0:if1):: prof, profK


      !
      ! Initialize variables
      !

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

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      al = Atom%damp(iterml,iz)
      aul = Atom%ldamp(itran,iz)
      at = (au + al + aul)/Dw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLl = Atom%rLval(iterml)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

!$omp parallel default(none) &
!$omp private(iL,el,rJl,f61,iU,rJu,f62,eu,Dfreq,ifreq,iU1,rJu1,f63) &
!$omp private(Kmin,Kmax,K,rK,iQ,iR,f64,tempR) &
!$omp shared(Atom,iterml,itermu,vpfil,prof,aprof,itran,Norma,if0) &
!$omp shared(if1,omega,vfacw,at,S,rLu,rLl,eps0,eps1,eps2,eps3,TB) &
!$omp shared(profk,rhs1,rhs2,rhs3,Dw,Flgsg,iz)


      !
      ! Compute emissivity
      !

      ! For each Jl
      do iL=1,Atom%nJ(iterml)

        ! Get eigenvalue lower level
        el = Atom%FSfreq(iL,iterml)/Dw

        ! Jl
        rJl = Atom%rJval(iL,iterml)

        ! Jl factor
        f61 = 2d0*rJl+1d0

        ! For each Ju
        do iU=1,Atom%nJ(itermu)

          ! Ju
          rJu = Atom%rJval(iU,itermu)

          ! 6-j
          f62 = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

          if (abs(f62).lt.TINYJS) cycle

          f62 = f62*f61*sqrt(2d0*rJu+1d0)

          ! Get eigenvalue upper level
          eu = Atom%FSfreq(iU,itermu)/Dw

          !
          ! Compute profile
          !

          ! If in file
          if (vpfil) then

!$omp workshare
            prof = aprof(:,Atom%i_Vind(itran)%indNB(iL,iU))
!$omp end workshare

          ! If stored
          else if (Norma%VRAM) then

!$omp workshare
            prof = Norma%prof(iL,iU,1,1)%cp
!$omp end workshare

          ! Not stored
          else

            ! Shift term
            Dfreq = eu - el

            ! For each frequency
!$omp do
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreq - omega(ifreq)*vfacw,at,prof(ifreq))

            end do ! frequencies
!$omp end do

            ! Normalize
!$omp workshare
            prof = dcmplx(dble(prof)*Norma%Norm(iL,iU,1,1), &
                          dimag(prof))
!$omp end workshare

          end if ! Storing

          ! For each Ju'
          do iU1=1,Atom%nJ(itermu)

            ! Ju'
            rJu1 = Atom%rJval(iU1,itermu)

            ! 6-j
            f63 = fun6j(rLu,rLl,1d0,rJl,rJu1,S,Flgsg)

            if (abs(f63).lt.TINYJS) cycle

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

                ! 6-j
                f64 = fun6j(1d0,1d0,rK,rJu,rJu1,rJl,Flgsg)

!$omp workshare
                profk = f64*f63*prof*Atom%crho(iR,iz)
!$omp end workshare

!$omp workshare
                ! Absorptivity
                eps0 = eps0 + dble(TB(0,iQ,K)*profK)
                eps1 = eps1 + dble(TB(1,iQ,K)*profK)
                eps2 = eps2 + dble(TB(2,iQ,K)*profK)
                eps3 = eps3 + dble(TB(3,iQ,K)*profK)

                ! Dispersion
                rhs1 = rhs1 + dimag(TB(1,iQ,K)*profK)
                rhs2 = rhs2 + dimag(TB(2,iQ,K)*profK)
                rhs3 = rhs3 + dimag(TB(3,iQ,K)*profK)
!$omp end workshare

              end do ! Q
            end do ! K
          end do ! iU1
        end do ! iU
      end do ! iL

      ! Common part for the two coefficients
      tempR = 1d3*sqrt3*IPI41*(2d0*rLu+1d0)* &
              Atom%Ecoeff(itermu,iterml)/Dw

      ! Final values
!$omp workshare
      eps0 = tempR*eps0
      eps1 = tempR*eps1
      eps2 = tempR*eps2
      eps3 = tempR*eps3
      rhs1 = tempR*rhs1
      rhs2 = tempR*rhs2
      rhs3 = tempR*rhs3
!$omp end workshare

!$omp end parallel

      return

      end subroutine emissNB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the inner summation on rhoKQ(u)\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!            iterm(integer): Term to sum over\n
      !!              iJ1(integer): J index of the level to sum
      !!                            over\n
      !!              iJ2(integer): J' index of the level to sum
      !!                            over\n
      !!              rJ1(integer): J of the level to sum over\n
      !!              rJ2(integer): J' of the level to sum over\n
      !!              rM1(integer): M of the level to sum over\n
      !!              rM2(integer): M' of the level to sum over\n
      !!              summ(dcmplx): Summation result
      subroutine emiss1(Atom,Flgsg,iz,iterm,iJ1,iJ2,rJ1,rJ2,rM1, &
                        rM2,summ)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz, iterm, iJ1, iJ2
      double precision, intent(in):: rJ1, rJ2, rM1, rM2
      complex(kind=8), intent(out):: summ

      ! Local

      integer:: K,iQ,Kmin,Kmax,iR

      double precision:: rK,Q,fKf3j

      ! Difference between magnetic momentums
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

      !> Computes the second order emission coefficient.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!                vx(dfloat): Velocity vector along X\n
      !!                vy(dfloat): Velocity vector along Y\n
      !!                vz(dfloat): Velocity vector along Z\n
      !!          omega(dfloat): Frequency array\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!    Fin(Frequencyc2_class): Structure with the input frequency
      !!                            information\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!            jtran(integer): Output transition index\n
      !!           itermu(integer): Upper term of output transition\n
      !!           itermf(integer): Lower term of output transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!               DwT(dfloat): Thermal part of Doppler width\n
      !!                Dw(dfloat): Doppler width output transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!               vmi(dfloat): Microturbulent velocity\n
      !!      TBout(dcmplx(:,:,:)): Geometrical tensor in the output
      !!                            direction in the magnetic field
      !!                            reference frame\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!     JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!          eps20(dfloat(:)): Intensity emissivity\n
      !!          eps21(dfloat(:)): Q emissivity\n
      !!          eps22(dfloat(:)): U emissivity\n
      !!          eps23(dfloat(:)): V emissivity
      subroutine emiss2ord_AA(Atom,Geom,vx,vy,vz,omega,Red,Fin, &
                              Flgsg,Norma,jtran,itermu,itermf,iz, &
                              if0,if1,DwT,Dw,vfac,Bfield,vmi,TBout, &
                              Stokes,JKQa,JKQ,JKQC,aprof, &
                              eps20,eps21,eps22,eps23)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      type(Frequencyc2_class), intent(inout):: Fin
      type(Nindex_class), intent(in):: Norma
      type(Redc2_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: jtran,itermu,itermf,iz,if0,if1
      double precision, intent(in):: DwT,Dw,vfac,vmi, vx, vy, vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(if0:if1), intent(out):: eps20
      double precision, dimension(if0:if1), intent(out):: eps21
      double precision, dimension(if0:if1), intent(out):: eps22
      double precision, dimension(if0:if1), intent(out):: eps23
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(in):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TBout

      ! Local

      logical:: PRDc,integrate,LPRAM,lNCHLT,cohIn,conj
#if _OPENMP
      logical:: ldo
#endif

      integer:: i,itran,ktran,iterml,ith1
      integer:: ifreq,iifreq,jfreq,iti,ios,iran
      integer:: K,iQQ,K1,iPP,iq,iq1,ip,ip1
      integer:: nMl,nMu,nMf,iMl,iMl1,iMu,iMu1,iMf
      integer:: iJl,iJl1,iJlb,iJlb1,iJu,iJu1,iJu2,iJu3,iJf,iJf1
      integer:: indF,indU,indU1,indL,indL1,icom
      integer:: jjfreq,kkfreq,kwfreq0,nmfreq
#ifdef _OPENMP
      integer:: tid
#endif
      integer:: iL,iL1,kL,kL1,kLb,kLb1
      integer:: iU,iU1,kU,kU1,kU2,kU3
      integer:: mF,kF,kF1

      double precision:: omegao,omegai
      double precision:: rLl,rLu,rLf,S,rJl,rJl1,rJlb,rJlb1
      double precision:: rJu,rJu1,rJu2,rJu3,rJf,rJf1
      double precision:: rJlmax,rJumax,rJfmax
      double precision:: rMl,rMl1,rMu,rMu1,rMf
      double precision:: el,el1,eu,eu1,ef,wlf
      double precision:: rK,QQ,rK1,PP,q,q1,p,p1
      double precision:: al,au,af,auf,aul,Dw1,hau
      double precision:: at,Dfreqw,vfacw,sig
      double precision:: Norme0,Norme1
      double precision:: ftmp,f1tmp,daux,rep,imp
     !double precision:: dNorme2
      double precision:: Cl,Cl1,Clb,Clb1,Cu,Cu1,Cu2,Cu3,Cf,Cf1
      double precision:: CC,CC1,CC2,CC3
#ifdef _OPENMP
      double precision, dimension(if0:if1):: leps20
      double precision, dimension(if0:if1):: leps21
      double precision, dimension(if0:if1):: leps22
      double precision, dimension(if0:if1):: leps23
#endif

      complex(kind=8):: tmpK,hanleden,rhoc,prof,y0
      complex(kind=8):: Norme2
      complex(kind=8), dimension(if0:if1):: PRD,CRD,CRD0
      complex(kind=8), dimension(0:3,if0:if1):: tmp
      complex(kind=8), dimension(-2:2,0:2):: Jrad
      complex(kind=8), dimension(:,:,:), allocatable, target:: JKQinMV
      complex(kind=8), dimension(:,:,:), allocatable, target:: JradC
      complex(kind=8), dimension(:), allocatable, target:: Warr2


      ! Pointers
      type(Frequencyd_class), pointer:: p_frec
      type(Redd_class), pointer :: p_red
      integer, pointer:: p_mfreq
      complex(kind=8), dimension(:), pointer:: p_warr2
      complex(kind=8), dimension(:), pointer:: p_JKQ
      complex(kind=8), dimension(:), pointer:: p_JKQC


      ! Routine name
      urou = 'emiss2ord_AA'

      ! Initialize pointers
      nullify(p_frec)
      nullify(p_red)
      nullify(p_mfreq)
      nullify(p_warr2)
      nullify(p_JKQ)
      nullify(p_JKQC)

      !
      ! Initializations
      !

      !
      ! Construct mean intensity if needed
      !

      ! If dynamics
      if (dyn) then

        ! Get JKQ in comoving frame
        call getJKQstar(Fin,Geom,iz,DwT,vx,vy,vz,omega, &
                        Flgsg,Stokes,JKQa,JKQinMV)

        ! If ad-hoc asymmetries, JradC needs to be rotated
        if (size(JKQa).ge.10.and.Bfield%Bstrength(iz).gt.TINYB) &
          call fieldB_alt(JKQinMV,Fin%ggf1-Fin%ggf0+1, &
                          Flgsg,Bfield%Btheta(iz), &
                          Bfield%Bphi(iz),1)

      ! No dynamics
      else

        ! Allocate and copy
        allocate(JKQinMV(-2:2,0:2,Fin%ggf0:Fin%ggf1))
        JKQinMV = JKQC(:,:,Fin%ggf0:Fin%ggf1)

        ! Rotate only if there is some magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) &
          call fieldB_alt(JKQinMV,Fin%ggf1-Fin%ggf0+1, &
                          Flgsg,Bfield%Btheta(iz), &
                          Bfield%Bphi(iz),1)

      end if ! Dynamics


      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      at = (au + af + auf)/Dw
      hau = 2d0*(au+auf)/Dw

      ! Units normalization factor for CRD
      Norme0 = (1d5/Dw)*.5d0*sqrt(IPI)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

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

      ! Trano index
      ktran = Atom%itrano(jtran)


      !
      ! Initialize the emission coefficient
      !
!$omp parallel workshare default(none) &
!$omp shared(eps20,eps21,eps22,eps23)
      eps20 = 0d0
      eps21 = 0d0
      eps22 = 0d0
      eps23 = 0d0
!$omp end parallel workshare


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible lower terms
      do i=1,Atom%nMulti-1

        ! If there is no transition or this term is larger
        ! than the upper term of the output transition, skip
        if(i.ge.itermu.or.Atom%irad(i,itermu).eq.0) cycle

        ! Store the input lower term index
        iterml = i

        ! Get index of input transition
        itran = Atom%irad(iterml,itermu)

        ! Get index of input transition in structure
        ios = -1
        do iti=1,Atom%trano(ktran)%nt
          if (Atom%trano(ktran)%ind(iti).eq.itran) then
            ios = 1
            exit
          end if
        end do
        if (ios.lt.0) cycle

        ! If PRAM, point to the redistribution subblock
        if (PRAM) then

          p_red => Red%trani(iti)
          LPRAM = PRAM.and.p_red%RAM

        ! If not, nothing stored
        else

          LPRAM = .False.

        end if

        ! Define if non-coherent lower term
        if (NCHLT.and.allocated(Atom%NCHLT)) then
          lNCHLT = Atom%NCHLT(iz,itran)
        else
          lNCHLT = .False.
        end if

        ! Point to input transition
        p_frec => Fin%trani(iti)

        ! Predict size of interpolation block
        nmfreq = sum(p_frec%mfreq)

        ! Get input radiation field
        call getJKQin(p_frec,Fin,nmfreq,omega,JradC,JKQinMV)

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

        !
        ! Initialize the emission coefficient
        !
#ifdef _OPENMP
        leps20 = 0d0
        leps21 = 0d0
        leps22 = 0d0
        leps23 = 0d0
#endif

! dnorme2 removed from list, commented in declarations
!$omp parallel if (omp) default(none) &
!$omp private(tid,ldo,kwfreq0,iMf,rMf,mF,ef,indF,iMu,rMu) &
!$omp private(q,iq,iU,eu,indU,CRD0,Dfreqw,Norme1,ifreq,prof,iMu1) &
!$omp private(rMu1,q1,QQ,iq1,iQQ,iU1,eu1,indU1,CRD,hanleden,iMl,rMl) &
!$omp private(rMl,p,ip,iL,el,indL,iMl1,rMl1,p1,PP,ip1,iPP,iL1,el1) &
!$omp private(indL1,tmp,icom,PRDc,kkfreq,Warr2,jjfreq,iifreq,iran) &
!$omp private(ifreq,omegao,p_mfreq,jfreq,omegai,ith1,rep,imp,K1,rK1) &
!$omp private(integrate,K,rK,ftmp,f1tmp,tmpK,kLb,Clb,iJlb,rJlb,KLb1) &
!$omp private(Clb1,iJlb1,rJlb1,rhoc,Kl,Cl,iJl,rJl,kU2,Cu2,iJu2,rJu2) &
!$omp private(CC2,kL1,Cl1,iJl1,rJl1,kU3,Cu3,iJu2,rJu3,CC3,kF,Cf,iJf) &
!$omp private(rJf,kU,Cu,iJu,rJu,CC,kF1,Cf1,iJf1,rJf1,kU1,Cu1,iJu1) &
!$omp private(rJu1,CC1,conj,sig,PRD,wlf,cohIn,p_JKQC,p_warr2,y0) &
!$omp private(Norme2,p_JKQ,daux) &
!$omp shared(Atom,nMf,rJfmax,itermf,iz,nMu,rJumax,itermu,omp,lNCHLT) &
!$omp shared(ktran,iti,vpfil,Norme0,aprof,Norma,Dw,if0,if1,vfacw,at) &
!$omp shared(nMl,rJlmax,iterml,cZero,LPRAM,p_red,omega,vfac,jtran) &
!$omp shared(p_frec,Geom,Dw1,au,al,af,aul,auf,IPI42,Flgsg,TINYJS) &
!$omp shared(TINYER,TINYCO,dyn,JKQinMV,nmfreq,JradC,Fin,Jrad,TBout) &
!$omp shared(rLl) &
!$omp reduction(+: leps20,leps21,leps22,leps23)

#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
        ldo = .True.
#endif

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
            indF = Atom%i_Wind(itermf)%ind(mF,iMf)

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
                indU = Atom%i_Wind(itermu)%ind(iU,iMu)
#ifdef _OPENMP
                ! If Openmp
                if (omp) then

                  ! If NCHLT
                  if (lNCHLT) then

                    ! Check if CRD is needed here
                    ldo = .not. &
                         (iU.lt.Atom%omp_2c(ktran)%mnnU(tid,iti).or. &
                          iU.gt.Atom%omp_2c(ktran)%mxnU(tid,iti))

                  ! If CHLT
                  else

                    ! Check if CRD is needed here
                    ldo = .not. &
                          (iU.lt.Atom%omp_2c(ktran)%mnU(tid,iti).or. &
                           iU.gt.Atom%omp_2c(ktran)%mxU(tid,iti))

                  end if ! NCHLT
                end if ! OpenMP

                ! If doing CRD profile
                if (ldo) then
#endif
                ! If in file
                if (vpfil) then

                  CRD0(if0:if1) = Norme0*conjg(aprof(:, &
                               Atom%i_Vind(jtran)%ind(mF,iMf,iU,iMu)))

                ! If stored
                else if (Norma%VRAM) then

                  ! Flat spectrum contribution
                  CRD0(if0:if1) = Norme0* &
                          conjg(Norma%prof(mF,iU,iMf,iMu)%cp(if0:if1))

                ! If not stored
                else

                  ! Shift term
                  Dfreqw  = (eu - ef + Atom%Dfreq(jtran))/Dw

                  ! Normalization factors
                  Norme1  = Norma%Norm(mF,iU,iMf,iMu)*Norme0


                  !
                  ! Flat contribution. Implicit branching
                  !

                  ! For each frequency
                  do ifreq=if0,if1

                    ! Calculate profile u-f
                    call voigt(Dfreqw - omega(ifreq)*vfacw, &
                               at,prof)

                    ! Normalize
                    prof = dcmplx(dble(prof)*Norme1, &
                                  dimag(prof)*Norme0)

                    ! Flat spectrum contribution
                    CRD0(ifreq) = conjg(prof)

                  end do ! frequencies

                end if ! Storing Voigt
#ifdef _OPENMP
                end if ! Computing CRD
#endif

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
                    indU1 = Atom%i_Wind(itermu)%ind(iU1,iMu1)
#ifdef _OPENMP
                    ! If Openmp
                    if (omp) then

                      ! If NCHLT
                      if (lNCHLT) then

                        ! Check if CRD is needed here
                        ldo = .not. &
                             (iU1.lt. &
                              Atom%omp_2c(ktran)%mnnU1(tid,iti).or. &
                              iU1.gt. &
                              Atom%omp_2c(ktran)%mxnU1(tid,iti))

                      ! If CHLT
                      else

                        ! Check if CRD is needed here
                        ldo = .not. &
                             (iU1.lt. &
                              Atom%omp_2c(ktran)%mnU1(tid,iti).or. &
                              iU1.gt. &
                              Atom%omp_2c(ktran)%mxU1(tid,iti))

                      end if ! NCHLT
                    end if ! OpenMP

                    ! If doing CRD profile
                    if (ldo) then
#endif
                    ! If in file
                    if (vpfil) then

                      CRD(if0:if1) = CRD0(if0:if1) + Norme0*aprof(:, &
                              Atom%i_Vind(jtran)%ind(mF,iMf,iU1,iMu1))

                    ! If stored
                    else if (Norma%VRAM) then
                      ! Flat spectrum contribution

                      CRD(if0:if1) = CRD0(if0:if1) + Norme0* &
                               Norma%prof(mF,iU1,iMf,iMu1)%cp(if0:if1)

                    ! If not stored
                    else

                      ! Shift term
                      Dfreqw = (eu1 - ef + Atom%Dfreq(jtran))/Dw

                      ! Normalization factors
                      Norme1 = Norma%Norm(mF,iU1,iMf,iMu1)*Norme0


                      !
                      ! Flat contribution. Implicit branching
                      !

                      ! For each frequency
                      do ifreq=if0,if1

                        ! Calculate profile u'-f
                        call voigt(Dfreqw - omega(ifreq)*vfacw, &
                                   at,prof)

                        ! Normalize
                        prof = dcmplx(dble(prof)*Norme1, &
                                      dimag(prof)*Norme0)

                        ! Flat spectrum contribution
                        CRD(ifreq) = CRD0(ifreq) + prof

                      end do ! frequencies

                    end if ! Storing Voigt
#ifdef _OPENMP
                    end if ! Computing CRD
#endif

                    ! Hanle factor
                    ! TODO ATTENTION TO THIS
                    hanleden = dcmplx(hau,(eu-eu1)/Dw)

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
                      if(abs(ip).gt.1) cycle

                      ! For each mu_l
                      do iL=1,Atom%nblk(iMl,iterml)

                        ! Get eigenvalue of input lower level
                        el = Atom%eval(iL,iMl,iterml,iz)

                        ! Get indexes
                        indL = Atom%i_Wind(iterml)%ind(iL,iMl)

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
                          if(abs(ip1).gt.1) cycle

                          ! If non-coherent lower term, skip if not
                          ! diagonal
                          if (lNCHLT) then
                            if (iMl.ne.iMl1) cycle
                          end if

                          ! For each mu_l'
                          do iL1 = 1,Atom%nblk(iMl1,iterml)

                            ! If non-coherent lower term, skip if not
                            ! diagonal
                            if (lNCHLT) then
                              if (iL.ne.iL1) cycle
                            end if

                            ! Get eigenvalue of input lower' level
                            el1 = Atom%eval(iL1,iMl1,iterml,iz)

                            ! Get indexes
                            indL1 = Atom%i_Wind(iterml)%ind(iL1,iMl1)

                            ! Initialize temporal variable
                            tmp = cZero

        !
        ! Reset indexing
        !

        ! Get the component index
        icom = Atom%trano(ktran)%trani(iti)% &
                    Wind(indL1,indL,indF,indU1,indU)
#ifdef _OPENMP
        ! Check if doing this component if in omp
        if (omp) then

          ! If NCHLT
          if (lNCHLT) then

            ! If out of own ranges
            if (icom.lt.Atom%omp_2c(ktran)%nif0(tid,iti)) then

              ! Add to index and skip
              kwfreq0 = kwfreq0 + nmfreq
              cycle

            end if ! Out by left

            ! If out by right
            if (icom.gt.Atom%omp_2c(ktran)%nif1(tid,iti)) exit

          ! If CHLT
          else

            ! If out of own ranges
            if (icom.lt.Atom%omp_2c(ktran)%if0(tid,iti)) then

              ! Add to index
              kwfreq0 = kwfreq0 + nmfreq
              cycle

            end if ! Out by left

            ! If out by right
            if (icom.gt.Atom%omp_2c(ktran)%if1(tid,iti)) exit

          end if ! NCHLT
        end if ! OMP
#endif

        if (LPRAM) then
          PRDc = p_red%iPPRD(icom)
        else
          PRDc = .True.
        end if

        ! Initialize frequency
        kkfreq = 0


        !
        ! Create array of Wfunc
        !

        if (PRDc.and.nmfreq.gt.0) then

          ! Initialize array
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

          ! Initialize frequency index
          jjfreq = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fin%nran
            do ifreq=Fin%if0(iran),Fin%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! Get output frequency
              omegao = omega(ifreq)*vfac - Atom%Dfreq(jtran)

              ! Point to dimension
              p_mfreq => p_frec%mfreq(iifreq)

              ! Skip coherent
              if (p_mfreq.lt.1) cycle

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq + 1
                kkfreq = kkfreq + 1

                ! Get input frequency
                omegai = p_frec%omega(jjfreq) - Atom%Dfreq(itran)

                ! Initialize Warr2
                Warr2(kkfreq) = cZero

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

          ! If storing
          if (LPRAM) then

            p_red%iPPRD(icom) = .False.
            do jfreq=kkfreq-nmfreq+1,kkfreq
              rep = dble(Warr2(jfreq))
              imp = dimag(Warr2(jfreq))
              if (rep.le.1e-30) rep = 0d0
              if (abs(imp).le.1e-30) imp = 0d0
              p_red%PWarr2(kwfreq0+jfreq) = &
                                          cmplx(real(rep),real(imp))
            end do

          end if ! If storing
        end if ! Initialized

        ! For each K'
        do K1=abs(iPP),2

          ! Get real value
          rK1 = dble(K1)

          ! Initialize flag for K1
          integrate = .True.

          ! For each K
          do K=abs(iQQ),2

            ! Get real value
            rK = dble(K)

            ! Racah algebra
            ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

            ! If forbidden (3j-sym=0), skip
            if(abs(ftmp).lt.TINYJS) cycle

            ! Racah algebra
            f1tmp = fun3j(1d0,1d0,rK1,-p,p1,-PP,Flgsg)

            ! If forbidden (3j-sym=0), skip
            if(abs(f1tmp).lt.TINYJS) cycle

            ! Add the rest of the factor
            ftmp = ftmp*Flgsg%sg(iq)*sqrt(2d0*rK+1d0)
            f1tmp = f1tmp*Flgsg%sg(ip1)*sqrt(2d0*rK1+1d0)

            ! Reset temporal variable
            tmpK = cZero

            ! For each Jlb
            do kLb=1,Atom%nblk(iMl,iterml)

              ! Get eigenvector for lower b level
              Clb = Atom%evec(kLb,iL,iMl,iterml,iz)

              ! If coefficient too small, skip
              if(abs(Clb).lt.TINYEV) cycle

              ! Get J level index
              iJlb = Atom%iJval(kLb,iMl,iterml)

              ! Get angular momentum
              rJlb = Atom%rJval(iJlb,iterml)

              ! For each Jlb'
              do kLb1=1,Atom%nblk(iMl1,iterml)

                ! Get eigenvector for lower b level
                Clb1 = Atom%evec(kLb1,iL1,iMl1,iterml,iz)

                ! If coefficient too small, skip
                if(abs(Clb1).lt.TINYEV) cycle

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
                if(abs(rhoc).lt.TINYER) cycle

                ! Add coefficients to rhoKQ
                rhoc = Flgsg%sg(nint(rJlb-rMl))*Clb*Clb1*rhoc

                ! For each Jl
                do kL=1,Atom%nblk(iMl,iterml)

                  ! Get eigenvector for lower level
                  Cl = Atom%evec(kL,iL,iMl,iterml,iz)

                  ! If coefficient too small, skip
                  if(abs(Cl).lt.TINYEV) cycle

                  ! Get J level index
                  iJl = Atom%iJval(kL,iMl,iterml)

                  ! Get angular momentum
                  rJl = Atom%rJval(iJl,iterml)

                  ! For each Ju''
                  do kU2=1,Atom%nblk(iMu,itermu)

                    ! Get eigenvector for upper'' level
                    Cu2 = Atom%evec(kU2,iU,iMu,itermu,iz)

                    ! If coefficient too small, skip
                    if(abs(Cu2).lt.TINYEV) cycle

                    ! Get J level index
                    iJu2 = Atom%iJval(kU2,iMu,itermu)

                    ! Get angular momentum
                    rJu2 = Atom%rJval(iJu2,itermu)

                    ! Add coefficient to dipole strength
                    CC2 = Cl*Cu2*Atom%rdip(itran)% &
                                      rdip(ip,iMu,iMl,iJu2,iJl)

                    ! If coefficient too small, skip
                    if(abs(CC2).lt.TINYCO) cycle

                    ! For each Jl'
                    do kL1=1,Atom%nblk(iMl1,iterml)

                      ! Get eigenvector for lower' level
                      Cl1 = Atom%evec(kL1,iL1,iMl1,iterml,iz)

                      ! If coefficient too small, skip
                      if(abs(Cl1).lt.TINYEV) cycle

                      ! Get J level index
                      iJl1 = Atom%iJval(kL1,iMl1,iterml)

                      ! Get angular momentum
                      rJl1 = Atom%rJval(iJl1,iterml)

                      ! For each Ju'''
                      do kU3=1,Atom%nblk(iMu1,itermu)

                        ! Get eigenvector for upper''' level
                        Cu3 = Atom%evec(kU3,iU1,iMu1,itermu,iz)

                        ! If coefficient too small, skip
                        if(abs(Cu3).lt.TINYEV) cycle

                        ! Get J level index
                        iJu3 = Atom%iJval(kU3,iMu1,itermu)

                        ! Get angular momentum
                        rJu3 = Atom%rJval(iJu3,itermu)

                        ! KCUT
!                       if (abs(rJu3-rJu2).gt.rK1) cycle

                        ! Add coefficients to dipole strength
                        CC3 = Cl1*Cu3* &
                              Atom%rdip(itran)% &
                                   rdip(ip1,iMu1,iMl1,iJu3,iJl1)

                        ! If coefficient too small, skip
                        if(abs(CC3).lt.TINYCO) cycle

        !
        ! Reset identation
        !

        ! For each Jf
        do kF=1,Atom%nblk(iMf,itermf)

          ! Get eigenvector for final level
          Cf = Atom%evec(kF,mF,iMf,itermf,iz)

          ! If coefficient too small, skip
          if(abs(Cf).lt.TINYEV) cycle

          ! Get J level index
          iJf = Atom%iJval(kF,iMf,itermf)

          ! Get angular momentum
          rJf = Atom%rJval(iJf,itermf)

          ! For each Ju
          do kU=1,Atom%nblk(iMu,itermu)

            ! Get eigenvector for upper level
            Cu = Atom%evec(kU,iU,iMu,itermu,iz)

            ! If coefficient too small, skip
            if(abs(Cu).lt.TINYEV) cycle

            ! Get J level index
            iJu = Atom%iJval(kU,iMu,itermu)

            ! Get angular momentum
            rJu = Atom%rJval(iJu,itermu)

            ! Add coefficients to dipole strength
            CC = Cf*Cu*Atom%rdip(jtran)% &
                            rdip(iq,iMu,iMf,iJu,iJf)

            ! If coefficient too small, skip
            if(abs(CC).lt.TINYCO) cycle

            ! For each Jf'
            do kF1=1,Atom%nblk(iMf,itermf)

              ! Get eigenvector for final' level
              Cf1 = Atom%evec(kF1,mF,iMf,itermf,iz)

              ! If coefficient too small, skip
              if(abs(Cf1).lt.TINYEV) cycle

              ! Get J level index
              iJf1 = Atom%iJval(kF1,iMf,itermf)

              ! Get angular momentum
              rJf1 = Atom%rJval(iJf1,itermf)

              ! For each Ju'
              do kU1=1,Atom%nblk(iMu1,itermu)

                ! Get eigenvector for upper' level
                Cu1 = Atom%evec(kU1,iU1,iMu1,itermu,iz)

                ! If coefficient too small, skip
                if(abs(Cu1).lt.TINYEV) cycle

                ! Get J level index
                iJu1 = Atom%iJval(kU1,iMu1,itermu)

                ! Get angular momentum
                rJu1 = Atom%rJval(iJu1,itermu)

                ! Add coefficients to dipole strength
                CC1 = Cf1*Cu1*Atom%rdip(jtran)% &
                                   rdip(iq1,iMu1,iMf,iJu1,iJf1)

                ! If coefficient big enough, add contribution to
                ! temporal variable
                if(abs(CC1).gt.TINYCO) &
                  tmpK = f1tmp*ftmp*CC*CC1*CC2*CC3*rhoc + tmpK

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


            !
            ! Integral over input frequencies
            !

            if (integrate) then

              !
              ! Check if conjugated
              conj = iPP.lt.0
              sig = Flgsg%sg(iPP)

              ! Initialize 2nd order part
              PRD = cZero

              ! Difference between l and f energies
              wlf = el - ef

              ! Check if static and coherent
              cohIn = dyn.or.abs(wlf).gt.0d0

              ! If coherent
              if (minval(p_frec%mfreq).lt.1) then

                ! Just point
                if (conj) then
                  p_JKQC(Fin%ggf0:Fin%ggf1) => JKQinMV(-iPP,K,:)
                else
                  p_JKQC(Fin%ggf0:Fin%ggf1) => JKQinMV(iPP,K,:)
                end if
              end if

              ! If storing Warr
              if (LPRAM) then
                if (nmfreq.gt.0) then
                  allocate(p_warr2(nmfreq))
                  p_warr2 = dcmplx(p_red% &
                            Pwarr2(kwfreq0+1:kwfreq0+nmfreq))
                end if
              ! If not storing
              else
                if (allocated(Warr2)) p_warr2 => Warr2
              end if

              ! Initialize frequency indexes
              jjfreq = 0
              kkfreq = 0

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

                    ! Full coherent
                    if (cohIn) then

                      ! Input frequency
                      omegai = omega(ifreq)*vfac - wlf

                      ! Get JKQ
                      y0 = getJKQinnu(omega(Fin%ggf0:Fin%ggf1), &
                                      p_JKQC, &
                                      ifreq-Fin%ggf0+1, &
                                      Fin%ggf1-Fin%ggf0+1,omegai)

                      ! Fully coherent contribution
                      if (conj) then
                        PRD(ifreq) = sig*CRD(ifreq)*conjg(y0)
                      else
                        PRD(ifreq) = CRD(ifreq)*y0
                      end if

                    ! Interpolate
                    else

                      ! Fully coherent contribution
                      if (conj) then
                        PRD(ifreq) = sig*CRD(ifreq)* &
                                     conjg(p_JKQC(ifreq))
                      else
                        PRD(ifreq) = CRD(ifreq)*p_JKQC(ifreq)
                      end if

                    end if

                    ! Skip rest
                    cycle

                  end if ! Coherent wing

                  ! If conjugate
                  if (conj) then

                    ! Point to positive Q
                    p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,-iPP,K1)

                    ! Integrate
                    PRD(ifreq) = sig*sum(conjg(p_JKQ)* &
                             p_frec%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                             p_warr2(kkfreq+1:kkfreq+p_mfreq))

                  ! Not conjugate
                  else

                    ! Point to positive Q
                    p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,iPP,K1)

                    ! Integrate
                    PRD(ifreq) = sum(p_JKQ* &
                             p_frec%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                             p_warr2(kkfreq+1:kkfreq+p_mfreq))

                  end if

                  ! Initialize
                  Norme2 = sum(p_warr2(kkfreq+1:kkfreq+p_mfreq)* &
                               p_frec%W_freq(jjfreq+1:jjfreq+p_mfreq))

                  ! Normalize to the first order profile
                  PRD(ifreq) = PRD(ifreq)*CRD(ifreq)/Norme2

                  ! This old normalization (only real part) was
                  ! problematic for non-axial cases with
                  ! angle-dependent
                 !! Normalize real part to the first order profile
                 !dNorme2 = dble(Norme2)
                 !if (dNorme2.gt.0d0) then
                 !  rep = dble(PRD(ifreq))*dble(CRD(ifreq))/ &
                 !                                     dNorme2
                 !  PRD(ifreq) = dcmplx(rep,dimag(PRD(ifreq)))
                 !else
                 !  PRD(ifreq) = dcmplx(0d0,dimag(PRD(ifreq)))
                 !end if

                  ! Update indexes
                  jjfreq = jjfreq + p_mfreq
                  kkfreq = kkfreq + p_mfreq

                end do ! Output frequencies
              end do ! Output frequencies ranges

              ! Clean p_warr2
              if (LPRAM.and.nmfreq.gt.0) deallocate(p_warr2)
              nullify(p_warr2)

              ! Clean JKQ
              nullify(p_JKQ)
              nullify(p_JKQC)

              !
              ! Flat spectrum contribution. Implicit branching ratio
              !

              ! For each output frequency
              do iran=1,Fin%nran
                do ifreq=Fin%if0(iran),Fin%if1(iran)

                  ! Substract the flat spectrum part due to just
                  ! radiative excitation
                  PRD(ifreq) = PRD(ifreq) - CRD(ifreq)*Jrad(iPP,K1)

                end do ! Output frequencies
              end do ! Output frequencies ranges

              integrate = .False.

            end if ! If flagged

            !
            ! Compute the main part of emiss2ord
            !

            ! For each output frequency
            do iran=1,Fin%nran
              do ifreq=Fin%if0(iran),Fin%if1(iran)

                ! Add TKQ to the PRD contribution and accumulate
                tmp(0,ifreq) = tmpK*TBout(0,-iQQ,K)*PRD(ifreq) + &
                               tmp(0,ifreq)
                tmp(1,ifreq) = tmpK*TBout(1,-iQQ,K)*PRD(ifreq) + &
                               tmp(1,ifreq)
                tmp(2,ifreq) = tmpK*TBout(2,-iQQ,K)*PRD(ifreq) + &
                               tmp(2,ifreq)
                tmp(3,ifreq) = tmpK*TBout(3,-iQQ,K)*PRD(ifreq) + &
                               tmp(3,ifreq)

              end do ! Output frequencies
            end do ! Output frequencies ranges

          end do ! K1
        end do ! K

        ! Update redistribution initial index
        kwfreq0 = kwfreq0 + nmfreq

        ! Apply hanle factor and Einstein coefficient
        daux = (2d0*rLl+1d0)*Atom%Ecoeff(iterml,itermu)
#ifdef _OPENMP
        leps20 = leps20 + dble(tmp(0,:)/hanleden)*daux
        leps21 = leps21 + dble(tmp(1,:)/hanleden)*daux
        leps22 = leps22 + dble(tmp(2,:)/hanleden)*daux
        leps23 = leps23 + dble(tmp(3,:)/hanleden)*daux
#else
        eps20 = eps20 + dble(tmp(0,:)/hanleden)*daux
        eps21 = eps21 + dble(tmp(1,:)/hanleden)*daux
        eps22 = eps22 + dble(tmp(2,:)/hanleden)*daux
        eps23 = eps23 + dble(tmp(3,:)/hanleden)*daux
#endif
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
!$omp end parallel

#ifdef _OPENMP
!$omp parallel workshare default(none) &
!$omp shared(eps20,eps21,eps22,eps23,leps20,leps21,leps22,leps23)
        eps20 = eps20 + leps20
        eps21 = eps21 + leps21
        eps22 = eps22 + leps22
        eps23 = eps23 + leps23
!$omp end parallel workshare
#endif
      end do ! Terms

      ! Apply common factor
      daux = 3d0*.5d0*IPI42*(2d0*rLu+1d0)* &
             Atom%Ecoeff(itermu,itermf)*1d-10/(c*Dw)
!$omp parallel workshare default(none) &
!$omp shared(eps20,eps21,eps22,eps23) shared(daux)
      eps20 = eps20*daux
      eps21 = eps21*daux
      eps22 = eps22*daux
      eps23 = eps23*daux
!$omp end parallel workshare

      ! Clean pointers
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_warr2)) nullify(p_warr2)

      return

      end subroutine emiss2ord_AA

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
      !!          omega(dfloat): Frequency array\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!    Fin(Frequencyc2_class): Structure with the input frequency
      !!                            information\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!             jdir(integer): Output direction in scattering
      !!                            indexing\n
      !!            jtran(integer): Output transition index\n
      !!           itermu(integer): Upper term of output transition\n
      !!           itermf(integer): Lower term of output transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!               DwT(dfloat): Thermal part of Doppler width\n
      !!                Dw(dfloat): Doppler width output transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!               vmi(dfloat): Microturbulent velocity\n
      !!      TBout(dcmplx(:,:,:)): Geometrical tensor in the output
      !!                            direction in the magnetic field
      !!                            reference frame\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!          eps20(dfloat(:)): Intensity emissivity\n
      !!          eps21(dfloat(:)): Q emissivity\n
      !!          eps22(dfloat(:)): U emissivity\n
      !!          eps23(dfloat(:)): V emissivity
      subroutine emiss2ord_AD(Atom,Geom,emerging,vx,vy,vz,omega,Red, &
                              Fin,Flgsg,Norma,jdir,jtran,itermu, &
                              itermf,iz,if0,if1,DwT,Dw,vfac,vmi, &
                              TBout,Stokes,JKQ,aprof, &
                              eps20,eps21,eps22,eps23)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Frequencyc2_class), intent(inout):: Fin
      type(Nindex_class), intent(in):: Norma
      type(Redc2_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: emerging
      integer, intent(in):: jtran,itermu,itermf,iz,if0,if1,jdir
      double precision, intent(in):: DwT,Dw,vfac,vmi, vx, vy, vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(if0:if1), intent(out):: eps20
      double precision, dimension(if0:if1), intent(out):: eps21
      double precision, dimension(if0:if1), intent(out):: eps22
      double precision, dimension(if0:if1), intent(out):: eps23
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TBout

      ! Local

      logical:: PRDc,integrate,LPRAM,lNCHLT,cohIn
#if _OPENMP
      logical:: ldo
#endif
      logical, dimension(Geom%nScatt):: skip_scatt

      integer:: i,itran,ktran,iterml,ith1,iph1
      integer:: ifreq,iifreq,jfreq,iti,ios,iran
      integer:: K,iQQ,K1,iPP,iq,iq1,ip,ip1
      integer:: nMl,nMu,nMf,iMl,iMl1,iMu,iMu1,iMf
      integer:: iJl,iJl1,iJlb,iJlb1,iJu,iJu1,iJu2,iJu3,iJf,iJf1
      integer:: indF,indU,indU1,indL,indL1,icom
      integer:: jjfreq,jjfreq0,kkfreq,kkfreq0,kkfreq0b,kwfreq0
      integer:: llfreq0,nmfreq,ish1,nfs,nskip,nScatt,stype
#ifdef _OPENMP
      integer:: tid
#endif
      integer:: iL,iL1,kL,kL1,kLb,kLb1
      integer:: iU,iU1,kU,kU1,kU2,kU3
      integer:: mF,kF,kF1
      integer, dimension(:), allocatable:: i_scatt

      double precision:: omegao,omegai,vfac1
      double precision:: rLl,rLu,rLf,S,rJl,rJl1,rJlb,rJlb1
      double precision:: rJu,rJu1,rJu2,rJu3,rJf,rJf1
      double precision:: rJlmax,rJumax,rJfmax
      double precision:: rMl,rMl1,rMu,rMu1,rMf
      double precision:: el,el1,eu,eu1,ef,wlf
      double precision:: rK,QQ,rK1,PP,q,q1,p,p1
      double precision:: al,au,af,auf,aul,Dw1,hau
      double precision:: at,Dfreqw,vfacw
      double precision:: Norme0,Norme1
      double precision:: ftmp,f1tmp,daux,rep,imp
     !double precision:: dNorme2
      double precision:: Cl,Cl1,Clb,Clb1,Cu,Cu1,Cu2,Cu3,Cf,Cf1
      double precision:: CC,CC1,CC2,CC3,cost,sint,cosc,sinc
      double precision, dimension(0:3):: StokesM
      double precision, dimension(:,:), allocatable:: Stokesin
#ifdef _OPENMP
      double precision, dimension(if0:if1):: leps20
      double precision, dimension(if0:if1):: leps21
      double precision, dimension(if0:if1):: leps22
      double precision, dimension(if0:if1):: leps23
#endif

      complex(kind=8):: intgr,tmpK,hanleden,rhoc,prof
      complex(kind=8):: Norme2,PRDin
      complex(kind=8), dimension(if0:if1):: PRD,CRD,CRD0
      complex(kind=8), dimension(0:3,if0:if1):: tmp
      complex(kind=8), dimension(-2:2,0:2):: Jrad
      complex(kind=8), dimension(:), allocatable, target:: Warr2
      complex(kind=8), dimension(:), allocatable:: Warr2xW
      complex(kind=8), dimension(:), allocatable:: intergrin


      ! Pointers
      type(Frequencyd_class), pointer:: p_frec
      type(Redd_class), pointer :: p_red
      integer, pointer:: p_mfreq
      complex(kind=8), dimension(:), pointer:: p_warr2


      ! Routine name
      urou = 'emiss2ord_AD'

      ! Initialize pointers
      nullify(p_frec)
      nullify(p_red)
      nullify(p_mfreq)
      nullify(p_warr2)


      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      at = (au + af + auf)/Dw
      hau = 2d0*(au+auf)/Dw

      ! Units normalization factor for CRD
      Norme0 = (1d5/Dw)*.5d0*sqrt(IPI)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

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

      ! Trano index
      ktran = Atom%itrano(jtran)


      !
      ! Initialize the emission coefficient
      !
!$omp parallel workshare default(none) &
!$omp shared(eps20,eps21,eps22,eps23)
      eps20 = 0d0
      eps21 = 0d0
      eps22 = 0d0
      eps23 = 0d0
!$omp end parallel workshare

      ! If there are frequencies
#ifdef _OPENMP
#else
      if (Fin%mxfreq.gt.0) then
        allocate(Warr2xW(Fin%mxfreq))
        allocate(intergrin(Fin%mxfreq))
      end if
#endif


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible lower terms
      do i=1,Atom%nMulti-1

        ! If there is no transition or this term is larger
        ! than the upper term of the output transition, skip
        if(i.ge.itermu.or.Atom%irad(i,itermu).eq.0) cycle

        ! Store the input lower term index
        iterml = i

        ! Get index of input transition
        itran = Atom%irad(iterml,itermu)

        ! Get index of input transition in structure
        ios = -1
        do iti=1,Atom%trano(ktran)%nt
          if (Atom%trano(ktran)%ind(iti).eq.itran) then
            ios = 1
            exit
          end if
        end do
        if (ios.lt.0) cycle

        ! Check if forward
        if (jtran.eq.itran.and. &
            Geom%V_CScatt(1).ge.1d0) then
          nfs = 1
        else
          nfs = 0
        end if

        ! If PRAM, point to the redistribution subblock
        if (PRAM) then

          p_red => Red%trani(iti)
          LPRAM = PRAM.and.p_red%RAM

        ! If not, nothing stored
        else

          LPRAM = .False.

        end if

        ! Define if non-coherent lower term
        if (NCHLT.and.allocated(Atom%NCHLT)) then
          lNCHLT = Atom%NCHLT(iz,itran)
        else
          lNCHLT = .False.
        end if

        ! If not storing or dynamic in quadrature
        if ((.not.LPRAM.or.dyn).and..not.emerging) then

          ! Copy restricted
          skip_scatt = Geom%skip_ksc
          nskip = Geom%nskip
          nScatt = Geom%nScatt - nskip
          i_scatt = Geom%k_scatt

        ! Storing, static, or LOS
        else

          ! Copy total
          skip_scatt = Geom%skip_jsc
          nskip = 0
          nScatt = Geom%nScatt
          i_scatt = Geom%j_scatt

        end if

        ! Point to input transition
        p_frec => Fin%trani(iti)

        ! Get interpolated intensity
        call getStkin(Geom,p_frec,Fin,omega,vx,vy,vz,jdir, &
                      nfs,Stokesin,Stokes)

        ! Get the 'flat' JKQ for this input transition
        JRad = JKQ(:,:,itran)

        ! Redistribution size
        nmfreq = sum(p_frec%mfreq)*(nScatt-nfs)

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

        !
        ! Initialize the emission coefficient
        !
#ifdef _OPENMP
        leps20 = 0d0
        leps21 = 0d0
        leps22 = 0d0
        leps23 = 0d0
#endif

! dnorme2 removed from list, commented in declarations
!$omp parallel if (omp) default(none) &
!$omp private(tid,ldo,kwfreq0,iMf,rMf,mF,ef,indF,iMu,rMu,q,iq,iU,eu) &
!$omp private(indU,CRD0,Dfreqw,Norme1,ifreq,prof,iMu1,rMu1,q1,QQ) &
!$omp private(iq1,iQQ,iU1,eu1,indU1,CRD,hanleden,iMl,rMl,p,ip,iL,el) &
!$omp private(indL,iMl1,rMl1,p1,PP,ip1,iPP,iL1,el1,indL1,tmp,icom) &
!$omp private(PRDc,Warr2,jjfreq0,kkfreq0,iifreq,iran,ifreq,p_mfreq) &
!$omp private(omegao,ish1,stype,jfreq,jjfreq,kkfreq,omegai,rep,imp) &
!$omp private(K1,rK1,integrate,K,rK,ftmp,f1tmp,tmpK,kLb,Clb,iJlb) &
!$omp private(rJlb,kLb1,Clb1,iJlb1,rJlb1,rhoc,kL,Cl,iJl,rJl,kU2,Cu2) &
!$omp private(iJu2,rJu2,CC2,kL1,Cl1,iJl1,rJl1,kU3,Cu3,iJu3,rJu3,CC3) &
!$omp private(kF,Cf,iJf,rJf,kU,Cu,iJu,rJu,CC,kF1,Cf1,iJf1,rJf1,kU1) &
!$omp private(Cu1,iJu1,rJu1,CC1,PRD,wlf,cohIn,p_warr2,llfreq0,ith1) &
!$omp private(iph1,cost,sint,cosc,sinc,vfac1,StokesM,intgr,kkfreq0b) &
!$omp private(Warr2xW,intergrin,Norme2,PRDin,daux) &
!$omp shared(Atom,nMf,rJfmax,itermf,iz,nMu,rJumax,itermu,omp,lNCHLT) &
!$omp shared(ktran,iti,vpfil,if0,if1,Norme0,aprof,jtran,Norma,Dw) &
!$omp shared(omega,vfacw,at,rJlmax,iterml,nMl,cZero,nmfreq,LPRAM) &
!$omp shared(p_red,p_frec,vfac,Geom,skip_scatt,itran,Dw1,au,al,af) &
!$omp shared(auf,aul,IPI42,Flgsg,TINYJS,TINYEV,TINYER,TINYCO,dyn) &
!$omp shared(axial,Fin,Stokesin,nScatt,nfs,Jrad,TBout,rLl,jdir) &
!$omp reduction(+: leps20,leps21,leps22,leps23)

      ! If there are frequencies
#ifdef _OPENMP
      if (Fin%mxfreq.gt.0) then
        allocate(Warr2xW(Fin%mxfreq))
        allocate(intergrin(Fin%mxfreq))
      end if
#endif

#ifdef _OPENMP
        tid = omp_get_thread_num() + 1
        ldo = .True.
#endif

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
            indF = Atom%i_Wind(itermf)%ind(mF,iMf)

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
                indU = Atom%i_Wind(itermu)%ind(iU,iMu)
#ifdef _OPENMP
                ! If Openmp
                if (omp) then

                  ! If NCHLT
                  if (lNCHLT) then

                    ! Check if CRD is needed here
                    ldo = .not. &
                         (iU.lt.Atom%omp_2c(ktran)%mnnU(tid,iti).or. &
                          iU.gt.Atom%omp_2c(ktran)%mxnU(tid,iti))

                  ! If CHLT
                  else

                    ! Check if CRD is needed here
                    ldo = .not. &
                         (iU.lt.Atom%omp_2c(ktran)%mnU(tid,iti).or. &
                          iU.gt.Atom%omp_2c(ktran)%mxU(tid,iti))

                  end if ! NCHLT
                end if ! OpenMP

                ! If doing CRD profile
                if (ldo) then
#endif
                ! If in file
                if (vpfil) then

                  CRD0(if0:if1) = Norme0*conjg(aprof(:, &
                               Atom%i_Vind(jtran)%ind(mF,iMf,iU,iMu)))

                ! If stored
                else if (Norma%VRAM) then

                  ! Flat spectrum contribution
                  CRD0(if0:if1) = Norme0* &
                          conjg(Norma%prof(mF,iU,iMf,iMu)%cp(if0:if1))

                ! If not stored
                else

                  ! Shift term
                  Dfreqw  = (eu - ef + Atom%Dfreq(jtran))/Dw

                  ! Normalization factors
                  Norme1  = Norma%Norm(mF,iU,iMf,iMu)*Norme0


                  !
                  ! Flat contribution. Implicit branching
                  !

                  ! For each frequency
                  do ifreq=if0,if1

                    ! Calculate profile u-f
                    call voigt(Dfreqw - omega(ifreq)*vfacw, &
                               at,prof)

                    ! Normalize
                    prof = dcmplx(dble(prof)*Norme1, &
                                  dimag(prof)*Norme0)

                    ! Flat spectrum contribution
                    CRD0(ifreq) = conjg(prof)

                  end do ! frequencies

                end if ! Storing Voigt
#ifdef _OPENMP
                end if ! Computing CRD
#endif

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
                    indU1 = Atom%i_Wind(itermu)%ind(iU1,iMu1)
#ifdef _OPENMP
                    ! If Openmp
                    if (omp) then

                      ! If NCHLT
                      if (lNCHLT) then

                        ! Check if CRD is needed here
                        ldo = .not. &
                             (iU1.lt. &
                              Atom%omp_2c(ktran)%mnnU1(tid,iti).or. &
                              iU1.gt. &
                              Atom%omp_2c(ktran)%mxnU1(tid,iti))

                      ! If CHLT
                      else

                        ! Check if CRD is needed here
                        ldo = .not. &
                             (iU1.lt. &
                              Atom%omp_2c(ktran)%mnU1(tid,iti).or. &
                              iU1.gt. &
                              Atom%omp_2c(ktran)%mxU1(tid,iti))

                      end if ! NCHLT
                    end if ! OpenMP

                    ! If doing CRD profile
                    if (ldo) then
#endif
                    ! If in file
                    if (vpfil) then

                      CRD(if0:if1) = CRD0(if0:if1) + Norme0*aprof(:, &
                              Atom%i_Vind(jtran)%ind(mF,iMf,iU1,iMu1))

                    ! If stored
                    else if (Norma%VRAM) then
                      ! Flat spectrum contribution

                      CRD(if0:if1) = CRD0(if0:if1) + Norme0* &
                               Norma%prof(mF,iU1,iMf,iMu1)%cp(if0:if1)

                    ! If not stored
                    else

                      ! Shift term
                      Dfreqw = (eu1 - ef + Atom%Dfreq(jtran))/Dw

                      ! Normalization factors
                      Norme1 = Norma%Norm(mF,iU1,iMf,iMu1)*Norme0


                      !
                      ! Flat contribution. Implicit branching
                      !

                      ! For each frequency
                      do ifreq=if0,if1

                        ! Calculate profile u'-f
                        call voigt(Dfreqw - omega(ifreq)*vfacw, &
                                   at,prof)

                        ! Normalize
                        prof = dcmplx(dble(prof)*Norme1, &
                                      dimag(prof)*Norme0)

                        ! Flat spectrum contribution
                        CRD(ifreq) = CRD0(ifreq) + prof

                      end do ! frequencies

                    end if ! Storing Voigt
#ifdef _OPENMP
                    end if ! Computing CRD
#endif

                    ! Hanle factor
                    ! TODO ATTENTION TO THIS
                    hanleden = dcmplx(hau,(eu-eu1)/Dw)

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
                      if(abs(ip).gt.1) cycle

                      ! For each mu_l
                      do iL=1,Atom%nblk(iMl,iterml)

                        ! Get eigenvalue of input lower level
                        el = Atom%eval(iL,iMl,iterml,iz)

                        ! Get indexes
                        indL = Atom%i_Wind(iterml)%ind(iL,iMl)

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
                          if(abs(ip1).gt.1) cycle

                          ! If non-coherent lower term, skip if not
                          ! diagonal
                          if (lNCHLT) then
                            if (iMl.ne.iMl1) cycle
                          end if

                          ! For each mu_l'
                          do iL1 = 1,Atom%nblk(iMl1,iterml)

                            ! If non-coherent lower term, skip if not
                            ! diagonal
                            if (lNCHLT) then
                              if (iL.ne.iL1) cycle
                            end if

                            ! Get eigenvalue of input lower' level
                            el1 = Atom%eval(iL1,iMl1,iterml,iz)

                            ! Get indexes
                            indL1 = Atom%i_Wind(iterml)%ind(iL1,iMl1)

                            ! Initialize temporal variable
                            tmp = cZero

        !
        ! Reset indexing
        !

        ! Get the component index
        icom = Atom%trano(ktran)%trani(iti)% &
                    Wind(indL1,indL,indF,indU1,indU)
#ifdef _OPENMP
        ! Check if doing this component if in omp
        if (omp) then

          ! If NCHLT
          if (lNCHLT) then

            ! If out of own ranges
            if (icom.lt.Atom%omp_2c(ktran)%nif0(tid,iti)) then

              ! Add to index, and skip
              kwfreq0 = kwfreq0 + nmfreq
              cycle

            end if ! Out by left

            ! If out by right
            if (icom.gt.Atom%omp_2c(ktran)%nif1(tid,iti)) exit

          ! If CHLT
          else

            ! If out of own ranges
            if (icom.lt.Atom%omp_2c(ktran)%if0(tid,iti)) then

              ! Add to index, and skip
              kwfreq0 = kwfreq0 + nmfreq
              cycle

            end if ! Out by left

            ! If out by right
            if (icom.gt.Atom%omp_2c(ktran)%if1(tid,iti)) exit

          end if ! NCHLT
        end if ! OMP
#endif

        if (LPRAM) then
          PRDc = p_red%iPPRD(icom)
        else
          PRDc = .True.
        end if


        !
        ! Create array of Wfunc
        !

        if (PRDc.and.nmfreq.gt.0) then

          ! Initialize array
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

          ! Initialize frequency index
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
              omegao = omega(ifreq)*vfac - Atom%Dfreq(jtran)

              ! For each scattering angle
              do ish1=1,Geom%nScatt

                ! Non-present scattering angle
                if (skip_scatt(ish1)) cycle

                ! Check forward scattering two-terms
                if (itran.eq.jtran.and. &
                    Geom%V_CScatt(ish1).ge.1d0) cycle

                ! Check backward
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
                  omegai = p_frec%omega(jjfreq) - &
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

          ! If storing
          if (LPRAM) then

            p_red%iPPRD(icom) = .False.
            do jfreq = kkfreq-nmfreq+1,kkfreq
              rep = dble(Warr2(jfreq))
              imp = dimag(Warr2(jfreq))
              if (rep.le.1e-30) rep = 0d0
              if (abs(imp).le.1e-30) imp = 0d0
              p_red%Pwarr2(kwfreq0+jfreq) = &
                                          cmplx(real(rep),real(imp))
            end do

          end if ! If storing
        end if ! Initialized

        ! For each K'
        do K1=abs(iPP),2

          ! Get real value
          rK1 = dble(K1)

          ! Initialize flag for K1
          integrate = .True.

          ! For each K
          do K=abs(iQQ),2

            ! Get real value
            rK = dble(K)

            ! Racah algebra
            ftmp = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

            ! If forbidden (3j-sym=0), skip
            if(abs(ftmp).lt.TINYJS) cycle

            ! Racah algebra
            f1tmp = fun3j(1d0,1d0,rK1,-p,p1,-PP,Flgsg)

            ! If forbidden (3j-sym=0), skip
            if(abs(f1tmp).lt.TINYJS) cycle

            ! Add the rest of the factor
            ftmp = ftmp*Flgsg%sg(iq)*sqrt(2d0*rK+1d0)
            f1tmp = f1tmp*Flgsg%sg(ip1)*sqrt(2d0*rK1+1d0)

            ! Reset temporal variable
            tmpK = cZero

            ! For each Jlb
            do kLb=1,Atom%nblk(iMl,iterml)

              ! Get eigenvector for lower b level
              Clb = Atom%evec(kLb,iL,iMl,iterml,iz)

              ! If coefficient too small, skip
              if(abs(Clb).lt.TINYEV) cycle

              ! Get J level index
              iJlb = Atom%iJval(kLb,iMl,iterml)

              ! Get angular momentum
              rJlb = Atom%rJval(iJlb,iterml)

              ! For each Jlb'
              do kLb1=1,Atom%nblk(iMl1,iterml)

                ! Get eigenvector for lower b level
                Clb1 = Atom%evec(kLb1,iL1,iMl1,iterml,iz)

                ! If coefficient too small, skip
                if(abs(Clb1).lt.TINYEV) cycle

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
                if(abs(rhoc).lt.TINYER) cycle

                ! Add coefficients to rhoKQ
                rhoc = Flgsg%sg(nint(rJlb-rMl))*Clb*Clb1*rhoc

                ! For each Jl
                do kL=1,Atom%nblk(iMl,iterml)

                  ! Get eigenvector for lower level
                  Cl = Atom%evec(kL,iL,iMl,iterml,iz)

                  ! If coefficient too small, skip
                  if(abs(Cl).lt.TINYEV) cycle

                  ! Get J level index
                  iJl = Atom%iJval(kL,iMl,iterml)

                  ! Get angular momentum
                  rJl = Atom%rJval(iJl,iterml)

                  ! For each Ju''
                  do kU2=1,Atom%nblk(iMu,itermu)

                    ! Get eigenvector for upper'' level
                    Cu2 = Atom%evec(kU2,iU,iMu,itermu,iz)

                    ! If coefficient too small, skip
                    if(abs(Cu2).lt.TINYEV) cycle

                    ! Get J level index
                    iJu2 = Atom%iJval(kU2,iMu,itermu)

                    ! Get angular momentum
                    rJu2 = Atom%rJval(iJu2,itermu)

                    ! Add coefficient to dipole strength
                    CC2 = Cl*Cu2*Atom%rdip(itran)% &
                                      rdip(ip,iMu,iMl,iJu2,iJl)

                    ! If coefficient too small, skip
                    if(abs(CC2).lt.TINYCO) cycle

                    ! For each Jl'
                    do kL1=1,Atom%nblk(iMl1,iterml)

                      ! Get eigenvector for lower' level
                      Cl1 = Atom%evec(kL1,iL1,iMl1,iterml,iz)

                      ! If coefficient too small, skip
                      if(abs(Cl1).lt.TINYEV) cycle

                      ! Get J level index
                      iJl1 = Atom%iJval(kL1,iMl1,iterml)

                      ! Get angular momentum
                      rJl1 = Atom%rJval(iJl1,iterml)

                      ! For each Ju'''
                      do kU3=1,Atom%nblk(iMu1,itermu)

                        ! Get eigenvector for upper''' level
                        Cu3 = Atom%evec(kU3,iU1,iMu1,itermu,iz)

                        ! If coefficient too small, skip
                        if(abs(Cu3).lt.TINYEV) cycle

                        ! Get J level index
                        iJu3 = Atom%iJval(kU3,iMu1,itermu)

                        ! Get angular momentum
                        rJu3 = Atom%rJval(iJu3,itermu)

                        ! KCUT
!                       if (abs(rJu3-rJu2).gt.rK1) cycle

                        ! Add coefficients to dipole strength
                        CC3 = Cl1*Cu3* &
                              Atom%rdip(itran)% &
                                   rdip(ip1,iMu1,iMl1,iJu3,iJl1)

                        ! If coefficient too small, skip
                        if(abs(CC3).lt.TINYCO) cycle

        !
        ! Reset identation
        !

        ! For each Jf
        do kF=1,Atom%nblk(iMf,itermf)

          ! Get eigenvector for final level
          Cf = Atom%evec(kF,mF,iMf,itermf,iz)

          ! If coefficient too small, skip
          if(abs(Cf).lt.TINYEV) cycle

          ! Get J level index
          iJf = Atom%iJval(kF,iMf,itermf)

          ! Get angular momentum
          rJf = Atom%rJval(iJf,itermf)

          ! For each Ju
          do kU=1,Atom%nblk(iMu,itermu)

            ! Get eigenvector for upper level
            Cu = Atom%evec(kU,iU,iMu,itermu,iz)

            ! If coefficient too small, skip
            if(abs(Cu).lt.TINYEV) cycle

            ! Get J level index
            iJu = Atom%iJval(kU,iMu,itermu)

            ! Get angular momentum
            rJu = Atom%rJval(iJu,itermu)

            ! Add coefficients to dipole strength
            CC = Cf*Cu*Atom%rdip(jtran)% &
                            rdip(iq,iMu,iMf,iJu,iJf)

            ! If coefficient too small, skip
            if(abs(CC).lt.TINYCO) cycle

            ! For each Jf'
            do kF1=1,Atom%nblk(iMf,itermf)

              ! Get eigenvector for final' level
              Cf1 = Atom%evec(kF1,mF,iMf,itermf,iz)

              ! If coefficient too small, skip
              if(abs(Cf1).lt.TINYEV) cycle

              ! Get J level index
              iJf1 = Atom%iJval(kF1,iMf,itermf)

              ! Get angular momentum
              rJf1 = Atom%rJval(iJf1,itermf)

              ! For each Ju'
              do kU1=1,Atom%nblk(iMu1,itermu)

                ! Get eigenvector for upper' level
                Cu1 = Atom%evec(kU1,iU1,iMu1,itermu,iz)

                ! If coefficient too small, skip
                if(abs(Cu1).lt.TINYEV) cycle

                ! Get J level index
                iJu1 = Atom%iJval(kU1,iMu1,itermu)

                ! Get angular momentum
                rJu1 = Atom%rJval(iJu1,itermu)

                ! Add coefficients to dipole strength
                CC1 = Cf1*Cu1*Atom%rdip(jtran)% &
                                   rdip(iq1,iMu1,iMf,iJu1,iJf1)

                ! If coefficient big enough, add contribution to
                ! temporal variable
                if(abs(CC1).gt.TINYCO) &
                  tmpK = f1tmp*ftmp*CC*CC1*CC2*CC3*rhoc + tmpK

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


            !
            ! Integral over input frequencies
            !

            if (integrate) then

              ! Initialize 2nd order part
              PRD = cZero

              ! Difference between l and f energies
              wlf = el - ef

              ! Check if static and coherent
              cohIn = dyn.or.abs(wlf).gt.0d0

              ! If storing Warr2
              if (LPRAM) then
                if (nmfreq.gt.0) then
                  allocate(p_warr2(nmfreq))
                  p_warr2 = dcmplx(p_red% &
                                   Pwarr2(kwfreq0+1:kwfreq0+nmfreq))
                end if
              ! If not storing
              else
                if (allocated(Warr2)) p_warr2 => Warr2
              end if

              ! If axial symmetry
              if (axial) then

                ! Initialize frequency indexes
                jjfreq0 = 0
                llfreq0 = 0
                kkfreq0 = 0

                ! For each output frequency
                iifreq = 0
                do iran=1,Fin%nran
                  do ifreq=Fin%if0(iran),Fin%if1(iran)

                    ! Advance index
                    iifreq = iifreq + 1

                    ! Point to dimension
                    p_mfreq => p_frec%mfreq(iifreq)

                    ! For each polar direction
                    do ith1=1,Geom%nTh

                      ! For each azimuthal direction
                      do iph1=1,Geom%nPh2

                        ! Scattering index
                        ish1 = Geom%i_scatt(iph1,ith1,jdir)

                        ! Special treatment if forward for two
                        ! terms
                        if ((jtran.eq.itran.and. &
                             Geom%V_CScatt(ish1).ge.1d0).or. &
                            (p_mfreq.lt.1)) then

                          ! Interpolate
                          if (cohIn) then

                            ! Input frequency
                            omegai = omega(ifreq)*vfac - wlf

                            ! If there are dynamics
                            if (dyn) then

                              ! Get director cosine
                              cost = Geom%V_mu(ith1)

                              ! Calculate Doppler shift factor
                              vfac1 = 1d0 - vz*cost

                              ! We will be using the inverse
                              vfac1 = 1d0/vfac1

                              ! Shift
                              omegai = omegai*vfac1

                            end if

                            ! Interpol
                            StokesM = getStkinnu(omega, &
                                                 Stokes(:,:,1,ith1), &
                                                 ifreq,omegai)

                          ! Fully coherent
                          else

                            StokesM = Stokes(:,ifreq,1,ith1)

                          end if

                          ! Sum over Stokes parameters,
                          ! integrand
                          intgr = sum(StokesM* &
                                      Geom%TB(:,iPP,K1, &
                                              iph1,ith1,iz))

                          PRD(ifreq) = PRD(ifreq) + &
                                       intgr*Geom%W_mu(ith1)* &
                                       Geom%W_mux2(iph1)* &
                                       CRD(ifreq)

                        ! Non-forward 2-term scattering
                        else

                          !
                          ! Find initial index for kkfreq

                          ! Scattering index
                          ish1 = i_scatt(ish1)

                          ! Shift in indexes
                          kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                          ! Multiply Warr2 and weights
                          Warr2xW(1:p_mfreq) = &
                              p_warr2(kkfreq0b+1:kkfreq0b+p_mfreq)* &
                              p_frec%W_freq(llfreq0+1:llfreq0+p_mfreq)

                          ! Sum Stokes
                          intergrin(1:p_mfreq) = &
                              Stokesin(0,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(0,iPP,K1,iph1,ith1,iz) + &
                              Stokesin(1,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(1,iPP,K1,iph1,ith1,iz) + &
                              Stokesin(2,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(2,iPP,K1,iph1,ith1,iz) + &
                              Stokesin(3,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(3,iPP,K1,iph1,ith1,iz)

                          ! Compute norm
                          Norme2 = sum(Warr2xW(1:p_mfreq))

                          ! Integrate
                          PRDin = sum(Warr2xW(1:p_mfreq)* &
                                      intergrin(1:p_mfreq))

                          ! Normalize to the first order profile
                          ! and add the directional weights
                          PRD(ifreq) = PRD(ifreq) + &
                                       PRDin*CRD(ifreq)* &
                                       Geom%W_mu(ith1)* &
                                       Geom%W_mux2(iph1)/Norme2

                          ! This old normalization (only real part)
                          ! was problematic for non-axial cases with
                          ! angle-dependent
                         !! Normalize real part to the first order
                         !! profile and add the directional
                         !! weights
                         !dNorme2 = dble(Norme2)
                         !if (dNorme2.gt.0d0) then
                         !  rep = dble(PRDin)* &
                         !        dble(CRD(ifreq))/dNorme2
                         !  imp = dimag(PRDin)
                         !  PRD(ifreq) = PRD(ifreq) + &
                         !               dcmplx(rep,imp)* &
                         !               Geom%W_mu(ith1)* &
                         !               Geom%W_mux2(iph1)
                         !else
                         !  PRD(ifreq) = PRD(ifreq) + &
                         !               dcmplx(0d0, &
                         !                      dimag(PRDin)* &
                         !                       Geom%W_mu(ith1)* &
                         !                       Geom%W_mux2(iph1))
                         !end if ! Positive norm

                        end if ! Type of scattering

                      end do ! azimuthal nodes

                      ! Update jjfreq
                      jjfreq0 = jjfreq0 + p_mfreq

                    end do ! polar nodes

                    ! Update jjfreq
                    llfreq0 = llfreq0 + p_mfreq
                    kkfreq0 = kkfreq0 + p_mfreq*(nScatt - nfs)

                  end do ! output frequencies
                end do ! output frequencies ranges

              ! If not axial symmetric
              else

                ! Initialize indexes
                jjfreq0 = 0
                llfreq0 = 0
                kkfreq0 = 0

                ! For each output frequency
                iifreq = 0
                do iran=1,Fin%nran
                  do ifreq=Fin%if0(iran),Fin%if1(iran)

                    ! Advance index
                    iifreq = iifreq + 1

                    ! Point to dimension
                    p_mfreq => p_frec%mfreq(iifreq)

                    ! For each polar direction
                    do ith1=1,Geom%nTh

                      ! For each azimuthal direction
                      do iph1=1,Geom%nPh2

                        ! Scattering index
                        ish1 = Geom%i_scatt(iph1,ith1,jdir)

                        ! Special treatment if forward for two
                        ! terms
                        if ((jtran.eq.itran.and. &
                             Geom%V_CScatt(ish1).ge.1d0).or. &
                            (p_mfreq.lt.1)) then

                          ! Interpolate
                          if (cohIn) then

                            ! Input frequency
                            omegai = omega(ifreq)*vfac - wlf

                            ! If there are dynamics
                            if (dyn) then

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
                                            ifreq,omegai)

                          ! Full coherent
                          else

                            StokesM = Stokes(:,ifreq,iph1,ith1)

                          end if

                          ! Sum over Stokes parameters,
                          ! integrand
                          intgr = sum(StokesM* &
                                      Geom%TB(:,iPP,K1, &
                                              iph1,ith1,iz))

                          PRD(ifreq) = PRD(ifreq) + CRD(ifreq)* &
                                       intgr*Geom%W_mu(ith1)* &
                                       Geom%W_mux2(iph1)

                        ! Non-forward 2-term scattering
                        else

                          ! Scattering index
                          ish1 = i_scatt(ish1)

                          ! Shift in indexes
                          kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                          ! Multiply Warr2 and weights
                          Warr2xW(1:p_mfreq) = &
                              p_warr2(kkfreq0b+1:kkfreq0b+p_mfreq)* &
                              p_frec%W_freq(llfreq0+1:llfreq0+p_mfreq)

                          ! Sum Stokes
                          intergrin(1:p_mfreq) = &
                              Stokesin(0,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(0,iPP,K1,iph1,ith1,iz) + &
                              Stokesin(1,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(1,iPP,K1,iph1,ith1,iz) + &
                              Stokesin(2,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(2,iPP,K1,iph1,ith1,iz) + &
                              Stokesin(3,jjfreq0+1:jjfreq0+p_mfreq)* &
                              Geom%TB(3,iPP,K1,iph1,ith1,iz)

                          ! Compute norm
                          Norme2 = sum(Warr2xW(1:p_mfreq))

                          ! Integrate
                          PRDin = sum(Warr2xW(1:p_mfreq)* &
                                      intergrin(1:p_mfreq))

                          ! Normalize to the first order profile
                          ! and add the directional weights
                          PRD(ifreq) = PRD(ifreq) + &
                                       PRDin*CRD(ifreq)* &
                                       Geom%W_mu(ith1)* &
                                       Geom%W_mux2(iph1)/Norme2

                          ! This old normalization (only real part)
                          ! was problematic for non-axial cases with
                          ! angle-dependent
                         !! Normalize real part to the first order
                         !! profile and add the directional weights
                         !dNorme2 = dble(Norme2)
                         !if (dNorme2.gt.0d0) then
                         !  rep = dble(PRDin)* &
                         !        dble(CRD(ifreq))/dNorme2
                         !  imp = dimag(PRDin)
                         !  PRD(ifreq) = PRD(ifreq) + &
                         !               dcmplx(rep,imp)* &
                         !                Geom%W_mu(ith1)* &
                         !                Geom%W_mux2(iph1)
                         !else
                         !  PRD(ifreq) = PRD(ifreq) + &
                         !               dcmplx(0d0, &
                         !                      dimag(PRDin)* &
                         !                       Geom%W_mu(ith1)* &
                         !                       Geom%W_mux2(iph1))
                         !end if ! Valid norm

                          ! Update indexes
                          jjfreq0 = jjfreq0 + p_mfreq

                        end if ! Type of scattering

                      end do ! azimuthal nodes
                    end do ! polar nodes

                    ! Advance
                    llfreq0 = llfreq0 + p_mfreq
                    kkfreq0 = kkfreq0 + p_mfreq*(nScatt - nfs)

                  end do ! output frequencies
                end do ! output frequencies ranges

              end if ! Axial symmetry

              ! Clean p_warr2
              if (LPRAM.and.nmfreq.gt.0) deallocate(p_warr2)
              nullify(p_warr2)

              !
              ! Flat spectrum contribution. Implicit branching ratio
              !

              ! For each output frequency
              do iran=1,Fin%nran
                do ifreq=Fin%if0(iran),Fin%if1(iran)

                  ! Substract the flat spectrum part due to just
                  ! radiative excitation
                  PRD(ifreq) = PRD(ifreq) - CRD(ifreq)*Jrad(iPP,K1)

                end do ! Output frequencies
              end do ! Output frequencies ranges

              integrate = .False.

            end if ! If flagged

            !
            ! Compute the main part of emiss2ord
            !

            ! For each output frequency
            do iran=1,Fin%nran
              do ifreq=Fin%if0(iran),Fin%if1(iran)

                ! Add TKQ to the PRD contribution and accumulate
                tmp(0,ifreq) = tmpK*TBout(0,-iQQ,K)*PRD(ifreq) + &
                               tmp(0,ifreq)
                tmp(1,ifreq) = tmpK*TBout(1,-iQQ,K)*PRD(ifreq) + &
                               tmp(1,ifreq)
                tmp(2,ifreq) = tmpK*TBout(2,-iQQ,K)*PRD(ifreq) + &
                               tmp(2,ifreq)
                tmp(3,ifreq) = tmpK*TBout(3,-iQQ,K)*PRD(ifreq) + &
                               tmp(3,ifreq)

              end do ! Output frequencies
            end do ! Output frequencies ranges

          end do ! K1
        end do ! K

        ! Update redistribution initial index
        kwfreq0 = kwfreq0 + nmfreq

        ! Apply hanle factor and Einstein coefficient
        daux = (2d0*rLl+1d0)*Atom%Ecoeff(iterml,itermu)
#ifdef _OPENMP
        leps20 = leps20 + dble(tmp(0,:)/hanleden)*daux
        leps21 = leps21 + dble(tmp(1,:)/hanleden)*daux
        leps22 = leps22 + dble(tmp(2,:)/hanleden)*daux
        leps23 = leps23 + dble(tmp(3,:)/hanleden)*daux
#else
        eps20 = eps20 + dble(tmp(0,:)/hanleden)*daux
        eps21 = eps21 + dble(tmp(1,:)/hanleden)*daux
        eps22 = eps22 + dble(tmp(2,:)/hanleden)*daux
        eps23 = eps23 + dble(tmp(3,:)/hanleden)*daux
#endif
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

#ifdef _OPENMP
      if (allocated(Warr2xW)) deallocate(Warr2xW,intergrin)
#endif

!$omp end parallel

#ifdef _OPENMP
!$omp parallel workshare default(none) &
!$omp shared(eps20,eps21,eps22,eps23,leps20,leps21,leps22,leps23)
        eps20 = eps20 + leps20
        eps21 = eps21 + leps21
        eps22 = eps22 + leps22
        eps23 = eps23 + leps23
!$omp end parallel workshare
#endif

      end do ! Terms

      ! Apply common factor
      daux = 3d0*.5d0*IPI42*(2d0*rLu+1d0)* &
             Atom%Ecoeff(itermu,itermf)*1d-10/(c*Dw)
!$omp parallel workshare default(none) &
!$omp shared(eps20,eps21,eps22,eps23) shared(daux)
      eps20 = eps20*daux
      eps21 = eps21*daux
      eps22 = eps22*daux
      eps23 = eps23*daux
!$omp end parallel workshare

      ! Clean pointers
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_warr2)) nullify(p_warr2)

      return

      end subroutine emiss2ord_AD

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the second order emission coefficient in absence
      !! of magnetic fields under the angle-averaged approximation.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!                vx(dfloat): Velocity vector along X\n
      !!                vy(dfloat): Velocity vector along Y\n
      !!                vz(dfloat): Velocity vector along Z\n
      !!          omega(dfloat(:)): Frequency array\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!    Fin(Frequencyc2_class): Structure with the input frequency
      !!                            information\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!            jtran(integer): Output transition index\n
      !!           itermu(integer): Upper term of output transition\n
      !!           itermf(integer): Lower term of output transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!               DwT(dfloat): Thermal part of Doppler width\n
      !!                Dw(dfloat): Doppler width output transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!               vmi(dfloat): Microturbulent velocity\n
      !!      TBout(dcmplx(:,:,:)): Geometrical tensor in the output
      !!                            direction\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!     JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!          eps20(dfloat(:)): Intensity emissivity\n
      !!          eps21(dfloat(:)): Q emissivity\n
      !!          eps22(dfloat(:)): U emissivity\n
      !!          eps23(dfloat(:)): V emissivity
      subroutine emiss2ordNB_AA(Atom,Geom,vx,vy,vz,omega,Red,Fin, &
                                Flgsg,Norma,jtran,itermu,itermf,iz, &
                                if0,if1,DwT,Dw,vfac,vmi,TBout, &
                                Stokes,JKQa,JKQ,JKQC,aprof, &
                                eps20,eps21,eps22,eps23)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Frequencyc2_class), intent(inout):: Fin
      type(Nindex_class), intent(in):: Norma
      type(Redc2_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: jtran,itermu,itermf,iz,if0,if1
      double precision, intent(in):: DwT,Dw,vfac,vmi,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(if0:if1), intent(out):: eps20
      double precision, dimension(if0:if1), intent(out):: eps21
      double precision, dimension(if0:if1), intent(out):: eps22
      double precision, dimension(if0:if1), intent(out):: eps23
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), &
                                             target, intent(in):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TBout

      ! Local

      logical:: PRDc,LPRAM,cohIn,conj

      integer:: ith1,i,iterml,itran,iti,icom,ios
      integer:: mF,iU,iU1,iL,iL1,iifreq,iran,ifreq,jfreq,iR
      integer:: K,iQ,K1,iQ1,iQl,Kmin,Kmax,Kl,ktran
      integer:: jjfreq,jjfreq0,nmfreq
      integer:: kwfreq0,indF,indU,indU1,indL,indL1
#ifdef _OPENMP
      integer:: tid
#endif

      double precision:: rLu,rLl,rLf,S,rJu,rJu1,rJl,rJl1,rJf
      double precision:: rK,Q,rK1,Q1,Ql,rKl,wlf
      double precision:: eu,eu1,el,el1,ef,au,af,al,auf,aul,at
      double precision:: f61,f62,f63,f64,f0tmp,ftmp,f1tmp,fKfj
      double precision:: vfacw,Dw1,Dfreqw,Norme0,Norme1,rep,imp
      double precision:: omegao,omegai,daux,sig
     !double precision:: dNorme2

      complex(kind=8):: hanleden,prof,rhoc
      complex(kind=8):: Norme2,y0
      complex(kind=8), dimension(if0:if1):: PRD,CRD0,CRD
      complex(kind=8), dimension(0:3,if0:if1):: tmp
      complex(kind=8), dimension(-2:2,0:2):: Jrad
      complex(kind=8), dimension(:), allocatable, target:: Warr2
      complex(kind=8), dimension(:,:,:), allocatable, target:: JKQinMV
      complex(kind=8), dimension(:,:,:), allocatable, target:: JradC

      ! Pointers
      type(Frequencyd_class), pointer:: p_frec
      type(Redd_class), pointer:: p_red
      integer, pointer:: p_mfreq
      complex(kind=8), dimension(:), pointer:: p_warr2
      complex(kind=8), dimension(:), pointer:: p_JKQ
      complex(kind=8), dimension(:), pointer:: p_JKQC


      ! Routine name
      urou = 'emiss2ordNB_AA'

      ! Initialize pointers
      nullify(p_frec)
      nullify(p_red)
      nullify(p_mfreq)
      nullify(p_warr2)
      nullify(p_JKQ)
      nullify(p_JKQC)


      !
      ! Construct mean intensity if needed
      !

      ! If dynamics
      if (dyn) then

        ! Get JKQ in comoving frame
        call getJKQstar(Fin,Geom,iz,DwT,vx,vy,vz,omega, &
                        Flgsg,Stokes,JKQa,JKQinMV)

      end if ! Dynamics


      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      at = (au + af + auf)/Dw

      ! Units normalization factor for CRD profile
      Norme0 = (1d5/Dw)*.5d0*sqrt(IPI)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLf = Atom%rLval(itermf)

      ! Trano index
      ktran = Atom%itrano(jtran)

      !
      ! Initialize the emission coefficient
      !
      eps20 = 0d0
      eps21 = 0d0
      eps22 = 0d0
      eps23 = 0d0
      PRDc = .True.


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible lower terms
      do i=1,Atom%nMulti-1

        ! If there is no transition or this term is larger
        ! than the upper term of the output transition, skip
        if(i.ge.itermu.or.Atom%irad(i,itermu).eq.0) cycle

        ! Store the input lower term index
        iterml = i

        ! Get index of input transition
        itran = Atom%irad(iterml,itermu)

        ! Get index of input transition in structure
        ios = -1
        do iti=1,Atom%trano(ktran)%nt
          if (Atom%trano(ktran)%ind(iti).eq.itran) then
            ios = 1
            exit
          end if
        end do
        if (ios.lt.0) cycle

        ! If PRAM, point to the redistribution subblock
        if (PRAM) then

          p_red => Red%trani(iti)
          LPRAM = PRAM.and.p_red%RAM

        ! If not, nothing stored
        else

          LPRAM = .False.

        end if

        ! Point to input transition
        p_frec => Fin%trani(iti)

        ! Predict size of interpolation block
        nmfreq = sum(p_frec%mfreq)

        ! If dynamic
        if (dyn) then

          ! Get input radiation field
          call getJKQin(p_frec,Fin,nmfreq,omega,JradC,JKQinMV)

        else

          ! Get input radiation field
          call getJKQin(p_frec,Fin,nmfreq,omega,JradC, &
                        JKQC(:,:,Fin%ggf0:Fin%ggf1))

        end if

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

          ! Get indexes
          indF = Atom%irho(itermf)%irho_ij(mF)

          ! Get Jf
          rJf = Atom%rJval(mF,itermf)

          ! For each Ju
          do iU=1,Atom%nJ(itermu)

            ! Get eigenvalue upper level
            eu = Atom%FSfreq(iU,itermu) - Atom%TRfreq(itermu)

            ! Get indexes
            indU = Atom%irho(itermu)%irho_ij(iU)

            ! Get Ju
            rJu = Atom%rJval(iU,itermu)

            f61 = fun6j(rJu,rJf,1d0,rLf,rLu,S,Flgsg)

            if (abs(f61).lt.TINYJS) cycle

            f61 = f61*(2d0*rJu+1d0)*(2d0*rJf+1d0)


            !
            ! Flat contribution. Implicit branching
            !

            ! If in file
            if (vpfil) then

!$omp parallel workshare default(none) &
!$omp shared(CRD0,if0,if1,Norme0,aprof,Atom,jtran,mF,iU)
              CRD0(if0:if1) = Norme0*conjg(aprof(:, &
                                     Atom%i_Vind(jtran)%indNB(mF,iU)))
!$omp end parallel workshare

            ! If stored
            else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(CRD0,if0,if1,Norme0,Norma,mF,iU)
              CRD0(if0:if1) = Norme0* &
                              conjg(Norma%prof(mF,iU,1,1)%cp(if0:if1))
!$omp end parallel workshare

            ! If not stored
            else

              ! Shift term
              Dfreqw  = (eu  - ef + Atom%Dfreq(jtran))/Dw

              ! Normalization factors
              Norme1  = Norma%Norm(mF,iU,1,1)*Norme0

              ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,Dfreqw,omega,vfacw,at,Norme1,Norme0,CRD0)
              do ifreq=if0,if1

                ! Calculate profile u-f
                call voigt(Dfreqw - omega(ifreq)*vfacw,at,prof)

                ! Normalize
                prof = dcmplx(dble(prof)*Norme1,dimag(prof)*Norme0)

                ! Flat spectrum contribution
                CRD0(ifreq) = conjg(prof)

              end do ! frequencies
!$omp end parallel do

            end if ! Storing Voigt

            ! For each Ju'
            do iU1=1,Atom%nJ(itermu)

              ! Get eigenvalue upper' level
              eu1 = Atom%FSfreq(iU1,itermu) - Atom%TRfreq(itermu)

              ! Get indexes
              indU1 = Atom%irho(itermu)%irho_ij(iU1)

              ! Get Ju'
              rJu1 = Atom%rJval(iU1,itermu)

              f62 = fun6j(rJu1,rJf,1d0,rLf,rLu,S,Flgsg)

              if (abs(f62).lt.TINYJS) cycle

              f62 = f62*(2d0*rJu1+1d0)*Flgsg%sg(nint(rJu1+rJf+1d0))

              !
              ! Flat contribution. Implicit branching
              !

              ! If in file
              if (vpfil) then

!$omp parallel workshare default(none) &
!$omp shared(CRD,CRD0,if0,if1,Norme0,aprof,Atom,jtran,mF,iU1)
                CRD(if0:if1) = CRD0(if0:if1) + Norme0*aprof(:, &
                                     Atom%i_Vind(jtran)%indNB(mF,iU1))
!$omp end parallel workshare

              ! If stored
              else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(CRD,CRD0,if0,if1,Norme0,Norma,mF,iU1)
                CRD(if0:if1) = CRD0(if0:if1) + Norme0*&
                               Norma%prof(mF,iU1,1,1)%cp(if0:if1)
!$omp end parallel workshare

              ! If not stored
              else

                ! Shift term
                Dfreqw = (eu1 - ef + Atom%Dfreq(jtran))/Dw

                ! Normalization factor
                Norme1 = Norma%Norm(mF,iU1,1,1)*Norme0

                ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,Dfreqw,omega,vfacw,at,Norme1,Norme0,CRD,CRD0)
                do ifreq=if0,if1

                  ! Calculate profile u'-f
                  call voigt(Dfreqw - omega(ifreq)*vfacw,at,prof)

                  ! Normalize
                  prof = dcmplx(dble(prof)*Norme1,dimag(prof)*Norme0)

                  ! Flat spectrum contribution
                  CRD(ifreq) = CRD0(ifreq) + prof

                end do ! frequencies
!$omp end parallel do

              end if ! Storing Voigt

              !
              ! Continue with the 2nd order emissivity
              !

              ! For each Jl
              do iL=1,Atom%nJ(iterml)

                ! Get eigenvalue of input lower level
                el = Atom%FSfreq(iL,iterml) - Atom%TRfreq(iterml)

                ! Get indexes
                indL = Atom%irho(iterml)%irho_ij(iL)

                ! Get Jl
                rJl = Atom%rJval(iL,iterml)

                f63 = fun6j(rJu,rJl,1d0,rLl,rLu,S,Flgsg)

                if (abs(f63).lt.TINYJS) cycle

                f63 = f63*sqrt(2d0*rJl+1d0)

                ! For each Jl'
                do iL1=1,Atom%nJ(iterml)

                  ! Get eigenvalue of input lower' level
                  el1 = Atom%FSfreq(iL1,iterml) - Atom%TRfreq(iterml)

                  ! Get indexes
                  indL1 = Atom%irho(iterml)%irho_ij(iL1)

                  ! Get Jl1
                  rJl1 = Atom%rJval(iL1,iterml)

                  f64 = fun6j(rJu1,rJl1,1d0,rLl,rLu,S,Flgsg)

                  if (abs(f64).lt.TINYJS) cycle

                  f64 = f64*sqrt(2d0*rJl1+1d0)

                  ! Hanle factor
                  ! TODO ATTENTION TO THIS
                  hanleden = dcmplx(2d0*(au+auf),eu-eu1)/Dw

                  ! Initialize temporal variable
                  tmp = cZero

        !
        ! Reset indexing
        !

        ! Get the component index
        icom = Atom%trano(ktran)%trani(iti)% &
                    WindNB(indL1,indL,indF,indU1,indU)

        if (LPRAM) then
          PRDc = p_red%iPPRD(icom)
        else
          PRDc = .True.
        end if


        !
        ! Create array of Wfunc
        !

        if (PRDc.and.nmfreq.gt.0) then

          ! Initialize array
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
!$omp private(tid,Warr2,jjfreq0,iifreq,iran,ifreq,p_mfreq,omegao) &
!$omp private(jfreq,jjfreq,omegai,ith1,red,imp) &
!$omp shared(Atom,nmfreq,p_frec,omp,Fin,vfac,Geom,el,el1,eu,eu1,ef) &
!$omp shared(al,au,af,aul,auf,IPI42,LPRAM,p_red) &

#ifdef _OPENMP
          tid = omp_get_thread_num() + 1
#endif

          ! Initialize Warr2
          if (nmfreq.gt.0) then
!$omp workshare
            Warr2 = cZero
!$omp end workshare
          end if

          ! Initialize frequency index
          jjfreq0 = 0

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
                  jjfreq0 = jjfreq0 + p_mfreq
                  cycle
                end if
                ! If out of range above
                if (iifreq.gt.Fin%oif1(tid)) exit
              end if
#endif
              ! Get output frequency
              omegao = omega(ifreq)*vfac - Atom%Dfreq(jtran)

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq

                ! Get input frequency
                omegai = p_frec%omega(jjfreq) - Atom%Dfreq(itran)

                ! For each direction in the integral AA quadrature
                do ith1=1,Geom%nThAA

                  ! Add the contribution to the angular integral
                  ! of the redistribution function
                  Warr2(jjfreq) = Warr2(jjfreq) + &
                                 Geom%W_muAA(ith1)* &
                                 Wfunc(omegai,omegao, &
                                       Dw,Dw1,el,el1,eu,eu1,ef, &
                                              al,au,af,aul,auf, &
                                       Geom%V_muAA(ith1), &
                                       Geom%V_siAA(ith1),0)*IPI42

                end do ! Direction quadrature
              end do ! input frequencies

              ! Update
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! output frequencies
          end do ! output frequency ranges
!$omp barrier
!$omp flush (Warr2)

          ! If storing
          if (LPRAM) then
!$omp single
            p_red%iPPRD(icom) = .False.
!$omp end single
!$omp do
            do jfreq=1,nmfreq
              rep = dble(Warr2(jfreq))
              imp = dimag(Warr2(jfreq))
              if (rep.le.1e-30) rep = 0d0
              if (abs(imp).le.1e-30) imp = 0d0
              p_red%PWarr2(kwfreq0+jfreq) = &
                                          cmplx(real(rep),real(imp))
            end do
!$omp end do
          end if ! If storing
!$omp end parallel
        end if ! Initialized

        f0tmp = f61*f62*f63*f64

        ! For each K
        do K=0,Krad

          ! Get real value
          rK = dble(K)

          ! Racah algebra
          ftmp = f0tmp*sqrt(2d0*rK+1d0)

          do iQ=-K,K

            Q = dble(iQ)

            !
            ! Check if conjugated
            conj = iQ.lt.0
            sig = Flgsg%sg(iQ)


      !
      ! Reset identation
      !

      !
      ! Integral over input frequencies
      !

      ! Difference between l and f energies
      wlf = el - ef

      ! Check if static and coherent
      cohIn = dyn.or.abs(wlf).gt.0d0

      ! If coherent
      if (minval(p_frec%mfreq).lt.1) then

        ! If dynamic
        if (dyn) then

          ! Just point
          if (conj) then
            p_JKQC(Fin%ggf0:Fin%ggf1) => JKQinMV(-iQ,K,:)
          else
            p_JKQC(Fin%ggf0:Fin%ggf1) => JKQinMV(iQ,K,:)
          end if

        ! If static
        else

          ! Just point
          if (conj) then
            p_JKQC(Fin%ggf0:Fin%ggf1) => JKQC(-iQ,K,Fin%ggf0:Fin%ggf1)
          else
            p_JKQC(Fin%ggf0:Fin%ggf1) => JKQC(iQ,K,Fin%ggf0:Fin%ggf1)
          end if

        end if
      end if

! dnorme2 removed from list, commented in declarations
!$omp parallel default(none) &
!$omp private(p_warr2,tid,jjfreq,iifreq,iran,ifreq,p_mfreq,omegai) &
!$omp private(y0,p_JKQ,Norme2) &
!$omp shared(LPRAM,nmfreq,p_red,kwfreq0,Warr2,Fin,p_frec,cohIn) &
!$omp shared(omega,vfac,wlf,p_JKQC,PRD,CRD,conj,JradC,iQ,K,sig) &

      ! If storing Warr
      if (LPRAM) then
        if (nmfreq.gt.0) then
!$omp single
          allocate(p_warr2(nmfreq))
!$omp end single
!$omp workshare
          p_warr2 = dcmplx(p_red%Pwarr2(kwfreq0+1:kwfreq0+nmfreq))
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

      ! Initialize frequency indexes
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

            ! Interpolate
            if (cohIn) then

              ! Input frequency
              omegai = omega(ifreq)*vfac - wlf

              ! Get JKQ
              y0 = getJKQinnu(omega(Fin%ggf0:Fin%ggf1), &
                              p_JKQC, &
                              ifreq-Fin%ggf0+1, &
                              Fin%ggf1-Fin%ggf0+1,omegai)

              ! Fully coherent contribution
              if (conj) then
                PRD(ifreq) = sig*CRD(ifreq)*conjg(y0)
              else
                PRD(ifreq) = CRD(ifreq)*y0
              end if

            ! Full coherent
            else

              ! Fully coherent contribution
              if (conj) then
                PRD(ifreq) = sig*CRD(ifreq)*conjg(p_JKQC(ifreq))
              else
                PRD(ifreq) = CRD(ifreq)*p_JKQC(ifreq)
              end if

            end if

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

          ! If conjugate
          if (conj) then

            ! Point to positive Q
            p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,-iQ,K)

            ! Integrate
            PRD(ifreq) = sig*sum(conjg(p_JKQ)* &
                             p_frec%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                             p_warr2(jjfreq+1:jjfreq+p_mfreq))

          ! Not conjugate
          else

            ! Point
            p_JKQ => JradC(jjfreq+1:jjfreq+p_mfreq,iQ,K)

            ! Integrate
            PRD(ifreq) = sum(p_JKQ* &
                             p_frec%W_freq(jjfreq+1:jjfreq+p_mfreq)* &
                             p_warr2(jjfreq+1:jjfreq+p_mfreq))

          end if

          ! Compute norm
          Norme2 = sum(p_warr2(jjfreq+1:jjfreq+p_mfreq)* &
                       p_frec%W_freq(jjfreq+1:jjfreq+p_mfreq))

          ! Normalize to the first order profile
          PRD(ifreq) = PRD(ifreq)*CRD(ifreq)/Norme2

          ! This old normalization (only real part) was problematic
          ! for non-axial cases with angle-dependent
         !! Normalize real part to the first order profile
         !dNorme2 = dble(Norme2)
         !if (dNorme2.gt.0d0) then
         !  rep = dble(PRD(ifreq))*dble(CRD(ifreq))/dNorme2
         !  PRD(ifreq) = dcmplx(rep,dimag(PRD(ifreq)))
         !else
         !  PRD(ifreq) = dcmplx(0d0,dimag(PRD(ifreq)))
         !end if

          ! Update index
          jjfreq = jjfreq + p_mfreq

        end do ! Output frequencies
      end do ! Output frequencies ranges
!$omp end parallel

      ! Clean p_warr2
      if (LPRAM.and.nmfreq.gt.0) deallocate(p_warr2)
      nullify(p_warr2)

      ! Clean JKQ
      nullify(p_JKQ)
      nullify(p_JKQC)

      !
      ! Flat spectrum contribution. Implicit branching ratio
      !

      ! For each output frequency
      do iran=1,Fin%nran
        do ifreq=Fin%if0(iran),Fin%if1(iran)

          ! Substract the flat spectrum part due to just
          ! radiative excitation
          PRD(ifreq) = PRD(ifreq) - CRD(ifreq)*Jrad(iQ,K)

        end do ! Output frequencies
      end do ! Output frequencies ranges

      ! For each K'
      do K1=0,2

        ! Get real value
        rK1 = dble(K1)

        ! Racah algebra
        f1tmp = fun6j(rK1,rJu,rJu1,rJf,1d0,1d0,Flgsg)

        ! If forbidden (6j-sym=0), skip
        if(abs(f1tmp).lt.TINYJS) cycle

        f1tmp = sqrt(2d0*rK1+1d0)*f1tmp

        do iQ1=-K1,K1

          Q1 = dble(iQ1)

          iQl = iQ + iQ1

          Ql = dble(iQl)

          ! Determine the limits in K
          Kmin = max(abs(iQl),nint(abs(rK-rK1)), &
                     nint(abs(rJl-rJl1)))
          Kmax = min(nint(rJl+rJl1),nint(rK+rK1),Kcut)

          rhoc = cZero

          do Kl=Kmin,Kmax

            ! Check sum K + Kl
            if((Kl+K).gt.Kcut) cycle

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

            if (abs(fKfj).lt.TINYJS) cycle

            ! Accumulate in the sum
            rhoc = fKfj*Atom%crho(iR,iz) + rhoc

          end do ! Kl

          ! If no rhoKQ, skip
          if(abs(rhoc).lt.TINYER) cycle

          rhoc = rhoc*ftmp*f1tmp

          ! For each output frequency
          do iran=1,Fin%nran
            do ifreq=Fin%if0(iran),Fin%if1(iran)

              tmp(:,ifreq) = tmp(:,ifreq) + &
                             rhoc*TBout(:,iQ1,K1)*PRD(ifreq)

            end do ! Output frequencies
          end do ! Output frequency ranges

        end do !Q'
      end do !K'

          end do ! Q
        end do ! K

        ! Update redistribution initial index
        kwfreq0 = kwfreq0 + nmfreq

        ! Apply hanle factor and Einstein coefficient
        daux = (2d0*rLl+1d0)*Atom%Ecoeff(iterml,itermu)
        eps20 = eps20 + dble(tmp(0,:)/hanleden)*daux
        eps21 = eps21 + dble(tmp(1,:)/hanleden)*daux
        eps22 = eps22 + dble(tmp(2,:)/hanleden)*daux
        eps23 = eps23 + dble(tmp(3,:)/hanleden)*daux

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
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_warr2)) nullify(p_warr2)

      return

      end subroutine emiss2ordNB_AA

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the second order emission coefficient in absence
      !! of magnetic fields.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!         emerging(logical): Indicates if emergence solution\n
      !!                vx(dfloat): Velocity vector along X\n
      !!                vy(dfloat): Velocity vector along Y\n
      !!                vz(dfloat): Velocity vector along Z\n
      !!          omega(dfloat(:)): Frequency array\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!    Fin(Frequencyc2_class): Structure with the input frequency
      !!                            information\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!       Norma(Nindex_class): Normalization factors for Voigt\n
      !!                            profiles or Voigt profiles\n
      !!             jdir(integer): Output direction in scattering
      !!                            indexing\n
      !!            jtran(integer): Output transition index\n
      !!           itermu(integer): Upper term of output transition\n
      !!           itermf(integer): Lower term of output transition\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!               DwT(dfloat): Thermal part of Doppler width\n
      !!                Dw(dfloat): Doppler width output transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!               vmi(dfloat): Microturbulent velocity\n
      !!      TBout(dcmplx(:,:,:)): Geometrical tensor in the output
      !!                            direction\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!        aprof(dcmplx(:,:)): Absorption profiles\n
      !!          eps20(dfloat(:)): Intensity emissivity\n
      !!          eps21(dfloat(:)): Q emissivity\n
      !!          eps22(dfloat(:)): U emissivity\n
      !!          eps23(dfloat(:)): V emissivity
      subroutine emiss2ordNB_AD(Atom,Geom,emerging,vx,vy,vz,omega, &
                                Red,Fin,Flgsg,Norma,jdir,jtran, &
                                itermu,itermf,iz,if0,if1,DwT,Dw, &
                                vfac,vmi,TBout,Stokes,JKQ,aprof, &
                                eps20,eps21,eps22,eps23)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Frequencyc2_class), intent(inout):: Fin
      type(Nindex_class), intent(in):: Norma
      type(Redc2_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      logical, intent(in):: emerging
      integer, intent(in):: jtran,itermu,itermf,iz,if0,if1,jdir
      double precision, intent(in):: DwT,Dw,vfac,vmi,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(if0:if1), intent(out):: eps20
      double precision, dimension(if0:if1), intent(out):: eps21
      double precision, dimension(if0:if1), intent(out):: eps22
      double precision, dimension(if0:if1), intent(out):: eps23
      complex(kind=8), dimension(:,:), intent(in):: aprof
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TBout

      ! Local

      logical:: PRDc,LPRAM,cohIn
      logical, dimension(Geom%nScatt):: skip_scatt

      integer:: ith1,iph1,i,iterml,itran,iti,icom,ios
      integer:: mF,iU,iU1,iL,iL1,iifreq,iran,ifreq,jfreq,iR
      integer:: K,iQ,K1,iQ1,iQl,Kmin,Kmax,Kl,ktran
      integer:: jjfreq,jjfreq0,kkfreq,kkfreq0,llfreq0,nmfreq
      integer:: kwfreq0,kkfreq0b,indF,indU,indU1,indL,indL1
      integer:: ish1,nfs,nskip,nScatt,stype
#ifdef _OPENMP
      integer:: mmfreq,iidir,tid
#endif
      integer, dimension(:), allocatable:: i_scatt

      double precision:: rLu,rLl,rLf,S,rJu,rJu1,rJl,rJl1,rJf
      double precision:: rK,Q,rK1,Q1,Ql,rKl,wlf,vfac1
      double precision:: eu,eu1,el,el1,ef,au,af,al,auf,aul,at
      double precision:: f61,f62,f63,f64,f0tmp,ftmp,f1tmp,fKfj
      double precision:: vfacw,Dw1,Dfreqw,Norme0,Norme1,rep,imp
      double precision:: omegao,omegai,daux,cost,sint,cosc,sinc
     !double precision:: dNorme2
      double precision, dimension(0:3):: StokesM
      double precision, dimension(:,:), allocatable:: Stokesin

      complex(kind=8):: hanleden,prof,rhoc,intgr
      complex(kind=8):: Norme2,PRDin
      complex(kind=8), dimension(if0:if1):: PRD,CRD0,CRD
      complex(kind=8), dimension(0:3,if0:if1):: tmp
      complex(kind=8), dimension(-2:2,0:2):: Jrad
      complex(kind=8), dimension(:), allocatable, target:: Warr2
      complex(kind=8), dimension(:), allocatable:: Warr2xW
      complex(kind=8), dimension(:), allocatable:: intergrin
#ifdef _OPENMP
      complex(kind=8), dimension(:,:), allocatable:: PRDdir
#endif

      ! Pointers
      type(Frequencyd_class), pointer:: p_frec
      type(Redd_class), pointer:: p_red
      integer, pointer:: p_mfreq
      complex(kind=8), dimension(:), pointer:: p_warr2
      complex(kind=8), dimension(:), pointer:: p_JKQ


      ! Routine name
      urou = 'emiss2ordNB_AD'

      ! Initialize pointers
      nullify(p_frec)
      nullify(p_red)
      nullify(p_mfreq)
      nullify(p_warr2)
      nullify(p_JKQ)

#ifdef _OPENMP
      ! Allocate auxiliar variable to deal with OpenMP
      allocate(PRDdir(Geom%nTh*Geom%nPh2,if0:if1))
#endif


      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      at = (au + af + auf)/Dw

      ! Units normalization factor for CRD profile
      Norme0 = (1d5/Dw)*.5d0*sqrt(IPI)

      ! Doppler shift in doppler units
      vfacw = vfac/Dw

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLf = Atom%rLval(itermf)

      ! Trano index
      ktran = Atom%itrano(jtran)

      !
      ! Initialize the emission coefficient
      !
      eps20 = 0d0
      eps21 = 0d0
      eps22 = 0d0
      eps23 = 0d0
      PRDc = .True.

#ifdef _OPENMP
#else
      ! If there are frequencies
      if (Fin%mxfreq.gt.0) then
        allocate(Warr2xW(Fin%mxfreq))
        allocate(intergrin(Fin%mxfreq))
      end if
#endif


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible lower terms
      do i=1,Atom%nMulti-1

        ! If there is no transition or this term is larger
        ! than the upper term of the output transition, skip
        if(i.ge.itermu.or.Atom%irad(i,itermu).eq.0) cycle

        ! Store the input lower term index
        iterml = i

        ! Get index of input transition
        itran = Atom%irad(iterml,itermu)

        ! Get index of input transition in structure
        ios = -1
        do iti=1,Atom%trano(ktran)%nt
          if (Atom%trano(ktran)%ind(iti).eq.itran) then
            ios = 1
            exit
          end if
        end do
        if (ios.lt.0) cycle

        ! Check if forward
        if (jtran.eq.itran.and. &
            Geom%V_CScatt(1).ge.1d0) then
          nfs = 1
        else
          nfs = 0
        end if

        ! If PRAM, point to the redistribution subblock
        if (PRAM) then

          p_red => Red%trani(iti)
          LPRAM = PRAM.and.p_red%RAM

        ! If not, nothing stored
        else

          LPRAM = .False.

        end if

        ! If not storing or dynamic in quadrature
        if ((.not.LPRAM.or.dyn).and..not.emerging) then

          ! Copy restricted
          skip_scatt = Geom%skip_ksc
          nskip = Geom%nskip
          nScatt = Geom%nScatt - nskip
          i_scatt = Geom%k_scatt

        ! Storing, static, or LOS
        else

          ! Copy total
          skip_scatt = Geom%skip_jsc
          nskip = 0
          nScatt = Geom%nScatt
          i_scatt = Geom%j_scatt

        end if

        ! Point to input transition
        p_frec => Fin%trani(iti)

        ! Get interpolated intensity
        call getStkin(Geom,p_frec,Fin,omega,vx,vy,vz,jdir, &
                      nfs,Stokesin,Stokes)

        ! Get the 'flat' JKQ for this input transition
        JRad = JKQ(:,:,itran)

        ! Redistribution size
        nmfreq = sum(p_frec%mfreq)*(nScatt-nfs)

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

          ! Get indexes
          indF = Atom%irho(itermf)%irho_ij(mF)

          ! Get Jf
          rJf = Atom%rJval(mF,itermf)

          ! For each Ju
          do iU=1,Atom%nJ(itermu)

            ! Get eigenvalue upper level
            eu = Atom%FSfreq(iU,itermu) - Atom%TRfreq(itermu)

            ! Get indexes
            indU = Atom%irho(itermu)%irho_ij(iU)

            ! Get Ju
            rJu = Atom%rJval(iU,itermu)

            f61 = fun6j(rJu,rJf,1d0,rLf,rLu,S,Flgsg)

            if (abs(f61).lt.TINYJS) cycle

            f61 = f61*(2d0*rJu+1d0)*(2d0*rJf+1d0)


            !
            ! Flat contribution. Implicit branching
            !

            ! If in file
            if (vpfil) then

!$omp parallel workshare default(none) &
!$omp shared(CRD0,if0,if1,Norme0,aprof,Atom,jtran,mF,iU)
              CRD0(if0:if1) = Norme0*conjg(aprof(:, &
                                     Atom%i_Vind(jtran)%indNB(mF,iU)))
!$omp end parallel workshare

            ! If stored
            else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(CRD0,if0,if1,Norme0,Norma,mF,iU)
              CRD0(if0:if1) = Norme0* &
                              conjg(Norma%prof(mF,iU,1,1)%cp(if0:if1))
!$omp end parallel workshare

            ! If not stored
            else

              ! Shift term
              Dfreqw  = (eu  - ef + Atom%Dfreq(jtran))/Dw

              ! Normalization factors
              Norme1  = Norma%Norm(mF,iU,1,1)*Norme0

              ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,Dfreqw,omega,vfacw,at,Norme1,Norme0,CRD0)
              do ifreq=if0,if1

                ! Calculate profile u-f
                call voigt(Dfreqw - omega(ifreq)*vfacw,at,prof)

                ! Normalize
                prof = dcmplx(dble(prof)*Norme1,dimag(prof)*Norme0)

                ! Flat spectrum contribution
                CRD0(ifreq) = conjg(prof)

              end do ! frequencies
!$omp end parallel do

            end if ! Storing Voigt

            ! For each Ju'
            do iU1=1,Atom%nJ(itermu)

              ! Get eigenvalue upper' level
              eu1 = Atom%FSfreq(iU1,itermu) - Atom%TRfreq(itermu)

              ! Get indexes
              indU1 = Atom%irho(itermu)%irho_ij(iU1)

              ! Get Ju'
              rJu1 = Atom%rJval(iU1,itermu)

              f62 = fun6j(rJu1,rJf,1d0,rLf,rLu,S,Flgsg)

              if (abs(f62).lt.TINYJS) cycle

              f62 = f62*(2d0*rJu1+1d0)*Flgsg%sg(nint(rJu1+rJf+1d0))

              !
              ! Flat contribution. Implicit branching
              !

              ! If in file
              if (vpfil) then

!$omp parallel workshare default(none) &
!$omp shared(CRD,CRD0,if0,if1,Norme0,aprof,Atom,jtran,mF,iU1)
                CRD(if0:if1) = CRD0(if0:if1) + Norme0*aprof(:, &
                                     Atom%i_Vind(jtran)%indNB(mF,iU1))
!$omp end parallel workshare

              ! If stored
              else if (Norma%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(CRD,CRD0,if0,if1,Norme0,Norma,mF,iU1)
                CRD(if0:if1) = CRD0(if0:if1) + Norme0*&
                               Norma%prof(mF,iU1,1,1)%cp(if0:if1)
!$omp end parallel workshare

              ! If not stored
              else

                ! Shift term
                Dfreqw = (eu1 - ef + Atom%Dfreq(jtran))/Dw

                ! Normalization factor
                Norme1 = Norma%Norm(mF,iU1,1,1)*Norme0

                ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,Dfreqw,omega,vfacw,at,Norme1,Norme0,CRD,CRD0)
                do ifreq=if0,if1

                  ! Calculate profile u'-f
                  call voigt(Dfreqw - omega(ifreq)*vfacw,at,prof)

                  ! Normalize
                  prof = dcmplx(dble(prof)*Norme1,dimag(prof)*Norme0)

                  ! Flat spectrum contribution
                  CRD(ifreq) = CRD0(ifreq) + prof

                end do ! frequencies
!$omp end parallel do

              end if ! Storing Voigt

              !
              ! Continue with the 2nd order emissivity
              !

              ! For each Jl
              do iL=1,Atom%nJ(iterml)

                ! Get eigenvalue of input lower level
                el = Atom%FSfreq(iL,iterml) - Atom%TRfreq(iterml)

                ! Get indexes
                indL = Atom%irho(iterml)%irho_ij(iL)

                ! Get Jl
                rJl = Atom%rJval(iL,iterml)

                f63 = fun6j(rJu,rJl,1d0,rLl,rLu,S,Flgsg)

                if (abs(f63).lt.TINYJS) cycle

                f63 = f63*sqrt(2d0*rJl+1d0)

                ! For each Jl'
                do iL1=1,Atom%nJ(iterml)

                  ! Get eigenvalue of input lower' level
                  el1 = Atom%FSfreq(iL1,iterml) - Atom%TRfreq(iterml)

                  ! Get indexes
                  indL1 = Atom%irho(iterml)%irho_ij(iL1)

                  ! Get Jl1
                  rJl1 = Atom%rJval(iL1,iterml)

                  f64 = fun6j(rJu1,rJl1,1d0,rLl,rLu,S,Flgsg)

                  if (abs(f64).lt.TINYJS) cycle

                  f64 = f64*sqrt(2d0*rJl1+1d0)

                  ! Hanle factor
                  ! TODO ATTENTION TO THIS
                  hanleden = dcmplx(2d0*(au+auf),eu-eu1)/Dw

                  ! Initialize temporal variable
                  tmp = cZero

        !
        ! Reset indexing
        !

        ! Get the component index
        icom = Atom%trano(ktran)%trani(iti)% &
                    WindNB(indL1,indL,indF,indU1,indU)

        if (LPRAM) then
          PRDc = p_red%iPPRD(icom)
        else
          PRDc = .True.
        end if


        !
        ! Create array of Wfunc
        !

        if (PRDc.and.nmfreq.gt.0) then

          ! Initialize array
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
!$omp private(tid,llfreq,jjfreq0,kkfreq0,iifreq,iran,ifreqp_mfreq) &
!$omp private(omegao,ish1,stype,omp,jfreq,jjfreq,kkfreq,omegai) &
!$omp shared(Atom,Fin,p_frec,omega,vfac,Geom,skip_scatt,itran,jtran) &
!$omp shared(Warr2,Dw,Dw1,el,el1,eu,eu1,ef,al,au,af,aul,auf,IPI42) &
!$omp shared(LPRAM,p_red)

#ifdef _OPENMP
          tid = omp_get_thread_num() + 1
          llfreq = 0
#endif
          ! Initialize frequency index
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
              omegao = omega(ifreq)*vfac - Atom%Dfreq(jtran)

              ! For each scattering angle
              do ish1=1,Geom%nScatt

                ! Non-present scattering angle
                if (skip_scatt(ish1)) cycle

                ! Check forward scattering two-terms
                if (itran.eq.jtran.and. &
                    Geom%V_CScatt(ish1).ge.1d0) cycle

                ! Check backward
                if (Geom%V_SScatt(ish1).le.0d0) then
                  stype = 1
                else
                  stype = 0
                end if

                ! For each input frequency
!$omp do
                do jfreq=1,p_mfreq

                  ! Advance indexes
                  jjfreq = jjfreq0 + jfreq
                  kkfreq = kkfreq0 + jfreq

                  ! Get input frequency
                  omegai = p_frec%omega(jjfreq) - &
                           Atom%Dfreq(itran)

                  ! Calculate redistribution function
                  Warr2(kkfreq) = Wfunc(omegai,omegao, &
                                        Dw,Dw1,el,el1,eu,eu1,ef, &
                                        al,au,af,aul,auf, &
                                        Geom%V_CScatt(ish1), &
                                        Geom%V_SScatt(ish1), &
                                        stype)*IPI42

                end do ! input frequencies
!$omp end do nowait

                ! Update kkfreq0
                kkfreq0 = kkfreq0 + p_mfreq

              end do  ! Scattering angles

              ! Update jjfreq0
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! output frequencies
          end do ! output frequencies ranges
!$omp barrier
!$omp flush (Warr2)

          ! If storing
          if (LPRAM) then
!$omp single
            p_red%iPPRD(icom) = .False.
!$omp end single
!$omp do
            do jfreq=1,nmfreq
              rep = dble(Warr2(jfreq))
              imp = dimag(Warr2(jfreq))
              if (rep.le.1e-30) rep = 0d0
              if (abs(imp).le.1e-30) imp = 0d0
              p_red%Pwarr2(kwfreq0+jfreq) = &
                                          cmplx(real(rep),real(imp))
            end do
!$omp end do
          end if ! If storing
!$omp end parallel
        end if ! Initialized

        f0tmp = f61*f62*f63*f64

        ! For each K
        do K=0,Krad

          ! Get real value
          rK = dble(K)

          ! Racah algebra
          ftmp = f0tmp*sqrt(2d0*rK+1d0)

          do iQ=-K,K

            Q = dble(iQ)

      !
      ! Reset identation
      !

      !
      ! Integral over input frequencies
      !

      ! Difference between l and f energies
      wlf = el - ef

      ! Check if static and coherent
      cohIn = dyn.or.abs(wlf).gt.0d0

! dnorme2 removed from list, commented in declarations
!$omp parallel default(none) &
!$omp private(tid,mmfreq,PRDdir,PRD,p_warr2,jjfreq0,llfreq0,kkfreq0) &
!$omp private(ifreq,iran,ifreq,p_mfreq,iidir,ith1,iph1,ish1,) &
!$omp private(omegai,cost,sint,cosc,sinc,vfac1,StokesM,intgr) &
!$omp private(kkfreq0b,Warr2xW,intergrin,Norme2,PRDin) &
!$omp shared(cZero,LPRAM,nmfreq,p_red,Warr2,axial,Fin,p_frec,Geom) &
!$omp shared(jdir,itran,jtran,omp,cohIn,omega,vfac,wlf,dyn,Stokes) &
!$omp shared(CRD,i_scatt,nfs,Stokesin)

#ifdef _OPENMP

      tid = omp_get_thread_num() + 1
      mmfreq = 0

      ! If there are frequencies
      if (Fin%mxfreq.gt.0) then
        allocate(Warr2xW(Fin%mxfreq))
        allocate(intergrin(Fin%mxfreq))
      end if

!$omp workshare
      ! Initialize 2nd order part
      PRDdir = cZero
      PRD = cZero
!$omp end workshare

#else
      ! Initialize 2nd order part
      PRD = cZero
#endif

      ! If storing Warr2
      if (LPRAM) then
        if (nmfreq.gt.0) then
!$omp single
          allocate(p_warr2(nmfreq))
!$omp end single
!$omp workshare
          p_warr2 = dcmplx(p_red%Pwarr2(kwfreq0+1:kwfreq0+nmfreq))
!$omp end workshare
        end if
      ! If not storing
      else
!$omp single
        if (allocated(Warr2)) p_warr2 => Warr2
!$omp end single
      end if

      ! If axial symmetry
      if (axial) then

        ! Initialize frequency indexes
        jjfreq0 = 0
        llfreq0 = 0
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

              ! For each azimuthal direction
              do iph1=1,Geom%nPh2
#ifdef _OPENMP
                ! Advance direction and thread index
                iidir = iidir + 1
                mmfreq = mmfreq + 1
#endif
                ! Scattering index
                ish1 = Geom%i_scatt(iph1,ith1,jdir)

                ! Special treatment if forward for two terms
                if ((jtran.eq.itran.and. &
                     Geom%V_CScatt(ish1).ge.1d0).or. &
                    (p_mfreq.lt.1)) then
#ifdef _OPENMP
                  ! If multi-threading
                  if (omp) then
                    ! If out of range below
                    if (mmfreq.lt.Fin%oif0(tid)) cycle
                    ! If out of range above
                    if (mmfreq.gt.Fin%oif1(tid)) exit
                  end if
#endif
                  ! Interpolate
                  if (cohIn) then

                    ! Input frequency
                    omegai = omega(ifreq)*vfac - wlf

                    ! If there are dynamics
                    if (dyn) then

                      ! Get director cosine
                      cost = Geom%V_mu(ith1)

                      ! Calculate Doppler shift factor
                      vfac1 = 1d0 - vz*cost

                      ! We will be using the inverse
                      vfac1 = 1d0/vfac1

                      ! Shift
                      omegai = omegai*vfac1

                    end if

                    ! Get Stokes
                    StokesM = getStkinnu(omega, &
                                         Stokes(:,:,1,ith1), &
                                         ifreq,omegai)

                  ! Fully coherent
                  else

                    StokesM = Stokes(:,ifreq,1,ith1)

                  end if

                  ! Sum over Stokes parameters, integrand
                  intgr = sum(Geom%TS(:,iQ,K,iph1,ith1)* &
                              StokesM)
#ifdef _OPENMP
                  PRDdir(iidir,ifreq) = PRDdir(iidir,ifreq) + &
                               CRD(ifreq)* &
                               intgr*Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)
#else
                  PRD(ifreq) = PRD(ifreq) + CRD(ifreq)* &
                               intgr*Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)
#endif

                ! Non-forward 2-term scattering
                else
#ifdef _OPENMP
                  ! If multi-threading
                  if (omp) then
                    ! If out of range below
                    if (mmfreq.lt.Fin%oif0(tid)) cycle
                    ! If out of range above
                    if (mmfreq.gt.Fin%oif1(tid)) exit
                  end if
#endif
                  !
                  ! Find initial index for kkfreq

                  ! Scattering index
                  ish1 = i_scatt(ish1)

                  ! Shift in indexes
                  kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                  ! Multiply Warr2 and weights
                  Warr2xW(1:p_mfreq) = p_warr2(kkfreq0b+1: &
                                               kkfreq0b+p_mfreq)* &
                                       p_frec%W_freq(llfreq0+1: &
                                                     llfreq0+p_mfreq)

                  ! Sum Stokes
                  intergrin(1:p_mfreq) = &
                        Stokesin(0,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(0,iQ,K,iph1,ith1) + &
                        Stokesin(1,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(1,iQ,K,iph1,ith1) + &
                        Stokesin(2,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(2,iQ,K,iph1,ith1) + &
                        Stokesin(3,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(3,iQ,K,iph1,ith1)

                  ! Compute norm
                  Norme2 = sum(Warr2xW(1:p_mfreq))

                  ! Integrate
                  PRDin = sum(Warr2xW(1:p_mfreq)* &
                              intergrin(1:p_mfreq))

                  
                  ! Normalize to the first order profile
                  ! and add the directional weights
#ifdef _OPENMP
                  PRDdir(iidir,ifreq) = PRDdir(iidir,ifreq) + &
                               PRDin*CRD(ifreq)* &
                               Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)/Norme2
#else
                  PRD(ifreq) = PRD(ifreq) + &
                               PRDin*CRD(ifreq)*&
                               Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)/Norme2
#endif
                  ! This old normalization (only real part) was
                  ! problematic for non-axial cases with
                  ! angle-dependent
                 !! Normalize real part to the first order
                 !! profile and add the directional
                 !! weights
                 !dNorme2 = dble(Norme2)
                 !if (dNorme2.gt.0d0) then
                 !  rep = dble(PRDin)*dble(CRD(ifreq))/dNorme2
                 !  imp = dimag(PRDin)
!ifdef _OPENMP
                 !  PRDdir(iidir,ifreq) = PRDdir(iidir,ifreq) + &
                 !               dcmplx(rep,imp)*Geom%W_mu(ith1)* &
                 !                               Geom%W_mux2(iph1)
!else
                 !  PRD(ifreq) = PRD(ifreq) + &
                 !               dcmplx(rep,imp)*Geom%W_mu(ith1)* &
                 !                               Geom%W_mux2(iph1)
!endif
                 !else
!ifdef _OPENMP
                 !  PRDdir(iidir,ifreq) = PRDdir(iidir,ifreq) + &
                 !               dcmplx(0d0,dimag(PRDin)* &
                 !                          Geom%W_mu(ith1)* &
                 !                          Geom%W_mux2(iph1))
!else
                 !  PRD(ifreq) = PRD(ifreq) + &
                 !               dcmplx(0d0,dimag(PRDin)* &
                 !                          Geom%W_mu(ith1)* &
                 !                          Geom%W_mux2(iph1))
!endif
                 !end if ! Positive norm

                end if ! Type of scattering

              end do ! azimuthal nodes

              ! Update jjfreq
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! polar nodes

            ! Update llfreq and kkfreq
            llfreq0 = llfreq0 + p_mfreq
            kkfreq0 = kkfreq0 + p_mfreq*(nScatt - nfs)

          end do ! output frequencies
        end do ! output frequencies ranges

      ! If not axial symmetric
      else

        ! Initialize indexes
        jjfreq0 = 0
        llfreq0 = 0
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
                mmfreq = mmfreq + 1
#endif
                ! Scattering index
                ish1 = Geom%i_scatt(iph1,ith1,jdir)

                ! Special treatment if forward for two
                ! terms
                if ((jtran.eq.itran.and. &
                     Geom%V_CScatt(ish1).ge.1d0).or. &
                    (p_mfreq.lt.1)) then
#ifdef _OPENMP
                  ! If multi-threading
                  if (omp) then
                    ! If out of range below
                    if (mmfreq.lt.Fin%oif0(tid)) cycle
                    ! If out of range above
                    if (mmfreq.gt.Fin%oif1(tid)) exit
                  end if
#endif
                  ! Interpolate
                  if (cohIn) then

                    ! Input frequency
                    omegai = omega(ifreq)*vfac - wlf

                    ! If there are dynamics
                    if (dyn) then

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

                    StokesM = getStkinnu(omega, &
                                         Stokes(:,:,iph1,ith1), &
                                         ifreq,omegai)

                  ! Full coherent
                  else

                    StokesM = Stokes(:,ifreq,iph1,ith1)

                  end if

                  ! Sum over Stokes parameters, integrand
                  intgr = sum(Geom%TS(:,iQ,K,iph1,ith1)* &
                              StokesM)

#ifdef _OPENMP
                  PRDdir(iidir,ifreq) = PRDdir(iidir,ifreq) + &
                               CRD(ifreq)* &
                               intgr*Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)
#else
                  PRD(ifreq) = PRD(ifreq) + CRD(ifreq)* &
                               intgr*Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)
#endif

                ! Non-forward 2-term scattering
                else
#ifdef _OPENMP
                  ! If multi-threading
                  if (omp) then
                    ! If out of range below
                    if (mmfreq.lt.Fin%oif0(tid)) then
                      jjfreq0 = jjfreq0 + p_mfreq
                      cycle
                    end if
                    ! If out of range above
                    if (mmfreq.gt.Fin%oif1(tid)) exit
                  end if
#endif
                  ! Scattering index
                  ish1 = i_scatt(ish1)

                  ! Shift in indexes
                  kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                  ! Multiply Warr2 and weights
                  Warr2xW(1:p_mfreq) = p_warr2(kkfreq0b+1: &
                                               kkfreq0b+p_mfreq)* &
                                       p_frec%W_freq(llfreq0+1: &
                                                     llfreq0+p_mfreq)

                  ! Sum Stokes
                  intergrin(1:p_mfreq) = &
                        Stokesin(0,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(0,iQ,K,iph1,ith1) + &
                        Stokesin(1,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(1,iQ,K,iph1,ith1) + &
                        Stokesin(2,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(2,iQ,K,iph1,ith1) + &
                        Stokesin(3,jjfreq0+1:jjfreq0+p_mfreq)* &
                        Geom%TS(3,iQ,K,iph1,ith1)

                  ! Compute norm
                  Norme2 = sum(Warr2xW(1:p_mfreq))

                  ! Integrate
                  PRDin = sum(Warr2xW(1:p_mfreq)* &
                              intergrin(1:p_mfreq))

                  ! Normalize to the first order profile
                  ! and add the directional weights
                  PRD(ifreq) = PRD(ifreq) + &
                               PRDin*CRD(ifreq)* &
                               Geom%W_mu(ith1)* &
                               Geom%W_mux2(iph1)/Norme2

                  ! This old normalization (only real part) was
                  ! problematic for non-axial cases with
                  ! angle-dependent
                 !! Normalize real part to the first order profile
                 !! and add the directional weights
                 !dNorme2 = dble(Norme2)
                 !if (dNorme2.gt.0d0) then
                 !  rep = dble(PRDin)*dble(CRD(ifreq))/dNorme2
                 !  imp = dimag(PRDin)
                 !  PRD(ifreq) = PRD(ifreq) + &
                 !               dcmplx(rep,imp)*Geom%W_mu(ith1)* &
                 !                               Geom%W_mux2(iph1)
                 !else
                 !  PRD(ifreq) = PRD(ifreq) + &
                 !               dcmplx(0d0,dimag(PRDin)* &
                 !                          Geom%W_mu(ith1)* &
                 !                          Geom%W_mux2(iph1))
                 !end if ! Valid norm

                  ! Update indexes
                  jjfreq0 = jjfreq0 + p_mfreq

                end if ! Type of scattering

              end do ! azimuthal nodes
            end do ! polar nodes

            ! Advance
            llfreq0 = llfreq0 + p_mfreq
            kkfreq0 = kkfreq0 + p_mfreq*(nScatt - nfs)

          end do ! output frequencies
        end do ! output frequencies ranges

      end if ! Axial symmetry
!$omp end parallel

#ifdef _OPENMP
      ! for each output frequency
      do iran=1,Fin%nran
        do ifreq=Fin%if0(iran),Fin%if1(iran)

          ! get proper prd quantity
          PRD(ifreq) = PRD(ifreq) + sum(PRDdir(:,ifreq))

        end do ! output frequencies
      end do ! output frequency range
#endif

      ! Clean p_warr2
      if (LPRAM.and.nmfreq.gt.0) deallocate(p_warr2)
      nullify(p_warr2)

      !
      ! Flat spectrum contribution. Implicit branching ratio
      !

      ! For each output frequency
      do iran=1,Fin%nran
        do ifreq=Fin%if0(iran),Fin%if1(iran)

          ! Substract the flat spectrum part due to just
          ! radiative excitation
          PRD(ifreq) = PRD(ifreq) - CRD(ifreq)*Jrad(iQ,K)

        end do ! Output frequencies
      end do ! Output frequencies ranges

      ! For each K'
      do K1=0,2

        ! Get real value
        rK1 = dble(K1)

        ! Racah algebra
        f1tmp = fun6j(rK1,rJu,rJu1,rJf,1d0,1d0,Flgsg)

        ! If forbidden (6j-sym=0), skip
        if(abs(f1tmp).lt.TINYJS) cycle

        f1tmp = sqrt(2d0*rK1+1d0)*f1tmp

        do iQ1=-K1,K1

          Q1 = dble(iQ1)

          iQl = iQ + iQ1

          Ql = dble(iQl)

          ! Determine the limits in K
          Kmin = max(abs(iQl),nint(abs(rK-rK1)), &
                     nint(abs(rJl-rJl1)))
          Kmax = min(nint(rJl+rJl1),nint(rK+rK1),Kcut)

          rhoc = cZero

          do Kl=Kmin,Kmax

            ! Check sum K + Kl
            if((Kl+K).gt.Kcut) cycle

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

            if (abs(fKfj).lt.TINYJS) cycle

            ! Accumulate in the sum
            rhoc = fKfj*Atom%crho(iR,iz) + rhoc

          end do ! Kl

          ! If no rhoKQ, skip
          if(abs(rhoc).lt.TINYER) cycle

          rhoc = rhoc*ftmp*f1tmp

          ! For each output frequency
          do iran=1,Fin%nran
            do ifreq=Fin%if0(iran),Fin%if1(iran)

              tmp(:,ifreq) = tmp(:,ifreq) + &
                             rhoc*TBout(:,iQ1,K1)*PRD(ifreq)

            end do ! Output frequencies
          end do ! Output frequency ranges

        end do !Q'
      end do !K'

          end do ! Q
        end do ! K

        ! Update redistribution initial index
        kwfreq0 = kwfreq0 + nmfreq

        ! Apply hanle factor and Einstein coefficient
        daux = (2d0*rLl+1d0)*Atom%Ecoeff(iterml,itermu)
        eps20 = eps20 + dble(tmp(0,:)/hanleden)*daux
        eps21 = eps21 + dble(tmp(1,:)/hanleden)*daux
        eps22 = eps22 + dble(tmp(2,:)/hanleden)*daux
        eps23 = eps23 + dble(tmp(3,:)/hanleden)*daux

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
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_warr2)) nullify(p_warr2)
      if (associated(p_JKQ)) nullify(p_JKQ)
      if (allocated(Warr2xW)) deallocate(Warr2xW)
      if (allocated(intergrin)) deallocate(intergrin)

      return

      end subroutine emiss2ordNB_AD

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient for a LTE line.\n
      !!       line(LTEline_class): Structure with the LTE line data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!      aprof(LTEprof_class): Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!         Bstrength(dfloat): Magnetic field strength\n
      !!                pE(dfloat): Unit transformation factor\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity
      subroutine absorbLTE(line,TB,omega,Flgsg,iz,if0,if1,aprof,Dw, &
                           vfac,Bstrength,pE, &
                           eta0,eta1,eta2,eta3,rha1,rha2,rha3)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(LTEprof_class), intent(in):: aprof
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw, pE, vfac, Bstrength
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0,eta1
      double precision, dimension(if0:if1), intent(out):: eta2,eta3
      double precision, dimension(if0:if1), intent(out):: rha1,rha2
      double precision, dimension(if0:if1), intent(out):: rha3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      integer:: ifreq,iMu,iMl,iq

      double precision:: at,feta,Dfreqw,dnubw,vfacw
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


      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE/Dw

      ! Line quantities

      ! Damping parameter
      at = line%damp(iz)/Dw

      ! Energy
      Dfreqw = (line%Eu - line%El)/Dw

      ! Shift
      vfacw = vfac/Dw

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
          tK0 = TB(0,0,0)*fun3j(1d0,1d0,0d0,q,-q,0d0,Flgsg)

          ! K = 1
          tK3 = sqrt3*TB(3,0,1)*fun3j(1d0,1d0,1d0,q,-q,0d0,Flgsg)

          ! K = 2
          tK2 = sqrt5*fun3j(1d0,1d0,2d0,q,-q,0d0,Flgsg)
          tK0 = tK0 + tK2*TB(0,0,2)
          tK1 = tK2*TB(1,0,2)
          tK2 = tK2*TB(2,0,2)

          ! Scale
          dtK0 = dble(tK0)*ftmp
          dtK1 = dble(tK1)*ftmp
          dtK2 = dble(tK2)*ftmp
          dtK3 = dble(tK3)*ftmp

          ! If stored in RAM
          if (aprof%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(aprof,eta0,eta1,eta2,eta3,rha1,rha2,rha3) &
!$omp shared(dtK0,dtK1,dtK2,dtK3)
            eta0 = eta0 + dtK0*dble(aprof%comp(iMl,iMu)%cp)
            eta1 = eta1 + dtK1*dble(aprof%comp(iMl,iMu)%cp)
            eta2 = eta2 + dtK2*dble(aprof%comp(iMl,iMu)%cp)
            eta3 = eta3 + dtK3*dble(aprof%comp(iMl,iMu)%cp)
            rha1 = rha1 + dtK1*dimag(aprof%comp(iMl,iMu)%cp)
            rha2 = rha2 + dtK2*dimag(aprof%comp(iMl,iMu)%cp)
            rha3 = rha3 + dtK3*dimag(aprof%comp(iMl,iMu)%cp)
!$omp end parallel workshare

          ! Not stored
          else

            ! Magnetic shift
            dnubw = B2LK*Bstrength* &
                    (line%gu*rMu - line%gl*rMl)/Dw

            ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feta,Dfreqw,dnubw,omega,vfacw,at) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3) &
!$omp shared(dtK0,dtK1,dtK2,dtK3)
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreqw + dnubw - omega(ifreq)*vfacw,at,prof)

              eta0(ifreq) = eta0(ifreq) + dtK0*dble(prof)
              eta1(ifreq) = eta1(ifreq) + dtK1*dble(prof)
              eta2(ifreq) = eta2(ifreq) + dtK2*dble(prof)
              eta3(ifreq) = eta3(ifreq) + dtK3*dble(prof)
              rha1(ifreq) = rha1(ifreq) + dtK1*dimag(prof)
              rha2(ifreq) = rha2(ifreq) + dtK2*dimag(prof)
              rha3(ifreq) = rha3(ifreq) + dtK3*dimag(prof)

            end do ! frequencies
!$omp end parallel do

          end if ! Type of profile calculation

        end do ! Ml
      end do ! Mu

      ! Add units
!$omp parallel workshare default(none) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3,feta)
      eta0 = eta0*feta
      eta1 = eta1*feta
      eta2 = eta2*feta
      eta3 = eta3*feta
      rha1 = rha1*feta
      rha2 = rha2*feta
      rha3 = rha3*feta
!$omp end workshare

      return

      end subroutine absorbLTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorption coefficient for a LTE line.\n
      !!       line(LTEline_class): Structure with the LTE line data\n
      !!         TB(dcmplx(:,:,:)): Geometry tensors in magnetic field
      !!                            reference frame\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!              if0(integer): First frequency index for this
      !!                            transition\n
      !!              if1(integer): Last frequency index for this
      !!                            transition\n
      !!      aprof(LTEprof_class): Voigt profiles\n
      !!                Dw(dfloat): Doppler width of transition\n
      !!              vfac(dfloat): Doppler shift factor\n
      !!         Bstrength(dfloat): Magnetic field strength\n
      !!                pE(dfloat): Unit transformation factor\n
      !!           eta0(dfloat(:)): Intensity absorptivity\n
      !!           eta1(dfloat(:)): Q absorptivity\n
      !!           eta2(dfloat(:)): U absorptivity\n
      !!           eta3(dfloat(:)): V absorptivity\n
      !!           rha1(dfloat(:)): Q dichroic absorptivity\n
      !!           rha2(dfloat(:)): U dichroic absorptivity\n
      !!           rha3(dfloat(:)): V dichroic absorptivity
      !!           eps0(dfloat(:)): Intensity emissivity\n
      !!           eps1(dfloat(:)): Q emissivity\n
      !!           eps2(dfloat(:)): U emissivity\n
      !!           eps3(dfloat(:)): V emissivity\n
      !!           rhs1(dfloat(:)): Q 'dichroic' emissivity\n
      !!           rhs2(dfloat(:)): U 'dichroic' emissivity\n
      !!           rhs3(dfloat(:)): V 'dichroic' emissivity
      subroutine rt1ordLTE(line,TB,omega,Flgsg,iz,if0,if1,aprof,Dw, &
                           vfac,Bstrength,pE, &
                           eta0,eta1,eta2,eta3,rha1,rha2,rha3, &
                           eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(LTEprof_class), intent(in):: aprof
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw, pE, vfac, Bstrength
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta0,eta1
      double precision, dimension(if0:if1), intent(out):: eta2,eta3
      double precision, dimension(if0:if1), intent(out):: eps0,eps1
      double precision, dimension(if0:if1), intent(out):: eps2,eps3
      double precision, dimension(if0:if1), intent(out):: rha1,rha2
      double precision, dimension(if0:if1), intent(out):: rhs1,rhs2
      double precision, dimension(if0:if1), intent(out):: rha3,rhs3
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TB

      ! Local

      integer:: ifreq,iMu,iMl,iq

      double precision:: at,feta,feps,Dfreqw,dnubw,vfacw
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


      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE/Dw

      ! Emissivity factor
      feps = line%nu(iz)*1d3*IPI41*line%Aul/Dw

      ! Line quantities

      ! Damping parameter
      at = line%damp(iz)/Dw

      ! Energy
      Dfreqw = (line%Eu - line%El)/Dw

      ! Shift
      vfacw = vfac/Dw

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
          tK0 = TB(0,0,0)*fun3j(1d0,1d0,0d0,q,-q,0d0,Flgsg)

          ! K = 1
          tK3 = sqrt3*TB(3,0,1)*fun3j(1d0,1d0,1d0,q,-q,0d0,Flgsg)

          ! K = 2
          tK2 = sqrt5*fun3j(1d0,1d0,2d0,q,-q,0d0,Flgsg)
          tK0 = tK0 + tK2*TB(0,0,2)
          tK1 = tK2*TB(1,0,2)
          tK2 = tK2*TB(2,0,2)

          ! Scale
          dtK0 = dble(tK0)*ftmp
          dtK1 = dble(tK1)*ftmp
          dtK2 = dble(tK2)*ftmp
          dtK3 = dble(tK3)*ftmp

          ! If stored in RAM
          if (aprof%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(aprof,eta0,eta1,eta2,eta3,rha1,rha2,rha3) &
!$omp shared(dtK0,dtK1,dtK2,dtK3)
            eta0 = eta0 + dtK0*dble(aprof%comp(iMl,iMu)%cp)
            eta1 = eta1 + dtK1*dble(aprof%comp(iMl,iMu)%cp)
            eta2 = eta2 + dtK2*dble(aprof%comp(iMl,iMu)%cp)
            eta3 = eta3 + dtK3*dble(aprof%comp(iMl,iMu)%cp)
            rha1 = rha1 + dtK1*dimag(aprof%comp(iMl,iMu)%cp)
            rha2 = rha2 + dtK2*dimag(aprof%comp(iMl,iMu)%cp)
            rha3 = rha3 + dtK3*dimag(aprof%comp(iMl,iMu)%cp)
!$omp end parallel workshare

          ! Not stored
          else

            ! Magnetic shift
            dnubw = B2LK*Bstrength* &
                    (line%gu*rMu - line%gl*rMl)/Dw

            ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,prof) &
!$omp shared(if0,if1,feta,Dfreqw,dnubw,omega,vfacw,at) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3) &
!$omp shared(dtK0,dtK1,dtK2,dtK3)
            do ifreq=if0,if1

              ! Calculate profile
              call voigt(Dfreqw + dnubw - omega(ifreq)*vfacw,at,prof)

              eta0(ifreq) = eta0(ifreq) + dtK0*dble(prof)
              eta1(ifreq) = eta1(ifreq) + dtK1*dble(prof)
              eta2(ifreq) = eta2(ifreq) + dtK2*dble(prof)
              eta3(ifreq) = eta3(ifreq) + dtK3*dble(prof)
              rha1(ifreq) = rha1(ifreq) + dtK1*dimag(prof)
              rha2(ifreq) = rha2(ifreq) + dtK2*dimag(prof)
              rha3(ifreq) = rha3(ifreq) + dtK3*dimag(prof)

            end do ! frequencies
!$omp end parallel do

          end if ! Type of profile calculation

        end do ! Ml
      end do ! Mu

      ! Add units
!$omp parallel workshare default(none) &
!$omp shared(eta0,eta1,eta2,eta3,rha1,rha2,rha3) &
!$omp shared(eps0,eps1,eps2,eps3,rhs1,rhs2,rhs3) &
!$omp shared(feta,feps)
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
!$omp end workshare

      return

      end subroutine rt1ordLTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the absorptivity due to photoionization\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!       itran(integer): Transition index\n
      !!     ilevell(integer): Lower level index\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!       eta(dfloat(:)): Absorptivity
      subroutine photoabs(Atom,itran,ilevell,iz,if0,if1,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran, ilevell, iz, if0, if1
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq,iterml,iJl,iR

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

      end subroutine photoabs

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
      !!       eta(dfloat(:)): Stimulated emissivity
      subroutine photoeps(Atom,omega,T,ne,itran,ilevelu,iz,if0,if1, &
                          eps,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran, ilevelu, iz, if0, if1
      double precision, intent(in):: T, ne
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta, eps

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

      omega3 = omega(if0:if1)
      omega3 = omega3*omega3*omega3

      ! For each frequency
!$omp parallel do default(none) &
!$omp private(ifreq,exu,pE) &
!$omp shared(eta,eps,if0,if1,c0,c1,omega,omega3,Atom,itran,rhou)
      do ifreq=if0,if1

        exu = c0*omega(ifreq)
        exu = diexp(exu)

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu*rhou

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies
!$omp end parallel do

      return

      end subroutine photoeps

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the emissivity due to recombination with
      !! precomputed quantities\n
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
      !!       eta(dfloat(:)): Stimulated emissivity
      subroutine photoepsS(Atom,omega3,exu,T,ne,itran, &
                           ilevelu,iz,if0,if1,eps,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran, ilevelu, iz, if0, if1
      double precision, intent(in):: T, ne
      double precision, dimension(if0:if1), intent(in):: omega3
      double precision, dimension(if0:if1), intent(in):: exu
      double precision, dimension(if0:if1), intent(out):: eta, eps

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
!$omp parallel do default(none) &
!$omp private(ifreq,pE) &
!$omp shared(eta,eps,if0,if1,c1,omega3,exu,Atom,itran,rhou)
      do ifreq=if0,if1

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu(ifreq)*rhou

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies
!$omp end parallel do

      return

      end subroutine photoepsS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Gets the JKQC for AA and dynamic\n
      !!    Fin(Frequencyc2_class): Structure with the input frequency
      !!                            information\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!               iz(integer): Height index\n
      !!               DwT(dfloat): Thermal part of Doppler width\n
      !!                vx(dfloat): Velocity vector along X\n
      !!                vy(dfloat): Velocity vector along Y\n
      !!                vz(dfloat): Velocity vector along Z\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!     JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!    JRadC(dcomplex(:,:,:)): JKQ comoving frame
      subroutine getJKQstar(Fin,Geom,iz,DwT,vx,vy,vz,omega, &
                            Flgsg,Stokes,JKQa,JradC)

      ! I/O
      type(Frequencyc2_class), intent(in):: Fin
      type(Geometry_class), intent(in):: Geom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz
      double precision, intent(in):: DwT,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(:,:,:), allocatable:: JradC

      ! Local
      logical:: shift,asym
      integer:: ifreq,ith1,K,iph1,iQ,jz
      double precision:: vfac1,omegao,cost,sint,cosc,sinc
      double precision, dimension(0:3):: StokesM

      ! Pointers

      ! If velocity is below threshold
      shift = (vx*vx + vy*vy + vz*vz)*1d6*c.ge.vrfrac*DwT
      asym = size(JKQa).gt.10

      ! Allocate
      allocate(JradC(-2:2,0:2,Fin%ggf0:Fin%ggf1))

      ! For each frequency, get mean intensity
!$omp parallel default(none) &
!$omp private(ifreq,ith1,cost,vfac1,omegao,StokesM,K,sint,iph1) &
!$omp private(cosc,sinc,iQ,jz) &
!$omp shared(Fin,Geom,shift,axial,vx,vy,vz,omega,Stokes,JradC,iz) &
!$omp shared(JKQa,asym,force_asym,Rz0)

      ! Initialize
      vfac1 = 1d0
!$omp workshare
      JradC = cZero
!$omp end workshare

      !
      ! If there are ad-hoc asymmetries
      !
      if (asym) then

        ! If forcing the ad-hoc asymmetry only anisotropy in vertical
        ! frame
        if (force_asym.or.axial) then
!$omp do
          ! Relevant frequencies
          do ifreq=Fin%ggf0,Fin%ggf1

            ! For each polar direction
            do ith1=1,Geom%nTh

              ! Get director cosines
              if (shift) cost = Geom%V_mu(ith1)

              ! If axial
              if (axial) then

                ! Calculate Doppler shift factor
                if (shift) then

                  ! Interpolate shift
                  vfac1 = 1d0 + vz*cost
                  omegao = omega(ifreq)*vfac1
                  StokesM = getStkinnu(omega, &
                                       Stokes(:,:,1,ith1), &
                                       ifreq,omegao)

                ! No shift
                else

                  ! Value
                  StokesM = Stokes(:,ifreq,1,ith1)

                end if ! shift

                ! Add to JKQ
                do K=0,Krad
                  JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                     Geom%W_mu(ith1)* &
                                   sum(StokesM* &
                                       Geom%TS(:,0,K,1,ith1))
                end do

              ! Non-axial
              else

                if (shift) sint = sqrt(1d0 - cost*cost)

                ! For each azimuth
                do iph1=1,Geom%nPh

                  ! Calculate Doppler shift factor
                  if (shift) then

                    ! Interpolate shift
                    cosc = Geom%v_mux(iph1)
                    sinc = Geom%v_muy(iph1)* &
                           sqrt(1d0 - cosc*cosc)
                    vfac1 = 1d0 + vx*sint*cosc + &
                                  vy*sint*sinc + &
                                  vz*cost
                    omegao = omega(ifreq)*vfac1
                    StokesM = getStkinnu(omega, &
                                         Stokes(:,:,iph1,ith1), &
                                         ifreq,omegao)

                  ! No shift
                  else

                    ! Value
                    StokesM = Stokes(:,ifreq,iph1,ith1)

                  end if ! shift

                  ! Add to JKQ
                  do K=0,Krad
                    JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                       Geom%W_mu(ith1)* &
                                       Geom%W_mux(iph1)* &
                                       sum(StokesM* &
                                           Geom%TS(:,0,K,iph1,ith1))
                  end do
                end do

              end if

            end do ! Polar directions
          end do ! Relevant frequencies
!$omp end do
        ! Additive ad-hoc asymmetry
        else
!$omp do
          ! Relevant frequencies
          do ifreq=Fin%ggf0,Fin%ggf1

            ! For each polar direction
            do ith1=1,Geom%nTh

              ! Compute z shift
              if (shift) then
                cost = Geom%V_mu(ith1)
                sint = sqrt(1d0 - cost*cost)
              end if

              ! For each azimuth
              do iph1=1,Geom%nPh

                ! Calculate Doppler shift factor
                if (shift) then

                  ! Interpolate shift
                  cosc = Geom%v_mux(iph1)
                  sinc = Geom%v_muy(iph1)* &
                         sqrt(1d0 - cosc*cosc)
                  vfac1 = 1d0 + vx*sint*cosc + &
                                vy*sint*sinc + &
                                vz*cost
                  omegao = omega(ifreq)*vfac1
                  StokesM = getStkinnu(omega, &
                                       Stokes(:,:,iph1,ith1), &
                                       ifreq,omegao)

                ! No shift
                else

                  ! Value
                  StokesM = Stokes(:,ifreq,iph1,ith1)

                end if ! Shift

                ! Add to JKQ
                do K=0,Krad
                  do iQ=0,K
                  JradC(iQ,K,ifreq) = JradC(iQ,K,ifreq) + &
                                        Geom%W_mu(ith1)* &
                                        Geom%W_mux2(iph1)* &
                                        sum(StokesM* &
                                     Geom%TS(:,iQ,K,iph1,ith1))
                  end do
                end do
              end do ! Azimuth
            end do ! Polar
          end do ! Frequencies
!$omp end do
        end if ! Type of ad-hoc asymmetry

        ! Current height for JKQa
        jz = iz - Rz0 + 1

        ! Now add ad-hoc assymetry
!$omp do
        ! Relevant frequencies
        do ifreq=Fin%ggf0,Fin%ggf1

        JradC(0:2,1:2,ifreq) = JradC(0:2,1:2,ifreq) + &
                               JKQa(3:5,:,jz)*dble(JradC(0,0,ifreq))

        end do
!$omp end do

      ! No Ad-hoc asymmetries
      else

      ! For each frequency, get JKQ
!$omp do
        do ifreq=Fin%ggf0,Fin%ggf1

          ! For each polar direction
          do ith1=1,Geom%nTh

            ! Get director cosines
            if (shift) cost = Geom%V_mu(ith1)

            ! If axial symmetric
            if (axial) then

              ! Calculate Doppler shift factor
              if (shift) then

                ! Interpolate shift
                vfac1 = 1d0 + vz*cost
                omegao = omega(ifreq)*vfac1
                StokesM = getStkinnu(omega, &
                                     Stokes(:,:,1,ith1), &
                                     ifreq,omegao)

              ! No shift
              else

                ! Value
                StokesM = Stokes(:,ifreq,1,ith1)

              end if ! shift

              ! Add to JKQ
              do K=0,Krad
                JradC(0,K,ifreq) = JradC(0,K,ifreq) + &
                                   Geom%W_mu(ith1)* &
                                   sum(StokesM* &
                                       Geom%TB(:,0,K,1,ith1,iz))
              end do

            ! Non axial symmetric
            else

              if (shift) sint = sqrt(1d0 - cost*cost)

              ! For each azimuth
              do iph1=1,Geom%nPh

                ! Calculate Doppler shift factor
                if (shift) then

                  ! Interpolate shift
                  cosc = Geom%v_mux(iph1)
                  sinc = Geom%v_muy(iph1)* &
                         sqrt(1d0 - cosc*cosc)
                  vfac1 = 1d0 + vx*sint*cosc + &
                                vy*sint*sinc + &
                                vz*cost
                  omegao = omega(ifreq)*vfac1
                  StokesM = getStkinnu(omega, &
                                       Stokes(:,:,iph1,ith1), &
                                       ifreq,omegao)


                ! No shift
                else

                  ! Value
                  StokesM = Stokes(:,ifreq,iph1,ith1)

                end if ! Shift

                ! Add to JKQ
                do K=0,Krad
                  do iQ=0,K
                    JradC(iQ,K,ifreq) = JradC(iQ,K,ifreq) + &
                                          Geom%W_mu(ith1)* &
                                          Geom%W_mux2(iph1)* &
                                          sum(StokesM* &
                                       Geom%TB(:,iQ,K,iph1,ith1,iz))
                  end do
                end do

              end do ! Azimuth

            end if ! Axial symmetry

          end do ! Polar
        enddo ! Frequencies
!$omp end do
      end if ! Are there ad-hoc asymmetries
!$omp end parallel


      ! If not axial 
      if (.not.axial) then

        ! Complete negative Q
        do K=1,Krad
          do iQ=0,K
            JradC(-iQ,K,:) = Flgsg%sg(iQ)*conjg(JradC(iQ,K,:))
          end do
        end do

      end if

      end subroutine getJKQstar

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolates Stokes parameters to the requested frequency\n
      !!     omega(dfloat(:)): Frequency array\n
      !!  Stokes(dfloat(:,:)): Stokes parameters\n
      !!       ifreq(integer): Frequency index of the output frequency
      !!                       associated to the requested input
      !!                       frequency\n
      !!            x(dfloat): Input frequency to interpolate into
      function getStkinnu(omega,Stokes,ifreq,x)

      ! I/O
      integer, intent(in):: ifreq
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
        if (x.ge.omega(nfreq)-TINYO) then

          getStkinnu = Stokes(:,nfreq)

          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,nfreq-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getStkinnu = Stokes(:,jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq).and. &
                    x.lt.omega(jfreq+1)) then

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
        if (x.le.omega(1)+TINYO) then

          getStkinnu = Stokes(:,1)

          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,2,-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getStkinnu = Stokes(:,jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq-1).and. &
                    x.lt.omega(jfreq)) then

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

      !> Interpolates Stokes parameters to the input frequency axis\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!  p_frec(Frequencyd_class): Input frequency data for output
      !!                            index and input transition\n
      !!    Fin(Frequencyc2_class): Input frequency data for output
      !!                            index\n
      !!          omega(dfloat(:)): Frequency array\n
      !!                vx(dfloat): Velocity vector along X\n
      !!                vy(dfloat): Velocity vector along Y\n
      !!                vz(dfloat): Velocity vector along Z\n
      !!             jdir(integer): Output direction index\n
      !!              nfs(integer): Number of forward scattering\n
      !!        Stkin(double(:,:)): Interpolated Stokes\n
      !!      Stk(double(:,:,:,:)): Original Stokes
      subroutine getStkin(Geom,p_frec,Fin,omega,vx,vy,vz,jdir, &
                          nfs,Stkin,Stk)

      ! I/O
      type(Geometry_class), intent(in):: Geom
      type(Frequencyd_class), intent(in), pointer:: p_frec
      type(Frequencyc2_class), intent(in):: Fin
      integer, intent(in):: jdir,nfs
      double precision, intent(in):: vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:,:,:,:), intent(in):: Stk
      double precision, dimension(:,:), &
                        allocatable, intent(out):: Stkin

      ! Local
      integer:: ith1,iph1,ish1
      integer:: nmfreq,iran,ifreq,iifreq,lifreq,ibfreq
      integer:: jfreq,jjfreq0,jjfreq,kkfreq0,kkfreq,nblock

      double precision:: dx,vfac1
      double precision:: cost,sint,cosc,sinc
      double precision, dimension(4):: y0,dy

      integer, pointer:: p_mfreq


      ! Nullify
      nullify(p_mfreq)

      ! Get J size
      nmfreq = 0

      ! Run over all output frequencies
      iifreq = 0
      do iran=1,Fin%nran
        do ifreq=Fin%if0(iran),Fin%if1(iran)

          ! Advance index
          iifreq = iifreq + 1

          ! Input frequency number
          p_mfreq => p_frec%mfreq(iifreq)

          ! Add size
          if (p_mfreq.gt.0) nmfreq = nmfreq + p_mfreq

        end do
      end do

      ! If no size, return
      if (nmfreq.lt.1) return

      ! Allocate intensity
      if (axiali) then
        allocate(Stkin(0:3,nmfreq*Geom%nTh))
      else
        allocate(Stkin(0:3,nmfreq*(Geom%nPh2*Geom%nTh-nfs)))
      end if


      ! If interpolation data
      if (p_frec%RAM) then

        ! Axial
        if (axiali) then

!$omp parallel default(none) &
!$omp private(iifreq,iran,ifreq,jfreq,p_mfreq,y0,dx,dy) &
!$omp private(ith1,jjfreq,jjfreq0) &
!$omp shared(Fin,Geom,p_frec,dyn)

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

              ! If coherent wing
              if (p_mfreq.lt.1) cycle

              ! For each polar direction
              do ith1=1,Geom%nTh

                ! For each input frequency
!$omp do
                do jfreq=1,p_mfreq

                  ! Advance index
                  jjfreq = jjfreq0 + jfreq
                  kkfreq = kkfreq0 + jfreq

                  ! Initial
                  y0 = Stk(:,p_frec%index1(jjfreq),1,ith1)

                  ! Slope
                  dy = Stk(:,p_frec%index2(jjfreq),1,ith1) - y0

                  ! Linear interpolation
                  Stkin(:,kkfreq) = dy*p_frec%dx(jjfreq) + y0

                end do ! Input frequencies
!$omp end do
                ! Update jjfreq
                if (dyn) jjfreq0 = jjfreq0 + p_mfreq
                kkfreq0 = kkfreq0 + p_mfreq

              end do ! polar nodes

              ! Update jjfreq
              if (.not.dyn) jjfreq0 = jjfreq0 + p_mfreq

            end do ! output frequencies
          end do ! output frequencies ranges

        ! Not axial
        else

!$omp parallel default(none) &
!$omp private(iifreq,iran,ifreq,jfreq,p_mfreq,y0,dx,dy) &
!$omp private(ith1,iph1,ish1,jjfreq,jjfreq0) &
!$omp shared(Fin,Geom,nfs,p_frec,dyn,jdir)

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

              ! For each polar direction
              do ith1=1,Geom%nTh

                ! For each azimuthal direction
                do iph1=1,Geom%nPh2

                  ! Scattering index
                  ish1 = Geom%i_scatt(iph1,ith1,jdir)

                  ! Special treatment if forward for two terms
                  if ((nfs.eq.1.and.Geom%V_CScatt(ish1).ge.1d0).or. &
                      (p_mfreq.lt.1)) cycle

                  ! For each input frequency
!$omp do
                  do jfreq=1,p_mfreq

                    ! Advance indexes
                    jjfreq = jjfreq0 + jfreq
                    kkfreq = kkfreq0 + jfreq

                    ! Linear interpolation
                    y0 = Stk(:,p_frec%index1(jjfreq),iph1,ith1)
                    dy = Stk(:,p_frec%index2(jjfreq),iph1,ith1) - y0
                    Stkin(:,kkfreq) = dy*p_frec%dx(jjfreq) + y0

                  end do ! Input frequencies
!$omp end do
                  ! Update indexes
                  if (dyn) jjfreq0 = jjfreq0 + p_mfreq
                  kkfreq0 = kkfreq0 + p_mfreq

                end do ! azimuthal nodes
              end do ! polar nodes

              ! Advance
              if (.not.dyn) jjfreq0 = jjfreq0 + p_mfreq

            end do ! output frequencies
          end do ! output frequencies ranges
!$omp end parallel

        end if ! Axial symmetry

      ! No interpolation data
      else

!$omp parallel default(none) &
!$omp private(iifreq,lifreq,iran,ifreq,jfreq,p_mfreq) &
!$omp private(ith1,iph1,cost,sint,cosc,sinc,vfac1) &
!$omp private(jjfreq,kkfreq,jjfreq0,kkfreq0,ibfreq) &
!$omp shared(dyn,axiali,nbth,nbph,vx,vy,vz,omega,Fin,Geom) &
!$omp shared(ffjtran,ffitran,p_frec) &
!$omp shared(nfreq)

        ! If axial
        if (axiali) then

          ! If dynamic
          if (dyn) then

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
                      Stkin(:,kkfreq) = Stk(:,1,1,ith1)

                    ! If out of range, take the value at the
                    ! boundary
                    else if (p_frec%omega(jjfreq)*vfac1.ge. &
                             (omega(nfreq) - TINYO)) then

                      ! We are in the last frequency
                      lifreq = nfreq

                      ! The index to take is nfreq
                      Stkin(:,kkfreq) = Stk(:,nfreq,1,ith1)

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
                          Stkin(:,kkfreq) = Stk(:,lifreq,1,ith1)

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
                          y0 = Stk(:,lifreq,1,ith1)

                          ! The second index is the upper
                          dy = Stk(:,lifreq+1,1,ith1) - y0

                          ! Inverse of the distance
                          ! between the two outputs
                          dx = (p_frec%omega(jjfreq)*vfac1 - &
                                omega(lifreq))/ &
                               (omega(lifreq+1) - omega(lifreq))

                          ! Interpolate
                          Stkin(:,kkfreq) = dx*dy + y0

                          exit

                        end if ! Check output frequency

                      end do ! Run output frequencies

                    end if ! Check if out of limits

                  end do ! Run input frequencies
!$omp end do
                  ! Update indexes
                  kkfreq0 = kkfreq0 + p_mfreq

                end do ! Input polar

                ! Update indexes
                jjfreq0 = jjfreq0 + p_mfreq

              end do ! Output frequencies
            end do ! Output frequency ranges
!$omp end parallel

          ! If static
          else

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

                ! Skip empty
                if (p_mfreq.lt.1) cycle

                ! Reset the search frequency
                lifreq = Fin%ggf0
!$omp do
                ! For each input frequency
                do jfreq=1,p_mfreq

                  ! Advance indexes
                  jjfreq = jjfreq0 + jfreq

                  ! If out of range, take the value at the
                  ! boundary
                  if (p_frec%omega(jjfreq).le. &
                      omega(1)+TINYO) then

                    ! We are still looking in the first one
                    lifreq = 1

                    ! Inclinations
                    do ith1=1,Geom%nTh

                      ! Output index
                      kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                      ! The index to take is 1
                      Stkin(:,kkfreq) = Stk(:,1,1,ith1)

                    end do

                  ! If out of range, take the value at the
                  ! boundary
                  else if (p_frec%omega(jjfreq).ge. &
                           (omega(nfreq) - TINYO)) then

                    ! We are in the last frequency
                    lifreq = nfreq

                    ! Inclinations
                    do ith1=1,Geom%nTh

                      ! Output index
                      kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                      ! The index to take is nfreq
                      Stkin(:,kkfreq) = Stk(:,nfreq,1,ith1)

                    end do

                  ! If within the boundaries
                  else

                    ! Search between the last found frequency and
                    ! all but the boundary
                    do ibfreq=lifreq,nfreq-1

                      ! If this exact frequency is in output
                      if (abs(p_frec%omega(jjfreq)- &
                              omega(ibfreq)).lt.TINYO) then

                        ! We are in the found frequency
                        lifreq = ibfreq

                        ! Inclinations
                        do ith1=1,Geom%nTh

                          ! Output index
                          kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                          ! The index to take is nfreq
                          Stkin(:,kkfreq) = Stk(:,lifreq,1,ith1)

                        end do

                        exit

                      ! If the input is between this output and
                      ! the next
                      else if(p_frec%omega(jjfreq).ge. &
                              omega(ibfreq).and. &
                              p_frec%omega(jjfreq).lt. &
                              omega(ibfreq+1)) then

                        ! We found it in the index of the lower
                        lifreq = ibfreq

                        ! Inverse of the distance
                        ! between the two outputs
                        dx = (p_frec%omega(jjfreq) - &
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
                          Stkin(:,kkfreq) = dx*dy + y0

                        end do

                        exit

                      end if ! Check output frequency

                    end do ! Run output frequencies

                  end if ! Check if out of limits

                end do ! Run input frequencies
!$omp end do
                ! Update frequency index
                jjfreq0 = jjfreq0 + p_mfreq

                ! Update Stokes index
                kkfreq0 = kkfreq0 + p_mfreq*Geom%nTh

              end do ! Output frequencies
            end do ! Output frequency ranges
!$omp end parallel

          end if ! Not dynamic

        ! Non-axially symmetric
        else

          ! If dynamic
          if (dyn) then

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

                ! Skip empty
                if (p_mfreq.lt.1) cycle

                ! For each input direction
                do ith1=1,Geom%nth
                  do iph1=1,Geom%nph2

                    ! If angle-dependent, check backward Rayleigh
                    ! scattering
                    if (nfs.eq.1.and. &
                        Geom%V_CScatt(Geom% &
                                i_scatt(iph1,ith1,jdir)).ge.1d0) &
                      cycle

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
                        Stkin(:,kkfreq) = Stk(:,1,iph1,ith1)

                      ! If out of range, take the value at the
                      ! boundary
                      else if (p_frec%omega(jjfreq)*vfac1.ge. &
                               (omega(nfreq) - TINYO)) then

                        ! We are in the last frequency
                        lifreq = nfreq

                        ! The index to take is nfreq
                        Stkin(:,kkfreq) = Stk(:,nfreq,iph1,ith1)

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
                            Stkin(:,kkfreq) = Stk(:,lifreq,iph1,ith1)

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
                            y0 = Stk(:,lifreq,iph1,ith1)

                            ! The second index is the upper
                            dy = Stk(:,lifreq+1,iph1,ith1) - y0

                            ! Inverse of the distance
                            ! between the two outputs
                            dx = (p_frec%omega(jjfreq)*vfac1 - &
                                  omega(lifreq))/ &
                                 (omega(lifreq+1) - omega(lifreq))

                            ! Interpolate
                            Stkin(:,kkfreq) = dx*dy + y0

                            exit

                          end if ! Check output frequency

                        end do ! Run output frequencies

                      end if ! Check if out of limits

                    end do ! Run input frequencies
!$omp end do
                    ! Update index
                    kkfreq0 = kkfreq0 + p_mfreq

                  end do ! Input azimuth
                end do ! Input polar

                ! Update index
                jjfreq0 = jjfreq0 + p_mfreq

              end do ! Output frequencies
            end do ! Output frequency ranges
!$omp end parallel

          ! If static
          else

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

                ! Skip empty
                if (p_mfreq.lt.1) cycle

                ! Reset the search frequency
                lifreq = Fin%ggf0

                ! For each input frequency
!$omp do
                do jfreq=1,p_mfreq

                  ! Advance indexes
                  jjfreq = jjfreq0 + jfreq

                  ! If out of range, take the value at the
                  ! boundary
                  if (p_frec%omega(jjfreq).le. &
                      omega(1)+TINYO) then

                    ! We are still looking in the first one
                    lifreq = 1

                    ! Block counter
                    nblock = -1

                    ! Directions
                    do ith1=1,Geom%nTh
                      do iph1=1,Geom%nPh

                        ! Skip forward scattering two-term
                        if (nfs.eq.1.and. &
                            Geom%V_CScatt(Geom% &
                                    i_scatt(iph1,ith1,jdir)).ge.1d0) &
                          cycle

                        ! Add block
                        nblock = nblock + 1

                        ! Output index
                        kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                        ! The index to take is 1
                        Stkin(:,kkfreq) = Stk(:,1,iph1,ith1)

                      end do
                    end do

                  ! If out of range, take the value at the
                  ! boundary
                  else if (p_frec%omega(jjfreq).ge. &
                           (omega(nfreq) - TINYO)) then

                    ! We are in the last frequency
                    lifreq = nfreq

                    ! Block counter
                    nblock = -1

                    ! Directions
                    do ith1=1,Geom%nTh
                      do iph1=1,Geom%nPh

                        ! Skip forward scattering two-term
                        if (nfs.eq.1.and. &
                            Geom%V_CScatt(Geom% &
                                    i_scatt(iph1,ith1,jdir)).ge.1d0) &
                          cycle

                        ! Add block
                        nblock = nblock + 1

                        ! Output index
                        kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                        ! The index to take is nfreq
                        Stkin(:,kkfreq) = Stk(:,nfreq,iph1,ith1)

                      end do
                    end do

                  ! If within the boundaries
                  else

                    ! Search between the last found frequency and
                    ! all but the boundary
                    do ibfreq=lifreq,nfreq-1

                      ! If this exact frequency is in output
                      if (abs(p_frec%omega(jjfreq)- &
                              omega(ibfreq)).lt.TINYO) then

                        ! We are in the found frequency
                        lifreq = ibfreq

                        ! Block counter
                        nblock = -1

                        ! Directions
                        do ith1=1,Geom%nTh
                          do iph1=1,Geom%nPh

                            ! Skip forward scattering two-term
                            if (nfs.eq.1.and. &
                                Geom%V_CScatt(Geom% &
                                    i_scatt(iph1,ith1,jdir)).ge.1d0) &
                              cycle

                            ! Add block
                            nblock = nblock + 1

                            ! Output index
                            kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                            ! This frequency gives us the value
                            Stkin(:,kkfreq) = Stk(:,lifreq,iph1,ith1)

                          end do
                        end do

                        exit

                      ! If the input is between this output and
                      ! the next
                      else if(p_frec%omega(jjfreq).ge. &
                              omega(ibfreq).and. &
                              p_frec%omega(jjfreq).lt. &
                              omega(ibfreq+1)) then

                        ! We found it in the index of the lower
                        lifreq = ibfreq

                        ! Inverse of the distance
                        ! between the two outputs
                        dx = (p_frec%omega(jjfreq) - &
                              omega(lifreq))/ &
                             (omega(lifreq+1) - omega(lifreq))

                        ! Block counter
                        nblock = -1

                        ! Directions
                        do ith1=1,Geom%nTh
                          do iph1=1,Geom%nPh

                            ! Skip forward scattering two-term
                            if (nfs.eq.1.and. &
                                Geom%V_CScatt(Geom% &
                                    i_scatt(iph1,ith1,jdir)).ge.1d0) &
                              cycle

                            ! Add block
                            nblock = nblock + 1

                            ! The first index is the lower
                            y0 = Stk(:,lifreq,iph1,ith1)

                            ! The second index is the upper
                            dy = Stk(:,lifreq+1,iph1,ith1) - y0

                            ! Output index
                            kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                            ! Interpolate
                            Stkin(:,kkfreq) = dx*dy + y0

                          end do
                        end do

                        exit

                      end if ! Check output frequency

                    end do ! Run output frequencies

                  end if ! Check if out of limits

                end do ! Run input frequencies
!$omp end do
                ! Update frequency index
                jjfreq0 = jjfreq0 + p_mfreq

                ! Update Stokes index
                kkfreq0 = kkfreq0 + p_mfreq*(Geom%nTh*Geom%nPh-nfs)

              end do ! Output frequencies
            end do ! Output frequency ranges
!$omp end parallel

          end if ! Static
        end if ! Non-axial
      end if ! Interpolation data

      ! Nullify
      nullify(p_mfreq)

      end subroutine getStkin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolates JKQ to the requested frequency\n
      !!  omega(dfloat(:)): Frequency array\n
      !!  JKQ(dcomplex(:)): Radiation field tensor\n
      !!    kfreq(integer): Frequency index of the output frequency
      !!                    associated to the requested input
      !!                    frequency (shifted to limited vector)\n
      !!    nfreq(integer): Size of input vectors\n
      !!                    associated to the requested input
      !!                    frequency (shifted to limited vector)\n
      !!         x(dfloat): Input frequency to interpolate into
      function getJKQinnu(omega,JKQ,kfreq,nfreq,x)

      ! I/O
      integer, intent(in):: kfreq,nfreq
      double precision, intent(in):: x
      double precision, dimension(:), intent(in):: omega
      complex(kind=8), dimension(:), intent(in):: JKQ

      complex(kind=8):: getJKQinnu

      ! Local
      integer:: jfreq,ifreq
      double precision:: dxs
      complex(kind=8):: dys


      ! Check limits
      if (kfreq.lt.1) then
        ifreq = 1
      else if (kfreq.gt.nfreq) then
        ifreq = nfreq
      else
        ifreq = kfreq
      end if

      ! Initialize as "equal"
      getJKQinnu = JKQ(ifreq)

      ! If omegai > omega(ifreq)
      if (x.gt.omega(ifreq)) then

        ! If out of right boundary
        if (x.ge.omega(nfreq)-TINYO) then

          getJKQinnu = JKQ(nfreq)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,nfreq-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getJKQinnu = JKQ(jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq).and. &
                    x.lt.omega(jfreq+1)) then

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
        if (x.le.omega(1)+TINYO) then

          getJKQinnu = JKQ(1)

          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,2,-1

            ! If this exact frequency is in
            ! output
            if (abs(x - omega(jfreq)).lt.TINYO) then

              getJKQinnu = JKQ(jfreq)

              return

            ! If the input is between this
            ! output and the next
            else if(x.ge.omega(jfreq-1).and. &
                    x.lt.omega(jfreq)) then

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

      !> Interpolates JKQ to the input frequency axis\n
      !!  p_frec(Frequencyd_class): Input frequency data for output
      !!                            index and input transition\n
      !!    Fin(Frequencyc2_class): Input frequency data for output
      !!                            index\n
      !!           nmfreq(integer): Size of frequency space\n
      !!          omega(dfloat(:)): Frequency array\n
      !!        Jin(dcmplx(:,:,:)): Interpolated mean intensity\n
      !!        JKQ(dcmplx(:,:,:)): Original mean intensity
      subroutine getJKQin(p_frec,Fin,nmfreq,omega,Jin,JKQ)

      ! I/O
      type(Frequencyd_class), intent(in), pointer:: p_frec
      type(Frequencyc2_class), intent(in):: Fin
      integer, intent(in):: nmfreq
      double precision, dimension(:), intent(in):: omega
      complex(kind=8), dimension(:,:,:), &
                                        allocatable, intent(out):: Jin
      complex(kind=8), dimension(-2:2,0:2,Fin%ggf0:Fin%ggf1), &
                                                      intent(in):: JKQ

      ! Local
      integer:: lifreq,ibfreq
      integer:: iran,ifreq,iifreq,jfreq,jjfreq0,jjfreq

      double precision:: dx

      complex(kind=8), dimension(0:2,0:2):: y0, dy

      integer, pointer:: p_mfreq


      ! Nullify
      nullify(p_mfreq)

      ! If no size, return
      if (nmfreq.lt.1) return

      ! Allocate J
      allocate(Jin(nmfreq,0:2,0:2))

      ! If interpolation data
      if (p_frec%RAM) then

!$omp parallel default(none) &
!$omp private(jjfreq,iifreq,iran,ifreq,p_mfreq,dy,y0) &
!$omp shared(Jin,JKQ,Fin,p_frec)

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
            if (p_mfreq.lt.1) cycle

            ! Input frequencies
!$omp do
            do jfreq=jjfreq+1,jjfreq+p_mfreq

              ! Interpolate
              y0 = JKQ(0:2,0:2,p_frec%index1(jfreq))

              ! Linear interpolation
              dy = JKQ(0:2,0:2,p_frec%index2(jfreq)) - y0

              ! Add J00 interpolated
              Jin(jfreq,0:2,0:2) = dy*p_frec%dx(jfreq) + y0

            end do
!$omp end do

            ! Update jjfreq
            jjfreq = jjfreq + p_mfreq

          end do ! output frequencies
        end do ! output frequencies ranges
!$omp end parallel

      ! No interpolation data
      else

!$omp parallel default(none) &
!$omp private(iifreq,lifreq,iran,ifreq,jfreq,p_mfreq) &
!$omp private(jjfreq,jjfreq0,ibfreq) &
!$omp shared(p_frec,omega,Fin,nfreq,Jin,JKQ)

        ! Initialize index
        jjfreq0 = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Fin%nran
          do ifreq=Fin%if0(iran),Fin%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Input frequency number
            p_mfreq => p_frec%mfreq(iifreq)

            ! Reset the search frequency
            lifreq = Fin%ggf0

            ! For each input frequency
!$omp do
            do jfreq=1,p_mfreq

              ! Advance indexes
              jjfreq = jjfreq0 + jfreq

              ! If out of range, take the value at the
              ! boundary
              if (p_frec%omega(jjfreq).le.omega(1)+TINYO) then

                ! We are still looking in the first one
                lifreq = 1

                ! The index to take is 1
                Jin(jjfreq,0:2,0:2) = JKQ(0:2,0:2,1)

              ! If out of range, take the value at the
              ! boundary
              else if (p_frec%omega(jjfreq).ge. &
                       (omega(nfreq) - TINYO)) then

                ! We are in the last frequency
                lifreq = nfreq

                ! The index to take is nfreq
                Jin(jjfreq,0:2,0:2) = JKQ(0:2,0:2,nfreq)

              ! If within the boundaries
              else

                ! Search between the last found frequency and
                ! all but the boundary
                do ibfreq=lifreq,nfreq-1

                  ! If this exact frequency is in output
                  if (abs(p_frec%omega(jjfreq) - &
                          omega(ibfreq)).lt.TINYO) then

                    ! We are in the found frequency
                    lifreq = ibfreq

                    ! This frequency gives us the value
                    Jin(jjfreq,0:2,0:2) = JKQ(0:2,0:2,lifreq)

                    exit

                  ! If the input is between this output and
                  ! the next
                  else if(p_frec%omega(jjfreq).ge. &
                          omega(ibfreq).and. &
                          p_frec%omega(jjfreq).lt. &
                          omega(ibfreq+1)) then

                    ! We found it in the index of the lower
                    lifreq = ibfreq

                    ! Inverse of the distance
                    ! between the two outputs
                    dx = (p_frec%omega(jjfreq) - omega(lifreq))/ &
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
!$omp end do
            ! Update index in general
            jjfreq0 = jjfreq0 + p_mfreq

          end do ! Output frequencies
        end do ! Output frequency ranges
!$omp end parallel
      end if ! Interpolation data

      ! Free pointers
      nullify(p_mfreq)

      end subroutine getJKQin

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeffaux_mod
