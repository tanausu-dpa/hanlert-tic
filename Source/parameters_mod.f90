      !> Physical constants and unit transformations
      module parameters_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     12/20/2017
!  Last version:
!     08/07/2023 V3.0.6
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.6 - Added ccons and sqrt5 (TdPA)
!
!     07/03/2023:    V3.0.5 - Added TINIDP (TdPA)
!
!     03/15/2023:    V3.0.4 - Added TINIPT and TINYREG (TdPA)
!
!     03/08/2023:    V3.0.3 - Added sg2fw, fw2sg, and TINYSVDS (TdPA)
!
!     11/24/2022:    V3.0.2 - Added isfrac, bbdis, and bfdis (TdPA)
!
!     10/25/2022:    V3.0.1 - Added TINYF and TINYANG (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added TINYSP (TdPA)
!
!     10/26/2020:   V1.1.12 - Added rydH and fnn and added d0 to the
!                             mhem definition (TdPA)
!
!     05/11/2020:   V1.1.11 - Added Wbiggauss and WbiggaussI (TdPA)
!
!     03/05/2020:   V1.1.10 - Added TINYR00 for rho00 correction in
!                             Initpopu (TdPA)
!
!     09/26/2019:    V1.1.9 - Added Avog, Avogadro's number (TdPA)
!
!     09/13/2019:    V1.1.8 - Added vrfrac, the fraction of the
!                             Doppler width above which emiss2ord for
!                             the angle-average and dynamic case takes
!                             into account Doppler shifts (TdPA)
!
!     07/23/2019:    V1.1.7 - Changed resol and resolin values because
!                             they changed units (TdPA)
!
!     07/19/2019:    V1.1.6 - Added BADNORM and jump (TdPA)
!
!     05/08/2019:    V1.1.5 - Added TINYPRO and TINYWAR (TdPA)
!
!     04/08/2019:    V1.1.4 - Added TINYFRC and kdif (TdPA)
!
!     03/22/2019:    V1.1.3 - Added CfA and KDT (TdPA)
!
!     03/18/2019:    V1.1.2 - Added VTINY, TINYR0, vbigexp, vbigexpv,
!                             and wien_limit (TdPA)
!
!     03/12/2019:    V1.1.1 - Added vbigexp and vbigexpv (TdPA)
!
!     02/20/2019:    V1.1.0 - Changed the TINY numerical quantities
!                             to specific for different uses (TdPA)
!                           - Added a smallexp to use only with
!                             exponential and removed small (TdPA)
!
!     05/16/2018:    V1.0.1 - Added some PI related quantities (TdPA)
!
!     12/20/2017:    V1.0.0 - First version with header (TdPA)
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
!    General constants and parameters
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Physical parameters
      double precision, parameter:: c = 2.99792458d-1               !> Lightspeed[10^9 m s^-1]
      double precision, parameter:: dopp = 4.301428d-7              !> Constant part of mean velocity
      double precision, parameter:: B2L = 8.7940981d-2              !> B[G] to Larmor[10^8 s^-1]
      double precision, parameter:: B2Lk = 4.6686437d-10            !> B[G] to Larmor [10^5 cm^-1]
      double precision, parameter:: c2 = 1.438769d1                 !> h*c*10^-4/kb for [10^5 cm^-1]
      double precision, parameter:: k2f = 1.88365156731d8           !> [10^5 cm^-1] to [10^8 s^-1]
      double precision, parameter:: convF = 1.9864474d-11           !> h*c for [10^5 cm^-1]
      double precision, parameter:: ryd = 1.09737405d0              !> Rydberg energy [10^5 cm^-1]
      double precision, parameter:: rydH = 1.09677576d0             !> Rydberg energy H [10^5 cm^-1]
      double precision, parameter:: amu = 1.660539040d-27           !> Atomic mass unit [kg]
      double precision, parameter:: Avog = 6.022140857d23           !> Avogadro's number [mol^-1]
      double precision, parameter:: me = 9.10938356d-31             !> electron mass [kg]
      double precision, parameter:: mem = 5.4857990932872017d-4     !> electron mass [AMU]
      double precision, parameter:: mhm = 1.00794d0                 !> Hydrogen mass [AMU]
      double precision, parameter:: mhem = 4.002602d0               !> Helium mass [AMU]
      double precision, parameter:: armass = 28d0                   !> Mean atomic weight [AMU]
      double precision, parameter:: kb = 1.38064852d-23             !> Boltzmann Constant [J K^-1]
      double precision, parameter:: qel = 1.60217662d-19            !> Electron charge [C]
      double precision, parameter:: eps0 = 8.854187817d-12          !> Vacuum permittivity [F m^-1]
      double precision, parameter:: hplanck = 6.62607004d-34        !> Planck constant [J s]
      double precision, parameter:: a4pieps0 = 6.67d-31             !> Hydrogen polarizability divided by pi4eps0 [m^3]
      double precision, parameter:: rb = 5.2917721067d-11           !> Bohr radius [m]
      double precision, parameter:: ktoev = 5.03974756d3            !> T[eV] to ktoev/T[K]
      double precision, parameter:: cSaha = 2.0706661673662705d-16  !> (h^2/2/pi/me/kb)^(3/2)/2 [cm^3]
      double precision, parameter:: PI = 3.14159265358979323846d0   !> pi
      double precision, parameter:: PI2 = 6.283185307179586231996d0 !> 2*pi
      double precision, parameter:: PI4 = 1.256637061435917246399d1 !> 4*pi
      double precision, parameter:: PI41 = 2.22733119873268279321d1 !> 4*pi*sqrt(PI)
      double precision, parameter:: PI42 = 3.94784176043574319692d1 !> 4*pi*pi
      double precision, parameter:: IPI = 3.18309886183790691216d-1 !> 1/PI
      double precision, parameter:: IPI2 = 1.5915494309189534560d-1 !> 1/2*pi
      double precision, parameter:: IPI4 = 7.9577471545947672804d-2 !> 1/4*pi
      double precision, parameter:: IPI41 = 4.489678053129164681d-2 !> 1/4*pi*sqrt(pi)
      double precision, parameter:: IPI42 = 2.533029591058444385d-2 !> 1/4*pi*pi
      double precision, parameter:: pi4eps0 = 4d0*PI*eps0           !> 4*pi*Vacuum permittivity [F m^-1]
      double precision, parameter:: fnn = 3.9205610341105217d0      !> Factor for the oscillator strength with gaunt factor
      double precision, parameter:: sg2fw = 2.35482004503094933d0   !> Transformation from sigma to FWHM in a gaussian
      double precision, parameter:: fw2sg = 4.2466090014400953d-1   !> Transformation from FWHM to sigma in a gaussian
      double precision, parameter:: ccons = 8.629112601874905d-6    !> Constant that goes into the oscillator strength collisional rate equation [cgs]

      ! Mathematical shortcut
      double precision, parameter:: sqrt2 = 1.4142135623730950488d0 !> square root of 2
      double precision, parameter:: sqrt3 = 1.7320508075688772935d0 !> square root of 3
      double precision, parameter:: sqrt5 = 2.2360679774997898051d0 !> square root of 5

      ! Code 'hardcoded' parameters
      double precision, parameter:: VTINY = 1d-299                  !> A really small number
      double precision, parameter:: TINYB = 1d-100                  !> What magnetic field is negligible
      double precision, parameter:: TINYA = 1d-6                    !> Difference between angles to be different
      double precision, parameter:: TINYF = 1d-6                    !> Small float (CLE radius)
      double precision, parameter:: TINYSP = 1d-8                   !> Single precision
      double precision, parameter:: TINYDP = 1d-16                  !> Double precision
      double precision, parameter:: TINYO = 1d-15                   !> Difference between frequencies to be different
      double precision, parameter:: TINYANG = 1d-20                 !> When an angle is very small
      double precision, parameter:: TINYJS = 1d-20                  !> When a J-symbol is 0
      double precision, parameter:: TINYN = 1d-15                   !> When the normalization is zero
      double precision, parameter:: TINYI = 1d-100                  !> When the intensity is 0
      double precision, parameter:: TINYS = 1d-100                  !> When the Stokes is 0
      double precision, parameter:: TINYCI = 1d-100                 !> When the intensity contribution function is 0
      double precision, parameter:: TINYCS = 1d-100                 !> When the contribution function is 0
      double precision, parameter:: TINYT = 1d-20                   !> When two optical depths are the same
      double precision, parameter:: TINYM = 1d-15                   !> When two mues are the same
      double precision, parameter:: TINYR = 1d-15                   !> When a rhoKQ is too small (relative to rho00)
      double precision, parameter:: TINYMRCJ = 1d-6                 !> When a J00 is too small to be accounter for in MRC
      double precision, parameter:: TINYMRCR = 1d-6                 !> When a rhoKQ is too small (relative to rho00) to be accounter for in MRC
      double precision, parameter:: TINYMRC0 = 1d-15                !> When a rho00 is too small to do MRC
      double precision, parameter:: TINYEV = 1d-20                  !> When a eigenvector is too small
      double precision, parameter:: TINYCO = 1d-20                  !> When a coefficient is too small in emiss2ord
      double precision, parameter:: TINYER = 1d-20                  !> When the addition of rhoKQ in emiss2ord is too small
      double precision, parameter:: TINYR0 = 1d-99                  !> When a rho00 in SEE is too small
      double precision, parameter:: TINYFRC = 1d-24                 !> Small fraction of ion
      double precision, parameter:: TINYPRO = 1d-32                 !> Small value for profile
      double precision, parameter:: TINYWAR = 1d-32                 !> Small value for second order redistribution profile
      double precision, parameter:: TINYR00 = 1d-64                 !> Small value for rho00 correction in initpopu
      double precision, parameter:: TINYSVDS = 1d-20                !> Small value in SVD solution
      double precision, parameter:: TINYPT = 1d-10                  !> Small penalty for inversion regularization
      double precision, parameter:: TINYREG = 1d-5                  !> Small regularization limit in the inversion
      double precision, parameter:: BADNORM = 0.95d0                !> Number to indicate bad normalization of profile
      double precision, parameter:: vrfrac = 1d-1                   !> Fraction of Thermal Doppler width below which fully AA for redistribution
      double precision, parameter:: isfrac = 1d-1                   !> Fraction that the limits of a line are extended when reading input spectra
      double precision, parameter:: bbdis = 1d0                     !> Distance for an input frequency in input spectro to be close to an actual frequency for bound-bound transitions
      double precision, parameter:: bfdis = 10d0                    !> Distance for an input frequency in input spectro to be close to an actual frequency for bound-free transitions
      double precision, parameter:: resol = 1d-6                    !> Maximum resolution of output frequencies [nm]
      double precision, parameter:: resolin = 1d-6                  !> Maximum resolution of intput frequencies [nm]
      double precision, parameter:: vacuum = 1d-50                  !> Small number that avoids division by zero when no opacity
      double precision, parameter:: WbiggaussI = 26d0               !> Argument for which, in WfuncI, the gaussian is assumed to be zero
      double precision, parameter:: Wbiggauss = 26d0                !> Argument for which, in Wfunc, the gaussian is assumed to be zero
      double precision, parameter:: smallexp = 1d-7                 !> When x in exp(-x) is small to use second order series
      double precision, parameter:: bigexp = 3d2                    !> When x in exp(-x) is large enough to consider exp(-x) = 0
      double precision, parameter:: vbigexp = 7d2                   !> When x in exp(x) is close to produce overflow
      double precision, parameter:: vbigexpv = 8.22d307             !> Value of exp(x) when x > vbigexp
      double precision, parameter:: wien_limit = 37d0               !> Value of exponential to apply Wien
      double precision, parameter:: KDT = 20.0                      !> Distance from line center to consider a Kurucz line in a certain wavelength
      double precision, parameter:: kdif = 0.02d0                   !> Allowed fractional difference between energies to be considered the same in Kurucz lines
      double precision, parameter:: jump = 20d0                     !> Maximum distance in nm between frequency nodes to consider two different ranges

      ! Convertions
      double precision, parameter:: RAD = 1.8d2/PI                  !> rad to deg
      double precision, parameter:: ergtoev = 6.241509126d11        !> [erg] to [eV]
      double precision, parameter:: fktoev = 1.2398419739d1         !> [10^5 cm^-1] to [eV]
      double precision, parameter:: fktoJ = 1.986445824d-18         !> [10^5 cm^-1] to [J]
      double precision, parameter:: CfA = 6.6702516615362047d1      !> Conversion between flu and Aul [10^8 s^-1], with frequency in 10^5 cm^-1

      ! Complex basic numbers
      complex(kind=8), parameter:: cZero = dcmplx(0d0, 0d0)         !> Zero in complex
      complex(kind=8), parameter:: cOne  = dcmplx(1d0, 0d0)         !> 1 in complex
      complex(kind=8), parameter:: cImag = dcmplx(0d0, 1d0)         !> i in complex
      complex(kind=8), parameter:: cTINY = dcmplx(1d-22, 0d0)       !> small real in complex

      end module parameters_mod
