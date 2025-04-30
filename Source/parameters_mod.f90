      !> Physical constants and unit transformations
      module parameters_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/12/2017
!  Last version:
!     30/05/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     30/05/2025:    V4.0.1 - Reduced TINYR to 1d-17 (TdPA)
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
!    General constants and parameters
!
!#####################################################################
!#####################################################################
!#####################################################################

      !
      ! Physical parameters
      !

      !> Lightspeed[10^9 m s^-1]
      double precision, parameter:: c = 2.99792458d-1

      !> Constant part of mean velocity
      double precision, parameter:: dopp = 4.301428d-7

      !> B[G] to Larmor [10^8 s^-1]
      double precision, parameter:: B2L = 8.7940981d-2

      !> B[G] to Larmor [10^5 cm^-1]
      double precision, parameter:: B2Lk = 4.6686437d-10

      !> h*c*10^-4/kb for [10^5 cm^-1]
      double precision, parameter:: c2 = 1.438769d1

      !> [10^5 cm^-1] to [10^8 s^-1]
      double precision, parameter:: k2f = 1.88365156731d8

      !> h*c for [10^5 cm^-1]
      double precision, parameter:: convF = 1.9864474d-11

      !> Rydberg energy [10^5 cm^-1]
      double precision, parameter:: ryd = 1.09737405d0

      !> Rydberg energy H [10^5 cm^-1]
      double precision, parameter:: rydH = 1.09677576d0

      !> Atomic mass unit [kg]
      double precision, parameter:: amu = 1.660539040d-27

      !> Avogadro's number [mol^-1]
      double precision, parameter:: Avog = 6.022140857d23

      !> electron mass [kg]
      double precision, parameter:: me = 9.10938356d-31

      !> electron mass [AMU]
      double precision, parameter:: mem = 5.4857990932872017d-4

      !> Hydrogen mass [AMU]
      double precision, parameter:: mhm = 1.00794d0

      !> Helium mass [AMU]
      double precision, parameter:: mhem = 4.002602d0

      !> Mean atomic weight [AMU]
      double precision, parameter:: armass = 28d0

      !> Boltzmann Constant [J K^-1]
      double precision, parameter:: kb = 1.38064852d-23

      !> Electron charge [C]
      double precision, parameter:: qel = 1.60217662d-19

      !> Vacuum permittivity [F m^-1]
      double precision, parameter:: eps0 = 8.854187817d-12

      !> Planck constant [J s]
      double precision, parameter:: hplanck = 6.62607004d-34

      !> Hydrogen polarizability divided by pi4eps0 [m^3]
      double precision, parameter:: a4pieps0 = 6.67d-31

      !> Bohr radius [m]
      double precision, parameter:: rb = 5.2917721067d-11

      !> T[eV] to ktoev/T[K]
      double precision, parameter:: ktoev = 5.03974756d3

      !> (h^2/2/pi/me/kb)^(3/2)/2 [cm^3]
      double precision, parameter:: cSaha = 2.0706661673662705d-16

      !> pi
      double precision, parameter:: PI = 3.14159265358979323846d0

      !> 2*pi
      double precision, parameter:: PI2 = 6.283185307179586231996d0

      !> 4*pi
      double precision, parameter:: PI4 = 1.256637061435917246399d1

      !> 4*pi*sqrt(PI)
      double precision, parameter:: PI41 = 2.22733119873268279321d1

      !> 4*pi*pi
      double precision, parameter:: PI42 = 3.94784176043574319692d1

      !> 1/PI
      double precision, parameter:: IPI = 3.18309886183790691216d-1

      !> 1/2*pi
      double precision, parameter:: IPI2 = 1.5915494309189534560d-1

      !> 1/4*pi
      double precision, parameter:: IPI4 = 7.9577471545947672804d-2

      !> 1/4*pi*sqrt(pi)
      double precision, parameter:: IPI41 = 4.489678053129164681d-2

      !> 1/4*pi*pi
      double precision, parameter:: IPI42 = 2.533029591058444385d-2

      !> 4*pi*Vacuum permittivity [F m^-1]
      double precision, parameter:: pi4eps0 = 4d0*PI*eps0

      !> Factor for the oscillator strength with gaunt factor
      double precision, parameter:: fnn = 3.9205610341105217d0

      !> Transformation from sigma to FWHM in a gaussian
      double precision, parameter:: sg2fw = 2.35482004503094933d0

      !> Transformation from FWHM to sigma in a gaussian
      double precision, parameter:: fw2sg = 4.2466090014400953d-1

      !> Constant that goes into the oscillator strength collisional
      !! rate equation [cgs]
      double precision, parameter:: ccons = 8.629112601874905d-6


      !
      ! Mathematical shortcut
      !

      !> square root of 2
      double precision, parameter:: sqrt2 = 1.4142135623730950488d0

      !> square root of 3
      double precision, parameter:: sqrt3 = 1.7320508075688772935d0

      !> square root of 5
      double precision, parameter:: sqrt5 = 2.2360679774997898051d0


      !
      ! Code 'hardcoded' parameters
      !

      !> A really small number
      double precision, parameter:: VTINY = 1d-299

      !> What magnetic field is negligible
      double precision, parameter:: TINYB = 1d-100

      !> Difference between angles to be different
      double precision, parameter:: TINYA = 1d-6

      !> Small float (CLE radius)
      double precision, parameter:: TINYF = 1d-6

      !> Single precision
      double precision, parameter:: TINYSP = 1d-8

      !> Double precision
      double precision, parameter:: TINYDP = 1d-16

      !> Difference between frequencies to be different
      double precision, parameter:: TINYO = 1d-15

      !> When an angle is very small
      double precision, parameter:: TINYANG = 1d-20

      !> When a J-symbol is 0
      double precision, parameter:: TINYJS = 1d-20

      !> When the normalization is zero
      double precision, parameter:: TINYN = 1d-15

      !> When the intensity is 0
      double precision, parameter:: TINYI = 1d-100

      !> When the Stokes parameter is 0
      double precision, parameter:: TINYS = 1d-100

      !> When the intensity contribution function is 0
      double precision, parameter:: TINYCI = 1d-100

      !> When the contribution function is 0
      double precision, parameter:: TINYCS = 1d-100

      !> When two optical depths are the same
      double precision, parameter:: TINYT = 1d-20

      !> When two cosines are the same
      double precision, parameter:: TINYM = 1d-15

      !> When a rhoKQ is too small (relative to rho00)
      double precision, parameter:: TINYR = 1d-17

      !> When a J00 is too small to be accounter for in MRC
      double precision, parameter:: TINYMRCJ = 1d-6

      !> When a rhoKQ is too small (relative to rho00) to be
      !! accounter for in MRC
      double precision, parameter:: TINYMRCR = 1d-6

      !> When a rho00 is too small to do MRC
      double precision, parameter:: TINYMRC0 = 1d-15

      !> When a eigenvector is too small
      double precision, parameter:: TINYEV = 1d-20

      !> When a coefficient is too small in emiss2ord
      double precision, parameter:: TINYCO = 1d-20

      !> When the addition of rhoKQ in emiss2ord is too small
      double precision, parameter:: TINYER = 1d-20

      !> When a rho00 in SEE is too small
      double precision, parameter:: TINYR0 = 1d-99

      !> Small fraction of ion
      double precision, parameter:: TINYFRC = 1d-24

      !> Small value for profile
      double precision, parameter:: TINYPRO = 1d-32

      !> Small value for second order redistribution profile
      double precision, parameter:: TINYWAR = 1d-32

      !> Small value for rho00 correction in initpopu
      double precision, parameter:: TINYR00 = 1d-64

      !> Small value in SVD solution
      double precision, parameter:: TINYSVDS = 1d-20

      !> Small penalty for inversion regularization
      double precision, parameter:: TINYPT = 1d-10

      !> Small regularization limit in the inversion
      double precision, parameter:: TINYREG = 1d-5

      !> Small velocity amplitude [Code units]
      double precision, parameter:: TINYVEL = 1d-12/c

      !> Number to indicate bad normalization of profile
      double precision, parameter:: BADNORM = 0.95d0

      !> Fraction of Thermal Doppler width below which fully AA for
      !! redistribution
      double precision, parameter:: vrfrac = 1d-1

      !> Fraction that the limits of a line are extended when reading
      !! input spectra
      double precision, parameter:: isfrac = 1d-1

      !> Distance for an input frequency in input spectro to be close
      !! to an actual frequency for bound-bound transitions
      double precision, parameter:: bbdis = 1d0

      !> Distance for an input frequency in input spectro to be close
      !! to an actual frequency for bound-free transitions
      double precision, parameter:: bfdis = 10d0

      !> Maximum resolution of output frequencies [nm]
      double precision, parameter:: resol = 1d-6

      !> Maximum resolution of intput frequencies [nm]
      double precision, parameter:: resolin = 1d-6

      !> Small number that avoids division by zero when no opacity
      double precision, parameter:: vacuum = 1d-50

      !> Argument for which, in WfuncI, the gaussian is assumed to be
      !! zero
      double precision, parameter:: WbiggaussI = 26d0

      !> Argument for which, in Wfunc, the gaussian is assumed to be
      !! zero
      double precision, parameter:: Wbiggauss = 26d0

      !> When x in exp(-x) is small to use second order series
      double precision, parameter:: smallexp = 1d-7

      !> When x in exp(-x) is large enough to consider exp(-x) = 0
      double precision, parameter:: bigexp = 3d2

      !> When x in exp(x) is close to produce overflow
      double precision, parameter:: vbigexp = 7d2

      !> Value of exp(x) when x > vbigexp
      double precision, parameter:: vbigexpv = 8.22d307

      !> Value of exponential to apply Wien
      double precision, parameter:: wien_limit = 37d0

      !> Distance from line center to consider a Kurucz line in a
      !! certain wavelength
      double precision, parameter:: KDT = 20.0

      !> Allowed fractional difference between energies to be
      !! considered the same in Kurucz lines
      double precision, parameter:: kdif = 0.02d0

      !> Maximum distance in nm between frequency nodes to consider
      !! two different ranges
      double precision, parameter:: jump = 20d0


      !
      ! Convertions
      !

      !> rad to deg
      double precision, parameter:: RAD = 1.8d2/PI

      !> [erg] to [eV]
      double precision, parameter:: ergtoev = 6.241509126d11

      !> [10^5 cm^-1] to [eV]
      double precision, parameter:: fktoev = 1.2398419739d1

      !> [10^5 cm^-1] to [J]
      double precision, parameter:: fktoJ = 1.986445824d-18

      !> Conversion between flu and Aul [10^8 s^-1], with frequency
      !! in 10^5 cm^-1
      double precision, parameter:: CfA = 6.6702516615362047d1


      !
      ! Complex basic numbers
      !

      !> Zero in complex
      complex(kind=8), parameter:: cZero = dcmplx(0d0, 0d0)

      !> 1 in complex
      complex(kind=8), parameter:: cOne  = dcmplx(1d0, 0d0)

      !> i in complex
      complex(kind=8), parameter:: cImag = dcmplx(0d0, 1d0)

      !> small real in complex
      complex(kind=8), parameter:: cTINY = dcmplx(1d-22, 0d0)

      end module parameters_mod
