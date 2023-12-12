RUN_MODE
--------

  * MANDATORY

  * Formats:
    
    - string: 1DS, 1.5DS, CLE, inversion; default: 1DS

  * Description: Indicate the type of problem to solve, a single 1D solution (1DS), a 1.5D solution of a 3D model (1.5DS), a 1.5D formal solution in an optically thin plasma (CLE), or pixel by pixel inversion (inversion).

ATMO_INPUT
----------

  * MANDATORY

  * Formats:
    
    - string: file path, HCC (only INV), HCP (only INV); default: HCC (only INV)

  * Description: Path, absolute or relative to the running directory, of the model atmosphere file. In inversion mode this specifies the initial model atmosphere and HCC or HCP can be specified to indicate the C or P models of Fontenla et al. (1993).

ATMO_SCALE
----------

  * OPTIONAL

  * Formats:
    
    - string: Height, tau; default: Height; optional: float (default 500)

  * Description: Type of stratification scale of the input model atmosphere (if 15D) or the initial 3D model atmosphere (if INV). Optional: reference wavelength for the optical depth calculation, in nanometers.

ATMO_CHAR
---------

  * OPTIONAL

  * Formats:
    
    - string: NENH, NE, PE, RHOE, PG, RHO; default: NENH

  * Description: Variable which variable the model atmosphere is providing to characterize the stratification, together with the temperature, both electron and Hydrogen number densities (NENH), electron number densities (NE), electron pressure in dyn cm^-2 (PE), electron density in g cm^-3 (RHOE), gas pressure in dyn cm^-2 (PG), or total density in g cm^-3 (RHO).

ATOM_INPUT
----------

  * MANDATORY, ADDITIVE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the atomic model whose density matrix elements are to be calculated, in principle, in non-LTE. Optional: Unique label to identify this atom in other input options.

ATOM_BACK
---------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the atomic model to include in the calculation as a background opacity source. Their populations will be calculated in LTE if no population file is provided (see ATOM_POPU).

ATOM_POPU
---------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the initial population (if the label points to a model in ATOM_FILE) or the fixed populations (if the label points to a model in ATOM_BACK) for the model atom specified by the label.

ATOM_FIX_POP
------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will kept their initial populations fixed when solving the iterative problem. These populations are calculated in LTE unless ATOM_POPU is specified for the same model atom.

ATOM_FIX_POP_LTERM
------------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will kept the initial populations of the ground level/term fixed when solving the iterative problem. These populations are calculated in LTE unless ATOM_POPU is specified for the same model atom.

ATOM_ZERO_ION
-------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will ignore the transitions with the last level/term in the model atom (of a different ionization stage than the rest of levels/terms), effectively keeping its initial ionization balance.

ATOM_NO_WAVE
------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will not contribute with frequency nodes when building the output frequency/wavelength axis.

MOLECULE_INPUT
--------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the molecular data that will be used to compute the chemical equilibrium and the contribution of some molecules to the background opacity.

OPACITY_FUDGE
-------------

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the fudge factors for the background continuum opacity.

PARFUN
------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with partition function and ionization potential tabulations.

ABUND
-----

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with atomic abundances.

BARK_SP
-------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with the tabulation of the Barklem broadening parameters for s-p transitions.

BARK_PD
-------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with the tabulation of the Barklem broadening parameters for p-d transitions.

BARK_DF
-------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with the tabulation of the Barklem broadening parameters for d-f transitions.

IGNORE_BB
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No

  * Description: Neglect (Yes) or include (No) the bound-bround transitions of the ATOM_BACK atomic models in the formal solution.

KURUCZ
------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with a list of transitions in the Kurucz's database format to be included in the background opacity.

LTE_LINE
--------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - entry
    - string: file path

  * Description: Entry with the data of an atomic transition to be included under the assumption of LTE. The format can be the same from the Kurucz's database (fill with zeros to fulfill the restricted size requirements) or the specific format of the code. | Path, absolute or relative to the running directory, to the file with a list of atomic transitions to be included under the assumption of LTE. The format of each entry in the file can be the same from the Kurucz's database (fill with zeros to fulfill the restricted size requirements) or the specific format of the code.

WAVELENGTHS
-----------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to a file with a list of wavelengths to include in the frequency/wavelength axis in the formal solution.

ASYMM_INPUT
-----------

  * OPTIONAL, ADDITIVE, ADVANCED

  * Formats:
    
    - integer*2 float*2 (only 1D)
    - string: file path

  * Description: Multipolas component idenfifiers K and Q, and real and imaginary values of the JKQ tensor to be considered as an ad-hoc (homogeneous) asymmetry in the formal solution. | Path, absolute or relative to the running directory, to a file with JKQ tensors to be considered as an ad-hoc assymetry in the formal solution.

FORCE_SYMM
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the input asymmetry must be added to the one properly generated in the formal solution (No) or must completely overwrite the actual radiation field tensors (Yes).

STIM
----

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: If stimulated emission must be accounted for.

OUT_FOLDER
----------

  * OPTIONAL

  * Formats:
    
    - string: folder path; default: Outputs/Default

  * Description: Path, absolute or relative to the running directory, to the folder to store the outputs. This folder does not have to exist, but must be possible to create it (the path up to the second to last folder must exist).

POLAR_NODES
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 4

  * Description: Number of nodes per hemisphere in the polar (gaussian) angular quadrature for the solution of the radiative transfer problem.

AXIAL_NODES
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 2

  * Description: Number of nodes per octant in the azimuthal (trapezoidal) angular quadrature for the solution of the radiative transfer problem. A negative value indicates axial symmetry.

POLARI_NODES
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: POLAR_NODES

  * Description: Number of nodes per hemisphere in the polar (gaussian) angular quadrature for the solution of the only intensity radiative transfer problem.

AXIALI_NODES
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: AXIAL_NODES

  * Description: Number of nodes per octant in the azimuthal (trapezoidal) angular quadrature for the solution of the only intensity radiative transfer problem. A negative value indicates axial symmetry.

POLAR_LOS
---------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: List of floats with the cosine of the heliocentric angle for which to compute the Stokes parameters formal solution.

AXIAL_LOS
---------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: List of floats with the azimuthal angle, in degrees, for which to compute the Stokes parameters formal solution.

MODE
----

  * OPTIONAL

  * Formats:
    
    - string: solve, read, both

  * Description: Specify what to do in the formal solution, solving from standard initialization (solve), read a previous solution and continue iterating (both), or read a previous solution and just perform the last formal solutions (read). WARNING: When using both with an only intensity solution file, it will be assumed that such a solution is fully converged.

FORCE
-----

  * OPTIONAL

  * Formats:
    
    - string: intensity, polarization, all, none; default: none

  * Description: Force the code to solve only the intensity problem (intensity), to solve the polarization problem skipping the intensity problem (polarization), to solve both intensity and polarization problem (all), or let the code decide from the available inputs (none).

RESTRICT_TAUC
-------------

  * OPTIONAL

  * Formats:
    
    - float*2

  * Description: Lower and upper limit in the decimal logarithm of the optical depth where to solve the radiative transfer problem. Nodes with values of the optical depths outside the specified range will be neglected.

RESTRICT_HEIGHT
---------------

  * OPTIONAL

  * Formats:
    
    - float*2

  * Description: Lower and upper limit in height in kilometers where to solve the radiative transfer problem. Nodes with values of the height outside the specified range will be neglected.

ZEEMAN_MODE
-----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: full, nozeeman, nosplit, linear; default: full

  * Description: Type of diagonalization of the atomic Hamiltonian, full diagonalization (full), neglect any magnetic splitting (nozeeman), neglect the splitting but not the level coupling (nosplit), or assume linear regime (linear).

P_CORR
------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Add one fake iteration between the solution of the only intensity and full Stokes radiative transfer problems to adjust the radiation field tensors of multi-term atomic models.

FCOL_TRANSFER
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Allow forbidden collisions to transfer polarization (Yes).

MIT_OFF
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No, none; default: none

  * Description: Neglect (Yes) or force the inclusion (No) of frequencies corresponding to magnetically induced transitions in multi-term model atoms. By default, they are only included if there is a magnetic field.

MIT_NODE
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float

  * Description: Factor for which to divide the number of nodes in the frequency/wavelength axis to allocate for the MIT with respect to the permitted transitions of a multiplet.

SKIP_SFILE
----------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the writing of the solution file (which allows initializing other runs) should be skipped (Yes) or not (No).

SOLUTION_INPUT
--------------

  * OPTIONAL

  * Formats:
    
    - string: file path (only 1D)
    - string: folder path (only 15D)

  * Description: Path, absolute or relative to the running directory, to the Solution file of a previous run to be read to initialize this run for appropriate modes. | Path, absolute or relative to the running directory, to the Solution folder of a previous run to be read to initialize this run for appropriate modes.

SOLUTION_KEEPI
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Keep a solution file with the solution of the intensity problem even if solving the polarized problem afterwards.

SOLUTION_KEEPS
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force the solution file to store the full Stokes parameters instead of the radiation field tensors independently of the parameters which decide it automatically.

K_CUT
-----

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum density matrix multipolar component to consider in the polarized problem. Negative means no restriction.

K_CUTAB
-------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the K_CUT must be accounted for in the combination of multipolar components which appears when partial frequency redistribution is taken into account.

K_CUT_TERM
----------

  * OPTIONAL, ADDITIVE, ADVANCED

  * Formats:
    
    - 3*integers; 4*integers

  * Description: If three integers, label for the (active) model atom, term index, and maximum rank of the density matrix allowed for that term. If four integers, label for the (active) model atom, initial and final indexes for a range of terms, and maximum rank of the density matrix allowed for those terms.

K_RAD
-----

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum radiation field tensor multipolar component to consider in the polarized problem. Negative means no restriction.

MEMOJ
-----

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Specify if memoization is to be used to avoid calculating any J-symbol more than once.

PIRAM
-----

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM some quantities used in photoionization to avoid recalculating them over and over.

VOI_TYPE
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: H82, A67, HAW78, GAUSS; default: H82

  * Description: Algorithm to compute the Voigt profiles. If solving for polarization, only H82 is allowed.

VOI_IRAM
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the Voigt profiles in the only intensity problem.

VOI_IFIL
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in a file the Voigt profiles in the only intensity problem.

LTE_VOI_IRAM
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the complex Voigt profiles for LTE transitions.

VOI_PRAM
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in RAM the complex Voigt profiles.

VOI_PFIL
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in a file the complex Voigt profiles.

RAM_LIMIT
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum amount of Megabytes which can be allocated in the form of Voigt profiles, photoionization pre-calculated quantities, interpolation precalculated quantities (partial frequency redistribution), or redistribution functions. Negative means no limit (NOT RECOMMENDED for complex problems).

RAMAN
-----

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Take into account Raman scattering (sometimes referred as cross redistribution, XRD).

NO_COH_L_TERM
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Asumme a non-coherent lower term when computing the radiative transfer coefficients.

RED_COHW
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float / string: No; Default: No

  * Description: Doppler widths from the line center from where to assume that the scattering is fully coherent.

REDI_COHW
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float / string: No; Default: RED_COHW

  * Description: Doppler widths from the line center from where to assume that the scattering is fully coherent in only intensity problems.

RED_MOD
-------

  * OPTIONAL

  * Formats:
    
    - string: AA, AD; default: AA

  * Description: Type of redistribution function, angle-average (AA) or angle-dependent (AD)

RED_AAINT
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force angle-average redistribution function for the only intensity problem, regardless of the RED_MOD input.

RED_IRAM
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the redistribution function of the only intensity problem.

RED_PRAM
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in RAM the complex redistribution function.

INT_IRAM
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the interpolation data for the partial frequency redistribution in the only intensity problem.

INT_PRAM
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the interpolation data for the partial frequency redistribution in the polarized problem.

RED_NODE
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 4

  * Description: Number of nodes per hemisphere in the scattering (gaussian) angular quadrature to compute the angle-averaged redistribution function.

REDI_NODE
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: RED_NODE

  * Description: Number of nodes per hemisphere in the scattering (gaussian) angular quadrature to compute the angle-averaged redistribution function in the only intensity problem.

RED_RANG
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3.5

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral.

RED_RESO
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3.5

  * Description: Distance in Doppler width from a transition resonance to move the search of input frequencies to include the full transition resonance.

RED_NEGL
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e4

  * Description: Distance in Doppler width from a frequency to the closest transition resonance in order to neglect partial frequency redistribution at that frequency.

RED_VLAR
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 7

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral.

RED_FSTP
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0.5

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances.

RED_MSTP
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 4

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances.

RED_CORE
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Doppler width distance to the closest transition resonance to consider a frequency as in a line core.

RED_RANG_CORE
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_RANG

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral, for frequencies in a line core.

RED_VLAR_CORE
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_VLAR

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral, for frequencies in a line core.

RED_FSTP_CORE
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_FSTP

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances, for frequencies in a line core.

RED_MSTP_CORE
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_MSTP

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances, for frequencies in a line core.

REDI_RANG
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_RANG

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral in the only intensity problem.

REDI_RESO
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_RESO

  * Description: Distance in Doppler width from a transition resonance to move the search of input frequencies to include the full transition resonance in the only intensity problem

REDI_NEGL
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_NEGL

  * Description: Distance in Doppler width from a frequency to the closest transition resonance in order to neglect partial frequency redistribution at that frequency in the only intensity problem.

REDI_VLAR
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_VLAR

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral in the only intensity problem.

REDI_FSTP
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_FSTP

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances in the only intensity problem.

REDI_MSTP
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_MSTP

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances in the only intensity problem.

REDI_CORE
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: RED_CORE

  * Description: Doppler width distance to the closest transition resonance to consider a frequency as in a line core in the only intensity problem.

REDI_RANG_CORE
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: REDI_RANG

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral, for frequencies in a line core in the only intensity problem.

REDI_VLAR_CORE
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: REDI_VLAR

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral, for frequencies in a line core in the only intensity problem.

REDI_FSTP_CORE
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: REDI_FSTP

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances, for frequencies in a line core in the only intensity problem.

REDI_MSTP_CORE
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: REDI_MSTP

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances, for frequencies in a line core in the only intensity problem.

DOP_WIDTH
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 2.5e3
    - string: Max, Min

  * Description: Doppler width in m s^-1 to consider when building the frequency axis to transform between Doppler widths and actual frequencies. | Take the maximum (Max) or minimum (Min) temperature to calculate the Doppler width used to build the frequency axis to transform between Doppler widths and actual frequencies.

MIN_T
-----

  * OPTIONAL

  * Formats:
    
    - float; default: 3e3

  * Description: Minimum temperature in kelvin expected in the 3D model atmosphere.

MAX_V
-----

  * OPTIONAL

  * Formats:
    
    - float; default: 1e1

  * Description: Minimum velocity in km s^-1 expected in the 3D model atmosphere.

RT_GROUP_N
----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 1

  * Description: Number of processes to solve the forward or inversion problem for each pixel of the input.

UNMAGNETIZED
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Assume no magnetic field in the input model atmosphere.

STATIC
------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Assume no velocity in the input model atmosphere.

ITER_MIN
--------

  * OPTIONAL

  * Formats:
    
    - integer; default: 1

  * Description: Index of the first iteration of a formal solution.

ITERI_MIN
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: ITER_MIN

  * Description: Index of the first iteration of an only intensity formal solution.

ITER_MAX
--------

  * OPTIONAL

  * Formats:
    
    - integer; default: 500

  * Description: Maximum iteration index allowed for in a formal solution.

ITERI_MAX
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: ITER_MAX

  * Description: Maximum iteration index allowed for in an only intensity formal solution.

ITER_2ORD
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Account for partial frequency redistribution effects.

ITER_NB
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If solving polarization with magnetic field, solve first the non-magnetic problem. It is most beneficial when the intensity problem is axial (and has been configured as axial). Otherwise, it is not necessarily faster, but problem dependent.

ITER_MRC_I
----------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-5

  * Description: Maximum relative change of the populations to consider that they have converged.

ITERI_MRC_I
-----------

  * OPTIONAL

  * Formats:
    
    - float; default: ITER_MRC_I

  * Description: Maximum relative change of the populations to consider that they have converged in the only intensity solution.

ITER_MRC_P
----------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-3

  * Description: Maximum relative change of the non-population density matrix elements to consider that they have converged.

ITER_J
------

  * OPTIONAL

  * Formats:
    
    - integer; default: 5

  * Description: Number of preliminar no-line iterations to perform to relax the initial radiation field.

ITER_PRD
--------

  * OPTIONAL

  * Formats:
    
    - integer; default: 4

  * Description: Maximum number of only-radiation iterations to perform when there is partial frequency redistribution in the only intensity problem.

ITER_MRC_R
----------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-3

  * Description: Maximum relative change of the intensity or mean intensity in only-radiation iterations to consider that it is converged.

NG_ACC
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Accelerate the convergence with Ng's algorithm.

NG_ORDER
--------

  * OPTIONAL

  * Formats:
    
    - integer; default: 3

  * Description: Order of the Ng acceleration.

NG_DELAY
--------

  * OPTIONAL

  * Formats:
    
    - integer; default: 20

  * Description: Iteration index from which data for Ng acceleration starts to accumulate.

NGI_ACC
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: NG_ACC

  * Description: Accelerate the convergence with Ng's algorithm in the only intensity problem.

NGI_ORDER
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: NG_ORDER

  * Description: Order of the Ng acceleration in the only intensity problem.

NGI_DELAY
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: NG_DELAY

  * Description: Iteration index from which data for Ng acceleration starts to accumulate in the only intensity problem.

PRD_DELAY
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 1

  * Description: Iteration index from which to start accounting for partial frequency redistribution in the only intensity problem.

ALI_DELAY
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 1

  * Description: Iteration index from which to start using the accelerated lambda iteration algorithm in the only intensity problem.

BCAST_MODE
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: BROADCAST, ALTSEND; default: BROADCAST

  * Description: Type of algorithm for message passing interface broadcasting.

ALLOW_NPHYS_STK
---------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: -1

  * Description: Iteration index from which finding a non-physical Stokes parameter results in termination of the formal solution. Negative means from the beginning.

ALLOW_NPHYS_RHO
---------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: -1

  * Description: Iteration index from which finding a non-physical density matrix results in termination of the formal solution. Negative means from the beginning.

SOLUTION_BOX
------------

  * OPTIONAL

  * Formats:
    
    - integer*4; default: -1 -1 -1 -1

  * Description: Indicate the initial x index, final x index, initial y index, and final y index, respectively, of the pixels to solve in a 3D model or data file. Negative numbers are wildcards (automatically adjusted to the relevant size of the input).

CONTRIBUTION
------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Calculate and save the contribution function for the last formal solutions.

TAU1
----

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Calculate and save the height where the optical depth is equal to one at each frequency for the last formal solutions.

KEEP_BACK
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the calculated background opacity, scattering coefficient, and emissivity.

KEEP_DAMP
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the damping parameter characteristic of the Voigt profile for the ATOM_FILE model atoms.

KEEP_COLS
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the rates of inelastic collisions between every term and level in the ATOM_FILE model atoms.

KEEP_ATMO
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Create a file with all the atmospheric quantities (original or derived) of the atmospheric model.

KEEP_POP
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the populations of the ATOM_FILE model atoms.

KEEP_DEP
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the departure coefficients of the ATOM_FILE model atoms.

KEEP_RHOKQ
----------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the density matrix elements of the ATOM_FILE model atoms.

KEEP_JKQ
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the line integrated radiation field tensors of the ATOM_FILE model atoms.

KEEP_STOKES_QUAD
----------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the Stokes parameters in the model atomsphere for the quadrature directions.

KEEP_JKQNU
----------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the frequency dependent radiation field tensors.

KEEP_MRC
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save the maximum relative change in a file.

KEEP_MPI_LOG
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Keep a log with the distribution of tasks among CPU.

KEEP_MPI_DETAIL_LOG
-------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Keep a more detailed log (written in binary) with the distribution of tasks and weights among CPU.

LIM_STK
-------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for Stokes parameters. Frequencies not included in the specified ranges will not be included in the output.

LIM_CTR
-------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for the contribution function. Frequencies not included in the specified ranges will not be included in the output.

LIM_TAU
-------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for the height where the optical depth is equal to one. Frequencies not included in the specified ranges will not be included in the output.

LIM_COLS_TT
-----------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string + integer*2

  * Description: Unique identifier of the model atom and doublet of indexes specifying the two terms whose collisional rates between them will be included in the output for collisional rates. Transitions not included in the specified ranges will not be included in the output.

LIM_COLS_LL
-----------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string + integer*2

  * Description: Unique identifier of the model atom and doublet of indexes specifying the two levels whose collisional rates between them will be included in the output for collisional rates. Transitions not included in the specified ranges will not be included in the output.

LIM_DAMP
--------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string + integer

  * Description: Unique identifier of the model atom and index specifying the bound-bound transition whose damping parameter will be included in the output for damping parameters. Transitions not included in the specified ranges will not be included in the output.

LIM_BACK
--------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for the background opacity quantities. Frequencies not included in the specified ranges will not be included in the output.

LIM_POP
-------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string + integer

  * Description: Unique identifier of the model atom and index of the level whose population will be included in the output for populations and departure coefficients. Levels not included in the specified ranges will not be included in the output.

PROTECT_H
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic Hydrogen number density.

CHEM_PROTECT_ALL
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic number density of any atom.

MPIDETAIL
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of a MPI detail file from a previous run.

OPERFORM
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of a performance detail file from a previous run.

MPISTKIP
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No

  * Description: Skip the first iteration when trying to optimize the distribution of frequencies between processes.

VERBOSE
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Indicate if the verbosity output of the code must go to the command line (Yes) or to a verbosity file to be created in the output directory (No).

