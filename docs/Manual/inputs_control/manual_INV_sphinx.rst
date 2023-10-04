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

DATA_FILE
---------

  * MANDATORY

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the data to invert.

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
    
    - float: No

  * Description: Doppler widths from the line center from where to assume that the scattering is fully coherent.

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

STATIC_INT
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Assume no velocity in the input model atmosphere when solving the only intensity problem.

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

VERBOSE_INV_LV
--------------

  * OPTIONAL

  * Formats:
    
    - integer: 0, 1, 2, 3; default: 0

  * Description: Level of verbosity, decides what messages go into which file.

VERBOSE_INV_SHUTUP
------------------

  * OPTIONAL

  * Formats:
    
    - integer: 0, 1, 2, 3; default: 3

  * Description: Largest level of verbosity in the inversion that will be outputted at all.

TYPE_INVERSION
--------------

  * OPTIONAL

  * Formats:
    
    - string: Thermal, Magnetic, All, Sequential, Sequential-Magnetic; default: Thermal

  * Description: Type of inversion, i.e., without magnetic field (thermal), only the magnetic field (magnetic), both (all), both but with a first convergence without magnetic field (sequential), first everythin but magnetic field and then only the magnetic field (sequential-magnetic). Note: only the first word will be registered when reading the input, so do not introduce spaces.

AUTO_WEIGHT
-----------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Let the code choose the weights automatically (NOT RECOMMENDED).

WEIGHT
------

  * MANDATORY, ADDITIVE

  * Formats:
    
    - float*4
    - float*6

  * Description: Weights for Stokes I, Q, U, and V in the merit function (only one entry). Ignored if AUTO_WEIGHT=Yes. Not mandatory if WEIGHT_FILE is specified. | Wavelength range in nm for which these weights must be used, and weights for Stokes I, Q, U, and V in the merit function. Ignored if AUTO_WEIGHT=Yes. Not mandatory if WEIGHT_FILE is specified. Note that the wavelength range does not truncate, but look for the closest wavelength in the data. Moreover, they are intended to be for distinguishing between different spectral ranges; in order to customize weights within a spectral line, use a file to specify the weights.

WEIGHT_FILE
-----------

  * MANDATORY, ADDITIVE

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of a file with the wavelengths ranges and weights for Stokes I, Q, U, and V to be used in the inversion. Ignored if AUTO_WEIGHT=Yes. Not mandatory if WEIGHT is specified.

INV_INIT
--------

  * MANDATORY

  * Formats:
    
    - string: Init, file path

  * Description: Path, absolute or relative to the running directory, of a file with the result of a previous inversion. Init means starting from scratch.

CENTERED_DERIVATIVE
-------------------

  * MANDATORY

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the derivative to compute the numerical response function must be centered, i.e., computing the perturbation with both signs.

ITER_MAX_INV
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 15

  * Description: Maximum iteration index allowed in the inversion.

NODES_B_METHOD
--------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the magnetic field strength or longitudinal component (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_BT_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the magnetic field inclination or transversal component (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_BP_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the magnetic field azimuth (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_T_METHOD
--------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the temperature field azimuth (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VX_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the velocity x component or inclination (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VY_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the velocity y component or azimuth (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VZ_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the velocity vertical or line of sight component (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VT_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the micro-turbulent velocity (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_PG_METHOD
---------------

  * OPTIONAL

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the gas pressure (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J21R_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the real part of the J21 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J21I_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the imaginary part of the J21 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J22R_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the real part of the J22 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J22I_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the imaginary part of the J22 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

INTERPOLATION
-------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: linear, quadratic bezier, cubir bezier; default: cubic bezier

  * Description: Type of interpolation to generate the atmospheric model stratification from the node values.

NODES_B_NUM
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the magnetic field strength or longitudinal component. Mutually exclusive with NODES_B_LOCATION.

NODES_BT_NUM
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the magnetic field inclination or transversal component. Mutually exclusive with NODES_BT_LOCATION.

NODES_BP_NUM
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the magnetic field azimuth. Mutually exclusive with NODES_BP_LOCATION.

NODES_F_NUM
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the diffuse light factor. Mutually exclusive with NODES_F_LOCATION.

NODES_T_NUM
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the temperature. Mutually exclusive with NODES_T_LOCATION.

NODES_VX_NUM
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the velocity x component or transversal component. Mutually exclusive with NODES_VX_LOCATION.

NODES_VY_NUM
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the velocity y component or azimuth. Mutually exclusive with NODES_VY_LOCATION.

NODES_VZ_NUM
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the velocity vertical or longitudinal component. Mutually exclusive with NODES_VZ_LOCATION.

NODES_PG_NUM
------------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the gas pressure. Mutually exclusive with NODES_PG_LOCATION.

NODES_J21R_NUM
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the real part of the J21 radiation field tensor. Mutually exclusive with NODES_J21R_LOCATION.

NODES_J21I_NUM
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the imaginary part of the J21 radiation field tensor. Mutually exclusive with NODES_J21I_LOCATION.

NODES_J22R_NUM
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the real part of the J22 radiation field tensor. Mutually exclusive with NODES_J22R_LOCATION.

NODES_J22I_NUM
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the imaginary part of the J22 radiation field tensor. Mutually exclusive with NODES_J22I_LOCATION.

NODES_B_LOCATION
----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the magnetic field strength or longitudinal component. Mutually exclusive with NODES_B_NUM.

NODES_BT_LOCATION
-----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the magnetic field inclination or transversal component. Mutually exclusive with NODES_BT_NUM.

NODES_BP_LOCATION
-----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the magnetic field azimuth. Mutually exclusive with NODES_BP_NUM.

NODES_F_LOCATION
----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the diffuse light factor. Mutually exclusive with NODES_F_NUM.

NODES_T_LOCATION
----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the temperature. Mutually exclusive with NODES_T_NUM.

NODES_VX_LOCATION
-----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the velocity x component or transversal component. Mutually exclusive with NODES_VX_NUM.

NODES_VY_LOCATION
-----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the velocity y component or azimuth. Mutually exclusive with NODES_VY_NUM.

NODES_VZ_LOCATION
-----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the velocity vertical or longitudinal component. Mutually exclusive with NODES_VZ_NUM.

NODES_PG_LOCATION
-----------------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the gas pressure. Mutually exclusive with NODES_PG_NUM.

NODES_J21R_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the real part of the J21 radiation field tensor. Mutually exclusive with NODES_J21R_NUM.

NODES_J21I_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the imaginary part of the J21 radiation field tensor. Mutually exclusive with NODES_J21I_NUM.

NODES_J22R_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the real part of the J22 radiation field tensor. Mutually exclusive with NODES_J22R_NUM.

NODES_J22I_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the imaginary part of the J22 radiation field tensor. Mutually exclusive with NODES_J22I_NUM.

BTYPE
-----

  * OPTIONAL

  * Formats:
    
    - string: Vertical, LOS; default: Vertical

  * Description: Invert the magnetic field strength, inclinarion, and azimuth in the vertical reference frame (vertical), or the longitudinal and transversal component, and the azimuth in the line of sight reference frame (LOS).

VTYPE
-----

  * OPTIONAL

  * Formats:
    
    - string: Vertical, LOS; default: Vertical

  * Description: Invert the velocity x, y, and vertical components in the vertical reference frame (vertical), or the longitudinal and transversal component, and the azimuth in the line of sight reference frame (LOS).

FIX_B
-----

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the magnetic field strength or longitudinal component.

FIX_BT
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the magnetic field inclination or transversal component.

FIX_BP
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the magnetic field azimuth.

FIX_F
-----

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the difusse light factor.

FIX_T
-----

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the temperature.

FIX_VX
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the velocity x or transversal component.

FIX_VY
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the velocity y component or azimuth.

FIX_VZ
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the velocity vertical or longitudinal component.

FIX_VT
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the micro-turbulent velocity.

FIX_PG
------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the gas pressure.

FIX_J21R
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the real part of the J21 radiation field tensor.

FIX_J21I
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the imaginary part of the J21 radiation field tensor.

FIX_J22R
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the real part of the J22 radiation field tensor.

FIX_J22I
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the imaginary part of the J22 radiation field tensor.

POSITION_CORRECTION
-------------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Correct the positions of the nodes to coincide with the closest node of the forward solver stratification.

REGUL_B
-------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the magnetic field strength or longitudinal component.

REGUL_BT
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the magnetic field inclination or longitudinal component.

REGUL_BP
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the magnetic field azimuth.

REGUL_F
-------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the diffuse light factor.

REGUL_T
-------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the temperature.

REGUL_VX
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the velocity x or transversal component.

REGUL_VY
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the velocity y component or azimuth.

REGUL_VZ
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the velocity vertical or longitudinal component.

REGUL_VT
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the micro-turbulent velocity.

REGUL_PG
--------

  * OPTIONAL

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the gas pressure.

REGUL_J21R
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the real part of the J21 radiation field tensor.

REGUL_J21I
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the imaginary part of the J21 radiation field tensor.

REGUL_J22R
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the real part of the J22 radiation field tensor.

REGUL_J22I
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the imaginary part of the J22 radiation field tensor.

REGUL_LIMITS
------------

  * OPTIONAL

  * Formats:
    
    - float; default: 0.1

  * Description: Upper limit of the relative weight that the regularization can have with respect to the merit function.

THRH_CHI2
---------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-4

  * Description: Threshold on the value of the merit function to stop the inversion.

INV_MRC
-------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-2

  * Description: Maximum relative change of the merit function to consider it converged.

SVD_TYPE
--------

  * OPTIONAL

  * Formats:
    
    - string: traditional, sir; default: sir

  * Description: Use the traditional singular value decomposition (traditional) or the one with the same modifications implemented in the SIR code (sir).

THRH_SVD
--------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-4

  * Description: Threshold for the singular value decomposition.

PG_TYPE
-------

  * OPTIONAL

  * Formats:
    
    - string: Hydrostatic equilibrium, stratified; default: hydrostatic equilibrium

  * Description: How to deal with the gas pressure stratification. The stratified option is NOT RECOMMENDED.

PG_BOUND
--------

  * OPTIONAL

  * Formats:
    
    - float; default: -1

  * Description: Initial or fixed value of the gas pressure at the top boundary, in dyn cm^-2, if PG_TYPE = hydrostatic. A negative value means that the gas pressure at the topmost node of the initial model atmosphere will be taken instead.

DIFFUSE_LIGHT
-------------

  * OPTIONAL

  * Formats:
    
    - float; default: -1

  * Description: Initial or fixed value of the diffuse light factor.

ATMO_NODES
----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes in the model atmosphere in the forward solution. A 0 value indicates that the stratification of the initial model atmosphere will be kept.

MAX_SVD_STEP
------------

  * OPTIONAL

  * Formats:
    ; float (1.0)

  * Description: Maximum step allowed in the singular value decomposition.

INV_ERROR
---------

  * OPTIONAL

  * Formats:
    
    - string: Hessian, RF, worst, recycle; default: Hessian

  * Description: Algorithm to estimate the error in the inversion. Recycle uses the last Hessian, Hessian recalculates the Hessian, RF is slightly different, using directly the response functions, and worst returns the worst between Hessian and RF.

PSF_FWHM
--------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float
    - float*3
    - float*2 + string: file path

  * Description: The FWHM of a gaussian spectral PSF in nanometers (only one entry). | Wavelength range in nm for which the specified FWHM must be used, and FWHM of a gaussian spectral PSF in nanometers. | Wavelength range in nm for which the specified PSF must be used and path, absolute or relative to the running directory, of a file with the spectral PSF. Note that the wavelength range does not truncate, but look for the closest wavelength in the data. Moreover, they are intended to be for distinguishing between different spectral ranges.

BOUNDS_B
--------

  * OPTIONAL

  * Formats:
    
    - float*2; default: 0 3000 (only BTYPE=vertical); default: -3000 3000 (only BTYPE=los)

  * Description: Minimum and maximum values that the magnetic field strengths or longitudinal component can take in the inversion.

BOUNDS_BT
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: 0 pi (only BTYPE=vertical); default: 0 3000 (only BTYPE=los)

  * Description: Minimum and maximum values that the magnetic field inclination or transversal component can take in the inversion.

BOUNDS_BP
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: 0 2pi

  * Description: Minimum and maximum values that the magnetic field azimuth can take in the inversion.

BOUNDS_F
--------

  * OPTIONAL

  * Formats:
    
    - float*2; default: 0 0.95

  * Description: Minimum and maximum values that the diffuse light factor can take in the inversion.

BOUNDS_VX
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: -20 20 (only vtype=vertical); default: 0 20 (only vtype=los)

  * Description: Minimum and maximum values that the velocity x or transversal component can take in the inversion.

BOUNDS_VY
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: -20 20 (only vtype=vertical); default: 0 2pi (only vtype=los)

  * Description: Minimum and maximum values that the velocity y component or azimuth can take in the inversion.

BOUNDS_VZ
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: -20 20

  * Description: Minimum and maximum values that the velocity vertical or longitudinal component can take in the inversion.

BOUNDS_VT
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: 0 40

  * Description: Minimum and maximum values that the micro-turbulent velocity can take in the inversion.

BOUNDS_PG
---------

  * OPTIONAL

  * Formats:
    
    - float*2; default: 0.1 15

  * Description: Minimum and maximum values that the gas pressure can take in the inversion.

BOUNDS_J21R
-----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the real part of the J21 radiation field tensor can take in the inversion.

BOUNDS_J21I
-----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the imaginary part of the J21 radiation field tensor can take in the inversion.

BOUNDS_J22R
-----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the real part of the J22 radiation field tensor can take in the inversion.

BOUNDS_J22I
-----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the imaginary part of the J22 radiation field tensor can take in the inversion.

EBOUNDS_B
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the magnetic field strengths or longitudinal component can take in the inversion.

EBOUNDS_BT
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the magnetic field inclination or transversal component can take in the inversion.

EBOUNDS_BP
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the magnetic field azimuth can take in the inversion.

EBOUNDS_VX
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the velocity x or transversal component can take in the inversion.

EBOUNDS_VY
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the velocity y component or azimuth can take in the inversion.

EBOUNDS_VZ
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the velocity vertical or longitudinal component can take in the inversion.

EBOUNDS_VT
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the micro-turbulent velocity can take in the inversion.

EBOUNDS_PG
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the gas pressure can take in the inversion.

EBOUNDS_J21R
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the real part of the J21 radiation field tensor can take in the inversion.

EBOUNDS_J21I
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the imaginary part of the J21 radiation field tensor can take in the inversion.

EBOUNDS_J22R
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the real part of the J22 radiation field tensor can take in the inversion.

EBOUNDS_J22I
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the imaginary part of the J22 radiation field tensor can take in the inversion.

SCALE_B
-------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 100

  * Description: Scale factor for the magnetic field strengths or longitudinal component can take in the inversion.

SCALE_BT
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1.8 (only BTYPE=vertical); default: 100 (only BTYPE=los)

  * Description: Scale factor for the magnetic field inclination or transversal component can take in the inversion.

SCALE_BP
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1.8

  * Description: Scale factor for the magnetic field azimuth can take in the inversion.

SCALE_F
-------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1

  * Description: Scale factor for the diffuse light factor can take in the inversion.

SCALE_VX
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 10

  * Description: Scale factor for the velocity x or transversal component can take in the inversion.

SCALE_VY
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 10 (only vtype=vertical); default: 1.8 (only vtype=los)

  * Description: Scale factor for the velocity y component or azimuth can take in the inversion.

SCALE_VZ
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 10

  * Description: Scale factor for the velocity vertical or longitudinal component can take in the inversion.

SCALE_VT
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 10

  * Description: Scale factor for the micro-turbulent velocity can take in the inversion.

SCALE_PG
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 2

  * Description: Scale factor for the gas pressure can take in the inversion.

SCALE_J21R
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the real part of the J21 radiation field tensor can take in the inversion.

SCALE_J21I
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the imaginary part of the J21 radiation field tensor can take in the inversion.

SCALE_J22R
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the real part of the J22 radiation field tensor can take in the inversion.

SCALE_J22I
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the imaginary part of the J22 radiation field tensor can take in the inversion.

PERTURB_B
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1

  * Description: Value of the perturbation to compute the response function for the magnetic field strengths or longitudinal component can take in the inversion.

PERTURB_BT
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-2 (only BTYPE=vertical); default: 2 (only BTYPE=los)

  * Description: Value of the perturbation to compute the response function for the magnetic field inclination or transversal component can take in the inversion.

PERTURB_BP
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3e-2

  * Description: Value of the perturbation to compute the response function for the magnetic field azimuth can take in the inversion.

PERTURB_F
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-2

  * Description: Value of the perturbation to compute the response function for the diffuse light factor can take in the inversion.

PERTURB_VX
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3e-1

  * Description: Value of the perturbation to compute the response function for the velocity x or transversal component can take in the inversion.

PERTURB_VY
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3e-1 (only vtype=vertical); default: 3e-2 (only vtype=los)

  * Description: Value of the perturbation to compute the response function for the velocity y component or azimuth can take in the inversion.

PERTURB_VZ
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3e-1

  * Description: Value of the perturbation to compute the response function for the velocity vertical or longitudinal component can take in the inversion.

PERTURB_VT
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 3e-1

  * Description: Value of the perturbation to compute the response function for the micro-turbulent velocity can take in the inversion.

PERTURB_PG
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 5e-2

  * Description: Value of the perturbation to compute the response function for the gas pressure can take in the inversion.

PERTURB_J21R
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the real part of the J21 radiation field tensor can take in the inversion.

PERTURB_J21I
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the imaginary part of the J21 radiation field tensor can take in the inversion.

PERTURB_J22R
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the real part of the J22 radiation field tensor can take in the inversion.

PERTURB_J22I
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the imaginary part of the J22 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_B
-----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the magnetic field strengths or longitudinal component can take in the inversion.

MIN_REL_PERTURB_BT
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the magnetic field inclination or transversal component can take in the inversion.

MIN_REL_PERTURB_BP
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the magnetic field azimuth can take in the inversion.

MIN_REL_PERTURB_F
-----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the diffuse light factor can take in the inversion.

MIN_REL_PERTURB_VX
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the velocity x or transversal component can take in the inversion.

MIN_REL_PERTURB_VY
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the velocity y component or azimuth can take in the inversion.

MIN_REL_PERTURB_VZ
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the velocity vertical or longitudinal component can take in the inversion.

MIN_REL_PERTURB_VT
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the micro-turbulent velocity can take in the inversion.

MIN_REL_PERTURB_PG
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the gas pressure can take in the inversion.

MIN_REL_PERTURB_J21R
--------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the real part of the J21 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_J21I
--------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the imaginary part of the J21 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_J22R
--------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the real part of the J22 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_J22I
--------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0

  * Description: Minimum relative perturbation to compute the response function for the imaginary part of the J22 radiation field tensor can take in the inversion.

INI_BPOS
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0.17 (only BTYPE=vertical); default: 10 (only BTYPE=los)

  * Description: Initial magnetic field inclination or transversal component when initializing from an inversion result with a very small value of this quantity.

INI_BAZI
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 0.5

  * Description: Initial magnetic field azimuth when initializing from an inversion result with a very small value of this quantity.

INI_VPOS
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 5

  * Description: Initial velocity x or transversal component when initializing form an inversion result with a very small value of this quantity.

INI_VAZI
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 1.5

  * Description: Initial velocity y or azimuth when initializing form an inversion result with a very small value of this quantity.

INV_FRACTION
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fit the fractional polarization (Yes) instead of the Stokes profiles (No).

INV_TAU_RANG
------------

  * OPTIONAL

  * Formats:
    
    - float*2; default: -8 1

  * Description: Logarithmic optical depth of the extreme nodes of the model atmosphere in the forward solution. This has no effect if ATMO_NODES is not specified.

ATMO_STRAT
----------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float*3

  * Description: Triplets of number, with the first two specifying a range in logarithmic optical depth in which the sampling step of the atmospheric model for the forward solver must be reduced by a factor given by the third number.

BROYDEN_LM
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Use Broyden's algorithm in the profile fitting.

LM_LAMBDA_RANG
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float*2; default: 1e-5 1e3

  * Description: Boundary limits for the lambda coefficient in the Levenmberg-Marquardt algorithm.

LM_LAMBDA_ACCEPT
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 5

  * Description: Factor dividing the lambda coefficient in the Levenberg-Marquardt algorithm when the previous lambda value was accepted.

LM_LAMBDA_REJECT
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 5

  * Description: Factor multiplying the lambda coefficient in the Levenberg-Marquardt algorithm when the previous lambda value was rejected.

INV_B_PROJECTION
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Make sure to correctly project the negative polarity when the magnetic field is in the vertical reference frame but you want to use it as its longitudinal component.

RF_INITSOL
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Initialize the calculation of the response function with the solution of the reference profile.

INV_NEGL_SIGMA
--------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Neglect the sigma (noise) in the observation for the merit function.

KEEP_RF
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Include the numerical response functions in the solution file.

FORCE_OBS_FREQ
--------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force the wavelengths of the data to be present in the frequency axis of the forward solver.

