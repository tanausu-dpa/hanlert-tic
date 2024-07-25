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

BFIELD_INPUT
------------

  * OPTIONAL

  * Formats:
    
    - string: file path; default: None (only 1D)
    - float*3

  * Description: Path, absolute or relative to the running directory, of the file with the stratification of the magnetic field vector. | Assume an homogeneous magnetic field with the strength (in Gauss), inclination (in degrees), and azimuth (in degrees) specified by such three floats, respectively.

OPACITY_FUDGE
-------------

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the fudge factors for the background continuum opacity.

ABUND
-----

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with atomic abundances.

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

  * Description: List of floats with the cosine of the heliocentric angle for which to compute the Stokes parameters formal solution. If both POLAR_LOS and AXIAL_LOS are absent, no formal solution is carried out and outputs for specific lines of sight are not generated.

AXIAL_LOS
---------

  * OPTIONAL

  * Formats:
    
    - floats

  * Description: List of floats with the azimuthal angle, in degrees, for which to compute the Stokes parameters formal solution. If both POLAR_LOS and AXIAL_LOS are absent, no formal solution is carried out and outputs for specific lines of sight are not generated.

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

RESTRICT_TAUC_STRICT
--------------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the ranges in RESTRICT_TAUC. If not, once the restriction is identified, the selected nodes will move one step to the extrema.

RESTRICT_TAUC
-------------

  * OPTIONAL

  * Formats:
    
    - float*2

  * Description: Lower and upper limit in the decimal logarithm of the optical depth where to solve the radiative transfer problem. Nodes with values of the optical depths outside the specified range will be neglected.

RESTRICT_HEIGHT_STRICT
----------------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the ranges in RESTRICT_HEIGHT. If not, once the restriction is identified, the selected nodes will move one step to the extrema.

RESTRICT_HEIGHT
---------------

  * OPTIONAL

  * Formats:
    
    - float*2

  * Description: Lower and upper limit in height in kilometers where to solve the radiative transfer problem. Nodes with values of the height outside the specified range will be neglected.

MIT_OFF
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No, none; default: none

  * Description: Neglect (Yes) or force the inclusion (No) of frequencies corresponding to magnetically induced transitions in multi-term model atoms. By default, they are only included if there is a magnetic field.

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

SOLUTION_BACKUP
---------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Rename a Solution file existing in the output directory (Yes) to avoid overwritting.

K_CUT
-----

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum density matrix multipolar component to consider in the polarized problem. Negative means no restriction.

K_RAD
-----

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum radiation field tensor multipolar component to consider in the polarized problem. Negative means no restriction.

VOI_IRAM
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the Voigt profiles in the only intensity problem.

VOI_PRAM
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in RAM the complex Voigt profiles.

RAM_LIMIT
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum amount of Megabytes which can be allocated in the form of Voigt profiles, photoionization pre-calculated quantities, interpolation precalculated quantities (partial frequency redistribution), or redistribution functions. Negative means no limit (NOT RECOMMENDED for complex problems).

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

FORCE_MICRO
-----------

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Constant microturbulence, in kilometers per second, to force in the model atmosphere.

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

ITER_MRC_J
----------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-4

  * Description: Maximum relative change of the radiation field to consider it converged in the only radiation field initial iterations.

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

STORE_STEP
----------

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: If positive, a solution file will be created when the iteration index in the polarized problem is a multiple of this input.

STOREI_STEP
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: If positive, a solution file will be created when the iteration index in the only intensity problem is a multiple of this input.

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

KEEP_QEL
--------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the elastic rates characteristic of each transition for the ATOM_FILE model atoms.

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

VERBOSE
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Indicate if the verbosity output of the code must go to the command line (Yes) or to a verbosity file to be created in the output directory (No).

