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

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the atomic model whose density matrix elements are to be calculated, in principle, in non-LTE. This keyword in MANDATORY in CLE mode. In other models, at least one entry of ATOM_INPUT or LTE_LINE must exist. Optional: Unique label to identify this atom in other input options.

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

  * Description: The model atom with the corresponding label will keep their initial populations fixed when solving the iterative problem. These populations are calculated in LTE unless ATOM_POPU is specified for the same model atom.

ATOM_FIX_POP_LTERM
------------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will keep the initial populations of the ground level/term fixed when solving the iterative problem. These populations are calculated in LTE unless ATOM_POPU is specified for the same model atom.

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

  * Description: Entry with the data of an atomic transition to be included under the assumption of LTE. The format can be the same from the SIR code (include the full SIR line after 'LTE_LINE =', including the line index and the equal sign, even though the index is not used), from the Kurucz's database (fill with zeros to fulfill the restricted size requirements), or the specific format of the code. | Path, absolute or relative to the running directory, to the file with a list of atomic transitions to be included under the assumption of LTE. The format of each entry in the file can be the same from the SIR code (include the full SIR line, including the line index and the equal sign, even though the index is not used), from the Kurucz's database (fill with zeros to fulfill the restricted size requirements) or the specific format of the code. While optional, at least one entry of ATOM_INPUT or LTE_LINE must exist.

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

TWO_STEP_INTENSITY_V
--------------------

  * OPTIONAL

  * Formats:
    ; string: Yes/No

  * Description: Solve the intensity problem in two steps, by first converging the problem without any velocity and then switching on the velocity field. WARNING: It does not work if a Solution file needs to be loaded.

RESTRICT_T_BOT_STRICT
---------------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the height in RESTRICT_T_BOT. If not, once the restriction is identified, the selected node will move one step to the extrema.

RESTRICT_T_BOT
--------------

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Lower limit to the temperature from the bottom of the model to solve the radiative transfer problem. Contiguous nodes with values of the temperature larger than the specified toward the bottom will be neglected.

RESTRICT_T_UP_STRICT
--------------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the height in RESTRICT_T_UP. If not, once the restriction is identified, the selected node will move one step to the extrema.

RESTRICT_T_UP
-------------

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Lower limit to the temperature from the top of the model to solve the radiative transfer problem. Contiguous nodes with values of the temperature larger than the specified toward the top will be neglected.

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

ANISOTROPY_FOCUS
----------------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: When computing the maximum relative change to decide for convergence, consider only the K=0 and K=2 multipoles.

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

  * Description: Maximum amount of Megabytes which can be allocated per CPU, used to limit the amount of Voigt and redistribution profiles to store in RAM. The counting is not perfect, so be conservative if there is a RAM limit.

RED_RESTRICT_HEIGHT
-------------------

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Upper limit in the decimal logarithm of the optical depth where to consider partial frequency distribution effects.

RED_RESTRICT_TAUC
-----------------

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Lower limit in height in kilometers where to consider partial frequency distribution effects.

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

TWO_STEP_AD
-----------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No

  * Description: Solve the PRD-AD problem in two steps, first converging the PRD-AA problem and then switching PRD-AD on. WARNING: It does not work if a Solution file needs to be loaded.

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

MIN_T
-----

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Minimum temperature in kelvin expected in the 3D model atmosphere. If not specified, the model will be explored to get it.

MAX_T
-----

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Maximum temperature in kelvin expected in the 3D model atmosphere. If not specified, the model will be explored to get it.

MAX_V
-----

  * OPTIONAL

  * Formats:
    
    - float

  * Description: Maximum velocity in km s^-1 expected in the 3D model atmosphere. If not specified, the model will be explored to get it.

RT_GROUP_N
----------

  * OPTIONAL

  * Formats:
    
    - integer; default: 1

  * Description: Number of processes to solve the forward or inversion problem for each pixel of the input.

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

  * Description: Maximum iteration index allowed for in a self-consistent solution.

ITERI_MAX
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: ITER_MAX

  * Description: Maximum iteration index allowed for in an only intensity self-consistent solution.

ITERAD_MAX
----------

  * OPTIONAL

  * Formats:
    
    - integer; default: ITER_MAX

  * Description: Maximum iteration index allowed for in the self-consistent solution for the PRD-AD phase if TWO_STEP_AD is activated.

ITERIAD_MAX
-----------

  * OPTIONAL

  * Formats:
    
    - integer; default: ITERI_MAX

  * Description: Maximum iteration index allowed for in an only intensity self-consistent solution for the PRD-AD phase if TWO_STEP_AD is activated.

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

ITER_MRC_ADI
------------

  * OPTIONAL

  * Formats:
    
    - float; default: ITER_MRC_I

  * Description: Maximum relative change of the populations to consider that they have converged in the PRD-AD phase if TWO_STEP_AD is activated.

ITERI_MRC_ADI
-------------

  * OPTIONAL

  * Formats:
    
    - float; default: ITERI_MRC_I

  * Description: Maximum relative change of the populations to consider that they have converged in the only intensity solution in the PRD-AD phase if TWO_STEP_AD is activated.

ITER_MRC_ADP
------------

  * OPTIONAL

  * Formats:
    
    - float; default: ITER_MRC_P

  * Description: Maximum relative change of the non-population density matrix elements to consider that they have converged in the PRD-AD phase if TWO_STEP_AD is activated.

ITER_J
------

  * OPTIONAL

  * Formats:
    
    - integer; default: 5

  * Description: Number of preliminar no-line iterations to perform to relax the initial radiation field.

ITERI_PRD
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: 4

  * Description: Maximum number of only-radiation iterations to perform when there is partial frequency redistribution in the only intensity problem.

ITERI_MRC_R
-----------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-3

  * Description: Maximum relative change of the mean intensity in only-radiation iterations to consider that it is converged in the intensity problem.

ITER_PRD
--------

  * OPTIONAL

  * Formats:
    
    - integer; default: 1

  * Description: Maximum number of only-radiation iterations to perform when there is partial frequency redistribution.

ITER_MRC_R
----------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-2

  * Description: Maximum relative change of the mean intensity in only-radiation iterations to consider that it is converged.

ITER_MRC_P_R
------------

  * OPTIONAL

  * Formats:
    
    - float; default: 1e-1

  * Description: Maximum relative change of the radiation field tensors (not mean intensity) in only-radiation iterations to consider that it is converged.

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

SOLUTION_BOX
------------

  * OPTIONAL

  * Formats:
    
    - integer*4; default: -1 -1 -1 -1

  * Description: Indicate the initial x index, final x index, initial y index, and final y index, respectively, of the pixels to solve in a 3D model or data file. Negative numbers are wildcards (automatically adjusted to the relevant size of the input).

EXCLUDE_PIXEL
-------------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - integer*2

  * Description: Pair of X and Y pixel coordinate for a pixel that must be excluded from the calculations. Note that, for the inversion model, this is not equivalent to setting a mask with INV_MASK, as the pixel will be completely skipped and nothing will be written.

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

LIM_QEL
-------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string + integer

  * Description: Unique identifier of the model atom and index specifying the bound-bound transition whose elastic rate will be included in the output for elastic rates. Transitions not included in the specified ranges will not be included in the output.

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

VERBOSE
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Indicate if the verbosity output of the code must go to the command line (Yes) or to a verbosity file to be created in the output directory (No).

