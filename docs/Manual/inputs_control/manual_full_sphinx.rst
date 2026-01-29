RUN_MODE
--------

  * MANDATORY

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: 1DS, 1.5DS, CLE, inversion; default: 1DS

  * Description: Indicate the type of problem to solve, a single 1D solution (1DS), a 1.5D solution of a 3D model (1.5DS), a 1.5D formal solution in an optically thin plasma (CLE), or pixel by pixel inversion (inversion).

ATMO_INPUT
----------

  * MANDATORY

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path, HCC (only INV), HCP (only INV), HCA (only INV), HCF (only INV), HCX (only INV); default: HCC (only INV)

  * Description: Path, absolute or relative to the running directory, of the model atmosphere file. In inversion mode this specifies the initial model atmosphere and HCC, HCP, HCA, HCF, and HCX can be specified to indicate the C, P, A, or F models of Fontenla et al. (1993) or the MCO (a.k.a. FALX) model of Ayres (1986), respectively.

DATA_FILE
---------

  * MANDATORY

  * Modes: INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the data to invert.

ATMO_SCALE
----------

  * OPTIONAL

  * Modes: 15D, INV

  * Formats:
    
    - string: Height, tau; default: Height; optional: float (default 500)

  * Description: Type of stratification scale of the input model atmosphere (if 15D) or the initial 3D model atmosphere (if INV). Optional: reference wavelength for the optical depth calculation, in nanometers.

ATMO_CHAR
---------

  * OPTIONAL

  * Modes: 15D, CLE, INV

  * Formats:
    
    - string: NENH, NE, PE, RHOE, PG, RHO; default: NENH

  * Description: Variable which variable the model atmosphere is providing to characterize the stratification, together with the temperature, both electron and Hydrogen number densities (NENH), electron number densities (NE), electron pressure in dyn cm^-2 (PE), electron density in g cm^-3 (RHOE), gas pressure in dyn cm^-2 (PG), or total density in g cm^-3 (RHO).

RESPECT_ALT_SCALE
-----------------

  * OPTIONAL, ADVANCED

  * Modes: 15D

  * Formats:
    
    - string: Yes/No

  * Description: If the primary provided stratification is geometrical height (optical depth), but the secondary optical depth (geometrical height) axis is also provided in the atmospheric model, skip the internal calculation. This is particularly useful when using keywords to limit the domain of the radiation transfer problem relying in such secondary axis.

ATOM_INPUT
----------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the atomic model whose density matrix elements are to be calculated, in principle, in non-LTE. This keyword in MANDATORY in CLE mode. In other models, at least one entry of ATOM_INPUT or LTE_LINE must exist. Optional: Unique label to identify this atom in other input options.

ATOM_BACK
---------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the atomic model to include in the calculation as a background opacity source. Their populations will be calculated in LTE if no population file is provided (see ATOM_POPU).

ATOM_POPU
---------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the initial population (if the label points to a model in ATOM_FILE) or the fixed populations (if the label points to a model in ATOM_BACK) for the model atom specified by the label.

ATOM_ION
--------

  * OPTIONAL, ADDITIVE

  * Modes: CLE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the ionization fraction for the model atom specified by the label.

ATOM_FIX_POP
------------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will keep their initial populations fixed when solving the iterative problem. These populations are calculated in LTE unless ATOM_POPU is specified for the same model atom.

ATOM_FIX_POP_LTERM
------------------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will keep the initial populations of the ground level/term fixed when solving the iterative problem. These populations are calculated in LTE unless ATOM_POPU is specified for the same model atom.

ATOM_ZERO_ION
-------------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will ignore the transitions with the last level/term in the model atom (of a different ionization stage than the rest of levels/terms), effectively keeping its initial ionization balance.

ATOM_NO_WAVE
------------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: a single word and unique label

  * Description: The model atom with the corresponding label will not contribute with frequency nodes when building the output frequency/wavelength axis.

MOLECULE_INPUT
--------------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the molecular data that will be used to compute the chemical equilibrium and the contribution of some molecules to the background opacity.

BFIELD_INPUT
------------

  * OPTIONAL

  * Modes: 1D, INV

  * Formats:
    
    - string: file path; default: None (only 1D)
    - float*3

  * Description: Path, absolute or relative to the running directory, of the file with the stratification of the magnetic field vector. | Assume an homogeneous magnetic field with the strength (in Gauss), inclination (in degrees), and azimuth (in degrees) specified by such three floats, respectively.

OPACITY_FUDGE
-------------

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the fudge factors for the background continuum opacity.

SPECT_INPUT
-----------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the spectra for the incoming radiation.

USE_ALLEN
---------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string: Yes/No

  * Description: Use Allen's tabulation to determine the intensity of the star.

CLV
---

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string: None/Allen

  * Description: Uses Allen's tabulation to determine the intensity CLV or neglect CLV. This is only used if there is no SPECT_INPUT data for a needed frequency.

FLAT_CLE_IN
-----------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string: Yes/No

  * Description: Assume that the radiation is exactly flat when computing the JKQ (only when no input spectra is provided).

CHIANTI_PATH
------------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string: folder path

  * Description: Path, absolute or relative to the running directory, to the folder with the CHIANTI database data.

T_RAD
-----

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - float; default: 5e3

  * Description: Effective temperature  in kelvin of the black-body radiation for the boundary condition if SPECT_INPUT is not specified or for frequencies not included in such file.

R_STAR
------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - float; default: 6.957e10

  * Description: Stellar radius in cm.

NEGLECT_CONTINUUM
-----------------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string Yes, No; default: No

  * Description: Neglect continuum contribution to the radiative transfer coefficients.

PARFUN
------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with partition function and ionization potential tabulations.

ABUND
-----

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with atomic abundances.

BARK_SP
-------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with the tabulation of the Barklem broadening parameters for s-p transitions.

BARK_PD
-------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with the tabulation of the Barklem broadening parameters for p-d transitions.

BARK_DF
-------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with the tabulation of the Barklem broadening parameters for d-f transitions.

IGNORE_BB
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No

  * Description: Neglect (Yes) or include (No) the bound-bround transitions of the ATOM_BACK atomic models in the formal solution.

KURUCZ
------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with a list of transitions in the Kurucz's database format to be included in the background opacity.

LTE_LINE
--------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, INV

  * Formats:
    
    - entry
    - string: file path

  * Description: Entry with the data of an atomic transition to be included under the assumption of LTE. The format can be the same from the SIR code (include the full SIR line after 'LTE_LINE =', including the line index and the equal sign, even though the index is not used), from the Kurucz's database (fill with zeros to fulfill the restricted size requirements), or the specific format of the code. | Path, absolute or relative to the running directory, to the file with a list of atomic transitions to be included under the assumption of LTE. The format of each entry in the file can be the same from the SIR code (include the full SIR line, including the line index and the equal sign, even though the index is not used), from the Kurucz's database (fill with zeros to fulfill the restricted size requirements) or the specific format of the code. While optional, at least one entry of ATOM_INPUT or LTE_LINE must exist.

WAVELENGTHS
-----------

  * OPTIONAL, ADDITIVE

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to a file with a list of wavelengths to include in the frequency/wavelength axis in the formal solution.

ASYMM_INPUT
-----------

  * OPTIONAL, ADDITIVE, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer*2 float*2 (only 1D)
    - string: file path

  * Description: Multipolas component idenfifiers K and Q, and real and imaginary values of the JKQ tensor to be considered as an ad-hoc (homogeneous) asymmetry in the formal solution. | Path, absolute or relative to the running directory, to a file with JKQ tensors to be considered as an ad-hoc assymetry in the formal solution.

FORCE_SYMM
----------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the input asymmetry must be added to the one properly generated in the formal solution (No) or must completely overwrite the actual radiation field tensors (Yes).

STIM
----

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: If stimulated emission must be accounted for.

OUT_FOLDER
----------

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: folder path; default: Outputs/Default

  * Description: Path, absolute or relative to the running directory, to the folder to store the outputs. This folder does not have to exist, but must be possible to create it (the path up to the second to last folder must exist).

POLAR_NODES
-----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 4

  * Description: Number of nodes per hemisphere in the polar (gaussian) angular quadrature for the solution of the radiative transfer problem.

AXIAL_NODES
-----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 2

  * Description: Number of nodes per octant in the azimuthal (trapezoidal) angular quadrature for the solution of the radiative transfer problem. A negative value indicates axial symmetry.

POLARI_NODES
------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: POLAR_NODES

  * Description: Number of nodes per hemisphere in the polar (gaussian) angular quadrature for the solution of the only intensity radiative transfer problem.

AXIALI_NODES
------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: AXIAL_NODES

  * Description: Number of nodes per octant in the azimuthal (trapezoidal) angular quadrature for the solution of the only intensity radiative transfer problem. A negative value indicates axial symmetry.

POLAR_LOS
---------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - floats

  * Description: List of floats with the cosine of the heliocentric angle for which to compute the Stokes parameters formal solution. If both POLAR_LOS and AXIAL_LOS are absent, no formal solution is carried out and outputs for specific lines of sight are not generated.

AXIAL_LOS
---------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - floats

  * Description: List of floats with the azimuthal angle, in degrees, for which to compute the Stokes parameters formal solution. If both POLAR_LOS and AXIAL_LOS are absent, no formal solution is carried out and outputs for specific lines of sight are not generated.

MODE
----

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: solve, read, both

  * Description: Specify what to do in the formal solution, solving from standard initialization (solve), read a previous solution and continue iterating (both), or read a previous solution and just perform the last formal solutions (read). WARNING: When using both with an only intensity solution file, it will be assumed that such a solution is fully converged.

FORCE
-----

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: intensity, polarization, all, none; default: none

  * Description: Force the code to solve only the intensity problem (intensity), to solve the polarization problem skipping the intensity problem (polarization), to solve both intensity and polarization problem (all), or let the code decide from the available inputs (none).

TWO_STEP_INTENSITY_V
--------------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    ; string: Yes/No

  * Description: Solve the intensity problem in two steps, by first converging the problem without any velocity and then switching on the velocity field. WARNING: It does not work if a Solution file needs to be loaded.

RESTRICT_T_BOT_STRICT
---------------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the height in RESTRICT_T_BOT. If not, once the restriction is identified, the selected node will move one step to the extrema.

RESTRICT_T_BOT
--------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - float

  * Description: Lower limit to the temperature from the bottom of the model to solve the radiative transfer problem. Contiguous nodes with values of the temperature larger than the specified toward the bottom will be neglected.

RESTRICT_T_UP_STRICT
--------------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the height in RESTRICT_T_UP. If not, once the restriction is identified, the selected node will move one step to the extrema.

RESTRICT_T_UP
-------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - float

  * Description: Lower limit to the temperature from the top of the model to solve the radiative transfer problem. Contiguous nodes with values of the temperature larger than the specified toward the top will be neglected.

RESTRICT_TAUC_STRICT
--------------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the ranges in RESTRICT_TAUC. If not, once the restriction is identified, the selected nodes will move one step to the extrema.

RESTRICT_TAUC
-------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - float*2

  * Description: Lower and upper limit in the decimal logarithm of the optical depth where to solve the radiative transfer problem. Nodes with values of the optical depths outside the specified range will be neglected.

RESTRICT_HEIGHT_STRICT
----------------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No

  * Description: If yes, the nodes will be truncated to not include the ranges in RESTRICT_HEIGHT. If not, once the restriction is identified, the selected nodes will move one step to the extrema.

RESTRICT_HEIGHT
---------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - float*2

  * Description: Lower and upper limit in height in kilometers where to solve the radiative transfer problem. Nodes with values of the height outside the specified range will be neglected.

ZEEMAN_MODE
-----------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: full, nozeeman, nosplit, linear; default: full

  * Description: Type of diagonalization of the atomic Hamiltonian, full diagonalization (full), neglect any magnetic splitting (nozeeman), neglect the splitting but not the level coupling (nosplit), or assume linear regime (linear).

P_CORR
------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Add one fake iteration between the solution of the only intensity and full Stokes radiative transfer problems to adjust the radiation field tensors of multi-term atomic models.

FCOL_TRANSFER
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Allow forbidden collisions to transfer polarization (Yes).

MIT_OFF
-------

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No, none; default: none

  * Description: Neglect (Yes) or force the inclusion (No) of frequencies corresponding to magnetically induced transitions in multi-term model atoms. By default, they are only included if there is a magnetic field.

MIT_NODE
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float

  * Description: Factor for which to divide the number of nodes in the frequency/wavelength axis to allocate for the MIT with respect to the permitted transitions of a multiplet.

SKIP_SFILE
----------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the writing of the solution file (which allows initializing other runs) should be skipped (Yes) or not (No).

SOLUTION_INPUT
--------------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: file path (only 1D)
    - string: folder path (only 15D)

  * Description: Path, absolute or relative to the running directory, to the Solution file of a previous run to be read to initialize this run for appropriate modes. | Path, absolute or relative to the running directory, to the Solution folder of a previous run to be read to initialize this run for appropriate modes.

SOLUTION_BACKUP
---------------

  * OPTIONAL

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Rename a Solution file existing in the output directory (Yes) to avoid overwritting.

SOLUTION_KEEPI
--------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Keep a solution file with the solution of the intensity problem even if solving the polarized problem afterwards.

SOLUTION_KEEPS
--------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force the solution file to store the full Stokes parameters instead of the radiation field tensors independently of the parameters which decide it automatically.

ANISOTROPY_FOCUS
----------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: When computing the maximum relative change to decide for convergence, consider only the K=0 and K=2 multipoles.

K_CUT
-----

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: -1

  * Description: Maximum density matrix multipolar component to consider in the polarized problem. Negative means no restriction.

K_CUTAB
-------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the K_CUT must be accounted for in the combination of multipolar components which appears when partial frequency redistribution is taken into account.

K_CUT_TERM
----------

  * OPTIONAL, ADDITIVE, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - 3*integers; 4*integers

  * Description: If three integers, label for the (active) model atom, term index, and maximum rank of the density matrix allowed for that term. If four integers, label for the (active) model atom, initial and final indexes for a range of terms, and maximum rank of the density matrix allowed for those terms.

K_RAD
-----

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: -1

  * Description: Maximum radiation field tensor multipolar component to consider in the polarized problem. Negative means no restriction.

MEMOJ
-----

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Specify if memoization is to be used to avoid calculating any J-symbol more than once.

PIRAM
-----

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM some quantities used in photoionization to avoid recalculating them over and over.

VOI_TYPE
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: H82, A67, HAW78, GAUSS; default: H82

  * Description: Algorithm to compute the Voigt profiles. If solving for polarization, only H82 is allowed.

VOI_IRAM
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the Voigt profiles in the only intensity problem.

LTE_VOI_IRAM
------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the complex Voigt profiles for LTE transitions.

VOI_PRAM
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in RAM the complex Voigt profiles.

RAM_LIMIT
---------

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: -1

  * Description: Maximum amount of Megabytes which can be allocated per CPU, used to limit the amount of Voigt and redistribution profiles to store in RAM. The counting is not perfect, so be conservative if there is a RAM limit.

RAM_REPORT
----------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Create a file with a report of the RAM used by different big consumers of memory.

RAMAN
-----

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Take into account Raman scattering (sometimes referred as cross redistribution, XRD).

NO_COH_L_TERM
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Asumme a non-coherent lower term when computing the radiative transfer coefficients.

RED_RESTRICT_HEIGHT
-------------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float

  * Description: Upper limit in the decimal logarithm of the optical depth where to consider partial frequency distribution effects.

RED_RESTRICT_TAUC
-----------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float

  * Description: Lower limit in height in kilometers where to consider partial frequency distribution effects.

RED_COHW
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float / string: No; Default: No

  * Description: Doppler widths from the line center from where to assume that the scattering is fully coherent.

REDI_COHW
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float / string: No; Default: RED_COHW

  * Description: Doppler widths from the line center from where to assume that the scattering is fully coherent in only intensity problems.

RED_INT_MODE
------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Linear, Splines; default: Splines

  * Description: Type of interpolation of the second order emissivity when transforming from the comoving to the observer's reference frame.

RED_MOD
-------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: AA, AD; default: AA

  * Description: Type of redistribution function, angle-average (AA) or angle-dependent (AD)

RED_AAINT
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force angle-average redistribution function for the only intensity problem, regardless of the RED_MOD input.

TWO_STEP_AD
-----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No

  * Description: Solve the PRD-AD problem in two steps, first converging the PRD-AA problem and then switching PRD-AD on. WARNING: It does not work if a Solution file needs to be loaded.

RED_IRAM
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: To store in RAM the redistribution function of the only intensity problem.

RED_PRAM
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: To store in RAM the complex redistribution function.

RED_NODE
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: 4

  * Description: Number of nodes per hemisphere in the scattering (gaussian) angular quadrature to compute the angle-averaged redistribution function.

REDI_NODE
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: RED_NODE

  * Description: Number of nodes per hemisphere in the scattering (gaussian) angular quadrature to compute the angle-averaged redistribution function in the only intensity problem.

RED_RANG
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 3.5

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral.

RED_RESO
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 3.5

  * Description: Distance in Doppler width from a transition resonance to move the search of input frequencies to include the full transition resonance.

RED_NEGL
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 1e4

  * Description: Distance in Doppler width from a frequency to the closest transition resonance in order to neglect partial frequency redistribution at that frequency. NOTE: The Doppler width in this variable is the one used to define the output frequency axis, and not the one corresponding to the atmospheric properties at each location.

RED_VLAR
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 7

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral.

RED_FSTP
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 0.5

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances.

RED_MSTP
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 4

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances.

RED_CORE
--------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 0

  * Description: Doppler width distance to the closest transition resonance to consider a frequency as in a line core.

RED_RANG_CORE
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_RANG

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral, for frequencies in a line core.

RED_VLAR_CORE
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_VLAR

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral, for frequencies in a line core.

RED_FSTP_CORE
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_FSTP

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances, for frequencies in a line core.

RED_MSTP_CORE
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_MSTP

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances, for frequencies in a line core.

REDI_RANG
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_RANG

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral in the only intensity problem.

REDI_RESO
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_RESO

  * Description: Distance in Doppler width from a transition resonance to move the search of input frequencies to include the full transition resonance in the only intensity problem

REDI_NEGL
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_NEGL

  * Description: Distance in Doppler width from a frequency to the closest transition resonance in order to neglect partial frequency redistribution at that frequency in the only intensity problem.

REDI_VLAR
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_VLAR

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral in the only intensity problem.

REDI_FSTP
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_FSTP

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances in the only intensity problem.

REDI_MSTP
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_MSTP

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances in the only intensity problem.

REDI_CORE
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: RED_CORE

  * Description: Doppler width distance to the closest transition resonance to consider a frequency as in a line core in the only intensity problem.

REDI_RANG_CORE
--------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: REDI_RANG

  * Description: Number of Doppler widths around transition and redistribution resonances to find the limits of the input frequencies in the redistribution integral, for frequencies in a line core in the only intensity problem.

REDI_VLAR_CORE
--------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: REDI_VLAR

  * Description: Additional Doppler width distance to include beyond the extremes of the input frequency ranges in the redistribution integral, for frequencies in a line core in the only intensity problem.

REDI_FSTP_CORE
--------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: REDI_FSTP

  * Description: Sampling in Doppler widths of the input frequency axis for the redistribution integral in a RED_RANG range from the transition and redistribution resonances, for frequencies in a line core in the only intensity problem.

REDI_MSTP_CORE
--------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: REDI_MSTP

  * Description: Multiplicative factor to apply to RED_FST to determine the sampling in Doppler widths of the input frequency axis for the redistribution integral for frequencies between RED_RANG and RED_VLAR Doppler widths from the transition and redistribution resonances, for frequencies in a line core in the only intensity problem.

DOP_WIDTH
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - float; default: 2.5e3
    - string: Max, Min

  * Description: Doppler width in m s^-1 to consider when building the frequency axis to transform between Doppler widths and actual frequencies. | Take the maximum (Max) or minimum (Min) temperature to calculate the Doppler width used to build the frequency axis to transform between Doppler widths and actual frequencies.

FORCE_MICRO
-----------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - float

  * Description: Constant microturbulence, in kilometers per second, to force in the model atmosphere.

MIN_T
-----

  * OPTIONAL

  * Modes: 15D, CLE, INV

  * Formats:
    
    - float

  * Description: Minimum temperature in kelvin expected in the 3D model atmosphere. If not specified, the model will be explored to get it.

MAX_T
-----

  * OPTIONAL

  * Modes: 15D, CLE, INV

  * Formats:
    
    - float

  * Description: Maximum temperature in kelvin expected in the 3D model atmosphere. If not specified, the model will be explored to get it.

MAX_V
-----

  * OPTIONAL

  * Modes: 15D, CLE, INV

  * Formats:
    
    - float

  * Description: Maximum velocity in km s^-1 expected in the 3D model atmosphere. If not specified, the model will be explored to get it.

RT_GROUP_N
----------

  * OPTIONAL

  * Modes: 15D, CLE, INV

  * Formats:
    
    - integer; default: 1

  * Description: Number of processes to solve the forward or inversion problem for each pixel of the input.

UNMAGNETIZED
------------

  * OPTIONAL, ADVANCED

  * Modes: 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Assume no magnetic field in the input model atmosphere.

STATIC
------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Assume no velocity in the input model atmosphere.

STATIC_INT
----------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: STATIC

  * Description: Assume no velocity in the input model atmosphere when solving the only intensity problem.

SKIP_DISK
---------

  * OPTIONAL

  * Modes: CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not calculate the formal solution is the line of sight intersects the stellar disk.

INIT_J_BB
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Include bound-bound transitions in the only radiation initial iterations.

ITER_MIN
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 1

  * Description: Index of the first iteration of a formal solution.

ITERI_MIN
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: ITER_MIN

  * Description: Index of the first iteration of an only intensity formal solution.

ITER_MAX
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 500

  * Description: Maximum iteration index allowed for in a self-consistent solution.

ITERI_MAX
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: ITER_MAX

  * Description: Maximum iteration index allowed for in an only intensity self-consistent solution.

ITERAD_MAX
----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: ITER_MAX

  * Description: Maximum iteration index allowed for in the self-consistent solution for the PRD-AD phase if TWO_STEP_AD is activated.

ITERIAD_MAX
-----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: ITERI_MAX

  * Description: Maximum iteration index allowed for in an only intensity self-consistent solution for the PRD-AD phase if TWO_STEP_AD is activated.

ITER_2ORD
---------

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Account for partial frequency redistribution effects.

TWO_STEP_B
----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If solving polarization with magnetic field, solve first the non-magnetic problem. It is most beneficial when the intensity problem is axial (and has been configured as axial). Otherwise, it is not necessarily faster, but problem dependent.

ITER_MRC_J
----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Maximum relative change of the radiation field to consider it converged in the only radiation field initial iterations.

ITER_MRC_I
----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: 1e-5

  * Description: Maximum relative change of the populations to consider that they have converged.

ITERI_MRC_I
-----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: ITER_MRC_I

  * Description: Maximum relative change of the populations to consider that they have converged in the only intensity solution.

ITER_MRC_P
----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: 1e-3

  * Description: Maximum relative change of the non-population density matrix elements to consider that they have converged.

ITER_MRC_ADI
------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: ITER_MRC_I

  * Description: Maximum relative change of the populations to consider that they have converged in the PRD-AD phase if TWO_STEP_AD is activated.

ITERI_MRC_ADI
-------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: ITERI_MRC_I

  * Description: Maximum relative change of the populations to consider that they have converged in the only intensity solution in the PRD-AD phase if TWO_STEP_AD is activated.

ITER_MRC_ADP
------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: ITER_MRC_P

  * Description: Maximum relative change of the non-population density matrix elements to consider that they have converged in the PRD-AD phase if TWO_STEP_AD is activated.

ITER_J
------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 5

  * Description: Number of preliminar no-line iterations to perform to relax the initial radiation field.

ITERI_PRD
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 4

  * Description: Maximum number of only-radiation iterations to perform when there is partial frequency redistribution in the only intensity problem.

ITERI_MRC_R
-----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: 1e-3

  * Description: Maximum relative change of the mean intensity in only-radiation iterations to consider that it is converged in the intensity problem.

ITER_PRD
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 1

  * Description: Maximum number of only-radiation iterations to perform when there is partial frequency redistribution.

ITER_MRC_R
----------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Maximum relative change of the mean intensity in only-radiation iterations to consider that it is converged.

ITER_MRC_P_R
------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - float; default: 1e-1

  * Description: Maximum relative change of the radiation field tensors (not mean intensity) in only-radiation iterations to consider that it is converged.

NG_ACC
------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Accelerate the convergence with Ng's algorithm.

NG_ORDER
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 3

  * Description: Order of the Ng acceleration.

NG_DELAY
--------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 20

  * Description: Iteration index from which data for Ng acceleration starts to accumulate.

NGI_ACC
-------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: NG_ACC

  * Description: Accelerate the convergence with Ng's algorithm in the only intensity problem.

NGI_ORDER
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: NG_ORDER

  * Description: Order of the Ng acceleration in the only intensity problem.

NGI_DELAY
---------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: NG_DELAY

  * Description: Iteration index from which data for Ng acceleration starts to accumulate in the only intensity problem.

PRD_DELAY
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 1

  * Description: Iteration index from which to start accounting for partial frequency redistribution in the only intensity problem.

ALI_PHOTO
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Compute and use the lambda operator for photoionization transitions.

ALI_DELAY
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - integer; default: 1

  * Description: Iteration index from which to start using the accelerated lambda iteration algorithm in the only intensity problem.

ALI_FORCE
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force more iterations with ALI if converged with delayed ALI iterations.

ALI_ALLOW_OFF
-------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Can switch off ALI if the SEE return negative populations.

APPEND_MRC
----------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Append the maximum relative change of each iteration into the existing file, if present.

APPENDI_MRC
-----------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Append the maximum relative change of each iteration into the existing file, if present, in the only intensity problem.

ALLOW_NPHYS_STK
---------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: -1

  * Description: Iteration index from which finding a non-physical Stokes parameter results in termination of the formal solution. Negative means from the beginning.

ALLOW_NPHYS_RHO
---------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: -1

  * Description: Iteration index from which finding a non-physical density matrix results in termination of the formal solution. Negative means from the beginning.

ALLOW_NPHYS_POP
---------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - integer; default: -1

  * Description: Iteration index from which finding a non-physical populations in the only intensity problem results in termination of the formal solution. Negative means from the beginning.

SOLUTION_BOX
------------

  * OPTIONAL

  * Modes: 15D, INV

  * Formats:
    
    - integer*4; default: -1 -1 -1 -1

  * Description: Indicate the initial x index, final x index, initial y index, and final y index, respectively, of the pixels to solve in a 3D model or data file. Negative numbers are wildcards (automatically adjusted to the relevant size of the input).

EXCLUDE_PIXEL
-------------

  * OPTIONAL, ADDITIVE

  * Modes: 15D, INV

  * Formats:
    
    - integer*2

  * Description: Pair of X and Y pixel coordinate for a pixel that must be excluded from the calculations. Note that, for the inversion model, this is not equivalent to setting a mask with INV_MASK, as the pixel will be completely skipped and nothing will be written.

STORE_STEP
----------

  * OPTIONAL

  * Modes: 1D

  * Formats:
    
    - integer; default: -1

  * Description: If positive, a solution file will be created when the iteration index in the polarized problem is a multiple of this input.

STOREI_STEP
-----------

  * OPTIONAL

  * Modes: 1D

  * Formats:
    
    - integer; default: -1

  * Description: If positive, a solution file will be created when the iteration index in the only intensity problem is a multiple of this input.

CONTRIBUTION
------------

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Calculate and save the contribution function for the last formal solutions.

TAU1
----

  * OPTIONAL

  * Modes: 1D, 15D, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Calculate and save the height where the optical depth is equal to one at each frequency for the last formal solutions.

KEEP_BACK
---------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the calculated background opacity, scattering coefficient, and emissivity.

KEEP_DAMP
---------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the damping parameter characteristic of the Voigt profile for the ATOM_FILE model atoms.

KEEP_QEL
--------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the elastic rates characteristic of each transition for the ATOM_FILE model atoms.

KEEP_APARAM
-----------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Create an ASCII file with the parameters that would be necessary to specify in the model atom at each height to get the same collisional broadening as with the current option in ATOM_FILE

KEEP_COLS
---------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the rates of inelastic collisions between every term and level in the ATOM_FILE model atoms.

KEEP_ATMO
---------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Create a file with all the atmospheric quantities (original or derived) of the atmospheric model.

KEEP_POP
--------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the populations of the ATOM_FILE model atoms.

KEEP_DEP
--------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the departure coefficients of the ATOM_FILE model atoms.

KEEP_RHOKQ
----------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the density matrix elements of the ATOM_FILE model atoms.

KEEP_JKQ
--------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the line integrated radiation field tensors of the ATOM_FILE model atoms.

KEEP_STOKES_QUAD
----------------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: Yes (only 1D), No (only 15D and CLE)

  * Description: Save in a file the Stokes parameters in the model atomsphere for the quadrature directions.

KEEP_JKQNU
----------

  * OPTIONAL

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save in a file the frequency dependent radiation field tensors.

KEEP_MRC
--------

  * OPTIONAL

  * Modes: 1D, 15D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save the maximum relative change in a file.

KEEP_COL_LOG
------------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Save the log of existing and missing collisional rates in the active atoms.

KEEP_MPI_LOG
------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Keep a log with the distribution of tasks among CPU.

KEEP_MPI_DETAIL_LOG
-------------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Keep a more detailed log (written in binary) with the distribution of tasks and weights among CPU.

LIM_STK
-------

  * OPTIONAL, ADDITIVE

  * Modes: 15D, CLE

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for Stokes parameters. Frequencies not included in the specified ranges will not be included in the output.

LIM_CTR
-------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for the contribution function. Frequencies not included in the specified ranges will not be included in the output.

LIM_TAU
-------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for the height where the optical depth is equal to one. Frequencies not included in the specified ranges will not be included in the output.

LIM_COLS_TT
-----------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - string + integer*2

  * Description: Unique identifier of the model atom and doublet of indexes specifying the two terms whose collisional rates between them will be included in the output for collisional rates. Transitions not included in the specified ranges will not be included in the output.

LIM_COLS_LL
-----------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - string + integer*2

  * Description: Unique identifier of the model atom and doublet of indexes specifying the two levels whose collisional rates between them will be included in the output for collisional rates. Transitions not included in the specified ranges will not be included in the output.

LIM_DAMP
--------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - string + integer

  * Description: Unique identifier of the model atom and index specifying the bound-bound transition whose damping parameter will be included in the output for damping parameters. Transitions not included in the specified ranges will not be included in the output.

LIM_QEL
-------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - string + integer

  * Description: Unique identifier of the model atom and index specifying the bound-bound transition whose elastic rate will be included in the output for elastic rates. Transitions not included in the specified ranges will not be included in the output.

LIM_BACK
--------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for the background opacity quantities. Frequencies not included in the specified ranges will not be included in the output.

LIM_POP
-------

  * OPTIONAL, ADDITIVE

  * Modes: 15D

  * Formats:
    
    - string + integer

  * Description: Unique identifier of the model atom and index of the level whose population will be included in the output for populations and departure coefficients. Levels not included in the specified ranges will not be included in the output.

REDO_NE
-------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Init, fin, both, No; No

  * Description: Recalculate the electron density before iterating (init), after iterating (fin), both (both), or not at all (no).

UPDATE_ATMOS
------------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, standard, ne, pe, rhoe, pg, rho, no; No

  * Description: Create a file with the updated atmospheric model with the same main variable than the input (yes), with electron and Hydrogen number densities (standard), with electron number density (ne), with electron pressure (Pe), with electron density (rhoe), with gas pressure (pg), or with total density (rho).

PROTECT_H
---------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic Hydrogen number density.

PROTECT_HM
----------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic Hydrogen minus number density.

CHEM_PROTECT_ALL
----------------

  * OPTIONAL, ADVANCED

  * Modes: 1D, 15D, CLE

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic number density of any atom.

WRITE_PERFORMANCE
-----------------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Create a file with the timing for some stops within the source.

WRITE_MPI_PERFORMANCE
---------------------

  * OPTIONAL, ADVANCED

  * Modes: 1D

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Create a file with the timing for some stops within the source, deeper and more verbose than WRITE_PERFORMANCE.

VERBOSE
-------

  * OPTIONAL

  * Modes: 1D, 15D, CLE, INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Indicate if the verbosity output of the code must go to the command line (Yes) or to a verbosity file to be created in the output directory (No).

VERBOSE_INV_LV
--------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer: 0, 1, 2, 3; default: 0

  * Description: Level of verbosity, decides what messages go into which file.

VERBOSE_INV_SHUTUP
------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer: 0, 1, 2, 3; default: 3

  * Description: Largest level of verbosity in the inversion that will be outputted at all.

TYPE_INVERSION
--------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Thermal, Magnetic, All, Sequential, Sequential-Magnetic, Sequential-Full; default: Thermal

  * Description: Type of inversion, i.e., without magnetic field (thermal), only the magnetic field (magnetic), both (all), both but with a first convergence without magnetic field (sequential), first everything but magnetic field and then only the magnetic field (sequential-magnetic), and the same with a last cycle with everything (sequential-full). Note: only the first word will be registered when reading the input, so do not introduce spaces.

AUTO_WEIGHT
-----------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Let the code choose the weights automatically (NOT RECOMMENDED).

WEIGHT
------

  * MANDATORY

  * Modes: INV

  * Formats:
    
    - float*4
    - float*6

  * Description: Weights for Stokes I, Q, U, and V in the merit function (only one entry). Ignored if AUTO_WEIGHT=Yes. Not mandatory if WEIGHT_FILE is specified. | Wavelength range in nm for which these weights must be used, and weights for Stokes I, Q, U, and V in the merit function. Ignored if AUTO_WEIGHT=Yes. Not mandatory if WEIGHT_FILE is specified. Note that the wavelength range does not truncate, but look for the closest wavelength in the data. Moreover, they are intended to be for distinguishing between different spectral ranges; in order to customize weights within a spectral line, use a file to specify the weights.

WEIGHT_FILE
-----------

  * MANDATORY

  * Modes: INV

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of a file with the wavelengths ranges and weights for Stokes I, Q, U, and V to be used in the inversion. Ignored if AUTO_WEIGHT=Yes. Not mandatory if WEIGHT is specified.

WEIGHT_FACTOR
-------------

  * OPTIONAL, ADDITIVE

  * Modes: INV

  * Formats:
    
    - string + float*3

  * Description: Each entry of this keyword specifies a multiplicative factor for the weights of a given Stokes parameter between two wavelengths. The string must be I, Q, U, or V, the first two floats specify the wavelength range of the weights to enhance in nanometers, and the last one is a non-negative multiplicative factor. If there are several ranges sharing at least one wavelength for the same Stokes parameter, the inversion will be aborted before starting. Ignored if AUTO_WEIGHT=Yes.

SIGMA_FACTOR
------------

  * OPTIONAL, ADDITIVE

  * Modes: INV

  * Formats:
    
    - string + float*3

  * Description: Each entry of this keyword specifies a multiplicative factor for the sigma of a given Stokes parameter between two wavelengths. The string must be I, Q, U, or V, the first two floats specify the wavelength range of the weights to enhance in nanometers, and the last one is a non-negative multiplicative factor. If there are several ranges sharing at least one wavelength for the same Stokes parameter, the inversion will be aborted before starting. Ignored if there is no sigma information in the data.

INV_INIT
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Init, file path

  * Description: Path, absolute or relative to the running directory, of a file with the result of a previous inversion. Init means starting from scratch.

INV_MASK
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, file path

  * Description: Path, absolute or relative to the running directory, of a file with a mask to determine pixels for which only the error will be calculated. This is only used when restarting from a previous file. The error will not be fully consistent with the model, as while it is kept unchanged, internally it is interpolated with splines from a (in general) new set of nodes.

CENTERED_DERIVATIVE
-------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: If the derivative to compute the numerical response function must be centered, i.e., computing the perturbation with both signs.

ITER_MAX_INV
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 15

  * Description: Maximum iteration index allowed in the inversion.

NODES_B_METHOD
--------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the magnetic field strength or longitudinal component (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_BT_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the magnetic field inclination or transversal component (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_BP_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the magnetic field azimuth (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_T_METHOD
--------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the temperature field azimuth (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VX_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the velocity x component or inclination (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VY_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the velocity y component or azimuth (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VZ_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the velocity vertical or line of sight component (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_VT_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the micro-turbulent velocity (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_PG_METHOD
---------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the gas pressure (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J21R_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the real part of the J21 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J21I_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the imaginary part of the J21 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J22R_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the real part of the J22 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

NODES_J22I_METHOD
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: value, value fix first, value fix last, value fix extremes, correction, correction fix first, correction fix last, correction fix extremes; default: value

  * Description: What should the nodes represent for the imaginary part of the J22 component of the radiation field tensors (value or correction) and if one or both extremes of the nodes in the atmosphere must remain fixed.

INTERPOLATION
-------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: linear, quadratic bezier, cubir bezier; default: cubic bezier

  * Description: Type of interpolation to generate the atmospheric model stratification from the node values.

NODES_B_NUM
-----------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the magnetic field strength or longitudinal component. Mutually exclusive with NODES_B_LOCATION.

NODES_BT_NUM
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the magnetic field inclination or transversal component. Mutually exclusive with NODES_BT_LOCATION.

NODES_BP_NUM
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the magnetic field azimuth. Mutually exclusive with NODES_BP_LOCATION.

NODES_F_NUM
-----------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the diffuse light factor. Mutually exclusive with NODES_F_LOCATION.

NODES_T_NUM
-----------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the temperature. Mutually exclusive with NODES_T_LOCATION.

NODES_VX_NUM
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the velocity x component or transversal component. Mutually exclusive with NODES_VX_LOCATION.

NODES_VY_NUM
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the velocity y component or azimuth. Mutually exclusive with NODES_VY_LOCATION.

NODES_VZ_NUM
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the velocity vertical or longitudinal component. Mutually exclusive with NODES_VZ_LOCATION.

NODES_PG_NUM
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the gas pressure. Mutually exclusive with NODES_PG_LOCATION.

NODES_J21R_NUM
--------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the real part of the J21 radiation field tensor. Mutually exclusive with NODES_J21R_LOCATION.

NODES_J21I_NUM
--------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the imaginary part of the J21 radiation field tensor. Mutually exclusive with NODES_J21I_LOCATION.

NODES_J22R_NUM
--------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the real part of the J22 radiation field tensor. Mutually exclusive with NODES_J22R_LOCATION.

NODES_J22I_NUM
--------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes for the imaginary part of the J22 radiation field tensor. Mutually exclusive with NODES_J22I_LOCATION.

NODES_B_LOCATION
----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the magnetic field strength or longitudinal component. Mutually exclusive with NODES_B_NUM.

NODES_BT_LOCATION
-----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the magnetic field inclination or transversal component. Mutually exclusive with NODES_BT_NUM.

NODES_BP_LOCATION
-----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the magnetic field azimuth. Mutually exclusive with NODES_BP_NUM.

NODES_F_LOCATION
----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the diffuse light factor. Mutually exclusive with NODES_F_NUM.

NODES_T_LOCATION
----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the temperature. Mutually exclusive with NODES_T_NUM.

NODES_VX_LOCATION
-----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the velocity x component or transversal component. Mutually exclusive with NODES_VX_NUM.

NODES_VY_LOCATION
-----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the velocity y component or azimuth. Mutually exclusive with NODES_VY_NUM.

NODES_VZ_LOCATION
-----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the velocity vertical or longitudinal component. Mutually exclusive with NODES_VZ_NUM.

NODES_PG_LOCATION
-----------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the gas pressure. Mutually exclusive with NODES_PG_NUM.

NODES_J21R_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the real part of the J21 radiation field tensor. Mutually exclusive with NODES_J21R_NUM.

NODES_J21I_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the imaginary part of the J21 radiation field tensor. Mutually exclusive with NODES_J21I_NUM.

NODES_J22R_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the real part of the J22 radiation field tensor. Mutually exclusive with NODES_J22R_NUM.

NODES_J22I_LOCATION
-------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - floats

  * Description: Positions of the nodes, in the decimal logarithm of the optical depth, for the imaginary part of the J22 radiation field tensor. Mutually exclusive with NODES_J22I_NUM.

NODES_B_EXTRAPOLATION
---------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the magnetic field strength toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_BT_EXTRAPOLATION
----------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the magnetic field inclination or transversal component toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_BP_EXTRAPOLATION
----------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the magnetic field azimuth toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_T_EXTRAPOLATION
---------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the temperature toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_VX_EXTRAPOLATION
----------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the velocity x component or transversal component toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_VY_EXTRAPOLATION
----------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the velocity y component or azimuth toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_VZ_EXTRAPOLATION
----------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the velocity vertical or longitudinal component toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_PG_EXTRAPOLATION
----------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the gas pressure toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_J21R_EXTRAPOLATION
------------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the real part of the J21 radiation field tensor toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_J21I_EXTRAPOLATION
------------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the imaginary part of the J21 randiation field tensor toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_J22R_EXTRAPOLATION
------------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the real part of the J22 radiation field tensor toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

NODES_J22I_EXTRAPOLATION
------------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string*2: none, zero, constant, linear; default: none none

  * Description: Type of extrapolation for the imaginary part of the J22 randiation field tensor toward the top and the bottom, respectively, of the model atmosphere when the extrema nodes are not at the extremes of the model. The first string defined the top and the second the bottom. If only one string is specified, both are consider equal. "none" means no extrapolation (kept the same than in the initial model), "zero" means to always set the values to zero, "constant" means to extent the values of the extrema nodes, and linear means linear interpolation.

BTYPE
-----

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Vertical, LOS; default: Vertical

  * Description: Invert the magnetic field strength, inclinarion, and azimuth in the vertical reference frame (vertical), or the longitudinal and transversal component, and the azimuth in the line of sight reference frame (LOS).

VTYPE
-----

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Vertical, LOS; default: Vertical

  * Description: Invert the velocity x, y, and vertical components in the vertical reference frame (vertical), or the longitudinal and transversal component, and the azimuth in the line of sight reference frame (LOS).

FIX_B
-----

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the magnetic field strength or longitudinal component.

FIX_BT
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the magnetic field inclination or transversal component.

FIX_BP
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the magnetic field azimuth.

FIX_F
-----

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the difusse light factor.

FIX_T
-----

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the temperature.

FIX_VX
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the velocity x or transversal component.

FIX_VY
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the velocity y component or azimuth.

FIX_VZ
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the velocity vertical or longitudinal component.

FIX_VT
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the micro-turbulent velocity.

FIX_PG
------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the gas pressure.

FIX_J21R
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the real part of the J21 radiation field tensor.

FIX_J21I
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the imaginary part of the J21 radiation field tensor.

FIX_J22R
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the real part of the J22 radiation field tensor.

FIX_J22I
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fix the node values of the imaginary part of the J22 radiation field tensor.

POSITION_CORRECTION
-------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Correct the positions of the nodes to coincide with the closest node of the forward solver stratification.

REGUL_B
-------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the magnetic field strength or longitudinal component. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_BT
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the magnetic field inclination or longitudinal component. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_BP
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the magnetic field azimuth. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_F
-------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the diffuse light factor. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_T
-------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the temperature. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_VX
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the velocity x or transversal component. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_VY
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the velocity y component or azimuth. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_VZ
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the velocity vertical or longitudinal component. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_VT
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the micro-turbulent velocity. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_PG
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the gas pressure. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_J21R
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the real part of the J21 radiation field tensor. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_J21I
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the imaginary part of the J21 radiation field tensor. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_J22R
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the real part of the J22 radiation field tensor. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_J22I
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: None, mean, constant, first derivative, second derivative; + float; default: none

  * Description: Type of regulatization and associated weight for the imaginary part of the J22 radiation field tensor. If the type of regularization is constant, a second float can be added to choose the value for the comparison; if not provided, the value will be the average of the initial values of the nodes.

REGUL_LIMITS
------------

  * MANDATORY, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1.0

  * Description: Additional factor to the regularization penalties.

THRH_CHI2
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Threshold on the value of the merit function to stop the inversion.

INV_MRC
-------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Maximum relative change of the merit function to consider it converged.

SVD_TYPE
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: traditional, sir; default: sir

  * Description: Use the traditional singular value decomposition (traditional) or the one with the same modifications implemented in the SIR code (sir).

THRH_SVD
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Threshold for the singular value decomposition.

PG_TYPE
-------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Hydrostatic equilibrium, stratified; default: hydrostatic equilibrium

  * Description: How to deal with the gas pressure stratification. The stratified option is NOT RECOMMENDED.

PG_BOUND
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float; default: -1

  * Description: Initial or fixed value of the gas pressure at the top boundary, in dyn cm^-2, if PG_TYPE = hydrostatic. A negative value means that the gas pressure at the topmost node of the initial model atmosphere will be taken instead.

DIFFUSE_LIGHT
-------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float; default: -1

  * Description: Initial or fixed value of the diffuse light factor.

ATMO_NODES
----------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Number of nodes in the model atmosphere in the forward solution. A 0 value indicates that the stratification of the initial model atmosphere will be kept.

MAX_SVD_STEP
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    ; float (1.0)

  * Description: Maximum step allowed in the singular value decomposition.

INV_ERROR
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Hessian, RF, worst, recycle; default: Hessian

  * Description: Algorithm to estimate the error in the inversion. Recycle uses the last Hessian, Hessian recalculates the Hessian, RF is slightly different, using directly the response functions, and worst returns the worst between Hessian and RF.

PSF_FWHM
--------

  * OPTIONAL, ADDITIVE

  * Modes: INV

  * Formats:
    
    - float
    - float*3
    - float*2 + string: file path

  * Description: The FWHM of a gaussian spectral PSF in nanometers (only one entry). | Wavelength range in nm for which the specified FWHM must be used, and FWHM of a gaussian spectral PSF in nanometers. | Wavelength range in nm for which the specified PSF must be used and path, absolute or relative to the running directory, of a file with the spectral PSF. Note that the wavelength range does not truncate, but look for the closest wavelength in the data. Moreover, they are intended to be for distinguishing between different spectral ranges.

BOUNDS_B
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: 0 3000 (only BTYPE=vertical); default: -3000 3000 (only BTYPE=los)

  * Description: Minimum and maximum values that the magnetic field strengths or longitudinal component can take in the inversion.

BOUNDS_BT
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: 0 pi (only BTYPE=vertical); default: 0 3000 (only BTYPE=los)

  * Description: Minimum and maximum values that the magnetic field inclination or transversal component can take in the inversion.

BOUNDS_BP
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: 0 2pi

  * Description: Minimum and maximum values that the magnetic field azimuth can take in the inversion.

BOUNDS_F
--------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: 0 0.95

  * Description: Minimum and maximum values that the diffuse light factor can take in the inversion.

BOUNDS_VX
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: -20 20 (only vtype=vertical); default: 0 20 (only vtype=los)

  * Description: Minimum and maximum values that the velocity x or transversal component can take in the inversion.

BOUNDS_VY
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: -20 20 (only vtype=vertical); default: 0 2pi (only vtype=los)

  * Description: Minimum and maximum values that the velocity y component or azimuth can take in the inversion.

BOUNDS_VZ
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: -20 20

  * Description: Minimum and maximum values that the velocity vertical or longitudinal component can take in the inversion.

BOUNDS_VT
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: 0 40

  * Description: Minimum and maximum values that the micro-turbulent velocity can take in the inversion.

BOUNDS_PG
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: 0.1 15

  * Description: Minimum and maximum values that the gas pressure can take in the inversion.

BOUNDS_J21R
-----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the real part of the J21 radiation field tensor can take in the inversion.

BOUNDS_J21I
-----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the imaginary part of the J21 radiation field tensor can take in the inversion.

BOUNDS_J22R
-----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the real part of the J22 radiation field tensor can take in the inversion.

BOUNDS_J22I
-----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*2; default: -1 1

  * Description: Minimum and maximum values that the imaginary part of the J22 radiation field tensor can take in the inversion.

EBOUNDS_B
---------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the magnetic field strengths or longitudinal component can take in the inversion.

EBOUNDS_BT
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the magnetic field inclination or transversal component can take in the inversion.

EBOUNDS_BP
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the magnetic field azimuth can take in the inversion.

EBOUNDS_VX
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the velocity x or transversal component can take in the inversion.

EBOUNDS_VY
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the velocity y component or azimuth can take in the inversion.

EBOUNDS_VZ
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the velocity vertical or longitudinal component can take in the inversion.

EBOUNDS_VT
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the micro-turbulent velocity can take in the inversion.

EBOUNDS_PG
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the gas pressure can take in the inversion.

EBOUNDS_J21R
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the real part of the J21 radiation field tensor can take in the inversion.

EBOUNDS_J21I
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the imaginary part of the J21 radiation field tensor can take in the inversion.

EBOUNDS_J22R
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the real part of the J22 radiation field tensor can take in the inversion.

EBOUNDS_J22I
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*4

  * Description: Minimum and maximum in logarithmic optical depth where to force the minimum and maximum specified value that the imaginary part of the J22 radiation field tensor can take in the inversion.

SCALE_B
-------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 100

  * Description: Scale factor for the magnetic field strengths or longitudinal component can take in the inversion.

SCALE_BT
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1.8 (only BTYPE=vertical); default: 100 (only BTYPE=los)

  * Description: Scale factor for the magnetic field inclination or transversal component can take in the inversion.

SCALE_BP
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1.8

  * Description: Scale factor for the magnetic field azimuth can take in the inversion.

SCALE_F
-------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1

  * Description: Scale factor for the diffuse light factor can take in the inversion.

SCALE_VX
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 10

  * Description: Scale factor for the velocity x or transversal component can take in the inversion.

SCALE_VY
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 10 (only vtype=vertical); default: 1.8 (only vtype=los)

  * Description: Scale factor for the velocity y component or azimuth can take in the inversion.

SCALE_VZ
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 10

  * Description: Scale factor for the velocity vertical or longitudinal component can take in the inversion.

SCALE_VT
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 10

  * Description: Scale factor for the micro-turbulent velocity can take in the inversion.

SCALE_PG
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 2

  * Description: Scale factor for the gas pressure can take in the inversion.

SCALE_J21R
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the real part of the J21 radiation field tensor can take in the inversion.

SCALE_J21I
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the imaginary part of the J21 radiation field tensor can take in the inversion.

SCALE_J22R
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the real part of the J22 radiation field tensor can take in the inversion.

SCALE_J22I
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Scale factor for the imaginary part of the J22 radiation field tensor can take in the inversion.

PERTURB_B
---------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1

  * Description: Value of the perturbation to compute the response function for the magnetic field strengths or longitudinal component can take in the inversion.

PERTURB_BT
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2 (only BTYPE=vertical); default: 2 (only BTYPE=los)

  * Description: Value of the perturbation to compute the response function for the magnetic field inclination or transversal component can take in the inversion.

PERTURB_BP
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 3e-2

  * Description: Value of the perturbation to compute the response function for the magnetic field azimuth can take in the inversion.

PERTURB_F
---------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-2

  * Description: Value of the perturbation to compute the response function for the diffuse light factor can take in the inversion.

PERTURB_VX
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 3e-1

  * Description: Value of the perturbation to compute the response function for the velocity x or transversal component can take in the inversion.

PERTURB_VY
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 3e-1 (only vtype=vertical); default: 3e-2 (only vtype=los)

  * Description: Value of the perturbation to compute the response function for the velocity y component or azimuth can take in the inversion.

PERTURB_VZ
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 3e-1

  * Description: Value of the perturbation to compute the response function for the velocity vertical or longitudinal component can take in the inversion.

PERTURB_VT
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 3e-1

  * Description: Value of the perturbation to compute the response function for the micro-turbulent velocity can take in the inversion.

PERTURB_PG
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 5e-2

  * Description: Value of the perturbation to compute the response function for the gas pressure can take in the inversion.

PERTURB_J21R
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the real part of the J21 radiation field tensor can take in the inversion.

PERTURB_J21I
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the imaginary part of the J21 radiation field tensor can take in the inversion.

PERTURB_J22R
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the real part of the J22 radiation field tensor can take in the inversion.

PERTURB_J22I
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1e-4

  * Description: Value of the perturbation to compute the response function for the imaginary part of the J22 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_B
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the magnetic field strengths or longitudinal component can take in the inversion.

MIN_REL_PERTURB_BT
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the magnetic field inclination or transversal component can take in the inversion.

MIN_REL_PERTURB_BP
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the magnetic field azimuth can take in the inversion.

MIN_REL_PERTURB_F
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the diffuse light factor can take in the inversion.

MIN_REL_PERTURB_VX
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the velocity x or transversal component can take in the inversion.

MIN_REL_PERTURB_VY
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the velocity y component or azimuth can take in the inversion.

MIN_REL_PERTURB_VZ
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the velocity vertical or longitudinal component can take in the inversion.

MIN_REL_PERTURB_VT
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the micro-turbulent velocity can take in the inversion.

MIN_REL_PERTURB_PG
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the gas pressure can take in the inversion.

MIN_REL_PERTURB_J21R
--------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the real part of the J21 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_J21I
--------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the imaginary part of the J21 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_J22R
--------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the real part of the J22 radiation field tensor can take in the inversion.

MIN_REL_PERTURB_J22I
--------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.02

  * Description: Minimum relative perturbation to compute the response function for the imaginary part of the J22 radiation field tensor can take in the inversion.

INI_BPOS
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.17 (only BTYPE=vertical); default: 10 (only BTYPE=los)

  * Description: Initial magnetic field inclination or transversal component when initializing from an inversion result with a very small value of this quantity.

INI_BAZI
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.5

  * Description: Initial magnetic field azimuth when initializing from an inversion result with a very small value of this quantity.

INI_VPOS
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 5

  * Description: Initial velocity x or transversal component when initializing form an inversion result with a very small value of this quantity.

INI_VAZI
--------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 1.5

  * Description: Initial velocity y or azimuth when initializing form an inversion result with a very small value of this quantity.

GUESS_POLARITY
--------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2 or float*4; default: nothing

  * Description: Try to estimate the polarity (if two floats) and the polarity and strength of the longitudinal magnetic field component (if four floats). The first two floats indicate the wavelength range to consider for the estimation, in nanometers, followed by the effective Landé factor and line's wavelength to consider.

INV_FRACTION
------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Fit the fractional polarization (Yes) instead of the Stokes profiles (No).

INV_TAU_RANG
------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - float*2; default: -8 1

  * Description: Logarithmic optical depth of the extreme nodes of the model atmosphere in the forward solution. This has no effect if ATMO_NODES is not specified.

ATMO_STRAT
----------

  * OPTIONAL, ADDITIVE

  * Modes: INV

  * Formats:
    
    - float*3

  * Description: Triplets of number, with the first two specifying a range in logarithmic optical depth in which the sampling step of the atmospheric model for the forward solver must be reduced by a factor given by the third number.

BROYDEN_LM
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Use Broyden's algorithm in the profile fitting.

LM_METHOD
---------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Traditional, backtracking; default: backtracking

  * Description: Type of Levenberg-Marquardt method for the inversion minimization.

LM_BACKTRACKING_MODE
--------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: desperation or nothing; default: nothing

  * Description: If "desperation" mode is activated, if the backtracking method (see LM_METHOD) gets stuck, it will inspect if the lambda values in the trials pertain to a certain regime and will force more trials in a different one.

LM_LAM_BIG_TEST
---------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.051

  * Description: Value that the minimum lambda in trials must surpass to be considered in the big regime in the backtracking when LM_BACKTRACKING_MODE is desperate. Must be within the limits specified in LM_LAMBDA_RANG.

LM_LAM_SMALL_TEST
-----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 9.9

  * Description: Value that the maximum lambda in trials must not surpass to be considered in the small regime in the backtracking when LM_BACKTRACKING_MODE is desperate. Must be within the limits specified in LM_LAMBDA_RANG.

LM_LAM_BIG_PROVE
----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 100.

  * Description: Value for the lambda in the next trial when forcing the big lambda regime in the backtracking when LM_BACKTRACKING_MODE is desperate.

LM_LAM_SMALL_PROVE
------------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 0.1

  * Description: Value for the lambda in the next trial when forcing the small lambda regime in the backtracking when LM_BACKTRACKING_MODE is desperate.

LM_LAMBDA_TRACK
---------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - integer; default: 0

  * Description: Specifies the amount of lambda values stored and used to choose the initial lambda for the current iteration. 0 means that no data is stored and the lambda value starts the same always. When 1, the initial lambda is the best of the previous iteration. When 2, the new lambda is linearly extrapolated. When 3, the new lambda is extrapolated with cubic splines. When extrapolating, the limits in LM_LAMBDA_RANG are respected. When 3, if the derivative changes sign the new value cannot overshoot the existing values.

LM_LAMBDA_RANG
--------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float*2; default: 1e-5 1e3

  * Description: Boundary limits for the lambda coefficient in the Levenmberg-Marquardt algorithm.

LM_LAMBDA_ACCEPT
----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 5

  * Description: Factor dividing the lambda coefficient in the Levenberg-Marquardt algorithm when the previous lambda value was accepted.

ALLOW_REDUCED_MODE
------------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes/No; default: No

  * Description: Allow a last iteration in which the weights are changed to try to focus on the worst fit regions.

LM_LAMBDA_REJECT
----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - float; default: 5

  * Description: Factor multiplying the lambda coefficient in the Levenberg-Marquardt algorithm when the previous lambda value was rejected.

INV_B_PROJECTION
----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Make sure to correctly project the negative polarity when the magnetic field is in the vertical reference frame but you want to use it as its longitudinal component.

RF_INITSOL
----------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Initialize the calculation of the response function with the solution of the reference profile.

TRIAL_INITSOL
-------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: Yes

  * Description: Allow trials to start from the last solution if the thermal parameters are fixed.

TRIAL_INITSOL_TP
----------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Allow trials to start from the last solution if temperature and gas pressure are fixed. Requires TRIAL_INITSOL = Yes to activate.

INV_NEGL_SIGMA
--------------

  * OPTIONAL, ADVANCED

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Neglect the sigma (noise) in the observation for the merit function.

KEEP_RF
-------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Include the numerical response functions in the solution file.

STOREINV_STEP
-------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - integer; default: -1

  * Description: If positive, the inversion solution will be stored in the file when the iteration index in the polarized problem is a multiple of this input.

FORCE_OBS_FREQ
--------------

  * OPTIONAL

  * Modes: INV

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Force the wavelengths of the data to be present in the frequency axis of the forward solver.

