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

ATOM_ION
--------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path; optional: string (a single word and unique label) (default: index given by reading order)

  * Description: Path, absolute or relative to the running directory, of the file with the ionization fraction for the model atom specified by the label.

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

SPECT_INPUT
-----------

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the spectra for the incoming radiation.

CHIANTI_PATH
------------

  * OPTIONAL

  * Formats:
    
    - string: folder path

  * Description: Path, absolute or relative to the running directory, to the folder with the CHIANTI database data.

T_RAD
-----

  * OPTIONAL

  * Formats:
    
    - float; default: 5e3

  * Description: Effective temperature  in kelvin of the black-body radiation for the boundary condition if SPECT_INPUT is not specified or for frequencies not included in such file.

R_STAR
------

  * OPTIONAL

  * Formats:
    
    - float; default: 6.957e10

  * Description: Stellar radius in cm.

ABUND
-----

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with atomic abundances.

KURUCZ
------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, to the file with a list of transitions in the Kurucz's database format to be included in the background opacity.

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

FORCE
-----

  * OPTIONAL

  * Formats:
    
    - string: intensity, polarization, all, none; default: none

  * Description: Force the code to solve only the intensity problem (intensity), to solve the polarization problem skipping the intensity problem (polarization), to solve both intensity and polarization problem (all), or let the code decide from the available inputs (none).

MIT_OFF
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No, none; default: none

  * Description: Neglect (Yes) or force the inclusion (No) of frequencies corresponding to magnetically induced transitions in multi-term model atoms. By default, they are only included if there is a magnetic field.

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

RAM_LIMIT
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum amount of Megabytes which can be allocated in the form of Voigt profiles, photoionization pre-calculated quantities, interpolation precalculated quantities (partial frequency redistribution), or redistribution functions. Negative means no limit (NOT RECOMMENDED for complex problems).

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

SKIP_DISK
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not calculate the formal solution is the line of sight intersects the stellar disk.

ITER_2ORD
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Account for partial frequency redistribution effects.

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

LIM_STK
-------

  * OPTIONAL, ADDITIVE

  * Formats:
    
    - float*2

  * Description: Wavelengths in nanometers to delimit a spectral range to include in the output file for Stokes parameters. Frequencies not included in the specified ranges will not be included in the output.

VERBOSE
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Indicate if the verbosity output of the code must go to the command line (Yes) or to a verbosity file to be created in the output directory (No).

