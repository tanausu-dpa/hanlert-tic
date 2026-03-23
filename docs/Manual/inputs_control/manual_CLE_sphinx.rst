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
    
    - string: file path, HCC (only INV), HCP (only INV), HCA (only INV), HCF (only INV), HCX (only INV); default: HCC (only INV)

  * Description: Path, absolute or relative to the running directory, of the model atmosphere file. In inversion mode this specifies the initial model atmosphere and HCC, HCP, HCA, HCF, and HCX can be specified to indicate the C, P, A, or F models of Fontenla et al. (1993) or the MCO (a.k.a. FALX) model of Ayres (1986), respectively.

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

SPECT_INPUT
-----------

  * OPTIONAL

  * Formats:
    
    - string: file path

  * Description: Path, absolute or relative to the running directory, of the file with the spectra for the incoming radiation.

USE_ALLEN
---------

  * OPTIONAL

  * Formats:
    
    - string: Yes/No

  * Description: Use Allen's tabulation to determine the intensity of the star.

CLV
---

  * OPTIONAL

  * Formats:
    
    - string: None/Allen

  * Description: Uses Allen's tabulation to determine the intensity CLV or neglect CLV. This is only used if there is no SPECT_INPUT data for a needed frequency.

FLAT_CLE_IN
-----------

  * OPTIONAL

  * Formats:
    
    - string: Yes/No

  * Description: Assume that the radiation is exactly flat when computing the JKQ (only when no input spectra is provided).

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

NEGLECT_CONTINUUM
-----------------

  * OPTIONAL

  * Formats:
    
    - string Yes, No; default: No

  * Description: Neglect continuum contribution to the radiative transfer coefficients.

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

VOI_TYPE
--------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: H82, A67, HAW78, GAUSS; default: H82

  * Description: Algorithm to compute the Voigt profiles. If solving for polarization, only H82 is allowed.

RAM_LIMIT
---------

  * OPTIONAL

  * Formats:
    
    - integer; default: -1

  * Description: Maximum amount of Megabytes which can be allocated per CPU, used to limit the amount of Voigt and redistribution profiles to store in RAM. The counting is not perfect, so be conservative if there is a RAM limit.

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

RED_INT_MODE
------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Linear, Cubic Hermite, Splines; default: Splines

  * Description: Type of interpolation of the second order emissivity when transforming from the comoving to the observer's reference frame.

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

  * Description: Distance in Doppler width from a frequency to the closest transition resonance in order to neglect partial frequency redistribution at that frequency. NOTE: The Doppler width in this variable is the one used to define the output frequency axis, and not the one corresponding to the atmospheric properties at each location.

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

MPI_CV_W_MINIMAL
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 4

  * Description: Minimum number of simultaneous communications to hold in comoving2ord.

MPI_CV_W_NOMINAL
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: 8

  * Description: Expected number of simultaneous communications to hold in comoving2ord.

MPI_CV_W_MEM_LIMIT
------------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 50e6

  * Description: Limit of memory to keep in a shared MPI pool in comoving2ord.

DOP_WIDTH
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - float; default: 2.5e3
    - string: Max, Min

  * Description: Doppler width in m s^-1 to consider when building the frequency axis to transform between Doppler widths and actual frequencies. | Take the maximum (Max) or minimum (Min) temperature to calculate the Doppler width used to build the frequency axis to transform between Doppler widths and actual frequencies.

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

ALLOW_NPHYS_POP
---------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - integer; default: -1

  * Description: Iteration index from which finding a non-physical populations in the only intensity problem results in termination of the formal solution. Negative means from the beginning.

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

PROTECT_H
---------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic Hydrogen number density.

PROTECT_HM
----------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic Hydrogen minus number density.

CHEM_PROTECT_ALL
----------------

  * OPTIONAL, ADVANCED

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Do not let the chemical equilibrium to change the atomic number density of any atom.

VERBOSE
-------

  * OPTIONAL

  * Formats:
    
    - string: Yes, No; default: No

  * Description: Indicate if the verbosity output of the code must go to the command line (Yes) or to a verbosity file to be created in the output directory (No).

