.. _outputs:

*******
Outputs
*******

HanleRT-TIC always generates output files (upon request in the
:ref:`control file <inputs_control>`) that can contain data
on the solution of the radiative transfer problem (e.g., density matrix elements),
the calculation of emergent Stokes profiles, or information about the execution
(e.g., how the tasks were divided between the different processes).

In the Tools directory (see :ref:`folder tree <folder_tree>`) we provide a python
class, which can be accessed by importing *hanlertio.py*::

    import hanlertio
    IO = hanlertio.hanlertio_class(<path_to_output_file(string)>)
    print(IO.help()) # This will show the list of available methods for that file
    print(IO.help(<method(string)>)) # This will show the description of the method and the possible inputs


This class can read (or will be able to read) all the relevant binary output
files, providing methods to conveniently extract the desired information. Whenever
an output file can be read with this class, we will not describe the details of
the format.

.. _cache_file:

Cache file
==========

Although it is not an actual output, as it does not contain any information
regarding the physics of the problem, knowing what this file is and what
is for is critical for most of the applications of the code.

Whevener a run is started in *1.5D Synthesis*, *CLE*, or *Inversion* mode,
a *cache* file is created in the *Output directory*. This file registers
which pixels/columns in the model or the data have been successfully
computed. If a calculation has to be halted, a new call with the same
parameters will use this cache file and will skip everything that was
already successfully finished. Therefore, it is important that this file
is **removed** to start a calculation anew.

.. warning::
   The *cache* file is not size-flexible. That is, changing the input model
   or the size of the output (e.g., with the **SOLUTION_BOX** keyword
   in the :ref:`control file <input_file>`), requires to start anew and
   thus remove the *cache* file. Not doing so will result in an error
   in the best cases, and in unknown consequences in the worst cases.

.. _solution_file:

Solution File
=============

The **Solution** file contains the state of the radiative transfer problem
at the last performed iteration (it can be an intermediate step if
**STOREI_STEP** or **STORE_STEP** are specified in the
:ref:`control file <inputs_control>`). This file contains all the necessary
quantities to either restart the problem from the last performed iteration
or to calculate the emergent Stokes profiles corresponding to the physical
state represented by the file. The **SolutionI** file is a **Solution**
file with the only-intensity solution which is created if the option
**SOLUTION_KEEPI** is *Yes* in the :ref:`control file <inputs_control>`.
Both files are written in binary.

.. include:: solutioni_file.rst
.. include:: solution_file.rst

.. _result_file:

Result File
===========

The **Result** file contains the solution to the inversion problem.

.. include:: result_file.rst

.. _population_file:

Population
==========

The **Population** file can be generated for each active atom present
in the calculation (see **ATOM_INPUT** in
:ref:`control file <inputs_control>`) and gives the number densities of the
excitation balance of the atom.

.. include:: popu_file.rst


.. _departure_file:

Departure
=========

The **Departure** file is generated for each active atom present
in the calculation (see **ATOM_INPUT** in
:ref:`control file <inputs_control>`) and gives the departure coefficients
of the exceitation balance from the local thermodynamical equilibrium.

.. include:: depar_file.rst


.. _jout_file:

Jout
====

The **Jout** file contains the radiation transfer tensor at the
last performed iteration (it can be an intermediate step if
**STOREI_STEP** or **STORE_STEP** are specified in the
:ref:`control file <inputs_control>`). The **JoutI** file is a **Jout**
file which is created if the option **SOLUTION_KEEPI** is *Yes*
in the :ref:`control file <inputs_control>`.

.. include:: jout_file.rst


.. _rhoout_file:

Rhoout
======

The **Rhoout** file contains the density matrix tensor at the
last performed iteration (it can be an intermediate step if
**STOREI_STEP** or **STORE_STEP** are specified in the
:ref:`control file <inputs_control>`). The **RhooutI** file is a **Rhoout**
file which is created if the option **SOLUTION_KEEPI** is *Yes*
in the :ref:`control file <inputs_control>`.

.. include:: rhoout_file.rst


.. _stokesout_file:

Stokesout
=========

The **Stokesout** file contains the Stokes parameters for
the directions in the quadrature which are directed towards
the top boundary for the last performed iteration (it can be an
intermediate step if **STOREI_STEP** or **STORE_STEP** are
specified in the :ref:`control file <inputs_control>`). The
**StokesoutI** file is a **Stokesout** file which is created
if the option **SOLUTION_KEEPI** is *Yes* in the
:ref:`control file <inputs_control>`.

.. include:: stokesout_file.rst


.. _jkqnu_file:

JKQnu
=====

.. include:: jkqnu_file.rst


.. _stokesi_file:

StokesI
=======

The **StokesI** file contains the emergent Stokes intensity for
the specified lines of sight (see **POLAR_LOS** and **AXIAL_LOS**
keywords in :ref:`control file <inputs_control>`), resulting from
the only intensity problem. The name of the output files has
the format **StokesI_X_Y**, with **X** and **Y** the indexes of the
polar and azimuthal angles of the line of sight in the order
they are specified, looping over both of them, with the polar
angle the slower loop.

.. include:: stokesi_file.rst


.. _stokes_file:

Stokes
======

The **Stokes** file contains the emergent Stokes parameters for
the specified lines of sight (see **POLAR_LOS** and **AXIAL_LOS**
keywords in :ref:`control file <inputs_control>`), resulting from
the general polarization problem. The name of the files has
the format **Stokes_X_Y**, with **X** and **Y** the indexes of the
polar and azimuthal angles of the line of sight in the order
they are specified, looping over both of them, with the polar
angle the slower loop. In inversion, if requested, it is named
just **Stokes**.

.. include:: stokes_file.rst


.. _contribution_file:

Contribution
============

The **Contribution** file contains the contribution function
for the four Stokes parameters for the specified lines of
sight (see **POLAR_LOS** and **AXIAL_LOS** keywords in
:ref:`control file <inputs_control>`) if the keyword
**CONTRIBUTION** is set to *Yes* in the
:ref:`control file <inputs_control>`.
The name of the files has the format **Contribution_X_Y**,
with **X** and **Y** the indexes of the polar and azimuthal
angles of the line of sight in the order they are specified,
looping over both of them, with the polar
angle the slower loop. In inversion, if requested, it is named
just **Contribution**.

.. include:: contribution_file.rst


.. _tau_file:

Tau
===

The **Tau** file contains the height where the optical depth
is equal to one for every frequency of the problem for
the specified lines of sight (see **POLAR_LOS** and
**AXIAL_LOS** keywords in
:ref:`control file <inputs_control>`) if the keyword
**TAU1** is set to *Yes* in the
:ref:`control file <inputs_control>`.
The name of the files has the format **Tau_X_Y**,
with **X** and **Y** the indexes of the polar and azimuthal
angles of the line of sight in the order they are specified,
looping over both of them, with the polar
angle the slower loop. In inversion, if requested, it is named
just **Tau**.

.. include:: tau_file.rst


.. _cols-tt_file:

Cols-TT
=======

The **Cols-TT** file contains the rates of inelastic collisions
between the atomic terms of the :ref:`model atom <atom_file>`.
This file is created only if the keyword **KEEP_COLS** is set
to *Yes* in the :ref:`control file <inputs_control>`.

.. include:: colstt_file.rst


.. _cols-ll_file:

Cols-LL
=======

The **Cols-LL** file contains the rates of inelastic collisions
between the energy levels of the :ref:`model atom <atom_file>`.
This file is created only if the keyword **KEEP_COLS** is set
to *Yes* in the :ref:`control file <inputs_control>`.

.. include:: colsll_file.rst


.. _damping_file:

Damping
=======

The **Damping** file contains the broadening parameter of the
bound-bound radiative transitions for all the active
(see **ATOM_INPUT** in :ref:`control file <inputs_control>`)
:ref:`atom models <atom_file>`.
This file is created only if the keyword **KEEP_DAMP** is set
to *Yes* in the :ref:`control file <inputs_control>`.

.. include:: damping_file.rst


.. _background_file:

Background
==========

The **Background** file contains the absorptivity, scattering
coefficient, and emissivity of the background continuum.
This file is created only if the keyword **KEEP_BACK** is set
to *Yes* in the :ref:`control file <inputs_control>`. It is
written in binary.

.. include:: background_file.rst


.. _aparam_file:

avdwparam & astkparam
=====================

The **avdwparam** and **astkparam** files contain the values
that one should put in the broadening columns for the
bound-bound transitions in a :ref:`model atom <atom_file>`
in order to get the broadening values closest to the current
computed values but using the *param* type for the Van der Waals
broadening and the parametric input (negative number) for the
Stark broadening (see :ref:`model atom <atom_file>` for further
details). These files are created if the keyword **KEEP_APARAM**
is set to *Yes* in the :ref:`control file <inputs_control>`.
Both files are written in ASCII format.

.. include:: avdwparam_file.rst
.. include:: astkparam_file.rst


.. _atmosdat_file:

atmos.dat
=========

The **atmos.dat** file contains all the physical quantities of the
:ref:`model atmosphere <atmo_file>`, including those provided
in the input model **ATMO_INPUT** in the :ref:`control file <inputs_control>`
and the ones computed in runtime. This file is created if the
keyword **KEEP_ATMO** is set to *Yes* in the
:ref:`control file <inputs_control>`. It is written in ASCII
format and includes a header describing the content and units of the
different columns.

.. include:: atmosdat_file.rst


.. _atmoatmos_file:

atmo.atmos
==========

The **atmo.atmos** file is an :ref:`atmospheric model <atmo_file>`
ready to be used as an input atmosphere. It is created is the keyword
**UPDATE_ATMOS** is set to any of the possible options except *No*.
It is a file with the :ref:`input model atmosphere <atmo_file>`
format for the 1D or 1.5D synthesis.


.. _mrc_file:

MRCI & MRC
==========

The **MRCI** and **MRC** files contain information about the convergence
of the iterative problem in the only intensity and the polarized cases,
respectively. They are written in ASCII format for the 1D synthesis and
in binary for the 1.5D synthesis.

.. include:: mrci_file.rst
.. include:: mrc_file.rst


.. _col_file:

COL
===

The **COL** file is a log of the defined collisions for each active
atom (**ATOM_INPUT** in :ref:`control file <inputs_control>`) and it is written
in ASCII format

.. include:: col_file.rst


.. _mpiinfo_file:

MPI_info
========

The **MPI_info** file contains information about the workload distribution
between the processes. It is written in ASCII format.

.. include:: mpiinfo_file.rst


.. _ram_file:

RAM file
========

The **RAM** file contains information about the memory that each CPU
allocates to store to store background quantities, the radiation field
tensors, the Stokes parameters, the exponentials for the bound-free
radiative transfer coefficients (**PIRAM** in
:ref:`control file <inputs_control>`), Voigt profiles (**VOI_IRAM/VOI_PRAM**
in :ref:`control file <inputs_control>`), redistribution
functions (**RED_IRAM/RED_PRAM** in :ref:`control file <inputs_control>`)
and interpolation data for the redistribution (**INT_IRAM/INT_PRAM**
in :ref:`control file <inputs_control>`). It is written in ASCII format.

.. include:: ram_file.rst


.. _timer_global_file:

Global Timer File
=================

This is an ASCII file with the time, in seconds, at which the code
reaches certain points in the run. Because this file is mostly for
debugging purposes, these points will not be specified here. Instead,
you can grep 'report_time' in the *"Source"* folder.


.. _timer_mpi_file:

MPI Timer File
==============

This is an ASCII file with the time, in seconds, at which each process
sends its data to the master. The file has the default name
*"time_mpiI_ID"* or *"time_mpi_ID"* for the intensity and polarization
formal solutions, respectively, and where ID is a number identifying
the run that is randomly generated each time.

The first line of the file is written by the master when entering the
formal solution, and the rest follows the format:


::

    cpu iteration [prd_iteration] time

where *"cpu"* is the process index, *"iteration"* is the iteration,
*"prd_iteration"* is the PRD iteration only present in the intensity formal
solution, and *"time"* is the time in seconds from the begining of the
run.

.. _units_file:

Units
=====

.. include:: units.rst
