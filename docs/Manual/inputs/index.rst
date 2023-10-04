.. _inputs:

******
Inputs
******

Running HanleRT-TIC requires the set up of several input files:
* :ref:`Control file <input_file>` with the information and configuration of the run.
* :ref:`Atomic files <atom_file>` with the model atoms to consider.
* :ref:`Atmospheric file <atmo_file>` with the model atmosphere (not required for the inversion).
* :ref:`Data file <data_file>` with the spectroscopic or spectro-polarimetric data to invert (only required for the inversion).
 
There are other files that are either used internally or that can be specified
in the :ref:`control file <input_file>`. In the following, we describe the format
of all these inputs.

.. _input_file:

Control File
============

The control or input file contains the list of parameters that fully configure
the HanleRT-TIC run. By default, the code looks for a file called
**INPUT** within the folder from where the command is being executed. Nevertheless,
the path to the input file can be specified in the command line after the
executable (see :ref:`running`).

The input file is in ASCII and allows for comments. Everything after a ``!``
symbol will be ignored until a line break happens. The configuration
is expected in the format ``KEYWORD = VALUE``. The allowed keywords are listed just below.
Some of the keywords are additive, that is, several lines can be included
with the same keyword to have more of the same. For example, multiple lines
with the ``ATOM_INPUT`` keyword will be interpreted as a request to include
all of the specified model atoms.

Some of the ``KEYWORDS`` are exclusive to different types of runs. We provide the
following lists of ``KEYWORDS``.

* :ref:`Full list <control_full>`.
* :ref:`1D synthesis (basic only) <control_1D_basic>`.
* :ref:`1D synthesis (full) <control_1D>`.
* :ref:`1.5D synthesis (basic only) <control_15D_basic>`.
* :ref:`1.5D synthesis (full) <control_15D>`.
* :ref:`CLE synthesis (basic only) <control_CLE_basic>`.
* :ref:`CLE synthesis (full) <control_CLE>`.
* :ref:`Inversion (basic only) <control_INV_basic>`.
* :ref:`Inversion (full) <control_INV>`.

.. note::
   A note on multiple cycles. Usually, inversion codes allows to configure
   multiple cycles for a single run changing, e.g., the number of nodes
   for each variable. In HanleRT-TIC we have opted for not allowing this
   kind of run (with the exception of the *sequential* and *sequential-magnetic*
   types of inversions, see the **TYPE_INVERSION** keyword in the control
   file). The inversions can get really expensive from the computational
   point of view, and we prefer to avoid the possibility of losing enormous
   amounts of computing time due to wrong or not optimal configuration
   for series of cycles. Intead, any result file can be used to initialize
   a new cycle for the same `data file <data_file>`. We encourage the user
   to check the result of each cycle before launching the next one.


.. _atmo_file:

Atmospheric Model File
======================

The format is different depending on the **RUN_MODE** in the :ref:`control file <input_file>`.

1D
--

We have chosen to keep the compatibility with the atmospheric model format
of the MULTI and RH codes
(`1986UppOR..33.....C <https://ui.adsabs.harvard.edu/abs/1986UppOR..33.....C/abstract>`_,
`2001ApJ...557..389U <https://iopscience.iop.org/article/10.1086/321659>`_,
respectively), but with some extensions to allow for a wider variety of inputs.
The atmospheric file is in ASCII and allows for comments,
everything after a ``!`` or a ``*`` symbol will be ignored until a line change
happens (the latter for compatibility reasons).

.. include:: atmo1d_file.rst

1.5D
----

The atmospheric file for the *1.5D Synthesis* is written in binary.

.. include:: atmo15d_file.rst


CLE
---

The atmospheric file for the *CLE Synthesis* is written in binary.

.. include:: atmocle_file.rst

.. _data_file:

Data File
=========

The data file for the inversion can be specified in *.fits.* or
binary formats.

.. include:: data_file.rst

.. _atom_file:

Atomic Model File
=================

The atomic files contains all the needed atomic information, in ASCII format,
to model its emergent spectrum in a given model atmosphere. The current version
of HanleRT-TIC admits two types of model atoms, namely multi-level and
multi-term models. The atomic file allows for comments, everything after a
"!" symbol will be ignored until a line change happens.
In order to explain the atomic format we will look at one
example for each of the allowed types.

.. include:: atom_file.rst


.. _mol_file:

Molecule Model File
===================

The molecule files contains all the needed molecular information
to solve the chemical equilibrium equations. It is in ASCII and allows for
comments, everything after a "!" will be ignored until a line change
happens. The format has been chosen in such a way that the conversion
from the molecular modes of the RH code
(`2001ApJ...557..389U <https://iopscience.iop.org/article/10.1086/321659>`_)
is straightforward.

.. include:: mol_file.rst


.. _bfield_file:

Magnetic Field File
===================

The magnetic field file contains the stratification of the magnetic field
vector. This file is only valid for the 1D synthesis or inversion
problems (will be the initialization in the latter). It is in ASCII format
and allows for comments, everything after a "!" will be ignored until a
line change happens.

.. include:: bfield_file.rst

.. _kur_file:

Kurucz File
===========

.. include:: kur_file.rst

.. _lte_file:

LTE line File
=============

.. include:: lte_file.rst

.. _ion_file:

Ion File
========

.. include:: ion_file.rst

.. _spec_file:

Spectral File
=============

.. include:: spec_file.rst


.. _fudge_file:

Background Fudge File
=====================

The background fudge file contains the background opacity fudge factors
to enhance the background continuum opacity at the indicated wavelengths.
It is an ASCII file and allows for comments, everything after a "!" symbol
will be ignored until a line change happens. The format has been chosen
in such a way that the convertion from the fudge files of the RH code
(`2001ApJ...557..389U <https://iopscience.iop.org/article/10.1086/321659>`_)
is almost straightforward.

.. include:: fudge_file.rst

.. note::

   A fudge file (*fudge.dat*) with the factors from Bruls 1992
   (*Formation of diagnostic lines in the solar spectrum*,
   PhD Thesis, Utrecht University) is included in the :ref:`Resources
   folder <folder_tree>`.


.. _partition_function_file:

Partition Function File
=======================

.. include:: partition.rst

.. include:: partition_function_file.rst


.. _abundances_file:

Abundances File
===============

This file allows the user to change the abundances of some or all atomic
species. Abundances that are not specified as inputs take their default
value. The abundance file is written in binary.

.. include:: abundances_file.rst


.. _barklem_file:

Barklem File
============

The Barklem files contain the information necessary for the
calculation of the Van der Waals contribution of the broadening
of atomic spectral lines if the *barklem* type is selected
(see :ref:`atomic model <atom_file>`). The format of the default
files is different to that of the input files. The default files
are written in ASCII, while the input files are written in binary.

Default files
-------------

.. include:: barklem.rst

Input files
-----------

.. include:: barklem_file.rst


.. _asymm_file:

Ad-hoc Asymmetries File
=======================

.. include:: asymm_file.rst



.. _wavelength_file:

Wavelength File
===============

The wavelength file contains a list of wavelengths that must be
included in the solution of the radiative transfer problem and
on the emergent Stokes profiles. This file has the following
binary format.

.. include:: wavelength_file.rst


.. _other_file:

Harcoded Quantities
===================

.. include:: hardcoded.rst

