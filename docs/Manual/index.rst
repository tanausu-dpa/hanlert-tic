.. HanleRT-TIC documentation master file, created by
   sphinx-quickstart on Fri Sep 11 12:15:21 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. rst format CheatSheet:
   http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html#comments-and-aliases 

HanleRT-TIC
===========

HanleRT-TIC is a numerical code for solving the problem of the generation and
transfer of polarized radiation in one-dimensional (1D and 1.5D) model
atmospheres for arbitrary atomic models. HanleRT-TIC also allows for the
solution of the inversion problem with its Tenerife Inversion Code (TIC)
module.

The code can take into account atomic and scattering polarization, quantum
interference among energy levels pertaining to a given atomic term, partial
frequency redistribution effects due to the photon coherence in scattering
processes, the impact of dynamics, and the impact of magnetic fields via the
joint action of the Hanle and Zeeman effects (incomplete Paschen-Back regime).
More details about the physics in this ratiative transfer code can be found
in :ref:`physics`.

Apart from solving the radiation transfer problem in optically thick
plane-parallel model atmospheres, HanleRT-TIC allows for the solution of the
coronal line emission (CLE) problem, in which the radiation field is assumed
to be dominated by the underlying stellar disk.

HanleRT-TIC is written in standard Fortran 2008, parallelized with the OpenMPI
(`<https://www.open-mpi.org/>`_) and OpenMP [1]_ (`<https://www.openmp.org/>`_)
libraries. The code also has some parsing routines written in python to
allow for more flexible input formats.

The code has been made publicly available within the framework of the POLMAG
project funded by an Advanced Grant of the European Research Council (see the
`POLMAG webpage <http://research.iac.es/proyecto/polmag/>`_).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction/index
   installation/index
   running/index
   inputs/index
   inputs_control/index
   outputs/index
   physics/index
   parallelization/index
   issues/index
   publications/index


.. [1] While at some point there was support for OpenMP, significant structural changes have rendered the OpenMP inneficient at best (with PRD), and wrong at worst. Support for OpenMP will come back in a future version.
