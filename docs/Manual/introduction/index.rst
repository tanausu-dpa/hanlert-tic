.. _intro:

*************
Introduction
*************

General information
===================

What is HanleRT
---------------

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

************
Terms of use
************

You can use HanleRT-TIC freely in your scientific research.

In case you publicly present your results obtained with the help of the
of the synthesis capability of the code (HanleRT), we ask you to quote the paper
`del Pino Alemán et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016ApJ...830L..24D/abstract>`_.
If using the inversion capabilities (HanleRT-TIC), we ask you to also quote
`Li et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022ApJ...933..145L/abstract>`_.
In both cases, please, include a link to the code repository
`https://gitlab.com/TdPA/hanlert-tic <https://gitlab.com/TdPA/hanlert-tic>`_.

Authors
=======

| Tanausú del Pino Alemán 
| Hao Li
| Roberto Casini

Contributors
============

| Ricky Egeland
| John Dennis

.. [1] While at some point there was support for OpenMP, significant structural changes have rendered the OpenMP inneficient at best (with PRD), and wrong at worst. Support for OpenMP will come back in a future version.

