.. _installation:

************
Installation
************

Getting the code
================

The code is available on a GitLab hosted repository:
`https://gitlab.com/TdPA/hanlert-tic <https://gitlab.com/TdPA/hanlert-tic>`_.
In order to clone it, just do:: 

    git clone https://gitlab.com/TdPA/hanlert-tic.git

or using SSH::

    git clone git@gitlab.com:TdPA/hanlert-tic.git

.. _folder_tree:

Folder Tree
===========

The directory structure of a clean version of HanleRT-TIC is as follows:

::

    .
    |-- data
    |  |-- Atmos
    |  |-- Atoms
    |  |   |-- Background
    |  |   |-- CaII
    |  |   |-- MgII
    |  |   `-- SrI
    |  |-- Data
    |  `-- Molecules
    |-- data-user
    |-- docs
    |   `-- Manual
    |-- objects
    |-- examples
    |  |-- 1D
    |  |-- 15D
    |  |-- CLE
    |  `-- Inversion
    |-- Outputs
    |-- Resources
    |-- Source
    `-- Tools

* The **data/Atmos** directory contains the C and P models of
  `Fontela et al. (1993) <https://ui.adsabs.harvard.edu/abs/1993ApJ...406..319F/abstract>`_
  to be used for 1D synthesis.

* The **data/Atoms** directory contains some sample atomic models.

  * The **Background** directory contains the typical atoms to include
    in the calculation of background continuum opacity, scattering,
    and emission coefficient.
  * The **CaII** directory contains some sample model atoms suited for
    the modeling of the Ca II H, K, and infrared triplet lines.
  * The **MgII** directory contains some sample model atoms suited for
    the modeling of the Mg II h and k lines.
  * The **SrI** directory contains some sample model atoms suited for
    the modeling of the Sr I line at 460.7 nm.

* The **data/Data** directory contains an example *data* file for
  the inversion example.

* The **data/Molecules** directory contains some sample molecule models 
  typically included in the calculation of the chemical equilibrium and
  of background continuum opacity, scattering, and emission coefficient.

* The **data-user** directory is not tracked in the respository. It is
  recommended to store your atomic, molecular, or atmospheric models in here.

* The **docs** directory contains:
  
  * The **Manual** directory contains the *"source code"* for this User
    Manual.

* The **examples** directory contain some relatively light-weight examples
  that the user can run to check that the code is working correctly and
  to learn how to run the code for different modes.

* The **objects** directory is using to hold the objects and modules
  during compilation. Its contain is not tracked by git.

* The **Outputs** directory can be used to run and store results, as
  its content is not tracked by git.

* The **Resources** directory contains some auxiliary files, including
  abundances, partition functions, and broadening parameters.

* The **Source** directory contains the source code.

* The **Tools** directory contains auxiliar routines to facilitate
  dealing with the output of the code.

.. _dependencies:

Dependencies
============

In order to compile and run the code, the following software and libraries
are required:

* A C compiler.
* A fortran compiler.
* Python 2.7 or Python 3.*.
* An MPI library with fortran support.
* The cfitsio library to handle fits files (optional, only for fits support).
* The lapack and blas libraries (or equivalents, such as Intel's mkl or openblas).
* The OpenMP library if you choose to compile and run the code with OpenMP.

.. note::
   HanleRT-TIC has been mostly tested with the *gcc* compiler, *OpenMPI* version 4,
   and *openblas*. Please, report about incompatibilities that you find with other
   software combinations.

There are several ways to acquire and install all these dependencies.

* In supercomputers, all these libraries are readily available. For example, in
  `LaPalma supercomputer <https://www.iac.es/en/science-and-technology/technology/technical-facilities/lapalma-supercomputer>`_
  it is enough to load the modules (assuming that fits support is desired)::

    module load gnu openmpi/gnu scalapack/2.0.2 python cfitsio

  and the environment would be ready for compilation.

.. warning::
   Remember that when using a queue system such as *slurm*, which is the case
   for any supercomputer, the same modules need to be loaded in the *slurm*
   script.

* In a computer where the user has administrator rights, the user may opt for
  installing the dependencies in the system. While this is not the method I
  would personally recommend, it is the most straightforward, albeit also
  the one with more risk of breaking with updates or breaking depencencies.

* A more elegant way to set-up the dependencies environment is using a
  package and environment manager such as `Spack <https://spack.io/>`_.
  A basic guide to install the dependencies using Spack is described in
  :ref:`Spack installation <spack_installation>`

.. note::
   If should be relatively straightforward to install the dependencies with
   other managers such as *conda*, but we have not tested such installation
   and thus we do not include instructions.

.. _spack_installation:

Environment installation with Spack
===================================

**Date**: 29 September 2023

Getting Spack
-------------

From `https://spack.io/ <https://spack.io/>`_, 
*Spack is a package manager for supercomputers, Linux, and macOS. It makes installing scientific software easy. Spack isn’t tied to a particular language; you can build a software stack in Python or R, link to libraries written in C, C++, or Fortran, and easily swap compilers or target specific microarchitectures. Learn more here* (`<https://spack.io/about/>`_).
It is thus a very convenient tool to install and manage dependencies for specific projects without the
need to modify the actual system and global dependencies.

This installation guide is thought for a computer without queue system. It is possible to also use *Spack* in a system using
*slurm* for which the user does not have administration privileges, but it can become system specific and thus it is not
included in here. However, *Spack* has its own `documentation <https://spack.readthedocs.io/en/latest/>`_ which can be
useful for more complicated cases.

*Spack* has its own dependencies, which are listed in this `link <https://spack.readthedocs.io/en/latest/getting_started.html>`_, and
are easily satisfied. In the provided link, the *Spack* team gives instructions on setting up the dependencies for several
operative systems.

The first step is to clone the *Spack* repository::

    git clone -c feature.manyFiles=true https://github.com/spack/spack.git

This will create a *Spack* directory. In order to activate *Spack* in a terminal, write the following command,
(this can be configured as an *alias* for easiness of use):

  * In bash/zsh/sh::

      . <path_to_spack_directory>/share/spack/setup-env.sh

  * In tcsh/csh::

      source <path_to_spack_directory>/share/spack/setup-env.csh

  * In fish::

      . <path_to_spack_directory>/share/spack/setup-env.fish

Now all *Spack* commands are available. The first time using *Spack*, it is recomended to run::

   spack spec zlib

in order to set-up *clingo* (see `documentation <https://spack.readthedocs.io/en/latest/>`_). Check the corresponding section in
the documentation if unable to bootstrap *clingo* from pre-built binaries.

The command::

  spack compilers

shows the list of compilers avaiable to *Spack*. If an existing compiler is missing, you can run::

  spack compiler find

either without arguments, for autodetection, or with a path where to look for, if it is known.

Installing a compiler with Spack
--------------------------------

If the user wants to use one of the compilers already available in the system, this section can be
skipped. Here we describe how to install a compiler with *Spack*. We install *gcc V9.5.0* as an example::

   spack install gcc@9.5.0
   spack load gcc@9.5.0
   spack compiler find
   spack unload gcc@9.5.0

Installing the dependencies
---------------------------

We will use a *yaml* file to configure the environment. If using a *<system_compiler>*, the file can be as
simple as::

    spack:
      specs: [openmpi, openblas, cfitsio]
      concretizer:
        unify: true
      view: true
      packages:
        all:
          compiler: [<system_compiler>]

The user can add *openmp* to the list of *specs* if planning to compile the *OpenMP* [1]_ support.

If using a compiler installed with *Spack*, the *yaml* file is slightly more complicated. First, we
need the path to the compiler, e.g.::

    spack compiler info gcc@9.5.0
    gcc@9.5.0:
        paths:
            cc = /spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.5.0/gcc-9.5.0-zo4ue572rugd5zxbip7k25a4zc4sr2ts/bin/gcc
            cxx = /spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.5.0/gcc-9.5.0-zo4ue572rugd5zxbip7k25a4zc4sr2ts/bin/g++
            f77 = /spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.5.0/gcc-9.5.0-zo4ue572rugd5zxbip7k25a4zc4sr2ts/bin/gfortran
            fc = /spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.5.0/gcc-9.5.0-zo4ue572rugd5zxbip7k25a4zc4sr2ts/bin/gfortran

we need the full path up to one level above the *bin* directory (identical for the four entries). For this example, the
*yaml* file would read::

  spack:
    specs: [gcc@9.5.0, openmpi%gcc@9.5.0, openblas, cfitsio]
    concretizer:
      unify: true
    view: true
    packages:
      all:
        compiler: [gcc@9.5.0]
      gcc:
        buildable: False
        externals:
        - spec: gcc@9.5.0
          prefix: /spack/opt/spack/linux-ubuntu18.04-skylake_avx512/gcc-7.5.0/gcc-9.5.0-zo4ue572rugd5zxbip7k25a4zc4sr2ts

where we have used the path from before in the *prefix* field.

Now everything is ready to install the environment. The commands::

  spack env create hanlert <path_to_yaml_file>
  spack env activate -p hanlert
  spack concretize -f
  spack install
  spack env deactivate

will create an environment named *hanlert* with the dependencies listed in the *specs* in the *yaml* file, using the
compiler specified in the *compiler* field of the same file.

.. warning::
   Be aware that the installation process can take a really long time. If more than one CPU is avaiable, consider adding
   "-jN", with "N" the number of CPU, to parallelize the installation.

Now, before compiling or running the code, it is enough to write (remember to activate *Spack* first)::

  spack env activate -p hanlert

.. _compilation:

Compilation
===========

The compilation of HanleRT-TIC involved to stage, *configuration* and actual
*compilation*.

.. _configuration:

1. **Configuration**

In the root directory of HanleRT-TIC there is a file named *configure*. This file
contains a script to generate the Makefile. The list of options can be accessed
by running ``./configure --help``. The most crucial fields to consider are:

* --lflags: The external libraries (lapack and blas, mkl, or openblas, and cfitsio)
  must be specified here in the Makefile language. For example,
  ``--lflags='-L/apps/SCALAPACK/2.0.2/lib -lreflapack -lrefblas -L/apps/CFITSIO/3.420/lib -lcfitsio'`` are the linking flags needed in the LaPalma supercomputer.
* --flags: The compilation flags. They can be customized, but the configure
  allows for some hard-coded configurations. If the computer does not have
  a list of preferred or specific flags, it is recommended to run the configuration
  with with the ``--medium3`` option.

.. note::
   If the dependencies were installed in a *Spack* environment, the linking options can be
   semi-automatized, using::

       --lflags='-L<spack_root_directory_path>/var/spack/environments/hanlert/.spack-env/view/lib -Wl,-rpath,<spack_root_directory_path>/var/spack/environments/hanlert/.spack-env/view/lib -lopenblas -lcfitsio'

   However, it then becomes necessary to execute the ``ld`` command on the compiled executable
   before running.

.. note::
   If the fits support is not needed, the *cfitsio* library is not necessary. To be able
   to compile the code without fits support, add the ``--nofits`` option when running ``./configure``.

.. warning::
   If dealing with an old version of the OpenMPI library (older than version 3.0),
   it is required to run the configuration with the ``--oldmpi`` option or the code
   will not compile.

2. **Compilation**

Once the configuration has been run successfully, just run ``make`` in the root
directory to start the compilation. The objects and modules will be stored in the
**objects** directory.

.. note::
   By default, the executable will be created in the root directory and will be
   named **hanlert**, unless you have specified a different path and name during
   the configuration using the ``--prefix=<path>`` option.

.. note::
   Run ``make clean`` in order to remove the compiled objects, modules, and
   excutable, or run ``make clean_obj`` to remove the compiled objects and
   modules.

.. [1] While at some point there was support for OpenMP, significant structural changes have rendered the OpenMP inneficient at best (with PRD), and wrong at worst. Support for OpenMP will come back in a future version.
