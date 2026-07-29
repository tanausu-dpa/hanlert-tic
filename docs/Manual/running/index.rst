.. _running:

****************
Running the code
****************

Executing HanleRT-TIC
=====================

Running HanleRT-TIC (after a successful :ref:`installation`) consist on a
simple run command:

.. code-block:: bash

   mpirun {-np N} <executable> [control_file]

* The ``-np N`` indicate that the code is to be executed with **N** processes
  (CPUs or hyperthreads). While it is a required argument when directly running
  in the command line, it may be optional when using a queue system such as
  **slurm**.

* ``<executable>`` is the path to the executable file (called **hanlert** and
  located in the root directory by default). Note that the code can be run
  from any folder, and the output can be directed to any path specified in the
  :ref:`control file <input_file>`.

* ``control_file`` is an optional argument to specify the path to the 
  :ref:`control file <input_file>`. If not specified, HanleRT-TIC will assume
  that there is a :ref:`control file <input_file>` named **INPUT** in the
  running directory.

.. note::
   HanleRT-TIC needs to access the **Source** and **Resources** folders in
   running time. These folders can be automatically found is the running
   directory is up to four levels inside the root directory. You can run
   from a deeper directory or from outside the root directory if you add
   the files **spath** and **sreso** with a single plain text line indicating
   the path to the **Source** and **Resources** folders, respectively.

.. warning::
   Some of the algebra libraries that can be use to fulfill the :ref:`dependencies <dependencies>`
   come with automatic parallelization which can conclict with the code's parallelization,
   signifncantly slowing down the running time. We have found this behavior with the
   *openblas* library in particular. For this specific library, you must set the environment
   variable ``export OPENBLAS_NUM_THREADS=1`` before running HanleRT-TIC. In case that the
   code is running much slower than expected, check that the algebra library is working in
   serial mode.

.. note::
   Remember that there is a :ref:`faq` section that can solve some of the most typical
   doubts regarding printed messages or the progress of the calculations.

Preparing the inputs
====================

Most of the HanleRT-TIC inputs and outputs are unformatted binary files. This
allows for very efficient parallel access in the minimum amount of storage
space without adding further dependencies and coding effort. However,
this means that the inputs (outputs) must be written (read) very carefully
taking into account the byte count.

While :ref:`inputs` and :ref:`outputs` provide a detailed description of the
structure of all these files, we understand the added burden on the user.
For these reasons we have included tools in the **Tools** folder to help
with this endeavor:

* **Tools/hanlertmodlib.py** contains python classes to manage the creation
  of the following intput files:

  * The class *invi_class* allows creating a :ref:`data_file` for
    the inversion.
  * The class *atmo15D_class* allows creating an :ref:`atmo_file` with
    a 3D model for 1.5D synthesis.
  * The class *atmo1D_class* allows creating an :ref:`atmo_file` with
    a 1D model for 1D synthesis. While this file is among the few ASCII
    exceptions, the class can still be convenient.
  * The class *weight_class* allows creating a :ref:`inversion_weights_file` with
    customized wavelength dependent weights for the inversion.

* **Tools/hanlertmodhelper_inversion.py** shows a template for the application
  of the *invi_class* in **Tools/hanlertmodlib.py**.

* **Tools/hanlertmodhelper_synthesis15D.py** shows a template for the application
  of the *atmo15D_class* in **Tools/hanlertmodlib.py**.

* **Tools/hanlertmodhelper_synthesis1D.py** shows a template for the application
  of the *atmo1D_class* in **Tools/hanlertmodlib.py**.

* **Tools/notebook.ipynb** is a notebook with very simple example applications
  to build the abovementioned files.
