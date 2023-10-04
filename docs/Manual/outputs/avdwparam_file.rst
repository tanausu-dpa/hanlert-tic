The **avdwparam** file is created for each atom and its real name is
**XX.avdwparam**, with **XX** the symbol of the element (e.g., H, Mg, Ca, etc.).

The first line of the file describes its use:

::

  This file contains the parameters that you have to put in the model atom with the parameter option to mimic the current Van dew Waals approximation

The second line specifies to what atom the file corresponds:
 

::

  Atom: CA


The third line shows the legend of the data coming below:

::

  Transition                             approxim.             A_H             B_H            A_He            B_He


Then, a line for each bound-bound radiative transition in the
:ref:`atomic model <atom_file>` follows, with the following
information for each of them:

1. The index of the bound-bound radiative transition.
2. The indexes of the upper and lower terms or levels (depending
   on the type of :ref:`atomic model <atom_file>`).
3. The approximation for the Van der Waals broadening that was
   used to compute the values.
4. The four parameters that should be written in the :ref:`atomic
   model <atom_file>` in order to get the closest value for the
   Van der Waals broadening of the bound-bound radiative transition
   using the parametric approach (*param*; see bound-bound radiative
   transitions in :ref:`atomic model <atom_file>` for further
   details).

::

         1          2-->           1     barklem  8.13865304E-02  3.80878279E-01  7.55558175E-01  3.00000000E-01


