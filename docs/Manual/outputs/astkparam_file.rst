The **astkparam** file is created for each atom and its real name is
**XX.astkparam**, with **XX** the symbol of the element (e.g., H, Mg, Ca, etc.).

The first line of the file describes its use:

::

  This file contains the parameters that you have to put in the model atom with the parameter option to mimic the current Stark approximation

The second line specifies to what atom the file corresponds:
 

::

  Atom: CA


The third line shows the legend of the data coming below:

::

  Transition                                         C


Then, a line for each bound-bound radiative transition in the
:ref:`atomic model <atom_file>` follows, with the following
information for each of them:

1. The index of the bound-bound radiative transition.
2. The indexes of the upper and lower terms or levels (depending
   on the type of :ref:`atomic model <atom_file>`).
3. The parameters that should be written in the :ref:`atomic
   model <atom_file>` in order to get the closest value for the
   Stark broadening of the bound-bound radiative transition using
   the parametric approach (inputting a negative number; see
   bound-bound radiative transitions in
   :ref:`atomic model <atom_file>` for further details).

::

                  1          2-->           1  2.52165562E-07


