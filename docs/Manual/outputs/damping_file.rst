This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **damping** file is the following:

    1. The two characters *da*.
    2. Integer with the number of active atoms (see **ATOM_INPUT** in the
       :ref:`control file <input_file>`).
    3. Integer with the number of height nodes.
    4. For each atom:

       * Integer with the number of bound-bound radiative transitions.
       * For each bound-bounda radiative transition, double precision array
         with the damping parameter for each height node.

    .. code-block:: none

       |-----------------------------+------------------------------------------------------|
       | 2 * (char)                  | The characters "da"                                  |   
       | 1 * (int)                   | Number of active atoms in the problem [NA]           |   
       | 1 * (int)                   | Number of height nodes [NZ]                          |   
       |------------------------------------------------------------------------------------|
       | For each atom:              |                                                      |   
       |   1 * (int)                 | Number of bound-bound radiative transitions [NT]     |   
       |                             |                                                      |   
       |   NT * NZ * (double)        | Damping parameter for each bound-bound radiative     |
       |                             |   transition (slowest) and height node (fastest)     |   
       |                             |                                                      |   
       |------------------------------------------------------------------------------------|

