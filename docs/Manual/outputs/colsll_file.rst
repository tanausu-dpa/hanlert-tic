This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Cols-LL** is the following:

    1. The two characters *cl*.
    2. Integer with the number of active atoms (see **ATOM_INPUT** in the
       :ref:`control file <input_file>`).
    3. Integer with the number of height nodes.
    4. For each atom:

      * Integer with the number of energy levels.
      * For each level of origin and destination, double precision array with
        the rate of inelastic collisions from the origin to the destination
        level for each height node.

    .. code-block:: none

       |-----------------------------+------------------------------------------------------|
       | 2 * (char)                  | The characters "cl"                                  |   
       | 1 * (int)                   | Number of active atoms in the problem [NA]           |   
       | 1 * (int)                   | Number of height nodes [NZ]                          |   
       |------------------------------------------------------------------------------------|
       | For each atom:              |                                                      |   
       |   1 * (int)                 | Number of energy levels [NL]                         |   
       |                             |                                                      |   
       |   NL * NL * NZ * (double)   | Rate of inelastic collisions from the origin         |
       |                             |   (fastest) to the destination level for each height |
       |                             |   node (slowest)                                     |   
       |------------------------------------------------------------------------------------|
