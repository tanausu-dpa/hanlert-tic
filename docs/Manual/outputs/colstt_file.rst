This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Cols-TT** is the following:

    1. The two characters *ct*.
    2. Integer with the number of active atoms (see **ATOM_INPUT** in the
       :ref:`control file <input_file>`).
    3. Integer with the number of height nodes.
    4. For each atom:

      * Integer with the number of atomic terms.
      * For each term of origin and destination, double precision array with
        the rate of inelastic collisions from the origin to the destination
        term for each height node.

    .. code-block:: none

       |-----------------------------+------------------------------------------------------|
       | 2 * (char)                  | The characters "ct"                                  |   
       | 1 * (int)                   | Number of active atoms in the problem [NA]           |   
       | 1 * (int)                   | Number of height nodes [NZ]                          |   
       |------------------------------------------------------------------------------------|
       | For each atom:              |                                                      |   
       |   1 * (int)                 | Number of atomic terms [NT]                          |   
       |                             |                                                      |   
       |   NT * NT * NZ * (double)   | Rate of inelastic collisions from the origin         |
       |                             |   (fastest) to the destination term for each height  |
       |                             |   node (slowest)                                     |   
       |------------------------------------------------------------------------------------|

