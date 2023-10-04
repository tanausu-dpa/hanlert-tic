This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Jout** file is the following:

    1. The two characters *bj*.
    2. Integer which is 1 if the problem includes stimulated emission
       and 0 otherwise.
    3. Integer with the number of height nodes of the problem.
    4. Integer with the number of active atoms (see **ATOM_INPUT** in the
       :ref:`control file <input_file>`) of the problem.
    5. Integer with the maximum number of radiative bound-bound transition
       for the individual atoms.
    6. Double precision array with the height or tau values
       (see :ref:`atmospheric model <atmo_file>`) for each node.
    7. For each atom:

       * The number of bound-bound radiative transitions.
       * For each bound-bound radiative transition, value of the multipole
         :math:`K` up to 2, and all possible :math:`Q` values, double
         precision complex array with the :math:`J^K_Q`
         radiation field tensor for each height node for the given
         transition, integrated over the absorption profile.

    8. If and only if there is stimulated emission, for each atom, bound-bound
       radiative transition, value of the multipole :math:`K` up to 2, and
       all possible :math:`Q` values, double precision complex array with the
       :math:`J^K_Q` radiation field tensor for each height node for the given
       transition, integrated over the emission profile.

    .. code-block:: none

       |-----------------------------+------------------------------------------------------|
       | 2 * (char)                  | The characters "bj"                                  |   
       | 1 * (int)                   | 1 if stimulated emission, 0 otherwise [ST]           |   
       | 1 * (int)                   | Number of height nodes [NZ]                          |   
       | 1 * (int)                   | Number of active atoms in the problem [NA]           |   
       | 1 * (int)                   | Maximum number of radiative bound-bound transitions  |   
       |                             |   for an atom in the file                            |   
       | NZ * (double)               | Height or optical depth values for the height nodes  |   
       |------------------------------------------------------------------------------------|
       | For each atom:              |                                                      |   
       |    1 * (int)                | Number of bound-bound radiative transitions          |   
       |                             |                                                      |   
       |   For each bound-bound      |                                                      |   
       |    radiative transition:    |                                                      |   
       |   For all K up to 2:        |                                                      |   
       |   For all possible Q for    |                                                      |   
       |    for the given K:         |                                                      |   
       |                             |                                                      |   
       |      NZ * (double complex)  | Radiation field tensor for each height node          |   
       |                             |   integrated over the absorption profile             |   
       | If ST:                      |                                                      |   
       |   For each atom:            |                                                      |   
       |   For each bound-bound      |                                                      |   
       |     radiative transition:   |                                                      |   
       |   For all K up to 2:        |                                                      |   
       |   For all possible Q for    |                                                      |   
       |     for the given K:        |                                                      |   
       |                             |                                                      |   
       |       NZ * (double complex) | Radiation field tensor for each height node          |   
       |                             |   integrated over the emission profile               |   
       |                             |                                                      |   
       |------------------------------------------------------------------------------------|

