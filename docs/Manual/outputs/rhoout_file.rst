This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Rhoout** file is the following:

    1. The two characters *br*.
    2. Integer with the number of height nodes of the problem.
    3. Integer with the number of active atoms (see **ATOM_INPUT** in the
       :ref:`control file <input_file>`) of the problem.
    4. Double precision array with the height or tau values
       (see :ref:`atmospheric model <atmo_file>`) for each node.
    5. For each atom:

      * Double precision array with the total atom population for each height node.
      * Integer with the number of atomic terms.
      * For each term:

          * Integer with the number of energy levels within the term.
          * For each pair of energy levels within the term:

            * Two Integers with twice the total angular momentum of
              the two energy levels.
            * For all possible multipole :math:`K` values, for all possible
              multipole :math:`Q` values, for each height node:

                * Double complex with the value of the
                  :math:`\rho^K_Q(J,J')` density matrix element in
                  the vertical reference for the given height.
                * Integer with the internal code flag which
                  indicates if a value is too small to be
                  taken into account in the calculation of the
                  radiative transfer coefficients (value 1 if
                  very small, 0 otherwise).

    .. code-block:: none

       |------------------------------+------------------------------------------------------|
       | 2 * (char)                   | The characters "br"                                  |   
       | 1 * (int)                    | Number of height nodes [NZ]                          |   
       | 1 * (int)                    | Number of active atoms in the problem [NA]           |   
       | NZ * (double)                | Height or optical depth values for the height nodes  |   
       |-------------------------------------------------------------------------------------|
       | For each atom:               |                                                      |   
       |   NZ * (double)              | Total element number density population for each     |   
       |                              |   height node in cgs                                 |   
       |   1 * (int)                  | Number of atomic terms                               |   
       |                              |                                                      |   
       |   For each term:             |                                                      |   
       |                              |                                                      |   
       |     1 * (int)                | Number of energy levels within the term              |   
       |                              |                                                      |   
       |     For each energy level    |                                                      |   
       |      within the term:        |                                                      |   
       |     For each energy level    |                                                      |   
       |      within the term:        |                                                      |   
       |                              |                                                      |   
       |       1 * (int)              | Total angular momentum of the energy level of the    |   
       |                              |   slow loop multiplied by two                        |   
       |       1 * (int)              | Total angular momentum of the energy level of the    |   
       |                              |   fast loop multiplied by two                        |   
       |       For all possible K for |                                                      |   
       |        the given levels:     |                                                      |   
       |       For all possible Q for |                                                      |   
       |        the given K:          |                                                      |   
       |       For each height:       |                                                      |   
       |                              |                                                      |   
       |         1 * (double complex) | Density matrix element                               |   
       |         1 * (int)            | Flag to indicate no contribution to the radiative    |   
       |                              |   transfer coefficients due to small value           |   
       |                              |                                                      |   
       |-------------------------------------------------------------------------------------|

