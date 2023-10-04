**General Polarized Case**

The format of the **Solution** file for the general polarized problem
is the following:

1. The two characters *sp*.
2. Integer with the number of frequencies of the problem.
3. Integer with the number of height nodes of the problem.
4. Integer with the number of polar quadrature nodes of the problem.
5. Integer with the number of azimuthal quadrature nodes of the problem.
6. Integer with the number of active atoms (see **ATOM_INPUT** in the
   :ref:`control file <input_file>`) of the problem.
7. Integer which is 1 if the problem is axially symmetric and 0
   otherwise.
8. Integer which is 1 if the problem includes stimulated emission
   and 0 otherwise.
9. Integer which is 1 if the redistribution function is angle-averaged
   and 0 otherwise.
10. For each atom:

  * Double precision array with the total atom population for each height
    node.
  * For each term, pair of energy levels within the term, possible multipole
    :math:`K` values, and possible multipole :math:`Q` values, double
    precision complex array with the value of the :math:`\rho^K_Q(J,J')`
    density matrix element in the vertical reference frame for each height node.

11. For each atom, bound-bound radiative transition, multipole :math:`K` values
    up to 2, and possible multipole :math:`Q` values, double precision complex
    array with with the :math:`J^K_Q` radiation field tensor for each height
    node for the given transition, integrated over the absorption profile.
    
12. If and only if there is stimulated emission,
    for each atom, bound-bound radiative transition, multipole :math:`K` values
    up to 2, and possible multipole :math:`Q` values, double precision complex
    array with with the :math:`J^K_Q` radiation field tensor for each height
    node for the given transition, integrated over the emission profile.

13. If and only if the problem is angle-averaged and not dynamic or if there
    is no partial frequency redistribution in any atomic transition line, for
    each height, frequency, multipole :math:`K` values up to 2, and possible
    multipole :math:`Q` values, a double complex number with the value of the
    :math:`J^K_Q` radiation field tensor.
   
14. If and only if the condition for 14 are not fulfilled, for each height, polar
    direction, azimuthal direction, and frequency, four double precision numbers
    with the *I*, *Q*, *U*, and *V* Stokes parameters.

.. code-block:: none

   |-----------------------------+------------------------------------------------------|
   | 2 * (char)                  | The characters "sp"                                  |   
   | 1 * (int)                   | Number of frequency nodes [NF]                       |   
   | 1 * (int)                   | Number of height nodes [NZ]                          |   
   | 1 * (int)                   | Number of polar quadrature nodes [NT]                |   
   | 1 * (int)                   | Number of azimuthal quadrature nodes [NP]            |   
   | 1 * (int)                   | Number of active atoms in the problem [NA]           |   
   | 1 * (int)                   | 1 if axially simmetryc, 0 otherwise [AX]             |   
   | 1 * (int)                   | 1 if stimulated emission, 0 otherwise [ST]           |   
   | 1 * (int)                   | 1 if angle-averaged redistribution, 0 otherwise [AA] |
   |------------------------------------------------------------------------------------|
   | For each atom:              |                                                      |   
   |   NZ * (double)             | Total element number density population for each     |   
   |                             |   height node in cgs                                 |   
   |   For each term:            |                                                      |   
   |   For each energy level     |                                                      |   
   |     within the term:        |                                                      |   
   |   For each energy level     |                                                      |   
   |     within the term:        |                                                      |   
   |   For all possible K for    |                                                      |   
   |     the given levels:       |                                                      |   
   |   For all possible Q for    |                                                      |   
   |     the given K:            |                                                      |   
   |                             |                                                      |   
   |     NZ * (double complex)   | Density matrix element for each height node          |   
   |                             |                                                      |   
   | For each atom:              |                                                      |   
   | For each bound-bound        |                                                      |   
   |   radiative transition:     |                                                      |   
   | For all K up to 2:          |                                                      |   
   | For all possible Q for      |                                                      |   
   |   for the given K:          |                                                      |   
   |                             |                                                      |   
   |   NZ * (double complex)     | Radiation field tensor for each height node          |   
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
   | If (AA and not dynamic) or  |                                                      |   
   |   not partial frequency     |                                                      |   
   |   redistribution:           |                                                      |   
   |                             |                                                      |   
   |   For each height:          |                                                      |   
   |   For each frequency:       |                                                      |   
   |   For all K up to 2:        |                                                      |   
   |   For all possible Q for    |                                                      |   
   |     for the given K:        |                                                      |   
   |                             |                                                      |   
   |     1 * (double complex)    | Radiation field tensor                               |   
   |                             |                                                      |   
   | else:                       |                                                      |   
   |                             |                                                      |   
   |   For each height:          |                                                      |   
   |   For each polar direction: |                                                      |   
   |   For each azimuthal        |                                                      |   
   |     direction:              |                                                      |   
   |   For each frequency:       |                                                      |   
   |                             |                                                      |   
   |     4 * (double)            | I, Q, U, and V Stokes parameters                     |   
   |                             |                                                      |   
   |------------------------------------------------------------------------------------|

