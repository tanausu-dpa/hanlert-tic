**Intensity Only**

The format of the only intensity **Solution** file and the **SolutionI**
file is the following:

1. The two characters *si*.
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
  * For each energy level, double precision array with the value of the
    density matrix element :math:`K=Q=0` for each height node.

11. For each atom and bound-bound radiative transition within the
    atom, double precision array with the value of the mean radiation
    field integrated over the absorption profile for each height node.

12. If and only if there is stimulated emission, for each atom and
    bound-bound radiative transition within the atom, double precision
    array with the value of the mean radiation field integrated over the
    emission profile for each height node.

13. For each atom and bound-free radiative transition within the atom, 
    double precision array with the value of the integrated photoionization
    terms (see equation 2.9 in 
    `1991A&A...245..171R <https://ui.adsabs.harvard.edu/abs/1991A%26A...245..171R/abstract>`_
    for the photoionization case) for the photoionization and recombination
    rates in the statistical equilibrium equations.

14. If and only if the problem is angle-averaged and not dynamic or if there
    is no partial frequency redistribution in any atomic transition line, for
    each height, double precision array the value of the mean radiation
    field for each frequency node.

15. If and only if the condition for 14 are not fulfilled, for each height,
    polar direction, and azimuthal direction, double precision array with
    the value of the intensity for each frequency node.

.. code-block:: none

   |-----------------------------+------------------------------------------------------|
   | 2 * (char)                  | The characters "si"                                  |   
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
   |   For each energy level:    |                                                      |   
   |     NZ * (double complex)   | Density matrix element  K=Q=0 for each height node   |   
   |                             |                                                      |   
   | For each atom:              |                                                      |   
   | For each bound-bound        |                                                      |   
   |  radiative transition:      |                                                      |   
   |   NZ * (double)             | Mean radiation field for each height node integrated |   
   |                             |   over the absorption profile                        |   
   | If ST:                      |                                                      |   
   |   For each atom:            |                                                      |   
   |   For each bound-bound      |                                                      |   
   |    radiative transition:    |                                                      |   
   |     NZ * (double)           | Mean radiation field for each height node integrated |   
   |                             |   over the emission profile                          |   
   | For each atom:              |                                                      |   
   | For each bound-free         |                                                      |   
   |  radiative transition:      |                                                      |   
   |   NZ * (double)             | Integrated photoionization and recombination rates   |   
   |                             |   for the statistical equilibrium equations          |   
   |                             |                                                      |   
   | If (AA and not dynamic) or  |                                                      |   
   |   not partial frequency     |                                                      |   
   |   redistribution:           |                                                      |   
   |                             |                                                      |   
   |   For each height:          |                                                      |   
   |     NF * (double)           | Radiation field tensor for each frequency            |   
   |                             |                                                      |   
   | else:                       |                                                      |   
   |                             |                                                      |   
   |   For each height:          |                                                      |   
   |   For each polar direction: |                                                      |   
   |   For each azimuthal        |                                                      |   
   |     direction:              |                                                      |   
   |   For each frequency:       |                                                      |   
   |                             |                                                      |   
   |     1 * (double)            | Intensity Stokes parameter                           |   
   |                             |                                                      |   
   |------------------------------------------------------------------------------------|

