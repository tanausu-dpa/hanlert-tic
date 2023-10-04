The **MRC** file, for the 1D synthesis, has lines with sixteen columns with the following information.

1. The iteration number.
2. The maximum relative change of the density matrix element with :math:`K=Q=0` among
   all the active atoms (see **ATOM_INPUT** in :ref:`control file <input_file>`).
   See note in :ref:`MRCI <mrc_file>` for details on the calculation of this
   value.
3. The index of the atom with the maximum relative change for the :math:`K=Q=0` density
   matrix element.
4. The index of the atomic term within the atom with the maximum relative change
   for the :math:`K=Q=0` density matrix element.
5. The total angular momentum of the energy level within the term with the maximum
   relative change for the :math:`K=Q=0` density matrix element, multiplied by 2.
6. The index of the height node with the maximum relative change for the :math:`K=Q=0` density
   matrix element.
7. The value of the height (in kilometers) or optical depth (depending on the
   :ref:`model atmosphere <atmo_file>`) with the maximum relative change for the :math:`K=Q=0`
   density matrix element.
8. The maximum relative change of the density matrix element with :math:`K\neq0` among
   all the active atoms (see **ATOM_INPUT** in :ref:`control file <input_file>`).

   .. note::

     If :math:`\rho^K_Q(J,J')` and :math:`{\rho^K_Q}^{\dagger}(J,J')` are the values
     of the density matrix element with :math:`K\neq0` for the last and second to
     last iterations, the maximum relative change (MRC) is calculated as:

     .. math::

        {\rm MRC} = \frac{\|\rho^K_Q(J,J') - {\rho^K_Q}^{\dagger}(J,J')\|}
                         {\|{\rho^K_Q}^{\dagger}(J,J')\|} 
                    & \quad {\rm if} \, \mathbb{F}[{\rho^K_Q}^{\dagger}(J,J')] > 10^{-15} \\
        {\rm MRC} = \frac{\|\rho^K_Q(J,J') - {\rho^K_Q}^{\dagger}(J,J')\|}
                           {\|{\rho^K_Q}(J,J')\|}
                    & \quad {\rm if} \, \mathbb{F}[{\rho^K_Q}(J,J')] > 10^{-15} \\
        {\rm MRC} = \|\rho^K_Q(J,J') - {\rho^K_Q}^{\dagger}(J,J')\|  & \quad {\rm otherwise}

     where

     .. math::

        \mathbb{F}[\rho^K_Q(J,J')] = \frac{\rho^K_Q(J,J')}
                                            {\sqrt{\rho^0_0(J,J)\rho^0_0(J',J')}}


9. The index of the atom with the maximum relative change for the :math:`K\neq0`
   density matrix element.
10. The index of the atomic term within the atom with the maximum relative change
    for the :math:`K\neq0` density matrix element.
11. The two total angular momentum of the energy levels within the term with the
    maximum relative change for the :math:`K\neq0` density matrix element,
    multiplied by 2.
12. The index of the height node with the maximum relative change for the
    :math:`K\neq0` density matrix element.
13. The value of the height (in kilometers) or optical depth (depending on the
    :ref:`model atmosphere <atmo_file>`) with the maximum relative change for the
    :math:`K\neq0` density matrix element.
14. The :math:`K` multipole index with the maximum relative change for the :math:`K\neq0`
15. The :math:`Q` multipole index with the maximum relative change for the :math:`K\neq0`
    density matrix element.

::

          1    6.186934732854E-02          1          3     5            33   2017.000      1.000000000000E+00          1          2     3     3            1   2219.430  2  0
          2    1.271413634840E-02          1          3     3            39   1580.000      8.684012222462E-01          1          3     5     5           43   1180.000  2  0
          3    7.121413928096E-03          1          3     3            41   1378.000      4.680716802297E-01          1          3     5     5           43   1180.000  2  0
          4    4.754039847017E-03          1          3     3            41   1378.000      1.247613802823E+00          1          3     5     5           43   1180.000  2  0
          .     .                          .          .     .             .       .          .                          .          .     .     .            .       .     .  .
          .     .                          .          .     .             .       .          .                          .          .     .     .            .       .     .  .


In the 1.5D synthesis, the file is written in binary and can be read with *hanlertio_class*.
