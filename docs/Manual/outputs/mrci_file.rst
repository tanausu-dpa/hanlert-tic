The **MRCI** file, for the 1D synthesis, has lines with six columns with the following information.

1. The iteration number.
2. The maximum relative change of the density matrix element with :math:`K=Q=0` among
   all the active atoms (see **ATOM_INPUT** in :ref:`control file <input_file>`).

   .. note::

     If :math:`\rho^0_0` and :math:`{\rho^0_0}^{\dagger}` are the values of the
     density matrix element with :math:`K=Q=0` for the last and second to last iterations,
     the maximum relative change (MRC) is calculated as:

     .. math::

        {\rm MRC} = \frac{\|\rho^0_0 - {\rho^0_0}^{\dagger}\|}{{\rho^0_0}^{\dagger}} & \quad
                       {\rm if} \,
                       {\rho^0_0}^{\dagger} > 10^{-15}\\
        {\rm MRC} = \frac{\|\rho^0_0 - {\rho^0_0}^{\dagger}\|}{{\rho^0_0}} & \quad
                       {\rm if} \,
                       {\rho^0_0} > 10^{-15}\\
        {\rm MRC} = \|\rho^0_0 - {\rho^0_0}^{\dagger}\| & \quad {\rm otherwise}


3. The index of the atom with the maximum relative change.
4. The index of the energy level within the atom with the maximum relative change.
5. The index of the height node with the maximum relative change.
6. The value of the height (in kilometers) or optical depth (depending on the
   :ref:`model atmosphere <atmo_file>`) with the maximum relative change.

::

           1    3.375551308513E+05          1           1           22   2200.840
           2    5.235407659015E-01          1           3           24   2199.000
           3    4.143517690516E-01          1           3           23   2200.100
           4    2.725264843674E-01          1           3           21   2201.160
           .     .                          .           .            .       .
           .     .                          .           .            .       .

In the 1.5D synthesis, the file is written in binary and can be read with *hanlertio_class*.
