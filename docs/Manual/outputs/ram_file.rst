The real name of the **RAM** files is **RAM-X-Y.dat**, with **X** being **I**
or **P** for the only intensity and polarization problems, respectively, and
**Y** being **F** or **E**, for the iterative problem and the calculation of
the emergent Stokes parameters, respectively. If the keyword
**RAM_REPORT** is set to *Yes* in the :ref:`control file <input_file>`, one
**RAM** file will be created just before each of the relevant steps that
correspond to the run (iteration of the only intensity problem, calculation
of the emergent intensity, iteration of the polarized problem, calculation
of the emergent Stokes parameters).

Each file shows, for each process in increasing rank, how many MB are being
used to store each of the above mentioned quantities. The sum of them
(shown in the column *Total*) is the quantity which is limited by the
**RAM_LIMIT** keyword in the :ref:`control file <input_file>`.

::

    --------------------------------------------------------------------------------------------
    #  CPU | Background  | Photoioni.  |  Radiation  | Voigt prf.  | Redistrib.  |    Total    |
    --------------------------------------------------------------------------------------------
         0 |       0.000 |       0.182 |       2.191 |       0.000 |       0.000 |       2.373 |
         1 |       0.996 |       0.000 |       2.191 |       0.526 |     138.453 |     142.167 |
         2 |       0.250 |       0.000 |       2.191 |       0.176 |     146.017 |     148.634 |
         3 |       0.149 |       0.000 |       2.191 |       0.116 |      96.332 |      98.788 |
         4 |       0.782 |       0.184 |       2.191 |       0.114 |      96.316 |      99.588 |
    --------------------------------------------------------------------------------------------

