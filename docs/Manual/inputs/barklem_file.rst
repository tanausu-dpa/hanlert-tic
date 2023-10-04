The format of the Barklem files is the following:

1. Integer with the number of effective quantum numbers for the
   first level (*s*, *p*, and *d* for **BARK_SP**, **BARK_PD**,
   and **BARK_DF** in the :ref:`control file <input_file>`,
   respectively).
2. Integer with the number of effective quantum numbers for the
   second level (*p*, *d*, and *f* for **BARK_SP**, **BARK_PD**,
   and **BARK_DF** in the :ref:`control file <input_file>`,
   respectively).
3. Double precision vector with the effective quantum number of
   the first level (see 1.).
4. Double precision vector with the effective quantum number of
   the second level (see 2.).
5. Table with the cross sections for each combination of first
   and second level's effective quantum numbers (with the first
   level in the slow dimension).
6. Table with the velocity parameter for each combination of first
   and second level's effective quantum numbers (with the first
   level in the slow dimension).

.. code-block:: none
    
   |-------------------------------+------------------------------------------------------|
   |  1 * (int)                    | Number of 1st level effective quantum numbers [N1]   |
   |  1 * (int)                    | Number of 2nd level effective quantum numbers [N2]   |
   |--------------------------------------------------------------------------------------|
   | N1 * (double)                 | Effective quantum numbers 1st level                  |
   | N2 * (double)                 | Effective quantum numbers 2nd level                  |
   | N1 * N2 * (double)            | Cross section                                        |
   | N1 * N2 * (double)            | Velocity parameter                                   |
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|
