This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Stokesout** file is the following:

    1. The two characters *bo*.
    2. Integer with the number of frequencies of the problem.
    3. Double precision array with the frequency axis of the problem.
    4. Integer with half the number of polar quadrature nodes directed
       toward the top boundary.
    5. Integer with the number of azimuthal quadrature nodes of the problem.
    6. For each polar quadrature direction going upwards, and for each azimuthal
       direction:

        * Two double precision numbers with the polar and azimuthal angles
          of the given direction.
        * Double precision array with the
          intensity Stokes *I* values for each frequency.
        * Double precision array with the
          Stokes *Q* values for each frequency.
        * Double precision array with the
          Stokes *U* values for each frequency.
        * Double precision array with the
          Stokes *V* values for each frequency.

    .. code-block:: none
        
       |-------------------------------+------------------------------------------------------|
       | 2 * (char)                    | The characters "bo"                                  |   
       | 1 * (int)                     | Number of frequency nodes [NF]                       |   
       | NF * (double)                 | Frequency array                                      |
       | 1 * (int)                     | Number of polar quadrature nodes in the outwards     |   
       |                               |   direction [NT]                                     |   
       | 1 * (int)                     | Number of azimuthal quadrature nodes [NP]            |   
       |--------------------------------------------------------------------------------------|
       | For each polar direction:     |                                                      |   
       | For each azimuthal direction: |                                                      |   
       |                               |                                                      |   
       |   1 * (double)                | Polar angle in degrees                               |   
       |   1 * (double)                | Azimuthal angle in degrees                           |   
       |   NF * (double)               | Stokes I intensity at the top boundary for each      |   
       |                               |   frequency                                          |   
       |   NF * (double)               | Stokes Q polarization at the top boundary for each   |   
       |                               |   frequency                                          |   
       |   NF * (double)               | Stokes U polarization at the top boundary for each   |   
       |                               |   frequency                                          |   
       |   NF * (double)               | Stokes V polarization at the top boundary for each   |   
       |                               |   frequency                                          |   
       |                               |                                                      |   
       |--------------------------------------------------------------------------------------|

