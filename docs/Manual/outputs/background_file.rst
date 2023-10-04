This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **background** file is the following:

    1. The two characters *ba*.
    2. Integer with the number of frequencies.
    3. Double precision array with the frequency axis.
    4. Integer with the number of height nodes.
    5. Integer with the number of directions with different
       continuum contributions.
    6. For each height and direction:

       * Double precision array with the absorptivity for each frequency.
       * Double precision array with the scattering coefficient for each frequency.
       * Double precision array with the emissivity for each frequency.

    .. code-block:: none
        
       |-------------------------------+------------------------------------------------------|
       | 2 * (char)                    | The characters "ba"                                  |   
       | 1 * (int)                     | Number of frequency nodes [NF]                       |   
       | NF * (double)                 | Frequency array                                      |
       | 1 * (int)                     | Number of height nodes [NZ]                          |   
       | 1 * (int)                     | Number of directions with different background       |   
       |                               |   radiative transfer coefficients [ND]               |   
       |--------------------------------------------------------------------------------------|
       |                               |                                                      |   
       | For each height:              |                                                      |   
       |   For each direction:         |                                                      |   
       |     NF * (double)             | Absorptivity for each frequency                      |   
       |     NF * (double)             | Scattering coefficient for each frequency            |   
       |     NF * (double)             | Emissivity for each frequency                        |   
       |                               |                                                      |   
       |--------------------------------------------------------------------------------------|

