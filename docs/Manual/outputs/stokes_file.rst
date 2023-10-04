This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Stokes** file is the following:

    1. The two characters *be*.
    2. Integer with the number of frequencies of the problem.
    3. Double precision number with the polar angle of the given direction.
    4. Double precision number with the azimuthal angle of the given direction.
    5. Double precision array with the frequency axis.
    6. Double precision array with the intensity Stokes *I* values for each frequency.
    7. Double precision array with the Stokes *Q* values for each frequency.
    8. Double precision array with the Stokes *U* values for each frequency.
    9. Double precision array with the Stokes *V* values for each frequency.

    .. code-block:: none
        
       |-------------------------------+------------------------------------------------------|
       | 2 * (char)                    | The characters "be"                                  |   
       | 1 * (int)                     | Number of frequency nodes [NF]                       |   
       | 1 * (double)                  | Polar angle                                          |   
       | 1 * (double)                  | Azimuthal angle                                      |   
       | NF * (double)                 | Frequency array                                      |
       |--------------------------------------------------------------------------------------|
       | NF * (double)                 | Stokes I intensity emerging for each frequency       |   
       | NF * (double)                 | Stokes Q polarization emerging for each frequency    |   
       | NF * (double)                 | Stokes U polarization emerging for each frequency    |   
       | NF * (double)                 | Stokes V polarization emerging for each frequency    |   
       |--------------------------------------------------------------------------------------|

