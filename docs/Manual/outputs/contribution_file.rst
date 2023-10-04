This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Contribution** file is the following:

    1. The two characters *bc*.
    2. Integer with the number of frequencies of the problem.
    3. Integer with the number of heights of the problem.
    4. Double precision number with the polar angle of the given direction.
    5. Double precision number with the azimuthal angle of the given direction.
    6. Double precision array with the frequency axis.
    7. Double precision array with the height or optical depth axis.
    8. Double precision array with the contribution function values of
       the intensity Stokes *I* for each height and frequency (slow and fast axes,
       respectively).
    9. Double precision array with the contribution function values of
       the Stokes *Q* for each height and frequency (slow and fast axes,
       respectively).
    10. Double precision array with the contribution function values of
        the Stokes *U* for each height and frequency (slow and fast axes,
        respectively).
    11. Double precision array with the contribution function values of
        the Stokes *V* for each height and frequency (slow and fast axes,
        respectively).

    .. code-block:: none
        
       |-------------------------------+------------------------------------------------------|
       | 2 * (char)                    | The characters "bc"                                  |   
       | 1 * (int)                     | Number of frequency nodes [NF]                       |   
       | 1 * (int)                     | Number of height nodes [NZ]                          |   
       | 1 * (double)                  | Polar angle                                          |   
       | 1 * (double)                  | Azimuthal angle                                      |   
       | NF * (double)                 | Frequency array                                      |
       | NZ * (double)                 | Height or optical depth array                        |
       |--------------------------------------------------------------------------------------|
       | NZ * NF * (double)            | Stokes I intensity contribution function for each    |   
       |                               |   height and frequency (slow, fast)                  |   
       | NZ * NF * (double)            | Stokes Q contribution function for each              |   
       |                               |   height and frequency (slow, fast)                  |   
       | NZ * NF * (double)            | Stokes U contribution function for each              |   
       |                               |   height and frequency (slow, fast)                  |   
       | NZ * NF * (double)            | Stokes V contribution function for each              |   
       |                               |   height and frequency (slow, fast)                  |   
       |--------------------------------------------------------------------------------------|

