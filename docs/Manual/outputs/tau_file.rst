This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **Tau** file is the following:

    1. The two characters *bt*.
    2. Integer with the number of frequencies of the problem.
    3. Integer with the number of heights of the problem.
    4. Double precision number with the polar angle of the given direction.
    5. Double precision number with the azimuthal angle of the given direction.
    6. Double precision array with the frequency axis.
    7. Double precision array with the height or reference optical depth
       where the optical depth at each frequency is equal to one.

    .. code-block:: none
        
       |-------------------------------+------------------------------------------------------|
       | 2 * (char)                    | The characters "bt"                                  |   
       | 1 * (int)                     | Number of frequency nodes [NF]                       |   
       | 1 * (int)                     | Number of height nodes [NZ]                          |   
       | 1 * (double)                  | Polar angle                                          |   
       | 1 * (double)                  | Azimuthal angle                                      |   
       | NF * (double)                 | Frequency array                                      |
       |--------------------------------------------------------------------------------------|
       | NF * (double)                 | Height or reference optical depth where the optical  |
       |                               |  depth is one for each frequency                     |
       |--------------------------------------------------------------------------------------|

