This file is written in binary and can be read with *hanlertio_class*.

..
    The format of the **.pop** file is the following:

    1. The two characters *bp*.
    2. Integer with the number of height nodes of the problem.
    3. Integer with the number of energy levels.
    4. For each height, double precision array with the
       population of each energy level.

    .. code-block:: none

       |-----------------------------+------------------------------------------------------|
       | 2 * (char)                  | The characters "bp"                                  |   
       | 1 * (int)                   | Number of height nodes [NZ]                          |   
       | 1 * (int)                   | Number of energy levels [NL]                         |   
       |------------------------------------------------------------------------------------|
       | NZ * NL * (double)          | Atomic populations for each height (slow) and energy |   
       |                             |   level (fast)                                       |   
       |------------------------------------------------------------------------------------|

