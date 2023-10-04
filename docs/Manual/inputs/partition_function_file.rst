.. code-block:: none
    
   |-------------------------------+------------------------------------------------------|
   |  1 * (int)                    | Number of tabulated temperatures [NT]                |
   | NT * (double)                 | Tabulated temperatures                               |
   |--------------------------------------------------------------------------------------|
   | For each of the 99 elements   |                                                      |   
   |                               |                                                      |   
   |   2 * (char)                  | Atomic symbol                                        |   
   |   1 * (char)                  | An non-used character                                |   
   |   1 * (int)                   | Number of stages for this atomic species [NS]        |   
   |  NS * NT * (double)           | Partition function for each stage (slow) and temp.   |   
   |  NS * (double)                | Ionization energy                                    |   
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|
