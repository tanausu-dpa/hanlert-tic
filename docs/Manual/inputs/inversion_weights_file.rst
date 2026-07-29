1. A 4 bytes integer with the number of wavelengths within the file.
2. A 4 bytes integer with the number of stokes parameters within the file.
3. A double precision array (8 bytes real numbers) with the weights for each Stokes parameter and wavelength.

.. code-block:: none
    
   |-------------------------------+------------------------------------------------------|
   | 1 * (int)                     | Number of wavelengths [NL]                           |   
   | 1 * (int)                     | Number of Stokes parameters [NS]                     |   
   |--------------------------------------------------------------------------------------|
   | NL * (double)                 | Stokes I intensity weights for each wavelength       |   
   | NL * (double) [if NS >= 2]    | Stokes Q polarization weights for each wavelength    |   
   | NL * (double) [if NS >= 3]    | Stokes U polarization weights for each wavelength    |   
   | NL * (double) [if NS == 4]    | Stokes V polarization weights for each wavelength    |   
   |--------------------------------------------------------------------------------------|

