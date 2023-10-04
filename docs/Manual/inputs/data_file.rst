Fits
----

Not documented.

Binary
------

The format of the binary data file is the following:

1. The four characters *invi*.
2. Integer with the number of nodes along the *X* axis.
3. Integer with the number of nodes along the *Y* axis.
4. Integer with the number of wavelengths in the spectra.
5. An integer indicating if the file contains only intensity (0) or full Stokes (1).
6. An integer indicating if all the lines of sight are the same (0) or change pixel by pixel (1).
7. An integer indicating if the type of error for the Stokes parameters:

   * 0 if none specified.
   * 1 if a constant :math:`\sigma` value for all wavelengths and pixels.
   * 2 if a wavelength dependent :math:`\sigma` value common to all pixels.
   * 3 if a constant :math:`\sigma` value changing pixel by pixel.
   * 4 if a wavelength dependent :math:`\sigma` value different for each pixel.

8. An integer indicating if the type of diffuse light profile included in the data file.

   * 0 if none specified.
   * 1 if an only intensity profile common for all pixels.
   * 2 if a full Stokes profiles common for all pixels.
   * 3 if an only intensity profile changing for each pixel.
   * 4 if a full Stokes profiles changing for each pixel.

9. Wavelength axis of the spectra in nm.

10. Only if the line of sight is common for all pixels, its heliocentric and azimuthal angles, in radians.

11. Only if constant for all pixels, the :math:`\sigma` value (or profile) for the relevant Stokes parameters depending on the data file, in SI units.

12. Only if constant for all pixels, the diffuse light profile/s in SI units.

13. For each X-Y position:

    * If the line of sight changes for each pixel, its heliocentric and azimuthal angles, in radians.
    * Observed Stokes parameters in SI units.
    * If changing for each pixel, the :math:`\sigma` value (or profile) for the relevant Stokes parameters in SI units.
    * If changing for each pixel, the diffuse light profile/s in SI units.

.. code-block:: none
    
   |-------------------------------+------------------------------------------------------|
   | 4 * (char)                    | The characters "invi"                                |   
   | 1 * (int)                     | Number of nodes along the *X* axis                   |
   | 1 * (int)                     | Number of nodes along the *Y* axis                   |
   | 1 * (int)                     | Number of nodes along the wavelength axis [NW]       |
   | 1 * (int)                     | If polarization [P]                                  |
   |                               |    - 0 only intensity                                |
   |                               |    - 1 full Stokes                                   |
   | 1 * (int)                     | Type of line of sight [L]                            |
   |                               |    - 0 common for the full field of view             |
   |                               |    - 1 changing pixel by pixel                       |
   | 1 * (int)                     | Type of line of :math:`\sigma` [S]                   |
   |                               |    - 0 none                                          |
   |                               |    - 1 constant and common for all pixels            |
   |                               |    - 2 wavelength dependent common for all pixels    |
   |                               |    - 3 constant pixel wise                           |
   |                               |    - 4 wavelength dependent pixel wise               |
   | 1 * (int)                     | Type of diffuse light [D]                            |
   |                               |    - 0 none                                          |
   |                               |    - 1 only intensity common for all pixels          |
   |                               |    - 2 full Stokes common for all pixels             |
   |                               |    - 3 only intensity pixel wise                     |
   |                               |    - 4 full Stokes pixel wise                        |
   | NL * (double)                 | Wavelength axis in nm                                |
   | if L == 0:                    |                                                      |
   |   1 * (double)                | Heliocentric angle for the line of sight             |
   |   1 * (double)                | Azimuthal angle for the line of sight                |
   | if S == 1:                    |                                                      |
   |   1 * (double)                | :math:`\sigma` value for the intensity (SI units)    |
   | if S == 1 and P == 1:         |                                                      |
   |   3 * (double)                | :math:`\sigma` value for Q, U, and V (SI units)      |
   | if S == 2:                    |                                                      |
   |   NL * (double)               | :math:`\sigma` profile for the intensity (SI units)  |
   | if S == 2 and P == 1:         |                                                      |
   |   3 * NL * (double)           | :math:`\sigma` profile for Q, U, V (SI units)        |
   | if D == 1 v D == 2:           |                                                      |
   |   NL * (double)               | Intensity diffuse light profile (SI units)           |
   | if D == 2:                    |                                                      |
   |   3 * NL * (double)           | Q, U, and V diffuse light profiles (SI units)        |
   |--------------------------------------------------------------------------------------|
   | For each X position:          |                                                      |   
   | For each Y position:          |                                                      |   
   |                               |                                                      |   
   |   if L == 1:                  |                                                      |   
   |     1 * (double)              | Heliocentric angle for the line of sight             |
   |     1 * (double)              | Azimuthal angle for the line of sight                |
   |   NL * (double)               | Intensity profile (SI units)                         |
   |   if P == 1:                  |                                                      |   
   |     3 * NL * (double)         | Q, U, and V Stokes profiles (SI units)               |
   |   if S == 3:                  |                                                      |
   |     1 * (double)              | :math:`\sigma` value for the intensity (SI units)    |
   |   if S == 3 and P == 1:       |                                                      |
   |     3 * (double)              | :math:`\sigma` value for Q, U, and V (SI units)      |
   |   if S == 4:                  |                                                      |
   |     NL * (double)             | :math:`\sigma` profile for the intensity (SI units)  |
   |   if S == 4 and P == 1:       |                                                      |
   |     3 * NL * (double)         | :math:`\sigma` profile for Q, U, V (SI units)        |
   |   if D == 3 v D == 4:         |                                                      |
   |     NL * (double)             | Intensity diffuse light profile (SI units)           |
   |   if D == 4:                  |                                                      |
   |     3 * NL * (double)         | Q, U, and V diffuse light profiles (SI units)        |
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|
