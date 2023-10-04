The format of the **atmospheric** file is the following:

1. The four characters *2Dat*.
2. Integer with the byte size of the main variables (4 for floats or 8 for doubles).
3. Integer with the number of nodes along the *X* axis.
4. Integer with the number of nodes along the *Y* axis.
5. Integer with the number of nodes along the *Z* axis.
6. For each X-Y position, several double precision vectors, with the same
   dimension than the *Z* axis, with the following information (24 variables):

    * Geometrical height in km.
    * Optical depth scale.
    * Continuum opacity in cm:sup:`-1` at the reference wavelength.
    * Temperature in K.
    * Gas pressure in dyn cm\ :sup:`-1`.
    * Density in g cm\ :sup:`-3`.
    * Magnetic field component along the *X* axis in G.
    * Magnetic field component along the *Y* axis in G.
    * Magnetic field component along the *Z* axis in G.
    * Velocity along the *X* axis in km s\ :sup:`-1`.
    * Velocity along the *Y* axis in km s\ :sup:`-1`.
    * Velocity along the *Z* axis in km s\ :sup:`-1`.
    * Microturbulent velocity in km s\ :sup:`-1`.
    * Electron pressure in dyn cm\ :sup:`-1` or electron density
      in g cm\ :sup:`-1`.
    * Electron number density in cm\ :sup:`-3`.
    * Total hydrogen number density in cm\ :sup:`-3`.
    * Total atomic hydrogen number density in cm\ :sup:`-3`.
    * H\ :sup:`-` number density.
    * Number density of hydrogen in the ground state in cm\ :sup:`-3`
    * Number density of hydrogen in the first excited state in cm\ :sup:`-3`
    * Number density of hydrogen in the second excited state in cm\ :sup:`-3`
    * Number density of hydrogen in the third excited state in cm\ :sup:`-3`
    * Number density of hydrogen in the fourth excited state in cm\ :sup:`-3`
    * Proton Number density in cm\ :sup:`-3`

.. note::
   Not all columns work as input. What quantities are used as
   input among densities, number densities, and pressures is determined
   by the **ATMO_CHAR** keyword in the :ref:`control file <inputs_control>`.

.. code-block:: none
    
   |-------------------------------+------------------------------------------------------|
   | 4 * (char)                    | The characters "2Dat"                                |   
   | 1 * (int)                     | Byte size (bsize) of the variables (4 or 8)          |
   | 1 * (int)                     | Number of nodes along the *X* axis                   |
   | 1 * (int)                     | Number of nodes along the *Y* axis                   |
   | 1 * (int)                     | Number of nodes along the *Z* axis [NZ]              |
   |--------------------------------------------------------------------------------------|
   | For each X position:          |                                                      |   
   | For each Y position:          |                                                      |   
   |                               |                                                      |   
   |  NZ * (bsize)                 | Geometrical height scale                             |   
   |  NZ * (bsize)                 | Optical depth scale                                  |   
   |  NZ * (bsize)                 | Continuum opacity at the reference wavelength        |   
   |  NZ * (bsize)                 | Temperature                                          |   
   |  NZ * (bsize)                 | Gas pressure                                         |   
   |  NZ * (bsize)                 | Density                                              |   
   |  NZ * (bsize)                 | Magnetic field X component                           |   
   |  NZ * (bsize)                 | Magnetic field Y component                           |   
   |  NZ * (bsize)                 | Magnetic field Z component                           |   
   |  NZ * (bsize)                 | Velocity X component                                 |   
   |  NZ * (bsize)                 | Velocity Y component                                 |   
   |  NZ * (bsize)                 | Velocity Z component                                 |   
   |  NZ * (bsize)                 | Microturbulence                                      |   
   |  NZ * (bsize)                 | Electron pressure or density                         |   
   |  NZ * (bsize)                 | Electron number density                              |   
   |  NZ * (bsize)                 | Total hydrogen number density                        |   
   |  NZ * (bsize)                 | Atomic hydrogen number density                       |   
   |  NZ * (bsize)                 | H- number density                                    |   
   |  NZ * (bsize)                 | Hydrogen ground state number density                 |   
   |  NZ * (bsize)                 | Hydrogen first excitation state number density       |   
   |  NZ * (bsize)                 | Hydrogen second excitation state number density      |   
   |  NZ * (bsize)                 | Hydrogen third excitation state number density       |   
   |  NZ * (bsize)                 | Hydrogen fourth excitation state number density      |   
   |  NZ * (bsize)                 | Proton number density                                |   
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|

