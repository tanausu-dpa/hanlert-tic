The format of the **atmospheric** file is the following:

1. The three characters *CLE*.
2. Integer with the byte size of the main variables (4 for floats or 8 for doubles).
3. Integer indicating the type of atmospheric model.

    * 0 for lines of sight in cartesian grid.
    * 1 for a list of slabs.
    * 2 for lines of sight in arbitrary positions.

4. An integer indicating if the spatial coordinates are normalized to the stellar radius (negative means they are not normalized and in cm.

**Mode 0 (LOSs in cartesian grid)**

5. Integer with the number of nodes along the *X* axis.
6. Integer with the number of nodes along the *Y* axis.
7. Integer with the number of nodes along the *Z* axis.
8. *X* axis.
9. *Y* axis.
10. *Z* axis.
11. For each Y-Z position, several vectors, with the same
    dimension than the *X* axis, with the following information (22 variables):

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
   | 3 * (char)                    | The characters "CLE"                                 |   
   | 1 * (int)                     | Byte size (bsize) of the variables (4 or 8)          |
   | 1 * (int)                     | Integer with value 0 (indicating LOSs in cartesian)  |
   | 1 * (int)                     | Integer indicating if positions are normalized       |
   | 1 * (int)                     | Number of nodes along the *X* axis [NX]              |
   | 1 * (int)                     | Number of nodes along the *Y* axis [NY]              |
   | 1 * (int)                     | Number of nodes along the *Z* axis [NZ]              |
   | NX * (bsize)                  | Positions along the *X* axis                         |
   | NY * (bsize)                  | Positions along the *Y* axis                         |
   | NZ * (bsize)                  | Positions along the *Z* axis                         |
   |--------------------------------------------------------------------------------------|
   | For each Y position:          |                                                      |   
   | For each Z position:          |                                                      |   
   |                               |                                                      |   
   |  NX * (bsize)                 | Continuum opacity at the reference wavelength        |   
   |  NX * (bsize)                 | Temperature                                          |   
   |  NX * (bsize)                 | Gas pressure                                         |   
   |  NX * (bsize)                 | Density                                              |   
   |  NX * (bsize)                 | Magnetic field X component                           |   
   |  NX * (bsize)                 | Magnetic field Y component                           |   
   |  NX * (bsize)                 | Magnetic field Z component                           |   
   |  NX * (bsize)                 | Velocity X component                                 |   
   |  NX * (bsize)                 | Velocity Y component                                 |   
   |  NX * (bsize)                 | Velocity Z component                                 |   
   |  NX * (bsize)                 | Microturbulence                                      |   
   |  NX * (bsize)                 | Electron pressure or density                         |   
   |  NX * (bsize)                 | Electron number density                              |   
   |  NX * (bsize)                 | Total hydrogen number density                        |   
   |  NX * (bsize)                 | Atomic hydrogen number density                       |   
   |  NX * (bsize)                 | H- number density                                    |   
   |  NX * (bsize)                 | Hydrogen ground state number density                 |   
   |  NX * (bsize)                 | Hydrogen first excitation state number density       |   
   |  NX * (bsize)                 | Hydrogen second excitation state number density      |   
   |  NX * (bsize)                 | Hydrogen third excitation state number density       |   
   |  NX * (bsize)                 | Hydrogen fourth excitation state number density      |   
   |  NX * (bsize)                 | Proton number density                                |   
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|

**Mode 1 (list of slabs)**

5. Integer with the number of slabs.
6. For each slab, several numbers, with the following information (25 variables):

    * Height over the stellar surface.
    * Cosine of the angle between the slab axis and the line of sight.
    * Angle between the slab's reference azimuth and the line of sight.
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
   | 3 * (char)                    | The characters "CLE"                                 |   
   | 1 * (int)                     | Byte size (bsize) of the variables (4 or 8)          |
   | 1 * (int)                     | Integer with value 0 (indicating LOSs in cartesian)  |
   | 1 * (int)                     | Integer indicating if positions are normalized       |
   | 1 * (int)                     | Number of slabs                                      |
   |--------------------------------------------------------------------------------------|
   | For each slab:                |                                                      |   
   |                               |                                                      |   
   |  1 * (bsize)                  | Height over the stellar surface                      |
   |  1 * (bsize)                  | Cosine of the angle between the slab axis and LOS    |
   |  1 * (bsize)                  | Angle between the slab's reference azimuth and LOS   |
   |  1 * (bsize)                  | Continuum opacity at the reference wavelength        |   
   |  1 * (bsize)                  | Temperature                                          |   
   |  1 * (bsize)                  | Gas pressure                                         |   
   |  1 * (bsize)                  | Density                                              |   
   |  1 * (bsize)                  | Magnetic field X component                           |   
   |  1 * (bsize)                  | Magnetic field Y component                           |   
   |  1 * (bsize)                  | Magnetic field Z component                           |   
   |  1 * (bsize)                  | Velocity X component                                 |   
   |  1 * (bsize)                  | Velocity Y component                                 |   
   |  1 * (bsize)                  | Velocity Z component                                 |   
   |  1 * (bsize)                  | Microturbulence                                      |   
   |  1 * (bsize)                  | Electron pressure or density                         |   
   |  1 * (bsize)                  | Electron number density                              |   
   |  1 * (bsize)                  | Total hydrogen number density                        |   
   |  1 * (bsize)                  | Atomic hydrogen number density                       |   
   |  1 * (bsize)                  | H- number density                                    |   
   |  1 * (bsize)                  | Hydrogen ground state number density                 |   
   |  1 * (bsize)                  | Hydrogen first excitation state number density       |   
   |  1 * (bsize)                  | Hydrogen second excitation state number density      |   
   |  1 * (bsize)                  | Hydrogen third excitation state number density       |   
   |  1 * (bsize)                  | Hydrogen fourth excitation state number density      |   
   |  1 * (bsize)                  | Proton number density                                |   
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|

**Mode 2 (list of LOSs)**

5. Integer with the number of LOSs.
6. For each LOS:

    * Y coordinate.
    * Z coordinate.
    * Integer with the number of nodes along the *X* axis.
    * Positions along the *X* axis.
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
   | 3 * (char)                    | The characters "CLE"                                 |   
   | 1 * (int)                     | Byte size (bsize) of the variables (4 or 8)          |
   | 1 * (int)                     | Integer with value 0 (indicating LOSs in cartesian)  |
   | 1 * (int)                     | Integer indicating if positions are normalized       |
   | 1 * (int)                     | Number of lines of sight                             |
   | 1 * (int)                     | Number of nodes along the *Y* axis [NY]              |
   |--------------------------------------------------------------------------------------|
   | For each LOS:                 |                                                      |   
   |                               |                                                      |   
   |   1 * (bsize)                 | Y position                                           |
   |   1 * (bsize)                 | Z position                                           |
   |   1 * (int)                   | Number of nodes along the *X* axis [NX]              |
   |  NX * (bsize)                 | Continuum opacity at the reference wavelength        |   
   |  NX * (bsize)                 | Temperature                                          |   
   |  NX * (bsize)                 | Gas pressure                                         |   
   |  NX * (bsize)                 | Density                                              |   
   |  NX * (bsize)                 | Magnetic field X component                           |   
   |  NX * (bsize)                 | Magnetic field Y component                           |   
   |  NX * (bsize)                 | Magnetic field Z component                           |   
   |  NX * (bsize)                 | Velocity X component                                 |   
   |  NX * (bsize)                 | Velocity Y component                                 |   
   |  NX * (bsize)                 | Velocity Z component                                 |   
   |  NX * (bsize)                 | Microturbulence                                      |   
   |  NX * (bsize)                 | Electron pressure or density                         |   
   |  NX * (bsize)                 | Electron number density                              |   
   |  NX * (bsize)                 | Total hydrogen number density                        |   
   |  NX * (bsize)                 | Atomic hydrogen number density                       |   
   |  NX * (bsize)                 | H- number density                                    |   
   |  NX * (bsize)                 | Hydrogen ground state number density                 |   
   |  NX * (bsize)                 | Hydrogen first excitation state number density       |   
   |  NX * (bsize)                 | Hydrogen second excitation state number density      |   
   |  NX * (bsize)                 | Hydrogen third excitation state number density       |   
   |  NX * (bsize)                 | Hydrogen fourth excitation state number density      |   
   |  NX * (bsize)                 | Proton number density                                |   
   |                               |                                                      |   
   |--------------------------------------------------------------------------------------|

