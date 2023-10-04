The partition function file contains the information of the ionization
energy and partition function of the atomic elements. It is used
for the equation of state and the chemical equilibrium. This file
is written in binary.

The format of the file is the following.

1. A 4 bytes integer with the number of temperatures for which the
   partition functions are specified.

2. The list of double precision numbers (8 bytes reals) with the
   temperatures for which the partition functions are
   specified. The logarithm of the partition function is linearly
   interpolated to the temperatures of the
   :ref:`model atmosphere <atmo_file>`.

3. 99 blocks of data for different atomic elements with the following
   information.

   1. Two bytes with the characters of the atomic element symbol (H,
      He, Li, ...).

   2. An empty byte.

   3. A 4 bytes integer with the number of ions present for the
      atomic species.

   4. For each ion, a vector of double precision numbers (8 bytes
      reals) with the values of the partition function for the
      previously indicated temperatures.

   5. For each ion, a double precision number (8 bytes
      real) with the values of the ionization energy in cm\ :sup:`-1`.
