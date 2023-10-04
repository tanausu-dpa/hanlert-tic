The background fudge file consists of a series of lines with four
floating points numbers indicating:

1. The wavelength to which the fudge factors correspond in nanometers.

2. The factor that should multiply the H\ :sup:`-` opacity at the
   given wavelength.

3. The factor that should multiply the scattering coefficient :math:`\sigma`
   at the given wavelength.

4. The factor that should multiply the contribution of the bound-free
   transitions of metals (every atom but Hydrogen) to the opacity
   at the given wavelength.

::

     100.000         1.00          1.00          1.00
     152.499         1.00          1.00          1.00
     152.500         1.00          1.00          2.00
        .             .             .             .
        .             .             .             .
        .             .             .             .

.. note::

   The fudge factors are linearly interpolated into the frequency
   axis of the radiative transfer problem. Therefore, it is
   recommended to have a line with factors set to 1.0 before and
   after the lines with fudge factors not unity.
