1. The first line contains the name of the model. This field is not used
   by the code.

   ::

      FALC93_70

2. The second line indicates the type of scale of the stratification. Accepted
   string are:

   * *HEIGHT SCALE*. The stratification is to be expressed in geometrical
     height in kilometers.
   * *TAU SCALE*. The stratification is to be expressed in optical depth
     for the reference wavelength.

   ::

       HEIGHT SCALE

   .. note::

      The default reference wavelength is 500 nm. However, this reference
      can be changed by just adding a new wavelength (in nanometers) in the
      same line where the scale string is specified, e.g., *TAU SCALE 600*.
      Note that when specifying *HEIGHT SCALE*, the reference wavelength
      will be the one for which the optical depths will be calculated
      if **KEEP_ATMO** = *Yes* in the :ref:`control file <input_file>`.
   
3. The third line is a floating point number with the logarithm of the
   gravity constant. This field is **only** used during the inversion
   when assuming hydrostatic equilibrium.

   ::

       4.44

4. The fourth line is an integer with the number of height nodes of the
   model atmosphere.

   ::

      70

5. Then we find a block with a number of columns between five and seven
   (the last two are optional) with the following data:

   1. Height in kilometers (if *HEIGHT SCALE*) or optical depth at the
      reference wavelength (if *TAU SCALE*) for each atmospheric node,
      ordered from the top to the bottom boundary (decreasing height
      or increasing optical depth).
   2. The temperature in kelvin for each atmospheric node.
   3. The electron number density in cm\ :sup:`-3` for each atmospheric
      node.
   4. The vertical velocity (v\ :sub:`z`) in kilometers per second for
      each atmospheric node. Positive values indicate upward flow.
   5. The microturbulent velocity in kilometers per second for each
      atmospheric node.
   6. The horizontal velocity component in the direction of zero
      azimuth (v\ :sub:`x`) in kilometers per second for each
      atmospheric node. This column is optional.
   7. The horizontal velocity compoenent in the direction of
      :math:`90^{\circ}` azimuth (v\ :sub:`y`) in kilometers per second
      for each atmospheric node. This column is optional. Note that
      in order to specify this quantity, the sixth columns stops being
      optional and should be present, even if just with zero values.

   ::

        2219.430       102770.00       6.560E+09       0       11.900
        2217.880        98790.00       6.810E+09       0       11.800
        2216.430        94800.00       7.070E+09       0       11.700
        2215.080        90816.00       7.360E+09       0       11.600
            .                .          .              .         .
            .                .          .              .         .
            .                .          .              .         .

6. Then we find a block with between two and six columns.

   1. From the first to the second to last column, we have the Hydrogen
      number density for the first, second, etc., atomic energy levels
      (specified only by the principal quantum number) in cm\ :sup:`-3`
      for each atmospheric node.
   2. The last column is the proton number density in cm\ :sup:`-3` for
      each atmospheric node.

   ::

       1.110E+05       3.700E-02       0       0       0       5.480E+09
       1.790E+05       5.900E-02       0       0       0       5.700E+09
       2.460E+05       8.000E-02       0       0       0       5.940E+09
       3.440E+05       1.110E-01       0       0       0       6.210E+09
        .               .              .       .       .        .
        .               .              .       .       .        .
        .               .              .       .       .        .

There are other options for the model atmosphere that do not need the
specification of the electron and Hydrogen number densities. In order
to use one of these alternative inputs, the Hydrogen number densities
(the block in point 6. above) must be substituted by just one string
which will indicate the nature of the variable in the third column
of the first block of data (subpoint 3. in point 5. above). The
available strings are:

* *ne*. The input is the electron number density in cm\ :sup:`-3`. The
  difference with the standard format is that the Hydrogen number
  density is not required.
* *pe*. The input is the electron pressure in :math:`{\rm dyn}\cdot{\rm cm}^{-2}`.
* *rhoe*. The input is the electron density in :math:`{\rm g}\cdot{\rm cm}^{-3}`.
* *pg*. The input is the gas pressure in :math:`{\rm dyn}\cdot{\rm cm}^{-2}`.
* *rho*. The input is the density in :math:`{\rm g}\cdot{\rm cm}^{-3}`.

.. note::

   For atmospheric file formats other than the standard (electron and
   Hydrogen number densities), the code uses the equation of state in
   order to compute these number densities from the provided input. For
   more details see :ref:`equation of state <eq_state>`

