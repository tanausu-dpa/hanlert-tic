The first block of the file contains general information about the atom.

.. note::
   This file admits comments, any text following the symbol *"!"* is ignored.

::

    MG           ! Atomic label
    24.305       ! Atomic mass (AMU)
    7.58  1  0   ! Abundance  %isotope  [renorm]
    ML           ! Type of model

#. The first line indicates the symbol of the atom, magnesium in this example.
#. The second is the atomic mass in atomic mass units.
#. The third is the abundance of the element (MG in this example), the amount of such
   abundance that is due to this isotope, and if it is to be normalized.
#. The fourth is the type of model, "ML" for a multi-level atom and "MT" for a multi-term atom.

.. note::
   Different isotopes can be simultaneously taken into account by specifying an atomic
   model for each of them (see **ATOM_INPUT** in :ref:`Control File <inputs_control>`).
   The abundance of the element (including all isotopes) must be specified in all the
   atomic files (be sure it is the same), and each isotope will represent an abundance
   given by the fraction following the abundance. If the third parameter is specified,
   all the fractions will be normalized to 1 to ensure the element abundance is the
   specified one.

For a multi-term atom, only the last line changes.

::

    MG           ! Atomic label
    24.305       ! Atomic mass (AMU)
    7.58  1  0   ! Abundance  %isotope  [renorm]
    MT           ! Type of model

The second block in the atomic file contains some dimensionallity information. This is
common for both multi-level and multi-term atoms.

::

    ! Terms   Lines   Photoionizations    Elastic Col.   
          6       5                  5              5

#. The first number indicates the number of atomic levels (terms) for a multi-level (multi-term)
   atomic model
#. The second number indicates the number of bound-bound radiation transitions entries.
#. The third number indicates the number of bound-free radiation transitions entries.
#. The fourth number indicates the number of depolarizing collisional rates.

For our multi-term example the numbers are slightly different.

::

    ! Terms   Lines   Photoionizations    Elastic Col.   
          4       2                  5              5

The third block specifies the energy structure of the atom, and it has different structure for
the multi-level and multi-term atoms. For the former,

::

    ! Energy levels
    !L     S ion.-stage   term
    !J     E              term level
    0.0	 0.5   2          1
    0.5	 0.0       -1     1     1

    1.0	 0.5   2          2
    0.5	 35669.31  -1     2     2

    1.0	 0.5   2          3
    1.5	 35760.88  -1     3     3

    2.0  0.5   2          4
    2.5  71490.19  -1     4      4

    2.0  0.5   2          5
    1.5  71491.06  -1     5      5

    0.0  0.0   3          6
    0.0  121267.64 -1     6      6

We have two lines for each atomic level with the following information:

#. In the first line, orbital angular momentum, spin angular momentum, stage of the atom (1 for neutral,
   2 for single ionized, etc.), and the level index.
#. In the second line, total angular momentum, energy in cm\ :sup:`-1`, the Landé factor,
   and the level index twice.

.. note::
   The level index is repeated three times for format compatibility with the multi-term format below.
   The indexes are not actually used in the HanleRT-TIC fortran code, but are there for sanity check
   of the atomic model during the reading phase.

.. note::
   The -1 value in the Landé factors implies that the proper Landé factor will be computed by
   HanleRT-TIC in LS coupling.

.. note::
   More than two ionization stages and multiple levels in each one of them can be specified, given
   that all of them have at least one connection to another level via at one bounb-bound or
   bound-free transition.

For the multi-term atom instead

::

    ! Energy levels
    !L     S ion.-stage  term
    !J     E             term level
    0.0	 0.5   2         1
    0.5	 0.0             1     1

    1.0	 0.5   2         2
    0.5	 35669.31        2     2
    1.5	 35760.88        2     3

    2.0  0.5   2         3
    2.5  71490.19        3      4
    1.5  71491.06        3      5

    0.0  0.0   3         4
    0.0  121267.64       4      6

Each block corresponds to an atomic term. The blocks have a header with term data, followed by a
line for each of the levels pertaining to the same term.

#. In the first line, orbital angular momentum, spin angular momentum, stage of the atom (1 for
   neutral, 2 for single ionized, etc.), and the term index.
#. For each level, a line with the total angular momentum, the level energy in cm\ :sup:`-1`,
   the term index (same than in the first line) and the level index.

.. note::
   The indexes are not actually used in the HanleRT-TIC fortran code, but are there for sanity check
   of the atomic model during the reading phase.

.. warning::
   In the multi-term model atom, every term must be complete, that is, given a term with orbital and
   spin angular momentum *L* and *S*, there must be (2 *L* + 1)(2 *S* + 1) levels with total angular
   momentum *J* between *|L - S|* and *L + S*.

The following block is the bound-bound radiative transitions. For the multi-level atom,

::

    ! Transition Lines
    ! Upper_l Lower_l    Aul       Type     Arg1   Arg2   Arg3   Arg4   Stark    NFT    NFC    DT    DC  PRD
          2         1  2.5700d8  unsold        1      0      1      0     1e0    251    101  2000    15    1
          3         1  2.6000d8  unsold        1      0      1      0     1e0    251    101  2000    15    1
          5         2  4.0100d8  unsold        1      0      1      0     1e0    125     75    20     5    0
          4         3  4.7900d8  unsold        1      0      1      0     1e0    125     75    20     5    0
          5         3  0.7980d8  unsold        1      0      1      0     1e0    125     75    20     5    0

There is one line for each bound-bound transition with the following information:

1. Two integers with the level index of the upper and lower levels of the transition.
2. A floating point number with the Einstein coefficient for spontaneous emission in s\ :sup:`-1`.
3. Five fields that characterize the Van der Waals broadening of the transition. Fire types are allowed:

  * *unsold*. Computes the broadening with the Unsold's approximation (`1960ZA.....49..231T <https://ui.adsabs.harvard.edu/abs/1960ZA.....49..231T/abstract>`_). This option only uses arguments 1 and 3, which are enhancement factors for the H and He contributions, respectively.
  * *barklem*. Computes the broadening with the Anstee & O'Mara formalistm (`1995MNRAS.276..859A <https://academic.oup.com/mnras/article/276/3/859/1034639>`_, `1997MNRAS.290..102B <https://academic.oup.com/mnras/article/290/1/102/1107836>`_, `1998MNRAS.296.1057B <https://academic.oup.com/mnras/article/296/4/1057/1064804>`_). This option uses the four arguments. Argument 1 and 3 are the type of orbital of the upper and lower levels, respectively (i.e., *s*, *p*, *d*, or *f*). Arguments 2 and 4 are the energy, in cm\ :sup:`-1`, of the ionizing level for the upper and lower levels of the transition to be used to compute the effective principal quantum number, respectively (usually the energy of the ground level of the next ionization stage). If the transition does not admit this option (due to not being a neutral atom, because the type of orbital is not in the available list, or because the principal effective quantum number is not in the data table) the broadening is computed using the *unsold* approximation without enhancements.
  * *kurucz*. The broadening is constant with temperature and the constant is given in decimal logarithm in Argument 1, the only one it uses. This constant is just multiplied by the number density of neutral hydrogen in the ground level.
  * *cross*. Computes the broadening with the Anstee & O'Mara formalistm. but the parameters are directly specified. This option uses the first two arguments. Argument 1 is the :math:`\sigma` parameter, given in cm\ :sup:`2` multiplied by 10\ :sup:`14`'. Argument 2 is the :math:`\alpha` parameter.
  * *param*. Computes the broadening in a parametric way. This option uses the four arguments. Arguments 1 and 2 determine parameters for the Hydrogen contribution, while Arguments 3 and 4 determine parameters for the Helium contribution. The broadening is computed in the following way,

  .. math::

      10^{-8}\cdot Arg1\cdot N_H\cdot ((1 + \frac{m_H}{m_{atom}})*T)^{Arg2} + \\
      10^{-9}\cdot Arg3\cdot N_{He}\cdot ((1 + \frac{m_{He}}{m_{atom}})*T)^{Arg4}

.. note::

    where N\ :sub:`H` is the ground level hydrogen number density, m\ :sub:`H` is
    the hydrogen mass, m\ :sub:`atom` is the atom mass (magnesium in our
    example), m\ :sub:`He` is the helium mass, and T is the temperature.

4. Parameter for the quadratic Stark broadening. If it is a possitive number,
   it is treated as an enhancement factor for the expression of C\ :sub:`4`
   in
   `1960ZA.....49..231T <https://ui.adsabs.harvard.edu/abs/1960ZA.....49..231T/abstract>`_
   . If negative, the quadratic Stark contribution is computed as 10\ :sup:`6` times
   the absolute value of the parameter times the electron density in cm\ :sup:`-3`. To specify a parameter in logarithmic form, add a letter *L* as the first character in this field (in this case, only put a minus sign if the value of the logarithm is negative).

5. Two integers that control the number of frequency nodes to be assigned to the
   transition. The first (NFT) indicates the total number of frequencies, while
   the second (NFC) indicates how many of the total must be considered in the core
   of the transition line.

6. Two floating point numbers that control the range of frequencies of the transition.
   The first (DT) is the total range from the line center in Doppler widths that the
   line should cover, while the second (DC) indicates the range from the line center in
   Doppler width that should be considered in the core of the transition line.

.. note::

   The distintion between core and not-core is due to the different frequency resolution
   in the two regimes. NFC frequencies will cover DC Doppler widths (to each side of the
   line center) with regular spacing between nodes. NFT-NFC frequencies will cover
   the range between DC and DT Doppler widths (to each side of the line center)
   with logarithmic step size.

7. Integer which flags a transition as a partial frequency redistribution transition (1)
   or a complete frequency redistribution transition (0).

For the case of the multi-term atom, the format is the same, but all quantities
refer to atomic terms instead of atomic levels. For our example,

::

    ! Transition Lines
    ! Upper_t Lower_t    Aul      Type  Arg1   Arg2   Arg3   Arg4   Stark    NFT    NFC    DT    DC  PRD
          2         1  2.5900d8 unsold     1      0      1      0     1e0    251    101  2000    15    1
          3         2  4.7972d8 unsold     1      0      1      0     1e0    125     75    20     5    0

For a multi-term transition, there is an optional last parameter which flags if a transition
should build the frequency axis by considering each of the fine structure transitions between
J levels (default behaviour with flag set to 1) or only considering the term-term transition
(flag set to 0).

.. note::

   This option has uses for the simulation of hyperfine structure transitions with
   zero spin momentum (that can be simulated with the transformation L→ J, S→ I, J→ F),
   as the hyperfine structure transitions between F levels usually overlap and it is
   enough to consider the nodes derived from the global transition.

The next block codes the information about the depolarizing elastic collisions. Depolarizing
elastic collisional rates are always specified in a level wise manner and thus it does not
change between the two types of atomic models. Therefore, for both our examples,

::

    ! Elastic collisions
    !level  number_of_entries
    !    K  type_of_entry nz
    !       entry
         1  1
         1  FIT
            2.9704  0.3600  1.0

         2  1
         1  FIT
            5.8448  0.3640  1.0

         3  3
         1  FIT
            5.2791  0.3440  1.0
         2  FIT
            6.1487  0.3600  1.0
         3  FIT
            5.8910  0.3660  1.0

         4  5
         1  FIT
            26.1186 0.3490  1.0
         2  FIT
            31.7442 0.3523  1.0
         3  FIT
            33.6754 0.3590  1.0
         4  FIT
            35.9476 0.3600  1.0
         5  FIT
            33.8079 0.3500  1.0

         5  3
         1  FIT
            30.3626 0.3540  1.0
         2  FIT
            34.7140 0.3600  1.0
         3  FIT
            33.9266 0.3580  1.0


There is an entry for each level in which a rate is specified. Each entry can have several lines:

1. The first line specifies the level index and the number of sub-entries to expect.
2. Each sub-entry has a header line with the K multipole for which the elastic collisonal rate
   is being specified and the type of input. There is currently only one possible type to
   determine the data format,

  * *Fit*. A single data line with three coefficients *"a  b  c"*. The elastic collisional rates
    are then computed with the expression

  .. math::

      a\cdot 10^{-9}\cdot\left(\frac{T}{5000}\right)^{b}c^{\frac{T}{5000}} N_{HI} [{\rm cm}^{-3}]

  with N\ :sub:`HI` is the number density of neutral hydrogen and T the temprature.


.. note::

   There use to be an alternative option to input these rates, reason why specifying
   **FIT** became necessary. Such an alternative is not available anymore, but the
   specifier has remained.

Next comes the photoionization. Because they are specified in a level wise manner, they are
common for the two types of model. For our example,

::

    ! Photoionization entries
    ! Upper_l lower_l       type   NF
    ! Entry
    !   MG II 2P6 4S 2SE
            6       1   EXPLICIT   31
       85.6   1.8590E-23
       84.9   1.8820E-23
       84.3   1.9050E-23
       83.7   1.9260E-23
       83.1   1.9470E-23
       82.5   1.9670E-23
       79.6   2.0560E-23
       76.9   2.1290E-23
       74.4   2.1870E-23
       72.0   2.2330E-23
       69.8   2.2680E-23
       67.7   2.2930E-23
       65.8   2.3110E-23
       63.9   2.3230E-23
       62.2   2.3280E-23
       60.5   2.3290E-23
       59.0   2.3250E-23
       57.5   2.3180E-23
       56.1   2.3090E-23
       54.7   2.2960E-23
       53.4   2.2820E-23
       52.2   2.2650E-23
       51.0   2.2480E-23
       49.9   2.2290E-23
       48.9   2.2080E-23
       47.8   2.1870E-23
       43.3   2.0720E-23
       39.5   1.9480E-23
       36.4   1.8210E-23
       33.7   1.6940E-23
       31.4   1.5720E-23

    !   MG II 2P6 4S 2PO
            6       2   EXPLICIT   32
      123.2   5.7890E-23
      121.8   5.5480E-23
      120.5   5.3250E-23
      119.3   5.1180E-23
      118.0   4.9270E-23
      116.8   4.7500E-23
      111.1   4.0440E-23
      106.0   3.5650E-23
      101.3   3.2410E-23
       97.0   3.0230E-23
       93.0   2.8800E-23
       89.3   2.7880E-23
       86.0   2.7320E-23
       82.8   2.7000E-23
       79.9   2.6860E-23
       77.2   2.6820E-23
       74.7   2.6860E-23
       72.3   2.6950E-23
       70.1   2.7060E-23
       68.0   2.7180E-23
       66.0   2.7300E-23
       64.2   2.7420E-23
       62.4   2.7520E-23
       60.8   2.7600E-23
       59.2   2.7660E-23
       57.7   2.7700E-23
       51.2   2.7490E-23
       46.0   2.6630E-23
       41.8   2.5230E-23
       38.3   2.3470E-23
       35.3   2.1530E-23
       32.8   1.9570E-23

    !   MG II 2P6 4S 2PO
            6       3   EXPLICIT   32
      123.2   5.7890E-23
      121.8   5.5480E-23
      120.5   5.3250E-23
      119.3   5.1180E-23
      118.0   4.9270E-23
      116.8   4.7500E-23
      111.1   4.0440E-23
      106.0   3.5650E-23
      101.3   3.2410E-23
       97.0   3.0230E-23
       93.0   2.8800E-23
       89.3   2.7880E-23
       86.0   2.7320E-23
       82.8   2.7000E-23
       79.9   2.6860E-23
       77.2   2.6820E-23
       74.7   2.6860E-23
       72.3   2.6950E-23
       70.1   2.7060E-23
       68.0   2.7180E-23
       66.0   2.7300E-23
       64.2   2.7420E-23
       62.4   2.7520E-23
       60.8   2.7600E-23
       59.2   2.7660E-23
       57.7   2.7700E-23
       51.2   2.7490E-23
       46.0   2.6630E-23
       41.8   2.5230E-23
       38.3   2.3470E-23
       35.3   2.1530E-23
       32.8   1.9570E-23

    !   MG II 2P6 3D 2DE
            6       4   HYDROGENIC 15
       28.0   5.2850E-22

    !   MG II 2P6 3D 2DE
            6       5   HYDROGENIC 15
       28.0   5.2850E-22


Every photoionization entry has a header line which specifies the upper and
lower level of the transition (which should have different ionization stages),
the type of photoionization input and the number of frequency nodes dedicated
to the bound-free transition. There are two types of inputs,

  * *Hydrogenic*. Only one line of data is required after the header, with the
    minimum wavelength to consider for the transition (in nanometers) and the
    photoionization cross section at the ionization frequency in m\ :sup:`2`.
    The cross section for each wavelength is computed as hydrogenic by HanleRT-TIC
    (`1970stat.book.....M <https://ui.adsabs.harvard.edu/abs/1970stat.book.....M/abstract>`_).
  * *Explicit*. As many entries as the number of frequencies specified in the
    header with the explicit value of the photoionization cross sections in
    m\ :sup:`2` for different wavelength (in nanometers).

  .. warning::

     Photoionization entries in *Explicit* format must be ordered in wavelength
     (the order direction does not matter).

Finally, there is a block with the inellastic collisional rates. This section is also
common for both types of models, although some options only have an impact on the
multi-term type. For our example,

::

    ! Collisions
    ! Flag_for_collision_reading  Type(level/term wise)  *Temperatures
    !   type_flag upper  lower  *rates
    COLS LEVEL ION 1000 3000 5000 10000 20000 30000
    BE    2   1  6.8635e-07  4.0730e-07  3.2342e-07  2.4293e-07  1.8917e-07  1.6766e-07
    BE    3   1  6.8908e-07  4.0572e-07  3.2342e-07  2.6537e-07  1.8917e-07  1.6816e-07
    BE    4   1  7.3684e-08  4.4380e-08  3.5190e-08  2.5746e-08  1.8714e-08  1.5446e-08
    BE    5   1  7.3684e-08  4.4511e-08  3.5394e-08  2.5890e-08  1.8765e-08  1.5446e-08
    BE    3   2  3.7797e-07  2.1271e-07  1.6842e-07  1.3247e-07  1.0237e-07  8.2835e-08
    BE    4   2  1.1053e-07  6.7226e-08  5.2277e-08  3.5815e-08  2.2884e-08  1.7190e-08
    BE    5   2  5.3353e-07  3.2182e-07  2.5752e-07  1.9547e-07  1.5561e-07  1.4200e-07
    BE    4   3  7.2775e-07  4.4117e-07  3.4987e-07  2.6322e-07  2.0443e-07  1.8186e-07
    BE    5   3  3.0565e-07  1.8553e-07  1.4585e-07  1.0356e-07  7.2465e-08  5.9292e-08
    BE    5   4  1.1257e-06  6.4206e-07  5.1260e-07  3.6893e-07  2.3647e-07  1.7314e-07
    COLS LEVEL ION 3000 5000 7000 10000 30000 100000
    FE    1   6  3.0832e-34  5.0489e-24  1.2771e-19  2.6987e-16  5.2644e-11  5.6345e-09
    FE    2   6  1.6573e-26  2.8963e-19  3.9018e-16  9.1410e-14  5.8252e-10  1.8826e-08
    FE    3   6  1.7402e-26  2.9882e-19  3.9954e-16  9.3076e-14  5.8795e-10  1.8944e-08
    FE    4   6  1.4182e-18  2.5692e-14  1.8206e-12  4.6847e-11  9.6119e-09  9.3326e-08
    FE    5   6  1.4188e-18  2.5699e-14  1.8209e-12  4.6852e-11  9.6123e-09  9.3327e-08

The format of this section consist on lines which start with the *COLS* string, which
determines the nature of every entry that follows them either until another *COLS*
string is found or the end of the file is reached. We can distinguish then between
header lines (the ones starting with the *COLS* string) and the data lines which follows
them:

1. The header lines contain the following information:

   * The *COLS* string which identifies the headers.
   * A string which indicates the scope of the collisional rates that follow. There are
     two types:

     * *Level*. The collisions which follow are expressed in a level wise way.
     * *Term*. The collisions which follow are expressed in a term wise way. Of course,
       this type of entry only makes sense for the multi-term model atom.

   * A string which determines how the interpolation to the temperatures of the model
     atmosphere is to be performed. Valid strings are *neutral*, *ion*, *sneutral*,
     *sion*, *lneutral*, *lion*, *linear*, *spline*, *none*.

     * Bound-bound and charge transfer transitions with a string containing *neutral*
       are scaled with :math:`\frac{1}{\sqrt{T}}` before interpolation.
     * Bound-bound and charge transfer transitions with a string containing *ion*
       are scaled with :math:`\sqrt{T}` before interpolation.
     * Bound-free transitions with either *neutral* or *ion* strings are scaled
       with :math:`\frac{e^{\frac{\Delta E}{T}}}{\sqrt{T}}` before interpolation,
       with :math:`\Delta E` the energy difference between the involved levels.
     * Transitions without *neutral* or *ion* (i.e., *linear*, *spline*, and
       *none*), are not scaled before interpolation.
     * The *none* and *spline* string are equivalent.
     * If the string starts with the letter *l*, the interpolation will always
       be linear. Otherwise, the interpolation will be performed with cubic
       splines.
      
       .. note::
         In the case that during the splines interpolation a negative value results,
         the code will switch to linear interpolation for that collisional rate.

   * A list of temperature values (in kelvin) which indicates the temperatures for which
     the collisional rates are given in the entries which follow.

2. The data entries which have a string indicating the type of collisions, the indexes of
   the involved levels, and the collisional rate for each temperature specified in the
   header line, following the same order, in s\ :sup:`-1` cm\ :sup:`3`. The types of
   transition are the following:

    * *BE*. Inelastic collision with electrons resulting in a bound-bound transition.
    * *FE*. Inelastic collision with electrons resulting in a bound-free transition.
    * *BP*. Inelastic collision with protons resulting in a bound-bound transition.
    * *BH*. Inelastic collision with Hydrogen atoms in the ground level resulting
      in a bound-bound transition.
    * *FH*. Inelastic collision with Hydrogen atoms in the ground level resulting
      in a bound-free transition.
    * *C0* (C-zero). Charge transfer collisional rate with Hydrogen atoms in the ground level.

      .. math::

          I^+ + H^0 \rightarrow I^0 + H^+


    * *C+*. Charge transfer collisional rate with protons.

      .. math::

          I^0 + H^+ \rightarrow I^+ + H^0


    * *CO* (C-Capital O). Charge transfer collisional rate with Hydrogen atoms in the ground level.

      .. math::

          I^0 + H^0 \rightarrow I^+ + H^-


    * *C-*. Charge transfer collisional rate with Hydrogen minus.

      .. math::

          I^+ + H^- \rightarrow I^0 + H^0


    .. note::

        Collisional rates of bound-bound transitions must be dexcitation rates (C\ :sub:`ul`).
        Collisional rates of bound-free transition must be ionization rates (C\ :sub:`lu`).
        Collisional rates for charge transfer events are non-symmetric, the charge
        transfer with neutral Hydrogen is a recombination (C\ :sub:`ul`), and the charge
        transfer with protons is an ionization (C\ :sub:`lu`).









