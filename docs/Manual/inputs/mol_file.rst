The molecular model files have the following structure line by line.

1. The name of the molecule.

  ::

    H2

2. The mass of the molecule in atomic mass units.

  ::

    2.01588

3. The charge of the molecule.

  ::

    0

4. The list of atoms in the molecule, separated with a comma. In
   case that there are more than one type of an atom, write the
   number together and before the symbol.

  ::

    2H

  ::

    C, O

5. The dissociation energy in cm\ :sup:`-1`
   (:math:`\Delta E`).

  ::

    36117.506

6. The type of fit for the partition function (:math:`P`)
   and equilibrium coefficients (:math:`Q`). This define
   how the coefficients (:math:`C`, see points 8. and 9.)
   must be combined. Available strings are:

  * *KURUCZ_70*.

  .. math::

    P = \sum_{j=0}^{n-1} p_j,

  .. math::

    p_0 & = C_0 \\
    p_j & = p_{j-1}\cdot T + C_j

  .. math::

    Q = e^{\frac{\Delta E[J]}{k\cdot T} + \sum_{j=0}^{n-1} q_j -
            \frac{3}{2}\cdot\ln{T}\cdot({\rm N_{\rm atoms}} - 1 - {\rm charge})},

  .. math::

     q_0 & = C_0 \\
     q_j & = q_{j-1}\cdot T + C_j

  * *KURUCZ_85*.

  .. math::

    P = e^{\sum_{j=0}^{n-1} p_j},

  .. math::
     p_0 & = C_0 \\
     p_j & = p_{j-1}\cdot T\cdot 10^{-4} + C_j

  .. math::

    Q = e^{\frac{\Delta E[J]}{k\cdot T} + \sum_{j=0}^{n-1} q_j -
            \frac{3}{2}\cdot\ln{T}\cdot({\rm N_{\rm atoms}} - 1 - {\rm charge})},

  .. math::

     q_0 & = C_0 \\
     q_j & = q_{j-1}\cdot T\cdot10^{-4} + C_j

  * *SAUVAL_TATUM_84*. Only for neutral molecules.

  .. math::

    P = 10^{\sum_{j=0}^{n-1} p_j},

  .. math::

     p_0 & = C_0 \\
     p_j & = p_{j-1}*\log_{10}{\frac{5039.74756}{T}} + C_j

  .. math::

    Q = 10^{\left(\frac{\Delta E[eV]\cdot5039.74756}{T} - \sum_{j=0}^{n-1} q_j\right)
             \frac{k\cdot T\cdot10^4}
                  {0.01^{\left({\rm N_{\rm atoms}} - 1 - {\rm charge}\right)}}},

  .. math::

     q_0 & = C_0 \\
     q_j & = q_{j-1}*\log_{10}{\frac{5039.74756}{T}} + C_j

  * *IRWIN_81*. Only for neutral molecules.

  .. math::

    P = e^{\sum_{j=0}^{n-1} p_j},

  .. math::

    p_0 & = C_0 \\
    p_j & = p_{j-1}*\ln{T} + C_j

  .. math::

    Q = 10^{\left(\frac{\Delta E[eV]\cdot5039.74756}{T} - \sum_{j=0}^{n-1} q_j\right)
             \frac{k\cdot T\cdot10^4}
                  {0.01^{\left({\rm N_{\rm atoms}} - 1 - {\rm charge}\right)}}},

  .. math::

    q_0 & = C_0 \\
    q_j & = q_{j-1}*\log_{10}{\frac{5039.74756}{T}} + C_j

  * *TSUJI_73*. Only for neutral molecules.

  .. math::

    Q = 10^{\left(-\sum_{j=0}^{n-1} q_j\right)\cdot\left(k\cdot T\right)^2},

  .. math::

    q_0 & = C_0 \\
    q_j & = q_{j-1}*\frac{5039.74756}{T} + C_j

  ::

     KURUCZ_85

7. Minimum and maximum temperatures for the molecule presence.

  ::

     1.0E+3  9.0E+3

8. Number of coefficients for the partition function, followed by the coefficients.

  ::

    7   0.582145  16.3760  -49.4684  112.049  -149.953  106.531  -30.9791

9. Number of coefficients for the chemical equilibrium constant, followed
   by the coefficients.

  ::

    7  -46.4584  16.3660  -49.3992  111.822  -149.567  106.206  -30.8720

