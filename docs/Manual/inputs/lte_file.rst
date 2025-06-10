SIR
---

This is the format of the Atomic parameter file in the
`SIR code <https://link.springer.com/article/10.1023/A:1002191418914>`_. From SIR's manual:

*Valid atomic parameter files have eight columns as follows:*

::

   Line=Ion  Wavelength  E  Ext.Pot  log(gf)   Transition   α     𝜎
    2=FE 1    6301.5012 1.0  3.654  -0.75   5P 2.0- 5D 2.0 0.243 2.3520e-14
    3=FE 1    6302.4936 1.0  3.686  -1.236  5P 1.0- 5D 0.0 0.240 2.3976e-14
    4=FE 1    5576.0888 1.0  3.428  -0.910  7D 1.0- 7D 0.0 0.232 2.3912e-14

*The first column gives the index with which the line is identified in*
*the profile and the wavelength grid files. The index is separated by*
*a = sign from the atomic symbol of the element. To specify the atomic*
*element, capital or lower case letters may be used (...)*
*For iron, either FE or XX can be employed. The ionization stage*
*is specified by a number: 1 means neutral atom, and 2 singly ionized*
*atom. (...)*

*The second column specified the (laboratory) central wavelength of the*
*transition (in Å). The third column gives the enhancement*
*factor to the van der Waals coefficient* :math:`\Gamma_6` *. The fourth*
*and fifth columns give the excitation potential of the lower level*
*(in eV) and the logarithm of the multiplicity of the level times the*
*oscillator strength, respectively. The sixth column specifies the atomic*
*transition. (...) Finally, the last two columns specify the collisional*
*broadening parameters* :math:`\alpha` *and* :math:`\sigma` *resulting*
*from the quantum mechanical theory of Anstee, Barklem, and O'Mara.*
:math:`\sigma` *is expressed in cm*\ :sup:`2`. *If these parameters are*
*zero, the classical Unsöld (1955) formula is used (together with the*
*enhancement factor mentioned above) for the calculation of the damping*
*factor.*

When using the SIR format in the :ref:`control file <input_file>`, it must include the full line, i.e.,

::

    LTE_LINE = 2=FE 1    6301.5012 1.0  3.654    -0.75   5P 2.0- 5D 2.0 0.243 2.3520e-14

.. note::
   The text in ellipsis has not been transcribed from SIR's manual
   because it does not apply for its HanleRT-TIC implementation.

.. warning::
   When choosing this format, some parameters are set to a default value. See warning at the end of this input for further details.

Kurucz
------

See :ref:`Kurucz File <kur_file>`.

.. warning::
   When choosing this format, some parameters are set to a default value. See warning at the end of this input for further details.

Native
------

This format allows for the most freedom in specifying LTE lines for HanleRT-TIC.
It can have up to 25 fields, where 21 are mandatory. The parameters are as follows:

#. Atomic label (e.g., Fe, Cr, etc.) or index (26, 24, etc.) designating the atomic species.
#. Ionization stage in the roman numeral format, i.e., I, II, III, etc., or as an integer (where 1 means neutral).
#. Energy of the lower level in cm\ :sup:`-1`.
#. Lower level total angular momentum.
#. Lower level Landé factor.
#. Energy of the upper level in cm\ :sup:`-1`.
#. Upper level total angular momentum.
#. Upper level Landé factor.
#. Einstein coefficient for spontaneous emission in s\ :sup:`-1`.
#. Five fields that characterize the Van der Waals broadening of the transition. The first field specifies the type, followed by four arguments. Five types are allowed: *unsold*, *barklem*, *kurucz*, *cross*, and *param*. See :ref:`Atom file <atom_file>` for further details.
#. Quadratic Stark broadening parameter. It can work as an enhancement or as a straightforward parameter. See :ref:`Atom file <atom_file>` for further details.
#. Collisional oscillator strength. Oscillator strength :math:`\Omega` (as in, e.g., `Bely & van Regemorter (1970) <https://ui.adsabs.harvard.edu/link_gateway/1970ARA&A...8..329B/doi:10.1146/annurev.aa.08.090170.001553>`_). This quantity is only used as a contribution to the line broadening.
#. Radiative broadening. Radiative :math:`\Gamma` parameter given in decimal logarithm. This quantity is approximately given by the natural width of the level (i.e., the sum of the Einstein coefficients for spontaneous emission for all possible transitions toward levels of lower energy).
#. Total number of frequencies that the line will add to the frequency axis.
#. Amount of frequencies, within the total, that will be located in what is considered the line core.
#. Number of Doppler widths (see *DOP_WIDTH* in :ref:`Control file <input_file>`) from the line center that the line will take up in the frequency axis.
#. Number of Doppler widths (see *DOP_WIDTH* in :ref:`Control file <input_file>`) from the line center that the core of the line will take up in the frequency axis.
#. Type of line. It can be *e1*, *m1*, *e2*, *m2*, or *un* (undefined). This field currently has no effect.
#. Lower limit of continuum optical depth (at the wavelength specified in the :ref:`Atmospheric file <atmo_file>` for *1D* or in the :ref:`Control file <input_file>` otherwise) where to consider the line to be present.
#. Upper limit of the temperature, above the minimum in the column, where to consider the line to be present.
#. Flag to make the line to not add frequencies to the frequency axis. It is activated by writing either yes or true.

.. warning::
   When any of the formats does not include one of the parameters listed in the native format, they are set to a default quantity. This default is hard-coded in the *rinput.py* file.

 | These defaults are the following:
 |           Quadratic Stark broadening. :math:`\Gamma=0`.
 |           Collisional oscillator strength. :math:`\Omega=0`.
 |           Radiative broadening. :math:`\Gamma={\rm log}_{10} A_{u\ell}` (calculated from the relevant input and only considering the line by itself).
 |           Amount of frequencies. :math:`N=41`.
 |           Amount of core frequencies. :math:`N_{\rm c}=31`.
 |           Doppler widths. :math:`\Delta\nu_D=10`
 |           Core Doppler widths. :math:`\Delta\nu_D=5`
 |           Type of line. e1.
 |           Lower limit of continuum optical depth. :math:`{\rm log}_{10}\tau_{\rm c}=-99.` (unlimited).
 |           Upper limit of the temperature. :math:`T = 10^{90}` (unlimited).
 |           Flag to avoid adding wavelengths. False.

