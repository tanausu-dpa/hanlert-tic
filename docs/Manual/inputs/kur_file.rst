The Kurucz file contains the list of atomic lines to include in the
background continuum (in LTE). It is in ASCII format and follows the
formatting of the linelist in the
`Kurucz database <http://kurucz.harvard.edu/atoms.html>`_.
From the webpage:

::

  1                                                                             80
  +++++++++++^^^^^^^++++++^^^^^^^^^^^^+++++^++++++++++^^^^^^^^^^^^+++++^++++++++++
     800.7110  0.116 27.00   45924.947  3.5 (3F)5s e2F   33439.661  4.5 (3F)4p y2G
      wl(nm)  log gf elem      E(cm-1)   J   label        E'(cm-1)   J'   label'
                     code
                          [char*28 level descriptor  ][char*28 level descriptor  ]
  continuing
  81                                                                           160
  ^^^^^^++++++^^^^^^++++^^++^^^++++++^^^++++++^^^^^+++++^+^+^+^+++^^^^^+++++^^^^^^
    8.19 -5.38 -7.59K88  0 0 59-2.584 59 0.000  104  -77F6 -5 0    1140 1165     0
    log   log   log ref NLTE iso log iso  log     hyper  F F'    eveglande     iso
   Gamma Gamma Gamma   level hyper f iso frac   shift(mK)     ^    oddglande shift
    rad  stark  vdW    numbers                    E    E'     ^abc  (x1000)   (mA)
                                                           I*1^char*3
                                                             codes
  FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
  3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3.2I5,I6)


* 1 wavelength (nm)  air above 200 nm   F11.4
* 2 log gf  F7.3
* 3 element code = element number + charge/100.  F6.2
* 4 first energy level  in cm-1   F12.3 (if allowed, with same parity as ground state) (negative energies are predicted or extrapolated) (sources for observed Es and gLandes are given in B* or C* files)
* 5 J for first level   F5.1
*  blank for legibility   1X
* 6 label field for first level   A10
* 7 second energy level  in cm-1   F12.3 (if allowed, with parity opposite first level)
* 8 J for second level   F5.1
*  blank for legibility   1X
* 9 label field for second level   A10
* 10 log of radiative damping constant, Gamma Rad  F6.2 or F6.3
* 11 log of stark damping constant/electron number. Gamma Stark  F6.2 or F6.3
* 12 log of van der Waals damping constant/neutral hydrogen number, Gamma van der Waals   F6.2 or F6.3
* 13 reference that can be expanded in gfxxyy*.lab   A4; Kxx are computed by Kurucz in year xx.  fourth character:  blank = E1, M = M1, N = M2, Q = E2, O = E4
* 14 non-LTE level index for first level   I2
* 15 non-LTE level index for second level   I2
* 16 isotope number   I3
* 17 hyperfine component log fractional strength  F6.3
* 18 isotope number  (for diatomics there are two and no hyperfine)   I3
* 19 log isotopic abundance fraction   F6.3
* 20 hyperfine shift for first level in mK to be added to E  I5
* 21 hyperfine shift for second level in mK to be added to E'  I5
*   the symbol "F" for legibilty   1X
* 22 hyperfine F for the first level    I1
* 23 note on character of hyperfine data for first level: z none, ? guessed  A1
*   the symbol "-" for legibility    1X
* 24 hyperfine F' for the second level  I1
* 25 note on character of hyperfine data for second level: z none, ? guessed  A1
* 26 1-digit code, sometimes for line strength classes   I1
* 27 3-character code such as AUT for autoionizing    A3
* 28 lande g for the even level times 1000   I5
* 29 lande g for the odd level times 1000   I5
* 30 isotope shift of wavelength in mA

.. note::
   Not all quantities are used in the code, but all must be present because the exact formatting is expected
