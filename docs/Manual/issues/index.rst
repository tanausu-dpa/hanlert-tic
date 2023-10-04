.. _issues:

**********************
Known Issues and to do
**********************

To be (fully) documented.

A bit of history and justification
==================================

I started the development of HanleRT back on 2016, during my first postdoc and
in collaboration with Roberto Casini, with the objective to create a numerical
code capable of synthesizing the Stokes parameters of spectral lines with partial
frequency redistribution effects, implementing the (by the time) recently
published redistribution function for lambda-type transitions
(`Casini et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014ApJ...791...94C/abstract>`_).
This objective was achieved during the same year, resulting in the Letter that
we ask to reference when using HanleRT.

However, over the years, more and more functionalities have been added (with their
associated bugs) and an uncontable amount of bugs have been fixed. In 2018, within
the framework of the POLMAG project, Hao Li joined the team and developed the
Tenerife Inversion Code, which used HanleRT as a forward engine. A plethora of
changes were required in the synthesis code to improve the interaction with the
inversion code.

Later, in 2023, the Tenerife Inversion Code was fully integrated in what we now
call HanleRT-TIC. Significant structural changes were made, and I took the
opportunity to fuse into the main code another two branches that had remained
independent (and thus really difficult to maintain), the 1.5D and CLE modes.

However, the reduced man-power in the development and the reduced number of users
with access to the code during this development, together with the significant rewriting
of important sections of code and addition of potentially conflicting functionalities,
have made it really impossible to test every possible option and combination of
modes and inputs. Therefore, **it is to be expected that even the public version has
some broken things inside**.

For this reason, please, if you find any obvious (or subtle) bug, not-self-explained
warnings, installation/compilation/running errors, or strange behavior in execution
or in the results, contact *Tanausú del Pino Alemán* (e.g., through gitlab or via
email, which can be found in papers published later than 2018).

