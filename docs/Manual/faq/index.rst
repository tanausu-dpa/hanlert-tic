.. _faq:

**************************
Frequently Asked Questions
**************************

I updated the code and now it is broken!
----------------------------------------

Seems like while adding a feature, optimization or trying to solve some issue I broke something else. Please, go back to the previous version (you can check the version control with *git tree* or *git log*) and report the bug (by email or using the `issues tracker <https://gitlab.com/TdPA/hanlert-tic/-/issues>`_) so I can try solving it as soon as possible.

What is the reference direction for the polarization?
-----------------------------------------------------

In synthesis mode the reference direction of the polarization is set to the radial direction, i.e., the gamma angle in the geometrical tensors is set to zero. In order to transform to the more usual (when working with scattering polarization) parallel to the limb reference, both Q and U should be multiplied by -1.

The reference direction for the inversion file is the opposite, and the linear polarization is expected to be referred to the parallel to the limb direction.

How do you request reveral lines of sight in the same run?
----------------------------------------------------------

Most of the keywords that admit more than one element per entry expect spaces as separators. Therefore, for several heliocentric angles one would write, e.g., *POLAR_LOS = 0.1 0.5 1.0*.

Why does The code explode when I execute a run with polarization starting from a previous solution?
---------------------------------------------------------------------------------------------------

The most common reason is that there was some compatibility problem with the chosen solution file and the radiation field could not be initialized correctly. If this happens, the following message should be present in the verbosity: *Warning: Number of frequencies in solution file different than in system; ignoring Stokes and J^K_Q(nu)*.

Why is the number of frequencies of the solution file different if I am using the same atom?
--------------------------------------------------------------------------------------------

This happens when trying to run a magnetic case with a multi-term atom while using a solution file from a non-magnetic run.

Short answer: add the line *MIT_OFF = No* to the :ref:`control file <input_file>` in the non-magnetic run.

Long answer: when neglecting :math:`J`-state interference in the atomic Hamiltonian or when there is no magnetic field, the total angular momentum :math:`J` is a good quantum number. The selection rules for an electric dipole transition indicate that a transition is allowed if :math:`|\Delta J|\le1` and both :math:`J` are different from zero.

However, when there is a magnetic field, if accounting for :math:`J`-state interference (multi-term atom) the total angular momentum :math:`J` is no longer a good quantum number. This means that the selection rules do not apply and the transitions with :math:`|\Delta J|>1` or with both :math:`J` equal to zero are no longer forbidden. These are the magnetically induced transitions (MIT).

For the sake of efficiency, a run with a multi-term atom without magnetic field will neglect any MIT, as there will not be any contribution to the line. However, introducing a magnetic field activates the MIT and will be included in the calculation. Consequently, the number of frequencies between the non-magnetic and magnetic cases will be different and the solution files will not be compatible.

To solve this, if the non-magnetic run is going to be used as an initialization of a magnetic problem, the line *MIT_OFF = No* should be added to the :ref:`control file <input_file>`. In this way, the non-magnetic case will include the MIT and the solution file will account for the MIT frequencies.

My inversion run stopped prematurely and I was using the STOREINV_STEP flag, so there is some work done in some of the pixels. Is that data wasted?
---------------------------------------------------------------------------------------------------------------------------------------------------

You can actually try to create a new initialization file. In the Tools folder there is a python file named *invo_injector.py* just for this. You can check how it is used by running in a terminal *python invo_injector.py help*. To summarize, this tool goes through the original :ref:`Result file <result_file>` and injects any valid model in the new (partial) :ref:`Result file <result_file>`, where valid means anything with an optical depth scale making sense (which means that something was written). A third file is created that can be used as a new initialization file. You can use the same name for the new file than for any of the existing ones, in which case a backup of the repeated one will be created in the same path; in this way, you can create an injected file with the same name than the original so the :ref:`control file <input_file>` does not need to be changed to continue the run. Of course, do not remove the :ref:`cache file <cache_file>`, as these partially solved columns are still flagged as not done.
