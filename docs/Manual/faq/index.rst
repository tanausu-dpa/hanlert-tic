.. _faq:

**************************
Frequently Asked Questions
**************************

I updated the code and now it is broken!
----------------------------------------

Seems like while adding a feature, optimization, or trying to solve some issue, I broke something else. Please, go back to the previous version (you can check the version control with *git tree* or *git log*) and report the bug (by email or using the `issues tracker <https://gitlab.com/TdPA/hanlert-tic/-/issues>`_) so I can try solving it as soon as possible.

What is the reference direction for the polarization?
-----------------------------------------------------

In synthesis mode the reference direction of the polarization is set to the radial direction, i.e., the gamma angle in the geometrical tensors is set to zero. In order to transform to the more usual (when working with scattering polarization) parallel to the limb reference, both Q and U should be multiplied by -1.

The reference direction for the inversion file is the opposite, and the linear polarization is expected to be referred to the parallel to the limb direction.

How do you request reveral lines of sight in the same run?
----------------------------------------------------------

Most of the keywords that admit more than one element per entry expect spaces as separators. Therefore, for several heliocentric angles one would write, e.g., *POLAR_LOS = 0.1 0.5 1.0*.

How can you check the progress of an inversion/synthesis with more than one pixel?
----------------------------------------------------------------------------------

It is possible to check the progress in two ways. First, you can try opening any of the output files with the class in *Tools/hanlertio.py*; it allows reading partial files, so by requesting a *plane* slice you will visually see the progress. An easier way is using the tool *Tools/check_cache.py*. By running in your terminal::

    python <path_to_Tools_folder>/check_cache.py  <path_to_the_output_folder>/cache

the current progress will be printed.

Why does the code explode when I execute a run with polarization starting from a previous solution?
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

The Maximum Relative Change (MRC) suddendly decreased by orders of magnitude and went below the convergence criteria. Is that normal behavior?
----------------------------------------------------------------------------------------------------------------------------------------------

It is normal in the sense that there is nothing wrong necessarily, but the fulfilment of the convergence criteria is likely not real. In rare ocasions, an iteration can be such that there are changes in the radiation field, but the integrated radiation field tensors do not change too much, resulting in a very similar result in the Statistical Equilibrium Equations, which can satisfy the convergence criteria. If you are suspicious of this jump (which generally is the correct approach) just restart the run using the last solution file and, if the convergence was not real, the iterations will continue.

I received a warning saying that two wavelengths were considered to be in different ranges. What does it mean?
--------------------------------------------------------------------------------------------------------------

Because the code can consider bound-free and bounb-bound processes in very different parts of the spectrum, the wavelength (internally frequency) axis is not continuum in general, as it will be the concatenation of the necessary ranges. The code needs to perform integrals over frequencies, and not distinguishing between ranges would result in wrong weights in the integral. These messages are just notifying about the split in independent wavelength ranges.

I received a warning saying that two wavelengths were considered to be in different ranges, but that they were protected from splitting. What does it mean?
-----------------------------------------------------------------------------------------------------------------------------------------------------------

The criteria to consider if wavelengths are too far apart is controlled by the parameter *jump* in the source file *parameters_mod.f90*. However, specially with bound-free processes, it is possible that those wavelengths are part of the same range for a given process (e.g., the steps in the cross-section for a photoionization can be larger than the criteria to consider two wavelengths in different ranges). Because splitting those frequencies would make the integral over wavelengths wrong, they are checked and not split. The message is then innocuous.

I received an **IMPORTANT WARNING** saying that a transition can be shifted more than half of the total or the core width specifief. What is happening?
-------------------------------------------------------------------------------------------------------------------------------------------------------

In the :ref:`Atomic file <atom_file>`, each transition is assigned a range (in Doppler widths given in the :ref:`Control file <input_file>`) for the line itself, and for the line core. This warning is notifying that the maximum Doppler shift possible in the provided model atmosphere (or the maximum possible velocity in your inversion) can result in the line shifting more than half of the total width (so at least one fourth of the line can be out of range) or more than half of the core width (so the line center gets too close to the range defined as line core). It is **recommended** to change the line parameters in the atomic model accordingly, i.e., expanding the ranges and increasing the number of wavelength points if necessary.

I received a message saying that an inelastic collisional rate or cross-section was negative with Spline interpolation. Why?
----------------------------------------------------------------------------------------------------------------------------

That means that interpolating the tabulated values of the collisional rates or of the cross-section to the model's temperature resulted in a negative number somewhere. The code assumed that the reason was a problem with the splines and applied a different interpolation instead. If the code proceeded as normal, the new interpolation resulted in a positive number. While this warning is not a problem by itself, it is recommended to check what happens to the provided tabulations when applying splines in order to understand the magnitude of the problem and the goodness of the interpolations.

I received a warning saying that a transition from an atom had bad normalization in a number of heights, directions, and components. What is that about?
--------------------------------------------------------------------------------------------------------------------------------------------------------

The absorption profiles for each transition need to be normalized for every height, direction (in dynamic cases), and magnetic component to ensure that the energy is exactly conserved in the iterative process. Usually, the integral over the absorption profile is close to one, and the normalization is just taking care of small numerical errors. This message is warning that, for the provided parameters to discretize the wavelengths of the mentioned spectral line, the profile prior to the normalization integrated to something that was clearly different from unity (this is controled by the parameters *BADNORM* in the *parameters_mod.f90* source file). If this line is important, it is recommended to improve the assignment of ranges and wavelength numbers for that line. The number is just the count of how many of the individual profiles (for each height, direction, and magnetic component) integrated badly to unity before the normalization, giving an idea of how important it is (it may not be necessary to improve the parameters if this happens in a few points and they are expected to be out of the region of formation, for example).

I received a warning saying that a processor *reached RAM limit*. What is that about?
-------------------------------------------------------------------------------------

The amount of data in the form of redistribution functions can get really large values in the most general cases, making it impossible to fill everything in RAM. This message indicates that the *RAM_LIMIT* in the :ref:`Control file <input_file>` (which should definitely be set-up to a reasonable value) was reached and that not all the redistribution functions could be kept in RAM.

The linear polarization has some strange wiggles to one side of the line center. Is that real?
----------------------------------------------------------------------------------------------

Those wiggles are due to numerical noise and are often associated to a poor discretization of the wavelength in the redistribution integral when there is partial frequency redistribution. They can be minimized by improving the spectral resolution in the integral, controled by the parameters starting with *RED_* in the :ref:`Control file <input_file>`. 

My results for the intensity differ from those by the RH code
-------------------------------------------------------------

There are a lot of parameters that go into a non-LTE calculation. One important one are collisions. If you are comparing, e.g., calculations with the Ca II atom, the default model in both codes use a different set of inelastic collisions. Moreover, RH does not account for the impact of inelastic collisions in the line width, so there can be differences even if the collisions were the same. However, if you find a very clear difference between both codes and you are sure that the parameters are close enough, do not hesitate in using the contact information.

Was there not OpenMP support before?
------------------------------------

Yes, the code had support for the OpenMP library for parallelization for few years. However, the continuous development of new functionalities and changes in structure made it really difficult to maintain. Given the small to none gain in performance by using OpenMP or OpenMP+OpenMPI in comparison to just using OpenMPI, I decided to fully remove that support.
