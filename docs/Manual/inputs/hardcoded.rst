In the current version of HanleRT-TIC there are some quantities that, while
they could be read from a file, are hardcoded in the code. This
quantities are:

1. The full Hydrogen atom if it is not specified in the input
   (see **ATOM_INPUT** and **ATOM_BACK** in the
   :ref:`control file <input_file>`). This is hardcoded in
   the subroutine *AtomH* in the *ratom_mod.f90* module file.
   There are also hardcoded energies for such 6 levels atom
   in the *set_densities* and *set_Hdensities* subroutines
   in the *chemicaux_mod.f90* module file, in the
   *InitHpopu* subroutine in the *initpopuaux_mod.f90*
   module file, and in the *ReviseHatmo* subroutine in the
   *initpopu_mod.f90* module file.

2. The atomic abundances if the atom is not specified as
   an input in the :ref:`control file <input_file>` (see
   (**ATOM_INPUT** and **ATOM_BACK** fields) and no file
   is specified (**ABUND**). These hardcoded values are
   used in the equation of state and chemical equilibrium
   and they are tabulated in the :ref:`documentation
   <abundances_file>` and can be
   found in the *recallabund* and *recallabund_ind* subroutines
   in the *chemicaux_mod.f90* module file.

3. The atomic mass if the atom is not specified as
   an input in the :ref:`control file <input_file>` (see
   **ATOM_INPUT** and **ATOM_BACK** fields). These hardcoded
   values are used in the equation of state and chemical
   equilibrium for such non-loaded atoms and they can be
   found in the *recallmass_ind* subroutine in the
   *chemicaux_mod.f90* module file.
   
4. Molecular data used in the equation of state. In the last
   version of HanleRT, only the H2 and H2+ molecules are
   considered. The hardcoded functions can be found in
   the *moldata_ind* subroutine in the *chemicaux_mod.f90*
   module file.
