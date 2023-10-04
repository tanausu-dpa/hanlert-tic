In this section we describe what are the units of all the quantities that can be found
in the output files in alphabetical order.

================================ =========================================  ============================================  ===================================================================================================================================
Symbol                           Description                                         Files                                         Units (*hanlertio_class* or in file otherwise)               
================================ =========================================  ============================================  ===================================================================================================================================
:math:`b`                        departure coefficient                      :ref:`Departure file <departure_file>`        None                                                                  
:math:`\vec{C_i}`                contribution function (:math:`i=I,Q,U,V`)  :ref:`Contribution file <contribution_file>`  :math:`{\rm J}\cdot{\rm m^{-2}}\cdot{\rm s}^{-1}\cdot{\rm sr}^{-1}` 
:math:`C_{J\rightarrow J'}`      inelastic collisional rate between         :ref:`Cols-LL file <cols-tt_file>`            :math:`{\rm s}^{-1}`
                                 levels
:math:`C_{L\rightarrow L'}`      inelastic collisional rate between         :ref:`Cols-TT file <cols-tt_file>`            :math:`{\rm s}^{-1}`
                                 terms
:math:`I`                        intensity                                  :ref:`Solution file <solution_file>`          :math:`{\rm erg}\cdot{\rm cm^{-2}}\cdot{\rm s}^{-1}\cdot{\rm sr}^{-1}\cdot{\rm Hz}^{-1}\cdot c[{\rm cm}\cdot{\rm s}^{-1}]\cdot10^6`
:math:`I`                        intensity                                  :ref:`Stokesout file <stokesout_file>`        :math:`{\rm J}\cdot{\rm m^{-2}}\cdot{\rm s}^{-1}\cdot{\rm sr}^{-1}\cdot{\rm Hz}^{-1}`

                                                                            :ref:`StokesI file <stokesi_file>`

                                                                            :ref:`Stokes file <stokes_file>`
:math:`J`                        total angular momentum                     :ref:`Rhoout file <rhoout_file>`              None
:math:`\hat{J}^0_0(u\ell)`       integrated mean radiation field            :ref:`Solution file <solution_file>`          :math:`{\rm erg}\cdot{\rm cm^{-2}}\cdot{\rm s}^{-1}\cdot c[{\rm cm}\cdot{\rm s}^{-1}]\cdot10^6`
:math:`\hat{J}^0_0(u\ell)`       integrated mean radiation field            :ref:`Jout file <jout_file>`                  :math:`{\rm J}\cdot{\rm m^{-2}}\cdot{\rm s}^{-1}`
:math:`J^0_0(\omega)`            mean radiation field                       :ref:`Solution file <solution_file>`          :math:`{\rm erg}\cdot{\rm cm^{-2}}\cdot{\rm s}^{-1}\cdot{\rm Hz}^{1}\cdot c[{\rm cm}\cdot{\rm s}^{-1}]\cdot10^6`
:math:`\hat{J}^K_Q(u\ell)`       integrated radiation field tensors         :ref:`Solution file <solution_file>`          :math:`{\rm erg}\cdot{\rm cm^{-2}}\cdot{\rm s}^{-1}\cdot c[{\rm cm}\cdot{\rm s}^{-1}]\cdot10^6`
:math:`\hat{J}^K_Q(u\ell)`       integrated radiation field tensors         :ref:`Jout file <jout_file>`                  :math:`{\rm J}\cdot{\rm m^{-2}}\cdot{\rm s}^{-1}`
:math:`J^K_Q(\omega)`            radiation field tensors                    :ref:`Solution file <solution_file>`          :math:`{\rm erg}\cdot{\rm cm^{-2}}\cdot{\rm s}^{-1}\cdot{\rm Hz}^{-1}\cdot c[{\rm cm}\cdot{\rm s}^{-1}]\cdot10^6`
:math:`\hat{J}^{\rm P_1}(u\ell)` integrated photoionization rate            :ref:`Solution file <solution_file>`          :math:`10^8 {\rm s}^{-1}`
:math:`\hat{J}^{\rm P_2}(u\ell)` integrated recombination rate              :ref:`Solution file <solution_file>`          :math:`10^8 {\rm s}^{-1}`
:math:`N`                        number density                             :ref:`Solution file <solution_file>`          :math:`{\rm cm}^{-3}`

                                                                            :ref:`Population file <population_file>`

                                                                            :ref:`Rhoout file <rhoout_file>`
:math:`\vec{S}`                  Stokes parameters (:math:`I, Q, U, V`)     :ref:`Solution file <solution_file>`          :math:`{\rm erg}\cdot{\rm cm^{-2}}\cdot{\rm s}^{-1}\cdot{\rm sr}^{-1}\cdot{\rm Hz}^{-1}\cdot c[{\rm cm}\cdot{\rm s}^{-1}]\cdot10^6`

:math:`\vec{S}`                  Stokes parameters (:math:`I, Q, U, V`)     :ref:`Stokesout file <stokesout_file>`        :math:`{\rm J}\cdot{\rm m^{-2}}\cdot{\rm s}^{-1}\cdot{\rm sr}^{-1}\cdot{\rm Hz}^{-1}`

                                                                            :ref:`Stokes file <stokes_file>`
:math:`z`                        geometrical height axis                    :ref:`Jout file <jout_file>`                  :math:`{\rm cm}`

                                                                            :ref:`Rhoout file <rhoout_file>`

                                                                            :ref:`Contribution file <contribution_file>`
:math:`z_\tau`                   geometrical height where                   :ref:`Tau file <tau_file>`                    :math:`{\rm cm}`
                                 :math:`\tau(\omega) = 1`
:math:`\Gamma_{u\ell}/4\pi`      damping parameter                          :ref:`Damping file <damping_file>`            None
:math:`\epsilon(\omega)`         continuum emissivity                       :ref:`Background file <background_file>`      :math:`{\rm erg}\cdot{\rm cm}^{-2}\cdot{\rm s}^{-1}\cdot{\rm sr}^{-1}\cdot{\rm Hz}^{-1}`
:math:`\theta`                   polar angle                                :ref:`StokesI file <stokesi_file>`            :math:`{\rm degrees}`

                                                                            :ref:`Stokes file <stokes_file>`

                                                                            :ref:`Contribution file <contribution_file>`

                                                                            :ref:`Tau file <tau_file>`
:math:`\kappa(\omega)`           continuum absorptivity                     :ref:`Background file <background_file>`      :math:`{\rm cm}^{-1}`
:math:`\rho^K_Q(J,J')`           density matrix elements                    :ref:`Solution file <solution_file>`          None

                                                                            :ref:`Rhoout file <rhoout_file>`
:math:`\sigma(\omega)`           continuum scattering coefficient           :ref:`Background file <background_file>`      :math:`{\rm cm}^{-1}`
:math:`\tau`                     optical depth axis                         :ref:`Jout file <jout_file>`                  None

                                                                            :ref:`Rhoout file <rhoout_file>`

                                                                            :ref:`Contribution file <contribution_file>`
:math:`\tau_\tau`                optical depth where                        :ref:`Tau file <tau_file>`                    None
                                 :math:`\tau(\omega) = 1`
:math:`\chi`                     azimuthal angle                            :ref:`StokesI file <stokesi_file>`            :math:`{\rm degrees}`

                                                                            :ref:`Stokes file <stokes_file>`

                                                                            :ref:`Contribution file <contribution_file>`

                                                                            :ref:`Tau file <tau_file>`
:math:`\omega`                   frequency axis                             :ref:`Stokesout file <Stokesout_file>`        :math:`10^5\cdot{\rm cm}^{-1}`

                                                                            :ref:`StokesI file <stokesi_file>`

                                                                            :ref:`Stokes file <stokes_file>`

                                                                            :ref:`Contribution file <contribution_file>`

                                                                            :ref:`Tau file <tau_file>`

                                                                            :ref:`Background file <background_file>`
================================ =========================================  ============================================  ===================================================================================================================================
