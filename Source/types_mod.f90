      !> Class definitions
      module types_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Hao Li (IAC/NSSCC)
!     Roberto Casini (HAO)
!  Start:
!     17/04/2017
!  Last version:
!     07/04/2025 V4.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     07/04/2025:    V4.0.4 - Changed s_inv_atmo and s_inv_res in
!                             Input_class to double precision (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!#####################################################################
!#####################################################################
!
!  To do:
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  Derived data types or structures:
!
!  box_class
!    Box with a double and a pointer to another box
!  tmp_col_box_class
!    Box to temporarily chain non-symmetric collisional rates
!  elastic_entry_class
!    Type to store individual entries for depolarizing elastic
!    collisional rates
!  elastic_class
!    Type to temporarily store depolarizing elastic collisional rates
!  inelastic_class
!    Type to temperarily store inelastic collisional rates
!  Tbox_class
!    Box to keep temperature tabulations
!  iabox_class
!    Box to link arrays when building frequency axes
!  dbabox_class
!    Box to link arrays when building frequency axes
!  strarr_class
!    Type with a long string (size 500)
!  strarr2_class
!    Type with a short string (size 2)
!  Freqflag_class
!    Transition absence flags
!  tranoindexb_class
!    Indexes for components of the redistribution function for a
!    pair of output and input transitions
!  tranoindex_class
!    Indexes for magnetic components for a transition and list of
!    input transitions in PRD
!  tranoIindex_class
!    List of input transitions for each transition in PRD
!  Phot_class
!    Photoionization data
!  FST_class
!    Fine structure transition data
!  rdip_class
!    Electric dipole strength for the components of a given transition
!  rdipev_class
!    Electric dipole strength for the components of a given transition
!    in the energy eigenbasis
!  Jrho_class
!    Indexing of KQ components
!  irho_class
!    Indexing of levels, J-J', and jM
!  Atom_class
!    Atomic model
!  LTEline_class
!    LTE line model
!  ckurucz_class
!    Atomic data for a Kurucz line
!  kurucz_class
!    Kurucz lines data
!  catm_class
!    Atomic data for a molecule
!  Mol_class
!    Molecular model
!  spect_class
!    Input spectra for CLE
!  chianti_real_pointer_class
!    Pointer to a real one-dimensional array
!  chianti_ioneq_class
!    Ionization equilibrium and ion fraction data
!  chianti_class
!    CHIANTI data structure
!  pf_class
!    Partition function and ionization energy data for an atom
!  Atmo_class
!    Model atmosphere
!  Bfield_class
!    Magnetic field stratification
!  Coronapoint_class
!    Geometry data for a point in CLE
!  Geometry_class
!    Geometry quadrature and LOS data
!  IO_helper_class
!    Data to restrict output ranges
!  FWHM_helper_class
!    Data on the spectral point spread function for inversions
!  strnum_class
!    Ionization fraction input for CLE
!  Node_class
!    Nodes location, value, error, and limits for a variable
!  Input_class
!    Configuration data
!  fudge_class
!    Opacity fudge data
!  Continuum_class
!    Background opacities
!  MRC_class
!    Maximum relative change data
!  MPI_class
!    MPI data
!  scalar, a1D, a2D, a3D, a4D, a5D, a6D, a7D, a8D, a9D
!    memoization
!  Fctsg_class
!    Factorials, signs, and J-symbols
!  Frequency_class
!    Frequency data
!  Redc2_class
!    Redistribution functions
!  Redb2_class
!    Input transitions for redistribution functions
!  Redc_class
!    Input frequency axes and weights
!  Redb_class
!    Output frequency ranges for redistribution, input transitions
!    for frequency axes and weights, and second order RT coefficients
!  Prof_class
!    Voigt and Faraday-Voigt profiles, and Voigt normalization factors
!  Red_class
!    Redistribution input frequency data, redistribution function
!    data, and profile or normalization data
!  Rhoc_class
!    Temporal rhoKQ or population data
!
!  !!!!!!!!!!!
!  TIC classes
!  !!!!!!!!!!!
!
!  Nodes_class
!    Inversion node data
!  Solution_class
!    Frequency and synthetic Stokes parameters in the frequency range
!    of the inverted data
!  Solution_F_class
!    Solution of the self-consistent problem and the corresponding
!    emergent profiles, contribution function, and height for optical
!    depth equal to one
!  Stokes_class
!    Inversion Stokes parameters data
!  Regul_class
!    Regularization data
!  LMFIT_class
!     Data for the Levenberg–Marquardt
!
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod

      integer, parameter:: nvar_inv = 14

!#####################################################################

      !> Box with a double and a pointer to another box
      type box_class

        ! Number
        double precision:: val

        ! Pointer to next box
        type(box_class), pointer:: next

      end type box_class

!#####################################################################

      !> Box to temporarily chain non-symmetric collisional rates
      type tmp_col_box_class

        ! Initial and final level, forbidden flag
        integer:: ifrom,ito,flag

        ! Rates
        double precision, dimension(:), allocatable:: C

        ! Pointer to next box
        type(tmp_col_box_class), pointer:: next

      end type tmp_col_box_class

!#####################################################################

      !> Type to store individual entries for depolarizing elastic
      !! collisions
      type elastic_entry_class

        ! K, type, nz
        integer:: K,typo,nz

        ! Coefficients
        double precision:: a,b,c

        ! Explicit
        double precision, dimension(:), allocatable:: Coeff

      end type elastic_entry_class

!#####################################################################

      !> Type to temporarily store depolarizing elastic collisional
      !! rates
      type elastic_class

        ! Level, number of entries
        integer:: ilevel, nentry

        ! Substructure with data
        type(elastic_entry_class), dimension(:), allocatable:: datum

      end type elastic_class

!#####################################################################

      !> Type to temperarily store inelastic collisional rates
      type inelastic_class

        ! Type of collision, up, low, forbidden, index in Tbox
        integer:: col_type,up,low,forbid,ind

        ! Substructure with data
        double precision, dimension(:), allocatable:: Cul

      end type inelastic_class

!#####################################################################

      !> Box to keep temperature tabulations
      type Tbox_class

        ! Type of interpolation
        logical:: flin

        ! Number of indexes, index, type of collision
        integer:: nTmp,ind,col_type,nion

        ! Temperatures
        double precision, dimension(:), allocatable:: temp

        ! Pointer to next
        type(Tbox_class), pointer:: next

      end type Tbox_class

!#####################################################################

      !> Box to link arrays when building frequency axes
      type dbabox_class

        ! Integers
        integer:: ifreq,mfreq,nback,iph,ith

        ! Array
        double precision, dimension(:), allocatable:: A

        ! Pointer to next
        type(dbabox_class), pointer:: next,prev

      end type dbabox_class

!#####################################################################

      !> Type with a long string (size 500)
      type strarr_class

        ! String with 500 characters
        character(len=500):: str

      end type strarr_class

!#####################################################################

      !> Type with a short string (size 2)
      type strarr2_class

        ! String with 2 characters
        character(len=2):: s

      end type strarr2_class

!#####################################################################

      !> Transition absence flags
      type Freqflag_class

        ! OR() of the vector flag Vabsent
        logical:: absent

        ! OR() of the vector flag Vabsent for Master
        logical, dimension(:), allocatable:: Mabsent

        ! Vector with flag of presence of a line in that frequency
        logical, dimension(:), allocatable:: Vabsent

      end type Freqflag_class

!#####################################################################

      !> Indexes for components of the redistribution function for a
      !! pair of output and input transitions
      type tranoindexb_class

        ! Number of coherences between lower terms
        integer:: nchlt

        ! Indexing of the mu and M of a term, same without magnetic
        ! field (J ordering)
        integer, dimension(:,:,:,:,:), allocatable:: indB,indNB

      end type tranoindexb_class

!#####################################################################

      !> Indexes for magnetic components for a transition and list of
      !! input transitions in PRD
      type tranoindex_class

        ! Number of input transitions,
        integer:: nt,ncomB,ncomNB

        ! Indexes of transitions
        integer, dimension(:), allocatable:: indT

        ! Indexing for the mu and M of a term for the first order
        ! part, same without magnetic field (jM order)
        integer, dimension(:,:), allocatable:: indB,indNB

        ! Indexing of the mu and m of a term
        type(tranoindexb_class), dimension(:), allocatable:: trani

      end type tranoindex_class

!#####################################################################

      !> List of input transitions for each transition in PRD
      type tranoIindex_class

        ! Number of input transitions,
        integer:: nt

        ! Indexes of transitions
        integer, dimension(:), allocatable:: indT

      end type tranoIindex_class

!#####################################################################

      !> Photoionization data
      type Phot_class

        ! Non presence of the photoionization in the CPU
        logical:: absent

        ! Non presence of the photoionization in the CPU for Master
        logical, dimension(:), allocatable:: Mabsent

        ! Upper and lower levels
        integer:: ilevell,ilevelu

        ! Type of input, number of frequencies for the cross section,
        ! initial frequency index, final frequency index
        integer:: mode,nfreq,if0,if1

        ! initial frequency index from master, final frequency from
        ! master
        integer, dimension(:), allocatable:: Mif0,Mif1

        ! Photoionization edge, (2J_l+1)/(2J_u+1), frequency weights
        ! in the boundaries of the subdomain
        double precision:: edge,glu,W0,W1

        ! Proper variation of cross section with frequency, integral
        ! for T_E SEE rate
        double precision, dimension(:), allocatable:: alpha,TEI
        ! Frequencies in the input atomic model and corresponding
        ! cross sections
        double precision, dimension(:), allocatable:: infreq,inalpha
        ! Frequency weights in the boundaries of the subdomain for
        ! master
        double precision, dimension(:), allocatable:: MW0,MW1

      end type Phot_class

!#####################################################################

      !> Fine structure transition data
      type FST_class

        ! Upper and lower terms
        integer:: iterml,itermu

        ! Number of FS transitions
        integer:: nt

        ! Index of upper and lower J levels
        integer, dimension(:), allocatable:: ilevell,ilevelu

        ! Index of FS transition (equivalent to irad in Atom_class)
        integer, dimension(:,:), allocatable:: irad

        ! Aul and Blu for each FS transition within a transition
        ! between terms
        double precision, dimension(:,:), allocatable:: Aul,Blu

      end type FST_class

!#####################################################################

      !> Electric dipole strength for the components of a given
      !! transition
      type rdip_class

        ! Dipole strength matrix
        double precision, dimension(:,:,:,:,:), allocatable:: rdip

      end type rdip_class

!#####################################################################

     !> Electric dipole strength for the components of a given
     !! transition in the energy eigenbasis
      type rdipev_class

        ! Dipole strength matrix
        type(rdip_class), dimension(:), allocatable:: rdipev

      end type rdipev_class

!#####################################################################

      !> Structure with the level-level indexing
      type Jrho_class

        ! K,Q indexing
        integer, dimension(:,:), allocatable:: kq

      end type Jrho_class

!#####################################################################

      !> Indexing of levels, J-J', and jM
      type irho_class

        ! Level-level indexing
        type(Jrho_class), dimension(:,:), allocatable:: Jrho

        ! Indexing of rho by term
        integer, dimension(:), allocatable:: irho_ij

        ! Indexing of j,M by term
        integer, dimension(:,:), allocatable:: jM

      end type irho_class

!#####################################################################

      !> Atomic model
      type Atom_class

        ! Information of the presence of a transition at some
        ! frequency
        type(Freqflag_class), dimension(:), allocatable:: fflag

        ! FS transitions information
        type(FST_class), dimension(:), allocatable:: fst

        ! Photoionization information
        type(Phot_class), dimension(:), allocatable:: phot

        ! Indexing for transitions (out and in)
        type(tranoindex_class), dimension(:), allocatable:: trano

        ! Indexing for transitions (out and in) for intensity
        type(tranoIindex_class), dimension(:), allocatable:: tranoI

        ! Pointer to non-symmetric rates
        type(tmp_col_box_class), pointer:: Ccoeff_special

        ! Array of inelastic collisional data
        type(inelastic_class), dimension(:), allocatable:: inelas

        ! Array of elastic collisional data
        type(elastic_class), dimension(:), allocatable:: elas

        ! Boxes of temperature tables
        type(Tbox_class), pointer:: Tbox

        ! Dipole strength array
        type(rdip_class), dimension(:), allocatable:: rdip

        ! Dipole strength array in energy representation
        type(rdipev_class), dimension(:), allocatable:: rdipev

        ! Atom indexing
        type(irho_class), dimension(:), allocatable:: irho

        ! Name of the atomic species
        character(len=2):: Element

        ! Label
        character(len=10):: file_label

        ! Normalize the total relative factor for abundance or not,
        ! multilevel flag, keep populations fixed, make zero the
        ! last ion, fix populations for lower term
        logical:: anorm,ML,fixp,zero_ion,fixplt

        ! Input spectrum for b-b and b-f transitions
        logical, dimension(:), allocatable:: bbspecin,bfspecin

        ! Flag to identify automatic Hydrogen background atom
        logical:: cust=.False.

        ! Flag to protect agains molecular equilibrium
        logical:: mol_protect=.False.

        ! Flag for the second order emissivity for each transition,
        ! flag to split between components in the multiterm atom when
        ! building frequency axis
        logical, dimension(:), allocatable:: lemiss2,splitf

        ! Flag to nullify density matrix elements if they are small
        logical, dimension(:,:), allocatable:: rhonull

        ! Non-coherent lower term
        logical, dimension(:,:), allocatable:: NCHLT

        ! Number of terms, transitions, maximum value of J, maximum
        ! value of K, number of magnetic levels, number of levels,
        ! size of the SEE system, number of output frequencies,
        ! number of photoionizations, number of depolarizing
        ! collisions entries, number of inelastic collisions entries
        ! total number of levels and number of FS transitions, shift
        ! in transition index, shift in photoionization index, shift
        ! in fine structure transitions, number of transitions with
        ! redistribution
        integer:: nMulti,ntran,nJmax,nKmax,nMmax,NNN,ndim,nfreq, &
                  nphot,ngk,ncol,nlevel,nftran,tshift,pshift, &
                  tfshift,ntrano

        ! Sizes of profile files headers
        integer:: hvifil,hwifil,hwtfil

        ! Number of frequency points for the transition, number of
        ! frequency points for the line core, ionization stage, type
        ! of Van der Waals broadening, term index, sublevel index,
        ! number of frequencies of photoionization input, type of
        ! inelastic collision, lower and upper limits of transitions
        ! in frequency indexes, term transition given FS,
        ! indexing of transitions with redistribution, initial
        ! and final frequency for Master, initial and final PRD
        ! frequencies for Master, initial and final frequency indexes
        ! in the input spectra for each bb and bf transition,
        ! initial and final index for the line and photoionization
        ! transitions for the CLE integrals, term-wise K cut,
        ! line-wise K cut, true (uncut) limits for transitions not
        ! mutable, sublevel transition given FS
        integer, dimension(:), allocatable:: nfreqt,nfreqtc,nJ, &
                                             stage,broad_type,term, &
                                             sublevel,nfreqph, &
                                             col_type,if0,if1,ifst, &
                                             itrano,rif0,rif1,rif20, &
                                             rif21,sbif0,sbif1, &
                                             sfif0,sfif1,ilf0,ilf1, &
                                             ipf0,ipf1,Kcut,Krad, &
                                             tif0,tif1,ifstj

        ! Transition indexing matrix, collisional transition indexing
        ! matrix, photoionization indexing matrix initial frequency
        ! index from master, final frequency from master, initial
        ! index for comoving version, final index for comoving version
        integer, dimension(:,:), allocatable:: irad,icol,iphot,Mif0, &
                                               Mif1,CMif0,CMif1

        ! Indexing of the FS transitions
        integer, dimension(:,:), allocatable:: ifst_ij

        ! Block dimension (index for M,term) and flag of forbidden
        ! collisions between levels
        integer, dimension(:,:), allocatable:: nblk,fcflag

        ! J indexing for each individual magnetic level
        integer, dimension(:,:,:), allocatable:: iJval

        ! Atomic part of the thermal Doppler width, atomic mass,
        ! abundance, multiplicative factor for abundance
        double precision:: cDopp,rmass,abun,abun_mod = 1d0

        ! orbital angular momentum, spin angular momentum, term
        ! frequency, Doppler widths to include, Doppler widths for
        ! the core, Frequency of term transition, number density of
        ! the atom in cm-3, Stark broadening, degeneration (terms
        ! and H custom), lower and upper frequency limit weights
        double precision, dimension(:), allocatable:: rLval,Sval, &
                                                      TRfreq,Dwvl, &
                                                      Dwvlc,Dfreq,n, &
                                                      broad_stark, &
                                                      deg,W0,W1,gL

        ! Inverse lifetime, level frequency and J values, radiative
        ! transition rates matrix, arguments for the Van der Waals
        ! Broadening, collisional transition-wise damping,
        ! populations, lte populations, departure coefficient,
        ! lower and upper frequency limit weights for master, elastic
        ! rates
        double precision, dimension(:,:), allocatable:: damp,FSfreq, &
                                                        rJval, &
                                                        Ecoeff, &
                                                        broad_args, &
                                                        ldamp,popu, &
                                                        populte, &
                                                        depar,MW0, &
                                                        MW1,qel

        ! Inelastic collisional terms rates for allowed transitions
        ! between terms and transition rates between levels
        double precision, dimension(:,:,:), allocatable:: Ccoeff, &
                                                          CcoeffJ

        ! Eigenvalues of the diagonalization
        double precision, dimension(:,:,:,:), allocatable:: eval

        ! Elastic collisional rates
        double precision, dimension(:,:,:,:,:), allocatable:: gk

        ! Eigenvectors of the diagonalization
        double precision, dimension(:,:,:,:,:), allocatable:: evec

        ! Vector with rhoKQ for the current
        complex(kind=8), dimension(:,:), allocatable:: crho

      end type Atom_class

!#####################################################################

      !> LTE line model
      type LTEline_class

        ! Name of the atomic species
        character(len=2):: Element

        ! If there is a background model atom, if limited in tau, if
        ! limited in T, line is absent for a CPU
        logical:: is_passive,taulim_l,Tlim_l,absent

        ! Element, ion, number of Ml, number of Mu, maximum height
        ! index for presence, total number of frequencies, number
        ! of frequencies for the line core, initial frequency, final
        ! frequency, index of the passive model atom, type of Van
        ! der Waals broadening, type of line, number of magnetic
        ! components
        integer:: ele,stage,nMl,nMu,Rz0,nfreq,nfreqc,if0,if1,ia, &
                  broad_type,ltype,ncom

        ! Upper level energy, lower level energy, upper level angular
        ! momentum, lower level angular momentum, upper level Landé
        ! factor, lower level Landé factor, Einstein coefficient for
        ! e.e., resonance frequency, Einstein coefficient for
        ! absorption, optical depth limit, temperature limit, range in
        ! Doppler widths for the line, range in Doppler widths for the
        ! core, Stark broadening, collisional oscillator strength,
        ! radiative broadening, abundance
        double precision:: eu,el,Ju,Jl,gu,gl,Aul,Dfreq,Blu,taulim, &
                           Tlim,Dwvl,Dwvlc,broad_stark,f_c, &
                           broad_rad,abund

        ! Inverse lifetime, inverse lifetime upper level, arguments
        ! for the Van der Waals broadening, collisional
        ! transition-wise damping
        double precision, dimension(:), allocatable:: damp,broad_args

        ! Lower level population (fraction), upper level population
        ! (fraction), total population
        double precision, dimension(:), allocatable:: nl,nu,n

        ! Atomic part of the thermal Doppler width, atomic mass,
        double precision:: cDopp,rmass

      end type LTEline_class

!#####################################################################

      !> Atomic data for a Kurucz line
      type ckurucz_class

        ! Element, charge
        integer:: A,Z

        ! Einstein coefficient, resonance, radiative, Van der Waals,
        ! and Stark broadening parameters, statistical weights upper
        ! and lower level, lower and upper limit in frequency,
        ! abundance, lower level energy
        double precision:: Aul,Dfreq,Grad,Gvdw,Gstk,gu,gl,O0,O1, &
                           abun,Ei

      end type ckurucz_class

!#####################################################################

      !> Kurucz lines data
      type kurucz_class

        ! Number of lines
        integer:: ntran

        ! Line data
        type(ckurucz_class), dimension(:), allocatable:: tran

        ! Maximum frequency distance
        double precision:: MDomg

      end type kurucz_class

!#####################################################################

      !> Atomic data for a molecule
      type catm_class

        ! Atom ID
        character(len=2):: s

        ! In model
        logical:: inmod

        ! Molecules where it is present
        integer:: pnmol = 0

        ! Number of pf stages
        integer:: nstg

        ! Indexes of molecules where it is present, and how many time
        integer, dimension(:), allocatable:: imol,nmol

        ! Abundance
        double precision:: abun = 0

        ! Ionization eng, PF single height
        double precision, dimension(:), allocatable:: Eion,pfs

        ! PF, height dependent
        double precision, dimension(:,:), allocatable:: pf

      end type catm_class

!#####################################################################

      !> Molecular model
      type Mol_class

        ! Name of atoms in molecule
        type(strarr2_class), dimension(:), allocatable:: catom

        ! Name of the molecule
        character(len=:), allocatable :: Molecule

        ! Charge, number of species, total number of atoms, type
        ! of fit, number of pf coefficients, number of eqc
        ! coefficients
        integer:: Charge,nA,nAT,pffit,npfcoeff,neqcoeff

        ! Number of each atom in the molecule and position in the
        ! atom list for chemical equilibrium
        integer, dimension(:), allocatable:: natom,iatom

        ! Atomic part of the thermal Doppler width, mass, dissociation
        ! energy, minimum temperature, maximum temperature
        double precision:: cDopp,rmass,Den,Tmin,Tmax

        ! Coefficients for partition function, partition function,
        ! Coefficients for equililibrum constant, equilibrium
        ! constant, Population
        double precision, dimension(:), allocatable:: pfcoeff,pf, &
                                                      eqcoeff,eq,n

      end type Mol_class

!#####################################################################

      !> Input spectra for CLE
      type spect_class

        ! Axial
        logical:: axial,valid,pol

        ! Size of spectrum
        integer:: nfreq,nmu,nphi,nstk

        ! Wavenumber, cosine of polar angle, phi angle
        double precision, dimension(:), allocatable:: omega,mu,phi
        ! Stokes spectrum
        double precision, dimension(:,:,:,:), allocatable:: stokes

        ! Stokes spectrum interpolated to the directional quadrature
        ! in a point
        double precision, dimension(:,:,:,:), allocatable:: mustokes

      end type spect_class

!#####################################################################

      !> Pointer to a real one-dimensional array
      type chianti_real_pointer_class

        ! pointer
        real, dimension(:), pointer:: p

      end type chianti_real_pointer_class

!#####################################################################

      !> Ionization equilibrium and ion fraction data
      type chianti_ioneq_class

        ! Dimensions
        integer:: nI

        ! Array containing pointers to ioneq_data
        type(chianti_real_pointer_class), dimension(:), &
                                          allocatable:: stage

        ! Splines data
        double precision, dimension(:,:), allocatable:: b,c,d


      end type chianti_ioneq_class

!#####################################################################

      !> CHIANTI data structure
      type chianti_class

        ! Dimensions
        integer:: nT,nE,nioneq

        ! Full ioneq data
        real, dimension(:), pointer:: ioneq_data

        ! Pointer to temperature
        real, dimension(:), pointer:: ioneq_T

        ! Array containing pointers to ioneq_data
        type(chianti_ioneq_class), dimension(:), &
                                   allocatable:: ioneq

      end type chianti_class

!#####################################################################

      !> Partition function and ionization energy data for an atom
      type pf_class

        ! Atom ID
        character(len=2):: Element

        ! Number of ion stages
        integer:: nstg

        ! Ionization energy
        double precision, dimension(:), allocatable:: Ei

        ! Partition function
        double precision, dimension(:,:), allocatable:: pf

      end type pf_class

!#####################################################################

      !> Model atmosphere
      type Atmo_class

        ! Type of scale in input
        character(len=1) scal

        ! Atom-atmospheric data
        type(pf_class), dimension(:), allocatable:: ele

        ! If each group of pointers is allocated, and not pointing
        logical:: alloc_a, alloc_b

        ! Indexes for quick location of variables in CLE
        integer:: d0,di0,di1,it,imi,iv,ine,inh,ipe,ipg,irh,ib

        ! Number of height nodes, type of density input, number of
        ! temperatures in partition function table, number of elements
        ! in partition function, atmospheric CLE mode, atmospheric
        ! CLE normalization of spatial coordinates
        integer:: nZ,typo,NT,nele,mode,norm

        ! log of gravity acceleration at surface, frequency
        ! of tau scale, y position for CLE, z position for CLE,
        ! diffuse light fraction
        double precision:: logg,tfreq,ypos,zpos,f_diff

        ! Height, Temperature, microturbulence, velocity components,
        ! magnetic field components, auxiliar velocity components to
        ! cheat in intensity
        double precision, dimension(:), pointer:: z,T,vmi,vx,vy,vz, &
                                                  Bx,By,Bz,vxa,vya,vza

        ! Total hydrogen, hydrogen minus, gas pressure, density,
        ! electron pressure, alternative height scale, and atomic
        ! hydrogen number density, partition function temperature,
        ! abundance, input asymmetry radiation field tensors, to
        ! be used by the inversion, electron number density, los
        ! velocity, pos velocity, pos azimuth velocity
        double precision, dimension(:), allocatable:: nHT,nHm,Pg, &
                                                      rho,Pe,zalt, &
                                                      nHA,pT,abund, &
                                                      JKQin,ne,vlos, &
                                                      vpos,vphi

        ! Continuum absorption at reference frequency
        double precision, dimension(:), allocatable:: chi500

        ! Array with zeros
        double precision, dimension(:), pointer:: zeros

        ! Hydrogen density, helium density
        double precision, dimension(:,:), allocatable:: nh,nhe

      end type Atmo_class

!#####################################################################

      !> Magnetic field stratification
      type Bfield_class

        ! Module, polar angle and azimuth of B field
        double precision, dimension(:), allocatable:: Bstrength, &
                          Btheta,Bphi

        ! Longitudinal magnetic field, transversal magnetic field, and
        ! magnetic field azimuth in the POS (these are used in
        ! inversion)
        double precision, dimension(:), allocatable:: Blos,Bpos, &
                                                      Azimuth

      end type Bfield_class

!#####################################################################

      !> Geometry data for a point in CLE
      type Coronapoint_class

        ! Heliocentric angle from local vertical, azimuth angle from
        ! local vertical, cosine of a angle, sine of a angle, cosine
        ! of b angle, sine of b angle, cosine of gamma angle, sine of
        ! gamma angle
        double precision:: theta,phi,CA,SA,CB,SB,CY,SY

        ! True theta, phi, and gamma angles in vertical frame
        double precision, dimension(3):: geom

        ! CLV geometrical parameters given a height
        double precision, dimension(6):: CLV

      end type Coronapoint_class

!#####################################################################

      !> Geometry quadrature and LOS data
      type Geometry_class

        ! Flag for axial symmetry
        logical:: axial

        ! Number of polar nodes, number of real azimuthal nodes,
        ! number of azimuthal nodes for emiss2ord, number of emergent
        ! polar directions, number of emergent azimuthal directionsm
        ! number of nodes for AA integral, number of scattering
        ! angles, number of current output directions
        integer:: nTh,nPh,nPh2,nThLOS,nPhLOS,nThAA,nScatt,njdir

        ! Index of polar direction for running index, index of
        ! azimuth direction for running index
        integer, dimension(:), allocatable:: ithv,iphv

        ! Indexing of 2D directions
        integer, dimension(:,:), allocatable:: i_geom

        ! Indexing of scattering angles
        integer, dimension(:,:,:), allocatable:: i_scatt

        ! Gamma angle
        double precision:: gam

        ! Vector of cosines of polar angle nodes, vector of cosiones
        ! of azimuthal angle nodes, vector of sign of sinus of
        ! azimuthal angle nodes, vector of polar angles, vector of
        ! azimuthal angles, weights of polar integral, weights of
        ! RT azimuthal integral, weights of emiss2 azimuthal
        ! integral, cosines of polar angle of emergent directions,
        ! azimuthal angles of emergent directions, Vector of cosines
        ! for the AA integral, Vector of sines for the AA integral,
        ! weights for the AA integral, nodes for classical gaussian
        ! quadrature, weights for classical gaussian quadrature, polar
        ! angle on the disk for a given quadrature in a point above
        ! the surface, cosine of scattering angles, sine of scattering
        ! angles
        double precision, dimension(:), allocatable:: V_mu,V_mux, &
                                                      V_muy,V_theta, &
                                                      V_phi,W_mu, &
                                                      W_mux,W_mux2, &
                                                      L_mu,L_theta, &
                                                      L_phi,V_muAA, &
                                                      V_siAA,W_muAA, &
                                                      V_gauss, &
                                                      W_gauss, &
                                                      V_mu_disk, &
                                                      V_CScatt, &
                                                      V_SScatt

        ! TKQ geometrical tensor in the vectical reference frame
        ! quadrature and LOS
        complex(kind=8), dimension(:,:,:,:), pointer:: TS,TSo,TSL

        ! TKQ geometrical tensor in the magnetic reference frame
        ! quadrature and LOS
        complex(kind=8), dimension(:,:,:,:,:), pointer:: TB,TBo,TBL

      end type Geometry_class

!#####################################################################

      !> Data to restrict output ranges
      type IO_helper_class

        ! Buffer size
        integer:: buffer_size=1

        ! Header size, number of 'ranges', integer for header
        integer:: head_size,nran,nn,geom_size

        ! Secondary indexes, size per range
        integer, dimension(:), allocatable:: sindx,nbuff

        ! Indexes to control output
        integer, dimension(:,:), allocatable:: indx

        ! Doubles to control output
        double precision, dimension(:,:), allocatable:: doub

      end type IO_helper_class

!#####################################################################

      !> Data on the spectral point spread function for inversions
      type FWHM_helper_class

        ! If single gaussian value, if pending initialization
        logical:: gaussian,toinit

        ! Number of ranges (duplicated), number of wavelengths (if
        ! not gaussian)
        integer:: nn,nfreq

        ! Secondary indexes, size per range
        integer, dimension(:), allocatable:: sindx

        ! Indexes to control output, indexes for quick
        ! interpolation
        integer, dimension(:), allocatable:: indx

        ! Indexes for quick interpolation
        integer, dimension(:,:), allocatable:: indx1,indx2

        ! Doubles to control output, wavelength, kernel
        double precision, dimension(:), allocatable:: doub,wave,kernel

        ! Auxiliar for quick interpolation
        double precision, dimension(:,:), allocatable:: idx

      end type FWHM_helper_class

!#####################################################################

      !> Ionization fraction input for CLE
      type strnum_class

        ! String with 500 characters
        character(len=500):: str

        ! Type
        integer:: typ

        ! Value
        double precision:: val

      end type strnum_class

!#####################################################################

      !> Nodes location, value, error, and limits for a variable
      type Node_class

        ! Number of special boundary conditions
        integer:: nebound = 0

        ! Indexes in atmosphere for values of tau in nodes
        integer, dimension(:), allocatable:: Tau_Indx

        ! Normal boundary limits
        double precision, dimension(2):: Bounds

        ! Positions, values, errors
        double precision, dimension(:), allocatable:: H,Var,Errors

        ! Positions and limits for special boundary conditions
        double precision, dimension(:,:), allocatable:: ebound

      end type Node_class

!#####################################################################

      !> Configuration data
      type Input_class

        ! Structures to help with 1.5D outputs
        type(IO_helper_class):: lim_stk,lim_ctr,lim_tau,lim_cols_tt, &
                                lim_cols_ll,lim_damp,lim_back, &
                                lim_pop,lim_atmo,lim_qel

        ! Structure with FWHM info
        type(FWHM_helper_class), dimension(:), allocatable:: lim_fwhm

        ! LTElines data
        type(LTEline_class), dimension(:), allocatable:: LTEline

        ! Angle average, append MRC file, append MRC file intensity
        ! part, write contribution function file, write tau=1 heights
        ! file, store partial solutions, store partial solutions of
        ! intensity, correction of rho00 due to the change in J00
        ! going from non-magnetic multilevel to magnetic multiterm,
        ! storing intensity solution, apply NG acceleration, store
        ! background continuum in file, store damping in file, store
        ! inelastic collisions in file, numerical magnetic field,
        ! save parameters of damping, add bound-bound background
        ! transitions, memoization of J symbols, protect hydrogen
        ! from the equation of state, apply NG acceleration to
        ! intensity, if RAM use should be reported, protect all
        ! atoms in chemical equilibrium, measure performance in
        ! blocks, measure performance per CPU, keep a file with the
        ! frequency dependent JKQ, if using coherent wings
        ! approximation, if keeping solution files, keep populations,
        ! keep departure coefficients, keep output rhoKQ, keep
        ! output JKQ, keep stokes in quadrature, keep MRC, skip the
        ! disk in CLE, if input spectrum loaded, restrict in tau_c,
        ! restrict in z, angle-averaged forced in intensity problem,
        ! force intensity problem to be static, if there is a file for
        ! weights, keep collisions log, keep MPI log, keep, MPI
        ! detailed log, if polarization with magnetic field must be
        ! done in two steps, if excluded pixels, if truncating tau or
        ! height restriction, force ALI iterations, initialize
        ! radiation field with bound-bound transitions, keep elastic
        ! rates, add the continuum to RT coefficients in CLE, use the
        ! Allen quantities for incoming intensity in CLE, assume flat
        ! spectrum when computing input JKQ in CLE if no input
        ! spectra, if restricting height for redistribution, if
        ! restricting tau for redistribution, if considering only K=2
        ! for MRC, if consider ALI for photoionization transitions,
        ! allow switching off ALI if SEE gives negative populations
        logical:: AV,appendMRC,appendMRCI,out_contr,out_tau1,store, &
                  storeI,Pcorr,Raman,keepIsol,NG,keep_back, &
                  keep_damp,keep_cols,bfieldn,keep_aparam,addbb, &
                  keep_atmo,memo,protect_H,NGI,RAMreport, &
                  chem_protect_all,asym,g_perf,mpi_perf,keep_jkqnu, &
                  cohw,cohwi,keep_sol,keep_pop,keep_dep,keep_rhoKQ, &
                  keep_JKQ,keep_stokesQ,keep_MRC,skip_disk, &
                  lspect_input,rest_tau,rest_z,AVI,static_int, &
                  linv_weight,keep_coll,keep_mpil,keep_mpidl, &
                  two_step_pol,lexcl,rest_tau_strc,rest_z_strc, &
                  ALI_force,init_J_bb,keep_qel,add_cont_cle, &
                  use_allen,flat_cle_in,rest_z_red,rest_tau_red, &
                  anisotropy_only,ALI_photo,ALI_allow_off

        ! If asymmetry input
        logical, dimension(2):: lasym

        ! Type of Doppler width input to build frequency axis
        character(len=3) dws

        ! Mode of solution switch, force type of problem (I or Stokes)
        ! heights, vertical scale for model atmosphere
        character(len=1) mode,force,atm_scale

        ! Run ID
        character(len=9) ID

        ! output folder, atmospheric file, continuum file,
        ! magnetic field files, solution file, file with fudge
        ! factors, input file name, cache file, partition function
        ! file name, abundance file name, barklem file with sp
        ! data, barklem file with pd data, barklem file with df data,
        ! path to CHIANTI database, file with CLE input spectra,
        ! file with weights for the inversion
        character(len=500) folder,atmo,continuum,bfield,solution, &
                           fudge,source,resource,input,cache,pf, &
                           abund,bark_sp,bark_pd,bark_df, &
                           chianti_path,spect_input,inv_weight

        ! Name of atomic files, population files, background atom
        ! files, background atom population files, molecules,
        ! Kurucz line files, wavelength files, asymmetry files
        type(strarr_class), dimension(:), allocatable:: atom,popu, &
                                                        atomback, &
                                                        popuback, &
                                                        mol,kurucz, &
                                                        waves, &
                                                        asym_fil, &
                                                        fwhm_fil

        ! Ionization fraction data
        type(strnum_class), dimension(:), allocatable:: ionf

        ! Force no magnetic field, force no velocity, to force
        ! observed frequencies in synthesis axis
        logical:: unmagnetized,static,force_inv_freq

        ! Keep atomic populations fixed, zero out the last ion
        ! populations, skip this atom for the wavelength axis,
        ! fix populations for lower term
        logical, dimension(:), allocatable:: fixp,zero_ion, &
                                             skip_wave,fixplt

        ! Number of first iteration, number of maximum iteration,
        ! order of iteration (emissivity), number of steps between
        ! saving data, number of atoms, number of background atoms,
        ! number of molecules, number of first iteration in intensity,
        ! number of maximum iteration in intensity, number of steps
        ! between saving data in intensity, maximum number of
        ! internal PRD iterations in intensity, maximum number of
        ! radiation field initial iterations, order of the NG
        ! acceleration, delay in iterations before applying NG
        ! acceleration, first iteration to apply ALI to, number of
        ! Kurucz line files, maximum iteration to allow non physical
        ! quantities in stokes and rho, number of wavelength files,
        ! mode of Zeeman effect, update atmospheric model at the end
        ! of the calculation for intensity, recompute electron
        ! density, order of the NG acceleration for intensity, delay
        ! in iterations before applying NG acceleration for intensity,
        ! number of asymmetry input given by constant numbers, number
        ! of asymmetry inputs given by files, number of asymmetry
        ! inputs, if need to take care of magnetically induced
        ! transitions, mode of running (synthesis 1D, syn. 1.5D, or
        ! inversion), number of CPU groups to split columns in 1.5D or
        ! inversion, type of presurre/density scale, polar nodes,
        ! azimuthal nodes, intensity polar nodes, intensity azimuthal
        ! nodes, AA integral nodes, intensity AA integral nodes, LOS
        ! polar directions, LOS azimuthal directions, number of LTE
        ! lines, number of pixels to exclude, for how many iterations
        ! allow negative populations, type of interpolation for the
        ! PRD to transform into the observer reference frame, number
        ! of internal PRD iterations
        integer:: iter_min,iter_max,iter_ord,store_step,nA,nAb,nM, &
                  iteri_min,iteri_max,storei_step,iteri_prd,iter_j, &
                  NG_ord,NG_delay,ALI_delay,NK,allownphys_stk, &
                  allownphys_rho,NW,zeeman_mode,update_atmos, &
                  redo_ne,NGI_ord,NGI_delay,PRD_delay,nasym_num, &
                  nasym_fil,nasym,MIT_input,run_mode,rt_group_n, &
                  atmo_char,nTh,nPh,nThI,nPhI,nThAA,nThAAI,nThLOS, &
                  nPhLOS,nLTE,nexcl,allownphys_pop,PRD_int_mode, &
                  iter_prd

        ! Box to solve in 1.5D synthesis problem
        integer, dimension(:), allocatable:: sol_box

        ! Input for additional K cuts, excluded pixels
        integer, dimension(:,:), allocatable:: Kcut_input,excl

        ! Value of the Doppler width to build the frequency axis,
        ! factor for the nodes dedicated to magnetically induced
        ! transitions, doppler widths for coherent wings,
        ! minimum expected temperature (or minimum temperature
        ! depending on inputs) and maximum temperature,
        ! maximum expected velocity (or maximum velocity, depending
        ! on inputs), reference wavelength for tau in model
        ! atmosphere, effective temperature for CLE radiation, radius
        ! of the star for CLE, minimum tauc to consider, maximum tauc
        ! to consider, minimum height to consider, maximum height to
        ! consider, forced microturbulence, maximum tau continuum to
        ! calculate PRD, minimim height to calculate PRD
        double precision:: dw,MIT_node,dcohw,dcohwi,minT,maxT,maxV, &
                           omega_ref,T_rad,R_star,r0tc,r1tc,r0z,r1z, &
                           fvmicro,r1tc_prd,r1z_prd

        ! LOS polar mus, LOS azimuthal angles
        double precision, dimension(:), allocatable:: L_mu,L_phi

        ! Parameters for redistribution input frequency axis: rang,
        ! reso, negl, vlar, fstp, mstp, core, rang_core, fstp_core,
        ! mstp_core, for polarization and for intensity, respectively
        double precision, dimension(11):: red_pars,redi_pars

        ! MRC for rho00, MRC for rhoKQ with K!=0, MRC for rho00 in the
        ! intensity problem, MRC for J00 in the intensity problem,
        ! MRC for J00 initial iterations, MRC for J00, MRC for JKQ
        ! with K!=0
        double precision:: mrc_i,mrc_p,mrci_i,mrci_r,mrcj,mrc_r, &
                           mrc_p_r

        ! Numerical values for field
        double precision, dimension(3):: bfieldv

        ! Asymmetry numbers in input
        complex(kind=8), dimension(:,:), allocatable:: asym_num


!!!!!!!!!
!!!!!!!!! Inversion only inputs
!!!!!!!!!

        ! Nodes for the inversion variables
        type(Node_class), dimension(:), allocatable:: Node

        ! Filename with input Stokes profiles, filename of file to
        ! restore the inversion from, name of the output file,
        ! filename of file with mask
        character(len=500):: Filename_Ob,Inv_init,Output_file, &
                             Inv_mask

        ! Inversion ID
        character(len=9) IDv

        ! If Broyden method in LM, if the input is a fit file,
        ! neglect sigma avobe 3 (not sure what it is for), if
        ! automatic weights for Stokes, centered derivative, if
        ! correcting the node positions from the atmosphere, if
        ! the pressure is given at the boundary (hydrostatic
        ! equilibrium), if return fractional polarization, if
        ! project the magnetic field, if using previous solution
        ! for RF calculation keep the response functions, if JKQ
        ! (assymetries) must be in the output, if tracking the
        ! value of lambda between iterations in backtracking,
        ! if storing incomplete inversion results
        logical:: Broyden,FITSFILE,Sigma_neglect,auto_weight, &
                  centered,Pos_Correction,hydroeq,Fractional, &
                  Projection,Popuinit,Keep_RF,out_jkqa, &
                  l_Lam_track,storeinv

        ! Flag to modify variable in the inversion, flag for
        ! the regularization of each variable
        logical, dimension(:), allocatable:: Nodes_Flags,Nodes_Regul

        ! Number of inversion variables
        integer:: nvar = nvar_inv
        integer:: nvar_th = 9
        integer:: nvar_mg = 3
        integer:: nvar_as = 4
        integer:: nvar_g = 1

        ! Maximum number of iterations, type of inversion,
        ! type of error, nodes in the atmosphere during synthesis,
        ! type of LM method, indicate initialization, type of
        ! interpolation method, type of magnetic field vector,
        ! type of velocity vector type of SVD, index where the
        ! extension ends in filenames, number of the weights, type
        ! of input atmosphere, order of the lambda tracking between
        ! iterations in backtracking, Backtracking mode when stuck,
        ! number of steps between saving inversion data
        integer:: Num_Iter,Type_Inversion,Err_Type,Atmo_Input, &
                  LM_Method,Init_Thermal,Interpolation,btype,vtype, &
                  SVD_type,fits_index,Num_Weight,atmoin_type, &
                  Lam_track,LM_Back_Mode,storeinv_step

        ! Output file sizes (int)
        integer:: s_inv_h,s_inv_atmo_c,s_inv_res_h,s_inv_res_c, &
                  s_inv_RF_h,s_inv_RF_c

        ! Output file sizes (double)
        double precision:: s_inv_atmo,s_inv_res

        ! Units to direct the verbosity
        integer, dimension(3):: Unit_VB

        ! Type of node value, number of nodes, index of the
        ! regularization for each variable the regulatization
        ! for each variable
        integer, dimension(:), allocatable:: Node_Type,Num_nodes, &
                                             Indx_regul

        ! Threshold in chi2, threshold for the fractional chi2,
        ! ratio limit for the regulatizations with respect to the
        ! proper chi2, threshold for the SVD, boundary value for
        ! the pressure, maximum step allowed in SVD, initial Bpos or
        ! B theta if input too small, initial B azimuth if input too
        ! small, factor to decrease lambda when LM iteration
        ! accepted, factor to increase lambda when LM iteration
        ! rejected, diffuse light factor, initial vz or vpos if
        ! input too small, initial vy or azimuth if input too small,
        ! Factor to regularization
        double precision:: Threshold_chisq,Chisq_fraction,  &
                           Regul_Limit,Threshold_svd,Pg_bound, &
                           Max_Step,ini_Bpos,ini_Bazi,factoraccept, &
                           factorreject,f_diff,ini_vpos,ini_vazi, &
                           LM_lam_big_test,LM_lam_small_test, &
                           LM_lam_big_prove,LM_lam_small_prove, &
                           Regul_factor

        ! scale for each parameter, perturbation for each parameter,
        ! minimum relative perturbation
        double precision, dimension(:), allocatable:: Scal,Perturb, &
                                                      min_rel_Pert

        ! Tau ranges to consider, LM lambda ranges to consider
        double precision, dimension(2):: Tau_Range,Lam_Range

        ! Weight of the regularization function for each variable
        double precision, dimension(:), allocatable:: Regul_weight

        ! Made-up stratification from inputs
        double precision, dimension(:), allocatable:: Atmo_strat_done

        ! Weights for each Stokes
        double precision, dimension(:,:), allocatable:: Weight

        ! Data to specify atmospheric stratification modifications
        double precision, dimension(:,:), allocatable:: Atmo_strat

        ! Additional weight factor, additional sigma factor
        double precision, dimension(:,:), allocatable:: &
                                                      Weight_Factor, &
                                                      Sigma_Factor
      end type Input_class

!#####################################################################

      !> Opacity fudge data
      type fudge_class

        ! Number of frequencies with data
        integer:: nfreq_f

        ! Fudge factor data
        double precision, dimension(:,:), allocatable:: fudge_v

      end type fudge_class

!#####################################################################

      !> Background opacities
      type Continuum_class

        ! Continuum presence and Angle dependence
        logical:: d

        ! Number of directions
        integer:: ndir

        ! Continuum absorption, scattering and emissivity
        ! (eta, sig, eps)
        double precision, dimension(:,:,:,:), allocatable:: c

      end type Continuum_class

!#####################################################################

      !> Maximum relative change data
      type MRC_class

        ! First index:
        ! 1:ia, 2:iz, 3:iterm, 4:iJ, 5:iJ1, 6(1):K, 6(2):Q
        ! Second index (only for first 1:5):
        ! 1:population, 2:polarization
        integer, dimension(6,2):: indexes

        ! First index: 1:z, 2:MRC
        ! Second index: 1:population, 2:polarization
        double precision, dimension(2,2):: values

      end type MRC_class

!#####################################################################

      !> MPI data
      type MPI_class

        ! Are we doing MPI, are we splitting frequencies,
        ! use alternative solverI, use alternative solver, use
        ! alternative solverJ, use alternative solver JKQgen,
        ! if sending columns to slaves
        logical:: mpi,alternI,alternP,alternJ,alternJgen,mpi15d

        ! Number of processors, identifier, number of
        ! slave processors, global number of processors, global
        ! identifier, number of groups to send columns
        integer:: nproc,pid,nnd,gnproc,gpid,ngroup

        ! Integers for non-blocking transfer
#ifdef oldmpi
        integer:: &
#else
        type(MPI_request):: &
#endif
            request1,request2,request3,request4,request5,request6, &
            request7,request8,request9,request0,request11

        ! Maximum number of frequencies per processor,
        ! maximum number of profile size per processor, same for
        ! intensity, same for photoionizations
        integer:: nxfreq,nxtfreq,nxtfreqi,nxpfreq

        ! Number of frequencies for the processor and its limits,
        ! sizes of packages to transfer, list of processors to send in
        ! custom bcast, leaders of the slave groups, number of
        ! frequencies for the processor and its limits (for the input
        ! in CLE)
        integer, dimension(:), allocatable:: nf,if0,if1, &
                                             size1,size2,size3, &
                                             size4,size5,size6, &
                                             size7,size8,&
                                             sizei1,sizei2,sizei3, &
                                             sizei4,sizei5,sizei6, &
                                             sizei7,sizei8,sizei9, &
                                             sizei0,size10,&
                                             sizei10,sizei11, &
                                             sizei12,sizei13, &
                                             sizei14,sizei4b, &
                                             ltslave,inf,iif0,iif1, &
                                             CMf0,CMf1

      end type MPI_class

!#####################################################################

      !> Memoization

      type scalar
         double precision, pointer :: d
      end type scalar

      type a1D
         type(scalar), dimension(:), pointer :: d
      end type a1D

      type a2D
         type(a1D), dimension(:), pointer :: d
      end type a2D

      type a3D
         type(a2D), dimension(:), pointer :: d
      end type a3D

      type a4D
         type(a3D), dimension(:), pointer :: d
      end type a4D

      type a5D
         type(a4D), dimension(:), pointer :: d
      end type a5D

      type a6D
         type(a5D), dimension(:), pointer :: d
      end type a6D

      type a7D
         type(a6D), dimension(:), pointer :: d
      end type a7D

      type a8D
         type(a7D), dimension(:), pointer :: d
      end type a8D

      type a9D
         type(a8D), dimension(:), pointer :: d
      end type a9D

!#####################################################################

      !> Factorials, signs, and J-symbols
      type Fctsg_class

        ! Doing memoization with jagged arrays
        logical:: memo

        ! Memoization of J symbols
        type(a6D) :: J6
        type(a6D) :: J3
        type(a9D) :: J9

        ! Factorial and sign
        double precision, dimension(:), allocatable:: flg,sg

      end type Fctsg_class

!#####################################################################

      !> Frequency data
      type Frequency_class

        ! Minimum and maximum index with photoionizations and lines,
        ! Number of polar and azimuthal directions, total directions,
        ! total directions that are in the quadrature, size of
        ! frequency space for profile messages, same for intensity,
        ! same for photoionizations
        integer:: pif0,pif1,lif0,lif1,nth,nph,ndir,nqdir,ntfreq, &
                  ntfreqi,npfreq

        ! Minimum and maximum index with photoionizations for the
        ! master, and for the lines, weight for each frequency node
        ! for sharing tasks, size of frequency space for profile
        ! messages, same for intensity, weight for each frequency
        ! node for sharing task but neglecting PRD, mapping of
        ! output frequencies (CLE) into general axis
        integer, dimension(:), allocatable:: Mpif0,Mpif1,Mlif0, &
                                             Mlif1,IW_freq,Mntfreq, &
                                             Mntfreqi,Mnpfreq, &
                                             IW_freq_in,mapping

        ! Frequency, weight (output), frequency to the cube
        double precision, dimension(:), allocatable:: omega, &
                                                      W_freq, &
                                                      omega3, &
                                                      omega_ou, &
                                                      omega3_ou

        ! Exponentials for epsIphoto
        double precision, dimension(:,:), pointer:: exu

      end type Frequency_class

!#####################################################################

      !> Redistribution functions
      type Redc2_class

        ! If Wfunc2 has to be calculated, if Wfunc2 can be stored
        logical:: iIPRD = .True., RAM

        ! If storing redistribution
        logical, dimension(:), allocatable:: iPPRD

        ! Redistribution (real)
        real, dimension(:), allocatable:: IWarr2

        ! Redistribution (complex)
        complex(kind=4), dimension(:), allocatable:: PWarr2

      end type Redc2_class

!#####################################################################

      !> Structure with input transitions for redistribution
      type Redb2_class

        ! Input transitions
        type(Redc2_class), dimension(:), pointer:: trani

      end type Redb2_class

!#####################################################################

      !> Input frequency axes and weights
      type Redc_class

        ! Size of omega and W_freq, size of interpolation data
        integer:: osize,isize

        ! Input frequency size
        integer, dimension(:), allocatable:: mfreq

        ! Frequencies and weights
        double precision, dimension(:), allocatable:: omega,W_freq

      end type Redc_class

!#####################################################################

      !> Output frequency ranges for redistribution, input transitions
      !! for frequency axes and weights, and second order RT
      !! coefficients
      type Redb_class

        ! Max of input frequencies, number of ranges, number of
        ! output frequencies, true maximum frequency limits for
        ! 2ord, maximum frequency limits for 2ord, limits for
        ! interpolation
        integer:: mxfreq,nran,nfreq,Igf0,Igf1,tgf0,tgf1,gf0,gf1, &
                  ggf0,ggf1

        ! Number of frequencies per CPU, lower frequency limits
        ! for all CPU, upper frequency limits, for all CPU,
        ! indexes within ranges
        integer, dimension(:), allocatable:: nf,Mif0,Mif1,Rif0,Rif1, &
                                             if0,if1

        ! Input transitions
        type(Redc_class), dimension(:), pointer:: trani

        ! Emissivity and ALI correction
        double precision, dimension(:,:), allocatable:: eps20,eps21, &
                                                        eps22,eps23, &
                                                        rpf
      end type Redb_class

!#####################################################################

      !> Voigt and Faraday-Voigt profiles, and Voigt normalization
      !! factors
      type Prof_class

        ! Tell if there is a 1st order stored profile
        logical:: VRAM

        ! Normalization of the first order profiles
        double precision, dimension(:), allocatable:: Norm

        ! First order profile for intensity
        double precision, dimension(:), allocatable:: p

        ! First order profile for polarization
        complex(kind=8), dimension(:,:), allocatable:: cp

      end type Prof_class

!#####################################################################

      !> Redistribution input frequency data, redistribution function
      !! data, and profile or normalization data
      type Red_class

        ! Sizes
        integer:: ndzao, nzao, ndzaoA

        ! Indexing height-atom-transition for PRD data and
        ! indexing for direction-height-atom_transition for
        ! normalization data
        integer, dimension(:,:,:), allocatable:: izao,idzao

        ! Redistribution class for height, atom, transition
        type(Prof_class), dimension(:), pointer:: dzao,pzao

        ! Redistribution class for height, atom, transition
        type(Redb_class), dimension(:), pointer:: zao

        ! Redistribution class for height, atom, transition
        type(Redb2_class), dimension(:), pointer:: rzao

      end type Red_class

!#####################################################################

      !> Temporal rhoKQ or population data
      type Rhoc_class

        ! population previous step
        double precision, dimension(:,:), allocatable:: rho

        ! rhoKQ for previous step
        complex(kind=8), dimension(:,:), allocatable:: crho

      end type Rhoc_class

!#####################################################################
!#####################################################################
!#####################################################################
!
!            TIC CLASSES
!
!#####################################################################
!#####################################################################
!#####################################################################

      !> Inversion node data
      type Nodes_class

        ! Nodes for the inversion variables
        type(Node_class), dimension(nvar_inv):: Node

        ! If regularizing, if the pressure is given at boundary,
        ! if correcting the node positions from the atmosphere,
        ! if taking the gas pressure from the input model, if
        ! hydrostatic equilibrium must be imposed each time
        logical:: Regul_Flag,hydroeq,Pg_Inv,Pos_Correction,Pg_Auto, &
                  hydros

        ! Flag to modify variable in the inversion, flag for
        ! the regularization of each variable
        logical, dimension(nvar_inv):: Nodes_Flags,Nodes_Regul

        ! Index of inversion variables
        integer:: nvar = nvar_inv
        integer:: index_B = 1
        integer:: index_Bt = 2
        integer:: index_Bp = 3
        integer:: index_f = 4
        integer:: index_T = 5
        integer:: index_vx = 6
        integer:: index_vy = 7
        integer:: index_vz = 8
        integer:: index_vm = 9
        integer:: index_Pg = 10
        integer:: index_J21R = 11
        integer:: index_J21I = 12
        integer:: index_J22R = 13
        integer:: index_J22I = 14

        ! Type of node location, Type of Magnetic field vector, type
        ! of velocity field vector, lower index of variables to
        ! consider, upper index of variables to consider, total number
        ! of nodes, number of nodes for the inversion, number of
        ! magnetic parameter nodes for the inversion, number of
        ! thermal parameters nodes for the inversion, number of
        ! ad-hoc assymmetries nodes for the inversion, type of
        ! interpolation method, Type of inversion, number of global
        ! nodes for the inversion
        integer:: Node_Location_Type,Btype,vtype,Indx_b,Indx_e, &
                  Tot_Nodes,Num_Fit,Num_Mag,Num_Thermal, &
                  Num_Asymmetry,Interpolation,Nodes_Type,Num_glob

        ! Type of node value for each variable, number of nodes for
        ! each variable, number of nodes that can change for each
        ! variable, index of the regularization for each variable,
        ! final number the regulatization "points"
        integer, dimension(nvar_inv):: Node_Type,Num_Nodes,Num_Vary, &
                                       Indx_regul,Num_regul

        ! Indexes of the first and last nodes that can change for each
        ! variable
        integer, dimension(2,nvar_inv):: Node_Vary

        ! Index of nodes and parameters
        integer, dimension(:,:), allocatable:: Inf_Inv

        ! Threshold for the SVD, maximum step allowed in SVD, boundary
        ! value for the pressure, cosine of the heliocentric angle for
        ! the emergence, azimuth for emergence
        double precision:: Threshold_svd,Max_Step,Pg_Bound,mu,azimuth

        ! Weight of the regularization function for each variable,
        ! value of the perturbations to the model, scale value for
        ! each parameter, value for constant regularization, minumum
        ! relative perturbation
        double precision, dimension(nvar_inv):: Regul_weight, &
                                                Perturb,Scal,Const, &
                                                min_rel_Pert

      end type Nodes_class

!#####################################################################

      !> Frequency and synthetic Stokes parameters in the frequency
      !! range of the inverted data
      type Solution_class

        !
        !type(Sol_class):: Sol_Tmp, Sol_Min, Sol_Sav

        ! If can issue warning for wrong memory count
        logical:: warning

        ! If interpolating frequencies, if return fractional
        ! polarization, if project the magnetic field, if
        ! diffuse light
        logical:: Fre_Intp,Fractional,Projection,Diff_flag

        ! Number of wavelength points, number of wavelength ranges
        integer:: Num_Wavelength,Num_Range

        ! Wavelength ranges
        integer, dimension(:,:), allocatable:: Range

        ! Input frequency axis
        double precision, dimension(:), allocatable:: omega_input

        ! Emergent Stokes and diffuse light Stokes (for RF),
        double precision, dimension(:,:), allocatable:: Stokes_out, &
                                                        Stokes_diff

        ! Stokes scale
        double precision, dimension(:), allocatable:: Scal_Stokes

      end type Solution_class

!#####################################################################

      !> Solution of the self-consistent problem and the corresponding
      !! emergent profiles, contribution function, and height for
      !! optical depth equal to one
      type Solution_F_class

        ! Flag for Hanle and initialization
        logical:: keep_solution, no_initialized

        ! Arrays of solutions
        real, dimension(:,:,:), allocatable:: e_tau1
        real, dimension(:,:,:,:,:), allocatable:: e_Ctr
        double precision, dimension(:,:), allocatable:: i_J00
        double precision, dimension(:,:), allocatable:: i_J00C
        double precision, dimension(:,:,:), allocatable:: i_J00P
        double precision, dimension(:,:,:,:), allocatable:: e_Stk
        double precision, dimension(:,:,:,:), allocatable:: i_StkI
        double precision, dimension(:,:,:,:,:), allocatable:: i_Stk
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQ
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQS
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQC

        ! Arrays of best backtrack
        real, dimension(:,:,:), allocatable:: e_tau1_t
        real, dimension(:,:,:,:,:), allocatable:: e_Ctr_t
        double precision, dimension(:,:), allocatable:: i_J00_t
        double precision, dimension(:,:), allocatable:: i_J00C_t
        double precision, dimension(:,:,:), allocatable:: i_J00P_t
        double precision, dimension(:,:,:,:), allocatable:: e_Stk_t
        double precision, dimension(:,:,:,:), allocatable:: i_StkI_t
        double precision, dimension(:,:,:,:,:), allocatable:: i_Stk_t
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQ_t
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQS_t
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQC_t

        ! Arrays of best solutions
        real, dimension(:,:,:), allocatable:: e_tau1_b
        real, dimension(:,:,:,:,:), allocatable:: e_Ctr_b
        double precision, dimension(:,:), allocatable:: i_J00_b
        double precision, dimension(:,:), allocatable:: i_J00C_b
        double precision, dimension(:,:,:), allocatable:: i_J00P_b
        double precision, dimension(:,:,:,:), allocatable:: e_Stk_b
        double precision, dimension(:,:,:,:), allocatable:: i_StkI_b
        double precision, dimension(:,:,:,:,:), allocatable:: i_Stk_b
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQ_b
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQS_b
        complex(kind=8), dimension(:,:,:,:), allocatable:: i_JKQC_b

        ! Atom solution
        type(Rhoc_class), dimension(:), allocatable:: i_rhoes
        type(Rhoc_class), dimension(:), allocatable:: i_rhoes_t
        type(Rhoc_class), dimension(:), allocatable:: i_rhoes_b

      end type Solution_F_class

!#####################################################################

      !> Inversion Stokes parameters data
      type Stokes_class

        ! If the sigma is the same for all profile points, if
        ! automatic weights for Stokes, if sigmas read at the
        ! beginning, if diffuse light, if diffuse light read at
        ! the beginning
        logical:: Sigma_flag,Auto_Weight,Sigma_ct,Diff_flag,Diff_ct

        ! Number of wavelengths, type of sigma in the input, number
        ! of wavelength ranges, total number of data points,
        ! intensity data points
        integer:: Num_Wavelength,Indx_Sigma,Num_Range,Num_freedom, &
                  Num_freedomI

        ! Wavelength ranges
        integer, dimension(:,:), allocatable:: Range

        ! Cosine of the heliocentric angle for the observation, and
        ! azimuth of the LOS
        double precision:: mu,azimuth

        ! Scale for each Stokes,
        double precision, dimension(:,:), allocatable:: Scales

        ! Observed Stokes parameters, frequency dependent sigma,
        ! Input frequency dependent sigma, diffuse light, input
        ! diffuse light
        double precision, dimension(:,:), allocatable:: Stokes_Ob, &
                                                        Sigma_W, &
                                                        Sigma_in, &
                                                        Diff_in

        ! Weights for each Stokes
        double precision, dimension(:,:), allocatable:: Weight

      end type Stokes_class

!#####################################################################

      !> Regularization data
      type Regul_class

        ! Scale of the penalty, regulatization penalty
        double precision:: Ratio, Penalty

        ! Regularization vector
        double precision, dimension(:), allocatable:: Regul_F

        ! Regularization matrix
        double precision, dimension(:,:), allocatable:: Regul_H

      end type Regul_class

!#####################################################################

      !> Data for the Levenberg–Marquardt
      type LMFIT_class

        ! Regularizations
        type(Regul_class):: Rgl

        ! If we need to calculate weights, if we need to compute the
        ! Jacobian, if the inversion step was accepted
        logical:: Flag_weight,Flag_Jac,accepted

        ! Number of Jacobian elements
        integer:: Num

        ! Factor to reduce lambda when accepted, factor to enhance
        ! lambda when rejected, parameter for the Levenberg-Marquardt,
        ! current chi^2, older chi^2, initial chi^2
        double precision:: factoraccept,factorreject,Lambda,Chisq, &
                           Chisq_og,Chisq_0

        ! Limits for the Levenberg-Marquardt lambda parameter
        double precision, dimension(2):: Lambda_bounds

        ! Residual for the intensity, weight for the L2 of the
        ! intensity, Jacobian vector, old Jacobian vector, diagonal
        ! of the Hessian
        double precision, dimension(:), allocatable:: ResidualI, &
                                                      WeightI, &
                                                      Jacfvec, &
                                                      Jacfvec_og, &
                                                      Diag

        ! Residual for Stokes parameters, weight for the L2 of the
        ! Stokes parameters, Hessian, old Hessian, Jacobian for the
        ! intensity
        double precision, dimension(:,:), allocatable:: Residual, &
                                                        Weight, &
                                                        Hessian, &
                                                        Hessian_og, &
                                                        JacobianI

        ! Jacobian for the Stokes parameters
        double precision, dimension(:,:,:), allocatable:: Jacobian

      end type LMFIT_class

!#####################################################################

      end module types_mod
