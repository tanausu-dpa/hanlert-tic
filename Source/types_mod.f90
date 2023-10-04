      !> Class definitions
      module types_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Hao Li (IAC)
!     Roberto Casini (HAO)
!  Contributors:
!     John Dennis (NCAR)
!  Start:
!     04/17/2017
!  Last version:
!     09/29/2023 V3.0.22
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:   V3.0.22 - Added Kcut and Krad to Atom_class (TdPA)
!                           - Added Kcut_input, keep_raml, keep_mpil,
!                             and keep_mpidl to Input_class (TdPA)
!
!     09/25/2023:   V3.0.21 - Added file_label to Atom_class (TdPA)
!                           - Changed size of vfile and wfile in
!                             Atom_class (TdPA)
!
!     09/18/2023:   V3.0.20 - Removed verbosity in Input_class (TdPA)
!
!     08/24/2023:   V3.0.19 - Added FWHM_helper_class (TdPA)
!                           - Changed the type of lim_fwhm in the
!                             Input_class from IO_helper_class to
!                             FWHM_helper_class (TdPA)
!                           - Added fwhm_fil and force_inv_freq to
!                             Input_class (TdPA)
!
!     08/11/2023:   V3.0.18 - Added linv_weight and inv_weight to
!                             Input_class (TdPA)
!                           - Added Num_freedom and Num_freedomI from
!                             Stokes_class (TdPA)
!                           - Removed Num_Stok in Stokes_class (TdPA)
!
!     08/07/2023:   V3.0.17 - Added LTEprof_sub_class, LTEprof_class,
!                             and LTEline_class (TdPA)
!                           - Added LTEline and nLTE to the
!                             Input_class structure (TdPA)
!
!     07/06/2023:   V3.0.16 - Added Atmo_strat_done to the type
!                             Input_class (TdPA)
!
!     07/03/2023:   V3.0.15 - Updated the number of variables (TdPA)
!                           - Removed module variables about
!                             specific subsets of variables (TdPA)
!                           - Removed module variables about
!                           - Added fixpt to Atom_class (TdPA)
!                           - Added alloc_a, alloc_b, f_diff, vlos,
!                             vpos, vphi, and chi500 to
!                             Atmo_class (TdPA)
!                           - Changed ne, nh, and nhe in Atmo_class
!                             from pointers to allocatables (TdPA)
!                           - Removed H_max and H_min from
!                             Node_class (TdPA)
!                           - Added lim_fwhm, skip_wave, fixplt,
!                             Inv_init, hydroeq, out_jkqa, nvar_g,
!                             vtype, atmoin_type, s_inv_h, s_inv_atmo,
!                             s_inv_atmo_c, s_inv_res_h, s_inv_res_c,
!                             s_inv_res, s_inv_RF_h, s_inv_RF_c,
!                             f_diff, ini_vpos, ini_vazi, and
!                             min_rel_Pert to Input_class (TdPA)
!                           - Removed Restore_File, Pg_Flag, Pg_inv,
!                             Num_File, Ind, Initpixel, Restartpixel,
!                             and FWHM to Input_class (TdPA)
!                           - Removed chi500 from
!                             Continuum_class (TdPA)
!                           - Added rho to Rhoc_class (TdPA)
!                           - Removbed Pg_Flag, Pg_Inv, and H_min from
!                             Nodes_class (TdPA)
!                           - Added hydros, index_B, index_Bt,
!                             index_Bp, indef_f, index_T, index_vx,
!                             index_vy, index_vz, index_vm, index_Pg,
!                             index_J21R, index_J21I, index_J22R,
!                             index_J22I, vtype, Num_glob, azimuth,
!                             and min_rel_Pert to Nodes_class (TdPA)
!                           - Removed Stokes_Ref and FWHM from
!                             Solution_class (TdPA)
!                           - Added Diff_flag and Stokes_diff to
!                             Solution_class (TdPA)
!                           - Added Solution_F_class (TdPA)
!                           - Added Sigma_ct, Diff_flag, Diff_ct,
!                             azimuth, Sigma_in, and Diff_in to
!                             Stokes_class (TdPA)
!                           - Added chisq_0 to LMFIT_class (TdPA)
!                           - Undid most of 3.0.14 changes because
!                             merging was impossible and some of
!                             the functionalities were replicated
!                             with existing keywords (TdPA)
!
!     06/12/2023:   V3.0.14 - Added file_size_class (HL)
!                           - update Restartpixel for 2 dimensions(HL)
!                           - added Num_FWHM, filesize, region, ix,
!                             iy, nx, ny, indxmu, indxstk, len_pixel,
!                             and len_header to Input_class (HL)
!                           - removed projection, Num_File, Ind (HL)
!                           - FWHM changed to an allocable array (HL)
!                           - Removed redundant Sol_class (HL)
!
!     04/25/2023:   V3.0.13 - Added Keep_RF to Input_class (TdPA)
!
!     04/11/2023:   V3.0.12 - Node_class should be defined before
!                             Input_class (HL)
!                           - Remove the keywords Hanle_Effect and
!                             CRD_RF (HL)
!                           - Update the weights, scale, and range
!                             for multi-wavelength ranges (HL)
!
!     03/15/2023:   V3.0.11 - Removed Bfieldlos_class (TdPA)
!                           - Added Blos, Bpos, and Azimuth to
!                             Bfield_class (TdPA)
!                           - Removed Pop_Init, CRD_Rf, Fre_Intp, and
!                             Bounds from Input_class (TdPA)
!                           - Renamed Bpos to ini_Bpos and Azimuth to
!                             ini_Bazi in Input_class (TdPA)
!                           - Added factoraccept, factorreject, Scal,
!                             Perturb, Lam_Range, and Atmo_strat to
!                             Input_class (TdPA)
!                           - Added nebound and ebound to Node_class
!                             and moved Bounds from Nodes_class to
!                             Node_class (TdPA)
!                           - Removed Output_flag and Pop_Init from
!                             Solution_class (TdPA)
!
!     03/08/2023:   V3.0.10 - Added Bfieldlos_class, Node_class,
!                             Nodes_class, Sol_class, Solution_class,
!                             Stokes_class, Regul_class, and
!                             LMFIT_class (TdPA)
!                           - Added nvar_inv, nvar_th_inv,
!                             nvar_mg_inv, and nvar_as_inv global
!                             header parameters (TdPA)
!                           - Added JKQin in Atmo_class (TdPA)
!                           - Added inversion inputs to Input_class
!                             structure (TdPA)
!
!     02/14/2023:    V3.0.9 - Added vxa, vya, and vza to
!                             Atmo_class (TdPA)
!                           - Added AVI, static_int, nTh, nPh, nThI,
!                             nPhI, nThAA, nThAAI, nThLOS, nPhLOS,
!                             L_mu, L_phi, and redi_pars to
!                             Input_class (TdPA)
!
!     11/24/2022:    V3.0.8 - Added bbspecin, bfspecin, sbif0, sbif1,
!                             sfif0, sfif1, ilf0, ilf1, ipf0, and
!                             ipf1 to Atom_class (TdPA)
!                           - Added V_gauss, W_gauss, and V_mu_disk to
!                             Geometry_class (TdPA)
!                           - Added geom_size to the IO_helper_class
!                             structure (TdPA)
!                           - Added inf, iif0, and iif1 to the
!                             MPI_class (TdPA)
!                           - Added IW_freq_in, mapping, omega_ou, and
!                             omega3_ou to Frequency_class (TdPA)
!                           - Added spect_class,
!                             chianti_real_pointer_class,
!                             chianti_ioneq_class, chianti_class, and
!                             Coronapoint_class (TdPA)
!
!     11/10/2022:    V3.0.7 - Added zero_ion to Atom_class and to
!                             Input_class (TdPA)
!                           - Removed stm from Input_class (TdPA)
!
!     10/26/2022:    V3.0.6 - Added rdip_class and changed the storage
!                             structure from Atom%rdip to
!                             Atom%rdip%rdip (TdPA)
!                           - Added irho_class and Jrho_class and
!                             changed the storage structure of
!                             Atom%irho and Atom%irho_ij using them
!                             to save memory in big models (TdPA)
!
!     10/25/2022:    V3.0.5 - Added d0, di0, di1, it, imi, iv, ine,
!                             inh, ipe, ipg, irh, ib, mode, norm,
!                             ypos, and zpos to Atmo_class (TdPA)
!                           - Added strnum_class (TdPA)
!                           - Added skip_disk, lspect_input,
!                             rest_tau, rest_z, chianti_path,
!                             spect_input, ionf, sol_box, T_rad,
!                             R_star, r0tc, r1tc, r0z, and r1z to
!                             Input_class (TdPA)
!
!     07/27/2022:    V3.0.4 - Removed ierr from MPI_class (TdPA)
!                           - Changed type of request variables for
!                             newer MPI versions, this can be chosen
!                             in the configure script (TdPA)
!
!     07/18/2022:    V3.0.3 - Added IWskip, MPIdetail, and operform
!                             to Input_class (TdPA)
!
!     07/13/2022:    V3.0.2 - Added pf_class (TdPA)
!                           - Added ele, NT, nele, pT, and abund to
!                             Atmo_class (TdPA)
!                           - Added pf, abund, bark_sp, bark_pd, and
!                             bark_df to Input_class (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: I missed the variable lim_atmo
!                             in Input_class (TdPA)
!                           - In IO_helper_class, doub was the wrong
!                             type (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!                           - Added Bx, By, Bz, and zeros pointers in
!                             Atmo_class (TdPA)
!                           - Changed v(:,:) into vx(:), vy(:), and
!                             vz(:) pointers in Atmo_class (TdPA)
!                           - z, T, vmi, ne, nh, and nhe are pointers
!                             now in Atmo_class (TdPA)
!                           - Reduces in two dimenions TSL and TBL
!                             in Geometry_class (TdPA)
!                           - Added IO_helper_class (TdPA)
!                           - Added lim_stk, lim_ctr, lim_tau,
!                             lim_cols_tt, lim_cols_ll, lim_damp,
!                             lim_back, lim_pop, keep_sol, keep_pop,
!                             keep_dep, keep_rhoKQ, keep_JKQ,
!                             keep_stokesQ, keep_MRC, atm_scale,
!                             cache, unmagnetized, static, run_mode,
!                             rt_group_n, atmo_char, minT, maxT,
!                             maxV, and omega_ref to
!                             Input_class (TdPA)
!                           - Added fudge_class (TdPA)
!                           - Added mpi15d, gnproc, gpid, ngroup,
!                             and ltslave to MPI_class (TdPA)
!
!     06/21/2022:    V2.0.2 - Added cohw and dcohw to the type
!                             Input_class (TdPA)
!
!     04/07/2022:    V2.0.1 - Added keep_jkqnu into the type
!                             Input_class (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added omp_1c_class and omp_c2_class
!                             types containing how OpenMP splits
!                             the magnetic RT coefficients (TdPA)
!                           - Added omp_1c, omp_c2, and omp_comp_1ord
!                             to the Atom_class structure, and
!                             removed rho from it (TdPA)
!                           - Removed domain decomposition related
!                             variables from Input_class and
!                             MPI_class (TdPA)
!                           - Added can_write to Fctsg_class (TdPA)
!                           - Added oif0 and oif1 to the structure
!                             Frequencyc2_class (TdPA)
!
!     02/12/2021:    V1.2.5 - Added g_perf and mpi_perf into the
!                             type Input_class (TdPA)
!
!     02/04/2021:    V1.2.4 - Added MIT_input and MIT_node into the
!                             type Input_class (TdPA)
!
!     01/13/2021:    V1.2.3 - Added asym, lasym, asym_fil, nasym_fil,
!                             nasym, and asym_num into the type
!                             Input_class (TdPA)
!
!     11/12/2020:    V1.2.2 - Added chem_protect_all and PRD_delay to
!                             Input_class (TdPA)
!
!     10/26/2020:    V1.2.1 - Changed index1, index2, and dx into
!                             pointers to avoid copying (TdPA)
!
!     09/11/2020:    V1.2.0 - Added iabox_class and dbabox_class
!                             classes (TdPA)
!                           - Initializing abun_mod in Atom_class
!                             to 1d0 at definition (TdPA)
!                           - Added depar to Atom_class (TdPA)
!                           - Added RAMreport to Input_class (TdPA)
!                           - Added RRAM and BRAM to MPI_class (TdPA)
!                           - Completely changed the chain of
!                             Frequency*_class and Red*_class
!                             structures (TdPA)
!
!     07/22/2020:   V1.1.23 - Added NGI, NGI_ord, and NGI_delay to
!                             Input_class (TdPA)
!
!     06/26/2020:   V1.1.22 - Added NCHLT to Atom_class (TdPA)
!
!     06/08/2020:   V1.1.21 - The %index array in Frequencyf_class has
!                             been split into index1 and index2 to
!                             speed-up emiss2ord (JD)
!
!     03/05/2020:   V1.1.20 - Added mol_protect in Atom_class (TdPA)
!                           - Added nHA in Atmo_class (TdPA)
!                           - Added redo_ne, protect_H, and
!                             update_atmos to Input_class (TdPA)
!
!     12/17/2019:   V1.1.19 - Added zeeman_mode to Input_class (TdPA)
!
!     12/10/2019:   V1.1.18 - Added structures needed by memoization
!                             as classes and in Fctsg_class (TdPA)
!
!     11/19/2019:   V1.1.17 - Added Ftrano_class (TdPA)
!                             Added Ftrano, vfile, wfile, hwifil,
!                             hwtfil, rif20, rif21, f0size, f1size,
!                             and dztwsize in Atom_class (TdPA)
!                           - Added VRAM and WRAM in MPI_class (TdPA)
!
!     11/13/2019:   V1.1.16 - Added AtomVindex_class (TdPA)
!                           - Changed Atom%Normp into a pointer (TdPA)
!                           - Added i_Vind, vfile, Mncom, hvifil,
!                             rif0, rif1, zsize, dsize, tsize, tBsize,
!                             and fsize into Atom_class (TdPA)
!                           - Changed i_geom from three dimensions
!                             to two (TdPA)
!
!     10/18/2019:   V1.1.15 - Added fixp in Atmo_class and
!                             Input_class (TdPA)
!
!     10/03/2019:   V1.1.14 - Added itrano and ntrano in structure
!                             Atom_class (TdPA)
!
!     09/26/2019:   V1.1.13 - Added classes elastic_entry_class,
!                             elastic_class, inelastic_class, and
!                             Tbox_class (TdPA)
!                           - Added inelas, elas, and Tbox in
!                             Atom_class (TdPA)
!                           - Added typo, Pg, rho, Pe, and zalt in
!                             Atmo_class (TdPA)
!                           - Added keep_atmo in Input_class (TdPA)
!
!     09/13/2019:   V1.1.12 - Added tfreq in Atmos_class (TdPA)
!                           - Added waves and NW in Input_class (TdPA)
!                           - Added chi500 in Continuum_class (TdPA)
!                           - Added ggf0 and ggf1 in
!                             Frequencyc2_class (TdPA)
!
!     08/14/2019:   V1.1.11 - Added tmp_col_box_class type (TdPA)
!                           - Added Ccoeff_special variable of type
!                             tmp_col_box_class in Atom_class (TdPA)
!
!     08/09/2019:   V1.1.10 - Added addbb to Input_class (TdPA)
!
!     06/04/2019:    V1.1.9 - Added ilevell and ilevelu to the
!                             FST_class structure (TdPA)
!
!     06/03/2019:    V1.1.8 - Added splitf to Atom_class (TdPA)
!
!     05/31/2019:    V1.1.7 - Added nxtfreq, nxtfreqi, and nxpfreq to
!                             MPI_class (TdPA)
!                           - Added the variables ntfreq, ntfreqi,
!                             npfreq, Mntfreq, Mntfreqi, and Mnpfreq
!                             to Frequency_class (TdPA)
!
!     05/08/2019:    V1.1.6 - Added tshift, pshift, and tfshift in
!                             Atom_class (TdPA)
!                           - Added allownphys_stk and allownphys_rho
!                             in Input_class (TdPA)
!                           - Added njdir in Red_class (TdPA)
!
!     04/15/2019:    V1.1.5 - Added altbcast in input_class and
!                             mpi_class (TdPA)
!
!     03/22/2019:    V1.1.4 - Added ckurucz_class and kurucz_class
!                             for background lines (TdPA)
!
!     03/18/2019:    V1.1.3 - Added iterml and itermu to FST_class
!                             and ilevell and ilevelu to
!                             Phot_class (TdPA)
!
!     03/13/2019:    V1.1.2 - Added ALI_delay in Input_class (TdPA)
!
!     03/12/2019:    V1.1.1 - Added keep_aparam in Input_class (TdPA)
!                           - Removed iexu from Frequency_class (TdPA)
!
!     02/20/2019:    V1.1.0 - Removed unit for verbosity in Input and
!                             added numerical field for magnetic field
!                             in the input file in Input_class (TdPA)
!
!     11/06/2018:   V1.0.28 - Added keep_back, keep_damp, and
!                             keep_cols to Input_class (TdPA)
!
!     09/20/2018:   V1.0.27 - Added NG, NG_ord, and NG_delay to
!                             Input_class (TdPA)
!
!     09/06/2018:   V1.0.26 - Added keepIsol to Input_class (TdPA)
!
!     09/04/2018:   V1.0.25 - Added alternI, alternP, alternJ,
!                             alternJgen, and sizei4b into
!                             MPI_class (TdPA)
!                           - Added iexu to Frequency_class (TdPA)
!
!     08/06/2018:   V1.0.24 - Added Voigt_class to store Voigt
!                             profiles (TdPA)
!                           - Added prof and RAM to Nindex_class
!                             (TdPA)
!                           - Added RAM and PRAM to MPI_class (TdPA)
!                           - Added omega3 and exu to Frequency_class
!                             (TdPA)
!                           - Added RAM to Redd_class (TdPA)
!
!     12/05/2017:   V1.0.23 - Added Raman to Input (TdPA)
!
!     11/27/2017:   V1.0.22 - Added Pcorr to Input (TdPA)
!
!     10/30/2017:   V1.0.21 - Added ML to Atom (TdPA)
!
!     10/23/2017:   V1.0.20 - Warr2 is single precision now (TdPA)
!
!     10/11/2017:   V1.0.19 - Added gL in atom (TdPA)
!
!     10/03/2017:   V1.0.18 - Added WindNB for atomic indexing (TdPA)
!
!     09/27/2017:   V1.0.17 - Added nsend, recv, lsend to MPI, for
!                             the new broadcasting (TdPA)
!
!     09/15/2017:   V1.0.16 - Added resource in Input (TdPA)
!
!     09/14/2017:   V1.0.15 - Added source and ID in Input (TdPA)
!
!     09/08/2017:   V1.0.14 - fout in frequency and fout and trani in
!                             redistribution are now pointers instead
!                             of simple allocatables (TdPA)
!
!     08/24/2017:   V1.0.13 - red_pars in Input now has 11 elements
!                             instead of 6 (TdPA)
!
!     08/22/2017:   V1.0.12 - Added sizei11, sizei12, sizei14, and
!                             sizei14 in MPI (TdPA)
!
!     08/21/2017:   V1.0.11 - Added size10 in MPI (TdPA)
!
!     07/21/2017:   V1.0.10 - omega and W_freq one step higher in
!                             Frec (TdPA)
!
!     06/28/2017:    V1.0.9 - Created Red_class and its own tree
!                             similar to Frec, exclusively for Warr
!                             storage (TdPA)
!
!     06/23/2017:    V1.0.8 - Added IW_freq to Frec (TdPA)
!
!     06/22/2017:    V1.0.7 - Added request11 to MPI (TdPA)
!
!     06/19/2017:    V1.0.6 - Changed the structure tree for input
!                             frequencies. The position of Warr2 for
!                             intensity and added Warr2 for
!                             polarization (TdPA)
!                           - Added indexing variables in Atom for
!                             the ordering of Warr2 in the
!                             polarized mode (TdPA)
!                           - Changed iPRD to iIPRD (TdPA)
!
!     06/16/2017:    V1.0.5 - Changed the structure tree for input
!                             frequencies. Before, it was tran(in,out)
!                             and now trano(out)%trani(in) (TdPA)
!                           - Added indexing of input transitions in
!                             trani of the frequency structure tree
!                             (TdPA)
!                           - Changed if0, if1, Mif0, and Mif1 in Frec
!                             to pif0, pif1, Mpif0, Mpif1, and added
!                             lif0, lif1, Mlif0, Mlif1 (TdPA)
!
!     06/13/2017:    V1.0.4 - Added variables if0, if1, Mif0, and Mif1
!                             to Frec (TdPA)
!                           - Removed interp from the Frec structure
!                             tree (TdPA)
!
!     06/12/2017:    V1.0.3 - Added variables if0, if1, W0, W1, Mif0,
!                             Mif1, MW0, and MW1 in Atom (TdPA)
!                           - Added absent in Phot (TdPA)
!                           - Removed l from Cont (TdPA)
!                           - Added Warr2 and iPRD in the Frec
!                             structure tree (TdPA)
!
!     06/08/2017:    V1.0.2 - Added variables for J and PRD
!                             iterations (TdPA)
!
!     05/04/2017:    V1.0.1 - MPI status variables no longer
!                             used (TdPA)
!
!     04/17/2017:    V1.0.0 - First version (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!    This module contains the definitions of the derived data types
!  or structures:
!
!    strarr_class
!    strarr2_class
!    Freqflag_class
!    Voigt_clalss
!    Nindex_class
!    Atomindex_class
!    Warrindexb_class
!    Warrindexa_class
!    Phot_class
!    FST_class
!    rdip_class
!    Jrho_class
!    irho_class
!    Atom_class
!    LTEprof_sub_class
!    LTEprof_class
!    LTEline_class
!    catm_class
!    Mol_class
!    spect_class
!    chianti_real_pointer_class
!    chianti_ioneq_class
!    chianti_class
!    pf_class
!    Atmo_class
!    Bfield_class
!    Coronapoint_class
!    Geometry_class
!    IO_helper_class
!    FWHM_helper_class
!    strnum_class
!    Input_class
!    fudge_class
!    Continuum_class
!    MRC_class
!    MPI_class
!    Fctsg_class
!    Frequencyf_class
!    Frequencye_class
!    Frequencyd_class
!    Frequencyc2_class
!    Frequencyc1_class
!    Frequencyb_class
!    Frequencya_class
!    Frequency_class
!    Redg_class
!    Redf_class
!    Rede_class
!    Redd_class
!    Redc2_class
!    Redc1_class
!    Redb_class
!    Reda_class
!    Red_class
!    Rhoc_class
!
!    !!!!!!!!!!!
!    TIC classes
!    !!!!!!!!!!!
!
!    Node_class
!    Nodes_class
!    Solution_class
!    Solution_F_class
!    Stokes_class
!    Regul_class
!    LMFIT_class
!
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod

      integer, parameter:: nvar_inv = 14

!#####################################################################

      !> Structure to temporally contain collisional rates
      !> non-symmetric
      type tmp_col_box_class

        ! Initial and final level, forbidden flag
        integer:: ifrom, ito, flag

        ! Rates
        double precision, dimension(:), allocatable:: C

        ! Pointer to next box
        type(tmp_col_box_class), pointer:: next

      end type tmp_col_box_class

!#####################################################################

      !> Structure within elastic_class to store individual entries
      type elastic_entry_class

        ! K, type, nz
        integer:: K, typo, nz

        ! Coefficients
        double precision:: a, b, c

        ! Explicit
        double precision, dimension(:), allocatable:: Coeff

      end type elastic_entry_class

!#####################################################################

      !> Structure to temporally contain elastic collisional rates
      type elastic_class

        ! Level, number of entries
        integer:: ilevel, nentry

        ! Substructure with data
        type(elastic_entry_class), dimension(:), allocatable:: datum

      end type elastic_class

!#####################################################################

      !> Structure to temporally contain inelastic collisional rates
      type inelastic_class

        ! Type of collision, up, low, forbidden, index in Tbox
        integer:: col_type, up, low, forbid, ind

        ! Substructure with data
        double precision, dimension(:), allocatable:: Cul

      end type inelastic_class

!#####################################################################

      !> Boxes to arrocate the temperature tables
      type Tbox_class

        ! Type of interpolation
        logical:: flin

        ! Number of indexes, index, type of collision
        integer:: nTmp, ind, col_type, nion

        ! Temperatures
        double precision, dimension(:), allocatable:: temp

        ! Pointer to next
        type(Tbox_class), pointer:: next

      end type Tbox_class

!#####################################################################

      !> Boxes to link arrays for omegabuildin/I
      type iabox_class

        ! Integers
        integer:: ifreq, iph, ith, mfreq, nback

        ! Array
        integer, dimension(:), allocatable:: A

        ! Pointer to next
        type(iabox_class), pointer:: next,prev

      end type iabox_class

!#####################################################################

      !> Boxes to link arrays for omegabuildin/I
      type dbabox_class

        ! Integers
        integer:: ifreq, mfreq, nback, iph, ith

        ! Array
        double precision, dimension(:), allocatable:: A

        ! Pointer to next
        type(dbabox_class), pointer:: next,prev

      end type dbabox_class

!#####################################################################

      !> Structure to have array of strings of big size
      type strarr_class

        ! String with 500 characters
        character(len=500):: str

      end type strarr_class

!#####################################################################

      !> Structure to have array of strings of size 2
      type strarr2_class

        ! String with 2 characters
        character(len=2):: s

      end type strarr2_class

!#####################################################################

      !> Structure to old data about the association between
      !! frequencies and transitions
      type Freqflag_class

        ! OR() of the vector flag Vabsent
        logical:: absent

        ! OR() of the vector flag Vabsent for Master
        logical, dimension(:), allocatable:: Mabsent

        ! Vector with flag of presence of a line in that frequency
        logical, dimension(:), allocatable:: Vabsent

      end type Freqflag_class

!#####################################################################

      !> Structure to hold voigt profiles
      type Voigt_class

        ! First order profile for intensity
        double precision, dimension(:), allocatable:: p

        ! First order profile for polarization
        complex(kind=8), dimension(:), allocatable:: cp

      end type Voigt_class

!#####################################################################

      !> Structure to hold normalization values
      type Nindex_class

        ! Tells if there is stored profile
        logical:: VRAM

        ! Normalization factor for absorption profiles
        double precision, dimension(:,:,:,:), allocatable:: Norm

        ! First order profile
        type(Voigt_class), dimension(:,:,:,:), allocatable:: prof

      end type Nindex_class

!#####################################################################

      !> Structure to store indexes of the M within a term for the
      !! atom
      type Atomindex_class

        ! Indexing of the mu and m of a term
        integer, dimension(:,:), allocatable:: ind

      end type Atomindex_class

!#####################################################################

      !> Structure to store indexes M within a term for the Voigt
      !! profiles
      type AtomVindex_class

        ! Number of components
        integer:: ncom, ncomNB

        ! Indexing of the Ju and Jf of a transition
        integer, dimension(:,:), allocatable:: indNB

        ! Indexing of the mu and mf of a transition
        integer, dimension(:,:,:,:), allocatable:: ind

      end type AtomVindex_class

!#####################################################################

      !> Structure to store indexes M within a term for the
      !! redistribution
      type Warrindexb_class

        ! Indexing of the mu and m of a term, same without magnetic
        ! field (J ordering)
        integer, dimension(:,:,:,:,:), allocatable:: Wind,WindNB

      end type Warrindexb_class

!#####################################################################

      !> Structure for the indexing of components for the
      !! redistribution function
      type Warrindexa_class

        ! Number of input transitions,
        integer:: nt

        ! Indexes of transitions
        integer, dimension(:), allocatable:: ind

        ! Indexing of the mu and m of a term
        type(Warrindexb_class), dimension(:), allocatable:: trani

      end type Warrindexa_class

!#####################################################################

      !> Structure with all the photoionization data
      type Phot_class

        ! Non presence of the photoionization in the CPU
        logical:: absent

        ! Non presence of the photoionization in the CPU for Master
        logical, dimension(:), allocatable:: Mabsent

        ! Upper and lower levels
        integer:: ilevell,ilevelu

        ! Type of input, number of frequencies for the cross section,
        ! initial frequency index, final frequency index
        integer:: mode, nfreq, if0, if1

        ! initial frequency index from master, final frequency from
        ! master
        integer, dimension(:), allocatable:: Mif0, Mif1

        ! Photoionization edge, (2J_l+1)/(2J_u+1), frequency weights
        ! in the boundaries of the subdomain
        double precision:: edge, glu, W0, W1

        ! Proper variation of cross section with frequency, integral
        ! for T_E SEE rate
        double precision, dimension(:), allocatable:: alpha, TEI
        ! Frequencies in the input atomic model and corresponding
        ! cross sections
        double precision, dimension(:), allocatable:: infreq,inalpha
        ! Frequency weights in the boundaries of the subdomain for
        ! master
        double precision, dimension(:), allocatable:: MW0, MW1

      end type Phot_class

!#####################################################################

      !> Structure with the fine structure transition data (because
      !! in general the atom is multiterm)
      type FST_class

        ! Upper and lower terms
        integer:: iterml,itermu

        ! Number of FS transitions
        integer:: nt

        ! Index of upper and lower J levels
        integer, dimension(:), allocatable:: ilevell, ilevelu

        ! Index of FS transition (equivalent to irad in Atom_class)
        integer, dimension(:,:), allocatable:: irad

        ! Aul and Blu for each FS transition within a transition
        ! between terms
        double precision, dimension(:,:), allocatable:: Aul,Blu

      end type FST_class

!#####################################################################

      !> Structure with the output frequency data for when using
      !! redistribution files
      type Ftrano_class

        ! Number of ranges, real limits in frequency axis, total
        ! number of frequencies, maximum number of input frequencies,
        ! limits in frequency for interpolation
        integer:: nran, gf0, gf1, nfreq, mxfreq, ggf0, ggf1

        ! Limits for each range
        integer, dimension(:), allocatable:: if0, if1

      end type Ftrano_class

!#####################################################################

      !> Structure with the dipole strength for a given transition
      type rdip_class

        ! Dipole strength matrix
        double precision, dimension(:,:,:,:,:), allocatable:: rdip

      end type rdip_class

!#####################################################################

      !> Structure with the level-level indexing
      type Jrho_class

        ! K,Q indexing
        integer, dimension(:,:), allocatable:: kq

      end type Jrho_class

!#####################################################################

      !> Structure with the dipole strength for a given transition
      type irho_class

        ! Level-level indexing
        type(Jrho_class), dimension(:,:), allocatable:: Jrho

        ! Indexing of rho by term
        integer, dimension(:), allocatable:: irho_ij

      end type irho_class

#ifdef _OPENMP
!#####################################################################

      !> Structure with how to split the components in a given
      !! transition for first order RT coefficients
      type omp_1c_class

        ! Limits for each range
        integer, dimension(:), allocatable:: if0, if1

      end type omp_1c_class

!#####################################################################

      !> Structure with how to split the components in a given
      !! transition for emiss2ord
      type omp_2c_class

        ! Limits for each range without and with NCHLT
        integer, dimension(:,:), allocatable:: if0, if1, nif0, nif1

        ! Minimum and maximums for iU and iU1 to skip CRD profile
        ! calculations for CHLT
        integer, dimension(:,:), allocatable:: mnU,mxU,mnU1,mxU1

        ! Minimum and maximums for iU and iU1 to skip CRD profile
        ! calculations for NCHLT
        integer, dimension(:,:), allocatable:: mnnU,mxnU,mnnU1,mxnU1

      end type omp_2c_class
#endif

!#####################################################################

      !> Structure with the atomic data
      type Atom_class

        ! Information of the presence of a transition at some
        ! frequency
        type(Freqflag_class), dimension(:), allocatable:: fflag

        ! FS transitions information
        type(FST_class), dimension(:), allocatable:: fst

        ! Photoionization information
        type(Phot_class), dimension(:), allocatable:: phot

        ! Indexing for Warr2
        type(Warrindexa_class), dimension(:), allocatable:: trano

        ! Internal indexing for Voigt
        type(AtomVindex_class), dimension(:), allocatable:: i_Vind

        ! Internal indexing for Warr2
        type(Atomindex_class), dimension(:), allocatable:: i_Wind

        ! Normalization of the first order profiles
        type(Nindex_class), dimension(:,:,:), pointer:: Normp

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

        ! Atom indexing
        type(irho_class), dimension(:), allocatable:: irho

        ! PRD transition ranges
        type(Ftrano_class), dimension(:), allocatable:: Ftrano
#ifdef _OPENMP
        ! How to split between components for first order
        type(omp_1c_class), dimension(:), allocatable:: omp_1c

        ! How to split between components for emiss2ord
        type(omp_2c_class), dimension(:), allocatable:: omp_2c
#endif
        ! Name of the atomic species
        character(len=2):: Element

        ! Label
        character(len=10):: file_label

        ! Name of the Voigt file
        character(len=20):: vfile, wfile
#ifdef _OPENMP
        ! Component splitting in absorb and emiss
        logical, dimension(:), allocatable:: omp_comp_1ord
#endif
        ! Normalize the total relative factor for abundance or not,
        ! multilevel flag, keep populations fixed, make zero the
        ! last ion, fix populations for lower term
        logical:: anorm, ML, fixp, zero_ion, fixplt

        ! Input spectrum for b-b and b-f transitions
        logical, dimension(:), allocatable:: bbspecin, bfspecin

        ! Flag to identify automatic Hydrogen background atom
        logical:: cust=.False.

        ! Flag to protect agains molecular equilibrium
        logical:: mol_protect=.False.

        ! Flag for the second order emissivity for each transition,
        ! flag to split between components in the multiterm atom when
        ! building frequency axis
        logical, dimension(:), allocatable:: lemiss2, splitf

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
        ! redistribution, maximum number of components for profile
        ! in file
        integer:: nMulti, ntran, nJmax, nKmax, nMmax, NNN, ndim, &
                  nfreq, nphot, ngk, ncol, nlevel, nftran, tshift, &
                  pshift, tfshift, ntrano, Mncom

        ! Sizes of profile files headers
        integer:: hvifil, hwifil, hwtfil

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
        ! line-wise K cut
        integer, dimension(:), allocatable:: nfreqt, nfreqtc, nJ, &
                                             stage, broad_type, &
                                             term, sublevel, &
                                             nfreqph, col_type, &
                                             if0, if1, ifst, itrano, &
                                             rif0, rif1, rif20, &
                                             rif21, sbif0, sbif1, &
                                             sfif0, sfif1, ilf0, &
                                             ilf1, ipf0, ipf1, Kcut, &
                                             Krad

        ! Number of the atomic dimension for each running index in
        ! first order
        integer, dimension(:), allocatable:: nL,nU,nMu,nMl

        ! Transition indexing matrix, collisional transition indexing
        ! matrix, photoionization indexing matrix initial frequency
        ! index from master, final frequency from master
        integer, dimension(:,:), allocatable:: irad, icol, iphot, &
                                               Mif0, Mif1

        ! Indexing of the FS transitions
        integer, dimension(:,:), allocatable:: ifst_ij

        ! Block dimension (index for M,term) and flag of forbidden
        ! collisions between levels
        integer, dimension(:,:), allocatable:: nblk, fcflag

        ! J indexing for each individual magnetic level
        integer, dimension(:,:,:), allocatable:: iJval

        ! Atomic part of the thermal Doppler width, atomic mass,
        ! abundance, multiplicative factor for abundance
        double precision:: cDopp, rmass, abun, abun_mod = 1d0

        ! Sizes of voigt file chunks
        double precision, dimension(:), allocatable:: zsize, dsize, &
                                                      tsize, tBsize, &
                                                      f0size, f1size

        ! orbital angular momentum, spin angular momentum, term
        ! frequency, Doppler widths to include, Doppler widths for
        ! the core, Frequency of term transition, number density of
        ! the atom in cm-3, Stark broadening, degeneration (terms
        ! and H custom), lower and upper frequency limit weights
        double precision, dimension(:), allocatable:: rLval, &
                                                      Sval, &
                                                      TRfreq, &
                                                      Dwvl, Dwvlc, &
                                                      Dfreq, n, &
                                                      broad_stark, &
                                                      deg, W0, W1, gL

        ! Inverse lifetime, level frequency and J values, radiative
        ! transition rates matrix, arguments for the Van der Waals
        ! Broadening, collisional transition-wise damping,
        ! populations, lte populations, departure coefficient,
        ! lower and upper frequency limit weights for master
        double precision, dimension(:,:), allocatable:: damp, &
                          FSfreq, rJval, Ecoeff, Norm, NormS, &
                          broad_args, ldamp, popu, populte, depar, &
                          MW0, MW1

        ! Inelastic collisional terms rates for allowed transitions
        ! between terms and transition rates between levels
        double precision, dimension(:,:,:), allocatable:: Ccoeff, &
                                                          CcoeffJ

        ! Size of jump for redistribution
        double precision, dimension(:,:,:), allocatable:: dztwsize

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

      !> Structure to hold complex profiles
      type LTEprof_sub_class

        ! First order profile for polarization
        complex(kind=8), dimension(:), allocatable:: cp

      end type LTEprof_sub_class

!#####################################################################

      !> Structure to hold normalization values
      type LTEprof_class

        ! Tells if there is stored profile
        logical:: VRAM

        ! First order profile for intensity
        double precision, dimension(:), allocatable:: p

        ! First order profile for polarization
        type(LTEprof_sub_class), dimension(:,:), allocatable:: comp

      end type LTEprof_class

!#####################################################################

      !> Structure with LTE lines data
      type LTEline_class

        ! Profiles
        type(LTEprof_class), dimension(:,:), pointer:: prof

        ! Name of the atomic species
        character(len=2):: Element

        ! If there is a background model atom, if limited in tau, if
        ! limited in T, line is absent for a CPU
        logical:: is_passive, taulim_l, Tlim_l, absent

        ! Element, ion, number of Ml, number of Mu, maximum height
        ! index for presence, total number of frequencies, number
        ! of frequencies for the line core, initial frequency, final
        ! frequency, index of the passive model atom, type of Van
        ! der Waals broadening, type of line
        integer:: ele, stage, nMl, nMu, Rz0, nfreq, nfreqc, &
                  if0, if1, ia, broad_type, ltype

        ! Upper level energy, lower level energy, upper level angular
        ! momentum, lower level angular momentum, upper level Landé
        ! factor, lower level Landé factor, Einstein coefficient for
        ! e.e., resonance frequency, Einstein coefficient for
        ! absorption, optical depth limit, temperature limit, range in
        ! Doppler widths for the line, range in Doppler widths for the
        ! core, Stark broadening, collisional oscillator strength,
        ! radiative broadening, abundance
        double precision:: eu, el, Ju, Jl, gu, gl, Aul, Dfreq, &
                           Blu, taulim, Tlim, Dwvl, Dwvlc, &
                           broad_stark, f_c, broad_rad, abund

        ! Inverse lifetime, inverse lifetime upper level, arguments
        ! for the Van der Waals broadening, collisional
        ! transition-wise damping
        double precision, dimension(:), allocatable:: damp, &
                                                      broad_args

        ! Lower level population (fraction), upper level population
        ! (fraction), total population
        double precision, dimension(:), allocatable:: nl, nu, n

        ! Atomic part of the thermal Doppler width, atomic mass,
        double precision:: cDopp, rmass

      end type LTEline_class

!#####################################################################

      !> Structure with the atomic data for a Kurucz line
      type ckurucz_class

        ! Element, charge
        integer:: A, Z

        ! Einstein coefficient, resonance, radiative, Van der Waals,
        ! and Stark broadening parameters, statistical weights upper
        ! and lower level, lower and upper limit in frequency,
        ! abundance, lower level energy
        double precision:: Aul, Dfreq, Grad, Gvdw, Gstk, gu, gl, &
                           O0, O1, abun, Ei

      end type ckurucz_class

!#####################################################################

      !> Structure with Kurucz data
      type kurucz_class

        ! Number of lines
        integer:: ntran

        ! Line data
        type(ckurucz_class), dimension(:), allocatable:: tran

        ! Maximum frequency distance
        double precision:: MDomg

      end type kurucz_class

!#####################################################################

      !> Structure with the atomic data for a molecule
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
        integer, dimension(:), allocatable:: imol
        integer, dimension(:), allocatable:: nmol

        ! Abundance
        double precision:: abun = 0

        ! Ionization eng
        double precision, dimension(:), allocatable:: Eion

        ! PF
        double precision, dimension(:,:), allocatable:: pf

      end type catm_class

!#####################################################################

      !> Structure with the full molecular data
      type Mol_class

        ! Name of atoms in molecule
        type(strarr2_class), dimension(:), allocatable:: catom

        ! Name of the molecule
        character(len=:), allocatable :: Molecule

        ! Charge, number of species, total number of atoms, type
        ! of fit, number of pf coefficients, number of eqc
        ! coefficients
        integer:: Charge, nA, nAT, pffit, npfcoeff, neqcoeff

        ! Number of each atom in the molecule and position in the
        ! atom list for chemical equilibrium
        integer, dimension(:), allocatable:: natom, iatom

        ! Atomic part of the thermal Doppler width, mass, dissociation
        ! energy, minimum temperature, maximum temperature
        double precision:: cDopp, rmass, Den, Tmin, Tmax

        ! Coefficients for partition function, partition function,
        ! Coefficients for equililibrum constant, equilibrium constant
        ! Population
        double precision, dimension(:), allocatable:: pfcoeff, pf, &
                                                      eqcoeff, eq, n

      end type Mol_class

!#####################################################################

      !> Input spectrum for CLE
      type spect_class

        ! Axial
        logical:: axial, valid, pol

        ! Size of spectrum
        integer:: nfreq,nmu,nphi,nstk

        ! Wavenumber, cosine of polar angle, phi angle
        double precision, dimension(:), allocatable:: omega, &
                                                      mu, &
                                                      phi
        ! Stokes spectrum
        double precision, dimension(:,:,:,:), allocatable:: stokes

        ! Stokes spectrum interpolated to the directional quadrature
        ! in a point
        double precision, dimension(:,:,:,:), allocatable:: mustokes

      end type spect_class

!#####################################################################

      !> Structure with a pointer to 1D arrays
      type chianti_real_pointer_class

        ! pointer
        real, dimension(:), pointer:: p

      end type chianti_real_pointer_class

!#####################################################################

      !> CHIANTI ioneq ion fraction data
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
      !! energy
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

      !> Atomic-atmospheric data (partition function and ioniuzation
      !! energy
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

      !> Structure with all the thermodynamical quantities of the
      !! atmosphere
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
        integer:: nZ, typo, NT, nele, mode, norm

        ! log of gravity acceleration at surface, frequency
        ! of tau scale, y position for CLE, z position for CLE,
        ! diffuse light fraction
        double precision:: logg, tfreq, ypos, zpos, f_diff

        ! Height, Temperature, microturbulence, velocity components,
        ! magnetic field components, auxiliar velocity components to
        ! cheat in intensity
        double precision, dimension(:), pointer:: z, T, vmi, &
                                                  vx, vy, vz, &
                                                  Bx, By, Bz, &
                                                  vxa, vya, vza

        ! Total hydrogen, hydrogen minus, gas pressure, density,
        ! electron pressure, alternative height scale, and atomic
        ! hydrogen number density, partition function temperature,
        ! abundance, input asymmetry radiation field tensors, to
        ! be used by the inversion, electron number density, los
        ! velocity, pos velocity, pos azimuth velocity
        double precision, dimension(:), allocatable:: nHT, nHm, &
                                                      Pg, rho, Pe, &
                                                      zalt, nHA, pT, &
                                                      abund, JKQin, &
                                                      ne, vlos, &
                                                      vpos, vphi

        ! Continuum absorption at reference frequency
        double precision, dimension(:), allocatable:: chi500

        ! Array with zeros
        double precision, dimension(:), pointer:: zeros

        ! Hydrogen density, helium density
        double precision, dimension(:,:), allocatable:: nh, nhe

      end type Atmo_class

!#####################################################################

      !> Structure with the magnetic field data
      type Bfield_class

        ! Module, polar angle and azimuth of B field
        double precision, dimension(:), allocatable:: Bstrength, &
                          Btheta, Bphi

        ! Longitudinal magnetic field, transversal magnetic field, and
        ! magnetic field azimuth in the POS (these are used in
        ! inversion)
        double precision, dimension(:), allocatable:: Blos, Bpos, &
                                                      Azimuth

      end type Bfield_class

!#####################################################################

      !> Structure with the geometry data
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

      !> Structure with the geometry data
      type Geometry_class

        ! Flag for axial symmetry
        logical:: axial

        ! Number of polar nodes, number of real azimuthal nodes,
        ! number of azimuthal nodes for emiss2ord, number of emergent
        ! polar directions, number of emergent azimuthal directions
        ! Number of nodes for AA integral
        integer:: nTh, nPh, nPh2, nThLOS, nPhLOS, nThAA

        ! Indexing of 2D directions
        integer, dimension(:,:), allocatable:: i_geom

        ! Gamma angle
        double precision:: gam

        ! Vector of cosines of polar angle nodes, vector of cosiones
        ! of azimuthal angle nodes, vector of sign of sinus of
        ! azimuthal angle nodes, vector of polar angles, vector of
        ! azimuthal angles, weights of polar integral, weights of
        ! RT azimuthal integral, weights of emiss2 azimuthal
        ! integral, cosines of polar angle of emergent directions,
        ! azimuthal angles of emergent directions, Vector of angles
        ! for the AA integral, weights for the AA integral, nodes for
        ! classical gaussian quadrature, weights for classical
        ! gaussian quadrature, polar angle on the disk for a given
        ! quadrature in a point above the surface
        double precision, dimension(:), allocatable:: V_mu, V_mux, &
                          V_muy, V_theta, V_phi, W_mu, W_mux, &
                          W_mux2, L_mu, L_theta, L_phi, V_thetaAA, &
                          W_muAA, V_gauss, W_gauss, V_mu_disk

        ! TKQ geometrical tensor in the vectical reference frames
        ! for the emergent problem
        complex(kind=8), dimension(:,:,:), allocatable:: TSL

        ! TKQ geometrical tensor in the magnetic reference frame
        ! for the emergent problem
        complex(kind=8), dimension(:,:,:,:), allocatable:: TBL

        ! TKQ geometrical tensor in the vectical reference frames
        ! for the formal problem
        complex(kind=8), dimension(:,:,:,:,:), allocatable:: TS

        ! TKQ geometrical tensor in the magnetic reference frame
        ! for the formal problems
        complex(kind=8), dimension(:,:,:,:,:,:), allocatable:: TB

      end type Geometry_class

!#####################################################################

      !> Structure with the settings
      type IO_helper_class

        ! Buffer size
        integer:: buffer_size=1

        ! Header size, number of 'ranges', integer for header
        integer:: head_size, nran, nn, geom_size

        ! Secondary indexes, size per range
        integer, dimension(:), allocatable:: sindx,nbuff

        ! Indexes to control output
        integer, dimension(:,:), allocatable:: indx

        ! Doubles to control output
        double precision, dimension(:,:), allocatable:: doub

      end type IO_helper_class

!#####################################################################

      !> Structure with the settings
      type FWHM_helper_class

        ! If single gaussian value, if pending initialization
        logical:: gaussian, toinit

        ! Number of ranges (duplicated), number of wavelengths (if
        ! not gaussian)
        integer:: nn, nfreq

        ! Secondary indexes, size per range
        integer, dimension(:), allocatable:: sindx

        ! Indexes to control output, indexes for quick
        ! interpolation
        integer, dimension(:), allocatable:: indx

        ! Indexes for quick interpolation
        integer, dimension(:,:), allocatable:: indx1,indx2

        ! Doubles to control output, wavelength, kernel
        double precision, dimension(:), allocatable:: doub, wave, &
                                                      kernel

        ! Auxiliar for quick interpolation
        double precision, dimension(:,:), allocatable:: idx

      end type FWHM_helper_class

!#####################################################################

      !> Structure to have array of strings of big size
      type strnum_class

        ! String with 500 characters
        character(len=500):: str

        ! Type
        integer:: typ

        ! Value
        double precision:: val

      end type strnum_class

!#####################################################################

      !> Structure with the node values and locations.
      type Node_class

        ! Number of special boundary conditions
        integer:: nebound = 0

        ! Indexes in atmosphere for values of tau in nodes
        integer, dimension(:), allocatable:: Tau_Indx

        ! Normal boundary limits
        double precision, dimension(2):: Bounds

        ! Positions, values, errors
        double precision, dimension(:), allocatable:: H, Var, Errors

        ! Positions and limits for special boundary conditions
        double precision, dimension(:,:), allocatable:: ebound

      end type Node_class


!#####################################################################

      !> Structure with the settings
      type Input_class

        ! Structures to help with 1.5D outputs
        type(IO_helper_class):: lim_stk,lim_ctr,lim_tau, &
                                lim_cols_tt,lim_cols_ll, &
                                lim_damp,lim_back,lim_pop, &
                                lim_atmo

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
        ! save parameters of damping, alternative broadcasting, add
        ! bound-bound background transitions, memoization of J
        ! symbols, protect hydrogen from the equation of state, apply
        ! NG acceleration to intensity, if RAM use should be reported,
        ! protect all atoms in chemical equilibrium, measure
        ! performance in blocks, measure performance per CPU, keep a
        ! file with the frequency dependent JKQ, if using coherent
        ! wings approximation, if keeping solution files, keep
        ! populations, keep departure coefficients, keep output
        ! rhoKQ, keep output JKQ, keep stokes in quadrature,
        ! keep MRC, skip first iteration when reading performance
        ! data from previous run, skip the disk in CLE, if
        ! input spectrum loaded, restrict in tau_c, restrict in z,
        ! angle-averaged forced in intensity problem, force intensity
        ! problem to be static, if there is a file for weights,
        ! keep collisions log, keep MPI log, keep, MPI detailed log
        logical:: AV, appendMRC, appendMRCI, out_contr, out_tau1, &
                  store, storeI, Pcorr, Raman, keepIsol, &
                  NG, keep_back, keep_damp, keep_cols, bfieldn, &
                  keep_aparam, altbcast, addbb, keep_atmo, memo, &
                  protect_H, NGI, RAMreport, chem_protect_all, &
                  asym, g_perf, mpi_perf, keep_jkqnu, cohw, &
                  keep_sol, keep_pop, keep_dep, keep_rhoKQ, &
                  keep_JKQ, keep_stokesQ, keep_MRC, IWskip, &
                  skip_disk, lspect_input, rest_tau, rest_z, AVI, &
                  static_int, linv_weight, keep_coll, keep_mpil, &
                  keep_mpidl

        ! If asymmetry input
        logical, dimension(2):: lasym

        ! Type of Doppler width input to build frequency axis
        character(len=3) dws

        ! Mode of solution switch, force type of problem (I or Stokes)
        ! heights, vertical scale for model atmosphere
        character(len=1) mode, force, atm_scale

        ! Run ID
        character(len=9) ID

        ! output folder, atmospheric file, continuum file,
        ! magnetic field files, solution file, file with fudge
        ! factors, input file name, cache file, partition function
        ! file name, abundance file name, barklem file with sp
        ! data, barklem file with pd data, barklem file with df data,
        ! file with previous MPI details, file with MPI performance,
        ! path to CHIANTI database, file with CLE input spectra,
        ! file with weights for the inversion
        character(len=500) folder, atmo, continuum, bfield, &
                           solution, fudge, source, resource, input, &
                           cache, pf, abund, bark_sp, bark_pd, &
                           bark_df, MPIdetail, operform, &
                           chianti_path, spect_input, inv_weight

        ! Name of atomic files, population files, background atom
        ! files, background atom population files, molecules,
        ! Kurucz line files, wavelength files, asymmetry files
        type(strarr_class), dimension(:), allocatable:: atom, popu, &
                                              atomback, popuback, &
                                              mol, kurucz, waves, &
                                              asym_fil, fwhm_fil

        ! Ionization fraction data
        type(strnum_class), dimension(:), allocatable:: ionf

        ! Force no magnetic field, force no velocity, to force
        ! observed frequencies in synthesis axis
        logical:: unmagnetized, static, force_inv_freq

        ! Keep atomic populations fixed, zero out the last ion
        ! populations, skip this atom for the wavelength axis,
        ! fix populations for lower term
        logical, dimension(:), allocatable:: fixp, zero_ion, &
                                             skip_wave, fixplt

        ! Number of first iteration, maximum possible iterations,
        ! order of iteration (emissivity), number of steps between
        ! saving data, number of atoms, number of background atoms,
        ! number of molecules, number of first iteration in intensity,
        ! maximum possible iterations in intensity, order of iteration
        ! (emissivity) in intensity, number of steps between saving
        ! data in intensity, order of the NG acceleration, delay in
        ! iterations before applying NG acceleration, first iteration
        ! to apply ALI to, number of Kurucz line files, maximum
        ! iteration to allow non physical quantities in stokes and
        ! rho, number of wavelength files, mode of Zeeman effect,
        ! update atmospheric model at the end of the calculation for
        ! intensity, recompute electron density, order of the NG
        ! acceleration for intensity, delay in iterations before
        ! applying NG acceleration for intensity, number of asymmetry
        ! input given by constant numbers, number of asymmetry inputs
        ! given by files, number of asymmetry inputs, if need to take
        ! care of magnetically induced transitions, mode of running
        ! (synthesis 1D, syn. 1.5D, or inversion), number of CPU
        ! groups to split columns in 1.5D or inversion, type of
        ! presurre/density scale, polar nodes, azimuthal nodes,
        ! intensity polar nodes, intensity azimuthal nodes, AA
        ! integral nodes, intensity AA integral nodes, LOS polar
        ! directions, LOS azimuthal directions, number of LTE lines
        integer:: iter_min, iter_max, iter_ord, store_step, &
                  nA, nAb, nM, iteri_min, iteri_max, &
                  storei_step, iteri_prd, iter_j, NG_ord, NG_delay, &
                  ALI_delay, NK, allownphys_stk, allownphys_rho, NW, &
                  zeeman_mode, update_atmos, redo_ne, NGI_ord, &
                  NGI_delay, PRD_delay, nasym_num, nasym_fil, nasym, &
                  MIT_input, run_mode, rt_group_n, atmo_char, nTh, &
                  nPh, nThI, nPhI, nThAA, nThAAI, nThLOS, nPhLOS, nLTE

        ! Box to solve in 1.5D synthesis problem
        integer, dimension(:), allocatable:: sol_box

        ! Input for additional K cuts
        integer, dimension(:,:), allocatable:: Kcut_input

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
        ! consider
        double precision:: dw, MIT_node, dcohw, minT, maxT, maxV, &
                           omega_ref,T_rad,R_star,r0tc,r1tc,r0z,r1z

        ! LOS polar mus, LOS azimuthal angles
        double precision, dimension(:), allocatable:: L_mu,L_phi

        ! Parameters for redistribution input frequency axis: rang,
        ! reso, negl, vlar, fstp, mstp, core, rang_core, fstp_core,
        ! mstp_core, for polarization and for intensity, respectively
        double precision, dimension(11):: red_pars, redi_pars

        ! MRC for populations, MRC for rhoKQ with K!=0 and for rho00
        ! in the intensity problem
        double precision:: mrc_i, mrc_p, mrci_i, mrci_r

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
        ! restore the inversion from, name of the output file
        character(len=500):: Filename_Ob, Inv_init, &
                             Output_file

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
        ! (assymetries) must be in the output
        logical:: Broyden, FITSFILE, Sigma_neglect, auto_weight, &
                  centered, Pos_Correction, hydroeq, Fractional, &
                  Projection, Popuinit, Keep_RF, out_jkqa

        ! Flag to modify variable in the inversion, flag for
        ! the regularization of each variable
        logical, dimension(:), allocatable:: Nodes_Flags, Nodes_Regul

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
        ! of input atmosphere
        integer:: Num_Iter, Type_Inversion, Err_Type, &
                  Atmo_Input, LM_Method, Init_Thermal, &
                  Interpolation, btype, vtype, SVD_type, fits_index, &
                  Num_Weight, atmoin_type

        ! Output file sizes
        integer:: s_inv_h, s_inv_atmo, s_inv_atmo_c, &
                  s_inv_res_h, s_inv_res_c, s_inv_res, &
                  s_inv_RF_h, s_inv_RF_c


        ! Units to direct the verbosity
        integer, dimension(3):: Unit_VB

        ! Type of node value, number of nodes, index of the
        ! regularization for each variable the regulatization
        ! for each variable
        integer, dimension(:), allocatable:: Node_Type, Num_nodes, &
                                             Indx_regul

        ! Threshold in chi2, threshold for the fractional chi2,
        ! ratio limit for the regulatizations with respect to the
        ! proper chi2, threshold for the SVD, boundary value for
        ! the pressure, maximum step allowed in SVD, initial Bpos or
        ! B theta if input too small, initial B azimuth if input too
        ! small, factor to decrease lambda when LM iteration
        ! accepted, factor to increase lambda when LM iteration
        ! rejected, diffuse light factor, initial vz or vpos if
        ! input too small, initial vy or azimuth if input too small
        double precision:: Threshold_chisq, Chisq_fraction,  &
                           Regul_Limit, Threshold_svd, Pg_bound, &
                           Max_Step, ini_Bpos, ini_Bazi, &
                           factoraccept, factorreject, f_diff, &
                           ini_vpos, ini_vazi

        ! scale for each parameter, perturbation for each parameter,
        ! minimum relative perturbation
        double precision, dimension(:), allocatable:: Scal, Perturb, &
                                                      min_rel_Pert

        ! Tau ranges to consider, LM lambda ranges to consider
        double precision, dimension(2):: Tau_Range, Lam_Range

        ! Weight of the regularization function for each variable
        double precision, dimension(:), allocatable:: Regul_weight

        ! Weights for each Stokes
        double precision, dimension(:,:), allocatable:: Weight

        ! Data to specify atmospheric stratification modifications
        double precision, dimension(:,:), allocatable:: Atmo_strat

        ! Made-up stratification from inputs
        double precision, dimension(:), allocatable:: Atmo_strat_done

      end type Input_class

!#####################################################################

      !> Structure for the background opacities fudge factors
      type fudge_class

        ! Number of frequencies with data
        integer:: nfreq_f

        ! Fudge factor data
        double precision, dimension(:,:), allocatable:: fudge_v

      end type fudge_class

!#####################################################################

      !> Structure for the background opacities
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

      !> Structure with the maximum relative change information
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

      !> Structure with the MPI data
      type MPI_class

        ! Are we doing MPI, are we splitting frequencies,
        ! use alternative solverI, use alternative solver, use
        ! alternative solverJ, use alternative solver JKQgen,
        ! alternative broadcasting, if sending columns to slaves
        logical:: mpi, lf, alternI, alternP, alternJ, &
                  alternJgen, altbcast, mpi15d

        ! Number of processors, identifier, number of
        ! slave processors, global number of processors, global
        ! identifier, number of groups to send columns
        integer:: nproc, pid, nnd, gnproc, gpid, ngroup

        ! Processor to receive in custom bcast, number of bcast steps
        integer:: nsend, recv, steps

        ! Integers for non-blocking transfer
#ifdef oldmpi
        integer:: &
#else
        type(MPI_request):: &
#endif
                  request1, request2, request3, request4, request5, &
                  request6, request7, request8, request9, request0, &
                  request11

        ! Maximum number of frequencies per processor,
        ! maximum number of profile size per processor, same for
        ! intensity, same for photoionizations
        integer:: nxfreq, nxtfreq, nxtfreqi, nxpfreq

        ! Number of frequencies for the processor and its limits,
        ! sizes of packages to transfer, list of processors to send in
        ! custom bcast, leaders of the slave groups, number of
        ! frequencies for the processor and its limits (for the input
        ! in CLE)
        integer, dimension(:), allocatable:: nf, if0, if1, &
                                             size1, size2, size3, &
                                             size4, size5, size6, &
                                             size7, size8, &
                                             sizei1, sizei2, sizei3, &
                                             sizei4, sizei5, sizei6, &
                                             sizei7, sizei8, sizei9, &
                                             sizei0, size10, &
                                             sizei10, sizei11, &
                                             sizei12, sizei13, &
                                             sizei14, sizei4b, &
                                             lsend, &
                                             ltslave, &
                                             inf,iif0,iif1

        ! Domain decomposition isend requests
#ifdef oldmpi
        integer, dimension(:,:), allocatable:: requestA
#else
        type(MPI_request), dimension(:,:), allocatable:: requestA
#endif

        ! RAM used to store Voigt and Warr, RAM used for ionizations
        ! RAM used to store radiation field quantities, RAM used
        ! to store background quantities
        double precision:: RAM, PRAM, VRAM, WRAM, RRAM, BRAM

      end type MPI_class

!#####################################################################

      ! Types for memoization

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

      !> Structure with factorial and signs
      type Fctsg_class

        ! Doing memoization with jagged arrays
        logical:: memo
#ifdef _OPENMP
        ! Allow writing in memoization
        logical:: can_write = .True.
#endif

        ! Memoization of J symbols
        type(a6D) :: J6
        type(a6D) :: J3
        type(a9D) :: J9

        ! Factorial and sign
        double precision, dimension(:), allocatable:: flg, sg

      end type Fctsg_class

!#####################################################################

      !> Structure with output frequencies for redistribution
      type Frequencyd_class

        ! If interpolation can be stored
        logical:: RAM

        ! Size of omega and W_freq, size of interpolation data
        integer:: osize, isize

        ! Input frequency size
        integer, dimension(:), allocatable:: mfreq

        ! Frequencies and weights
        double precision, dimension(:), allocatable:: omega,W_freq

        ! Indexes for interpolation
        integer, dimension(:), pointer:: index1, index2

        ! Interpolation double
        double precision, dimension(:), pointer:: dx

      end type Frequencyd_class

!#####################################################################

      !> Structure with input transitions for redistribution
      type Frequencyc2_class

        ! Max of input frequencies, number of output frequencies,
        ! number of output frequency ranges, maximum frequency limits
        ! for 2ord, maximum frequency limits for interpolation in
        ! 2ord (for AA hybrid)
        integer:: mxfreq, nfreq, nran, gf0, gf1, ggf0, ggf1

        ! Indexes of transitions, frequency limits for 2ord
        integer, dimension(:), allocatable:: if0,if1

        ! Input by pair of transitions
        type(Frequencyd_class), dimension(:), pointer:: trani

#ifdef _OPENMP
        integer, dimension(:), allocatable:: oif0,oif1
#endif
      end type Frequencyc2_class

!#####################################################################

      !> Structure with the frequency tree data for the redistribution
      type Frequency_class

        ! Minimum and maximum index with photoionizations and lines,
        ! Number of polar and azimuthal directions, total directions,
        ! total directions that are in the quadrature, size of
        ! frequency space for profile messages, same for intensity,
        ! same for photoionizations, dimension of direction/height/
        ! atom/output transition array
        integer:: pif0, pif1, lif0, lif1, nth, nph, ndir, nqdir, &
                  ntfreq, ntfreqi, npfreq, ndzao

        ! Minimum and maximum index with photoionizations for the
        ! master, and for the lines, weight for each frequency node
        ! for sharing tasks, size of frequency space for profile
        ! messages, same for intensity, number of forward scattering
        ! directions, weight for each frequency node for sharing task
        ! but neglecting PRD, mapping of output frequencies (CLE) into
        ! general axis
        integer, dimension(:), allocatable:: Mpif0, Mpif1, &
                                             Mlif0, Mlif1, &
                                             IW_freq, Mntfreq, &
                                             Mntfreqi, Mnpfreq, &
                                             nfs,IW_freq_in, &
                                             mapping

        ! Type of scattering: -1: forward, 0: normal, 1: backward
        integer, dimension(:,:,:), allocatable:: stype

        ! Indexing
        integer, dimension(:,:,:,:), allocatable:: indx

        ! Frequency, weight (output), frequency to the cube
        double precision, dimension(:), allocatable:: omega, &
                                                      W_freq, &
                                                      omega3, &
                                                      omega_ou, &
                                                      omega3_ou

        ! Exponentials for epsIphoto
        double precision, dimension(:,:), pointer:: exu

        ! Input by height
        type(Frequencyc2_class), dimension(:), pointer:: dzao

      end type Frequency_class

!#####################################################################

      !> Structure with a flag and output frequencies for
      !! redistribution
      type Redd_class

        ! If Wfunc2 has to be calculated, if Wfunc2 can be stored
        logical:: iIPRD = .True., RAM

        ! If Wfunc2 has to be calculated
        logical, dimension(:), allocatable:: iPPRD

        ! Redistribution function for intensity
        real, dimension(:), allocatable:: IWarr2

        ! Redistribution function for polarization
        complex(kind=4), dimension(:), allocatable:: PWarr2

      end type Redd_class

!#####################################################################

      !> Structure with input transitions for redistribution
      type Redc2_class

        type(Redd_class), dimension(:), pointer:: trani

      end type Redc2_class

!#####################################################################

      !> Structure to store the redistribution function tree data
      type Red_class

        ! Number of polar and azimuthal directions, total directions,
        ! total directions in non-redistribution problem, dimensions
        ! of direction/height/atom/output transition array
        integer:: nth, nph, ndir, njdir, ndzao

        ! Indexing
        integer, dimension(:,:,:,:), allocatable:: indx

        ! Input by height, directions and atoms
        type(Redc2_class), dimension(:), pointer:: dzao

      end type Red_class

!#####################################################################

      !> Structure to store rhoKQ data temporarily
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

      !> Structure with the node information.
      type Nodes_class

        ! Nodes for the inversion variables
        type(Node_class), dimension(nvar_inv):: Node

        ! If regularizing, if the pressure is given at boundary,
        ! if correcting the node positions from the atmosphere,
        ! if taking the gas pressure from the input model, if
        ! hydrostatic equilibrium must be imposed each time
        logical:: Regul_Flag, hydroeq, Pg_Inv, &
                  Pos_Correction, Pg_Auto, hydros

        ! Flag to modify variable in the inversion, flag for
        ! the regularization of each variable
        logical, dimension(nvar_inv):: Nodes_Flags, Nodes_Regul

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
        integer:: Node_Location_Type, Btype, vtype, Indx_b, &
                  Indx_e, Tot_Nodes, Num_Fit, Num_Mag, Num_Thermal, &
                  Num_Asymmetry, Interpolation, Nodes_Type, Num_glob

        ! Type of node value for each variable, number of nodes for
        ! each variable, number of nodes that can change for each
        ! variable, index of the regularization for each variable,
        ! final number the regulatization "points"
        integer, dimension(nvar_inv):: Node_Type, Num_Nodes, &
                                       Num_Vary, Indx_regul, Num_regul

        ! Indexes of the first and last nodes that can change for each
        ! variable
        integer, dimension(2,nvar_inv):: Node_Vary

        ! Index of nodes and parameters
        integer, dimension(:,:), allocatable:: Inf_Inv

        ! Threshold for the SVD, maximum step allowed in SVD, boundary
        ! value for the pressure, cosine of the heliocentric angle for
        ! the emergence, azimuth for emergence
        double precision:: Threshold_svd, Max_Step, Pg_Bound, mu, &
                           azimuth

        ! Weight of the regularization function for each variable,
        ! value of the perturbations to the model, scale value for
        ! each parameter, value for constant regularization, minumum
        ! relative perturbation
        double precision, dimension(nvar_inv):: Regul_weight, &
                                                Perturb, Scal, &
                                                Const, min_rel_Pert

      end type Nodes_class

!#####################################################################

      !> Structure with the solution of the RT problem
      type Solution_class

        !
        !type(Sol_class):: Sol_Tmp, Sol_Min, Sol_Sav

        ! If interpolating frequencies, if return fractional
        ! polarization, if project the magnetic field, if
        ! diffuse light
        logical:: Fre_Intp, Fractional, Projection, Diff_flag

        ! Number of wavelength points, number of wavelength ranges
        integer:: Num_Wavelength, Num_Range

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

      !> Structure with the full solution of the RT problem
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

      !> Structure with the information of input Stokes profiles.
      type Stokes_class

        ! If the sigma is the same for all profile points, if
        ! automatic weights for Stokes, if sigmas read at the
        ! beginning, if diffuse light, if diffuse light read at
        ! the beginning
        logical:: Sigma_flag, Auto_Weight, Sigma_ct, Diff_flag, &
                  Diff_ct

        ! Number of wavelengths, type of sigma in the input, number
        ! of wavelength ranges, total number of data points,
        ! intensity data points
        integer:: Num_Wavelength, Indx_Sigma, Num_Range, &
                  Num_freedom, Num_freedomI

        ! Wavelength ranges
        integer, dimension(:,:), allocatable:: Range

        ! Cosine of the heliocentric angle for the observation, and
        ! azimuth of the LOS
        double precision:: mu, azimuth

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

      !> Structure with information for the regularization
      type Regul_class

        ! Scale of the penalty, regulatization penalty
        double precision:: Ratio, Penalty

        ! Regularization vector
        double precision, dimension(:), allocatable:: Regul_F

        ! Regularization matrix
        double precision, dimension(:,:), allocatable:: Regul_H

      end type Regul_class

!#####################################################################

      ! Struture for the LM step
      type LMFIT_class

        ! Regularizations
        type(Regul_class):: Rgl

        !
        logical:: Flag_weight, Flag_Jac, accepted

        !
        integer:: Num

        !
        double precision:: factoraccept, factorreject, Lambda, &
                           Chisq, Chisq_og, Chisq_0

        !
        double precision, dimension(2):: Lambda_bounds

        !
        double precision, dimension(:), allocatable:: ResidualI, &
                                        WeightI, Jacfvec, &
                                        Jacfvec_og, Diag

        !
        double precision, dimension(:,:), allocatable:: Residual, &
                                          Weight, Hessian, &
                                          Hessian_og, JacobianI

        !
        double precision, dimension(:,:,:), allocatable:: Jacobian

      end type LMFIT_class

!#####################################################################

      end module types_mod
