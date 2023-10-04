      !> Routines to generate frequency axis
      module omegabuild_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Contributors:
!     John Dennis (NCAR)
!  Start:
!     04/18/2017
!  Last version:
!     08/24/2023 V3.0.10
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/24/2023:   V3.0.10 - Added the possibility to force
!                             wavelengths in the inversion data in the
!                             synthesis axis. This is the first case
!                             of added frequencies with high priority
!                             being able to remove others within the
!                             resolution threshold (TdPA)
!
!     08/07/2023:    V3.0.8 - Added the contribution of LTE lines to
!                             the frequency axis (TdPA)
!
!     07/03/2023:    V3.0.7 - Implemented the possibility of a model
!                             atom not contributing to the frequency
!                             axis with nodes (TdPA)
!                           - Limited some verbosity to the great
!                             master (TdPA)
!                           - Removed a control call (TdPA)
!
!     02/14/2023:    V3.0.6 - Distinguish between the PRD problem
!                             being fully AA, or only the intensity
!                             solution (TdPA)
!                           - Distinguish between the problem
!                             being fully axial, or only the intensity
!                             solution (TdPA)
!                           - The parameters to build the frequency
!                             quadrature for the PRD integral is
!                             different for intensity and polarization
!                             in general (TdPA)
!
!     11/24/2022:    V3.0.5 - Create a second array of weights if
!                             running CLE (TdPA)
!                           - Added refitfrec (TdPA)
!                           - Moved the creation of an auxiliar array
!                             of frequency weights in frecresize
!                             inside its relevant block (TdPA)
!
!     10/26/2022:    V3.0.4 - Master deallocates Atom%fflag%Vabsent
!                             which it does not use (TdPA)
!                           - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!                           - Remove Vabsent because it was a stupid
!                             waste of RAM (TdPA)
!
!     10/25/2022:    V3.0.3 - Implemented the restriction of the
!                             height axis (TdPA)
!                           - Initialize pointers in omegabuildin
!                             and omegabuildinI (TdPA)
!                           - Changed call to cleanFrecandRed which
!                             requires a new argument (TdPA)
!                           - Removed unused resetWarr routine (TdPA)
!                           - Bugfix: Deallocate mfreq array in
!                             Frec%dzao%trani, memory leak (TdPA)
!                           - Bugfix: Added deallocation of omp
!                             indexes (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: The allocation of the Master
!                             quantities cannot be done in omegabuild
!                             because for the non-1D cases the Master
!                             is still unknown (TdPA)
!                           - Added the omegainitmaster routine (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Atmo and Bstrength are no longer
!                                inputs for omegabuild. Instead, maxB
!                                is a new input to indicate if there
!                                will be a magnetic field.
!                              o Atmo%v has changed to Atmo%vx,%vy,
!                                and %vz.
!                             (TdPA)
!
!     06/21/2022:    V2.1.0 - Modifications to account for coherent
!                             scattering in the observers frame for
!                             a selected Doppler width far from the
!                             line core (TdPA)
!                             NOTE: Limited testing (AA, static, and
!                             non-magnetic)
!
!     03/24/2021:    V2.0.2 - Bugfix: The indexes where indocrectly
!                             split among the threads in the routines
!                             omegabuildin/I (TdPA)
!                           - Bugfix: In omegabuildin, when using
!                             OpenMP one shared variable was being
!                             changed within a single construct before
!                             every thread was done with it. Added a
!                             barrier before the construct (TdPA)
!                           - Now allocating only the actually needed
!                             space for input transitions in Frec and
!                             Red for PRD (TdPA)
!                           - Added an extra OpenMP barrier in
!                             omegabuildin/I just in case (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Added OpenMP when building the input
!                             frequency axis (TdPA)
!
!     02/09/2021:   V1.16.1 - In omegabuild, making different the
!                             initial and final limits of b-b and b-f
!                             transitions when not present in a
!                             given process (TdPA)
!
!     02/04/2021:   V1.16.0 - Now take into account MIT when building
!                             the frequency axis, but depending on
!                             the user input (TdPA)
!
!     01/13/2021:   V1.15.6 - Added a check for magnetic field
!                             presence in check_nchlt (TdPA)
!                           - Added a missing call to verbose in
!                             check_nchlt (TdPA)
!
!     11/17/2020:   V1.15.5 - Makes sure that the amount of elements
!                             to allocate to store the PRD functions
!                             and the interpolation quantities is
!                             positive (TdPA)
!
!     11/12/2020:   V1.15.4 - Bugfix: Because index1, index2, and dx
!                             are now pointers, they must be
!                             nullified when allocating %trani in
!                             omegabuildin(I) (TdPA)
!
!     11/04/2020:   V1.15.3 - When locating the extremes in the line
!                             and photoionization indexes, the
!                             variables were being reset inside the
!                             atomic loop (TdPA)
!
!     10/26/2020:   V1.15.2 - Clean the Frec and Red structures using
!                             the specific subroutine (TdPA)
!                           - Changes because index1 and index2 are
!                             now pointers (TdPA)
!
!     09/15/2020:   V1.15.1 - Bugfix: Forgot to limit the transition
!                             loop of the Frec%indx array to the lines
!                             that are present in any given CPU (TdPA)
!                           - The mfreq array was allocated with the
!                             wrong size (TdPA)
!
!     09/11/2020:   V1.15.0 - Changed completely the Frec and Red
!                             structures (TdPA):
!                             .Output directions, heights,
!                              atoms, and output transitions have been
!                              combined in a single array dzao and an
!                              index has been created to map such
!                              dimensions.
!                             .Output frequencies, input directions,
!                              and output frequencies have been
!                              combined and the omega, W_freq,
!                              index1, index2, and dx have the total
!                              size. There is no indexing, we trust
!                              in the calling order. This applies to
!                              the redistribution functions as well,
!                              including the atomic levels dimension
!                              in the polarization case.
!                           - Storing the interpolation quantities is
!                             now optional and constraint by the RAM
!                             limit (TdPA)
!                           - Input frequencies and weights are taken
!                             into account to compute the RAM (TdPA)
!                           - Added an extra variable to communicate
!                             if the RAM limit was reached inside the
!                             routine (TdPA)
!                           - Added routine to clean the frequency
!                             and redistribution structures (TdPA)
!
!     07/31/2020:   V1.14.1 - Made sure that dimensions for future
!                             allocations (ntfreqi, ntfreq, and
!                             npfreq in Frec%) are not zero (TdPA)
!
!     07/10/2020:   V1.14.0 - Changes in omegabuildin to take into
!                             account the new NHCLT approximation,
!                             which is applied height by height
!                             depending on the magnetic field (TdPA)
!
!     07/01/2020:   V1.13.2 - Bugfix: Some CPU numbers could lead to
!                             an allocation problem when the range
!                             of output frequencies in the
!                             redistribution is limited through
!                             the options in the input file (TdPA)
!
!     06/26/2020:   V1.13.1 - Changed what check_nchlt does (TdPA)
!
!     06/05/2020:   V1.13.0 - In omegabuildin and omegabuildinI, now
!                             the interpolation index array has been
!                             split into two, one for the previous and
!                             one for the next point (JD)
!
!     06/02/2020:   V1.12.6 - In omegabuildin and omegabuildinI, now
!                             check if there were formal solutions
!                             when going in through emergence, in
!                             order to avoid double passings (TdPA)
!
!     06/01/2020:   V1.12.5 - Bugfix: There are instances in which
!                             Atom(ia)%i_Wind could be allocated
!                             twice. Solved the problem (TdPA)
!                           - Reduce the dimensionality of the
!                             redistribution matrix when the
!                             non-coherent lower term approximation is
!                             applied (TdPA)
!                           - Added check_nchlt, routine to print a
!                             message when the magnetic field is too
!                             small (but larger than zero) for the
!                             non-coherent lower term (TdPA)
!
!     11/19/2019:   V1.12.4 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - Changes of how the memory is traced
!                             among photoionizations, voigt, and
!                             redistribution (TdPA)
!                           - Now omegabuildin/I are called with the
!                             LOS parameter (TdPA)
!
!     11/12/2019:   V1.12.3 - Added definition of Atom%rif0 and
!                             Atom%rif1, so every CPU know where the
!                             lines really start and end (TdPA)
!                           - Bugfix: When calculating the space that
!                             the redistribution takes, it was
!                             considered double instead of double
!                             complex (TdPA)
!
!     10/31/2019:   V1.12.2 - Bugfix: The initializers for the Frec
!                             variables %ntfreq, %ntfreqi, and %npfreq
!                             were inside the atomic loop, and must
!                             be before (TdPA)
!
!     10/03/2019:   V1.12.1 - Now output transitions get indexed, so
!                             you do not have to be careful with
!                             transition ordering (TdPA)
!
!     09/13/2019:   V1.12.0 - Added the possibility to import extra
!                             wavelength nodes from files (TdPA)
!                           - Added determination of input frequencies
!                             for angle-averaged PRD with velocities.
!                             Had to add the computation of the total
!                             frequency range that emiss2ord uses
!                             for interpolation (TdPA)
!
!     07/23/2019:   V1.11.5 - The distance between frequencies is now
!                             checked in wavelength units (TdPA)
!                           - Bugfix: The fix in 1.11.4 was wrongly
!                             coded (TdPA)
!
!     07/19/2019:   V1.11.4 - Bugfix: The split in ranges for far
!                             frequencies could affect
!                             photoionizations of high energy levels.
!                             Frequencies in the range of a transition
!                             with only one FS component or photo-
!                             ionizations are protected from
!                             splitting. Added warning for splitting
!                             and for protected (TdPA)
!                           - The distance between frequencies to
!                             consider different ranges is now a
!                             parameter, jump (TdPA)
!
!     06/06/2019:   V1.11.3 - Now the memory that will be filled in
!                             Red% is predicted for each block before
!                             allocating (TdPA)
!
!     06/03/2019:   V1.11.2 - Not there is an option to not split
!                             between FS components when building the
!                             frequency axis (TdPA)
!
!     05/31/2019:   V1.11.1 - Handles the new variables in the
!                             Frec (frequency) structure to handle
!                             the profile and ratio variables.
!                           - Bugfix: Red%njdir was not initialized
!                             for angle-averaged or static runs (TdPA)
!                           - Left only a placeholder for
!                             omegainoptimize (TdPA)
!
!     05/08/2019:   V1.11.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!                           - Bugfix: There must be a index for
!                             directions to be used by Frec%stype that
!                             is limited to Red%njdir, before I got
!                             out of bounds for some cases (TdPA)
!                           - Introduced routines to optimize the
!                             input frequency ranges. But not using
!                             them (TdPA)
!
!     03/18/2019:   V1.10.1 - If a transition has only one frequency
!                             in its range, it is ignored (TdPA)
!
!     02/20/2019:   V1.10.0 - New verbosity (TdPA)
!                           - Using specific TINY parameters (TdPA)
!
!     02/14/2019:    V1.9.4 - Bugfix: Wrong formatting of warning
!                             message for big Doppler shifts (TdPA)
!
!     02/11/2019:    V1.9.3 - Bugfix: With velocities, the cosine
!                             and sine values were being reused and
!                             the resonances were completely out of
!                             place, resulting in zero integrals in
!                             the redistribution (TdPA)
!                           - Bugfix: Now the dx in the frequency
!                             structure contains the (x-x0) term too,
!                             it is convenient for the dynamic cases
!                             specially (TdPA)
!
!     02/08/2019:    V1.9.2 - Introduced warnings when the velocity
!                             shifts are too big for the ranges
!                             specified for the lines (TdPA)
!
!     09/21/2018:    V1.9.1 - Allocating stype if PRD even if AA or
!                             not dynamic, to avoid complains when
!                             doing RAM debugging (TdPA)
!
!     08/06/2018:    V1.9.0 - Allows a limit on the amount of RAM to
!                             allocate in order to store Wfunc2 in
!                             omegabuildin and omegabuildinI (TdPA)
!
!     08/03/2018:    V1.8.2 - Make omegabuildin coherent with the
!                             changes in emiss2ord (TdPA)
!
!     07/27/2018:    V1.8.1 - Skipping in the memory allocations the
!                             same combinations that will not be
!                             computed in emiss2ord (TdPA)
!
!     05/16/2018:    V1.8.0 - Introduced the definition of an index
!                             to flag forward and backward scattering
!                             when AD redistribution (TdPA)
!
!     12/05/2017:    V1.7.0 - Possibility to ignore l!=f in second
!                             order emissivity. (TdPA)
!
!     10/30/2017:    V1.6.1 - By popular demand, stored 1/4 into a
!                             parameter (what is used to determine
!                             if two quantum numbers are the
!                             same (TdPA)
!
!     10/03/2017:    V1.6.0 - Taken into account non-magnetic case for
!                             the allocation of Warr (TdPA)
!
!     09/22/2017:    V1.5.0 - Possibility to limit K (TdPA)
!
!     09/08/2017:    V1.4.0 - Introduced pointers in buildin to
!                             reduce verbosity and to improve the
!                             performance when compiling without
!                             optimization flags, important for the
!                             debugging mode (TdPA)
!
!     08/29/2017:    V1.3.0 - Introduced two ranges for construction
!                             of input frequency axis (TdPA)
!                           - Changed the weighting of frequencies
!                             with PRD depending on the region, wing
!                             or core (TdPA)
!
!     08/14/2017:    V1.2.4 - Bugfix: One absence check was missing
!                             in omegabuildinI (TdPA)
!
!     08/10/2017:    V1.2.3 - Changed nM variable to nMm, to avoid
!                             possible conflicts with the common
!                             molecule count variable (TdPA)
!                           - Added checks for absent lines in both
!                             omegabuildin and omegabuildinI to avoid
!                             entering in the loops when it is obvious
!                             that is not needed (TdPA)
!                           - Bugfix: Added a special treatment for
!                             single frequency processors, was not
!                             working as intended (TdPA)
!                           - Added an early scape when the frequency
!                             finder determines that the outputs of
!                             the process are out of range (TdPA)
!
!     08/09/2017:    V1.2.2 - Bugfix: There was a jdir that should be
!                             jbdir (TdPA)
!
!     08/01/2017:    V1.2.1 - Only weight a frequency as PRD if there
!                             is PRD computation (TdPA)
!
!     07/21/2017:    V1.2.0 - omega and Wfreq are one step higher in
!                             Frec (TdPA)
!
!     07/20/2017:    V1.1.3 - Changed the initial value of bf1 in
!                             omegabuild, so if no transition is
!                             found, it does not run the second loop
!                             to count frequencies (TdPA)
!                           - In omegabuild, always deallocate Red,
!                             even if not renewing it with the
!                             polarized one (TdPA)
!                           - If not storing Warr2, no need to
!                             allocate iPPRD (TdPA)
!
!     07/19/2017:    V1.1.2 - Bugfix: In omegabuildinI, the output
!                             frequencies were being compared with
!                             the input transition resonance (TdPA)
!                           - Bugfix: If line absent, no need for
!                             2nd order information (TdPA)
!                           - Bugfix: Typo in array element, a jtran
!                             that should be a jdir (TdPA)
!                           - Bugfix: mxfreq was reset for each
!                             input transition (therefore, was not
!                             the maximum) (TdPA)
!
!     06/29/2017:    V1.1.1 - Removed resetWarrI and simplified
!                             resetWarr (TdPA)
!                           - Bugfixing from compiler complains (TdPA)
!
!     06/28/2017:    V1.1.0 - Enormous amount of changes to take
!                             into account velocities and AD
!                             redistribution in omegabuildin and
!                             omegabuildinI (TdPA)
!
!     06/23/2017:    V1.0.9 - Omegabuild also build an array with
!                             weights for nodes (TdPA)
!
!     06/22/2017:    V1.0.8 - There were some remaining itran that
!                             had to be iti in omegabuildin (TdPA)
!
!     06/20/2017:    V1.0.7 - Added resetWarrI and resetWarr, that
!                             makes the Warr redistribution function
!                             not initialized (TdPA)
!
!     06/19/2017:    V1.0.6 - Changed the structure of the Warr2 that
!                             is allocated to store Warr2 for
!                             the intensity part (TdPA)
!                           - Added code to allocate space for Warr2
!                             in polarization mode (TdPA)
!
!     06/16/2017:    V1.0.5 - Changed RAM to IRAM (TdPA)
!                           - Changed how the input axis are build
!                             to accomodate for the new structure
!                             of the class tree (TdPA)
!                           - Accomodated the change of if0/1 to
!                             pif0/1 and added the correspondent code
!                             for lif0/1 (TdPA)
!
!     06/14/2017:    V1.0.4 - Bugfix: Frec%if0/1 where reset inside
!                             the photoionization loop, they have to
!                             be outside such loop (TdPA)
!
!     06/12/2017:    V1.0.3 - Limits for b-b transitions are stored,
!                             and weights for these limits (TdPA)
!                           - The limits are adjusted to each CPU in
!                             frecresize (TdPA)
!                           - Allocated space for Warr2 (TdPA)
!                           - omegabulidin/I also uses the b-b
!                             boundary information (TdPA)
!                           - There is no flag to decide to
!                             interpolate in emiss2
!
!     05/12/2017:    V1.0.2 - Implemented multilevel version of the
!                             input frequency axis build. The
!                             multi-term checks if this structure
!                             exists and deallocates it before (TdPA)
!
!     05/05/2017:    V1.0.1 - Extended limits must be double size
!                             than the standard (TdPA)
!                           - omegabuildin should not be called if
!                             the line is not flagged as PRD (TdPA)
!
!     04/18/2017:    V1.0.0 - First version (TdPA)
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
!  omegabuild:
!    Calculates the frequency axis and the weights
!
!  omegainitmaster:
!    Initializes frequency related variables needed by the RT master
!
!  omegabuildin:
!    Calculates the input frequency axis and the weights
!
!  freqresize:
!    Resize weights and absence vector for each processor domain
!
!  refitfrec::
!    Resize arrays and setup index limits for the task splitting
!  for the CLE RT
!
!  check_nchlt
!    Checks where the non-coherent lower term approximation can be
!    applied
!
!  cleanFrecandRed
!    Safely clean the Frec and Red structures
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use omp_mod
      use parameters_mod , only : c , TINYA, TINYB, TINYO, TINYWAR , &
                                  resol, resolin , pi , IPI42, jump
      use profile_mod
      use qsort_mod
      use types_mod

      integer, parameter:: ContW = 1
      integer, parameter:: PhotW = 2
      integer, parameter:: LLTEW = 5
      integer, parameter:: LCRDW = 10
      integer, parameter:: LCOHW = 20
      integer, parameter:: LPRDW = 500


      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Build the output frequency axis.\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!       Atom(Atom_class): Structure with the atomic data\n
      !!     Input(Input_class): Structure with settings data\n
      !!           maxB(dfloat): Maximum magnetic field strength\n
      !!            lp(logical): Bool that says if it is the
      !!                         polarization problem or not\n
      !!    obs_wave(dfloat(:)): In inversion mode, the data
      !!                         wavelengths, a dummy array
      !!                         otherwise
      subroutine omegabuild(Frec,Atom,Input,maxB,lp,obs_wave)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Input_class):: Input
      type(Frequency_class):: Frec
      logical, intent(in):: lp
      double precision, intent(in):: maxB
      double precision, dimension(:), intent(in):: obs_wave

      ! Local

      logical:: init, core, cohw
      logical:: Yfield
      logical:: warmed, MIT
      logical, dimension(:), allocatable:: protect
      logical, dimension(:), allocatable:: vMIT

      integer:: ia,itran,jtran,ktran,iterml,itermu,itermf
      integer:: ifreq,jfreq,cfreq,it,lnfreq,lnfreqc
      integer:: i,i1,iJl,iJu,nt,ios
      integer:: iL,iL1,iU,iU1,mF,iMf,iMu,iMu1,iMl,iMl1,nMf,nMl,nMu
      integer:: ip,iq,ip1,iq1
      integer, dimension(:), allocatable:: flag

      double precision, parameter:: dLSJ = .25d0

      double precision:: norm,norm1,v,dv,O0,O1,prdf,prdc,cmpf,mincf
      double precision:: S,rLf,rLu,rLl,rJf,rJu,rJu1,rJl,rJl1
      double precision:: rJfmax,rJumax,rJlmax,maxV,disp
      double precision:: rMf,rMu,rMu1,rMl,rMl1,q,q1,p,p1
      double precision:: prdcohw
      double precision, dimension(:), allocatable:: omegaaux
      double precision, dimension(:), allocatable:: omega
      double precision, dimension(:), allocatable:: nut
      double precision, dimension(:), allocatable:: DwT
      double precision, dimension(:), allocatable:: DwTL
      double precision, dimension(:), allocatable:: DwTmin
      double precision, dimension(:), allocatable:: compfact
      double precision, dimension(:,:), allocatable:: tmp_lim
      double precision, dimension(:,:), allocatable:: tmp_limL

      ! Routine name
      urou = 'omegabuild'


      !
      ! Read the number of wavelengths from the wavelength files
      ! and add it to the total frequencies
      !
      do i=1,Input%NW

        ! Open
        open (200,file=trim(Input%waves(i)%str), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', form='unformatted')

        ! Read size
        read(200,err=1100) jfreq

        ! Add size to the total
        nfreq = nfreq + jfreq

        ! Close
        close(200)

      end do ! For each wavelength file

      ! If forcing inversion data frequencies
      if (Input%force_inv_Freq) then

        ! Add size to the total
        nfreq = nfreq + size(obs_wave)

      end if

      ! Take the Doppler width from the input
      allocate(DwT(nA))

      ! Maximum Doppler width
      if(Input%dws.eq.'MAX')then

        do ia=1,nA
          DwT(ia) = Atom(ia)%cDopp*sqrt(Input%maxT)
        end do

      ! Minimum Doppler width
      else if(Input%dws.eq.'MIN')then

        do ia=1,nA
          DwT(ia) = Atom(ia)%cDopp*sqrt(Input%minT)
        end do

      ! Or fixed input Doppler width
      else if(Input%dws.eq.'NUM')then

        do ia=1,nA
          DwT(ia) = Input%dw*1d-9/c
        end do

      end if

      ! If LTE lines
      if (nLTEl.gt.0) then

        ! Take the Doppler width from the input
        allocate(DwTL(nLTEl))

        ! Maximum Doppler width
        if(Input%dws.eq.'MAX')then

          do ia=1,nLTEl
            DwTL(ia) = Input%LTEline(ia)%cDopp*sqrt(Input%maxT)
          end do

        ! Minimum Doppler width
        else if(Input%dws.eq.'MIN')then

          do ia=1,nLTEl
            DwTL(ia) = Input%LTEline(ia)%cDopp*sqrt(Input%minT)
          end do

        ! Or fixed input Doppler width
        else if(Input%dws.eq.'NUM')then

          do ia=1,nLTEl
            DwTL(ia) = Input%dw*1d-9/c
          end do

        end if
      end if ! LTE lines

      ! Compute minimum Doppler width if coherent wing
      if (Input%cohw) then

        allocate(DwTmin(nA))

        do ia=1,nA
          DwTmin(ia) = Atom(ia)%cDopp*sqrt(Input%minT)
        end do

      end if

      ! Correct nfreq if including MIT transitions
      if (Input%MIT_input.ge.0) then

        ! For each atom
        do ia=1,nA

          ! Check if need to include MIT transitions
          if (Atom(ia)%ML) then
            MIT = .False.
          else
            if (Input%MIT_input.gt.0) then
              MIT = .True.
            else if (Input%MIT_input.lt.0) then
              MIT = .False.
            else
              if (maxB.gt.TINYB) then
                MIT = .True.
              else
                MIT = .False.
              end if
            end if
          end if

          ! If not MIT, skip
          if (.not.MIT) cycle

          ! Skipping wavelengths?
          if (Input%skip_wave(ia)) cycle

          ! bound-bound transitions
          do itran=1,Atom(ia)%ntran

            ! If not splitting between components, skip
            if (.not.Atom(ia)%splitf(itran)) cycle

            ! Identify terms
            do i=1,Atom(ia)%nMulti-1
              do i1=i+1,Atom(ia)%nMulti
                if (Atom(ia)%irad(i,i1).eq.itran) then
                  iterml = i
                  itermu = i1
                end if
              end do
            end do

            ! Loop fine transitions
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If transition exist, skip
                if (Atom(ia)%fst(itran)%irad(iJu,iJl).ge.1) cycle

                ! If MIT
                if (MIT) then

                  lnfreq = nint(dble(Atom(ia)%NfreqT(itran))* &
                                Input%MIT_node)
                  if (mod(lnfreq,2).eq.0) lnfreq = lnfreq + 1
                  lnfreqc = nint(dble(Atom(ia)%NfreqTc(itran))* &
                                 Input%MIT_node)
                  if (mod(lnfreqc,2).eq.0) lnfreqc = lnfreqc + 1
                  if (lnfreqc.lt.0) lnfreqc = 3
                  if (lnfreq.lt.lnfreqc) lnfreq = lnfreqc + 2

                  ! Add frequency
                  nfreq = nfreq + lnfreq
                  Atom(ia)%nfreq = Atom(ia)%nfreq + lnfreq

                end if

              end do ! Lower J level
            end do ! Upper J level
          end do ! L-L transition
        end do ! Atom

      end if ! If possible MIT transitions

      ! Allocate temporal limit array
      cfreq = nfreq*2
      allocate(tmp_lim(2,nxtran))
      if (nLTEl.gt.0) allocate(tmp_limL(2,nLTEl))

      ! Set counter of frequencies to 0
      cfreq = 0

      ! Get maximum number of FS transitions
      nMu = 0
      do ia=1,nA
        nMl = maxval(Atom(ia)%nJ)
        if (nMl.gt.nMu) nMu = nMl
      end do
      nMu = nMu*nMu

      ! Define an initial vector
      allocate(omega(nfreq))
      ! Fine structure frequencies
      allocate(nut(nMu))
      ! Vector for MIT flag
      allocate(vMIT(nMU))

      ! Initialize quantities to check velocity shifts
      if (dyn) then
        maxV = Input%maxV
      end if

      ! Initialize flag
      warmed = .False.


      !
      ! Collect the frequencies
      !

      ! Take the frequencies of the transitions
      do ia=1,nA

        ! Make sure that there is memory to work with
        if(.not.allocated(omegaaux))then
          allocate(omegaaux(Atom(ia)%nfreq))
        else
          if(size(omegaaux).lt.Atom(ia)%nfreq)then
            deallocate(omegaaux)
            allocate(omegaaux(Atom(ia)%nfreq))
          end if
        end if

        ! Check if need to include MIT transitions
        if (Atom(ia)%ML) then
          MIT = .False.
        else
          if (Input%MIT_input.gt.0) then
            MIT = .True.
          else if (Input%MIT_input.lt.0) then
            MIT = .False.
          else
            if (maxB.gt.TINYB) then
              MIT = .True.
            else
              MIT = .False.
            end if
          end if
        end if

        ! bound-bound transitions
        do itran=1,Atom(ia)%ntran

          ! Apply atom shift
          jtran = itran + Atom(ia)%tshift

          ! Initialize the limits of the line
          tmp_lim(1,jtran) =  1D99
          tmp_lim(2,jtran) = -1D99

          ! Identify terms
          do i=1,Atom(ia)%nMulti-1
            do i1=i+1,Atom(ia)%nMulti
              if (Atom(ia)%irad(i,i1).eq.itran) then
                iterml = i
                itermu = i1
              end if
            end do
          end do

          ! If splitting between components
          if (Atom(ia)%splitf(itran)) then

            ! Count fine transitions
            nt = 0
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If transition does not exist
                if (Atom(ia)%fst(itran)%irad(iJu,iJl).lt.1) then

                  ! If MIT
                  if (MIT) then

                    nt = nt + 1
                    vMIT(nt) = .True.

                  ! If not MIT
                  else

                    ! Just cycle
                    cycle

                  end if

                ! Permitted transition
                else

                  nt = nt + 1
                  vMIT(nt) = .False.

                end if

                nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJl,iterml)
              end do
            end do

          ! Not splitting between components
          else

            nt = 1
            nut(1) = Atom(ia)%Dfreq(itran)
            vMIT(1) = .False.

          end if

          ! Add frequencies for each FS transition
          do it=1,nt

            ! If MIT
            if (vMIT(it)) then
              lnfreq = nint(dble(Atom(ia)%NfreqT(itran))* &
                            Input%MIT_node)
              if (mod(lnfreq,2).eq.0) lnfreq = lnfreq + 1
              lnfreqc = nint(dble(Atom(ia)%NfreqTc(itran))* &
                             Input%MIT_node)
              if (mod(lnfreqc,2).eq.0) lnfreqc = lnfreqc + 1
              if (lnfreqc.lt.0) lnfreqc = 3
              if (lnfreq.lt.lnfreqc) lnfreq = lnfreqc + 2
            ! If not MIT
            else
              lnfreq = Atom(ia)%NfreqT(itran)
              lnfreqc = Atom(ia)%NfreqTc(itran)
            end if

            ! Linear part of the axis
            omegaaux(1) = 0d0
            do ifreq=2,lnfreqc/2+1
              omegaaux(ifreq) = omegaaux(ifreq-1) + &
                                Atom(ia)%Dwvlc(itran)/ &
                                dble(lnfreqc/2)
            end do

            ! Logarithmic part
            v = log10(Atom(ia)%Dwvlc(itran))
            dv = log10(Atom(ia)%Dwvl(itran)) - v
            do ifreq=lnfreqc/2+2,lnfreq/2+1
              v = v + dV/dble(lnfreq/2 - lnfreqc/2)
              omegaaux(ifreq) = 1d1**v
            end do

            ! Build symmetric axis
            omegaaux(lnfreq/2+2:lnfreq) = omegaaux(2:lnfreq/2+1)
            omegaaux(1:lnfreq/2) = -omegaaux(lnfreq/2+1:2:-1)
            omegaaux(lnfreq/2+1) = 0d0

            ! Check limits
            O0 = minval(omegaaux(1:lnfreq))* &
                 Atom(ia)%Dfreq(itran)*DwT(ia) + nut(it)
            O1 = maxval(omegaaux(1:lnfreq))* &
                 Atom(ia)%Dfreq(itran)*DwT(ia) + nut(it)
            if (O0.lt.tmp_lim(1,jtran)) tmp_lim(1,jtran) = O0
            if (O1.gt.tmp_lim(2,jtran)) tmp_lim(2,jtran) = O1

            ! Skip wavelengths?
            if (Input%skip_wave(ia)) cycle

            ! Build real contribution
            do ifreq=1,lnfreq

              ! Add to frequency axis
              omega(ifreq + cfreq) = nut(it) + &
                                     omegaaux(ifreq)* &
                                     Atom(ia)%Dfreq(itran)*DwT(ia)

            end do ! Atomic frequencies

            ! Accumulate the defined frequencies
            cfreq = cfreq + lnfreq

          end do ! FS transition

          ! Check that the range is adecuated to the velocity imposed
          ! in this frame
          if (dyn.and.gpid.eq.0) then
            do it=1,nt
              disp = maxV*nut(it)/DwT(ia)/Atom(ia)%Dfreq(itran)
              if (2d0*disp.gt.Atom(ia)%Dwvl(itran)) then
                if (.not.warmed) then
                  umsg = '###'
                  call verbose
                  umsg = '### IMPORTANT WARNING'
                  call verbose
                  umsg = '###'
                  call verbose
                  warmed = .True.
                end if
                write(umsg,'(A,i4,",",i4,3A)') &
                  ' - Warning: transition ',itran,it,' in ', &
                  Atom(ia)%Element,' atom can be shifted more '// &
                  'than half of the total width specified'
                call verbose
              else if (disp.gt.Atom(ia)%Dwvlc(itran)) then
                if (.not.warmed) then
                  umsg = '###'
                  call verbose
                  umsg = '### IMPORTANT WARNING'
                  call verbose
                  umsg = '###'
                  call verbose
                  warmed = .True.
                end if
                write(umsg,'(A,i4,",",i4,3A)') &
                  ' - Warning: transition ',itran,it,' in ', &
                  Atom(ia)%Element,' atom can be shifted more '// &
                  'than the core width specified'
                call verbose
              else if (2d0*disp.gt.Atom(ia)%Dwvlc(itran)) then
                if (.not.warmed) then
                  umsg = '###'
                  call verbose
                  umsg = '### IMPORTANT WARNING'
                  call verbose
                  umsg = '###'
                  call verbose
                  warmed = .True.
                end if
                write(umsg,'(A,i4,",",i4,3A)') &
                  ' - Warning: transition ',itran,it,' in ', &
                  Atom(ia)%Element,' atom can be shifted more '// &
                  'than half of the core width specified'
                call verbose
              end if
            end do
          end if

        end do ! b-b Transition

        ! bound-free transitions
        do itran=1,Atom(ia)%nphot

          ! Skip?
          if (Input%skip_wave(ia)) cycle

          ! If it is explicit, just add the input frequencies to the
          ! axis
          if (Atom(ia)%phot(itran)%mode.eq.0) then

            do ifreq=1,Atom(ia)%phot(itran)%nfreq

              omega(ifreq+cfreq) = Atom(ia)%phot(itran)%infreq(ifreq)

            end do

          ! If it is hydrogenic, introduce a linear axis for the
          ! interval between the edge and the maximum frequency
          else

            omega(1 + cfreq) = Atom(ia)%phot(itran)%edge

            nt = (1 - Atom(ia)%phot(itran)%mode)* &
                 Atom(ia)%phot(itran)%nfreq + &
                 Atom(ia)%phot(itran)%mode
            dv = (Atom(ia)%phot(itran)%infreq(nt) - omega(1+cfreq))/ &
                 dble(Atom(ia)%phot(itran)%nfreq - 1)

            do ifreq=2,Atom(ia)%phot(itran)%nfreq

              omega(ifreq + cfreq) = omega(ifreq - 1 + cfreq) + dv

            end do

          end if

          ! Accumulate the defined frequencies
          cfreq = cfreq + Atom(ia)%phot(itran)%nfreq

        end do ! b-f Transiton

      end do ! Atom

      !
      ! LTE lines
      !

      ! Take the frequencies of the transitions
      do ia=1,nLTEl

        ! If no frequencies, do not bother
        if (Input%LTEline(ia)%nfreq.le.0) cycle

        ! Make sure that there is memory to work with
        if(.not.allocated(omegaaux))then
          allocate(omegaaux(Input%LTEline(ia)%nfreq))
        else
          if(size(omegaaux).lt.Input%LTEline(ia)%nfreq)then
            deallocate(omegaaux)
            allocate(omegaaux(Input%LTEline(ia)%nfreq))
          end if
        end if

        ! Initialize the limits of the line
        tmp_limL(1,ia) =  1D99
        tmp_limL(2,ia) = -1D99

        ! Linear part of the axis
        omegaaux(1) = 0d0
        do ifreq=2,Input%LTEline(ia)%nfreqc/2+1
          omegaaux(ifreq) = omegaaux(ifreq-1) + &
                            Input%LTEline(ia)%Dwvlc/ &
                            dble(Input%LTEline(ia)%nfreqc/2)
        end do

        ! Logarithmic part
        v = log10(Input%LTEline(ia)%Dwvlc)
        dv = log10(Input%LTEline(ia)%Dwvl) - v
        do ifreq=Input%LTEline(ia)%nfreqc/2+2, &
                 Input%LTEline(ia)%nfreq/2+1
          v = v + dV/dble(Input%LTEline(ia)%nfreq/2 - &
                          Input%LTEline(ia)%nfreqc/2)
          omegaaux(ifreq) = 1d1**v
        end do

        ! Build symmetric axis
        omegaaux(Input%LTEline(ia)%nfreq/2+2: &
                 Input%LTEline(ia)%nfreq) = &
                               omegaaux(2:Input%LTEline(ia)%nfreq/2+1)
        omegaaux(1:Input%LTEline(ia)%nfreq/2) = &
                           -omegaaux(Input%LTEline(ia)%nfreq/2+1:2:-1)
        omegaaux(Input%LTEline(ia)%nfreq/2+1) = 0d0

        ! Check limits
        tmp_limL(1,ia) = minval(omegaaux(1:Input%LTEline(ia)%nfreq))*&
                         Input%LTEline(ia)%Dfreq*DwTL(ia) + &
                         Input%LTEline(ia)%Dfreq
        tmp_limL(2,ia) = maxval(omegaaux(1:Input%LTEline(ia)%nfreq))*&
                           Input%LTEline(ia)%Dfreq*DwTL(ia) + &
                           Input%LTEline(ia)%Dfreq

        ! Build real contribution
        do ifreq=1,Input%LTEline(ia)%nfreq

          ! Add to frequency axis
          omega(ifreq + cfreq) = Input%LTEline(ia)%Dfreq + &
                                 omegaaux(ifreq)* &
                                 Input%LTEline(ia)%Dfreq*DwTL(ia)

        end do ! Atomic frequencies

        ! Accumulate the defined frequencies
        cfreq = cfreq + Input%LTEline(ia)%nfreq

      end do ! LTE lines

      !
      ! Wavelength files
      !
      do i=1,Input%NW

        ! Open file
        open (200,file=trim(Input%waves(i)%str), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', form='unformatted')

        ! Read size of array
        read(200,err=1100) jfreq

        ! Read wavelengths
        read(200,err=1100) omega(cfreq+1:cfreq+jfreq)

        ! Update last index
        cfreq = cfreq + jfreq

        ! Close file
        close(200)

      end do ! Wavelength files

      !
      ! Inversion data wavelengths
      !
      if (Input%force_inv_Freq) then

        ! Add frequencies
        omega(cfreq+1:cfreq+size(obs_wave)) = 1d2/obs_wave

        ! For each data wavelengths
        do i=1,size(obs_wave)

          ! Check if already in
          do ifreq=1,cfreq

            ! If a close frequency exist, negate it
            if(abs(1d2/omega(ifreq)-obs_wave(i)).lt.resol) &
              omega(ifreq) = -1d0

          end do ! Synthesis freqs.
        end do ! Data freqs.

        ! Update size
        cfreq = cfreq + size(obs_wave)

      end if ! Force inversion wavelengths


      !
      ! Check for duplicates
      !

      ! Allocate a flag of valid frequencies
      allocate(flag(cfreq))
      flag = 1

      ! For each frequency
      do ifreq=1,cfreq

        ! If it has been flagged, we already checked
        if (flag(ifreq).lt.1) cycle

        ! If nosense, flag
        if (omega(ifreq).le.0d0) flag(ifreq) = 0

        ! Check the following ones
        do jfreq=ifreq+1,cfreq

          ! If it has been flagged, we already checked
          if (flag(jfreq).lt.1) cycle

          ! If some of them are repeated, flag them to be removed
          if(abs(1d2/omega(ifreq)-1d2/omega(jfreq)).lt.resol) &
            flag(jfreq) = 0

        end do ! jfreq
      end do ! ifreq

      ! Reset the running real index
      jfreq = 0

      ! For each frequency in the vector
      do ifreq=1,cfreq

        ! If it is flagged correct, add to real vector
        if(flag(ifreq).gt..5)then
          jfreq = jfreq + 1
          omega(jfreq) = omega(ifreq)
        end if

      end do

      ! The number of frequencies is the number of admitted
      ! frequencies in the flag vector
      cfreq = sum(flag)
      nfreq = cfreq

      ! Deallocate the flag
      deallocate(flag)

      ! Allocate the true axis
      ! Frequencies
      allocate(Frec%omega(nfreq))
      ! Integration weights
      allocate(Frec%W_freq(nfreq))
      ! Node weights
      allocate(Frec%IW_freq(nfreq))
      Frec%IW_freq = ContW
      ! If CLE
      if (run_mode.eq.2) then
        allocate(Frec%IW_freq_in(nfreq))
        Frec%IW_freq_in = 0
      end if

      ! Take only the valid frequencies and deallocate the auxiliar
      Frec%omega = omega(1:nfreq)
      deallocate(omega)

      ! Order the frequencies in the axis
      call QsortC(Frec%omega)


      !
      ! Check the presence of transitions at each frequency
      !

      ! Reset the ranges that the transitions spawns
      Frec%lif0 = 100000000
      Frec%lif1 = -1
      Frec%pif0 = 100000000
      Frec%pif1 = -1

      ! Factor of PRD in the core with respect to the wings
      prdf = Input%red_pars(5)/Input%red_pars(10)

      ! The default state of Yfield is False
      Yfield = .False.

      ! If there is PRD and calculating polarization, check if
      ! there is magnetic field
      if (PRD.and.lp) then

        Yfield = maxB.gt.TINYB

      end if ! PRD

      ! If there is PRD, calculate factors for the components
      ! of each transition as output
      if (PRD) then

        allocate(compfact(nxtran))
        compfact = 0d0
        mincf = 1d99

        ! For each atom
        do ia=1,nA

          ! If we have to run M components
          if (Yfield) then

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! If not PRD, ignore
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! Apply atomic shift
              ktran = jtran + Atom(ia)%tshift

              ! Identify terms
              do i=1,Atom(ia)%nMulti-1
                do i1=i+1,Atom(ia)%nMulti
                  if (Atom(ia)%irad(i,i1).eq.jtran) then
                    itermf = i
                    itermu = i1
                  end if
                end do
              end do

              ! Spin
              S = Atom(ia)%Sval(itermu)

              ! Orbital angular momentum
              rLu = Atom(ia)%rLval(itermu)
              rLf = Atom(ia)%rLval(itermf)

              ! Determine the maximum angular momentum and the
              ! number of magnetic sublevels for that maximum
              ! momentum
              rJumax = rLu + S
              nMu = nint(2d0*rJumax+1d0)
              rJfmax = rLf + S
              nMf = nint(2d0*rJfmax+1d0)

              ! For all the possible lower terms
              do i=1,Atom(ia)%nMulti-1

                ! If there is no transition or this term is larger
                ! than the upper term of the output transition, skip
                if (i.ge.itermu.or.Atom(ia)%irad(i,itermu).eq.0) &
                  cycle

                ! Store the input lower term index
                iterml = i

                ! Get index of input transition
                itran = Atom(ia)%irad(iterml,itermu)

                ! Angular momentum input lower level
                rLl = Atom(ia)%rLval(iterml)

                ! Determine maximum value of J and number of
                ! magnetic sublevels for this maximum J
                rJlmax = rLl + S
                nMl = nint(2d0*rJlmax+1d0)


        !
        ! Reset identation
        !

        ! For each Mf
        do iMf=1,nMf

          ! Value of Mf
          rMf = -rJfmax + dble(iMf-1)

          ! For each mu_f
          do mF=1,Atom(ia)%nblk(iMf,itermf)

            ! For each Mu
            do iMu=1,nMu

              ! Value of Mu
              rMu = -rJumax + dble(iMu-1)

              ! Difference between M momentums, done integer
              q = rMu - rMf
              iq = nint(q)

              ! If not pi nor sigma, skip
              if(abs(iq).gt.1) cycle

              ! For each mu_u
              do iU=1,Atom(ia)%nblk(iMu,itermu)

                ! For each Mu'
                do iMu1=1,nMu

                  ! Value of Mu'
                  rMu1 = -rJumax + dble(iMu1-1)

                  ! Difference between M momentums
                  q1 = rMu1-rMf

                  ! Convert to integers
                  iq1 = nint(q1)

                  ! If not pi or sigma, skip
                  if(abs(iq1).gt.1) cycle

                  ! For each mu_u'
                  do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                    ! For each Ml
                    do iMl=1,nMl

                      ! Value of Ml
                      rMl = -rJlmax + dble(iMl-1)

                      ! Difference between M momentums, in integer
                      p = rMu-rMl
                      ip = nint(p)

                      ! If not pi nor sigma, skip
                      if(abs(ip).gt.1) cycle

                      ! For each mu_l
                      do iL=1,Atom(ia)%nblk(iMl,iterml)

                        ! For each Ml'
                        do iMl1=1,nMl

                          ! Value of Ml'
                          rMl1 = -rJlmax + dble(iMl1-1)

                          ! Difference between M momentums
                          p1 = rMu1-rMl1

                          ! Convert to integer
                          ip1 = nint(p1)

                          ! If not pi nor sigma, skip
                          if(abs(ip1).gt.1) cycle

                          ! For each mu_l'
                          do iL1 = 1,Atom(ia)%nblk(iMl1,iterml)

                            ! Add to the weight
                            compfact(ktran) = compfact(ktran) + 1d0

                          end do ! iL1
                        end do ! iMl1
                      end do ! iL
                    end do ! iMl
                  end do ! iU1
                end do ! iMu1
              end do ! iU
            end do ! iMu
          end do ! mF
        end do ! iMF

                !
                ! Recover identation
                !

              end do ! Lower terms

              if (compfact(ktran).lt.mincf) mincf = compfact(ktran)

            end do ! Output transitions

          ! No need to run over M components
          else

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! If not PRD, ignore
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! Apply atomic shift
              ktran = jtran + Atom(ia)%tshift

              ! Identify terms
              do i=1,Atom(ia)%nMulti-1
                do i1=i+1,Atom(ia)%nMulti
                  if (Atom(ia)%irad(i,i1).eq.jtran) then
                    itermf = i
                    itermu = i1
                  end if
                end do
              end do

              ! Spin
              S = Atom(ia)%Sval(itermu)

              ! Orbital angular momentum
              rLu = Atom(ia)%rLval(itermu)
              rLf = Atom(ia)%rLval(itermf)

              ! For all the possible lower terms
              do i=1,Atom(ia)%nMulti-1

                ! If there is no transition or this term is larger
                ! than the upper term of the output transition, skip
                if (i.ge.itermu.or.Atom(ia)%irad(i,itermu).eq.0) &
                  cycle

                ! Store the input lower term index
                iterml = i

                ! Get index of input transition
                itran = Atom(ia)%irad(iterml,itermu)

                ! Angular momentum input lower level
                rLl = Atom(ia)%rLval(iterml)


        !
        ! Reset identation
        !

        ! For each Jf
        do mF=1,Atom(ia)%nJ(itermf)

          ! Get Jf
          rJf = Atom(ia)%rJval(mF,itermf)

          ! For each Ju
          do iU=1,Atom(ia)%nJ(itermu)

            ! Get Ju
            rJu = Atom(ia)%rJval(iU,itermu)

            if (nint(abs(rJu-rJf)).gt.1.or.(rJu+rJf).lt.dLSJ) cycle

            ! For each Ju'
            do iU1=1,Atom(ia)%nJ(itermu)

              ! Get Ju'
              rJu1 = Atom(ia)%rJval(iU1,itermu)

              if (nint(abs(rJu1-rJf)).gt.1.or.(rJu1+rJf).lt.dLSJ) &
                cycle

              ! For each Jl
              do iL=1,Atom(ia)%nJ(iterml)

                ! Get Jl
                rJl = Atom(ia)%rJval(iL,iterml)

                if (nint(abs(rJu-rJl)).gt.1.or.(rJu+rJl).lt.dLSJ) &
                  cycle

                ! For each Jl'
                do iL1=1,Atom(ia)%nJ(iterml)

                  ! Get Jl1
                  rJl1 = Atom(ia)%rJval(iL1,iterml)

                  if (nint(abs(rJu1-rJl1)).gt.1.or. &
                      (rJu1+rJl1).lt.dLSJ) cycle

                  ! Add to the weight
                  compfact(ktran) = compfact(ktran) + 1d0

                end do ! Jl'
              end do ! Jl
            end do ! Ju'
          end do ! Ju
        end do ! Jf

                !
                ! Recover identation
                !

              end do ! Lower terms

              if (compfact(ktran).lt.mincf) mincf = compfact(ktran)

            end do ! Output transitions

          end if ! M or not M components

        end do ! Atoms

        compfact = compfact/mincf

      end if ! If PRD

      ! Initialize size of profile variable
      Frec%ntfreq = 0
      Frec%ntfreqi = 0
      Frec%npfreq = 0

      ! For each atom
      do ia=1,nA

        ! Distance to consider PRD in the core
        prdc = Input%red_pars(7)*DwT(ia)

        ! Distance to consider coherent scattering in wing
        if (Input%cohw) prdcohw = Input%dcohw*DwTmin(ia)

        ! For each b-b transitions
        do itran=1,Atom(ia)%ntran

          ! Apply atomic shift
          ktran = itran + Atom(ia)%tshift

          ! Identify terms
          do i=1,Atom(ia)%nMulti-1
            do i1=i+1,Atom(ia)%nMulti
              if (Atom(ia)%irad(i,i1).eq.itran) then
                iterml = i
                itermu = i1
              end if
            end do
          end do

          ! Initialize limits
          Atom(ia)%if0(itran) = 100000000
          Atom(ia)%if1(itran) = -1

          ! Check for each frequency if it is within limits
          do ifreq=1,nfreq

            ! Line present
            if (Frec%omega(ifreq).ge.tmp_lim(1,ktran).and. &
                Frec%omega(ifreq).le.tmp_lim(2,ktran)) then

              ! Add weight to the node
              if (Atom(ia)%lemiss2(itran).and.PRD) then

                ! Initialize flag
                core = .False.

                ! Coherent wing?
                if (Input%cohw) then

                  ! Initialize flag
                  cohw = .True.

                  do iJu=1,Atom(ia)%nJ(itermu)
                    do iJl=1,Atom(ia)%nJ(iterml)

                      if (Atom(ia)%fst(itran)%irad(iJu,iJl).lt.1) &
                        cycle

                      v = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJl,iterml)

                      if (abs(v - Frec%omega(ifreq)).lt.prdcohw) then
                        cohw = .False.
                        exit
                      end if

                    end do
                    if (.not.cohw) exit
                  end do

                ! No coherent wing
                else

                  cohw = .False.

                end if ! Coherent wings?

                ! If not coherent wing
                if (.not.cohw) then

                  cmpf = compfact(ktran)

                  do iJu=1,Atom(ia)%nJ(itermu)
                    do iJl=1,Atom(ia)%nJ(iterml)

                      if (Atom(ia)%fst(itran)%irad(iJu,iJl).lt.1) &
                        cycle

                      v = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJl,iterml)

                      if (abs(v - Frec%omega(ifreq)).lt.prdc) then
                        core = .True.
                        exit
                      end if

                    end do
                    if (core) exit
                  end do

                end if ! Non-coherent wings

                if (core) then

                    if (nint(prdf*cmpf).eq.1) &
                      Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) - 1

                    Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + &
                                        nint(cmpf*prdf* &
                                             (LPRDW + LCRDW + 1))

                else if (cohw) then

                    if (nint(cmpf).eq.1) &
                      Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) - 1

                    Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + &
                                        nint(cmpf*(LCOHW + LCRDW + 1))

                else

                    if (nint(cmpf).eq.1) &
                      Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) - 1

                    Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + &
                                        nint(cmpf*(LPRDW + LCRDW + 1))

                end if

              ! CRD
              else

                Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + LCRDW

              end if

              ! If CLE
              if (run_mode.eq.2) &
                Frec%IW_freq_in(ifreq) = Frec%IW_freq_in(ifreQ) + &
                                         LCRDW

              ! Update lower limit
              if (ifreq.lt.Atom(ia)%if0(itran)) &
                Atom(ia)%if0(itran) = ifreq

              ! Update upper limit
              if (ifreq.gt.Atom(ia)%if1(itran)) &
                Atom(ia)%if1(itran) = ifreq

            end if

          end do

          ! If only one frequency, remove the transition
          if (Atom(ia)%if1(itran).le.Atom(ia)%if0(itran)) then
            Atom(ia)%fflag(itran)%absent = .True.
            Atom(ia)%if0(itran) = 1
            Atom(ia)%if1(itran) = 0
            Atom(ia)%W0(itran) = 0d0
            Atom(ia)%W1(itran) = 0d0
          else
            Atom(ia)%fflag(itran)%absent = .False.
          end if

          ! Store real Master frequency (for slaves)
          ! They will be overwritten later in general
          Atom(ia)%rif0(itran) = Atom(ia)%if0(itran)
          Atom(ia)%rif1(itran) = Atom(ia)%if1(itran)

          ! If absent, skip
          if (Atom(ia)%fflag(itran)%absent) cycle

          ! Compute weights for the limits
          Atom(ia)%W0(itran) = .5d5* &
                               (Frec%omega(Atom(ia)%if0(itran)+1) - &
                                Frec%omega(Atom(ia)%if0(itran)))
          Atom(ia)%W1(itran) = .5d5* &
                               (Frec%omega(Atom(ia)%if1(itran)) - &
                                Frec%omega(Atom(ia)%if1(itran)-1))

          ! Limits for the range with lines
          if (Atom(ia)%if0(itran).lt.Frec%lif0) Frec%lif0 = &
                                                   Atom(ia)%if0(itran)
          if (Atom(ia)%if1(itran).gt.Frec%lif1) Frec%lif1 = &
                                                   Atom(ia)%if1(itran)

          ! Add frequencies to count of profiles
          Frec%ntfreq = Frec%ntfreq + Atom(ia)%if1(itran) - &
                        Atom(ia)%if0(itran) + 1
          Frec%ntfreqi = Frec%ntfreqi + (Atom(ia)%if1(itran) - &
                         Atom(ia)%if0(itran) + 1)* &
                         Atom(ia)%fst(itran)%nt

        end do

        ! For each b-f transitions
        do itran=1,Atom(ia)%nphot

          ! Collect the maximum frequency
          nt = (1 - Atom(ia)%phot(itran)%mode)* &
               Atom(ia)%phot(itran)%nfreq + &
               Atom(ia)%phot(itran)%mode
          v = Atom(ia)%phot(itran)%infreq(nt)

          Atom(ia)%phot(itran)%absent = .False.

          ! Reset the ranges that the transition spawns
          Atom(ia)%phot(itran)%if0 = -1
          Atom(ia)%phot(itran)%if1 = -1

          ! If explicit
          if (Atom(ia)%phot(itran)%mode.eq.0) then

            ! For each frequency
            do ifreq=1,nfreq

              ! If we are below range, keep searching
              if (Frec%omega(ifreq).lt. &
                  Atom(ia)%phot(itran)%infreq(1)) cycle

              ! If we are above maximum, identify the index and go out
              if (Frec%omega(ifreq).gt.v) then
                Atom(ia)%phot(itran)%if1 = ifreq - 1
                exit
              end if

              ! If we are above the minimum, identify the index
              if (Atom(ia)%phot(itran)%if0.lt.0) then
                if (Frec%omega(ifreq).ge. &
                    Atom(ia)%phot(itran)%infreq(1)) &
                  Atom(ia)%phot(itran)%if0 = ifreq
              end if

            end do ! ifreq

          ! If hydrogenic
          else

            ! For each frequency
            do ifreq=1,nfreq

              ! If we are bwlor the edge, keep searching
              if (Frec%omega(ifreq).lt.Atom(ia)%phot(itran)%edge) &
                cycle

              ! If we are above maximum, identify the index and go out
              if (Frec%omega(ifreq).gt.v) then
                Atom(ia)%phot(itran)%if1 = ifreq - 1
                exit
              end if

              ! If we are above the minimum, identify the index
              if (Atom(ia)%phot(itran)%if0.lt.0) then
                if (Frec%omega(ifreq).ge.Atom(ia)%phot(itran)%edge) &
                  Atom(ia)%phot(itran)%if0 = ifreq
              end if

            end do ! ifreq

          end if ! Type of b-f transition

          ! If you didn't find the maximum, that means it was the
          ! last frequency itself
          if (Atom(ia)%phot(itran)%if1.lt.0) &
            Atom(ia)%phot(itran)%if1 = nfreq

          ! Limits for the range with photoionizations
          if (Atom(ia)%phot(itran)%if0.lt.Frec%pif0) Frec%pif0 = &
                                              Atom(ia)%phot(itran)%if0
          if (Atom(ia)%phot(itran)%if1.gt.Frec%pif1) Frec%pif1 = &
                                              Atom(ia)%phot(itran)%if1

          ! Add frequencies to count of profiles
          Frec%npfreq = Frec%npfreq + &
                        (Atom(ia)%phot(itran)%if1 - &
                         Atom(ia)%phot(itran)%if0 + 1)

          ! Add weight to the node
          do ifreq=Atom(ia)%phot(itran)%if0,Atom(ia)%phot(itran)%if1
            Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + PhotW
          end do

          ! If CLE
          if (run_mode.eq.2) then
            ! Add weight to the node
            do ifreq=Atom(ia)%phot(itran)%if0,Atom(ia)%phot(itran)%if1
              Frec%IW_freq_in(ifreq) = Frec%IW_freq_in(ifreq) + PhotW
            end do
          end if

        end do ! nphot

      end do ! nA

      ! Deallocate compfact if necessary
      if (PRD) deallocate(compfact)

      ! For each LTE line
      do ia=1,nLTEl

        ! Initialize limits
        Input%LTEline(ia)%if0 = nfreq+1
        Input%LTEline(ia)%if1 = -1

        ! Check for each frequency if it is within limits
        do ifreq=1,nfreq

          ! Line present
          if (Frec%omega(ifreq).ge.tmp_limL(1,ia).and. &
              Frec%omega(ifreq).le.tmp_limL(2,ia)) then

            Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + LLTEW

            ! Update lower limit
            if (ifreq.lt.Input%LTEline(ia)%if0) &
              Input%LTEline(ia)%if0 = ifreq

            ! Update upper limit
            if (ifreq.gt.Input%LTEline(ia)%if1) &
              Input%LTEline(ia)%if1 = ifreq

          end if

        end do ! Frequencies

        ! If only one frequency, remove the transition
        if (Input%LTEline(ia)%if1.le.Input%LTEline(ia)%if0) then
          Input%LTEline(ia)%absent = .True.
          Input%LTEline(ia)%if0 = 1
          Input%LTEline(ia)%if1 = 0
        else
          Input%LTEline(ia)%absent = .False.
        end if

      end do ! LTE lines


      !
      ! Protect frequencies against jumps due to big range
      !

      ! Allocate protection
      allocate(protect(Nfreq))
      protect = .False.

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          ! If between terms, we cannot protect
          if (Atom(ia)%fst(itran)%nt.gt.1) then
            cycle
          end if

          ! Protect FS lines
          do ifreq=Atom(ia)%if0(itran)+1,Atom(ia)%if1(itran)
            protect(ifreq) = .True.
          end do ! Each frequency
        end do ! b-b transitions

        ! For each b-f transitions
        do itran=1,Atom(ia)%nphot

          ! Protect ionizations
          do ifreq=Atom(ia)%phot(itran)%if0+1, &
                   Atom(ia)%phot(itran)%if1
            protect(ifreq) = .True.
          end do ! Each frequency
        end do ! b-b transitions
      end do ! Atom


      !
      ! Define the integration weights
      !

      ! The first point is special in compound trapezoidal rule
      Frec%W_freq(1) = .5d0*(Frec%omega(2) - Frec%omega(1))

      ! Initialize the integral to normalize the weights
      norm1 = Frec%W_freq(1)

      ! The initial lower limit is the first point
      O0 = Frec%omega(1)

      ! This is the pointer to the first element of the current
      ! interval, we are pointing to the first element
      cfreq = 1

      ! Flag that says that the point 2 is not the initial point
      ! of the interval (because 1 is the initial point)
      init = .FALSE.

      ! For the rest of frequencies except the last
      do ifreq=2,Nfreq-1

        ! If ifreq is the initial point of an interval
        if(init)then

          ! The first point is special in compound trapezoidal rule
          Frec%W_freq(ifreq) = .5d0*(Frec%omega(ifreq+1) - &
                                       Frec%omega(ifreq))

          ! The next point cannot be a first point
          init = .FALSE.

          ! Initialize the integral to normalize the weights
          norm1 = Frec%W_freq(ifreq)

          ! Pointer is now in this frequency
          cfreq = ifreq

          ! And it is the beginning of the current interval
          O0 = Frec%omega(ifreq)

        ! If ifreq is not the initial point of an interval
        else

          if(abs(1d2/Frec%omega(ifreq+1) - &
             1d2/Frec%omega(ifreq)).gt.jump.and. &
             .not.protect(ifreq+1))then

            ! Notify of this jump
            if (gpid.eq.0) then
              write(umsg,'(A,1x,f10.3,1x,A,1x,f10.3,1x,A)') &
                  ' # Wavelengths ',1d2/Frec%omega(ifreq+1),'and', &
                  1d2/Frec%omega(ifreq),'are considered to be'// &
                  ' in different ranges.'
              call verbose
            end if

            ! The last point is special in compound trapezoidal rule
            Frec%W_freq(ifreq) = .5d0*(Frec%omega(ifreq) - &
                                         Frec%omega(ifreq-1))

            ! It is the end of the current interval
            O1 = Frec%omega(ifreq)

            ! Add to the integral
            norm1 = norm1 + Frec%W_freq(ifreq)

            ! We know that the integral must be
            ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
            norm = 1d5*(O1 - O0)

            ! The normalizing factor is thus
            norm = norm/norm1

            ! Normalize the weights of this interval
            do jfreq=cfreq,ifreq
              Frec%W_freq(jfreq) = Frec%W_freq(jfreq)*norm
            end do

            ! The next point is the first point of its interval
            init = .TRUE.

          ! If ifreq is not the last point of an interval
          else

            ! Check if ifreq is the last point of an interval
            if(abs(1d2/Frec%omega(ifreq+1) - &
                   1d2/Frec%omega(ifreq)).gt.jump.and. &
               protect(ifreq+1).and.gpid.eq.0) then

              write(umsg,'(A,1x,f10.3,1x,A,1x,f10.3,1x,A)') &
                  ' # Wavelengths ',1d2/Frec%omega(ifreq+1),'and', &
                  1d2/Frec%omega(ifreq),'are considered to be'// &
                  ' in different ranges, but they are '// &
                  'protected from the splitting.'
              call verbose
            end if

            ! Compound trapezoidal rule weight
            Frec%W_freq(ifreq) = .5d0*(Frec%omega(ifreq+1) - &
                                         Frec%omega(ifreq-1))

            ! Add to the integral
            norm1 = norm1 + Frec%W_freq(ifreq)

          end if ! Last point

        end if ! Initial point

      end do ! ifreq

      ! The last point is special in compound trapezoidal rule
      Frec%W_freq(Nfreq) = .5d0*(Frec%omega(Nfreq) - &
                                   Frec%omega(Nfreq-1))

      ! It is the end of the interval
      O1 = Frec%omega(Nfreq)

      ! Add to the integral
      norm1 = norm1 + Frec%W_freq(Nfreq)

      ! We know that the integral must be
      ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
      norm = 1d5*(O1 - O0)

      ! The normalizing factor is thus
      norm = norm/norm1

      ! Normalize the weights of this interval
      do ifreq=cfreq,nfreq
        Frec%W_freq(ifreq) = Frec%W_freq(ifreq)*norm
      end do

      ! Initialize ndzao
      Frec%ndzao = 0

      ! Deallocate arrays
      deallocate(tmp_lim)
      deallocate(protect)

      ! Check if everything is fine
      call control

      return

1000  write(umsg,'(A,1x,i2)') 'Error opening wave file',i
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  write(umsg,'(A,1x,i2)') 'Error reading wave file',i
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine omegabuild

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize Master limits\n
      !!    Atom(Atom_class): Structure with the atomic data
      subroutine omegainitmaster(Atom)

      ! I/O

      type(Atom_class), dimension(:):: Atom

      ! Local
      integer:: ia,itran

      ! Only the Máster with MPI
      if (pid.eq.0.and.nproc.gt.1) then

        ! For each atom
        do ia=1,nA

          ! For each b-b transition
          do itran=1,Atom(ia)%ntran

            ! Allocate the presence flag for master
            allocate(Atom(ia)%fflag(itran)%Mabsent(0:nproc-1))
            Atom(ia)%fflag(itran)%Mabsent = .False.

            ! If only one frequency, remove the transition
            if (Atom(ia)%if1(itran).le.Atom(ia)%if0(itran)) then
              ! Initialize master
              if (pid.eq.0.and.nproc.gt.1) then
                Atom(ia)%Mif0(itran,:) = Atom(ia)%if0(itran)
                Atom(ia)%Mif1(itran,:) = Atom(ia)%if1(itran)
                Atom(ia)%MW0(itran,:) = Atom(ia)%W0(itran)
                Atom(ia)%MW1(itran,:) = Atom(ia)%W1(itran)
              end if
            end if

            ! If absent
            if (Atom(ia)%fflag(itran)%absent) cycle

            ! Initialize master
            Atom(ia)%Mif0(itran,:) = Atom(ia)%if0(itran)
            Atom(ia)%Mif1(itran,:) = Atom(ia)%if1(itran)
            Atom(ia)%MW0(itran,:) = Atom(ia)%W0(itran)
            Atom(ia)%MW1(itran,:) = Atom(ia)%W1(itran)

          end do ! b-b transitions

          ! For each b-f tarnsition
          do itran=1,Atom(ia)%nphot

            ! Allocate the presence flag for master
            allocate(Atom(ia)%phot(itran)%Mabsent(0:nproc-1))
            Atom(ia)%phot(itran)%Mabsent = .False.

            ! Initialize master
            Atom(ia)%phot(itran)%Mif0 = Atom(ia)%phot(itran)%if0
            Atom(ia)%phot(itran)%Mif1 = Atom(ia)%phot(itran)%if1

          end do ! b-f transitions
        end do ! Atoms

      end if ! Máster

      end subroutine omegainitmaster

!#####################################################################
!#####################################################################
!#####################################################################

      !> Build the input frequency axis.\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bstrength(dfloat(:)): Magnetic field strength\n
      !!        Input(Input_class): Structure with settings data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!               lp(logical): If doing formal solution in this
      !!                            run\n
      !!            ofram(logical): Indicates if out of RAM\n
      !!              LOS(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine omegabuildin(Frec,Red,Atom,Atmo,Bstrength,Input, &
                              Geom,MPID,lp,ofram,LOS)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(Input_class):: Input
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class):: MPID
      logical, intent(in):: lp,LOS
      logical, intent(out):: ofram
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: lskip,skip,init,reset,nfound,core,Yfield,Nfield
      logical:: RAMOF, lNCHLT, cohw
      logical, dimension(:), allocatable:: WNCHLT

      integer:: ios,i,i1,i2,it,iz,ia,ip,ipp,ir,itran,jtran,ktran
      integer:: ifreq,jfreq,iifreq,lifreq,kfreq,cfreq,if0,if1,bf0,bf1
      integer:: iJl,iJu,iJf,itermf,itermu,iterml,iran,nran,ibfreq
      integer:: jdir,jbdir,jcdir,njdir,iYYF,iYNF,iNF,iDF,iDFR
      integer:: nr,nt,ntj,ni,nie,np,np0,npp,nti,iti,ith,iph,ith1,iph1
      integer:: minto,maxto,mint,maxt,nat,mina,maxa
      integer:: jjfreq,jjfreq0,kkfreq,kkfreq0,jufreq
      integer:: nMm,nMu,nMl,nMf,nblk,indx,indxf,nn
      integer:: iMf,mF,iMu,iU,iMu1,iU1,iMl,iL,iMl1,iL1
      integer:: iq,iq1,ip1,iQQ
      integer:: indU,indU1,indF,indL,indL1
      integer:: MiindU,MiindU1,MiindF,MiindL,MiindL1
      integer:: MindU,MindU1,MindF,MindL,MindL1
#ifdef _OPENMP
      integer:: tid
#endif
      integer, dimension(:), allocatable:: minti,maxti
      integer, dimension(:), allocatable:: flag,ithV,iphV
#ifdef _OPENMP
      integer, dimension(:), allocatable:: oif0,oif1,orif0,orif1,onf
#endif

      double precision:: DwT,Dw,Dw1,vph,vpl,vfac,vfac1
      double precision:: norm,norm1,O0,O1,dnlmin,dnlmax
      double precision:: red_rangW1,red_resoW
      double precision:: red_vlarW1,red_fstpW1
      double precision:: red_mstpW1,red_neglW
      double precision:: red_rangwW1,red_vlarwW1
      double precision:: red_fstpwW1,red_mstpwW1
      double precision:: red_rangcW1,red_vlarcW1
      double precision:: red_fstpcW1,red_mstpcW1
      double precision:: red_coreW,red_cohwW
      double precision:: rJumax,rJlmax,rJfmax
      double precision:: rJu,rJu1,rJl,rJl1,rJf
      double precision:: rMf,rMu,rMu1,rMl,rMl1
      double precision:: q,q1,p,p1,QQ,PP
      double precision:: ct,st,cc,sc,ThK,ct1,st1,cc1,sc1,SRAM
      double precision, dimension(:), allocatable:: dnl, nut
      double precision, dimension(:), allocatable:: vphv, vplv, vpr
      double precision, dimension(:), allocatable:: vphve, vplve
      double precision, dimension(:), allocatable:: vpp
      double precision, dimension(:), allocatable:: Wvpp

      ! Box
      type(dbabox_class), pointer:: bomega, bw_freq, bdaux

      ! Pointer
      type(Frequencyd_class), pointer:: p_frec

      ! Initialize
      nullify(bomega,bw_freq,bdaux)

      ! Routine name
      urou = 'omegabuildin'

      ! Initialize
      ofram = .False.

      ! Check if we don't have to reset
      if (AV.and..not.dyn.and.LOS) then
        if (lp) return
      end if

      ! If line of sight
      if ((.not.AV.or.dyn).and.LOS) then
        PRAM = .False.
        IRAM = .False.
      end if

      ! Inititlize RAM counters
      MPID%WRAM = 0d0

      ! Deallocate existing Frec and Red
      call cleanFrecandRed(Frec,Red,MPID)

      ! If the process is the master, it does not need this
      if (MPID%mpi.and.pid.eq.0) then
        call control
        return
      end if

      Yfield = .False.
      Nfield = .False.
      RAMOF = .False.

      ! Check magnetic field
      do iz=Rz0,Rz1
        if (Bstrength(iz).lt.TINYB) then
          Nfield = .True.
        else
          Yfield = .True.
        end if
        if (Yfield.and.Nfield) exit
      end do

      ! If angle dependent
      if (.not.AV) then

        !
        ! Output and input directions

        ! If for emergence
        if (LOS) then

          njdir = Geom%nPhLOS*Geom%nThLOS

        ! If quadrature
        else

          njdir = Geom%nPh*Geom%nTh

        end if

        ! If dynamic, axis changes
        if (dyn) then

          Frec%ndir = njdir
          Frec%nth = Geom%nTh

          ! Azimuth depends on symmetry
          if (axial) then
            Frec%nph = 1
          else
            Frec%nph = Geom%nPh
          end if

        ! If static, they don't
        else

          Frec%ndir = 1
          Frec%nth = 1
          Frec%nph = 1

        end if

        Red%ndir = Geom%nPh2*Geom%nTh
        Red%njdir = njdir
        Red%nth = Geom%nTh
        Red%nph = Geom%nPh2


        ! Allocate auxiliar quantities for Doppler shifts

        ! Index of polar direction of quadrature
        allocate(ithv(njdir))
        ! Index of azimuthal direction of quadrature
        allocate(iphv(njdir))

        ! Allocate type of scattering
        if (allocated(Frec%stype)) deallocate(Frec%stype)
        if (allocated(Frec%nfs)) deallocate(Frec%nfs)
        allocate(Frec%stype(Geom%nPh2,Geom%nTh,njdir))
        allocate(Frec%nfs(njdir))

        ! Initialize to normal scattering
        Frec%stype = 0
        Frec%nfs = 0


        !
        ! De-index the directions

        jdir = 0

        ! If lines of sight
        if (LOS) then

          do ith=1,Geom%nThLOS
            do iph=1,Geom%nPhLOS

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

              ! If angle-dependent
              if (.not.AV) then

                do ith1=1,Geom%nTh
                  do iph1=1,Geom%nPh2

                    ! Calculate scattering angle between the
                    ! quadrature direction and the LOS direction
                    ThK = atom2lab(Geom%L_theta(ith), &
                                   Geom%L_phi(iph), &
                                   Geom%V_theta(ith1), &
                                   Geom%V_phi(iph1))

                    ! Backward condition
                    if (abs(pi - ThK).le.TINYA) then

                      Frec%stype(iph1,ith1,jdir) = 1

                    ! Forward condition
                    else if (ThK.le.TINYA) then

                      Frec%stype(iph1,ith1,jdir) = -1
                      Frec%nfs(jdir) = Frec%nfs(jdir) + 1

                    end if ! Backward or forward scattering

                  end do ! Quadrature directions (input)
                end do

              end if ! AD redistribution

            end do ! LOS directions
          end do

        ! If quadrature
        else

          do ith=1,Geom%nTh
            do iph=1,Geom%nPh

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

              ! If angle-dependent
              if (.not.AV) then

                do ith1=1,Geom%nTh
                  do iph1=1,Geom%nPh2

                    ! Backward condition
                    if ((Geom%nTh - ith + 1).eq.ith1.and. &
                        int(iph + .5d0*Geom%V_muy(iph)* &
                            Geom%nPh2).eq.iph1) then

                      Frec%stype(iph1,ith1,jdir) = 1

                    ! Forward condition
                    else if (ith1.eq.ith.and. &
                             iph1.eq.iph) then

                      Frec%stype(iph1,ith1,jdir) = -1
                      Frec%nfs(jdir) = Frec%nfs(jdir) + 1

                    end if ! Backward or forward scattering

                  end do ! Quadrature directions (input)
                end do

              end if ! AD redistribution

            end do ! Quadrature directions
          end do

        end if ! LOS

        ! If static, initialize velocity factors
        if (.not.dyn) then

          vfac = 1d0
          vfac1 = 1d0

        end if

      ! If dynamic (and AA)
      else if (dyn) then

        !
        ! Output and input directions

        ! If lines of sight
        if (LOS) then
          njdir = Geom%nPhLOS*Geom%nThLOS
        ! If quadrature
        else
          njdir = Geom%nPh*Geom%nTh
        end if
        Frec%ndir = njdir
        Frec%nth = 1
        Frec%nph = 1

        Red%ndir = Geom%nPh*Geom%nTh
        Red%njdir = 1
        Red%nth = 1
        Red%nph = 1

        ! Allocate auxiliar quantities for Doppler shifts

        ! Index of polar direction of quadrature
        allocate(ithv(njdir))
        ! Index of azimuthal direction of quadrature
        allocate(iphv(njdir))

        ! Unneccessary type of scattering
        if (allocated(Frec%stype)) deallocate(Frec%stype)
        if (allocated(Frec%nfs)) deallocate(Frec%nfs)
        allocate(Frec%stype(1,1,1))
        allocate(Frec%nfs(1))

        ! Initialize to normal scattering
        Frec%stype = 0
        Frec%nfs = 0

        !
        ! De-index the directions

        jdir = 0

        ! LOS
        if (LOS) then

          do ith=1,Geom%nThLOS
            do iph=1,Geom%nPhLOS

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

            end do ! LOS directions
          end do

        ! Quadrature
        else

          do ith=1,Geom%nTh
            do iph=1,Geom%nPh

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

            end do ! Quadrature directions
          end do

        end if ! LOS

      else

        ! Output directions
        Frec%ndir = 1
        Red%ndir = 1

        ! Input directions
        Frec%nth = 1
        Frec%nph = 1
        Red%nth = 1
        Red%nph = 1
        Red%njdir = 1

        ! Doppler factors
        vfac = 1d0
        vfac1 = 1d0

        ! Allocate stype to avoid undefined
        if (PRD) then

          ! Allocate type of scattering
          if (allocated(Frec%stype)) deallocate(Frec%stype)
          if (allocated(Frec%nfs)) deallocate(Frec%nfs)
          allocate(Frec%stype(1,1,1))
          allocate(Frec%nfs(1))

          ! Initialize to normal scattering
          Frec%stype = 0
          Frec%nfs = 0

        end if

      end if

      !
      ! Count maximum index of PRD atom and transition

      ! Initialize atomic and transition index, and counter
      ! of real elements
      mina = 10000
      maxa = 0
      minto = 10000
      maxto = 0
      nat = 0

      ! For each atom
      do ia=1,nA
        ! For all transitions
        do jtran=1,Atom(ia)%ntran

          ! If PRD line
          if (Atom(ia)%lemiss2(jtran)) then
            ! If not absent in this CPU
            if (.not.Atom(ia)%fflag(jtran)%absent) then

              ! Add to counter
              nat = nat + 1

              ! Update limits
              if (jtran.lt.minto) minto = jtran
              if (jtran.gt.maxto) maxto = jtran
              if (ia.lt.mina) mina = ia
              if (ia.gt.maxa) maxa = ia

            end if ! Presence of line
          end if ! PRD line

        end do ! Transitions
      end do ! Atoms

      ! Allocate indexing array and first step of Frec and Red
      allocate(Frec%indx(minto:maxto,mina:maxa,Rz0:Rz1,Frec%ndir))
      Frec%ndzao = nat*Rnz*Frec%ndir
      allocate(Frec%dzao(Frec%ndzao))
      do indx=1,Frec%ndzao
        nullify(Frec%dzao(indx)%trani)
      end do
      if (PRAM) then
        allocate(Red%indx(minto:maxto,mina:maxa,Rz0:Rz1,Red%ndir))
        Red%ndzao = nat*Rnz*Red%ndir
        allocate(Red%dzao(Red%ndzao))
        do indx=1,Red%ndzao
          nullify(Red%dzao(indx)%trani)
        end do
      end if

      !
      ! Build index

      ! Initialize
      ip = 0

      ! Directions
      do jdir=1,Frec%ndir
        ! Height
        do iz=Rz0,Rz1
          ! Atom
          do ia=mina,maxa
            ! Transition
            do jtran=1,Atom(ia)%ntran

              ! If PRD line
              if (Atom(ia)%lemiss2(jtran).and. &
                  .not.Atom(ia)%fflag(jtran)%absent) then
                ip = ip + 1
                Frec%indx(jtran,ia,iz,jdir) = ip
              end if

            end do ! Transition
          end do ! Atom
        end do ! Height
      end do ! Directions

#ifdef _OPENMP
      ! If multiple threads
      if (omp) then

        ! Allocate limits for threads
        allocate(oif0(nthread),oif1(nthread),onf(nthread))

        ! Work per thread
        ios = Frec%ndzao/nthread

        ! Give this first stimation to each process
        do tid=1,nthread
          onf(tid) = ios
        end do

        ! Put the remaining heights in some of the threads if there
        ! are remaining ones
        if(ios*nthread.ne.Frec%ndzao)then

          ! Number of nodes to distribute
          ios = Frec%ndzao - ios*nthread

          ! Give them to the first threads
          do tid=1,ios
            onf(tid) = onf(tid) + 1
          end do

        end if

        ! Set first thread boundary
        oif0(1) = 1
        oif1(1) = onf(1)

        ! For each other thread
        do tid=2,nthread
          oif0(tid) = oif1(tid-1) + 1
          oif1(tid) = oif0(tid) + onf(tid) - 1
        end do

        ! Deallocate onf
        deallocate(onf)

      end if ! multithreaded
#endif
      ! If storing redistribution
      if (PRAM) then

        ! Initialize
        ip = 0

        ! Directions
        do jdir=1,Red%ndir
          ! Height
          do iz=Rz0,Rz1
            ! Atom
            do ia=mina,maxa
              ! Transition
              do jtran=1,Atom(ia)%ntran

                ! If PRD line
                if (Atom(ia)%lemiss2(jtran).and. &
                    .not.Atom(ia)%fflag(jtran)%absent) then
                  ip = ip + 1
                  Red%indx(jtran,ia,iz,jdir) = ip
                end if

              end do ! Transition
            end do ! Atom
          end do ! Height
        end do ! Directions
#ifdef _OPENMP
        ! If threading
        if (omp) then

          ! Allocate limits for threads
          allocate(orif0(nthread),orif1(nthread),onf(nthread))

          ! Work per thread
          ios = Red%ndzao/nthread

          ! Give this first stimation to each process
          do tid=1,nthread
            onf(tid) = ios
          end do

          ! Put the remaining heights in some of the threads if there
          ! are remaining ones
          if(ios*nthread.ne.Red%ndzao)then

            ! Number of nodes to distribute
            ios = Red%ndzao - ios*nthread

            ! Give them to the first threads
            do tid=1,ios
              onf(tid) = onf(tid) + 1
            end do

          end if

          ! Set first thread boundary
          orif0(1) = 1
          orif1(1) = onf(1)

          ! For each other thread
          do tid=2,nthread
            orif0(tid) = orif1(tid-1) + 1
            orif1(tid) = orif0(tid) + onf(tid) - 1
          end do

          ! Deallocate onf
          deallocate(onf)

        end if ! Multi-thread
#endif
      end if ! PRAM

!$omp parallel default(none) &
!$omp private(np0,vpp,flag,tid,ia,skip,jtran,minto,maxto,nt,itermu) &
!$omp private(itermf,iterml,itran,ktran,i,i1,nti,jdir,iz,indx,i2) &
!$omp private(iti,nr,iJl,iJf,dnl,dnlmax,dnlmin,iJu,ntj,nut,ni) &
!$omp private(vphv,vplv,vphve,vplve,vpr,ith1,iph1,ct,st,cc,sc) &
!$omp private(vfac,if0,if1,DwT,Dw1,Dw,red_resoW,red_neglW) &
!$omp private(red_coreW,red_rangwW1,red_vlarwW1,red_fstpwW1) &
!$omp private(red_mstpwW1,red_rangcW1,red_vlarcW1,red_fstpcW1) &
!$omp private(red_mstpcW1,nran,bf0,bf1,lskip,ifreq,np,bomega) &
!$omp private(bw_freq,iran,iifreq,core,it,red_rangW1,red_vlarW1) &
!$omp private(red_fstpW1,red_mstpW1,vph,vpl,ir,nie,reset,npp) &
!$omp private(nfound,ip,ipp,bdaux,Wvpp,norm1,O0,cfreq,init,O1) &
!$omp private(norm,nn,p_frec,SRAM,jjfreq0,kkfreq0,ith,iph,ct1,vfac1) &
!$omp private(st1,cc1,sc1,jjfreq,kkfreq,lifreq,jfreq,jufreq,ibfreq) &
!$omp private(nMm,nblk,rJumax,nMu,rJfmax,nMf,ios,rJlmax,nMl) &
!$omp private(MindU,MindU1,MindL,MindL1,MindF,iMf,rMf,mF,indF,iMu) &
!$omp private(rMu,q,iq,iU,indU,iMu1,rMu1,q1,QQ,iq1,iQQ,iU1,indU1) &
!$omp private(iMl,rMl,p,iL,indL,iMl1,rMl1,p1,PP,ip1,iL1) &
!$omp private(indL1,MiindU,MiindU1,MiindL,MiindL1,MiindF) &
!$omp private(rJf,rJu,rJu1,rJl,rJl1,jbdir,jcdir,indxf,iDf,iDFR) &
!$omp private(red_cohwW,cohw) &
!$omp shared(nA,Atom,minti,maxti,mint,maxt,Input,Frec,nZ,omp,oif0) &
!$omp shared(oif1,PRAM,Red,orif0,orif1,dyn,LOS,Geom,MPID,AV,axial) &
!$omp shared(TPRAM,Yfield,NCHLT,lNCHLT,iYNF,WNCHLT,iYYF,Nfield,iNF) &
!$omp shared(Bstrength,ofram,ithv,iphv,Atmo,RLIM,nfreq,Rz0,Rz1)

      ! Preliminar allocation of vpp, auxiliar vector to store
      ! frequencies
      np0 = 10000
      allocate(vpp(np0))
      allocate(flag(np0))

#ifdef _OPENMP
      vfac = 1d0
      vfac1 = 1d0
      tid = omp_get_thread_num() + 1
#endif

      ! For each atom
      do ia=1,nA

        ! Check that there is at least 1 PRD line in this atom
        skip = .True.

        do jtran=1,Atom(ia)%ntran
          if (Atom(ia)%lemiss2(jtran)) then
            if (.not.Atom(ia)%fflag(jtran)%absent) then
              skip = .False.
              exit
            end if
          end if
        end do

        ! There are no PRD lines for this atom
        if (skip) cycle

!$omp barrier

!$omp single
        ! Allocate trano
        if (allocated(Atom(ia)%trano)) deallocate(Atom(ia)%trano)

        ! Allocate itrano
        if (allocated(Atom(ia)%itrano)) deallocate(Atom(ia)%itrano)

        ! Reallocate minti and maxti
        if (allocated(minti)) deallocate(minti,maxti)
        allocate(minti(Atom(ia)%ntran))
        allocate(maxti(Atom(ia)%ntran))


        !
        ! Count the transition combinations
        !

        ! Reset index
        minto = Atom(ia)%ntran + 1
        maxto = 0
        minti = Atom(ia)%ntran + 1
        maxti = 0
        mint = Atom(ia)%nMulti + 1
        maxt = 0
        nt = 0

        ! For each upper term
        do itermu=2,Atom(ia)%nMulti

          ! For each final lower term
          do itermf=1,itermu-1

            jtran = Atom(ia)%irad(itermu,itermf)

            if (jtran.le.0) cycle
            if (.not.Atom(ia)%lemiss2(jtran)) cycle
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Update size
            nt = nt + 1

            if (jtran.lt.minto) minto = jtran
            if (jtran.gt.maxto) maxto = jtran

            ! For each other lower term
             do iterml=1,itermu-1

               itran = Atom(ia)%irad(itermu,iterml)

               if (itran.le.0) cycle

               if (.not.Input%Raman.and.itran.ne.jtran) cycle

               ! Look for the maximum and minimum indexes
               if (itermu.lt.mint) mint = itermu
               if (itermf.lt.mint) mint = itermf
               if (iterml.lt.mint) mint = iterml
               if (itermu.gt.maxt) maxt = itermu
               if (itermf.gt.maxt) maxt = itermf
               if (iterml.gt.maxt) maxt = iterml

               if (itran.lt.minti(jtran)) minti(jtran) = itran

               if (itran.gt.maxti(jtran)) maxti(jtran) = itran

            end do ! iterml
          end do ! itermf
        end do ! itermu

        ! Store amount of output transitions
        Atom(ia)%ntrano = nt

        ! Index trano
        allocate(Atom(ia)%itrano(minto:maxto))
        nt = 0

        ! For each upper term
        do itermu=2,Atom(ia)%nMulti

          ! For each final lower term
          do itermf=1,itermu-1

            jtran = Atom(ia)%irad(itermu,itermf)

            if (jtran.le.0) cycle
            if (.not.Atom(ia)%lemiss2(jtran)) cycle
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Update index
            nt = nt + 1
            Atom(ia)%itrano(jtran) = nt

          end do ! itermf
        end do ! itermu


        !
        ! Allocate needed transition structures
        !

        ! Atomic structure
        allocate(Atom(ia)%trano(Atom(ia)%ntrano))
!$omp end single

        ! For each output transition
        do jtran = 1,Atom(ia)%ntran

          ! If no PRD line, skip
          if (.not.Atom(ia)%lemiss2(jtran)) cycle
          if (Atom(ia)%fflag(jtran)%absent) cycle

          ! structure index
          ktran = Atom(ia)%itrano(jtran)

          ! Look for the terms
          itermu = -1
          do i=1,Atom(ia)%nMulti-1
            do i1=i+1,Atom(ia)%nMulti
              if (Atom(ia)%irad(i,i1).eq.jtran) then
                itermf = i
                itermu = i1
                exit
              end if
            end do
            if (itermu.ge.0) exit
          end do

          ! Number of input transitions
          nti = 0

          ! For each other transition that shares upper term
          do i2=1,Atom(ia)%nMulti-1

            if(i2.ge.itermu.or.Atom(ia)%irad(i2,itermu).eq.0) cycle
            iterml = i2

            ! Input transition
            itran = Atom(ia)%irad(iterml,itermu)

            ! Skip if no Raman and is not Rayleigh
            if (.not.Input%Raman.and.itran.ne.jtran) cycle

            ! Advance the continuous index
            nti = nti + 1

          end do ! Lower transitions

          ! Allocate atomic structure
!$omp single
          allocate(Atom(ia)%trano(ktran)%trani(nti))
          allocate(Atom(ia)%trano(ktran)%ind(nti))
          Atom(ia)%trano(ktran)%nt = nti

          ! Reset index
          Atom(ia)%trano(ktran)%ind = -1
!$omp end single

          ! For each output direction
          do jdir=1,Frec%ndir

            ! For each height
            do iz=Rz0,Rz1

              ! Get index
              indx = Frec%indx(jtran,ia,iz,jdir)
#ifdef _OPENMP
              ! If multi-thread
              if (omp) then
                ! Skip if not assigned
                if (indx.lt.oif0(tid)) cycle
                if (indx.gt.oif1(tid)) exit
              end if
#endif
              ! Allocate the input transition type
              if(.not.associated(Frec%dzao(indx)%trani)) then

                ! Allocate
                allocate(Frec%dzao(indx)%trani(nti))

                ! Nullify pointers
                do i2=1,nti
                  nullify(Frec%dzao(indx)%trani(i2)%index1)
                  nullify(Frec%dzao(indx)%trani(i2)%index2)
                  nullify(Frec%dzao(indx)%trani(i2)%dx)
                end do

              end if ! No input transitions associated yet

              ! Initialize the number of input frequencies
              Frec%dzao(indx)%mxfreq = 0

            end do ! heights
          end do ! output directions

          ! If storing Warr
          if (PRAM) then

            ! For each output direction
            do jdir=1,Red%ndir

              ! For each height
              do iz=Rz0,Rz1

                ! Get index
                indx = Red%indx(jtran,ia,iz,jdir)
#ifdef _OPENMP
                ! If multi-thread
                if (omp) then
                  ! Skip if not assigned
                  if (indx.lt.orif0(tid)) cycle
                  if (indx.gt.orif1(tid)) exit
                end if
#endif
                ! Allocate the input transition type
                if(.not.associated(Red%dzao(indx)%trani)) &
                  allocate(Red%dzao(indx)%trani(nti))

              end do ! heights
            end do ! output directions

          end if ! Storing PRAM
!$omp single
          ! Reset index input transition
          iti = 0

          ! For each other transition that shares upper term
          do i2=1,Atom(ia)%nMulti-1

            if(i2.ge.itermu.or.Atom(ia)%irad(i2,itermu).eq.0) cycle
            iterml = i2

            ! Input transition
            itran = Atom(ia)%irad(iterml,itermu)

            ! Skip if no Raman and is not Rayleigh
            if (.not.Input%Raman.and.itran.ne.jtran) cycle

            ! Advance the continuous index
            iti = iti + 1

            ! Store index of this input transition
            Atom(ia)%trano(ktran)%ind(iti) = itran

          end do ! Lower transitions
!$omp end single

          ! Reset index of input transition
          iti = 0

          ! For each other transition that shares upper term
          do i2=1,Atom(ia)%nMulti-1

            if(i2.ge.itermu.or.Atom(ia)%irad(i2,itermu).eq.0) cycle
            iterml = i2

            ! Input transition
            itran = Atom(ia)%irad(iterml,itermu)

            ! Skip if no Raman and is not Rayleigh
            if (.not.Input%Raman.and.itran.ne.jtran) cycle

            ! Advance the continuous index
            iti = iti + 1


            !
            ! Count the number of resonances
            !

            ! Reset the index
            nr = 0

            ! We have as many resonances as combinations between
            ! lower levels (the repeated ones will be removes later)
            do iJl=1,Atom(ia)%nJ(iterml)
              do iJf=1,Atom(ia)%nJ(itermf)
                nr = nr + 1
              end do
            end do

            ! Make sure that we have enough space to store resonances
            if(.not.allocated(dnl))then
              allocate(dnl(nr))
            else
              if(size(dnl).lt.nr)then
                deallocate(dnl)
                allocate(dnl(nr))
              end if
            end if

            ! Reset number of resonances
            nr = 0

            ! Reset ranges
            dnlmax = -1D99
            dnlmin = 1D99

            ! For each pair of lower levels for input and output
            do iJl=1,Atom(ia)%nJ(iterml)
              do iJf=1,Atom(ia)%nJ(itermf)

                ! Add the resonance
                nr = nr + 1
                dnl(nr) = Atom(ia)%FSfreq(iJf,itermf) - &
                          Atom(ia)%FSfreq(iJl,iterml)

                ! And check that we know what are the futhest ones
                if(dnl(nr).gt.dnlmax)dnlmax = dnl(nr)
                if(dnl(nr).lt.dnlmin)dnlmin = dnl(nr)

              end do ! iJl
            end do ! iJf

            !
            ! Count number of transitions
            !

            ! Reset the index
            nt = 0

            ! Count the FS transitions in the output line
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJf=1,Atom(ia)%nJ(itermf)
                if(abs(Atom(ia)%rJval(iJu,itermu) - &
                       Atom(ia)%rJval(iJf,itermf)).gt.1.or. &
                   Atom(ia)%rJval(iJu,itermu) + &
                   Atom(ia)%rJval(iJf,itermf).lt..4d0) cycle
                nt = nt + 1
              end do
            end do

            ! This is the number of transitions that are in the output
            ntj = nt

            ! Count the FS transitions in the input line
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJl=1,Atom(ia)%nJ(iterml)
                if(abs(Atom(ia)%rJval(iJu,itermu) - &
                       Atom(ia)%rJval(iJl,iterml)).gt.1.or. &
                   Atom(ia)%rJval(iJu,itermu) + &
                   Atom(ia)%rJval(iJl,iterml).lt..4d0)cycle
                nt = nt + 1
              end do
            end do

            ! Make sure that we have enough space to store frequencies
            if(.not.allocated(nut))then
              allocate(nut(nt))
            else
              if(size(nut).lt.nt)then
                deallocate(nut)
                allocate(nut(nt))
              end if
            end if

            ! The total number of limits to consider is the sum of
            ! lines and resonances
            ni = nr + nt


            ! Make sure that we have enough space to store ranges
            ! Upper ranges
            if(.not.allocated(vphv))then
              allocate(vphv(ni))
            else
              if(size(vphv).lt.ni)then
                deallocate(vphv)
                allocate(vphv(ni))
              end if
            end if
            ! Lower ranges
            if(.not.allocated(vplv))then
              allocate(vplv(ni))
            else
              if(size(vplv).lt.ni)then
                deallocate(vplv)
                allocate(vplv(ni))
              end if
            end if
            ! Extended upper range
            if(.not.allocated(vphve))then
              allocate(vphve(ni*2))
            else
              if(size(vphve).lt.ni*2)then
                deallocate(vphve)
                allocate(vphve(ni*2))
              end if
            end if
            ! Extended lower range
            if(.not.allocated(vplve))then
              allocate(vplve(ni*2))
            else
              if(size(vplve).lt.ni*2)then
                deallocate(vplve)
                allocate(vplve(ni*2))
              end if
            end if
            ! Resonances
            if(.not.allocated(vpr))then
              allocate(vpr(nr))
            else
              if(size(vpr).lt.nr)then
                deallocate(vpr)
                allocate(vpr(nr))
              end if
            end if

            !
            ! Store the frequencies of the FS transitions
            !

            ! Reset index
            nt = 0

            ! Output transition
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJf=1,Atom(ia)%nJ(itermf)
                if(abs(Atom(ia)%rJval(iJu,itermu) - &
                       Atom(ia)%rJval(iJf,itermf)).gt.1.or. &
                   Atom(ia)%rJval(iJu,itermu) + &
                   Atom(ia)%rJval(iJf,itermf).lt..4d0) cycle
                nt = nt + 1
                nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJf,itermf)
              end do
            end do

            ! Input transition
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJl=1,Atom(ia)%nJ(iterml)
                if(abs(Atom(ia)%rJval(iJu,itermu) - &
                       Atom(ia)%rJval(iJl,iterml)).gt.1.or. &
                   Atom(ia)%rJval(iJu,itermu) + &
                   Atom(ia)%rJval(iJl,iterml).lt..4d0)cycle
                nt = nt + 1
                nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJl,iterml)
              end do
            end do

            ! For each output direction
            do jdir=1,Frec%ndir

              if (dyn) then

                ! If line of sight
                if (LOS) then

                  ith1 = ithv(jdir)
                  iph1 = iphv(jdir)
                  ct = Geom%L_mu(ith1)
                  st = sqrt(1d0 - ct*ct)
                  cc = cos(Geom%L_phi(iph1))
                  sc = sin(Geom%L_phi(iph1))

                ! If quadrature
                else

                  ith1 = ithv(jdir)
                  iph1 = iphv(jdir)
                  ct = Geom%V_mu(ith1)
                  st = sqrt(1d0 - ct*ct)
                  cc = Geom%v_mux(iph1)
                  sc = Geom%v_muy(iph1)*sqrt(1d0 - cc*cc)

                end if
              end if

              ! For each height node
              do iz=Rz0,Rz1

                ! Get index
                indx = Frec%indx(jtran,ia,iz,jdir)
#ifdef _OPENMP
                ! Multi-thread
                if (omp) then
                  ! Skip if not assigned
                  if (indx.lt.oif0(tid)) cycle
                  if (indx.gt.oif1(tid)) exit
                end if
#endif
                ! Calculate Doppler shift factor
                if (dyn) &
                  vfac = 1d0 - atmo%vx(iz)*st*cc - &
                               atmo%vy(iz)*st*sc - &
                               atmo%vz(iz)*ct

                ! Store limits
                if0 = Atom(ia)%if0(jtran)
                if1 = Atom(ia)%if1(jtran)


                !
                ! Calculate the Doppler width of both input and output
                ! transitions
                !

                ! Thermal common part
                DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

                ! Input
                Dw1 = Atom(ia)%Dfreq(itran)*sqrt(DwT*DwT + &
                                             Atmo%vmi(iz)**2d0)
                ! Output
                Dw  = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                             Atmo%vmi(iz)**2d0)

                ! Transform the searching parameters from normalized
                ! to proper frequency units
                red_cohwW = Input%dcohw*Dw
                red_resoW = Input%red_pars(2)*Dw
                red_neglW = Input%red_pars(3)*Dw
                red_coreW = Input%red_pars(7)*Dw
                red_rangwW1 = Input%red_pars(1)*Dw1
                red_vlarwW1 = Input%red_pars(4)*Dw1
                red_fstpwW1 = Input%red_pars(5)*Dw1
                red_mstpwW1 = red_fstpwW1*Input%red_pars(6)
                red_rangcW1 = Input%red_pars(8)*Dw1
                red_vlarcW1 = Input%red_pars(9)*Dw1
                red_fstpcW1 = Input%red_pars(10)*Dw1
                red_mstpcW1 = red_fstpcW1*Input%red_pars(11)


                !
                ! Find limits for the second order output
                !
                if (.not.allocated(Frec%dzao(indx)%if0)) then

                  nran = 0
                  bf0 = -1
                  bf1 = -2

                  ! Only one output
                  if (if0.eq.if1) then

                    ! Check if we are close to a transition frequency
                    do it=1,ntj
                      if(abs(Frec%omega(if0)*vfac-nut(it)).lt. &
                         red_neglW)then
                        bf0 = if0
                        bf1 = if1
                        nran = 1
                        exit
                      end if
                    end do

                  ! More than one output
                  else

                    ! Reset logical variable
                    lskip = .True.

                    ! Look for the limits, for each output frequency
                    do ifreq=if0,if1

                      ! Initialize the flag
                      skip = .True.

                      ! Check if we are close to a transition
                      ! frequency
                      do it=1,ntj
                        if(abs(Frec%omega(ifreq)*vfac-nut(it)).lt. &
                           red_neglW)then
                          skip = .False.
                          exit
                        end if
                      end do

                      ! If we skip this frequency
                      if(skip.or.ifreq.eq.if1) then

                        if(.not.lskip) then
                          if (ifreq.eq.if1) then
                            bf1 = ifreq
                          else
                            bf1 = ifreq - 1
                          end if
                          nran = nran + 1
                        end if


                      ! If we cannot skip this frequency
                      else

                        ! If this is the first index for this
                        ! transition
                        if (bf0.lt.0) bf0 = ifreq

                      end if

                      lskip = skip

                    end do ! output frequencies

                  end if ! Number of outputs

                  ! Store in the array
                  Frec%dzao(indx)%nran = nran

                  if (Frec%dzao(indx)%nran.lt.1) cycle

                  ! Allocate the ranges
                  allocate(Frec%dzao(indx)%if0(nran))
                  Frec%dzao(indx)%if0 = bf0
                  allocate(Frec%dzao(indx)%if1(nran))
                  Frec%dzao(indx)%if1 = bf1

                  ! Only one output
                  if (if0.eq.if1) then

                    nran = 1
                    np = 1
                    Frec%dzao(indx)%if0(nran) = if0
                    Frec%dzao(indx)%if1(nran) = if0

                  ! More than one output
                  else

                    ! Only one output
                    if (bf1.eq.bf0) then

                      ! Initialize the flag
                      skip = .True.

                      ! Check if we are close to a transition
                      ! frequency
                      do it=1,ntj
                        if(abs(Frec%omega(ifreq)*vfac-nut(it)).lt. &
                           red_neglW)then
                          skip = .False.
                          exit
                        end if
                      end do

                      ! No transition
                      if (skip) then

                        Frec%dzao(indx)%nran = 0
                        deallocate(Frec%dzao(indx)%if0)
                        deallocate(Frec%dzao(indx)%if1)
                        cycle

                      ! Yes transition
                      else

                        Frec%dzao(indx)%if0(1) = bf0
                        Frec%dzao(indx)%if1(1) = bf0
                        np = 1
                        nran = 1

                      end if

                    ! Several outputs
                    else

                      ! Reset logical variable
                      lskip = .True.

                      nran = 0
                      np = 0

                      ! Look for the limits, for each output frequency
                      do ifreq=bf0,bf1

                        ! Initialize the flag
                        skip = .True.

                        ! Check if we are close to a transition
                        ! frequency
                        do it=1,ntj
                          if(abs(Frec%omega(ifreq)*vfac-nut(it)).lt. &
                             red_neglW)then
                            skip = .False.
                            exit
                          end if
                        end do

                        ! If we skip this frequency
                        if(skip.or.ifreq.eq.bf1) then

                          if(.not.lskip) then
                            if (ifreq.eq.bf1) then
                              Frec%dzao(indx)%if1(nran) = ifreq
                              np = np + 1
                            else
                              Frec%dzao(indx)%if1(nran) = ifreq - 1
                            end if
                          end if

                        ! If we cannot skip this frequency
                        else

                          np = np + 1

                          if (lskip) then
                            nran = nran + 1
                            Frec%dzao(indx)%if0(nran) = ifreq
                          end if
                        end if

                        lskip = skip

                      end do ! output frequencies

                    end if ! Number of outputs
                  end if ! Number of outputs

                  ! Set global limits
                  Frec%dzao(indx)%gf0 = minval(Frec%dzao(indx)%if0)
                  Frec%dzao(indx)%gf1 = maxval(Frec%dzao(indx)%if1)
                  Frec%dzao(indx)%ggf0 = 10000000
                  Frec%dzao(indx)%ggf1 = -1

                  ! Count frequencies
                  Frec%dzao(indx)%nfreq = np

                end if ! Already found limits

                ! Initialize sizes
                Frec%dzao(indx)%trani(iti)%osize = 0
                Frec%dzao(indx)%trani(iti)%isize = 0

                if (Frec%dzao(indx)%nran.lt.1) cycle

                ! Allocate input frequency size
                np = Frec%dzao(indx)%nfreq
                allocate(Frec%dzao(indx)%trani(iti)%mfreq(np))

                ! Initialize pointers and back dimension trace
                allocate(bomega, bw_freq)
                nullify(bomega%next)
                nullify(bomega%prev)
                nullify(bw_freq%next)
                nullify(bw_freq%prev)
                bomega%nback = 0
                bw_freq%nback = 0

                ! For each output frequency
                iifreq = 0
                do iran=1,Frec%dzao(indx)%nran
                  do ifreq=Frec%dzao(indx)%if0(iran), &
                           Frec%dzao(indx)%if1(iran)

                    iifreq = iifreq + 1

      !
      ! Reset indentation
      !


      ! Reset number freq
      Frec%dzao(indx)%trani(iti)%mfreq(iifreq) = 0

      !
      ! Find the limits specified by transitions
      !

      ! Coherent wings?
      if (Input%cohw) then

        ! Initialize flag
        cohw = .True.

        ! Check if close to any transition
        do it=1,ntj
          if (abs(nut(it) - Frec%omega(ifreq)*vfac).lt.red_cohwW) then
            cohw=.False.
            exit
          end if
        end do

      else

        cohw = .False.

      end if

      ! If Coherent wing
      if (cohw) then

        ! Definitely not a core frequency (if its wing!)
        core = .False.

        !
        ! Store the determined vector
        !

        ! Advance boxes
        if (.not.allocated(bomega%A)) then

          allocate(bomega%A(1))
          allocate(bw_freq%A(1))

        else

          ! bomega
          bdaux => bomega
          allocate(bomega%next)
          bomega => bomega%next
          bomega%prev => bdaux
          allocate(bomega%A(1))
          nullify(bomega%next)
          nullify(bdaux)
          bomega%nback = bomega%prev%nback + bomega%prev%mfreq

          ! bw_freq
          bdaux => bw_freq
          allocate(bw_freq%next)
          bw_freq => bw_freq%next
          bw_Freq%prev => bdaux
          allocate(bw_freq%A(1))
          nullify(bw_freq%next)
          nullify(bdaux)

        end if

        ! Store the fequency axis
        bomega%A = Frec%omega(ifreq)*vfac + dnlmin
        ! Store the dimension of the axis
        bomega%mfreq = 0
        Frec%dzao(indx)%trani(iti)%mfreq(iifreq) = 0
        ! Store indexes
        bomega%ifreq = iifreq

      ! Non-coherent wing
      else

        ! Reset logical
        core = .False.

        ! The output transition does not matter here, if the
        ! input is close to these, that means that they are also
        ! the input
        do it=1,ntj
          if (abs(nut(it) - Frec%omega(ifreq)*vfac).le.red_coreW) &
            core=.True.
          vphv(it) = -1D0
          vplv(it) = -1D0
        end do

        if (core) then
          red_rangW1 = red_rangcW1
          red_vlarW1 = red_vlarcW1
          red_fstpW1 = red_fstpcW1
          red_mstpW1 = red_mstpcW1
        else
          red_rangW1 = red_rangwW1
          red_vlarW1 = red_vlarwW1
          red_fstpW1 = red_fstpwW1
          red_mstpW1 = red_mstpwW1
        end if

        ! Take the ones in the input
        do it=ntj+1,nt
          vphv(it) = nut(it) + red_rangW1
          vplv(it) = nut(it) - red_rangW1
        end do


        !
        ! Find the limits specified by resonances
        !
        do ir=1,nr
          vpr(ir) = Frec%omega(ifreq)*vfac + dnl(ir)
          vphv(nt + ir) = vpr(ir) + red_rangW1
          vplv(nt + ir) = vpr(ir) - red_rangW1
        end do


        !
        ! Define the true limits:
        ! The resonances specify the true limits, but if
        ! they are close to a transition, we expand the limit
        !

        ! Take the furthest resonances
        vph = Frec%omega(ifreq)*vfac + dnlmax
        vpl = Frec%omega(ifreq)*vfac + dnlmin

        ! Check if we are close to an input transition
        do it=ntj+1,nt

          ! If we are, move the limit to the transition
          ! instead of the resonance
          if(abs(vph - nut(it)).lt.red_resoW.or. &
             abs(vpl - nut(it)).lt.red_resoW)then

            if(nut(it).lt.vpl)vpl=nut(it)
            if(nut(it).gt.vph)vph=nut(it)

          end if

        end do

        ! Now add the range from the parameters
        vph = vph + red_rangW1
        vpl = vpl - red_rangW1


        !
        ! Flag lines and resonances out of limits
        !

        ! The total number of resonances that we had
        ni = nr + nt

        ! We are going to change ni, so store it because
        ! we need the original value
        np = ni

        ! For each resonance and transition
        do ir=1,np

          ! Check if it is out of range
          if(vphv(ir).lt.vpl.or.vplv(ir).gt.vph)then

            vplv(ir) = vpl - 1
            ni = ni - 1

          ! It is not out of range, but the lower limit is
          ! out
          else if(vplv(ir).lt.vpl)then
            vplv(ir) = vpl

          ! It is not out of range, but the upper limit is
          ! out
          else if(vphv(ir).gt.vph)then
            vphv(ir) = vph
          end if

        end do


        !
        ! Check lines and resonances that spawns the same
        ! range
        !

        ! For each pair of transitions and resonances
        do ir=np,2,-1

          ! This is equivalent to be flagged out
          if(vplv(ir).lt.vpl)cycle

          ! The other index to make a pair
          do it=ir-1,1,-1

            ! This is equivalent to be flagged out
            if(vplv(it).lt.vpl)cycle

            ! If the ranges overlap, combine them into
            ! just one and flag the other out
            if((vphv(ir).ge.vplv(it).and.vphv(ir).le.vphv(it)).or. &
               (vplv(ir).ge.vplv(it).and.vplv(ir).le.vphv(it)))then
              vplv(it) = min(vplv(ir),vplv(it))
              vphv(it) = max(vphv(ir),vphv(it))
              vplv(ir) = vpl - 1
              ni = ni - 1
              exit
            end if

          end do ! it
        end do ! ir


        !
        ! Shift the individual limits moving the valid ones
        ! to the first part of the vector
        !

        ! If we have changed the number of ranges from the
        ! beginning
        if(ni.ne.np)then

          ! For each pair of resonances
          do ir=1,np-1

            ! If it is flagged out
            if(vplv(ir).lt.vpl)then

              ! Reset the added index
              ip = 1

              ! For all the resonances in front of this one
              do it=ir+1,np

                ! If it is flagged, skip it
                if(vplv(it).lt.vpl)then

                  ip = ip + 1

                ! If it is not flagged, move it to the
                ! position of the flagged one
                else

                  vplv(it-ip) = vplv(it)
                  vphv(it-ip) = vphv(it)
                  vplv(it) = vpl - 1
                  exit

                end if

              end do ! it
            end if ! ir flagged
          end do ! ir
        end if ! if something flagged


        !
        ! Order the individual limits
        !

        ! Lower limits
        call QsortC(vplv(1:ni))
        ! Upper limits
        call QsortC(vphv(1:ni))


        !
        ! Define extended limits
        !

        ! Reset the index of extended limits
        nie = 0

        ! For each normal limit, get an extended version
        do ir=1,ni
          nie = nie + 1
          vplve(nie) = vplv(ir) - red_vlarW1
          vphve(nie) = vplv(ir) - red_fstpW1
          nie = nie + 1
          vplve(nie) = vphv(ir) + red_fstpW1
          vphve(nie) = vphv(ir) + red_vlarW1
        end do

        !
        ! Check lines and resonances that spawns the same range
        !

        ! Store the number of original limits
        np = nie

        ! For each pair of transitions and resonances
        do ir=np,2,-1

          ! This is equivalent to be flagged out
          if(vplve(ir).lt.0d0)cycle

          ! The other index to make a pair
          do it=ir-1,1,-1

            ! If the ranges overlap, combine them into
            ! just one and flag the other out
            if(vplve(it).lt.0d0)cycle
            if((vphve(ir).ge.vplve(it).and. &
                vphve(ir).le.vphve(it)).or. &
               (vplve(ir).ge.vplve(it).and. &
                vplve(ir).le.vphve(it)))then
              vplve(it) = min(vplve(ir),vplve(it))
              vphve(it) = max(vphve(ir),vphve(it))
              vplve(ir) = -1d0
              nie = nie - 1
              exit
            end if
          end do ! it
        end do ! ir


        !
        ! Shift the individual limits moving the valid ones
        ! to the first part of the vector
        !

        ! If we have changed the number of ranges from the
        ! beginning
        if(nie.ne.np)then

          ! For each pair of resonances
          do ir=1,np-1

            ! If it is flagged out
            if(vplve(ir).lt.0d0)then

              ! Reset the added index
              ip = 1

              ! For all the resonances in front of this one
              do it=ir+1,np

                ! If it is flagged, skip it
                if(vplve(it).lt.0d0)then

                  ip = ip + 1

                ! If it is not flagged, move it to the
                ! position of the flagged one
                else

                  vplve(it-ip) = vplve(it)
                  vphve(it-ip) = vphve(it)
                  vplve(it) = -1d0
                  exit

                end if
              end do ! it
            end if ! ir flagged
          end do ! ir
        end if ! if something flagged


        !
        ! Order the individual limits
        !

        ! Lower limits
        call QsortC(vplve(1:nie))
        ! Upper limits
        call QsortC(vphve(1:nie))

        !
        ! Build the vector of input frequencies from the
        ! limits
        !

        ! Flag to reset in case we run out of space while
        ! building a vector
        reset = .False.

        ! Do until we are finished
        do while (.True.)

          ! If the flag is one, we need more space for the
          ! vector
          if(reset)then
            np0 = np0*2
            deallocate(vpp)
            allocate(vpp(np0))
            deallocate(flag)
            allocate(flag(np0))
            reset = .False.
          end if

          ! Reset index counter
          np = 0


          !
          ! Build for the short limits
          !

          ! For each short range
          do it=1,ni

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if(np.gt.np0)then
              reset = .True.
              exit
            end if

            ! Start with the lower limit
            vpp(np) = vplv(it)

            ! Do until finished
            do while(.True.)

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Next frequency is the previous plus the step
              vpp(np) = vpp(np-1) + red_fstpW1

              ! If we are over the range, we are done
              if(vpp(np).ge.vphv(it))then
                vpp(np) = vphv(it)
                exit
              end if

            end do

            ! If we have no space, we need to reset
            if (reset) exit

          end do

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          !
          ! Build for the extended limits
          !

          ! For each long range
          do it=1,nie

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if(np.gt.np0)then
              reset = .True.
              exit
            end if

            ! We start with the lower limit
            vpp(np) = vplve(it)

            ! Do until finished
            do while(.True.)

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Next frequency is the previous plus the step
              vpp(np) = vpp(np-1) + red_mstpW1

              ! If we are over the range, we are done
              if(vpp(np).ge.vphve(it))then
                vpp(np) = vphve(it)
                exit
              end if

            end do

            ! If we have no space, we need to allocate it
            ! above
            if (reset) exit

          end do

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          !
          ! Add the resonance frequencies to the vector
          !

          ! The number of frequencies we already have
          npp = np

          ! For each resonance
          do ir=1,nr

            ! Reset the flag
            nfound = .True.

            ! Run over the existing frequencies
            do ip=1,npp

              ! If the frequency is there, do not add it
              if (abs(1d2/vpp(ip) - 1d2/vpr(ir)).lt.resolin) then
                nfound = .False.
                exit
              end if

            end do

            ! If we did not find it
            if (nfound) then

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Add the frequency
              vpp(np) = vpr(ir)

            end if

            ! If we have no space, we need to allocate it
            ! above
            if (reset) exit

          end do

          if(reset)cycle

          exit ! If we get to this point, we have everything

        end do ! First do while


        !
        ! Check for duplicates
        !

        ! Reset the flag
        flag(1:np) = 1

        ! For each frequency
        do ip=1,np

          ! If it has been flagged, we already checked
          if (flag(ip).lt.1) cycle

          ! Check the following ones
          do ipp=ip+1,np

            ! If it has been flagged, we already checked
            if (flag(ipp).lt.1) cycle

            ! If some of them are repeated, flag them to be
            ! removed
            if(abs(1d2/vpp(ip)-1d2/vpp(ipp)).lt.resolin) &
              flag(ipp) = 0

          end do ! ipp
        end do ! ip

        ! Reset the running real index
        ipp = 0

        ! For each frequency in the vector
        do ip=1,np

          ! If it is flagged correct, add to real vector
          if(flag(ip).gt..5)then
            ipp = ipp + 1
            vpp(ipp) = vpp(ip)
          end if

        end do

        ! The number of frequencies is the last value of ipp
        np = ipp

        ! Order the frequency axis
        call QsortC(vpp(1:np))

        !
        ! Store the determined vector
        !

        ! Advance boxes
        if (.not.allocated(bomega%A)) then

          allocate(bomega%A(np))
          allocate(bw_freq%A(np))

        else

          ! bomega
          bdaux => bomega
          allocate(bomega%next)
          bomega => bomega%next
          bomega%prev => bdaux
          allocate(bomega%A(np))
          nullify(bomega%next)
          nullify(bdaux)
          bomega%nback = bomega%prev%nback + bomega%prev%mfreq

          ! bw_freq
          bdaux => bw_freq
          allocate(bw_freq%next)
          bw_freq => bw_freq%next
          bw_Freq%prev => bdaux
          allocate(bw_freq%A(np))
          nullify(bw_freq%next)
          nullify(bdaux)

        end if

        ! Check that we have enough space to work below
        if(.not.allocated(Wvpp))then
          allocate(Wvpp(np))
        else
          if(size(Wvpp).lt.np)then
            deallocate(Wvpp)
            allocate(Wvpp(np))
          end if
        end if

        ! Store the fequency axis
        bomega%A = vpp(1:np)
        ! Store the dimension of the axis
        bomega%mfreq = np
        Frec%dzao(indx)%trani(iti)%mfreq(iifreq) = np
        ! Store indexes
        bomega%ifreq = iifreq

        ! Update the maximum of input frequencies
        if (np.gt.Frec%dzao(indx)%mxfreq) Frec%dzao(indx)%mxfreq = np


        !
        ! Define the integration weights (same than
        ! omegabuild)
        !

        ! The first point is special in compound trapezoidal
        ! rule
        Wvpp(1) = .5d0*(vpp(2) - vpp(1))

        ! Initialize the integral to normalize the weights
        norm1 = Wvpp(1)

        ! The initial lower limit is the first point
        O0 = vpp(1)

        ! This is the pointer to the first element of the
        ! current interval, we are pointing to the first
        ! element
        cfreq = 1

        ! Flag that says that the point 2 is not the initial
        ! point of the interval (because 1 is the initial
        ! point)
        init = .FALSE.

        ! For the rest of frequencies except the last
        do jfreq=2,np-1

          ! If ifreq is the initial point of an interval
          if(init)then

            ! The first point is special in compound
            ! trapezoidal rule
            Wvpp(jfreq) = .5d0*(vpp(jfreq+1) - vpp(jfreq))

            ! The next point cannot be a first point
            init = .FALSE.

            ! Initialize the integral to normalize the
            ! weights
            norm1 = Wvpp(jfreq)

            ! Pointer is now in this frequency
            cfreq = jfreq

            ! And it is the beginning of the current interval
            O0 = vpp(jfreq)

          ! If jfreq is not the initial point of an interval
          else

            ! Check if ifreq is the last point of an interval
            if(abs(vpp(jfreq+1) - vpp(jfreq)).gt.red_neglW)then

              ! The last point is special in compound
              ! trapezoidal rule
              Wvpp(jfreq) = .5d0*(vpp(jfreq) - vpp(jfreq-1))

              ! It is the end of the current interval
              O1 = vpp(jfreq)

              ! Add to the integral
              norm1 = norm1 + Wvpp(jfreq)

              ! We know that the integral must be
              ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
              norm = 1d5*(O1 - O0)/norm1

              ! Normalize the weights of this interval
              do kfreq=cfreq,jfreq
                Wvpp(kfreq) = Wvpp(kfreq)*norm
              end do

              ! The next point is the first point of its
              ! interval
              init = .TRUE.

            ! If jfreq is not the last point of an interval
            else

              ! Compound trapezoidal rule weight
              Wvpp(jfreq) = .5d0*(vpp(jfreq+1) - vpp(jfreq-1))

              ! Add to the integral
              norm1 = norm1 + Wvpp(jfreq)

            endif ! Last point

          end if ! Initial point

        end do ! jfreq

        ! The last point is special in compound trapezoidal
        ! rule
        Wvpp(np) = .5d0*(vpp(np) - vpp(np-1))

        ! It is the end of the interval
        O1 = vpp(np)

        ! Add to the integral
        norm1 = norm1 + Wvpp(np)

        ! We know that the integral must be
        ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
        norm = 1d5*(O1 - O0)/norm1

        ! Normalize the weights of this interval
        do jfreq=cfreq,np
          Wvpp(jfreq) = Wvpp(jfreq)*norm
        end do

        ! Store the weights
        bw_freq%A = Wvpp(1:np)

      end if ! Coherent wing

                  end do ! output frequencies
                end do ! output frequency ranges

                !
                ! Properly store and index the data
                !

                ! Total dimension of omega and w_freq
                nn = bomega%nback + bomega%mfreq

                ! Allocate omega and W_freq
                allocate(Frec%dzao(indx)%trani(iti)%omega(nn))
                allocate(Frec%dzao(indx)%trani(iti)%w_freq(nn))

                ! Determine size
                Frec%dzao(indx)%trani(iti)%osize = nn

                ! Go backwards in the linked lists
                do while (.True.)

                  iifreq = bomega%ifreq
                  ip = bomega%nback + 1
                  ipp = ip + bomega%mfreq - 1
                  if (ipp.ge.ip) then
                    Frec%dzao(indx)%trani(iti)% &
                                    omega(ip:ipp) = bomega%A
                    Frec%dzao(indx)%trani(iti)% &
                                    W_freq(ip:ipp)= bw_freq%A
                  end if

                  ! Deallocate arrays
                  deallocate(bomega%A,bw_freq%A)

                  ! If last one, clean and quit
                  if (.not.associated(bomega%prev)) then
                    deallocate(bomega,bw_freq)
                    nullify(bomega,bw_freq)
                    exit
                  ! Not done with the list
                  else
                    bomega => bomega%prev
                    bw_freq => bw_freq%prev
                    nullify(bomega%next%prev,bw_freq%next%prev)
                    deallocate(bomega%next,bw_freq%next)
                    nullify(bomega%next,bw_freq%next)
                  end if

                end do ! Run backwards the frequency axes

                ! Update RAM
!$omp flush(MPID)
                MPID%RAM = MPID%RAM + 16d-6*dble(nn)
                MPID%WRAM = MPID%WRAM + 16d-6*dble(nn)
!$omp flush(MPID)

                ! End of array constructions

              end do ! Height
            end do ! Output directions
          end do ! Input transition
        end do ! Output transition

!$omp barrier

        !
        ! Allocate space for interpolation and define it or
        ! find the index limits
        !

        ! For each output direction
        do jdir=1,Frec%ndir

          !
          ! If coherent wings, compute vfac
          if (Input%cohw) then

            ! Actually dynamic
            if (dyn) then

              ! If line of sight
              if (LOS) then

                ith1 = ithv(jdir)
                iph1 = iphv(jdir)
                ct = Geom%L_mu(ith1)
                st = sqrt(1d0 - ct*ct)
                cc = cos(Geom%L_phi(iph1))
                sc = sin(Geom%L_phi(iph1))

              ! If quadrature
              else

                ith1 = ithv(jdir)
                iph1 = iphv(jdir)
                ct = Geom%V_mu(ith1)
                st = sqrt(1d0 - ct*ct)
                cc = Geom%v_mux(iph1)
                sc = Geom%v_muy(iph1)*sqrt(1d0 - cc*cc)

              end if

            else

              vfac = 1d0

            end if
          end if

          ! For each height
          do iz=Rz0,Rz1

            !
            ! If coherent wings and dynamic, compute vfac
            if (Input%cohw.and.dyn) &
              vfac = 1d0 - atmo%vx(iz)*st*cc - &
                           atmo%vy(iz)*st*sc - &
                           atmo%vz(iz)*ct

            ! For each output transition level
            do jtran=1,Atom(ia)%ntran

              ! If no PRD line, skip
              if (.not.Atom(ia)%lemiss2(jtran)) cycle
              if (Atom(ia)%fflag(jtran)%absent) cycle

              ! Look for the terms
              itermu = -1
              do i=1,Atom(ia)%nMulti-1
                do i1=i+1,Atom(ia)%nMulti
                  if (Atom(ia)%irad(i,i1).eq.jtran) then
                    itermf = i
                    itermu = i1
                    exit
                  end if
                end do
                if (itermu.ge.0) exit
              end do

              ! Get index
              indx = Frec%indx(jtran,ia,iz,jdir)
#ifdef _OPENMP
              ! Multi-thread
              if (omp) then
                ! Skip if not assigned
                if (indx.lt.oif0(tid)) cycle
                if (indx.gt.oif1(tid)) exit
              end if
#endif
              if (Frec%dzao(indx)%nran.lt.1) cycle

              ! structure index
              ktran = Atom(ia)%itrano(jtran)

              ! For all the possible lower terms
              do i=1,itermu-1

                ! If there is no transition or this term is larger
                ! than the upper term of the output transition, skip
                if(i.ge.itermu.or.Atom(ia)%irad(i,itermu).eq.0) cycle

                ! Store the input lower term index
                iterml = i

                ! Get index of input transition
                itran = Atom(ia)%irad(iterml,itermu)

                ! Get index of input transition in structure
                ios = -1
                do iti=1,Atom(ia)%trano(ktran)%nt
                  if (Atom(ia)%trano(ktran)%ind(iti).eq.itran) then
                    ios = 1
                    exit
                  end if
                end do
                if (ios.lt.0) cycle

        !
        ! Reset identation
        !

        ! Point to input frequency
        p_frec => Frec%dzao(indx)%trani(iti)

        ! Predict size of next block
        nn = sum(p_frec%mfreq)

        ! If angle-dependent
        if (.not.AV) then

          ! If dynamic extra dimensions, if static just frequencies
          if (dyn) then

            ! For axial problems
            if (axial) then

              ! Size is just polar
              nn = nn*Geom%nTh

            ! For non-axial problems
            else

              ! Skip backward rayleigh
              nn = nn*(Geom%nTh*Geom%nPh2 - Frec%nfs(jdir))

            end if ! Axial
          end if ! Dynamic
        end if ! AD

        ! Predict aditional frequency
        SRAM = 16d-6*dble(nn)

        ! If can store
        if (TPRAM) then
!$omp flush(MPID,ofram)
          ! If no more space
          if (floor(MPID%RAM+SRAM).gt.RLIM.or.SRAM.le.0d0) then
            ofram = .True.
            p_frec%RAM = .False.
          else
            allocate(p_frec%index1(nn))
            allocate(p_frec%index2(nn))
            allocate(p_frec%dx(nn))
            MPID%RAM = MPID%RAM + SRAM
            MPID%WRAM = MPID%WRAM + SRAM
            p_frec%RAM = .True.
          end if
!$omp flush(MPID,ofram)
        ! Cannot store
        else
          p_frec%RAM = .False.
        end if

        !
        ! If coherent wings, recover resonance limits
        if (Input%cohw) then

          ! Make sure that we have enough space to store resonances
          if(.not.allocated(dnl))then
            allocate(dnl(1))
          end if

          ! Reset ranges
          dnlmax = -1D99
          dnlmin = 1D99

          ! For each pair of lower levels for input and output
          do iJl=1,Atom(ia)%nJ(iterml)
            do iJf=1,Atom(ia)%nJ(itermf)

              ! Add the resonance
              dnl(1) = Atom(ia)%FSfreq(iJf,itermf) - &
                       Atom(ia)%FSfreq(iJl,iterml)

              ! And check that we know what are the futhest ones
              if(dnl(1).gt.dnlmax)dnlmax = dnl(1)
              if(dnl(1).lt.dnlmin)dnlmin = dnl(1)

            end do ! iJl
          end do ! iJf

        end if ! Coherent wings


        !
        ! Define interpolation

        ! Initialize index
        jjfreq0 = 0
        kkfreq0 = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Frec%dzao(indx)%nran
          do ifreq=Frec%dzao(indx)%if0(iran), &
                   Frec%dzao(indx)%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Input frequency number
            np = p_frec%mfreq(iifreq)

            ! For each input direction
            do ith=1,Frec%nth
              do iph=1,Frec%nph

                ! If dynamics and AD
                if (dyn.and..not.AV) then

                  ! For axial problems
                  if (axial) then

                    ! Automatically skip extra azimuths
                    if (iph.gt.1) cycle

                    ! Get director cosines
                    ct1 = Geom%V_mu(ith)

                    ! Calculate Doppler shift factor
                    vfac1 = 1d0 - atmo%vz(iz)*ct1

                    ! We will be using the inverse
                    vfac1 = 1d0/vfac1

                  ! For non-axial problems
                  else

                    ! If angle-dependent, check backward Rayleigh
                    ! scattering
                    if (jtran.eq.itran.and. &
                        Frec%stype(iph,ith,jdir).lt.0) cycle

                    ! Get director cosines
                    ct1 = Geom%V_mu(ith)
                    st1 = sqrt(1d0 - ct1*ct1)
                    cc1 = Geom%v_mux(iph)
                    sc1 = Geom%v_muy(iph)*sqrt(1d0 - cc1*cc1)

                    ! Calculate Doppler shift factor
                    vfac1 = 1d0 - atmo%vx(iz)*st1*cc1 - &
                                  atmo%vy(iz)*st1*sc1 - &
                                  atmo%vz(iz)*ct1

                    ! We will be using the inverse
                    vfac1 = 1d0/vfac1

                  end if ! Axial

                ! Not dynamic or AV
                else

                  ! No shift
                  vfac1 = 1d0

                  ! Only one direction
                  if (iph.gt.1.or.ith.gt.1) cycle

                end if ! Dynamics

                ! Reset indexes
                jjfreq = jjfreq0
                kkfreq = kkfreq0

                ! Skip empty
                if (np.lt.1) then

                  ! Left and right limits from resonance
                  lifreq = 1
                  jfreq = nfreq
                  O0 = Frec%omega(ifreq)*vfac + dnlmin
                  O1 = Frec%omega(ifreq)*vfac + dnlmax

                  !
                  ! Look for the indexes

                  ! Left
                  do while (.True.)

                    ! Check if already inside
                    if (Frec%omega(lifreq)*vfac1.ge.O0) exit

                    lifreq = lifreq + 1
                    if (lifreq.gt.nfreq) exit

                  end do

                  ! Right
                  do while (.True.)

                    ! Check if already inside
                    if (Frec%omega(jfreq)*vfac1.le.O1) exit

                    jfreq = jfreq - 1
                    if (jfreq.lt.1) exit

                  end do


                  ! Update global limits
                  if (lifreq.lt.Frec%dzao(indx)%ggf0) &
                    Frec%dzao(indx)%ggf0 = lifreq
                  if (jfreq.gt.Frec%dzao(indx)%ggf1) &
                    Frec%dzao(indx)%ggf1 = jfreq

                  ! Skip rest
                  cycle

                end if

                ! Add to size
                Frec%dzao(indx)%trani(iti)%isize = np + &
                                      Frec%dzao(indx)%trani(iti)%isize

                ! If storing
                if (p_frec%RAM) then

                  !
                  ! Reset identation
                  !

      ! Reset the search frequency
      lifreq = 1

      ! For each input frequency
      do jfreq=1,np

        ! Advance indexes
        jjfreq = jjfreq + 1
        kkfreq = kkfreq + 1

        ! If out of range, take the value at the
        ! boundary
        if (p_frec%omega(jjfreq)*vfac1.le.Frec%omega(1)+TINYO) then

          ! We are still looking in the first one
          lifreq = 1

          ! The index to take is 1
          p_frec%index1(kkfreq) = 1

          ! The index to take is 1
          p_frec%index2(kkfreq) = 1

          ! We do not need this number
          p_frec%dx(kkfreq) = 0d0

        ! If out of range, take the value at the boundary
        else if (p_frec%omega(jjfreq)*vfac1.ge. &
                 (Frec%omega(nfreq) - TINYO)) then

          ! We are in the last frequency
          lifreq = nfreq

          ! The index to take is nfreq
          p_frec%index1(kkfreq) = nfreq

          ! The index to take is nfreq
          p_frec%index2(kkfreq) = nfreq

          ! We do not need this number
          p_frec%dx(kkfreq) = 0d0

        ! If within the boundaries
        else

          ! Search between the last found frequency and
          ! all but the boundary
          do ibfreq=lifreq,nfreq-1

            ! If this exact frequency is in output
            if (abs(p_frec%omega(jjfreq)*vfac1 - &
                    Frec%omega(ibfreq)).lt.TINYO) then

              ! We are in the found frequency
              lifreq = ibfreq

              ! This frequency gives us the value
              p_frec%index1(kkfreq) = lifreq

              ! This frequency gives us the value
              p_frec%index2(kkfreq) = lifreq

              ! We do not need this number
              p_frec%dx(kkfreq) = 0d0

              exit

            ! If the input is between this output and
            ! the next
            else if(p_frec%omega(jjfreq)*vfac1.ge. &
                    Frec%omega(ibfreq).and. &
                    p_frec%omega(jjfreq)*vfac1.lt. &
                    Frec%omega(ibfreq+1)) then

              ! We found it in the index of the lower
              lifreq = ibfreq

              ! The first index is the lower
              p_frec%index1(kkfreq) = lifreq

              ! The second index is the upper
              p_frec%index2(kkfreq) = lifreq+1

              ! Store the inverse of the distance between
              ! the two outputs
              p_frec%dx(kkfreq) = &
                  (p_frec%omega(jjfreq)*vfac1 - Frec%omega(lifreq))/ &
                  (Frec%omega(lifreq+1) - Frec%omega(lifreq))

              exit

            end if ! Check output frequency

          end do ! Run output frequencies

        end if ! Check if out of limits

        ! Update global limits
        if (p_frec%index1(kkfreq).lt.Frec%dzao(indx)%ggf0) &
          Frec%dzao(indx)%ggf0 = p_frec%index1(kkfreq)
        if (p_frec%index2(kkfreq).gt.Frec%dzao(indx)%ggf1) &
          Frec%dzao(indx)%ggf1 = p_frec%index2(kkfreq)

      end do ! Run input frequencies

                  !
                  ! Restore identation
                  !

                ! Not storing
                else

                  !
                  ! Reset identation
                  !

       ! Compute jump
       lifreq = 1
       jufreq = np
       if (np.gt.1) jufreq = jufreq - 1

       ! First and last frequencies
       do jfreq=jjfreq+1,jjfreq+np,jufreq

         ! If out of range, take the value at the
         ! boundary
         if (p_frec%omega(jfreq)*vfac1.le.Frec%omega(1)+TINYO) then

           if (Frec%dzao(indx)%ggf0.gt.1) Frec%dzao(indx)%ggf0 = 1
           if (Frec%dzao(indx)%ggf1.lt.1) Frec%dzao(indx)%ggf1 = 1

         ! If out of range, take the value at the boundary
         else if (p_frec%omega(jfreq)*vfac1.ge. &
                  (Frec%omega(nfreq) - TINYO)) then

           ! We are in the last frequency
           if (Frec%dzao(indx)%ggf0.gt.nfreq) &
             Frec%dzao(indx)%ggf0 = nfreq
           if (Frec%dzao(indx)%ggf1.lt.nfreq) &
             Frec%dzao(indx)%ggf1 = nfreq

         ! If within the boundaries
         else

           ! Search between the last found frequency and
           ! all but the boundary
           do ibfreq=lifreq,nfreq-1

             ! If this exact frequency is in output
             if (abs(p_frec%omega(jfreq)*vfac1 - &
                     Frec%omega(ibfreq)).lt.TINYO) then

               ! Found frequency
               if (Frec%dzao(indx)%ggf0.gt.ibfreq) &
                 Frec%dzao(indx)%ggf0 = ibfreq
               if (Frec%dzao(indx)%ggf1.lt.ibfreq) &
                 Frec%dzao(indx)%ggf1 = ibfreq

               exit

             ! If the input is between this output and the next
             else if(p_frec%omega(jfreq)*vfac1.ge. &
                                            Frec%omega(ibfreq).and. &
                     p_frec%omega(jfreq)*vfac1.lt. &
                                            Frec%omega(ibfreq+1)) then

               ! Found frequency
               if (Frec%dzao(indx)%ggf0.gt.ibfreq) &
                 Frec%dzao(indx)%ggf0 = ibfreq
               if (Frec%dzao(indx)%ggf1.lt.ibfreq+1) &
                 Frec%dzao(indx)%ggf1 = ibfreq+1

               exit

             end if ! Check output frequency

           end do ! Run output frequencies

         end if ! Check if out of limits

       end do ! Run input frequencies

       ! Fake the advance of frequencies
       jjfreq = jjfreq + np
       kkfreq = kkfreq + np


                  !
                  ! Restore indentation
                  !

                end if ! Storing

                if (.not.AV.and..not.axial.and.dyn) &
                  kkfreq0 = kkfreq

              end do ! Input azimuth

              ! Update kkfreq
              if (.not.AV.and.dyn) kkfreq0 = kkfreq

            end do ! Input polar

            ! Update index
            jjfreq0 = jjfreq
            kkfreq0 = kkfreq

            end do ! Output frequencies
          end do ! Output frequency ranges

          ! Nullify local pointer
          nullify(p_frec)

                  !
                  ! Restore identation
                  !

              end do ! Lower intput term
            end do ! Output transition
          end do ! height nodes
        end do ! output directions

!$omp barrier

        !
        ! Build indexes for atom and get sizes
        !

!$omp single

        ! If there is at least one node with magnetic field and
        ! not already defined
        if (Yfield.and..not.allocated(Atom(ia)%i_Wind)) then

          ! Allocate internal index term range
          allocate(Atom(ia)%i_Wind(mint:maxt))

          ! For each term in the range where we need indexes
          do itermu=mint,maxt

            nMm = nint(2d0*(Atom(ia)%rLval(itermu) + &
                           Atom(ia)%Sval(itermu)) + 1d0)
            nblk = 0

            ! For each mu
            do i=1,nMm
              if (Atom(ia)%nblk(i,itermu).gt.nblk) &
                nblk = Atom(ia)%nblk(i,itermu)
            end do

            ! Allocate index size
            allocate(Atom(ia)%i_Wind(itermu)%ind(nblk,nMm))

            ! Index the components
            i = 0
            do i1=1,nMm
              do i2=1,Atom(ia)%nblk(i1,itermu)
                i = i+1
                Atom(ia)%i_Wind(itermu)%ind(i2,i1) = i
              end do
            end do

          end do ! For each term in the range of second order

        end if ! There is one node with magnetic field


        !
        ! Build the indexes for Warr2
        !

        ! If non-coherent lower term and we are storing
        ! redistribution with magnetic field
        if (Yfield.and.NCHLT.and.PRAM) then

          ! Check if for this atom there is any NCHLT
          ! height
          lNCHLT = .False.
          do iz=Rz0,Rz1
            if (Atom(ia)%NCHLT(iz,itran)) then
              lNCHLT = .True.
              exit
            end if
          end do

        else

          lNCHLT = .False.

        end if
!$omp end single

        ! For each upper term
        do itermu=2,Atom(ia)%nMulti

          ! For each final lower term
          do itermf=1,itermu-1

            jtran = Atom(ia)%irad(itermu,itermf)

            if (jtran.le.0) cycle
            if (.not.Atom(ia)%lemiss2(jtran)) cycle
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Look for the ktran
            ktran = Atom(ia)%itrano(jtran)

            rJumax = Atom(ia)%rLval(itermu) + &
                     Atom(ia)%Sval(itermu)
            nMu = nint(2d0*rJumax + 1d0)
            rJfmax = Atom(ia)%rLval(itermf) + &
                     Atom(ia)%Sval(itermf)
            nMf = nint(2d0*rJfmax + 1d0)

            ! For each other lower term
            do iterml=1,itermu-1

              itran = Atom(ia)%irad(itermu,iterml)

              if (itran.le.0) cycle

              ! Find the transition index
              ios = -1
              do iti=1,Atom(ia)%trano(ktran)%nt

                if (Atom(ia)%trano(ktran)%ind(iti).eq.itran) then
                  ios = 1
                  exit
                end if

              end do ! Input transitions
              if (ios.lt.0) cycle

              rJlmax = Atom(ia)%rLval(iterml) + &
                       Atom(ia)%Sval(iterml)
              nMl = nint(2d0*rJlmax + 1d0)
!$omp barrier
!$omp single
              ! If at least one height with magnetic feld
              if (Yfield) then

                do i1=1,2

                  if (i1.eq.1) then
                    MindU = 0
                    MindU1 = 0
                    MindL = 0
                    MindL1 = 0
                    MindF = 0
                  else

                    ! Allocate input transition part
                    allocate(Atom(ia)%trano(ktran)%trani(iti)% &
                      Wind(MindL1,MindL,MindF,MindU1,MindU))
                    ! Initialize to zero
                    Atom(ia)%trano(ktran)%trani(iti)%Wind = 0

                  end if

                  ! Reset indexing
                  i = 0

                  ! For each Mf
                  do iMf=1,nMf

                    ! Value of Mf
                    rMf = -rJfmax + dble(iMf-1)

                    ! For each mu_f
                    do mF=1,Atom(ia)%nblk(iMf,itermf)

                      indF = Atom(ia)%i_Wind(itermf)%ind(mF,iMf)

                      ! For each Mu
                      do iMu=1,nMu

                        ! Value of Mu
                        rMu = -rJumax + dble(iMu-1)

                        ! Difference between M momentums, done integer
                        q = rMu - rMf
                        iq = nint(q)

                        ! If not pi nor sigma, skip
                        if(abs(iq).gt.1) cycle

                        ! For each mu_u
                        do iU=1,Atom(ia)%nblk(iMu,itermu)

                          indU = Atom(ia)%i_Wind(itermu)%ind(iU,iMu)

                          ! For each Mu'
                          do iMu1=1,nMu

                            ! Value of Mu'
                            rMu1 = -rJumax + dble(iMu1-1)

                            ! Difference between M momentums
                            q1 = rMu1-rMf
                            QQ = q1-q

                            ! Convert to integers
                            iq1 = nint(q1)
                            iQQ = nint(QQ)

                            ! If not pi or sigma, skip
                            if(abs(iq1).gt.1) cycle

                            ! For each mu_u'
                            do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                              indU1 = Atom(ia)%i_Wind(itermu)% &
                                               ind(iU1,iMu1)

      !
      ! Reset indexing
      !

      ! For each Ml
      do iMl=1,nMl

        ! Value of Ml
        rMl = -rJlmax + dble(iMl-1)

        ! Difference between M momentums
        p = rMu-rMl
        ip = nint(p)

        ! If not pi nor sigma, skip
        if(abs(ip).gt.1) cycle

        ! For each mu_l
        do iL=1,Atom(ia)%nblk(iMl,iterml)

          indL = Atom(ia)%i_Wind(iterml)%ind(iL,iMl)

          ! For each Ml'
          do iMl1=1,nMl

            ! Value of Ml'
            rMl1 = -rJlmax + dble(iMl1-1)

            ! Difference between M momentums
            p1 = rMu1-rMl1
            PP = p1-p

            ! Convert to integer
            ip1 = nint(p1)
            iPP = nint(PP)

            ! If not pi nor sigma, skip
            if(abs(ip1).gt.1) cycle

            ! For each mu_l'
            do iL1=1,Atom(ia)%nblk(iMl1,iterml)

              indL1 = Atom(ia)%i_Wind(iterml)%ind(iL1,iMl1)

              ! First round
              if (i1.eq.1) then

                if (indU.gt.MindU) MindU = indU
                if (indU1.gt.MindU1) MindU1 = indU1
                if (indL.gt.MindL) MindL = indL
                if (indL1.gt.MindL1) MindL1 = indL1
                if (indF.gt.MindF) MindF = indF

              ! Second round
              else

                i = i+1
                Atom(ia)%trano(ktran)%trani(iti)% &
                  Wind(indL1,indL,indF,indU1,indU) = i

              end if

            end do ! iL1
          end do ! iMl1
        end do ! iL
      end do ! iMl
                          !
                          ! Recover indexing
                          !
                            end do ! iU1
                          end do ! iMu1
                        end do ! iU
                      end do ! iMu
                    end do ! mF
                  end do ! iMf
                end do ! i1 (round of counting)

                iYNF = i

                ! If non-coherent lower term and we are storing
                ! redistribution
                if (lNCHLT) then

                  ! Allocate
                  if (allocated(WNCHLT)) deallocate(WNCHLT)
                  allocate(WNCHLT(iYNF))
                  WNCHLT = .True.

                  ! Count elements
                  iYYF = 0

                  ! For each Mf
                  do iMf=1,nMf

                    ! Value of Mf
                    rMf = -rJfmax + dble(iMf-1)

                    ! For each mu_f
                    do mF=1,Atom(ia)%nblk(iMf,itermf)

                      indF = Atom(ia)%i_Wind(itermf)%ind(mF,iMf)

                      ! For each Mu
                      do iMu=1,nMu

                        ! Value of Mu
                        rMu = -rJumax + dble(iMu-1)

                        ! If not pi nor sigma, skip
                        if(abs(nint(rMu-rMf)).gt.1) cycle

                        ! For each mu_u
                        do iU=1,Atom(ia)%nblk(iMu,itermu)

                          indU = Atom(ia)%i_Wind(itermu)%ind(iU,iMu)

                          ! For each Mu'
                          do iMu1=1,nMu

                            ! Value of Mu'
                            rMu1 = -rJumax + dble(iMu1-1)

                            ! If not pi or sigma, skip
                            if(abs(nint(rMu1-rMf)).gt.1) cycle

                            ! For each mu_u'
                            do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                              indU1 = Atom(ia)%i_Wind(itermu)% &
                                               ind(iU1,iMu1)

                              ! For each Ml
                              do iMl=1,nMl

                                ! Value of Ml
                                rMl = -rJlmax + dble(iMl-1)

                                ! If not pi nor sigma, skip
                                if(abs(nint(rMu-rMl)).gt.1) cycle

                                ! For each mu_l
                                do iL=1,Atom(ia)%nblk(iMl,iterml)

                                  indL = Atom(ia)%i_Wind(iterml)% &
                                                  ind(iL,iMl)

                                  ! For each Ml'
                                  do iMl1=1,nMl

                                    ! NCHLT
                                    if (iMl.ne.iMl1) cycle

                                    ! Value of Ml'
                                    rMl1 = -rJlmax + dble(iMl1-1)

                                    ! If not pi nor sigma, skip
                                    if(abs(nint(rMu1-rMl1)).gt.1) &
                                      cycle

                                    ! For each mu_l'
                                    do iL1=1,Atom(ia)% &
                                                     nblk(iMl1,iterml)

                                      indL1 = Atom(ia)% &
                                              i_Wind(iterml)% &
                                              ind(iL1,iMl1)

                                      if (iL.ne.iL1) cycle

                                      iYYF = iYYF+1
                                      i1 = Atom(ia)%trano(ktran)% &
                                                    trani(iti)% &
                                                    Wind(indL1,indL,&
                                                         indF,indU1,&
                                                         indU)
                                      WNCHLT(i1) = .False.

                                    end do ! iL1
                                  end do ! iMl1
                                end do ! iL
                              end do ! iMl
                            end do ! iU1
                          end do ! iMu1
                        end do ! iU
                      end do ! iMu
                    end do ! mF
                  end do ! iMf

                end if ! NCHLT and storing redistribution
              end if ! Yes field

              ! If there is at least one point without magnetic field
              if (Nfield) then

                do i1=1,2

                  if (i1.eq.1) then
                    MiindU = 1000000
                    MiindU1 = 1000000
                    MiindL = 1000000
                    MiindL1 = 1000000
                    MiindF = 1000000
                    MindU = 0
                    MindU1 = 0
                    MindL = 0
                    MindL1 = 0
                    MindF = 0
                  else

                    ! Allocate input transition part
                    allocate(Atom(ia)%trano(ktran)%trani(iti)% &
                             WindNB(MiindL1:MindL1,MiindL:MindL, &
                                    MiindF:MindF,MiindU1:MindU1, &
                                    MiindU:MindU))

                  end if

                  ! Reset indexing
                  i = 0

                  ! For each Jf
                  do mF=1,Atom(ia)%nJ(itermf)

                    ! Get indexes
                    indF = Atom(ia)%irho(itermf)%irho_ij(mF)

                    ! Get Jf
                    rJf = Atom(ia)%rJval(mF,itermf)

                    ! For each Ju
                    do iU=1,Atom(ia)%nJ(itermu)

                      ! Get indexes
                      indU = Atom(ia)%irho(itermu)%irho_ij(iU)

                      ! Get Ju
                      rJu = Atom(ia)%rJval(iU,itermu)

                      if (abs(rJu-rJf).gt.1d0.or. &
                          rJu+rJf.lt..25) cycle

                      ! For each Ju'
                      do iU1=1,Atom(ia)%nJ(itermu)

                        ! Get indexes
                        indU1 = Atom(ia)%irho(itermu)%irho_ij(iU1)

                        ! Get Ju'
                        rJu1 = Atom(ia)%rJval(iU1,itermu)

                        if (abs(rJu1-rJf).gt.1d0.or. &
                            rJu1+rJf.lt..25) cycle

                        ! For each Jl
                        do iL=1,Atom(ia)%nJ(iterml)

                          ! Get indexes
                          indL = Atom(ia)%irho(iterml)%irho_ij(iL)

                          ! Get Jl
                          rJl = Atom(ia)%rJval(iL,iterml)

                          if (abs(rJu-rJl).gt.1d0.or. &
                              rJu+rJl.lt..25) cycle

                          ! For each Jl'
                          do iL1=1,Atom(ia)%nJ(iterml)

                            ! Get indexes
                            indL1 = Atom(ia)%irho(iterml)%irho_ij(iL1)

                            ! Get Jl1
                            rJl1 = Atom(ia)%rJval(iL1,iterml)

                            if (abs(rJu1-rJl1).gt.1d0.or. &
                                rJu1+rJl1.lt..25) cycle

                            ! First round
                            if (i1.eq.1) then

                              if (indU.lt.MiindU) MiindU = indU
                              if (indU1.lt.MiindU1) MiindU1 = indU1
                              if (indL.lt.MiindL) MiindL = indL
                              if (indL1.lt.MiindL1) MiindL1 = indL1
                              if (indF.lt.MiindF) MiindF = indF
                              if (indU.gt.MindU) MindU = indU
                              if (indU1.gt.MindU1) MindU1 = indU1
                              if (indL.gt.MindL) MindL = indL
                              if (indL1.gt.MindL1) MindL1 = indL1
                              if (indF.gt.MindF) MindF = indF

                            ! Second round
                            else

                              i = i+1
                              Atom(ia)%trano(ktran)%trani(iti)% &
                                WindNB(indL1,indL,indF,indU1,indU) = i

                            end if
                          end do ! Jl'
                        end do ! Jl
                      end do ! Ju'
                    end do ! Ju
                  end do ! Jf
                end do ! Round of counting
                iNF = i
              end if ! No B field

!$omp end single

              !
              ! Allocate space for Warr2
              !
              ! If Storing in RAM
              if (PRAM) then

                ! For each output direction in the quadrature
                do jdir=1,Red%ndir

                  jbdir = min(jdir,Frec%ndir)
                  jcdir = min(jdir,Red%njdir)

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Get index
                    indxf = Frec%indx(jtran,ia,iz,jbdir)
                    indx = Red%indx(jtran,ia,iz,jdir)
#ifdef _OPENMP
                    ! If multi-thread
                    if (omp) then
                      ! Skip if not assigned
                      if (indx.lt.orif0(tid)) cycle
                      if (indx.gt.orif1(tid)) exit
                    end if
#endif
                    ! Check not already allocated
                    if (allocated(Red%dzao(indx)% &
                                      trani(iti)%PWarr2)) cycle

                    if (Bstrength(iz).lt.TINYB) then

                      iDF = iNF
                      iDFR = iNF

                    else

                      iDF = iYNF

                      ! Check if NCHLT
                      if (lNCHLT) then

                        ! If NCHLT applies
                        if (Atom(ia)%NCHLT(iz,itran)) then
                          iDFR = iYYF
                        else
                          iDFR = iYNF
                        end if

                      else
                        iDFR = iYNF
                      end if

                    end if ! Magnetic field or not

                    ! Predict size of next block
                    nn = sum(Frec%dzao(indxf)%trani(iti)%mfreq)*iDFR
                    if (jtran.eq.itran) then
                      nn = nn*(Red%nth*Red%nph - Frec%nfs(jcdir))
                    else
                      nn = nn*Red%nth*Red%nph
                    end if
                    SRAM = 8d-6*dble(nn)

                    ! If no more space
!$omp flush(MPID,ofram)
                    if (floor(MPID%RAM+SRAM).gt.RLIM.or. &
                      SRAM.le.0d0) then
                      ofram = .True.
!$omp flush(ofram)
                      Red%dzao(indx)%trani(iti)%RAM = .False.
                      cycle
                    end if
                    MPID%RAM = MPID%RAM + SRAM
                    MPID%WRAM = MPID%WRAM + SRAM
!$omp flush(MPID)

                    Red%dzao(indx)%trani(iti)%RAM = .True.

                    ! Allocate initialize flag
                    allocate(Red%dzao(indx)%trani(iti)%iPPRD(iDF))
                    Red%dzao(indx)%trani(iti)%iPPRD = .True.

                    ! Allocate Warr2
                    allocate(Red%dzao(indx)%trani(iti)%PWarr2(nn))


                  end do ! heights
                end do ! output directions

              ! If not storing in RAM
              else

                ! Allocate dummy array
!$omp single
                if (.not.associated(Red%dzao)) then
                  allocate(Red%dzao(1))
                  nullify(Red%dzao(1)%trani)
                end if
!$omp end single

              end if ! Storing in RAM

            end do ! Input lower term
          end do ! Output lower term
        end do ! Upper term
      end do ! Atom
!$omp end parallel

#ifdef _OPENMP
      ! deallocate limits for threads
      if (omp) deallocate(oif0,oif1)
#endif

      ! Check if everything is fine
      call control

      return

      end subroutine omegabuildin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Build the input frequency axis for the multilevel intensity
      !! problem.\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Input(Input_class): Structure with settings data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!              lio(logical): If doing formal solution in this
      !!                            run\n
      !!            ofram(logical): Indicates if out of RAM\n
      !!              LOS(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine omegabuildinI(Frec,Red,Atom,Atmo,Input,Geom, &
                               MPID,lio,ofram,LOS)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Atmo_class):: Atmo
      type(Input_class):: Input
      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class):: MPID
      logical, intent(in):: lio,LOS
      logical, intent(out):: ofram

      ! Local

      logical:: skip,lskip,init,reset,nfound,core,cohw

      integer:: ios,i2,it,iz,ia,ip,ipp,ir,itran,jtran,iran,nran
      integer:: fitran,fjtran,ffitran,ffjtran,ffktran,ibfreq
      integer:: ifreq,jfreq,iifreq,lifreq,kfreq,cfreq,if0,if1
      integer:: jdir,jbdir,jcdir,njdir,ith,iph,ith1,iph1
      integer:: iJl,iJu,iJf,itermf,itermu,iterml,iti,bf0,bf1
      integer:: ni,nie,np,np0,npp,nt,nti,nat,mina,maxa
      integer:: jjfreq0,jjfreq,kkfreq0,kkfreq,jufreq
      integer:: minto,maxto,indx,indxf,nn
#ifdef _OPENMP
      integer:: tid
#endif
      integer, dimension(:), allocatable:: minti,maxti
      integer, dimension(:), allocatable:: flag,ithV,iphV
#ifdef _OPENMP
      integer, dimension(:), allocatable:: oif0,oif1,orif0,orif1,onf
#endif

      double precision:: DwT,Dw,Dw1,vph,vpl,vfac,vfac1
      double precision:: norm,norm1,O0,O1
      double precision:: red_rangW1,red_resoW
      double precision:: red_vlarW1,red_fstpW1
      double precision:: red_mstpW1,red_neglW
      double precision:: red_rangwW1,red_vlarwW1
      double precision:: red_fstpwW1,red_mstpwW1
      double precision:: red_rangcW1,red_vlarcW1
      double precision:: red_fstpcW1,red_mstpcW1
      double precision:: red_coreW,red_cohwW
      double precision:: dnl, nut, nutout, vpr
      double precision:: ct,st,cc,sc,ThK,ct1,st1,cc1,sc1,SRAM
      double precision, dimension(2):: vphv, vplv
      double precision, dimension(4):: vphve, vplve
      double precision, dimension(:), allocatable:: vpp
      double precision, dimension(:), allocatable:: Wvpp

      ! Box
      type(dbabox_class), pointer:: bomega, bw_freq, bdaux

      ! Pointer
      type(Frequencyd_class), pointer:: p_frec

      ! Initialize
      nullify(bomega,bw_freq,bdaux)

      ! Initialize
      ofram = .False.

      ! Check if we don't have to reset
      if (AVI.and..not.dyn.and.LOS) then
        if (lio) return
      end if

      ! If line of sight
      if ((.not.AVI.or.dyn).and.LOS) then
        PRAM = .False.
        IRAM = .False.
      end if

      ! Inititlize RAM counters
      MPID%WRAM = 0d0

      ! Deallocate existing Frec and Red
      call cleanFrecandRed(Frec,Red,MPID)

      ! If the process is the master, it does not need this, but
      ! it will participate if storing in a file
      if (MPID%mpi.and.pid.eq.0) then

        ! And leave
        call control
        return

      end if

      ! If angle dependent
      if (.not.AVI) then

        !
        ! Output and input directions

        ! If for emergence
        if (LOS) then

          njdir = Geom%nPhLOS*Geom%nThLOS

        ! If quadrature
        else

          njdir = Geom%nPh*Geom%nTh

        end if ! LOS

        ! If dynamic, axis changes
        if (dyn) then

          Frec%ndir = njdir
          Frec%nth = Geom%nTh

          ! Azimuth depends on symmetry
          if (axiali) then
            Frec%nph = 1
          else
            Frec%nph = Geom%nPh
          end if

        ! If static, they don't
        else

          Frec%ndir = 1
          Frec%nth = 1
          Frec%nph = 1

        end if

        Red%ndir = Geom%nPh2*Geom%nTh
        Red%njdir = njdir
        Red%nth = Geom%nTh
        Red%nph = Geom%nPh2

        ! Allocate auxiliar quantities for Doppler shifts

        ! Index of polar direction of quadrature
        allocate(ithv(njdir))
        ! Index of azimuthal direction of quadrature
        allocate(iphv(njdir))

        ! Allocate type of scattering
        if (allocated(Frec%stype)) deallocate(Frec%stype)
        if (allocated(Frec%nfs)) deallocate(Frec%nfs)
        allocate(Frec%stype(Geom%nPh2,Geom%nTh,njdir))
        allocate(Frec%nfs(njdir))

        ! Initialize to normal scattering
        Frec%stype = 0
        Frec%nfs = 0

        !
        ! De-index the directions

        jdir = 0

        ! If lines of sight
        if (LOS) then

          do ith=1,Geom%nThLOS
            do iph=1,Geom%nPhLOS

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

              ! If angle-dependent
              if (.not.AVI) then

                do ith1=1,Geom%nTh
                  do iph1=1,Geom%nPh2

                    ! Calculate scattering angle between the
                    ! quadrature direction and the LOS direction
                    ThK = atom2lab(Geom%L_theta(ith), &
                                   Geom%L_phi(iph), &
                                   Geom%V_theta(ith1), &
                                   Geom%V_phi(iph1))

                    ! Backward condition
                    if (abs(pi - ThK).le.TINYA) then

                      Frec%stype(iph1,ith1,jdir) = 1

                    ! Forward condition
                    else if (ThK.le.TINYA) then

                      Frec%stype(iph1,ith1,jdir) = -1
                      Frec%nfs(jdir) = Frec%nfs(jdir) + 1

                    end if ! Backward or forward scattering

                  end do ! Quadrature directions (input)
                end do

              end if ! AD redistribution

            end do ! LOS directions
          end do

        ! If quadrature
        else

          ! Quadrature
          do ith=1,Geom%nTh
            do iph=1,Geom%nPh

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

              ! If angle-dependent
              if (.not.AVI) then

                do ith1=1,Geom%nTh
                  do iph1=1,Geom%nPh2

                    ! Backward condition
                    if ((Geom%nTh - ith + 1).eq.ith1.and. &
                        int(iph + .5d0*Geom%V_muy(iph)* &
                            Geom%nPh2).eq.iph1) then

                      Frec%stype(iph1,ith1,jdir) = 1

                    ! Forward condition
                    else if (ith1.eq.ith.and. &
                             iph1.eq.iph) then

                      Frec%stype(iph1,ith1,jdir) = -1
                      Frec%nfs(jdir) = Frec%nfs(jdir) + 1

                    end if ! Backward or forward scattering

                  end do ! Quadrature directions (input)
                end do

              end if ! AD redistribution

            end do
          end do ! Quadrature directions

        end if ! LOS

        ! If static, initialize velocity factors
        if (.not.dyn) then

          vfac = 1d0
          vfac1 = 1d0

        end if

      ! If dynamic (and AA)
      else if (dyn) then

        !
        ! Output and input directions

        ! If lines of sight
        if (LOS) then
          njdir = Geom%nPhLOS*Geom%nThLOS
        ! If quadrature
        else
          njdir = Geom%nPh*Geom%nTh
        end if

        Frec%ndir = njdir
        Frec%nth = 1
        Frec%nph = 1

        Red%ndir = Geom%nPh*Geom%nTh
        Red%njdir = 1
        Red%nth = 1
        Red%nph = 1

        ! Allocate auxiliar quantities for Doppler shifts

        ! Index of polar direction of quadrature
        allocate(ithv(njdir))
        ! Index of azimuthal direction of quadrature
        allocate(iphv(njdir))

        ! Unneccessary type of scattering
        if (allocated(Frec%stype)) deallocate(Frec%stype)
        if (allocated(Frec%nfs)) deallocate(Frec%nfs)
        allocate(Frec%stype(1,1,1))
        allocate(Frec%nfs(1))

        ! Initialize to normal scattering
        Frec%stype = 0
        Frec%nfs = 0

        !
        ! De-index the directions

        jdir = 0

        ! If lines of sight
        if (LOS) then

          ! LOS
          do ith=1,Geom%nThLOS
            do iph=1,Geom%nPhLOS

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

            end do ! LOS directions
          end do

        ! If quadrature
        else

          ! Quadrature
          do ith=1,Geom%nTh
            do iph=1,Geom%nPh

              jdir = jdir + 1
              ithv(jdir) = ith
              iphv(jdir) = iph

            end do ! Quadrature directions
          end do

        end if ! LOS

      ! Angle-averaged and static
      else

        ! Output directions
        Frec%ndir = 1
        Red%ndir = 1

        ! Input directions
        Frec%nth = 1
        Frec%nph = 1
        Red%nth = 1
        Red%nph = 1
        Red%njdir = 1

        ! Doppler factors
        vfac = 1d0
        vfac1 = 1d0

        ! Allocate stype if PRD to avoid undefined
        if (PRD) then

          ! Allocate type of scattering
          if (allocated(Frec%stype)) deallocate(Frec%stype)
          if (allocated(Frec%nfs)) deallocate(Frec%nfs)
          allocate(Frec%stype(1,1,1))
          allocate(Frec%nfs(1))

          ! Initialize to normal scattering
          Frec%stype = 0
          Frec%nfs = 0

        end if ! PRD

      end if

      !
      ! Count maximum index of PRD atom and transition

      ! Initialize atomic and transition indexes, and counter
      ! of real elements
      mina = 10000
      maxa = 0
      minto = 10000
      maxto = 0
      nat = 0

      ! For each atom
      do ia=1,nA
        ! For all transitions
        do jtran=1,Atom(ia)%ntran

          ! If PRD line
          if (Atom(ia)%lemiss2(jtran)) then
            ! If not absent in this CPU
            if (.not.Atom(ia)%fflag(jtran)%absent) then

              ! Go by fine structure
              do it=1,Atom(ia)%fst(jtran)%nt

                ! Get index
                ffjtran = Atom(ia)%ifst_ij(it,jtran)

                ! Add to counter
                nat = nat + 1

                ! Update limits
                if (ffjtran.lt.minto) minto = ffjtran
                if (ffjtran.gt.maxto) maxto = ffjtran
                if (ia.lt.mina) mina = ia
                if (ia.gt.maxa) maxa = ia

              end do ! Fine structure

            end if ! Presence of line
          end if ! PRD line

        end do ! Transitions
      end do ! Atoms

      ! Allocate indexing array and first step of Frec and Red
      allocate(Frec%indx(minto:maxto,mina:maxa,Rz0:Rz1,Frec%ndir))
      Frec%ndzao = nat*Rnz*Frec%ndir
      allocate(Frec%dzao(Frec%ndzao))
      do indx=1,Frec%ndzao
        nullify(Frec%dzao(indx)%trani)
      end do
      if (IRAM) then
        allocate(Red%indx(minto:maxto,mina:maxa,Rz0:Rz1,Red%ndir))
        Red%ndzao = nat*Rnz*Red%ndir
        allocate(Red%dzao(Red%ndzao))
        do indx=1,Red%ndzao
          nullify(Red%dzao(indx)%trani)
        end do
      end if

      !
      ! Build index

      ! Initialize
      ip = 0

      ! Directions
      do jdir=1,Frec%ndir
        ! Height
        do iz=Rz0,Rz1
          ! Atom
          do ia=1,nA
            ! Transition
            do jtran=1,Atom(ia)%ntran

              ! If PRD line
              if (Atom(ia)%lemiss2(jtran)) then
                ! If not absent in this CPU
                if (.not.Atom(ia)%fflag(jtran)%absent) then

                  ! Fine structure
                  do it=1,Atom(ia)%fst(jtran)%nt

                    ffjtran = Atom(ia)%ifst_ij(it,jtran)
                    ip = ip + 1
                    Frec%indx(ffjtran,ia,iz,jdir) = ip

                  end do ! Fine structure

                end if ! Presence
              end if ! PRD line

            end do ! Transition
          end do ! Atom
        end do ! Height
      end do ! Directions

#ifdef _OPENMP
      ! If multiple threads
      if (omp) then

        ! Allocate limits for threads
        allocate(oif0(nthread),oif1(nthread),onf(nthread))

        ! Work per thread
        ios = Frec%ndzao/nthread

        ! Give this first stimation to each process
        do tid=1,nthread
          onf(tid) = ios
        end do

        ! Put the remaining heights in some of the threads if there
        ! are remaining ones
        if(ios*nthread.ne.Frec%ndzao)then

          ! Number of nodes to distribute
          ios = Frec%ndzao - ios*nthread

          ! Give them to the first threads
          do tid=1,ios
            onf(tid) = onf(tid) + 1
          end do

        end if

        ! Set first thread boundary
        oif0(1) = 1
        oif1(1) = onf(1)

        ! For each other thread
        do tid=2,nthread
          oif0(tid) = oif1(tid-1) + 1
          oif1(tid) = oif0(tid) + onf(tid) - 1
        end do

        ! Deallocate onf
        deallocate(onf)

      end if ! multithreaded
#endif
      ! If storing redistribution
      if (IRAM) then

        ! Initialize
        ip = 0

        ! Directions
        do jdir=1,Red%ndir
          ! Height
          do iz=Rz0,Rz1
            ! Atom
            do ia=1,nA
              ! Transition
              do jtran=1,Atom(ia)%ntran

                ! If PRD line
                if (Atom(ia)%lemiss2(jtran)) then
                  ! If not absent in this CPU
                  if (.not.Atom(ia)%fflag(jtran)%absent) then

                    ! Fine structure
                    do it=1,Atom(ia)%fst(jtran)%nt

                      ffjtran = Atom(ia)%ifst_ij(it,jtran)
                      ip = ip + 1
                      Red%indx(ffjtran,ia,iz,jdir) = ip

                    end do ! Fine structure

                  end if ! Presence
                end if ! PRD line

              end do ! Transition
            end do ! Atom
          end do ! Height
        end do ! Directions
#ifdef _OPENMP
        ! If threading
        if (omp) then

          ! Allocate limits for threads
          allocate(orif0(nthread),orif1(nthread),onf(nthread))

          ! Work per thread
          ios = Red%ndzao/nthread

          ! Give this first stimation to each process
          do tid=1,nthread
            onf(tid) = ios
          end do

          ! Put the remaining heights in some of the threads if there
          ! are remaining ones
          if(ios*nthread.ne.Red%ndzao)then

            ! Number of nodes to distribute
            ios = Red%ndzao - ios*nthread

            ! Give them to the first threads
            do tid=1,ios
              onf(tid) = onf(tid) + 1
            end do

          end if

          ! Set first thread boundary
          orif0(1) = 1
          orif1(1) = onf(1)

          ! For each other thread
          do tid=2,nthread
            orif0(tid) = orif1(tid-1) + 1
            orif1(tid) = orif0(tid) + onf(tid) - 1
          end do

          ! Deallocate onf
          deallocate(onf)

        end if ! Multi-thread
#endif
      end if ! IRAM

!$omp parallel default(none) &
!$omp private(ia,skip,jtran,minto,maxto,nt,itermu,iJu,itermf,iJf) &
!$omp private(fjtran,ffjtran,iterml,itran,iJl,fitran,ffitran,tid) &
!$omp private(ffktran,nti,iti,dnl,nutout,nut,jdir,iz,ith,jbdir) &
!$omp private(iph,ct,st,cc,sc,indx,vfac,vfac1,if0,if1,DwT,Dw1,Dw) &
!$omp private(red_resoW,red_neglW,red_coreW,red_rangwW1,red_vlarwW1) &
!$omp private(red_fstpwW1,red_mstpwW1,red_rangcW1,red_vlarcW1) &
!$omp private(red_fstpcW1,red_mstpcW1,nran,bf0,bf1,lskip,ifreq,np) &
!$omp private(bomega,bw_freq,iifreq,iran,core,red_rangW1,red_vlarW1) &
!$omp private(red_fstpW1,red_mstpW1,vphv,vplv,vpr,vph,vpl,ni,ir,O0) &
!$omp private(nie,vplve,vphve,it,reset,vpp,flag,np0,npp,nfound) &
!$omp private(ip,ipp,bdaux,Wvpp,norm1,cfreq,init,jfreq,O1) &
!$omp private(nn,p_frec,ios,SRAM,jjfreq0,kkfreq0,ct1,st1,cc1,sc1) &
!$omp private(jjfreq,kkfreq,lifreq,ibfreq,norm,jcdir,indxf) &
!$omp private(red_cohwW,cohw) &
!$omp shared(na,Atom,minti,maxti,Input,Frec,Red,nz,LOS,atmo,IRAM) &
!$omp shared(MPID,AVI,dyn,axiali,Geom,nfreq,jufreq,ithv,iphv,TIRAM) &
!$omp shared(oif0,oif1,orif0,orif1,RLIM,ofram,omp,Rz0,Rz1)

      ! Preliminar allocation of vpp, auxiliar vector to store
      ! frequencies
      np0 = 10000
      allocate(vpp(np0))
      allocate(flag(np0))

#ifdef _OPENMP
      vfac = 1d0
      vfac1 = 1d0
      tid = omp_get_thread_num() + 1
#endif

      ! For each atom
      do ia=1,nA

        ! Check that there is at least 1 PRD line in this atom
        skip = .True.

        do jtran=1,Atom(ia)%ntran
          if (Atom(ia)%lemiss2(jtran)) then
            if (.not.Atom(ia)%fflag(jtran)%absent) then
              skip = .False.
              exit
            end if
          end if
        end do

        ! There are no PRD lines for this atom
        if (skip) cycle

!$omp barrier
!$omp single
        ! Reallocate minti and maxti
        if (allocated(minti)) deallocate(minti,maxti)
        allocate(minti(Atom(ia)%nftran))
        allocate(maxti(Atom(ia)%nftran))

        !
        ! Count the transition combinations
        !

        ! Reset index
        minto = Atom(ia)%nftran + 1
        maxto = 0
        minti = Atom(ia)%nftran + 1
        maxti = 0
        nt = 0

        ! For each upper level
        do itermu=2,Atom(ia)%nMulti
          do iJu=1,Atom(ia)%nJ(itermu)

            ! For each final lower level
            do itermf=1,itermu-1
              do iJf=1,Atom(ia)%nJ(itermf)

                jtran = Atom(ia)%irad(itermu,itermf)

                if (jtran.le.0) cycle
                if (.not.Atom(ia)%lemiss2(jtran)) cycle
                if (Atom(ia)%fflag(jtran)%absent) cycle

                fjtran = Atom(ia)%fst(jtran)%irad(iJu,iJf)

                if (fjtran.le.0) cycle

                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                nt = nt + 1

                if (ffjtran.lt.minto) minto = ffjtran
                if (ffjtran.gt.maxto) maxto = ffjtran

                ! For each other lower level
                do iterml=1,itermu-1

                  itran = Atom(ia)%irad(itermu,iterml)

                  if (itran.le.0) cycle

                  if (.not.Input%Raman.and.itran.ne.jtran) cycle

                  do iJl=1,Atom(ia)%nJ(iterml)

                    fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

                    if (fitran.le.0) cycle

                    ffitran = Atom(ia)%ifst_ij(fitran,itran)

                    if (ffitran.lt.minti(ffjtran)) &
                      minti(ffjtran) = ffitran

                    if (ffitran.gt.maxti(ffjtran)) &
                      maxti(ffjtran) = ffitran

                  end do ! iJl
                end do ! iterml
              end do ! iJf
            end do ! itermf
          end do ! iJu
        end do ! itermu

        ! Store amount of output transitions
        Atom(ia)%ntrano = nt

        ! Index trano
        if (allocated(Atom(ia)%itrano)) deallocate(Atom(ia)%itrano)
        allocate(Atom(ia)%itrano(minto:maxto))
        nt = 0

        ! Index transitions

        ! For each upper level
        do itermu=2,Atom(ia)%nMulti
          do iJu=1,Atom(ia)%nJ(itermu)

            ! For each final lower level
            do itermf=1,itermu-1
              do iJf=1,Atom(ia)%nJ(itermf)

                jtran = Atom(ia)%irad(itermu,itermf)

                if (jtran.le.0) cycle
                if (.not.Atom(ia)%lemiss2(jtran)) cycle
                if (Atom(ia)%fflag(jtran)%absent) cycle

                fjtran = Atom(ia)%fst(jtran)%irad(iJu,iJf)

                if (fjtran.le.0) cycle

                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

                ! Store index
                nt = nt + 1
                Atom(ia)%itrano(ffjtran) = nt

              end do ! iJf
            end do ! itermf
          end do ! iJu
        end do ! itermu

        !
        ! Allocate needed transition structures
        !

        ! Atomic structure
        if (allocated(Atom(ia)%trano)) deallocate(Atom(ia)%trano)
        allocate(Atom(ia)%trano(Atom(ia)%ntrano))
!$omp end single

        ! For each upper level
        do itermu=2,Atom(ia)%nMulti
          do iJu=1,Atom(ia)%nJ(itermu)

            ! For each final lower level
            do itermf=1,itermu-1
              do iJf=1,Atom(ia)%nJ(itermf)

                jtran = Atom(ia)%irad(itermu,itermf)

                if (jtran.le.0) cycle
                if (.not.Atom(ia)%lemiss2(jtran)) cycle
                if (Atom(ia)%fflag(jtran)%absent) cycle

                fjtran = Atom(ia)%fst(jtran)%irad(iJu,iJf)

                if (fjtran.le.0) cycle

                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                ffktran = Atom(ia)%itrano(ffjtran)

                ! Number of input transitions
                nti = 0

                ! For each other lower level
                do iterml=1,itermu-1

                  itran = Atom(ia)%irad(itermu,iterml)

                  if (itran.le.0) cycle

                  ! Skip if no Raman and is not Rayleigh
                  if (.not.Input%Raman.and.itran.ne.jtran) cycle

                  do iJl=1,Atom(ia)%nJ(iterml)

                    fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

                    if (fitran.le.0) cycle

                    ! Advance the continuous index
                    nti = nti + 1

                  end do ! Lower level
                end do ! Lower term

                ! Allocate atomic structure
!$omp single
                allocate(Atom(ia)%trano(ffktran)%trani(nti))
                allocate(Atom(ia)%trano(ffktran)%ind(nti))
                Atom(ia)%trano(ffktran)%nt = nti

                ! Reset index
                Atom(ia)%trano(ffktran)%ind = -1
!$omp end single

                ! For each direction
                do jdir=1,Frec%ndir

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Get index
                    indx = Frec%indx(ffjtran,ia,iz,jdir)
#ifdef _OPENMP
                    ! If multi-thread
                    if (omp) then
                      ! Skip if not assigned
                      if (indx.lt.oif0(tid)) cycle
                      if (indx.gt.oif1(tid)) exit
                    end if
#endif
                    ! Allocate the input transition type
                    if(.not.associated(Frec%dzao(indx)%trani)) then

                      ! Allocate
                      allocate(Frec%dzao(indx)%trani(nti))

                      ! Nullify pointers
                      do i2=1,nti
                        nullify(Frec%dzao(indx)%trani(i2)%index1)
                        nullify(Frec%dzao(indx)%trani(i2)%index2)
                        nullify(Frec%dzao(indx)%trani(i2)%dx)
                      end do

                    end if ! No input transitions associated yet

                    ! Initialize the number of input frequencies
                    Frec%dzao(indx)%mxfreq = 0

                  end do ! heights
                end do ! output directions

                ! If storing Warr
                if (IRAM) then

                  ! For each direction
                  do jdir=1,Red%ndir

                    ! For each height
                    do iz=Rz0,Rz1

                      ! Get index
                      indx = Red%indx(ffjtran,ia,iz,jdir)
#ifdef _OPENMP
                      ! If multi-thread
                      if (omp) then
                        ! Skip if not assigned
                        if (indx.lt.orif0(tid)) cycle
                        if (indx.gt.orif1(tid)) exit
                      end if
#endif
                      ! Allocate the input transition type
                      if(.not.associated(Red%dzao(indx)%trani)) &
                        allocate(Red%dzao(indx)%trani(nti))

                    end do ! heights
                  end do ! output directions

                end if ! Storing Warr
!$omp single
                ! Reset index of input transition
                iti = 0

                ! For each other lower level
                do iterml=1,itermu-1

                  itran = Atom(ia)%irad(itermu,iterml)

                  if (itran.le.0) cycle

                  ! Skip if no Raman and is not Rayleigh
                  if (.not.Input%Raman.and.itran.ne.jtran) cycle

                  do iJl=1,Atom(ia)%nJ(iterml)

                    fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

                    if (fitran.le.0) cycle

                    ffitran = Atom(ia)%ifst_ij(fitran,itran)

                    ! Advance the continuous index
                    iti = iti + 1

                    ! Store index of this input transition
                    Atom(ia)%trano(ffktran)%ind(iti) = ffitran

                  end do ! Lower level
                end do ! Lower term
!$omp end single
                ! Reset index of input transition
                iti = 0


        !
        ! Reset identation
        !

        ! For each other lower level
        do iterml=1,itermu-1

          itran = Atom(ia)%irad(itermu,iterml)

          if (itran.le.0) cycle

          ! Skip if no Raman and is not Rayleigh
          if (.not.Input%Raman.and.itran.ne.jtran) cycle

          do iJl=1,Atom(ia)%nJ(iterml)

            fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

            if (fitran.le.0) cycle

            ffitran = Atom(ia)%ifst_ij(fitran,itran)

            ! Advance the continuous index
            iti = iti + 1

            !
            ! Store resonance
            !

            dnl = Atom(ia)%FSfreq(iJf,itermf) - &
                  Atom(ia)%FSfreq(iJl,iterml)


            !
            ! Store the frequencies of transitions
            !

            ! Output transition
            nutout = Atom(ia)%FSfreq(iJu,itermu) - &
                     Atom(ia)%FSfreq(iJf,itermf)

            ! Input transition
            nut = Atom(ia)%FSfreq(iJu,itermu) - &
                  Atom(ia)%FSfreq(iJl,iterml)

            ! For each output direction
            do jdir=1,Frec%ndir

              if (dyn) then

                ! If line of sight
                if (LOS) then

                  ith = ithv(jdir)
                  iph = iphv(jdir)
                  ct = Geom%L_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = cos(Geom%L_phi(iph))
                  sc = sin(Geom%L_phi(iph))

                ! If quadrature
                else

                  ith = ithv(jdir)
                  iph = iphv(jdir)
                  ct = Geom%V_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = Geom%v_mux(iph)
                  sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

                end if
              end if

              ! For each height node
              do iz=Rz0,Rz1

                ! Get index
                indx = Frec%indx(ffjtran,ia,iz,jdir)
#ifdef _OPENMP
                ! Multi-thread
                if (omp) then
                  ! Skip if not assigned
                  if (indx.lt.oif0(tid)) cycle
                  if (indx.gt.oif1(tid)) exit
                end if
#endif
                ! Calculate Doppler shift factor
                if (dyn) &
                  vfac = 1d0 - atmo%vx(iz)*st*cc - &
                               atmo%vy(iz)*st*sc - &
                               atmo%vz(iz)*ct

                ! Store limits
                if0 = Atom(ia)%if0(jtran)
                if1 = Atom(ia)%if1(jtran)


                !
                ! Calculate the Doppler width of both input and output
                ! transitions
                !

                ! Thermal common part
                DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

                ! Input
                Dw1 = nut*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

                ! Output
                Dw  = nutout*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

                ! Transform the searching parameters from normalized
                ! to proper frequency units
                red_cohwW = Input%dcohw*Dw
                red_resoW = Input%redi_pars(2)*Dw
                red_neglW = Input%redi_pars(3)*Dw
                red_coreW = Input%redi_pars(7)*Dw
                red_rangwW1 = Input%redi_pars(1)*Dw1
                red_vlarwW1 = Input%redi_pars(4)*Dw1
                red_fstpwW1 = Input%redi_pars(5)*Dw1
                red_mstpwW1 = red_fstpwW1*Input%redi_pars(6)
                red_rangcW1 = Input%redi_pars(8)*Dw1
                red_vlarcW1 = Input%redi_pars(9)*Dw1
                red_fstpcW1 = Input%redi_pars(10)*Dw1
                red_mstpcW1 = red_fstpcW1*Input%redi_pars(11)


                !
                ! Find limits for the second order output
                !
                if (.not.allocated(Frec%dzao(indx)%if0)) then

                  nran = 0
                  bf0 = -1
                  bf1 = -2

                  ! Only one output
                  if (if0.eq.if1) then

                    if(abs(Frec%omega(if0)*vfac-nutout).lt. &
                       red_neglW) then
                      bf0 = if0
                      bf1 = if1
                      nran = 1
                    end if

                  ! More than one output
                  else

                    ! Reset skip var
                    lskip = .True.

                    ! Look for the limits, for each output frequency
                    do ifreq=if0,if1

                      ! Initialize the flag
                      skip = .True.

                      ! Check if we are close to the transition
                      ! frequency
                      if(abs(Frec%omega(ifreq)*vfac-nutout).lt. &
                         red_neglW) &
                        skip = .False.

                      ! If we skip this frequency
                      if(skip.or.ifreq.eq.if1) then

                        if(.not.lskip) then

                          if (ifreq.eq.if1) then
                            bf1 = ifreq
                          else
                            bf1 = ifreq - 1
                          end if
                          nran = nran + 1

                        end if


                      ! If we cannot skip this frequency
                      else

                        ! If this is the first index for this
                        ! transition
                        if (bf0.lt.0) &
                          bf0 = ifreq

                      end if

                      lskip = skip

                    end do ! output frequencies

                  end if ! Only one output

                  ! Store in the array
                  Frec%dzao(indx)%nran = nran

                  if (Frec%dzao(indx)%nran.lt.1) cycle

                  ! Allocate the ranges
                  allocate(Frec%dzao(indx)%if0(nran))
                  Frec%dzao(indx)%if0 = bf0
                  allocate(Frec%dzao(indx)%if1(nran))
                  Frec%dzao(indx)%if1 = bf1

                  ! Only one output
                  if (if0.eq.if1) then

                    nran = 1
                    np = 1
                    Frec%dzao(indx)%if0(nran) = if0
                    Frec%dzao(indx)%if1(nran) = if0

                  ! More than one output
                  else

                    ! Only one output
                    if (bf1.eq.bf0) then

                      ! Check if we are close to a transition
                      ! frequency
                      if(abs(Frec%omega(ifreq)*vfac-nutout).lt. &
                         red_neglW) then
                        Frec%dzao(indx)%if0(1) = bf0
                        Frec%dzao(indx)%if1(1) = bf0
                        np = 1
                        nran = 1
                      else
                        Frec%dzao(indx)%nran = 0
                        deallocate(Frec%dzao(indx)%if0)
                        deallocate(Frec%dzao(indx)%if1)
                        cycle
                      end if

                    ! Multiple outputs
                    else

                      ! Reset logical variable
                      lskip = .True.

                      nran = 0
                      np = 0

                      ! Look for the limits, for each output frequency
                      do ifreq=bf0,bf1

                        ! Initialize the flag
                        skip = .True.

                        ! Check if we are close to a transition
                        ! frequency
                        if(abs(Frec%omega(ifreq)*vfac-nutout).lt. &
                           red_neglW) &
                          skip = .False.

                        ! If we skip this frequency
                        if(skip.or.ifreq.eq.bf1) then

                          if(.not.lskip) then
                            if (ifreq.eq.bf1) then
                              Frec%dzao(indx)%if1(nran) = ifreq
                              np = np + 1
                            else
                              Frec%dzao(indx)%if1(nran) = ifreq - 1
                            end if
                          end if

                        ! If we cannot skip this frequency
                        else

                          np = np + 1

                          if (lskip) then
                            nran = nran + 1
                            Frec%dzao(indx)%if0(nran) = ifreq
                          end if

                        end if

                        lskip = skip

                      end do ! output frequencies

                    end if ! Number of outputs
                  end if ! Number of outputs

                  ! Set global limits
                  Frec%dzao(indx)%gf0 = minval(Frec%dzao(indx)%if0)
                  Frec%dzao(indx)%gf1 = maxval(Frec%dzao(indx)%if1)

                  ! If storing interpolation data
                  Frec%dzao(indx)%ggf0 = 10000000
                  Frec%dzao(indx)%ggf1 = -1

                  ! Count frequencies
                  Frec%dzao(indx)%nfreq = np

                end if ! Already found limits

                if (Frec%dzao(indx)%nran.lt.1) cycle

                ! Allocate input frequency size
                np = Frec%dzao(indx)%nfreq
                allocate(Frec%dzao(indx)%trani(iti)%mfreq(np))

                ! Initialize pointers and back dimension trace
                allocate(bomega, bw_freq)
                nullify(bomega%next)
                nullify(bomega%prev)
                nullify(bw_freq%next)
                nullify(bw_freq%prev)
                bomega%nback = 0
                bw_freq%nback = 0

                ! Initialize sizes
                Frec%dzao(indx)%trani(iti)%osize = 0
                Frec%dzao(indx)%trani(iti)%isize = 0

                ! For each output frequency
                iifreq = 0
                do iran=1,Frec%dzao(indx)%nran
                  do ifreq=Frec%dzao(indx)%if0(iran), &
                           Frec%dzao(indx)%if1(iran)

                    iifreq = iifreq + 1

      !
      ! Reset indentation B
      !

      ! Reset number freq
      Frec%dzao(indx)%trani(iti)%mfreq(iifreq) = 0

      ! Coherent wings?
      if (Input%cohw) then

        ! Check distance
        if (abs(nut - Frec%omega(ifreq)*vfac).lt.red_cohwW) then
          cohw = .False.
        else
          cohw = .True.
        end if

      else

        cohw = .False.

      end if

      ! If coherent wings
      if (cohw) then

        ! Cannot be core then
        core = .False.

        !
        ! Store the determined vector
        !

        ! Advance boxes
        if (.not.allocated(bomega%A)) then

          allocate(bomega%A(1))
          allocate(bw_freq%A(1))

        else

          ! bomega
          bdaux => bomega
          allocate(bomega%next)
          bomega => bomega%next
          bomega%prev => bdaux
          allocate(bomega%A(1))
          nullify(bomega%next)
          nullify(bdaux)
          bomega%nback = bomega%prev%nback + bomega%prev%mfreq

          ! bw_freq
          bdaux => bw_freq
          allocate(bw_freq%next)
          bw_freq => bw_freq%next
          bw_Freq%prev => bdaux
          allocate(bw_freq%A(1))
          nullify(bw_freq%next)
          nullify(bdaux)

        end if

        ! Store the fequency axis
        bomega%A = Frec%omega(ifreq)*vfac
        ! Store the dimension of the axis
        bomega%mfreq = 0
        Frec%dzao(indx)%trani(iti)%mfreq(iifreq) = 0
        ! Store indexes
        bomega%ifreq = iifreq

      ! Non-coherent wings
      else

        if (abs(nut - Frec%omega(ifreq)*vfac).le.red_coreW) then
          core = .True.
        else
          core = .False.
        end if

        if (core) then
          red_rangW1 = red_rangcW1
          red_vlarW1 = red_vlarcW1
          red_fstpW1 = red_fstpcW1
          red_mstpW1 = red_mstpcW1
        else
          red_rangW1 = red_rangwW1
          red_vlarW1 = red_vlarwW1
          red_fstpW1 = red_fstpwW1
          red_mstpW1 = red_mstpwW1
        end if

        !
        ! Find the limits specified by transitions
        !

        ! Limits for input transition
        vphv(1) = nut + red_rangW1
        vplv(1) = nut - red_rangW1


        !
        ! Find the limits specified by the resonance
        !
        vpr = Frec%omega(ifreq)*vfac + dnl
        vphv(2) = vpr + red_rangW1
        vplv(2) = vpr - red_rangW1


        !
        ! Define the true limits:
        ! The resonance specify the true limits, but if
        ! they are close to a transition, we expand the limit
        !

        ! Take the resonance
        vph = vpr
        vpl = vpr


        ! If we are, move the limit to the transition
        ! instead of the resonance
        if (abs(vph - nut).lt.red_resoW) then

          if(nut.lt.vpl)vpl=nut
          if(nut.gt.vph)vph=nut

        end if

        ! Now add the range from the parameters
        vph = vph + red_rangW1
        vpl = vpl - red_rangW1


        !
        ! Flag lines and resonances out of limits
        !

        ! The total number of resonances that we had
        ni = 2

        ! Check if first one is out of limits
        if(vphv(1).lt.vpl.or.vplv(1).gt.vph)then

          ni = 1
          vphv(1) = vphv(2)
          vplv(1) = vplv(2)

          if(vphv(2).lt.vpl.or.vplv(2).gt.vph) &
            ni = 0

        end if


        ! For each resonance and transition
        do ir=1,ni

          ! The lower limit is out
          if(vplv(ir).lt.vpl) &
            vplv(ir) = vpl

          ! The upper limit is out
          if(vphv(ir).gt.vph) &
            vphv(ir) = vph

        end do


        !
        ! Check lines and resonances that spawns the same
        ! range
        !
        if (ni.gt.1) then

          ! If the ranges overlap, combine them into
          ! just one and flag the other out
          if ((vphv(2).ge.vplv(1).and.vphv(2).le.vphv(1)).or. &
              (vplv(2).ge.vplv(1).and.vplv(2).le.vphv(1))) then

            vplv(1) = min(vplv(1),vplv(2))
            vphv(1) = max(vphv(1),vphv(2))

            ni = 1

          else

            ! Order the individual limits
            if (vplv(1).gt.vplv(2)) then
              O0 = vplv(1)
              vplv(1) = vplv(2)
              vplv(2) = O0
            end if
            if (vphv(1).gt.vphv(2)) then
              O0 = vphv(1)
              vphv(1) = vphv(2)
              vphv(2) = O0
            end if

          end if

        end if ! Two ranges


        !
        ! Define extended limits
        !

        ! Reset the index of extended limits
        nie = 0

        ! For each normal limit, get an extended version
        do ir=1,ni
          nie = nie + 1
          vplve(nie) = vplv(ir) - red_vlarW1
          vphve(nie) = vplv(ir) - red_fstpW1
          nie = nie + 1
          vplve(nie) = vphv(ir) + red_fstpW1
          vphve(nie) = vphv(ir) + red_vlarW1
        end do


        !
        ! Check lines and resonances that spawns the same range
        !

        ! Store the number of original limits
        np = nie

        ! For each pair of transitions and resonances
        do ir=np,2,-1

          ! This is equivalent to be flagged out
          if(vplve(ir).lt.0d0)cycle

          ! The other index to make a pair
          do it=ir-1,1,-1

            ! If the ranges overlap, combine them into
            ! just one and flag the other out
            if(vplve(it).lt.0d0)cycle
            if((vphve(ir).ge.vplve(it).and. &
                vphve(ir).le.vphve(it)).or. &
               (vplve(ir).ge.vplve(it).and. &
                vplve(ir).le.vphve(it)))then
              vplve(it) = min(vplve(ir),vplve(it))
              vphve(it) = max(vphve(ir),vphve(it))
              vplve(ir) = -1d0
              nie = nie - 1
              exit
            end if

          end do ! it
        end do ! ir


        !
        ! Shift the individual limits moving the valid ones
        ! to the first part of the vector
        !

        ! If we have changed the number of ranges from the
        ! beginning
        if(nie.ne.np)then

          ! For each pair of resonances
          do ir=1,np-1

            ! If it is flagged out
            if(vplve(ir).lt.0d0)then

              ! Reset the added index
              ip = 1

              ! For all the resonances in front of this one
              do it=ir+1,np

                ! If it is flagged, skip it
                if(vplve(it).lt.0d0)then

                  ip = ip + 1

                ! If it is not flagged, move it to the
                ! position of the flagged one
                else

                  vplve(it-ip) = vplve(it)
                  vphve(it-ip) = vphve(it)
                  vplve(it) = -1d0
                  exit

                end if
              end do ! it
            end if ! ir flagged
          end do ! ir
        end if ! if something flagged


        !
        ! Order the individual limits
        !

        ! Lower limits
        call QsortC(vplve(1:nie))
        ! Upper limits
        call QsortC(vphve(1:nie))


        !
        ! Build the vector of input frequencies from the
        ! limits
        !

        ! Flag to reset in case we run out of space while
        ! building a vector
        reset = .False.

        ! Do until we are finished
        do while (.True.)

          ! If the flag is one, we need more space for the
          ! vector
          if(reset)then
            np0 = np0*2
            deallocate(vpp)
            allocate(vpp(np0))
            deallocate(flag)
            allocate(flag(np0))
            reset = .False.
          end if

          ! Reset index counter
          np = 0


          !
          ! Build for the short limits
          !

          ! For each short range
          do it=1,ni

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if(np.gt.np0)then
              reset = .True.
              exit
            end if

            ! Start with the lower limit
            vpp(np) = vplv(it)

            ! Do until finished
            do while(.True.)

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Next frequency is the previous plus the step
              vpp(np) = vpp(np-1) + red_fstpW1

              ! If we are over the range, we are done
              if(vpp(np).ge.vphv(it))then
                vpp(np) = vphv(it)
                exit
              end if

            end do

            ! If we have no space, we need to reset
            if (reset) exit

          end do

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          !
          ! Build for the extended limits
          !

          ! For each long range
          do it=1,nie

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if(np.gt.np0)then
              reset = .True.
              exit
            end if

            ! We start with the lower limit
            vpp(np) = vplve(it)

            ! Do until finished
            do while(.True.)

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Next frequency is the previous plus the step
              vpp(np) = vpp(np-1) + red_mstpW1

              ! If we are over the range, we are done
              if(vpp(np).ge.vphve(it))then
                vpp(np) = vphve(it)
                exit
              end if

            end do

            ! If we have no space, we need to allocate it
            ! above
            if (reset) exit

          end do

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          !
          ! Add the resonance frequencies to the vector
          !

          ! The number of frequencies we already have
          npp = np

          !
          ! Add the resonance
          !

          ! Reset the flag
          nfound = .True.

          ! Run over the existing frequencies
          do ip=1,npp

            ! If the frequency is there, do not add it
            if (abs(1d2/vpp(ip) - 1d2/vpr).lt.resolin) then
              nfound = .False.
              exit
            end if

          end do

          ! If we did not find it
          if (nfound) then

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if(np.gt.np0)then
              reset = .True.
            end if

            ! Add the frequency
            if (.not.reset) vpp(np) = vpr

          end if

          if(reset)cycle

          exit ! If we get to this point, we have everything

        end do ! First do while


        !
        ! Check for duplicates
        !

        ! Reset the flag
        flag(1:np) = 1

        ! For each frequency
        do ip=1,np

          ! If it has been flagged, we already checked
          if (flag(ip).lt.1) cycle

          ! Check the following ones
          do ipp=ip+1,np

            ! If it has been flagged, we already checked
            if (flag(ipp).lt.1) cycle

            ! If some of them are repeated, flag them to be removed
            if(abs(1d2/vpp(ip)-1d2/vpp(ipp)).lt.resolin) flag(ipp) = 0

          end do ! ipp
        end do ! ip

        ! Reset the running real index
        ipp = 0

        ! For each frequency in the vector
        do ip=1,np

          ! If it is flagged correct, add to real vector
          if(flag(ip).gt..5)then
            ipp = ipp + 1
            vpp(ipp) = vpp(ip)
          end if

        end do

        ! The number of frequencies is the last value of ipp
        np = ipp

        ! Order the frequency axis
        call QsortC(vpp(1:np))


        !
        ! Store the determined vector
        !

        ! Advance boxes
        if (.not.allocated(bomega%A)) then

          allocate(bomega%A(np))
          allocate(bw_freq%A(np))

        else

          ! bomega
          bdaux => bomega
          allocate(bomega%next)
          bomega => bomega%next
          bomega%prev => bdaux
          allocate(bomega%A(np))
          nullify(bomega%next)
          nullify(bdaux)
          bomega%nback = bomega%prev%nback + bomega%prev%mfreq

          ! bw_freq
          bdaux => bw_freq
          allocate(bw_freq%next)
          bw_freq => bw_freq%next
          bw_Freq%prev => bdaux
          allocate(bw_freq%A(np))
          nullify(bw_freq%next)
          nullify(bdaux)

        end if

        ! Check that we have enough space to work below
        if(.not.allocated(Wvpp))then
          allocate(Wvpp(np))
        else
          if(size(Wvpp).lt.np)then
            deallocate(Wvpp)
            allocate(Wvpp(np))
          end if
        end if

        ! Store the fequency axis
        bomega%A = vpp(1:np)
        ! Store the dimension of the axis
        bomega%mfreq = np
        Frec%dzao(indx)%trani(iti)%mfreq(iifreq) = np
        ! Store indexes
        bomega%ifreq = iifreq

        ! Update the maximum of input frequencies
        if (np.gt.Frec%dzao(indx)%mxfreq) Frec%dzao(indx)%mxfreq = np


        !
        ! Define the integration weights (same than
        ! omegabuild)
        !

        ! The first point is special in compound trapezoidal
        ! rule
        Wvpp(1) = .5d0*(vpp(2) - vpp(1))

        ! Initialize the integral to normalize the weights
        norm1 = Wvpp(1)

        ! The initial lower limit is the first point
        O0 = vpp(1)

        ! This is the pointer to the first element of the
        ! current interval, we are pointing to the first
        ! element
        cfreq = 1

        ! Flag that says that the point 2 is not the initial
        ! point of the interval (because 1 is the initial
        ! point)
        init = .FALSE.

        ! For the rest of frequencies except the last
        do jfreq=2,np-1

          ! If ifreq is the initial point of an interval
          if(init)then

            ! The first point is special in compound
            ! trapezoidal rule
            Wvpp(jfreq) = .5d0*(vpp(jfreq+1) - vpp(jfreq))

            ! The next point cannot be a first point
            init = .FALSE.

            ! Initialize the integral to normalize the
            ! weights
            norm1 = Wvpp(jfreq)

            ! Pointer is now in this frequency
            cfreq = jfreq

            ! And it is the beginning of the current interval
            O0 = vpp(jfreq)

          ! If jfreq is not the initial point of an interval
          else

            ! Check if ifreq is the last point of an interval
            if(abs(vpp(jfreq+1) - vpp(jfreq)).gt.red_neglW)then

              ! The last point is special in compound
              ! trapezoidal rule
              Wvpp(jfreq) = .5d0*(vpp(jfreq) - vpp(jfreq-1))

              ! It is the end of the current interval
              O1 = vpp(jfreq)

              ! Add to the integral
              norm1 = norm1 + Wvpp(jfreq)

              ! We know that the integral must be
              ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
              norm = 1d5*(O1 - O0)/norm1

              ! Normalize the weights of this interval
              do kfreq=cfreq,jfreq
                Wvpp(kfreq) = Wvpp(kfreq)*norm
              end do

              ! The next point is the first point of its
              ! interval
              init = .TRUE.

            ! If jfreq is not the last point of an interval
            else

              ! Compound trapezoidal rule weight
              Wvpp(jfreq) = .5d0*(vpp(jfreq+1) - vpp(jfreq-1))

              ! Add to the integral
              norm1 = norm1 + Wvpp(jfreq)

            endif ! Last point

          end if ! Initial point

        end do ! jfreq

        ! The last point is special in compound trapezoidal
        ! rule
        Wvpp(np) = .5d0*(vpp(np) - vpp(np-1))

        ! It is the end of the interval
        O1 = vpp(np)

        ! Add to the integral
        norm1 = norm1 + Wvpp(np)

        ! We know that the integral must be
        ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
        norm = 1d5*(O1 - O0)/norm1

        ! Normalize the weights of this interval
        do jfreq=cfreq,np
          Wvpp(jfreq) = Wvpp(jfreq)*norm
        end do

        ! Store the weights
        bw_freq%A = Wvpp(1:np)

      end if ! Coherent wings

                  end do ! Output frequencies
                end do ! output frequency ranges

                !
                ! Properly store and index the data
                !

                ! Total dimension of omega and w_freq
                nn = bomega%nback + bomega%mfreq

                ! Allocate omega and W_freq
                allocate(Frec%dzao(indx)%trani(iti)%omega(nn))
                allocate(Frec%dzao(indx)%trani(iti)%w_freq(nn))

                ! Determine size
                Frec%dzao(indx)%trani(iti)%osize = nn

                ! Go backwards in the linked lists
                do while (.True.)

                  iifreq = bomega%ifreq
                  ip = bomega%nback + 1
                  ipp = ip + bomega%mfreq - 1
                  if (ipp.ge.ip) then
                    Frec%dzao(indx)%trani(iti)%omega(ip:ipp) = &
                                                              bomega%A
                    Frec%dzao(indx)%trani(iti)%W_freq(ip:ipp) = &
                                                             bw_freq%A
                  end if

                  ! Deallocate arrays
                  deallocate(bomega%A,bw_freq%A)

                  ! If last one, clean and quit
                  if (.not.associated(bomega%prev)) then
                    deallocate(bomega,bw_freq)
                    nullify(bomega,bw_freq)
                    exit
                  ! Not done with the list
                  else
                    bomega => bomega%prev
                    bw_freq => bw_freq%prev
                    nullify(bomega%next%prev,bw_freq%next%prev)
                    deallocate(bomega%next,bw_freq%next)
                    nullify(bomega%next,bw_freq%next)
                  end if

                end do ! Run backwards the frequency axes

                ! Update RAM
!$omp flush(MPID)
                MPID%RAM = MPID%RAM + 16d-6*dble(nn)
                MPID%WRAM = MPID%WRAM + 16d-6*dble(nn)
!$omp flush(MPID)

                ! End of array constructions

              end do ! Heights
            end do !Output directions
          end do ! Lower input level
        end do ! Lower input term

              !
              ! Restore identation
              !

              end do ! Lower output level
            end do ! Lower output term
          end do ! upper level
        end do ! upper term

!$omp barrier

        !
        ! Allocate space for interpolation and define it or
        ! find the index limits
        !

        ! For each output direction
        do jdir=1,Frec%ndir

          ! For each height
          do iz=Rz0,Rz1

            ! For each upper level
            do itermu=2,Atom(ia)%nMulti
              do iJu=1,Atom(ia)%nJ(itermu)

                ! For each final lower level
                do itermf=1,itermu-1
                  do iJf=1,Atom(ia)%nJ(itermf)

                    jtran = Atom(ia)%irad(itermu,itermf)

                    if (jtran.le.0) cycle
                    if (.not.Atom(ia)%lemiss2(jtran)) cycle
                    if (Atom(ia)%fflag(jtran)%absent) cycle

                    fjtran = Atom(ia)%fst(jtran)%irad(iJu,iJf)

                    if (fjtran.le.0) cycle

                    ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                    ffktran = Atom(ia)%itrano(ffjtran)

                    ! Get index
                    indx = Frec%indx(ffjtran,ia,iz,jdir)
#ifdef _OPENMP
                    ! Multi-thread
                    if (omp) then
                      ! Skip if not assigned
                      if (indx.lt.oif0(tid)) cycle
                      if (indx.gt.oif1(tid)) exit
                    end if
#endif
                    if (Frec%dzao(indx)%nran.lt.1) cycle


        !
        ! Reset identation
        !

        ! For each other lower level
        do iterml=1,itermu-1
          do iJl=1,Atom(ia)%nJ(iterml)

            itran = Atom(ia)%irad(itermu,iterml)
            if (itran.le.0) cycle

            fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

            if (fitran.le.0) cycle

            ffitran = Atom(ia)%ifst_ij(fitran,itran)

            ! Find the transition index
            ios = -1
            do iti=1,Atom(ia)%trano(ffktran)%nt

              if (Atom(ia)%trano(ffktran)%ind(iti).eq.ffitran) then
                ios = 1
                exit
              end if

            end do ! Input transitions
            if (ios.lt.0) cycle

            ! Point to input frequency
            p_frec => Frec%dzao(indx)%trani(iti)

            ! Predict size of next block
            nn = sum(p_frec%mfreq)

            ! If angle-dependent
            if (.not.AVI) then

              ! If dynamic extra dimensions, if static just
              ! frequencies
              if (dyn) then

                ! For axial problems
                if (axiali) then

                  ! Size is just polar
                  nn = nn*Geom%nTh

                ! For non-axial problems
                else

                  ! Skip backward rayleigh
                  nn = nn*(Geom%nTh*Geom%nPh2 - Frec%nfs(jdir))

                end if ! Axial
              end if ! Dynamic
            end if ! AD

            ! Predict aditional frequency
            SRAM = 16d-6*dble(nn)

            ! If can store
            if (TIRAM) then
!$omp flush(MPID,ofram)
              ! If no more space
              if (floor(MPID%RAM+SRAM).gt.RLIM.or.SRAM.le.0d0) then
                ofram = .True.
                p_frec%RAM = .False.
              else
                allocate(p_frec%index1(nn))
                allocate(p_frec%index2(nn))
                allocate(p_frec%dx(nn))
                MPID%RAM = MPID%RAM + SRAM
                MPID%WRAM = MPID%WRAM + SRAM
                p_frec%RAM = .True.
              end if
!$omp flush(MPID,ofram)
            ! Cannot store
            else
              p_frec%RAM = .False.
            end if

            ! Store size
            Frec%dzao(indx)%trani(iti)%isize = nn

            !
            ! Define interpolation

            ! Initialize index
            jjfreq0 = 0
            kkfreq0 = 0

            ! For each output frequency
            iifreq = 0
            do iran=1,Frec%dzao(indx)%nran
              do ifreq=Frec%dzao(indx)%if0(iran), &
                       Frec%dzao(indx)%if1(iran)

                ! Advance index
                iifreq = iifreq + 1

                ! Input frequency number
                np = p_frec%mfreq(iifreq)

                ! For each input direction
                do ith=1,Frec%nth
                  do iph=1,Frec%nph

                    ! If dynamics and AD
                    if (dyn.and..not.AVI) then

                      ! For axial problems
                      if (axiali) then

                        ! Automatically skip extra azimuths
                        if (iph.gt.1) cycle

                        ! Get director cosines
                        ct1 = Geom%V_mu(ith)

                        ! Calculate Doppler shift factor
                        vfac1 = 1d0 - atmo%vz(iz)*ct1

                        ! We will be using the inverse
                        vfac1 = 1d0/vfac1

                      ! For non-axial problems
                      else

                        ! If angle-dependent, check backward Rayleigh
                        ! scattering
                        if (ffjtran.eq.ffitran.and. &
                            Frec%stype(iph,ith,jdir).lt.0) cycle

                        ! Get director cosines
                        ct1 = Geom%V_mu(ith)
                        st1 = sqrt(1d0 - ct1*ct1)
                        cc1 = Geom%v_mux(iph)
                        sc1 = Geom%v_muy(iph)*sqrt(1d0 - cc1*cc1)

                        ! Calculate Doppler shift factor
                        vfac1 = 1d0 - atmo%vx(iz)*st1*cc1 - &
                                      atmo%vy(iz)*st1*sc1 - &
                                      atmo%vz(iz)*ct1

                        ! We will be using the inverse
                        vfac1 = 1d0/vfac1

                      end if ! Axial

                    ! Not dynamic or AV
                    else

                      ! No shift
                      vfac1 = 1d0

                      ! Only one direction
                      if (iph.gt.1.or.ith.gt.1) cycle

                    end if ! Dynamics

                    ! Reset indexes
                    jjfreq = jjfreq0
                    kkfreq = kkfreq0

                    ! Skip empty
                    if (np.lt.1) then

                      ! Left and right limits from resonance
                      lifreq = 1
                      jfreq = nfreq
                      O0 = Frec%omega(ifreq)*vfac

                      !
                      ! Look for the indexes

                      ! Left
                      do while (.True.)

                        ! Check if already inside
                        if (Frec%omega(lifreq)*vfac1.ge.O0) exit

                        lifreq = lifreq + 1
                        if (lifreq.gt.nfreq) exit

                      end do

                      ! Right
                      do while (.True.)

                        ! Check if already inside
                        if (Frec%omega(jfreq)*vfac1.le.O0) exit

                        jfreq = jfreq - 1
                        if (jfreq.lt.1) exit

                      end do

                      ! Update global limits
                      if (lifreq.lt.Frec%dzao(indx)%ggf0) &
                        Frec%dzao(indx)%ggf0 = lifreq
                      if (jfreq.gt.Frec%dzao(indx)%ggf1) &
                        Frec%dzao(indx)%ggf1 = jfreq

                      ! Skip rest
                      cycle

                    end if ! Coherent wing

                    ! If storing
                    if (p_frec%RAM) then

                      !
                      ! Reset identation
                      !

      ! Reset the search frequency
      lifreq = 1

      ! For each input frequency
      do jfreq=1,np

        ! Advance indexes
        jjfreq = jjfreq + 1
        kkfreq = kkfreq + 1

        ! If out of range, take the value at the
        ! boundary
        if (p_frec%omega(jjfreq)*vfac1.le.Frec%omega(1)+TINYO) then

          ! We are still looking in the first one
          lifreq = 1

          ! The index to take is 1
          p_frec%index1(kkfreq) = 1

          ! The index to take is 1
          p_frec%index2(kkfreq) = 1

          ! We do not need this number
          p_frec%dx(kkfreq) = 0d0

        ! If out of range, take the value at the boundary
        else if (p_frec%omega(jjfreq)*vfac1.ge. &
                 (Frec%omega(nfreq) - TINYO)) then

          ! We are in the last frequency
          lifreq = nfreq

          ! The index to take is nfreq
          p_frec%index1(kkfreq) = nfreq

          ! The index to take is nfreq
          p_frec%index2(kkfreq) = nfreq

          ! We do not need this number
          p_frec%dx(kkfreq) = 0d0

        ! If within the boundaries
        else

          ! Search between the last found frequency and
          ! all but the boundary
          do ibfreq=lifreq,nfreq-1

            ! If this exact frequency is in output
            if (abs(p_frec%omega(jjfreq)*vfac1 - &
                    Frec%omega(ibfreq)).lt.TINYO) then

              ! We are in the found frequency
              lifreq = ibfreq

              ! This frequency gives us the value
              p_frec%index1(kkfreq) = lifreq

              ! This frequency gives us the value
              p_frec%index2(kkfreq) = lifreq

              ! We do not need this number
              p_frec%dx(kkfreq) = 0d0

              exit

            ! If the input is between this output and
            ! the next
            else if(p_frec%omega(jjfreq)*vfac1.ge. &
                    Frec%omega(ibfreq).and. &
                    p_frec%omega(jjfreq)*vfac1.lt. &
                    Frec%omega(ibfreq+1)) then

              ! We found it in the index of the lower
              lifreq = ibfreq

              ! The first index is the lower
              p_frec%index1(kkfreq) = lifreq

              ! The second index is the upper
              p_frec%index2(kkfreq) = lifreq+1

              ! Store the inverse of the distance between
              ! the two outputs
              p_frec%dx(kkfreq) = &
                  (p_frec%omega(jjfreq)*vfac1 - Frec%omega(lifreq))/ &
                  (Frec%omega(lifreq+1) - Frec%omega(lifreq))

              exit

            end if ! Check output frequency

          end do ! Run output frequencies

        end if ! Check if out of limits

        ! Update global limits
        if (p_frec%index1(kkfreq).lt.Frec%dzao(indx)%ggf0) &
          Frec%dzao(indx)%ggf0 = p_frec%index1(kkfreq)
        if (p_frec%index2(kkfreq).gt.Frec%dzao(indx)%ggf1) &
          Frec%dzao(indx)%ggf1 = p_frec%index2(kkfreq)

      end do ! Run input frequencies

                      !
                      ! Restore identation
                      !

                    ! Not storing
                    else

                      !
                      ! Reset identation
                      !

       ! Compute jump
       lifreq = 1
       jufreq = np
       if (np.gt.1) jufreq = jufreq - 1

       ! First and last frequencies
       do jfreq=jjfreq+1,jjfreq+np,jufreq

         ! If out of range, take the value at the
         ! boundary
         if (p_frec%omega(jfreq)*vfac1.le.Frec%omega(1)) then

           if (Frec%dzao(indx)%ggf0.gt.1) Frec%dzao(indx)%ggf0 = 1
           if (Frec%dzao(indx)%ggf1.lt.1) Frec%dzao(indx)%ggf1 = 1

         ! If out of range, take the value at the boundary
         else if (p_frec%omega(jfreq)*vfac1.ge. &
                  (Frec%omega(nfreq) - TINYO)) then

           ! We are in the last frequency
           if (Frec%dzao(indx)%ggf0.gt.nfreq) &
             Frec%dzao(indx)%ggf0 = nfreq
           if (Frec%dzao(indx)%ggf1.lt.nfreq) &
             Frec%dzao(indx)%ggf1 = nfreq

         ! If within the boundaries
         else

           ! Search between the last found frequency and
           ! all but the boundary
           do ibfreq=lifreq,nfreq-1

             ! If this exact frequency is in output
             if (abs(p_frec%omega(jfreq)*vfac1 - &
                     Frec%omega(ibfreq)).lt.TINYO) then

               ! Found frequency
               if (Frec%dzao(indx)%ggf0.gt.ibfreq) &
                 Frec%dzao(indx)%ggf0 = ibfreq
               if (Frec%dzao(indx)%ggf1.lt.ibfreq) &
                 Frec%dzao(indx)%ggf1 = ibfreq

               exit

             ! If the input is between this output and the next
             else if(p_frec%omega(jfreq)*vfac1.ge. &
                                            Frec%omega(ibfreq).and. &
                     p_frec%omega(jfreq)*vfac1.lt. &
                                            Frec%omega(ibfreq+1)) then

               ! Found frequency
               if (Frec%dzao(indx)%ggf0.gt.ibfreq) &
                 Frec%dzao(indx)%ggf0 = ibfreq
               if (Frec%dzao(indx)%ggf1.lt.ibfreq+1) &
                 Frec%dzao(indx)%ggf1 = ibfreq+1

               exit

             end if ! Check output frequency

           end do ! Run output frequencies

         end if ! Check if out of limits

       end do ! Run input frequencies

       ! Fake the advance of frequencies
       jjfreq = jjfreq + np
       kkfreq = kkfreq + np


                      !
                      ! Restore indentation
                      !

                    end if ! Storing

                    if (.not.AVI.and..not.axiali.and.dyn) &
                      kkfreq0 = kkfreq

                  end do ! Input azimuth

                  ! Update kkfreq
                  if (.not.AVI.and.dyn) kkfreq0 = kkfreq

                end do ! Input polar

                ! Update index
                jjfreq0 = jjfreq
                kkfreq0 = kkfreq

              end do ! Output frequencies
            end do ! Output frequency ranges

            ! Nullify local pointer
            nullify(p_frec)

          end do ! iJl
        end do ! iterml

                  !
                  ! Restore identation
                  !

                  end do ! Lower output level
                end do ! Lower output term
              end do ! upper level
            end do ! upper term
          end do ! height nodes
        end do ! output directions

!$omp barrier

        !
        ! Allocate space for Warr2
        !
        if (IRAM) then

          ! For each output direction
          do jdir=1,Red%ndir

            jbdir = min(jdir,Frec%ndir)
            jcdir = min(jdir,Red%njdir)

            ! For each height
            do iz=Rz0,Rz1

              ! For each upper level
              do itermu=2,Atom(ia)%nMulti
                do iJu=1,Atom(ia)%nJ(itermu)

                  ! For each final lower level
                  do itermf=1,itermu-1
                    do iJf=1,Atom(ia)%nJ(itermf)

                      jtran = Atom(ia)%irad(itermu,itermf)

                      if (jtran.le.0) cycle
                      if (.not.Atom(ia)%lemiss2(jtran)) cycle
                      if (Atom(ia)%fflag(jtran)%absent) cycle

                      fjtran = Atom(ia)%fst(jtran)%irad(iJu,iJf)

                      if (fjtran.le.0) cycle

                      ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                      ffktran = Atom(ia)%itrano(ffjtran)

                      ! Get index
                      indxf = Frec%indx(ffjtran,ia,iz,jbdir)
                      if (Frec%dzao(indxf)%nran.lt.1) cycle

                      ! Get index
                      indx = Red%indx(ffjtran,ia,iz,jdir)
#ifdef _OPENMP
                      ! If multi-thread
                      if (omp) then
                        ! Skip if not assigned
                        if (indx.lt.orif0(tid)) cycle
                        if (indx.gt.orif1(tid)) exit
                      end if
#endif

        !
        ! Reset identation
        !

        ! For each other lower level
        do iterml=1,itermu-1
          do iJl=1,Atom(ia)%nJ(iterml)

            itran = Atom(ia)%irad(itermu,iterml)
            if (itran.le.0) cycle

            fitran = Atom(ia)%fst(itran)%irad(iJu,iJl)

            if (fitran.le.0) cycle

            ffitran = Atom(ia)%ifst_ij(fitran,itran)

            ! Find the transition index
            ios = -1
            do iti=1,Atom(ia)%trano(ffktran)%nt

              if (Atom(ia)%trano(ffktran)%ind(iti).eq.ffitran) then
                ios = 1
                exit
              end if

            end do ! Input transitions
            if (ios.lt.0) cycle

            ! Predict size of next block
            nn = sum(Frec%dzao(indxf)%trani(iti)%mfreq)
            if (ffjtran.eq.ffitran) then
              nn = nn*(Red%nth*Red%nph - Frec%nfs(jcdir))
            else
              nn = nn*Red%nth*Red%nph
            end if
            SRAM = 4d-6*dble(nn)

            ! If no more space
!$omp flush(MPID,ofram)
            if (floor(MPID%RAM+SRAM).gt.RLIM.or. &
                SRAM.le.0d0) then
              ofram = .True.
!$omp flush(ofram)
              Red%dzao(indx)%trani(iti)%RAM = .False.
              cycle
            end if
            MPID%RAM = MPID%RAM + SRAM
            MPID%WRAM = MPID%WRAM + SRAM
!$omp flush(MPID)

            Red%dzao(indx)%trani(iti)%RAM = .True.
            allocate(Red%dzao(indx)%trani(iti)%IWarr2(nn))

          end do ! iJl
        end do ! iterml

                  !
                  ! Restore identation
                  !

                    end do ! Lower output level
                  end do ! Lower output term
                end do ! upper level
              end do ! upper term
            end do ! height nodes
          end do ! output directions

        ! If not storing in RAM
        else

          ! Allocate dummy array
!$omp single
          if (.not.associated(Red%dzao)) then
            allocate(Red%dzao(1))
            nullify(Red%dzao(1)%trani)
          end if
!$omp end single

        end if ! IRAM
      end do ! Atom
!$omp end parallel

#ifdef _OPENMP
      ! deallocate limits for threads
      if (omp) deallocate(oif0,oif1)
#endif

      ! Check if everything is fine
      call control

      return

      end subroutine omegabuildinI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resize some frequency dependent quantities and set CPU wise
      !! limits for the indexes where transitions are present.\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!        Input(Input_class): Structure with settings data\n
      !!           MPID(MPI_class): Structure with MPI data
      subroutine frecresize(Frec,Atom,Input,MPID)

      ! I/O

      type(Atom_class), dimension(:):: Atom
      type(Frequency_class):: Frec
      type(Input_class):: Input
      type(MPI_class), intent(in):: MPID

      ! Local

      integer:: ia,itran,iproc,if0,if1

      double precision, dimension(nfreq):: daux

      ! If no MPI
      if (.not.MPID%mpi) then
        if (Frec%ntfreqi.lt.1) Frec%ntfreqi = 1
        if (Frec%npfreq.lt.1) Frec%npfreq = 1
        if (Frec%ntfreq.lt.1) Frec%ntfreq = 1
        return
      end if


      ! Allocate size of profiles
      allocate(Frec%Mntfreq(0:nproc-1))
      allocate(Frec%Mntfreqi(0:nproc-1))
      allocate(Frec%Mnpfreq(0:nproc-1))


      !
      ! If master
      !
      if (pid.eq.0.and.MPID%mpi) then

        ! Limits of photoionizations
        if (pid.eq.0.and.nproc.gt.1) then
          allocate(Frec%Mlif0(0:nproc-1))
          allocate(Frec%Mlif1(0:nproc-1))
          allocate(Frec%Mpif0(0:nproc-1))
          allocate(Frec%Mpif1(0:nproc-1))
          Frec%Mlif0(0) = 0
          Frec%Mlif1(0) = 0
          Frec%Mpif0(0) = 0
          Frec%Mpif1(0) = 0
        end if

        ! Reset Mntfreq
        Frec%Mntfreq(0) = Frec%ntfreq
        Frec%Mntfreqi(0) = Frec%ntfreqi
        Frec%Mnpfreq(0) = Frec%npfreq
        Frec%Mntfreq(1:nproc-1) = 0
        Frec%Mntfreqi(1:nproc-1) = 0
        Frec%Mnpfreq(1:nproc-1) = 0

        ! For each CPU
        do iproc=1,nproc-1

          if0 = MPID%if0(iproc)
          if1 = MPID%if1(iproc)


          ! Initialize master limits
          Frec%Mlif0(iproc) = 100000000
          Frec%Mlif1(iproc) = -1
          Frec%Mpif0(iproc) = 100000000
          Frec%Mpif1(iproc) = -1

          ! For each atom
          do ia=1,nA

            ! For each b-b transition
            do itran=1,Atom(ia)%ntran

              ! Rearrange the limits taking into account the CPU
              ! limits

              ! The line is totally absent
              if (Atom(ia)%Mif0(itran,iproc).gt.if1.or. &
                  Atom(ia)%Mif1(itran,iproc).lt.if0) then
                Atom(ia)%Mif0(itran,iproc) = if1
                Atom(ia)%Mif1(itran,iproc) = if0-1
                Atom(ia)%MW0(itran,iproc) = 0d0
                Atom(ia)%MW1(itran,iproc) = 0d0
                Atom(ia)%fflag(itran)%Mabsent(iproc) = .True.

              ! The line is present
              else

                ! If the lower limit is out of range
                if (Atom(ia)%Mif0(itran,iproc).lt.if0) then
                  Atom(ia)%Mif0(itran,iproc) = if0
                  Atom(ia)%MW0(itran,iproc) = Frec%W_freq(if0)
                end if

                ! If the upper limit is out of range
                if (Atom(ia)%Mif1(itran,iproc).gt.if1) then
                  Atom(ia)%Mif1(itran,iproc) = if1
                  Atom(ia)%MW1(itran,iproc) = Frec%W_freq(if1)
                end if

                ! If there is only one frequency in this CPU
                if (Atom(ia)%Mif0(itran,iproc).eq. &
                    Atom(ia)%Mif1(itran,iproc)) &
                  Atom(ia)%MW1(itran,iproc) = 0d0

                ! Add frequencies to count
                Frec%Mntfreq(iproc) = Frec%Mntfreq(iproc) + 1 + &
                                      Atom(ia)%Mif1(itran,iproc) - &
                                      Atom(ia)%Mif0(itran,iproc)
                Frec%Mntfreqi(iproc) = Frec%Mntfreqi(iproc) + (1 + &
                                       Atom(ia)%Mif1(itran,iproc) - &
                                       Atom(ia)%Mif0(itran,iproc))* &
                                       Atom(ia)%fst(itran)%nt

              end if ! Line presence

              ! Handle the line range
              if (Atom(ia)%Mif0(itran,iproc).lt.Frec%Mlif0(iproc)) &
                Frec%Mlif0(iproc) = Atom(ia)%Mif0(itran,iproc)
              if (Atom(ia)%Mif1(itran,iproc).gt.Frec%Mlif1(iproc)) &
                Frec%Mlif1(iproc) = Atom(ia)%Mif1(itran,iproc)

            end do ! b-b Transition

            ! For each b-f transition
            do itran=1,Atom(ia)%nphot

              ! The line is totally absent
              if (Atom(ia)%phot(itran)%Mif0(iproc).gt.if1.or. &
                  Atom(ia)%phot(itran)%Mif1(iproc).lt.if0) then
                Atom(ia)%phot(itran)%Mif0(iproc) = if1
                Atom(ia)%phot(itran)%Mif1(iproc) = if0-1
                Atom(ia)%phot(itran)%MW0(iproc) = 0d0
                Atom(ia)%phot(itran)%MW1(iproc) = 0d0
                Atom(ia)%phot(itran)%Mabsent(iproc) = .True.

              ! The line is present
              else

                ! If the lower limit is out of range
                if (Atom(ia)%phot(itran)%Mif0(iproc).lt.if0) then
                  Atom(ia)%phot(itran)%Mif0(iproc) = if0
                  Atom(ia)%phot(itran)%MW0(iproc) = Frec%W_freq(if0)
                end if

                ! If the upper limit is out of range
                if (Atom(ia)%phot(itran)%Mif1(iproc).gt.if1) then
                  Atom(ia)%phot(itran)%Mif1(iproc) = if1
                  Atom(ia)%phot(itran)%MW1(iproc) = Frec%W_freq(if1)
                end if

                ! If there is only one frequency in this CPU
                if (Atom(ia)%phot(itran)%Mif0(iproc).eq. &
                    Atom(ia)%phot(itran)%Mif1(iproc)) &
                  Atom(ia)%phot(itran)%MW1(iproc) = 0d0

                ! Add frequencies to count
                Frec%Mnpfreq(iproc) = Frec%Mnpfreq(iproc) + 1 + &
                                  Atom(ia)%phot(itran)%Mif1(iproc) - &
                                  Atom(ia)%phot(itran)%Mif0(iproc)

              end if ! Line presence

              ! Handle the photoionization range
              if (Atom(ia)%phot(itran)%Mif0(iproc).lt. &
                  Frec%Mpif0(iproc)) &
                Frec%Mpif0(iproc) = Atom(ia)%phot(itran)%Mif0(iproc)
              if (Atom(ia)%phot(itran)%Mif1(iproc).gt. &
                  Frec%Mpif1(iproc)) &
                Frec%Mpif1(iproc) = Atom(ia)%phot(itran)%Mif1(iproc)

            end do ! b-f Transition

          end do ! Atom

          ! Control limits
          if (Frec%Mntfreqi(iproc).lt.1) Frec%Mntfreqi(iproc) = 1
          if (Frec%Mnpfreq(iproc).lt.1) Frec%Mnpfreq(iproc) = 1
          if (Frec%Mntfreq(iproc).lt.1) Frec%Mntfreq(iproc) = 1

        end do ! CPU

      !
      ! If slave
      !
      else

        ! Store MPI limits in short variables
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! Resize the true array
        if(size(Frec%W_freq).ne.MPID%nf(pid))then

          ! Store into a temporal array
          daux = Frec%W_freq

          ! Resize
          deallocate(Frec%W_freq)
          allocate(Frec%W_freq(if0:if1))

          ! Recover the data
          Frec%W_freq = daux(if0:if1)

          ! For each atom
          do ia=1,nA

            ! For each b-b transition
            do itran=1,Atom(ia)%ntran

              ! Skip if already absent
              if (Atom(ia)%fflag(itran)%absent) cycle

              ! If upper line limit above lower CPU limit
              if (Atom(ia)%if1(itran).lt.MPID%if0(pid)) then
                Atom(ia)%fflag(itran)%absent = .True.
                cycle
              end if

              ! If lower line limit above upper CPU limit
              if (Atom(ia)%if0(itran).gt.MPID%if1(pid)) then
                Atom(ia)%fflag(itran)%absent = .True.
                cycle
              end if

            end do ! b-b Transition
          end do ! Atom

        end if ! Need to resize

        ! Reset number of frequencies for profiles
        Frec%ntfreq = 0
        Frec%ntfreqi = 0
        Frec%npfreq = 0

        ! For each atom
        do ia=1,nA

          ! For each b-b transition
          do itran=1,Atom(ia)%ntran

            ! Rearrange the limits taking into account the CPU limits

            ! The line is totally absent
            if (Atom(ia)%fflag(itran)%absent) then
              Atom(ia)%if0(itran) = if1
              Atom(ia)%if1(itran) = if0-1
              Atom(ia)%W0(itran) = 0d0
              Atom(ia)%W1(itran) = 0d0

            ! The line is present
            else

              ! If the lower limit is out of range
              if (Atom(ia)%if0(itran).lt.if0) then
                Atom(ia)%if0(itran) = if0
                Atom(ia)%W0(itran) = Frec%W_freq(if0)
              end if

              ! If the upper limit is out of range
              if (Atom(ia)%if1(itran).gt.if1) then
                Atom(ia)%if1(itran) = if1
                Atom(ia)%W1(itran) = Frec%W_freq(if1)
              end if

              ! If there is only one frequency in this CPU
              if (Atom(ia)%if0(itran).eq.Atom(ia)%if1(itran)) &
                Atom(ia)%W1(itran) = 0d0

              Frec%ntfreq = Frec%ntfreq + 1 + &
                            Atom(ia)%if1(itran) - Atom(ia)%if0(itran)
              Frec%ntfreqi = Frec%ntfreqi + (1 + &
                              Atom(ia)%if1(itran) - &
                              Atom(ia)%if0(itran))* &
                             Atom(ia)%fst(itran)%nt

            end if ! Line presence

          end do ! b-b Transition

          ! For each b-f transition
          do itran=1,Atom(ia)%nphot

            ! The line is totally absent
            if (Atom(ia)%phot(itran)%if0.gt.if1.or. &
                Atom(ia)%phot(itran)%if1.lt.if0) then
              Atom(ia)%phot(itran)%absent = .True.
              Atom(ia)%phot(itran)%if0 = if1
              Atom(ia)%phot(itran)%if1 = if0-1
              Atom(ia)%phot(itran)%W0 = 0d0
              Atom(ia)%phot(itran)%W1 = 0d0

            ! The line is present
            else

              ! If the lower limit is out of range
              if (Atom(ia)%phot(itran)%if0.lt.if0) then
                Atom(ia)%phot(itran)%if0 = if0
                Atom(ia)%phot(itran)%W0 = Frec%W_freq(if0)
              end if

              ! If the upper limit is out of range
              if (Atom(ia)%phot(itran)%if1.gt.if1) then
                Atom(ia)%phot(itran)%if1 = if1
                Atom(ia)%phot(itran)%W1 = Frec%W_freq(if1)
              end if

              ! If there is only one frequency in this CPU
              if (Atom(ia)%phot(itran)%if0.eq. &
                  Atom(ia)%phot(itran)%if1) &
                Atom(ia)%phot(itran)%W1 = 0d0

              Frec%npfreq = Frec%npfreq + 1 + &
                            Atom(ia)%phot(itran)%if1 - &
                            Atom(ia)%phot(itran)%if0

            end if ! Line presence

          end do ! b-f Transition

        end do ! Atom

        ! Control limits
        if (Frec%ntfreqi.lt.1) Frec%ntfreqi = 1
        if (Frec%npfreq.lt.1) Frec%npfreq = 1
        if (Frec%ntfreq.lt.1) Frec%ntfreq = 1

        ! LTE lines
        do ia=1,nLTEl

          ! The line is totally absent
          if (Input%LTEline(ia)%absent) then

            ! Fake out of bound indexes
            Input%LTEline(ia)%if0 = if1
            Input%LTEline(ia)%if1 = if0-1

          ! The line is present
          else

            ! If the lower limit is out of range
            if (Input%LTEline(ia)%if0.lt.if0) &
              Input%LTEline(ia)%if0 = if0

            ! If the upper limit is out of range
            if (Input%LTEline(ia)%if1.gt.if1) &
              Input%LTEline(ia)%if1 = if1

          end if ! Line presence

        end do ! LTE lines

      end if ! Can resize

      ! Check if everything is fine
      call control

      return

      end subroutine frecresize

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resize some frequency dependent quantities, set CPU wise
      !! limits for the indexes where transitions are present and
      !! create output frequency axis for CLE.\n
      !!      Input(Input_class): Structure with settings data\n
      !!   Frec(Frequency_class): Structure with frequency data\n
      !!     Atom(Atom_class(:)): Structure with the atomic data\n
      !!         MPID(MPI_class): Structure with MPI data
      subroutine refitfrec(Input,Frec,Atom,MPID)

      ! I/O

      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Atom_class), dimension(:):: Atom
      type(Frequency_class), intent(inout):: Frec

      ! Local

      integer:: ia,iran,ifreq,jfreq,itran
      integer:: if0,if1,iif0,iif1,jf0,jf1,kf0,kf1
      integer, dimension(:), allocatable:: invmapping

      double precision, dimension(nfreq):: daux

      ! Store MPI limits in short variables
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)
      iif0 = MPID%iif0(pid)
      iif1 = MPID%iif1(pid)

      !
      ! Resize the true array
      !
      if(size(Frec%W_freq).ne.MPID%inf(pid))then

        ! Store weights into a temporal array
        daux = Frec%W_freq

        ! Resize weights
        deallocate(Frec%W_freq)
        allocate(Frec%W_freq(iif0:iif1))

        ! Recover the data
        Frec%W_freq = daux(iif0:iif1)

      end if ! Need to resize

      !
      ! Create indexing axis
      !

      ! Set omega3
      allocate(Frec%omega3(nfreq))
      Frec%omega3 = Frec%omega*Frec%omega*Frec%omega

      ! Allocate
      allocate(Frec%mapping(Input%lim_stk%nn))
      allocate(Frec%omega_ou(Input%lim_stk%nn))
      allocate(Frec%omega3_ou(Input%lim_stk%nn))

      ! Allocate "Master" line and photoionization limits
      allocate(Frec%Mlif0(pid:pid))
      allocate(Frec%Mlif1(pid:pid))
      allocate(Frec%Mpif0(pid:pid))
      allocate(Frec%Mpif1(pid:pid))
      ! And initialize
      Frec%Mlif0 = 10000000
      Frec%Mlif1 = -1
      Frec%Mpif0 = 10000000
      Frec%Mpif1 = -1
      Frec%lif0 = 10000000
      Frec%lif1 = -1
      Frec%pif0 = 10000000
      Frec%pif1 = -1


      ! For each atom
      do ia=1,nA

        ! Allocate input limits
        allocate(Atom(ia)%ilf0(Atom(ia)%ntran))
        allocate(Atom(ia)%ilf1(Atom(ia)%ntran))
        allocate(Atom(ia)%ipf0(Atom(ia)%nphot))
        allocate(Atom(ia)%ipf1(Atom(ia)%nphot))

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          ! Copy the limits
          Atom(ia)%ilf0(itran) = Atom(ia)%if0(itran)
          Atom(ia)%ilf1(itran) = Atom(ia)%if1(itran)

          ! Rearrange the limits taking into account the CPU limits

          ! The line is totally absent
          if (Atom(ia)%ilf1(itran).lt.iif0.or. &
              Atom(ia)%ilf0(itran).gt.iif1) then

            Atom(ia)%ilf0(itran) = iif1
            Atom(ia)%ilf1(itran) = iif0-1
            Atom(ia)%W0(itran) = 0d0
            Atom(ia)%W1(itran) = 0d0

          ! The line is present
          else

            ! If the lower limit is out of range
            if (Atom(ia)%ilf0(itran).lt.iif0) then
              Atom(ia)%ilf0(itran) = iif0
              Atom(ia)%W0(itran) = Frec%W_freq(iif0)
            end if

            ! If the upper limit is out of range
            if (Atom(ia)%ilf1(itran).gt.iif1) then
              Atom(ia)%ilf1(itran) = iif1
              Atom(ia)%W1(itran) = Frec%W_freq(iif1)
            end if

            ! If there is only one frequency in this CPU
            if (Atom(ia)%ilf0(itran).eq.Atom(ia)%ilf1(itran)) &
              Atom(ia)%W1(itran) = 0d0

          end if ! Line presence

          ! Handle the line range
          if (Atom(ia)%ilf0(itran).lt.Frec%Mlif0(pid)) &
            Frec%Mlif0(pid) = Atom(ia)%ilf0(itran)
          if (Atom(ia)%ilf1(itran).gt.Frec%Mlif1(pid)) &
            Frec%Mlif1(pid) = Atom(ia)%ilf1(itran)

        end do ! b-b Transition

        ! For each b-f transition
        do itran=1,Atom(ia)%nphot

          ! Copy
          Atom(ia)%ipf0(itran) = Atom(ia)%phot(itran)%if0
          Atom(ia)%ipf1(itran) = Atom(ia)%phot(itran)%if1

          ! The line is totally absent
          if (Atom(ia)%ipf0(itran).gt.iif1.or. &
              Atom(ia)%ipf1(itran).lt.iif0) then

            Atom(ia)%ipf0(itran) = iif1
            Atom(ia)%ipf1(itran) = iif0-1
            Atom(ia)%phot(itran)%W0 = 0d0
            Atom(ia)%phot(itran)%W1 = 0d0

          ! The line is present
          else

            ! If the lower limit is out of range
            if (Atom(ia)%ipf0(itran).lt.iif0) then
              Atom(ia)%ipf0(itran) = iif0
              Atom(ia)%phot(itran)%W0 = Frec%W_freq(if0)
            end if

            ! If the upper limit is out of range
            if (Atom(ia)%ipf1(itran).gt.iif1) then
              Atom(ia)%ipf1(itran) = iif1
              Atom(ia)%phot(itran)%W1 = Frec%W_freq(if1)
            end if

            ! If there is only one frequency in this CPU
            if (Atom(ia)%ipf0(itran).eq.Atom(ia)%ipf1(itran)) &
              Atom(ia)%phot(itran)%W1 = 0d0

            Frec%npfreq = Frec%npfreq + 1 + &
                          Atom(ia)%ipf1(itran) - &
                          Atom(ia)%ipf0(itran)

          end if ! Line presence

          ! Handle the photoionization range
          if (Atom(ia)%ipf0(itran).lt.Frec%Mpif0(pid)) &
            Frec%Mpif0(pid) = Atom(ia)%ipf0(itran)
          if (Atom(ia)%ipf1(itran).gt.Frec%Mpif1(pid)) &
            Frec%Mpif1(pid) = Atom(ia)%ipf1(itran)

        end do ! b-f Transition
      end do ! Atom

      ! Check limits
      if (Frec%Mlif0(pid).gt.Frec%Mlif1(pid)) then
        Frec%Mlif0(pid) = iif1
        Frec%Mlif1(pid) = iif0-1
      end if
      if (Frec%Mpif0(pid).gt.Frec%Mpif1(pid)) then
        Frec%Mpif0(pid) = iif1
        Frec%Mpif1(pid) = iif0-1
      end if

      ! If no limitations in output
      if (Input%lim_stk%nran.le.0) then

        ! Just copy the relevant data
        do ifreq=1,nfreq
          Frec%mapping(ifreq) = ifreq
        end do
        Frec%omega_ou = Frec%omega
        Frec%omega3_ou = Frec%omega3

        ! Copy in outputs
        do ia=1,nA

          ! b-b transitions
          do itran=1,Atom(ia)%ntran

            Atom(ia)%if0(itran) = Atom(ia)%ilf0(itran)
            Atom(ia)%if1(itran) = Atom(ia)%ilf1(itran)
            Atom(ia)%fflag(itran)%absent = &
                                      Atom(ia)%if1(itran).lt.if0.or. &
                                      Atom(ia)%if0(itran).gt.if1
          end do ! b-b

          ! b-f transitions
          do itran=1,Atom(ia)%nphot

            Atom(ia)%phot(itran)%if0 = Atom(ia)%ipf0(itran)
            Atom(ia)%phot(itran)%if1 = Atom(ia)%ipf1(itran)
            Atom(ia)%phot(itran)%absent = &
                                 Atom(ia)%phot(itran)%if1.lt.if0.or. &
                                 Atom(ia)%phot(itran)%if0.gt.if1

          end do ! b-f
        end do ! Atoms

        ! And copy limits
        Frec%lif0 = Frec%Mlif0(pid)
        Frec%lif1 = Frec%Mlif1(pid)
        Frec%pif0 = Frec%Mpif0(pid)
        Frec%pif1 = Frec%Mpif1(pid)

      ! If limiting the output axis
      else

        ! Allocate inverse mapping
        allocate(invmapping(nfreq))
        invmapping = -1

        !
        ! Initialize everything to absent
        !

        ! For each atom
        do ia=1,nA

          ! b-b transitions
          do itran=1,Atom(ia)%ntran
            Atom(ia)%fflag(itran)%absent = .True.
          end do

          ! b-f transitions
          do itran=1,Atom(ia)%nphot
            Atom(ia)%phot(itran)%absent = .True.
          end do
        end do

        ! Initialize rolling index
        jfreq = 0

        ! For each range
        do iran=1,Input%lim_stk%nran

          ! For each frequency in range
          do ifreq=Input%lim_stk%indx(1,iran), &
                   Input%lim_stk%indx(2,iran)

            ! Advance rolling index
            jfreq = jfreq + 1

            ! Copy correct frequency
            Frec%mapping(jfreq) = ifreq
            Frec%omega_ou(jfreq) = Frec%omega(ifreq)
            Frec%omega3_ou(jfreq) = Frec%omega3(ifreq)

            ! If below, skip
            if (jfreq.lt.MPID%if0(pid)) cycle
            ! If above, skip
            if (jfreq.gt.MPID%if1(pid)) cycle

            ! Inverse mapping if in this CPU
            invmapping(ifreq) = jfreq

          end do ! Frequency in range
        end do ! Ranges

        ! For each range (again)
        do iran=1,Input%lim_stk%nran

          ! Translate range limits
          iif0 = invmapping(Input%lim_stk%indx(1,iran))
          iif1 = invmapping(Input%lim_stk%indx(2,iran))


          ! If this range is out of mine, skip rest
          if (iif0.gt.MPID%if1(pid)) cycle
          if (iif1.lt.MPID%if0(pid)) cycle

          ! Get limits
          if0 = max(iif0,MPID%if0(pid))
          if1 = min(iif1,MPID%if1(pid))

          !
          ! Check lines present
          !

          ! For each atom
          do ia=1,nA

            ! b-b transition
            do itran=1,Atom(ia)%ntran

              ! If present, skip
              if (.not.Atom(ia)%fflag(itran)%absent) cycle

              ! Get limits
              jf0 = Atom(ia)%if0(itran)
              jf1 = Atom(ia)%if1(itran)

              ! Translate
              kf0 = invmapping(jf0)
              kf1 = invmapping(jf1)

              ! Find lower
              do while (kf0.lt.1)
                jf0 = jf0 + 1
                kf0 = invmapping(jf0)
                if (jf0.gt.nfreq) then
                  jf0 = -1
                  kf0 = -1
                  exit
                end if
              end do

              ! Bad?
              if (jf0.lt.1) cycle

              ! Find upper
              do while (kf1.lt.1)
                jf1 = jf1 - 1
                kf1 = invmapping(jf1)
                if (jf1.lt.1) then
                  jf1 = -1
                  kf1 = -1
                  exit
                end if
              end do

              ! Bad?
              if (jf1.lt.if0) cycle
              if (kf1.lt.if0) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.lt.jf0) cycle
              if (kf1.lt.kf0) cycle

              ! Save
              Atom(ia)%if0(itran) = kf0
              Atom(ia)%if1(itran) = kf1

              ! Is present here!
              Atom(ia)%fflag(itran)%absent = .False.

              ! Update totals
              if (Atom(ia)%if0(itran).lt.Frec%lif0) &
                Frec%lif0 = Atom(ia)%if0(itran)
              if (Atom(ia)%if1(itran).gt.Frec%lif1) &
                Frec%lif1 = Atom(ia)%if1(itran)

            end do ! b-b

            ! b-f transition
            do itran=1,Atom(ia)%nphot

              ! If present, skip
              if (.not.Atom(ia)%phot(itran)%absent) cycle

              ! Get limits
              jf0 = Atom(ia)%phot(itran)%if0
              jf1 = Atom(ia)%phot(itran)%if1

              ! Translate
              kf0 = invmapping(jf0)
              kf1 = invmapping(jf1)

              ! Find lower
              do while (kf0.lt.1)
                jf0 = jf0 + 1
                kf0 = invmapping(jf0)
                if (jf0.gt.nfreq) then
                  jf0 = -1
                  kf0 = -1
                  exit
                end if
              end do

              ! Bad?
              if (jf0.lt.1) cycle

              ! Find upper
              do while (kf1.lt.1)
                jf1 = jf1 - 1
                kf1 = invmapping(jf1)
                if (jf1.lt.1) then
                  jf1 = -1
                  kf1 = -1
                  exit
                end if
              end do

              ! Bad?
              if (jf1.lt.if0) cycle
              if (kf1.lt.if0) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.lt.jf0) cycle
              if (kf1.lt.kf0) cycle

              ! Translate limits
              Atom(ia)%phot(itran)%if0 = kf0
              Atom(ia)%phot(itran)%if1 = kf1

              ! Is present here!
              Atom(ia)%phot(itran)%absent = .False.

              ! Update totals
              if (Atom(ia)%phot(itran)%if0.lt.Frec%pif0) &
                Frec%pif0 = Atom(ia)%phot(itran)%if0
              if (Atom(ia)%phot(itran)%if1.gt.Frec%pif1) &
                Frec%pif1 = Atom(ia)%phot(itran)%if1

            end do ! b-f
          end do ! Atoms
        end do ! Ranges

      end if ! Limited output

      ! Allocate exu?
      if (Frec%pif1.ge.Frec%pif0) &
        allocate(Frec%exu(Frec%pif0:Frec%pif1,1))


      ! Check if everything is fine
      call control

      return

      end subroutine refitfrec

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks the non-coherent lower term approximation\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data
      subroutine check_nchlt(Atom,JKQ,Bfield)

      ! I/O
      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Bfield_class), intent(in):: Bfield
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                                                      intent(in):: JKQ

      ! Parameter
      double precision, parameter:: Bsat = 10d0

      ! Local

      logical:: skip

      integer:: ia,iz,jtran,itermu,itermf,iterml,itran,itran0,iJl
      integer:: itranmin,itranmax,litran
      integer, dimension(:), allocatable:: nmsg

      double precision:: Blu,gJ,rJ,S,rL,Bcrit,Bcrit0,efield


      ! Check that there is at least one height with non-zero
      ! magnetic field

      ! Initialize flag
      skip = .True.

      ! For each height
      do iz=Rz0,Rz1
        ! Check there is magnetic field
        if (Bfield%Bstrength(iz).gt.0d0) then
          skip = .False.
          exit
        end if
      end do ! Every height

      ! If no field, not necessary to check
      if (skip) return

      ! For each Atom
      do ia=1,nA

        ! Check that there is at least 1 PRD line in this atom
        skip = .True.

        do jtran=1,Atom(ia)%ntran
          if (Atom(ia)%lemiss2(jtran)) then
            skip = .False.
            exit
          end if
        end do

        ! There are no PRD lines for this atom
        if (skip) cycle

        !
        ! Allocate space for NCHLT approximation
        !

        ! Find limit in transition indexes
        itranmin = Atom(ia)%ntran + 1
        itranmax = -1

        ! Now go trough every pair of terms
        do itermu=2,Atom(ia)%nMulti
          do itermf=1,itermu-1

            ! Get transition
            jtran = Atom(ia)%irad(itermu,itermf)

            ! Check it exists and is PRD
            if (jtran.le.0) cycle
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! For each other lower term
            do iterml=1,itermu-1

              ! Check input transition
              itran = Atom(ia)%irad(itermu,iterml)

              ! That exists
              if (itran.le.0) cycle

              if (itranmin.gt.itran) itranmin = itran
              if (itranmax.lt.itran) itranmax = itran

            end do

          end do
        end do

        ! Allocate
        allocate(Atom(ia)%NCHLT(Rz0:Rz1,itranmin:itranmax))

        ! Initialize transition shift
        itran0 = Atom(ia)%tshift

        ! Allocate logical array to not repeat messages
        allocate(nmsg(Atom(ia)%nMulti))
        nmsg = 0

        ! Now go trough every pair of terms
        do itermu=2,Atom(ia)%nMulti
          do itermf=1,itermu-1

            ! Get transition
            jtran = Atom(ia)%irad(itermu,itermf)

            ! Check it exists and is PRD
            if (jtran.le.0) cycle
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! For each other lower term
            do iterml=1,itermu-1

              ! Check input transition
              litran = Atom(ia)%irad(itermu,iterml)

              ! That exists
              if (litran.le.0) cycle

              ! Get Blu
              Blu = Atom(ia)%Ecoeff(iterml,itermu)
              ! Get L and S
              rL = Atom(ia)%rLval(iterml)
              S = Atom(ia)%Sval(iterml)

              ! Get real itran
              itran = litran + itran0

              ! For each height
              do iz=Rz0,Rz1

                ! Critical field (Blu has factor 10^8)
                Bcrit0 = 1.137d1*Blu*dble(JKQ(0,0,itran,iz))*Bsat

                ! Initialize applicability
                skip = .True.

                ! Effective magnetic field
                efield = Bfield%Bstrength(iz)*sin(Bfield%Btheta(iz))

                ! For each level
                do iJl=1,Atom(ia)%nJ(iterml)

                  ! Get J
                  rJ = Atom(ia)%rJval(iJl,iterml)

                  ! Get Lande factor
                  if (Atom(ia)%ML) then
                    gJ = Atom(ia)%gL(iterml)
                  else
                    gJ = 1d0 + .5d0*(rJ*(rJ+1d0) + S*(S+1d0) - &
                                     rL*(rL+1d0))/rJ/(rJ+1d0)
                  end if

                  ! Quenching and Landé factor
                  if (gJ.gt.0d0) then
                    Bcrit = Bcrit0/gJ
                  else
                    Bcrit = 1d99
                  end if

                  ! Check field is big enough
                  if (efield.lt.Bcrit) then
                    skip = .False.
                    if (Bfield%Bstrength(iz).gt.0d0) &
                      nmsg(iterml) = nmsg(iterml) + 1
                    exit
                  end if

                end do ! Every lower level

                ! Store
                Atom(ia)%NCHLT(iz,litran) = skip

              end do ! Every height
            end do ! Every other lower term
          end do ! Every pair
        end do ! of terms

        ! Now, print the messages
        do iterml=1,Atom(ia)%nMulti

          ! If there are messages for this term
          if (nmsg(iterml).gt.0) then
            ! Multi-level
            if (Atom(ia)%ML) then
              write(umsg,'(A,A,1x,i3,1x,A,1x,i5,1x,A)') &
                ' # WARNING: The magnetic field is smaller than ', &
                'the saturation field of the level',iterml,'for ', &
                nmsg(iterml),'heights'
            ! Multi-term
            else
              write(umsg,'(A,A,1x,i3,1x,A,1x,i5,1x,A)') &
                ' # WARNING: The magnetic field is smaller than ', &
                'the saturation field of the term',iterml,'for ', &
                nmsg(iterml),'heights'
            end if ! Multi level or multi term
            call verbose
          end if ! Messages to print

        end do ! Messages to print

        ! Deallocate logical array of this atom
        deallocate(nmsg)

      end do ! Every atom

      return

      end subroutine check_nchlt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resets the redistribution function structure\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      subroutine cleanFrecandRed(Frec,Red,MPID)

      ! I/O

      type(Frequency_class):: Frec
      type(Red_class):: Red
      type(MPI_class):: MPID

      integer:: indx,jndx

      ! Remove counted memory
      MPID%RAM = MPID%RAM - MPID%WRAM
      MPID%WRAM = 0

      !
      ! Deallocate Frequency structure
      !

      ! Deallocate indexing
      if (allocated(Frec%indx)) deallocate(Frec%indx)

      ! If there is PRD data
      if (associated(Frec%dzao)) then

        ! For each index allocated
        do indx=1,Frec%ndzao

          ! If input transition data
          if (associated(Frec%dzao(indx)%trani)) then

            ! For each input transition
            do jndx=1,size(Frec%dzao(indx)%trani)

              ! Deallocate PRD info
              if (associated(Frec%dzao(indx)%trani(jndx)%index1)) then
                deallocate(Frec%dzao(indx)%trani(jndx)%index1)
                nullify(Frec%dzao(indx)%trani(jndx)%index1)
              end if
              if (associated(Frec%dzao(indx)%trani(jndx)%index2)) then
                deallocate(Frec%dzao(indx)%trani(jndx)%index2)
                nullify(Frec%dzao(indx)%trani(jndx)%index2)
              end if
              if (associated(Frec%dzao(indx)%trani(jndx)%dx)) then
                deallocate(Frec%dzao(indx)%trani(jndx)%dx)
                nullify(Frec%dzao(indx)%trani(jndx)%dx)
              end if
              if (allocated(Frec%dzao(indx)%trani(jndx)%omega)) &
                deallocate(Frec%dzao(indx)%trani(jndx)%omega)
              if (allocated(Frec%dzao(indx)%trani(jndx)%W_freq)) &
                deallocate(Frec%dzao(indx)%trani(jndx)%W_freq)
              if (allocated(Frec%dzao(indx)%trani(jndx)%mfreq)) &
                deallocate(Frec%dzao(indx)%trani(jndx)%mfreq)

            end do ! Input transitions

            ! Deallocate input transition and free the pointer
            deallocate(Frec%dzao(indx)%trani)
            nullify(Frec%dzao(indx)%trani)

          end if ! Input transition data
#ifdef _OPENMP
          if (allocated(Frec%dzao(indx)%oif0)) &
            deallocate(Frec%dzao(indx)%oif0, &
                       Frec%dzao(indx)%oif1)
#endif

        end do ! Indexes

        ! Deallocate and free first pointer
        deallocate(Frec%dzao)
        nullify(Frec%dzao)
        Frec%ndzao = 0

      end if ! There is PRD data

      ! Deallocate indexing
      if (allocated(Red%indx)) deallocate(Red%indx)

      ! If there is PRD data
      if (associated(Red%dzao)) then

        ! For each index allocated
        do indx=1,size(Red%dzao)

          ! If input transition data
          if (associated(Red%dzao(indx)%trani)) then

            ! For each input transition
            do jndx=1,size(Red%dzao(indx)%trani)

              ! Deallocate PRD info
              if (allocated(Red%dzao(indx)%trani(jndx)%iPPRD)) &
                deallocate(Red%dzao(indx)%trani(jndx)%iPPRD)
              if (allocated(Red%dzao(indx)%trani(jndx)%IWarr2)) &
                deallocate(Red%dzao(indx)%trani(jndx)%IWarr2)
              if (allocated(Red%dzao(indx)%trani(jndx)%PWarr2)) &
                deallocate(Red%dzao(indx)%trani(jndx)%PWarr2)

            end do ! Input transitions

            ! Deallocate input transition and free the pointer
            deallocate(Red%dzao(indx)%trani)
            nullify(Red%dzao(indx)%trani)

          end if ! Input transition data

        end do ! Indexes

        ! Deallocate and free first pointer
        deallocate(Red%dzao)
        nullify(Red%dzao)
        Red%ndzao = 0

      end if ! There is PRD data

      return

      end subroutine cleanFrecandRed

!#####################################################################
!#####################################################################
!#####################################################################

      end module omegabuild_mod
