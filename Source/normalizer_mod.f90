      !> Profile normalization
      module normalizer_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
!  Last version:
!     09/23/2024 V3.0.12
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/23/2024:   V3.0.12 - Changed the format in the writing of
!                             atmospheric models in ASCII to the
!                             recommended ratio between digits and
!                             decimals (TdPA)
!
!     05/14/2024:   V3.0.11 - Added a new argument to the
!                             normalization subroutine. If this
!                             argument is true, the routine returns
!                             after setting up the geometry indexing
!                             of the relevant structure (TdPA)
!
!     03/15/2024:   V3.0.10 - Save Voigt profile files in the output
!                             folder (TdPA)
!                           - Use array-stored indexes to run over
!                             the Normp pointer in the Atom class
!                             when freeing normalization (TdPA)
!                           - Use array-stored indexes to run over
!                             the prof pointer in the LTEline class
!                             when freeing profiles (TdPA)
!
!     10/16/2023:    V3.0.9 - Made LTElines allocatable to satisfy
!                             memory warnings (TdPA)
!                           - The bad normalization files are now
!                             written in the output directory (TdPA)
!                           - Use the lower and upper limits of
!                             Atom%Normp to run the loops (TdPA)
!                           - Check that %Norm is allocated before
!                             trying to deallocate (TdPA)
!
!     10/04/2023:    V3.0.8 - Bugfix: LTE lines were not working when
!                             not storing their Voigt profiles. Added
!                             a dummy prof structure of size 1 to
!                             point to in RT coeff whenever that is
!                             the case (TdPA)
!                           - The count of M components for LTE lines
!                             is now done elsewhere (TdPA)
!
!     09/26/2023:    V3.0.7 - Improved message for badnorm (TdPA)
!                           - Changed named of badnorm files to
!                             global ID, not local (TdPA)
!
!     09/25/2023:    V3.0.6 - Changed name of Voigt p. files (TdPA)
!
!     08/28/2023:    V3.0.5 - Changed check in allocation of
!                             normalization in atom class (TdPA)
!
!     08/07/2023:    V3.0.4 - Added normalization(), which contain the
!                             logic that used to be in the hanle
!                             routine, with the LTE lines added (TdPA)
!                           - Moved the geometry part of normalize and
!                             normalizeI to normalization(), avoiding
!                             repeating the same twice. This entails
!                             new arguments in both routines (TdPA)
!                           - Added getprof_LTE and getprofI_LTE
!                             routines (TdPA)
!
!     11/24/2022:    V3.0.3 - Set VRAM to False always in CLE (TdPA)
!
!     10/25/2022:    V3.0.2 - Implemented the restriction of the
!                             height axis (TdPA)
!                           - Added normalize_cle to just put
!                             everything to 1 in the CLE case (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!                           - funit is now global (TdPA)
!                           - MPID is no longer necessary in
!                             validatevoigtI (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Atmo%v has changed to Atmo%vx,%vy,
!                                and %vz.
!                              o Added returns for when we need to
!                                abort the run.
!                              o Changed the MPI communicator from
!                                MPI_COMM_WORLD to MPI_COMM_RT.
!                             (TdPA)
!                           - Fixed a typo in a warning message (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     02/09/2021:    V1.7.4 - Bugfix: Missing recall of quantum
!                             numbers when counting the expected
!                             memory to be allocated for profiles
!                             in Normalize in non-magnetic height
!                             nodes (TdPA)
!
!     02/04/2021:    V1.7.3 - Added the possibility of ignoring MIT
!                             for multi-term atoms with magnetic
!                             field. It completely ignores them if
!                             the magnetic resonance is off of the
!                             permitted resonances by at least
!                             MIT_Dw Doppler widths (TdPA)
!                           - Bugfix: There were some VIRAM in
!                             Normalize instead of VPRAM (TdPA)
!
!     09/15/2020:    V1.7.2 - Bugfix: In static magnetic cases there
!                             could be a problem with the dimensions
!                             of checkram (TdPA)
!
!     09/11/2020:    V1.7.1 - Added an extra variable to communicate
!                             if the RAM limit was reached inside the
!                             routine (TdPA)
!                           - Moved the indexing of directions to the
!                             beginning of the routine, because it
!                             was not called in certain situations
!                             for the LOS directions (TdPA)
!                           - Bugfix: There was a VIRAM instead of a
!                             VPRAM in normalize (TdPA)
!                           - Now the storage of the normalization
!                             values also counts for the RAM limit
!                             in the Voigt profiles (TdPA)
!                           - Bugfix: At some point in normalize,
!                             there were some loops in the magnetic
!                             sections where the upper boundaries
!                             were not correctly initialized (TdPA)
!                           - When deallocating profiles because they
!                             could not be allocated in another CPU,
!                             now the memory counter is updated (TdPA)
!
!     07/31/2020:    V1.7.0 - Bugfix: If we are storing Voigt profiles
!                             in RAM, but the CPUs run out of memory,
!                             this needs to be notified, because for
!                             a given height, direction and
!                             transition, it cannot happen that part
!                             of the profile is kept and others are
!                             not (TdPA)
!
!     06/05/2020:    V1.6.7 - size1 needs initialization after the
!                             change cutting the directional
!                             dimension (TdPA)
!
!     06/02/2020:    V1.6.6 - Avoid calling normalizations twice when
!                             is not needed (TdPA)
!                           - Limited the directional size of the
!                             normalization array when possible (TdPA)
!
!     03/18/2020:    V1.6.5 - Bugfix: With domain decomposition, the
!                             profile at the shared height node was
!                             only being computed in one of the CPU.
!                             Corrected it while avoiding double
!                             communication to the master (TdPA)
!
!     01/23/2020:    V1.6.4 - Added output to file when the
!                             normalization of profiles is not
!                             accurate. There is a logical variable
!                             that can switch it off upon compilation
!                             if it produces files too big (TdPA)
!                           - Bugfix: When validating a Voigt file
!                             for intensity, the frequency jumps were
!                             wrong (TdPA)
!
!     01/14/2020:    V1.6.3 - Bugfix: When adding file support I
!                             changed the meaning of id(3) and id(4)
!                             in normalize, but did not change the
!                             index used to check the magnetic field
!                             array element (TdPA)
!
!     12/10/2019:    V1.6.2 - Added digits in format of bad
!                             normalization warning (TdPA)
!
!     11/19/2019:    V1.6.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - Changed the memory initialization (TdPA)
!
!     11/12/2019:    V1.6.0 - Added the possibility to store the Voigt
!                             profiles in binary files (TdPA)
!                           - Added routines to check existing Voigt
!                             profile files (TdPA)
!                           - Prepared the normalization routines to
!                             work either for quadrature directions
!                             or LOS directions, depending on an input
!                             argument (TdPA)
!                           - Changed the dimension order of the
!                             Atom%Normp array (TdPA)
!                           - Bugfix: When computing the RAM taken by
!                             profiles, there were considered double
!                             instead of double complex (TdPA)
!
!     10/02/2019:    V1.5.3 - Bugfix: Fixed wrong formatting for
!                             warnings in normalizer (TdPA)
!
!     07/19/2019:    V1.5.2 - Added a check for bad normalizations.
!                             The threshold is a parameter called
!                             BADNORM (TdPA)
!
!     04/16/2019:    V1.5.1 - Implemented back the normal broadcasting
!                             from mpi libraries (TdPA)
!
!     02/20/2019:    V1.5.0 - New verbosity (TdPA)
!                           - Using specific TINY variables (TdPA)
!                           - Changed cpulimit for old way of
!                             communication (TdPA)
!                           - Checks success of communications when
!                             possible (TdPA)
!
!     09/05/2018:    V1.4.1 - Using the old way of communication for
!                             smaller amounts of CPU, and updated
!                             also the intensity subroutine (TdPA)
!
!     09/04/2018:    V1.4.0 - Had to change the communication
!                             algorithm because it started to fail
!                             in supercomputers due to excess of
!                             messages (TdPA)
!
!     08/09/2018:    V1.3.2 - Bugfix: Forgot to apply the full change
!                             in V1.3.1 to the normalizeI (TdPA)
!
!     08/08/2018:    V1.3.1 - Bugfix: Master does not need to manage
!                             normalizations or normalize profiles
!                             when they are stored. Removed that for
!                             the Master (TdPA)
!
!     08/04/2018:    V1.3.0 - Stores Voigt profiles if requested
!                             in the input (TdPA)
!
!     10/03/2017:    V1.2.0 - Implemented non-magnetic case (TdPA)
!                           - Forgot to identify the terms in one non
!                             magnetic part of normalize (TdPA)
!
!     09/28/2017:    V1.1.0 - Implemented send_tree algorithm to
!                             avoid the terrible scaling of the native
!                             mpi_bcast (TdPA)
!
!     07/19/2017:    V1.0.4 - Now using the proper weights in the
!                             boundaries (TdPA)
!
!     06/16/2017:    V1.0.3 - Using the transition limits instead
!                             of logical flags to check presence
!                             of a transition (TdPA)
!
!     05/04/2017:    V1.0.2 - Eliminated status calls (TdPA)
!                           - Changed domain decomposition isend by
!                             send (TdPA)
!                           - Avoiding to allocate Geom%i_geom
!                             multiple times (TdPA)
!
!     05/03/2017:    V1.0.1 - Made sure that the isend with domain
!                             decomposition waits properly every
!                             time to not mix buffers and overlap
!                             orders (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!  normalization:
!    manages the normalization and profile storage
!
!  normalize:
!    normalize profiles for term-term transitions taking into account
!  the magnetic component
!
!  getprof_LTE:
!    pre-calculate profiles for LTE lines accounting for potential
!  magnetic splitting
!
!  validatevoigt
!    validates an existing file with Voigt profiles
!
!  normalizeI:
!    normalize profiles for FS transitions not taking into account
!  the magnetic component
!
!  getprofI_LTE:
!    pre-calculate profiles for LTE lines
!
!  validatevoigtI
!    validates an existing file with Voigt (intensity) profiles
!
!  writebadbound
!    outputs into a file data about bad normalization of Voigt
!  profiles
!
!  normalize_cle
!    Initialize normalization structures and set them to 1 for the
!  CLE case
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use funnj_mod
      use parameters_mod , only : cZero, TINYB, TINYJS, TINYN, &
                                  sqrt3, PI, c, BADNORM , TINYO, B2LK
      use profile_mod
      use setmpi_mod
      use types_mod

      ! To output the normalization problems
      logical:: obadnorm = .True.

      ! CPU threshold to ask the master for permission
      integer:: cpulimit = 16

      ! Maximum minimum Doppler distance for MIT not to be sampled
      double precision:: MIT_Dw = 10d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manages the normalization of profiles\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      lines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bstrength(dfloat(:)): Magnetic field strength\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!        Input(Input_class): Structure with settings data\n
      !!     Frec(Frequency_class): Structure with frequency data
      !!            rlimw(logical): Write RAM limit message\n
      !!              lit(logical): If the code had to iterate\n
      !!     polarization(logical): Normalizing profiles for
      !!                            polarization problem\n
      !!              los(logical): Normalizing for LOS formal
      !!                            solutions\n
      !!        only_geom(logical): If only entering to set-up
      !!                            the geometry
      subroutine normalization(Atom,lines,Atmo,Bstrength,Geom,Frec, &
                               Input,Flgsg,MPID,rlimw,lit, &
                               polarization,los,only_geom)

      ! IO
      type(Atom_class), dimension(:):: Atom
      type(LTEline_class), dimension(:), allocatable:: lines
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(inout):: Geom
      type(Frequency_class), intent(in):: Frec
      type(Input_class):: Input
      type(Fctsg_class):: Flgsg
      type(MPI_class):: MPID
      logical, intent(in):: polarization, los, lit, only_geom
      logical, intent(inout):: rlimw
      double precision, dimension(:), intent(in):: Bstrength

      ! Local
      logical:: norm, exists, ofram

      integer:: ia, ios, jdir, ith, iph, njdir
      integer, dimension(:), allocatable:: ithv,iphv


      ! Polarization
      if (polarization) then

        !
        ! Initialize geometry
        !

        ! LOS
        if (los) then

          ! Directions
          njdir = Geom%nPhLOS*Geom%nThLOS

          ! Allocate indexing
          if (allocated(Geom%i_geom)) deallocate(Geom%i_geom)
          allocate(Geom%i_geom(Geom%nPhLOS,Geom%nThLOS))
          Geom%i_geom = 0

          ! Index of polar direction
          allocate(ithv(njdir))
          ! Index of azimuthal direction
          allocate(iphv(njdir))

          ! De-index the directions
          jdir = 0
          do ith=1,Geom%nThLOS
            do iph=1,Geom%nPhLOS
              jdir = jdir + 1
              Geom%i_geom(iph,ith) = jdir
              ithv(jdir) = ith
              iphv(jdir) = iph
            end do
          end do

        ! Quadrature
        else

          ! Directions
          njdir = Geom%nPh*Geom%nTh

          ! Allocate indexing
          if (allocated(Geom%i_geom)) deallocate(Geom%i_geom)
          allocate(Geom%i_geom(Geom%nPh,Geom%nTh))
          Geom%i_geom = 0

          ! Index of polar direction
          allocate(ithv(njdir))
          ! Index of azimuthal direction
          allocate(iphv(njdir))

          ! De-index the directions
          jdir = 0
          do ith=1,Geom%nTh
            do iph=1,Geom%nPh
              jdir = jdir + 1
              Geom%i_geom(iph,ith) = jdir
              ithv(jdir) = ith
              iphv(jdir) = iph
            end do
          end do

          ! If Voigt profiles in file, initialize variable
          if (vpfil) ios = 0

          ! Reset RAM usage
          MPID%VRAM = 0d0

        end if

        ! Only geometry?
        if (only_geom) return

        ! For each atom
        do ia=1,nA

          ! Initialie to true
          norm = .True.

          ! If using a file with Voigt profiles, check there is
          ! an existing one
          if (vpfil) then

            ! If LOS
            if (los) then

              ! Set filename
              Atom(ia)%vfile = trim(Input%Folder)//'voigt-P-L-'// &
                               trim(Atom(ia)%file_label)

            ! If quadrature
            else

              ! Set filename
              Atom(ia)%vfile = trim(Input%folder)//'voigt-P-G-'// &
                               trim(Atom(ia)%file_label)
            end if

            inquire(file=trim(Atom(ia)%vfile), exist=exists)

            ! If there is a file with such name, check it is valid
            if (exists) then

              call validatevoigt(Atom(ia),Bstrength,MPID,Frec, &
                                 njdir,Flgsg,norm)
            end if
          end if

          ! If we have to normalize
          if (norm) then

            if (gpid.eq.0) then
              if (los) then
                umsg = ' - Normalizing LOS profiles for '// &
                       Atom(ia)%Element
              else
                umsg = ' - Normalizing quadrature profiles for '// &
                       Atom(ia)%Element
              end if
              call verbose
            end if

            call normalize(Atom(ia),Atmo,Bstrength,Geom,MPID,Frec, &
                           Flgsg,Input%folder,njdir,ithv,iphv, &
                           Input%MIT_input,lit,ofram,los)

            if (laborted) return

            if (gpid.eq.0) then
              umsg = ' - Profiles normalized '
              call verbose
            end if

            ! If storing Voigt profiles
            if (VPRAM) then

             ! If CPU went above ram
              if (ofram.and.rlimw) then

                write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                    ' reached the limit of profile allocations.'
                call verbose

                rlimw = .False.

              end if

              call MPI_BARRIER(MPI_COMM_RT, ierr)

            end if ! Storing Voigt profiles
          end if ! Reading file with profiles

          ! If from file, update variable
          if (.not.los.and.vpfil) ios = max(ios,Atom(ia)%Mncom)

        end do ! Atoms

        ! If from file, store maximum components in first atom
        if (.not.los) Atom(1)%Mncom = ios

        ! If no LTE lines, leave
        if (Input%nLTE.lt.1) return

        ! If not storing
        if (.not.LVPRAM) then

          ! Allocate dummy and leave
          allocate(lines(1)%prof(1,1))
          lines(1)%prof(1,1)%VRAM = .False.
          return

        end if

        !
        ! Normalize LTE lines
        !
        do ia=1,Input%nLTE

          ! Get profiles
          call getprof_LTE(lines(ia),Atmo,Bstrength,Geom,MPID,Frec, &
                           njdir,ithv,iphv,lit,ofram,los)

          if (laborted) return

        end do ! LTE lines

        ! If storing Voigt profiles
        if (LVPRAM) then

          ! If CPU went above ram
          if (ofram.and.rlimw) then

            write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                ' reached the limit of profile allocations.'
            call verbose

            rlimw = .False.

          end if

          call MPI_BARRIER(MPI_COMM_RT, ierr)

        end if ! Storing Voigt profiles

      ! Intensity
      else

        !
        ! Initialize geometry
        !

        ! LOS
        if (los) then

          ! Total number of directions
          njdir = Geom%nPhLOS*Geom%nThLOS

          ! Allocate indexing
          if (allocated(Geom%i_geom)) deallocate(Geom%i_geom)
          allocate(Geom%i_geom(Geom%nPhLOS,Geom%nThLOS))
          Geom%i_geom = 0

          ! Index of polar direction
          allocate(ithv(njdir))
          ! Index of azimuthal direction
          allocate(iphv(njdir))

          ! De-index the directions
          jdir = 0
          do ith=1,Geom%nThLOS
            do iph=1,Geom%nPhLOS
              jdir = jdir + 1
              Geom%i_geom(iph,ith) = jdir
              ithv(jdir) = ith
              iphv(jdir) = iph
            end do
          end do

        ! Quadrature
        else

          ! Total number of directions
          njdir = Geom%nPh*Geom%nTh

          ! Allocate indexing
          if (allocated(Geom%i_geom)) deallocate(Geom%i_geom)
          allocate(Geom%i_geom(Geom%nPh,Geom%nTh))
          Geom%i_geom = 0

          ! Index of polar direction
          allocate(ithv(njdir))
          ! Index of azimuthal direction
          allocate(iphv(njdir))

          ! De-index the directions
          jdir = 0
          do ith=1,Geom%nTh
            do iph=1,Geom%nPh
              jdir = jdir + 1
              Geom%i_geom(iph,ith) = jdir
              ithv(jdir) = ith
              iphv(jdir) = iph
            end do
          end do

          ! Reset RAM usage
          MPID%VRAM = 0d0

        end if

        ! Only geometry?
        if (only_geom) return

        !
        ! Normalize FS profiles
        !

        ! For each Atom
        do ia=1,nA

          ! Initialie to true
          norm = .True.

          ! If using a file with Voigt profiles, check there is
          ! an existing one
          if (vifil) then

            ! If LOS normalization
            if (los) then

              ! Set filename
              Atom(ia)%vfile = trim(Input%folder)//'voigt-I-L-'// &
                               trim(Atom(ia)%file_label)
            ! Iterations
            else

              ! Set filename
              Atom(ia)%vfile = trim(Input%folder)//'voigt-I-G-'// &
                               trim(Atom(ia)%file_label)

            end if ! Type of solution

            ! Check if file exist
            inquire(file=trim(Atom(ia)%vfile), exist=exists)

            ! If there is a file with such name, check it is valid
            if (exists) then

              call validatevoigtI(Atom(ia),Frec,njdir,norm)

            end if
          end if

          ! If we have to normalize
          if (norm) then

            if (gpid.eq.0) then
              if (los) then
                umsg = ' - Normalizing LOS profiles for '// &
                       Atom(ia)%Element
              else
                umsg = ' - Normalizing quadrature profiles for '// &
                       Atom(ia)%Element
              end if
              call verbose
            end if

            ! Call normalization
            call normalizeI(Atom(ia),Atmo,Geom,MPID,Frec, &
                            Input%folder,njdir,ithv,iphv,lit, &
                            ofram,los)

            ! Control
            if (laborted) return

            if (gpid.eq.0) then
              umsg = ' - Profiles normalized '
              call verbose
            end if

            ! If storing Voigt profiles
            if (VIRAM) then

              ! If CPU went above ram
              if (ofram.and.rlimw) then

                write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                    ' reached the limit of profile allocations.'
                call verbose

                rlimw = .False.

              end if

              call MPI_BARRIER(MPI_COMM_RT, ierr)

            end if ! Storing Voigt profiles
          end if ! Reading file with profiles

        end do ! Atoms

        ! If no LTE lines, leave
        if (Input%nLTE.lt.1) return

        ! If not storing
        if (.not.LVIRAM) then

          ! Allocate dummy and leave
          allocate(lines(1)%prof(1,1))
          lines(1)%prof(1,1)%VRAM = .False.
          return

        end if

        !
        ! Normalize LTE lines
        !

        ! For each line
        do ia=1,Input%nLTE

          ! Skip if absent
          if (lines(ia)%absent) cycle

          ! Get profiles
          call getprofI_LTE(lines(ia),Atmo,Geom,MPID,Frec,njdir, &
                            ithv,iphv,lit,ofram,los)

          ! Control
          if (laborted) return

          ! If storing Voigt profiles
          if (LVIRAM) then

            ! If CPU went above ram
            if (ofram.and.rlimw) then

              write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                  ' reached the limit of profile allocations.'
              call verbose

              rlimw = .False.

            end if

            call MPI_BARRIER(MPI_COMM_RT, ierr)

          end if ! Storing Voigt profiles

        end do ! Atoms

      end if ! Intensity

      end subroutine normalization

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes normalization factors for all the Voigt profiles\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bstrength(dfloat(:)): Magnetic field strength\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!     Frec(Frequency_class): Structure with frequency data
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!    folder(character(500)): Output folder path\n
      !!            njdir(integer): Number of directions\n
      !!          ithv(integer(:)): Indexing of polar directions\n
      !!          ithv(integer(:)): Indexing of azimuth directions\n
      !!                            signs\n
      !!              MIT(integer): Input about MIT\n
      !!               lp(logical): If doing formal solution in this
      !!                            run\n
      !!            ofram(logical): Indicates if out of RAM\n
      !!              LOS(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine normalize(Atom,Atmo,Bstrength,Geom,MPID,Frec, &
                           Flgsg,folder,njdir,ithv,iphv,MIT,lp, &
                           ofram,LOS)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(inout):: Geom
      type(Frequency_class), intent(in):: Frec
      type(Fctsg_class):: Flgsg
      type(MPI_class):: MPID
      character(len=500), intent(in):: folder
      logical, intent(in):: lp, LOS
      logical, intent(out):: ofram
      integer, intent(in):: MIT, njdir
      integer, dimension(:), intent(in):: ithv,iphv
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      character(len=4):: record

      logical:: extracomm,field,checkMIT,found
      integer, dimension(:), allocatable:: outofbound

      integer:: i,i1,jtran,itermf,itermu
      integer:: idir,jdir,iz,ii,jj,nodir,nrdir
      integer:: ith,iph,ifreq,if0,if1,rif0,rif1,nfreqt
      integer:: ncom,ncomNB
      integer:: nMu,nMf,icom
      integer:: iMl,iMu,iMf
      integer:: iL,iU,mF,nt
      integer:: nL,nU,znjdir
      integer:: ios,istep,lcheckram
      integer, dimension(:,:,:), allocatable:: checkram

      double precision:: d1,W0,W1
      double precision:: rLu,rLf,S
      double precision:: rJumax,rJfmax,rJu,rJf
      double precision:: rMu,rMf,f62
      double precision:: el,eu
      double precision:: au,af,auf,Dfreqw
      double precision:: DwT,Dw,vfac,vfacw
      double precision:: ct,st,cc,sc
      double precision:: d2,d3,d1T,d2T,loffset,loffsetin
      double precision, dimension(:), allocatable:: nut

      complex(kind=8):: prof

      ! Buffers, counter, and sizes for MPI
      logical, dimension(:,:,:), allocatable:: tosend
      integer:: finished
      integer, dimension(5):: id, id1_b
      integer, dimension(0:nproc-1):: nbf1
      integer, dimension(:,:,:), allocatable:: size1
      double precision, dimension(:),allocatable:: buff1

      ! MPI offset type
      integer(kind=MPI_OFFSET_KIND):: offset


      ! Routine name
      urou = 'normalize'

      ! Check if there is magnetic field
      field = .False.

      ! Initialize
      ofram = .False.

      ! Check magnetic field strength
      do iz=1,nz
        if (Bstrength(iz).ge.TINYB) then
          field = .True.
          exit
        end if
      end do

      ! If there is magnetic field, the atom is multi-term
      ! and the MIT are switched off
      if (field.and..not.Atom%ML.and.MIT.lt.0) then

        ! We need to check later for MIT
        checkMIT = .True.

        !
        ! Get maximum number of components
        nt = 0

        ! For each transition
        do jtran=1,Atom%ntran

          ! Identify the terms
          itermf = Atom%fst(jtran)%iterml
          itermu = Atom%fst(jtran)%itermu

          ! Update nt
          if (Atom%nJ(itermf)*Atom%nJ(itermu).gt.nt) &
            nt = Atom%nJ(itermf)*Atom%nJ(itermu)

        end do ! Transitions

        ! Allocate vector for frequencies
        allocate(nut(nt))

      ! No need to check for MIT
      else

        ! No need to check
        checkMIT = .False.

      end if

      ! If LOS and not dynamic, if already iterated, maybe no need to
      ! repeat the normalization
      if (LOS.and..not.dyn.and.lp) then

        ! If no magnetic, copy the file if using it and get out
        if (.not.field) then

          ! If using files, just copy to new name
          if (vpfil.and.pid.eq.0) then

            ! Orginal
            open(200,file='voigt-P-G-'//trim(Atom%file_label), &
                 access='stream',status='old',action='read', &
                 iostat=ios)

            ! New
            open(300,file=trim(Atom%vfile),access='stream', &
                 status='unknown',action='write',iostat=ios)

            ! Copy records until finished
            do while (.True.)

              ! Read characters
              read(200, end=3000) record
              ! Write the same integer
              write(300) record

              cycle
3000          exit

            end do

            close(200)
            close(300)

          end if ! Copy vfile

          ! Everyone control and return
          call control
          return

        end if ! Non-magnetic
      end if ! LOS, not dynamic, previously normalized


      !
      ! Initialize comm flag
      !
      if (nproc.gt.cpulimit) then
        extracomm = .True.
      else
        extracomm = .False.
      end if

      ! Get real size of direction dimension
      if (dyn) then
        nodir = njdir
        nrdir = njdir
      else
        nodir = 1
        if (field) then
          nrdir = njdir
        else
          nrdir = 1
        end if
      end if


      !
      ! Allocate vector for norm (part I)
      !

      ! Check Normp is not allocated
      if (associated(Atom%Normp)) then
        do jdir=lbound(Atom%Normp,3),ubound(Atom%Normp,3)
          do iz=lbound(Atom%Normp,2),ubound(Atom%Normp,2)
            do jtran=lbound(Atom%Normp,1),ubound(Atom%Normp,1)
              if (allocated(Atom%Normp(jtran,iz,jdir)%prof)) &
                deallocate(Atom%Normp(jtran,iz,jdir)%prof)
              if (allocated(Atom%Normp(jtran,iz,jdir)%Norm)) &
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)
            end do
          end do
        end do
        deallocate(Atom%Normp)
        nullify(Atom%Normp)
      end if

      ! Structure with the norm for each component
      allocate(Atom%Normp(Atom%ntran,Rz0:Rz1,nrdir))
      ! Allocate sizes for subcomponents for each transition
      if (.not.allocated(Atom%nL)) allocate(Atom%nL(Atom%ntran))
      if (.not.allocated(Atom%nU)) allocate(Atom%nU(Atom%ntran))
      if (.not.allocated(Atom%nMl)) allocate(Atom%nMl(Atom%ntran))
      if (.not.allocated(Atom%nMu)) allocate(Atom%nMu(Atom%ntran))
      ! Allocate size for MPI
      allocate(size1(Atom%ntran,Rz0:Rz1,nrdir))
      size1 = 0
      ! Allocate bool for MPI
      allocate(tosend(Atom%ntran,Rz0:Rz1,nrdir))
      ! Initialize tosend
      tosend = .False.

      ! Allocate the check for the master
      if (pid.eq.0.and.MPID%mpi) then
        if (field) then
          allocate(checkram(Atom%ntran,Rz0:Rz1,nrdir))
        else
          allocate(checkram(Atom%ntran,Rz0:Rz1,nodir))
        end if
        if (VPRAM) checkram = 1
      end if


      !
      ! Count indexes
      !

      ! Reset index
      nbf1 = 0

      ! For each transition
      do jtran=1,Atom%ntran

        ! Identify the terms
        itermf = Atom%fst(jtran)%iterml
        itermu = Atom%fst(jtran)%itermu

        ! Get atomic quantities
        S = Atom%Sval(itermu)

        rLu = Atom%rLval(itermu)
        rJumax = rLu+S
        nMu = nint(2d0*rJumax+1d0)

        rLf = Atom%rLval(itermf)
        rJfmax = rLf + S
        nMf = nint(2d0*rJfmax+1d0)

        ! Fill two of the sizes arrays
        Atom%nMu(jtran) = nMu
        Atom%nMl(jtran) = nMf

        ! Initialize counters
        nL = 0
        nU = 0

        ! For each Mu
        do iMu=1,nMu

          rMu = -rJumax + dble(iMu-1)

          ! For each mu_u
          do iU=1,Atom%nblk(iMu,itermu)

            ! Update counter
            if (iU.gt.nU) nu = iMu

            ! For each Mf
            do iMf=1,nMf

              rMf = -rJfmax + dble(iMf-1)

              ! If not allowed, skip
              if (nint(abs(rMu-rMf)).gt.1) cycle

              ! Sum over mu_f
              do mF=1,Atom%nblk(iMf,itermf)

                ! Update counter
                if (mF.gt.nL) nL = mf

              end do ! iL
            end do ! Ml
          end do ! iU
        end do ! Mu

        ! Update variables
        Atom%nU(jtran) = nU
        Atom%nL(jtran) = nL

      end do ! Transitions


      !
      ! If using a file, prepare for it
      !
      if (vpfil) then

        ! Deallocate sizes
        if (allocated(Atom%zsize)) deallocate(Atom%zsize)
        if (allocated(Atom%dsize)) deallocate(Atom%dsize)
        if (allocated(Atom%tsize)) deallocate(Atom%tsize)
        if (allocated(Atom%tBsize)) deallocate(Atom%tBsize)
        if (allocated(Atom%f0size)) deallocate(Atom%f0size)
        if (allocated(Atom%f1size)) deallocate(Atom%f1size)
        if (allocated(Atom%i_Vind)) deallocate(Atom%i_Vind)

        ! Allocate sizes and indexes
        allocate(Atom%zsize(nz))
        allocate(Atom%dsize(nrdir))
        allocate(Atom%tsize(Atom%ntran))
        allocate(Atom%tBsize(Atom%ntran))
        allocate(Atom%f0size(Atom%ntran))
        allocate(Atom%f1size(Atom%ntran))
        allocate(Atom%i_Vind(Atom%ntran))

        ! Initialize total blocks of sizes and transition index
        d1T = 0d0
        d2T = 0d0
        Atom%zsize(1) = 0d0
        Atom%dsize(1) = 0d0
        Atom%tsize(1) = 0d0
        Atom%tBsize(1) = 0d0
        Atom%Mncom = 0
        ncom = 0
        ncomNB = 0

        ! For each transition
        do jtran=1,Atom%ntran

          ! Get real limits
          if0 = Atom%rif0(jtran)
          if1 = Atom%rif1(jtran)

          ! Get terms
          itermf = Atom%fst(jtran)%iterml
          itermu = Atom%fst(jtran)%itermu

          ! Get atomic quantities
          S = Atom%Sval(itermu)

          rLu = Atom%rLval(itermu)
          rJumax = rLu+S
          nMu = nint(2d0*rJumax+1d0)

          rLf = Atom%rLval(itermf)
          rJfmax = rLf + S
          nMf = nint(2d0*rJfmax+1d0)

          ! Allocate magnetic indexing
          allocate(Atom%i_Vind(jtran)% &
                        ind(maxval(Atom%nblk(1:nMf,itermf)),nMf, &
                            maxval(Atom%nblk(1:nMu,itermu)),nMu))

          ! Allocate non-magnetic indexing
          allocate(Atom%i_Vind(jtran)% &
                        indNB(Atom%nJ(itermf),Atom%nJ(itermu)))

          ! No magnetic field
          d1 = 0d0

          ! Reset index
          i = 0

          ! sum over Ju
          do iU=1,Atom%nJ(itermu)

            ! Get Ju
            rJu = Atom%rJval(iU,itermu)

            ! sum over Jl
            do mF=1,Atom%nJ(itermf)

              ! Get Jl
              rJf = Atom%rJval(mF,itermf)

              ! 6-j
              f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

              if (abs(f62).lt.TINYJS) cycle

              d1 = d1 + dble(if1 - if0 + 1)

              ! Store and advance index
              i = i + 1
              Atom%i_Vind(jtran)%indNB(mF,iU) = i

            end do ! Final levels
          end do ! Upper levels

          ! Store number of components
          Atom%i_Vind(jtran)%ncomNB = i
          if (i.gt.ncomNB) ncomNB = i

          ! Yes magnetic field
          d2 = 0d0

          ! Reset index
          i = 0

          ! sum over Mu
          do iMu=1,nMu

            rMu = -rJumax + dble(iMu-1)

            ! sum over mu_u
            do iU=1,Atom%nblk(iMu,itermu)

              ! sum over Ml
              do iMf=1,nMf

                rMf = -rJfmax + dble(iMf-1)

                if (nint(abs(rMu-rMf)).gt.1) cycle

                ! sum over mu_l
                do mF=1,Atom%nblk(iMf,itermf)

                  d2 = d2 + dble(if1 - if0 + 1)

                  ! Store and advance index
                  i = i + 1
                  Atom%i_Vind(jtran)%ind(mF,iMf,iU,iMu) = i

                end do ! mu_l
              end do ! Ml
            end do ! mu_u
          end do ! Mu

          ! Store number of components
          Atom%i_Vind(jtran)%ncom = i
          if (i.gt.ncomNB) ncom = i

          ! Add to the total size
          d1T = d1T + d1
          d2T = d2T + d2

          ! If not the last, add size to next
          if (jtran.lt.Atom%ntran) then
            Atom%tsize(jtran+1) = Atom%tsize(jtran) + d1
            Atom%tBsize(jtran+1) = Atom%tBsize(jtran) + d2
          end if

          ! deallocate indexes of magnetic components
          if (maxval(Bstrength).lt.TINYB) &
            deallocate(Atom%i_Vind(jtran)%ind)

          ! deallocate indexes of non magnetic components
          if (minval(Bstrength).ge.TINYB) &
            deallocate(Atom%i_Vind(jtran)%indNB)

        end do ! Transitions

        ! If no magnetic field
        if (maxval(Bstrength).lt.TINYB) ncom = 0
        ! If only magnetic field
        if (minval(Bstrength).ge.TINYB) ncomNB = 0
        ! Get maximum
        Atom%Mncom = max(ncom,ncomNB)

        ! Initialize
        d3 = 0d0

        ! For each height
        do iz=1,nz

          ! If not magnetic
          if (Bstrength(iz).lt.TINYB) then
            if (iz.lt.nz) Atom%zsize(iz+1) = Atom%zsize(iz) + d1T
            d3 = d3 + d1T
          ! If magnetic
          else
            if (iz.lt.nz) Atom%zsize(iz+1) = Atom%zsize(iz) + d2T
            d3 = d3 + d2T
          end if

        end do ! Heights

        ! For each direction
        do idir=1,nrdir-1
          Atom%dsize(idir+1) = Atom%dsize(idir) + d3
        end do

        Atom%dsize = Atom%dsize*16d0
        Atom%zsize = Atom%zsize*16d0
        Atom%tsize = Atom%tsize*16d0
        Atom%tBsize = Atom%tBsize*16d0

        ! Size of header
        Atom%hvifil = 4*4 + & ! Dimension integers
                      8*nrdir + & ! Directions sizes
                      8*nz + & ! Height sizes
                      4*2*Atom%ntran + & ! Line limits
                      8*2*Atom%ntran + & ! Line sizes
                      8*nfreq ! Frequencies

        !
        ! Only Master writes
        !
        if (pid.eq.0) then

          ! Open files
          open(200, file=trim(Atom%vfile), status='unknown', &
               iostat=ios, err=1000, access='stream', &
               action='write', form='unformatted')

          ! Write dimensions
          write(200, err=1100) nrdir
          write(200, err=1100) nz
          write(200, err=1100) nfreq
          write(200, err=1100) Atom%ntran
          write(200, err=1100) Atom%dsize
          write(200, err=1100) Atom%zsize

          ! For each transition
          do jtran=1,Atom%ntran

            write(200, err=1100) Atom%rif0(jtran)
            write(200, err=1100) Atom%rif1(jtran)
            write(200, err=1100) Atom%tsize(jtran)
            write(200, err=1100) Atom%tBsize(jtran)

          end do ! Transitions

          ! Write frequency
          write(200, err=1100) Frec%omega

          ! And close this
          close(200, err=1100)

        end if ! Master

        ! The master does not need the indexes
        if (pid.eq.0.and.MPID%mpi) deallocate(Atom%i_Vind)

        ! Control
        call control

      end if ! Voigt file


      !
      ! Allocate the norm array
      !

      ! For each height
      do iz=Rz0,Rz1

        ! No magnetic field
        if (Bstrength(iz).lt.TINYB) then

          ! For each direction
          do jdir=1,nodir

            ! For each transition
            do jtran=1,Atom%ntran

              ! Identify the terms
              itermf = Atom%fst(jtran)%iterml
              itermu = Atom%fst(jtran)%itermu

              ! Get atomic quantities
              S = Atom%Sval(itermu)
              rLu = Atom%rLval(itermu)
              rLf = Atom%rLval(itermf)

              ! Allocate
              allocate(Atom%Normp(jtran,iz,jdir)%Norm( &
                       Atom%nJ(itermf), Atom%nJ(itermu), 1, 1))
              Atom%Normp(jtran,iz,jdir)%Norm = 0d0

              ! Determine size of this array
              size1(jtran,iz,jdir) = Atom%nJ(itermf)*Atom%nJ(itermu)

              ! Skip unless slave
              if (MPID%mpi.and.pid.eq.0) cycle

              ! Count components

              jj = 0

              ! sum over Ju
              do iU=1,Atom%nJ(itermu)

                ! Get Ju
                rJu = Atom%rJval(iU,itermu)

                ! sum over Jl
                do mF=1,Atom%nJ(itermf)

                  ! Get Jl
                  rJf = Atom%rJval(mF,itermf)

                  ! 6-j
                  f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

                  if (abs(f62).lt.TINYJS) cycle

                  jj = jj + 1

                end do ! Final levels
              end do ! Upper levels

              ! Allocate profile itself if storing and it is present
              if (VPRAM.and..not.Atom%fflag(jtran)%absent) then

                ! Prediction
                d1 = 16d-6*dble((Atom%if1(jtran) - &
                                 Atom%if0(jtran) + 1)*jj)

                ! If no more space
                if (floor(MPID%RAM+d1).gt.RLIM) then

                  ! No stored
                  Atom%Normp(jtran,iz,jdir)%VRAM = .False.
                  ofram = .True.

                  ! Add normalization to RAM
                  d1 = 8d-6*dble(jj)
                  MPID%RAM = MPID%RAM + d1
                  MPID%VRAM = MPID%VRAM + d1

                ! If there is space
                else

                  ! Storing
                  Atom%Normp(jtran,iz,jdir)%VRAM = .True.

                  ! Allocate
                  allocate(Atom%Normp(jtran,iz,jdir)%prof( &
                             Atom%nJ(itermf),Atom%nJ(itermu),1,1))

                  ! sum over Ju
                  do iU=1,Atom%nJ(itermu)

                    ! Get Ju
                    rJu = Atom%rJval(iU,itermu)

                    ! sum over Jl
                    do mF=1,Atom%nJ(itermf)

                      ! Get Jl
                      rJf = Atom%rJval(mF,itermf)

                      ! 6-j
                      f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

                      if (abs(f62).lt.TINYJS) cycle

                      ! Allocate
                      allocate(Atom%Normp(jtran,iz,jdir)%prof( &
                                  mF,iU,1,1)%cp( &
                               Atom%if0(jtran):Atom%if1(jtran)))

                    end do ! Final levels
                  end do ! Upper levels

                  ! Update RAM
                  MPID%RAM = MPID%RAM + d1
                  MPID%VRAM = MPID%VRAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Atom%Normp(jtran,iz,jdir)%VRAM = .False.

                ! Add normalization to RAM
                d1 = 8d-6*dble(jj)
                MPID%RAM = MPID%RAM + d1
                MPID%VRAM = MPID%VRAM + d1

              end if ! Storing

            end do ! transitions
          end do ! directions

        ! Yes magnetic field
        else

          ! For each direction
          do jdir=1,nrdir

            ! For each transition
            do jtran=1,Atom%ntran

              ! Allocate
              allocate(Atom%Normp(jtran,iz,jdir)%Norm( &
                       Atom%nL(jtran), Atom%nU(jtran), &
                       Atom%nMl(jtran), Atom%nMu(jtran)))
              Atom%Normp(jtran,iz,jdir)%Norm = 0d0

              ! Determine size of this array
              size1(jtran,iz,jdir) = Atom%nL(jtran)*Atom%nU(jtran)* &
                                     Atom%nMl(jtran)*Atom%nMu(jtran)

              ! Skip if not slave
              if (MPID%mpi.and.pid.eq.0) cycle

              ! Identify the terms
              itermf = Atom%fst(jtran)%iterml
              itermu = Atom%fst(jtran)%itermu

              ! Get atomic quantities
              S = Atom%Sval(itermu)

              rLu = Atom%rLval(itermu)
              rJumax = rLu+S
              nMu = nint(2d0*rJumax+1d0)

              rLf = Atom%rLval(itermf)
              rJfmax = rLf + S
              nMf = nint(2d0*rJfmax+1d0)

              ! Get indexes
              if0 = Atom%if0(jtran)
              if1 = Atom%if1(jtran)

              ! Count components

              ! Initialize
              jj = 0

              ! sum over Mu
              do iMu=1,nMu

                rMu = -rJumax + dble(iMu-1)

                ! sum over mu_u
                do iU=1,Atom%nblk(iMu,itermu)

                  ! sum over Ml
                  do iMf=1,nMf

                    rMf = -rJfmax + dble(iMf-1)

                    if (nint(abs(rMu-rMf)).gt.1) cycle

                    ! sum over mu_l
                    do mF=1,Atom%nblk(iMf,itermf)

                      jj = jj + 1

                    end do ! mu_l
                  end do ! Ml
                end do ! mu_u
              end do ! Mu

              ! Allocate profile itself if storing and it is present
              if (VPRAM.and..not.Atom%fflag(jtran)%absent) then

                ! Prediction
                d1 = 16d-6*dble((Atom%if1(jtran) - &
                                 Atom%if0(jtran) + 1)*jj)

                ! If no more space
                if (floor(MPID%RAM+d1).gt.RLIM) then

                  ! No stored
                  Atom%Normp(jtran,iz,jdir)%VRAM = .False.
                  ofram = .True.
                  ! Add normalization to RAM
                  d1 = 8d-6*dble(jj)
                  MPID%RAM = MPID%RAM + d1
                  MPID%VRAM = MPID%VRAM + d1

                ! If there is space
                else

                  ! Storing
                  Atom%Normp(jtran,iz,jdir)%VRAM = .True.

                  ! Allocate
                  allocate(Atom%Normp(jtran,iz,jdir)%prof( &
                           Atom%nL(jtran), Atom%nU(jtran), &
                           Atom%nMl(jtran), Atom%nMu(jtran)))

                  ! sum over Mu
                  do iMu=1,nMu

                    rMu = -rJumax + dble(iMu-1)

                    ! sum over mu_u
                    do iU=1,Atom%nblk(iMu,itermu)

                      ! sum over Ml
                      do iMf=1,nMf

                        rMf = -rJfmax + dble(iMf-1)

                        if (nint(abs(rMu-rMf)).gt.1) cycle

                        ! sum over mu_l
                        do mF=1,Atom%nblk(iMf,itermf)

                          ! Allocate
                          allocate(Atom%Normp(jtran,iz,jdir)% &
                                   prof(mF,iU,iMf,iMu)%cp(if0:if1))

                        end do ! mu_l
                      end do ! Ml
                    end do ! mu_u
                  end do ! Mu

                  ! Update RAM
                  MPID%RAM = MPID%RAM + d1
                  MPID%VRAM = MPID%VRAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Atom%Normp(jtran,iz,jdir)%VRAM = .False.

                ! Add normalization to RAM
                d1 = 8d-6*dble(jj)
                MPID%RAM = MPID%RAM + d1
                MPID%VRAM = MPID%VRAM + d1

              end if ! Storing

            end do ! transitions
          end do ! directions

        end if ! No magnetic field

      end do ! heights


      ! Check the maximum size to transfer
      ios = maxval(size1)

      ! Control
      call control
      if (laborted) return

      ! Gather the maximum size that each processor is holding
      do while (.True.)
        call MPI_ALLGATHER(ios,1,MPI_INTEGER,nbf1(0),1, &
                           MPI_INTEGER,MPI_COMM_RT,ierr)
        if (ierr.eq.0) exit
      end do

      ! Allocate buffers
      allocate(buff1(nbf1(pid)))

      ! Control
      call control
      if (laborted) return


      !
      ! MASTER
      !

      if (MPID%mpi.and.pid.eq.0) then

        !
        ! Initialize finished
        !
        finished = 1

        !
        ! Calculate normalization
        !

        do while (finished.lt.nproc)

          ! Receive the informative package with indexes
          if (extracomm) then
            do while (.True.)
              call MPI_recv(id(1),5,MPI_INTEGER, &
                            MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do
          else
            call MPI_recv(id(1),5,MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
          end if

          ! If ending signal
          if (id(1).lt.0) then
            finished = finished + 1
            cycle
          end if

          if (extracomm) then
            do while (.True.)
              call MPI_SEND(id(1),1,MPI_INTEGER,id(1),id(1), &
                            MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do
          end if

          ! Receive the buffer with the integral
          if (extracomm) then
            do while (.True.)
              call MPI_recv(buff1(1), nbf1(id(1)), &
                            MPI_DOUBLE_PRECISION, id(1), &
                            id(1), MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do
          else
            call MPI_recv(buff1(1), nbf1(id(1)), &
                          MPI_DOUBLE_PRECISION, id(1), &
                          id(1), MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
          end if

          ! Store checkram
          if (VPRAM.and.id(5).lt.0) &
            checkram(id(2),id(3),id(4)) = -1

          ! Reset index
          jj = 0

          ! No magnetic field
          if (Bstrength(id(3)).lt.TINYB) then

            ! Identify the terms
            do i=1,Atom%nMulti-1
              do i1=i+1,Atom%nMulti
                if (Atom%irad(i,i1).eq.id(2)) then
                  itermf = i
                  itermu = i1
                end if
              end do
            end do

            ! Run over the components
            do iU=1,Atom%nJ(itermu)
              do iL=1,Atom%nJ(itermf)

                ! Advance indexes
                jj = jj + 1

                ! Accumulate the sub-integrals
                Atom%Normp(id(2),id(3),id(4))% &
                     Norm(iL,iU,1,1) = &
                      Atom%Normp(id(2),id(3),id(4))% &
                         Norm(iL,iU,1,1) + buff1(jj)

              end do ! iJl
            end do ! iJu

          ! Yes magnetic field
          else

            ! Run over the components
            do iMu=1,Atom%nMu(id(2))
              do iMl=1,Atom%nMl(id(2))
                do iU=1,Atom%nU(id(2))
                  do iL=1,Atom%nL(id(2))

                    ! Advance indexes
                    jj = jj + 1

                    ! Accumulate the sub-integrals
                    Atom%Normp(id(2),id(3),id(4))% &
                         Norm(iL,iU,iMl,iMu) = &
                      Atom%Normp(id(2),id(3),id(4))% &
                         Norm(iL,iU,iMl,iMu) + buff1(jj)

                  end do ! iL
                end do ! iU
              end do ! iMl
            end do ! iMu

          end if ! Magnetic field

        end do ! Communication to do


      !
      ! SLAVE OR SINGLE PROCESSOR
      !

      else

        !
        ! Calculate normalization
        !

        ! For each height
        do iz=Rz0,Rz1

          ! Thermal part of the Doppler width
          DwT = Atom%cDopp*sqrt(Atmo%T(iz))

          ! Select number of directions
          if (Bstrength(iz).ge.TINYB) then
            znjdir = nrdir
          else
            znjdir = nodir
          end if

          ! For each direction
          do jdir=1,znjdir

            ! Recover the indexes
            ith = ithv(jdir)
            iph = iphv(jdir)

            ! If emergent
            if (LOS) then

              ct = Geom%L_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = cos(Geom%L_phi(iph))
              sc = sin(Geom%L_phi(iph))

            else

              ct = Geom%V_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = Geom%v_mux(iph)
              sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

            end if

            ! Calculate Doppler shift factor

            vfac = 1d0

            if (dyn) &
              vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                           atmo%vz(iz)*ct

            ! For each transition
            do jtran=1,Atom%ntran

              ! If this process have frequencies for this line
              if (.not.Atom%fflag(jtran)%absent) then

                ! Flag to send
                tosend(jtran,iz,jdir) = .True.

                ! Output Doppler width
                Dw = Atom%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

                ! Find the term indexes for this transition
                do i=1,Atom%nMulti-1
                  do i1=i+1,Atom%nMulti
                    if (Atom%irad(i,i1).eq.jtran) then
                      itermf = i
                      itermu = i1
                    end if
                  end do
                end do

                ! Get contributions to damping parameter
                au = Atom%damp(itermu,iz)/Dw
                af = Atom%damp(itermf,iz)/Dw
                auf = Atom%ldamp(jtran,iz)/Dw

                ! Get atomic quantities
                S = Atom%Sval(itermu)

                rLu = Atom%rLval(itermu)
                rJumax = rLu+S
                nMu = nint(2d0*rJumax+1d0)

                rLf = Atom%rLval(itermf)
                rJfmax = rLf + S
                nMf = nint(2d0*rJfmax+1d0)

                ! Get indexes
                if0 = Atom%if0(jtran)
                if1 = Atom%if1(jtran)
                ! Get weights
                W0 = Atom%W0(jtran)
                W1 = Atom%W1(jtran)

                !
                ! Proper normalization
                !

                ! No magnetic field
                if (Bstrength(iz).lt.TINYB) then

                  ! sum over Ju
                  do iU=1,Atom%nJ(itermu)

                    eu = Atom%FSfreq(iU,itermu)/Dw

                    ! Get Ju
                    rJu = Atom%rJval(iU,itermu)

                    ! sum over Jl
                    do mF=1,Atom%nJ(itermf)

                      ! Get Jl
                      rJf = Atom%rJval(mF,itermf)

                      el = Atom%FSfreq(mF,itermf)/Dw

                      ! 6-j
                      f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

                      if (abs(f62).lt.TINYJS) cycle

                      !
                      ! Calculate profile
                      !

                      ! Common quantities
                      Dfreqw = eu - el
                      vfacw = vfac/Dw
                      d1 = 1d-5/(sqrt(PI)*Dw)

                      ! Boundaries

                      ! Lower
                      call voigt(Dfreqw - Frec%omega(if0)*vfacw, &
                                 au+af+auf,prof)

                      Atom%Normp(jtran,iz,jdir)% &
                           Norm(mF,iU,1,1) = dble(prof)*(W0*d1)

                      if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                        Atom%Normp(jtran,iz,jdir)%prof(mF,iU,1,1)% &
                        cp(if0) = prof

                      ! For each frequency
                      do ifreq=if0+1,if1-1

                        call voigt(Dfreqw - Frec%omega(ifreq)*vfacw, &
                                   au+af+auf,prof)

                        ! Add to the integral
                        Atom%Normp(jtran,iz,jdir)% &
                             Norm(mF,iU,1,1) = &
                             Atom%Normp(jtran,iz,jdir)% &
                                  Norm(mF,iU,1,1) + &
                                  dble(prof)*(Frec%W_freq(ifreq)*d1)

                        if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                          Atom%Normp(jtran,iz,jdir)%prof(mF,iU,1,1)% &
                          cp(ifreq) = prof

                      end do ! frequencies

                      ! Upper
                      call voigt(Dfreqw - Frec%omega(if1)*vfacw, &
                                 au+af+auf,prof)

                      Atom%Normp(jtran,iz,jdir)% &
                           Norm(mF,iU,1,1) = &
                           Atom%Normp(jtran,iz,jdir)% &
                                Norm(mF,iU,1,1) + dble(prof)*(W1*d1)

                      if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                        Atom%Normp(jtran,iz,jdir)%prof(mF,iU,1,1)% &
                        cp(if1) = prof

                    end do ! Jf
                  end do ! Ju

                ! Yes magnetic field
                else

                  ! If we need to check for MIT, store
                  ! resonances
                  if (checkMIT) then

                    ! Initialize size and index
                    nt = 0
                    ii = 0

                    ! For every J combination
                    do i=1,Atom%nJ(itermf)
                      do i1=1,Atom%nJ(itermu)

                        ! Skip if forbidden
                        if (Atom%fst(jtran)%irad(i1,i).lt.1) cycle

                        ! Add resonance
                        ii = ii + 1
                        nut(ii) = (Atom%FSfreq(i1,itermu) - &
                                   Atom%FSfreq(i,itermf))/Dw

                      end do ! Upper J level
                    end do ! Lower J level

                    ! Store size
                    nt = ii

                  end if ! If checking for MIT

                  ! sum over Mu
                  do iMu=1,nMu

                    rMu = -rJumax + dble(iMu-1)

                    ! sum over mu_u
                    do iU=1,Atom%nblk(iMu,itermu)

                      eu = Atom%eval(iU,iMu,itermu,iz)/Dw

                      ! sum over Ml
                      do iMf=1,nMf

                        rMf = -rJfmax + dble(iMf-1)

                        if (nint(abs(rMu-rMf)).gt.1) cycle

                        ! sum over mu_l
                        do mF=1,Atom%nblk(iMf,itermf)

                          el = Atom%eval(mF,iMf,itermf,iz)/Dw


                          !
                          ! Calculate profile
                          !

                          ! Common quantities
                          Dfreqw = eu - el + Atom%Dfreq(jtran)/Dw
                          vfacw = vfac/Dw
                          d1 = 1d-5/(sqrt(PI)*Dw)

                          !
                          ! If checking for MIT
                          if (checkMIT) then

                            ! Initialize found flag
                            found = .False.

                            ! Compare with all resonances
                            do ii=1,nt

                              ! Check if close
                              if (abs(Dfreqw - nut(ii)).lt.MIT_Dw) then
                                found = .True.
                                exit
                              end if

                            end do

                            ! If not found, we set to zero and skip
                            if (.not.found) then

                              Atom%Normp(jtran,iz,jdir)% &
                                         Norm(mF,iU,iMf,iMu) = 0d0
                              if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                                Atom%Normp(jtran,iz,jdir)% &
                                     prof(mF,iU,iMf,iMu)% &
                                     cp(if0:if1) = cZero

                              cycle

                            end if

                          end if ! Checking for MIT

                          ! Boundaries

                          ! Lower
                          call voigt(Dfreqw - Frec%omega(if0)*vfacw, &
                                     au+af+auf,prof)
                          Atom%Normp(jtran,iz,jdir)% &
                                               Norm(mF,iU,iMf,iMu) = &
                                                    dble(prof)*(W0*d1)

                          if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                            Atom%Normp(jtran,iz,jdir)% &
                                 prof(mF,iU,iMf,iMu)% &
                                 cp(if0) = prof

                          ! For each frequency
                          do ifreq=if0+1,if1-1

                            call voigt(Dfreqw - &
                                       Frec%omega(ifreq)*vfacw, &
                                       au+af+auf,prof)

                            ! Add to the integral
                            Atom%Normp(jtran,iz,jdir)% &
                                               Norm(mF,iU,iMf,iMu) = &
                               Atom%Normp(jtran,iz,jdir)% &
                                               Norm(mF,iU,iMf,iMu) + &
                               dble(prof)*(Frec%W_freq(ifreq)*d1)

                            if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                              Atom%Normp(jtran,iz,jdir)% &
                                   prof(mF,iU,iMf,iMu)% &
                                   cp(ifreq) = prof

                          end do ! frequencies

                          ! Upper
                          call voigt(Dfreqw - Frec%omega(if1)*vfacw, &
                                     au+af+auf,prof)
                          Atom%Normp(jtran,iz,jdir)% &
                                               Norm(mF,iU,iMf,iMu) = &
                               Atom%Normp(jtran,iz,jdir)% &
                                               Norm(mF,iU,iMf,iMu) + &
                                                    dble(prof)*(W1*d1)

                          if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                            Atom%Normp(jtran,iz,jdir)% &
                                 prof(mF,iU,iMf,iMu)% &
                                 cp(if1) = prof

                        end do ! iL
                      end do ! Ml
                    end do ! iU
                  end do ! Mu

                end if ! Magnetic field presence

              end if ! Line presence in processor test

              if (MPID%mpi.and.tosend(jtran,iz,jdir)) then

                !
                ! Share data with master
                !

                ! Check last send was received
                if (.not.extracomm) then
                  call MPI_WAIT(MPID%request1,MPI_STATUS_IGNORE, &
                                ierr)
                  call MPI_WAIT(MPID%request2,MPI_STATUS_IGNORE, &
                                ierr)
                end if

                ! Send the indexing data
                if (Atom%Normp(jtran,iz,jdir)%VRAM) then
                  id1_b = (/ pid, jtran, iz, jdir, 1 /)
                else
                  id1_b = (/ pid, jtran, iz, jdir,-1 /)
                end if
                if (extracomm) then
                  do while (.True.)
                    call MPI_SEND(id1_b(1),5,MPI_INTEGER,0,0, &
                                  MPI_COMM_RT,ierr)
                    if (ierr.eq.0) exit
                  end do
                else
                  call MPI_ISEND(id1_b(1),5,MPI_INTEGER,0,0, &
                                 MPI_COMM_RT,MPID%request1,ierr)
                end if

                ! Reorder the normalization into the buffer
                buff1(1:size1(jtran,iz,jdir)) = &
                           reshape(Atom%Normp(jtran,iz,jdir)%Norm, &
                                  (/ size1(jtran,iz,jdir) /))

                if (extracomm) then
                  do while (.True.)
                    call MPI_recv(id1_b(1),1,MPI_INTEGER, &
                                  0, pid, MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)
                    if (ierr.eq.0) exit
                  end do
                end if

                ! Send the actual normalization values
                if (extracomm) then
                  do while (.True.)
                    call MPI_SEND(buff1(1),nbf1(pid), &
                                  MPI_DOUBLE_PRECISION, 0, pid, &
                                  MPI_COMM_RT, ierr)
                    if (ierr.eq.0) exit
                  end do
                else
                  call MPI_ISEND(buff1(1),nbf1(pid), &
                                 MPI_DOUBLE_PRECISION, 0, pid, &
                                 MPI_COMM_RT, &
                                 MPID%request2, ierr)
                end if

              end if ! MPI

            end do ! output transition
          end do ! output direction
        end do ! height

        ! If MPI send finished signal
        if (MPID%mpi) then

          ! Check last send was received
          if (.not.extracomm) &
            call MPI_WAIT(MPID%request1,MPI_STATUS_IGNORE,ierr)

          ! Send the indexing data
          id = (/ -1, -1, -1, -1, -1 /)
          if (extracomm) then
            do while (.True.)
              call MPI_SEND(id(1),5,MPI_INTEGER,0,0, &
                             MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do
          else
            call MPI_ISEND(id(1),5,MPI_INTEGER,0,0, &
                            MPI_COMM_RT,MPID%request1,ierr)
          end if

        end if

      end if ! MPI

      ! If MPI
      if (MPID%mpi) then

        !
        ! Broadcast the results
        !

        ! Slaves allocate checkram
        if (pid.gt.0.and.VPRAM) then
          ! Allocate checkram
          allocate(checkram(Atom%ntran,Rz0:Rz1,nodir))
        end if

        ! For each height
        do iz=Rz0,Rz1

          ! Select number of directions
          if (Bstrength(iz).ge.TINYB) then
            znjdir = nrdir
          else
            znjdir = nodir
          end if

          ! For each direction
          do jdir=1,znjdir

            ! For each transition
            do jtran=1,Atom%ntran

              ! Alternative bcast
              if (MPID%altbcast) then

                ! If not master, receive first
                if (pid.ne.0) then

                  ! Receive Norm
                  call MPI_RECV(Atom%Normp(jtran,iz,jdir)% &
                                Norm(1,1,1,1), &
                                size1(jtran,iz,jdir), &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                end if ! No master

                ! For each send
                do istep=1,MPID%nsend

                  ! Send Norm
                  call MPI_ISEND(Atom%Normp(jtran,iz,jdir)% &
                                 Norm(1,1,1,1), &
                                 size1(jtran,iz,jdir), &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,1), ierr)

                end do ! Sends

                ! For each slave
                do istep=1,MPID%nsend

                  ! Wait for everyone to receive the radiation data
                  ! continuing
                  call MPI_WAIT(MPID%requestA(istep,1), &
                                MPI_STATUS_IGNORE,ierr)

                end do ! Sends

              ! Normal bcast
              else

                ! Share Norm
                call MPI_BCAST(Atom%Normp(jtran,iz,jdir)% &
                                    Norm(1,1,1,1), &
                               size1(jtran,iz,jdir), &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

              end if ! Type of bcast

            end do ! transitions
          end do ! directions
        end do ! heights

        ! If storing in RAM
        if (VPRAM) then

          ! Send checkram
          lcheckram = Atom%ntran*Rnz*nodir

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive Norm
              call MPI_RECV(checkram(1,Rz0,1),lcheckram, &
                            MPI_INTEGER,MPID%recv, pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No master

            ! For each send
            do istep=1,MPID%nsend

              ! Send Norm
              call MPI_ISEND(checkram(1,Rz0,1),lcheckram, &
                             MPI_INTEGER,MPID%lsend(istep), &
                             MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,1), ierr)

            end do ! Sends

            ! For each slave
            do istep=1,MPID%nsend

              ! Wait for everyone to receive the radiation data
              ! continuing
              call MPI_WAIT(MPID%requestA(istep,1), &
                            MPI_STATUS_IGNORE,ierr)

            end do ! Sends

          ! Normal bcast
          else

            ! Share Norm
            call MPI_BCAST(checkram(1,Rz0,1),lcheckram, &
                           MPI_INTEGER, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast

          ! Slaves deal with it
          if (pid.gt.0) then

            ! For each direction
            do jdir=1,nodir

              ! For each height
              do iz=Rz0,Rz1

                ! For each transition
                do jtran=1,Atom%ntran

                  ! If was saving but cannot
                  if (Atom%Normp(jtran,iz,jdir)%VRAM.and. &
                      checkram(jtran,iz,jdir).lt.0) then

                    ! Not storing now
                    Atom%Normp(jtran,iz,jdir)%VRAM = .False.

                    ! Identify the terms
                    itermf = Atom%fst(jtran)%iterml
                    itermu = Atom%fst(jtran)%itermu

                    ! No magnetic field
                    if (Bstrength(iz).lt.TINYB) then

                      ! sum over Ju
                      do iU=1,Atom%nJ(itermu)
                        ! sum over Jl
                        do mF=1,Atom%nJ(itermf)

                          ! If stored already, remove
                          if (allocated(Atom%Normp(jtran,iz,jdir)% &
                                        prof(mF,iU,1,1)%cp)) then
                            deallocate(Atom%Normp(jtran,iz,jdir)% &
                                        prof(mF,iU,1,1)%cp)
                            ! Update RAM
                            d1 = 16d-6*dble(Atom%if1(jtran) - &
                                            Atom%if0(jtran) + 1)
                            MPID%RAM = MPID%RAM - d1
                            MPID%VRAM = MPID%VRAM - d1
                            d1 = 8d-6
                            MPID%RAM = MPID%RAM + d1
                            MPID%VRAM = MPID%VRAM + d1

                          end if ! Allocated

                        end do ! Jf
                      end do ! Ju

                    ! Yes magnetic field
                    else

                      rLu = Atom%rLval(itermu)
                      rJumax = rLu+S
                      nMu = nint(2d0*rJumax+1d0)

                      rLf = Atom%rLval(itermf)
                      rJfmax = rLf + S
                      nMf = nint(2d0*rJfmax+1d0)

                      ! sum over Mu
                      do iMu=1,nMu
                        ! sum over mu_u
                        do iU=1,Atom%nblk(iMu,itermu)
                          ! sum over Ml
                          do iMf=1,nMf
                            if (nint(abs(rMu-rMf)).gt.1) cycle
                            ! sum over mu_l
                            do mF=1,Atom%nblk(iMf,itermf)

                              ! If stored already, remove
                              if (allocated(Atom% &
                                            Normp(jtran,iz,jdir)% &
                                            prof(mF,iU,iMf,iMu)% &
                                            cp)) then
                                deallocate(Atom% &
                                           Normp(jtran,iz,jdir)% &
                                           prof(mF,iU,iMf,iMu)%cp)
                                ! Update RAM
                                d1 = 16d-6*dble(Atom%if1(jtran) - &
                                                Atom%if0(jtran) + 1)
                                MPID%RAM = MPID%RAM - d1
                                MPID%VRAM = MPID%VRAM - d1
                                d1 = 8d-6
                                MPID%RAM = MPID%RAM + d1
                                MPID%VRAM = MPID%VRAM + d1

                              end if ! Allocated

                            end do ! iL
                          end do ! Ml
                        end do ! iU
                      end do ! Mu

                    end if ! Magnetic field presence

                    ! And deallocate prof
                    deallocate(Atom%Normp(jtran,iz,jdir)%prof)

                  end if ! Was storing and not now

                end do ! transition
              end do ! height
            end do ! direction

          end if ! Slaves
        end if ! If saving in RAM
      end if ! MPI

      !
      ! Calculate multiplicative factor (instead of division factor)
      !

      ! If MPI, master does not need this
      if (MPID%mpi.and.pid.eq.0) then

        do jdir=1,size(Atom%Normp,3)
          do iz=Rz0,Rz1
            do jtran=1,Atom%ntran
              if (allocated(Atom%Normp(jtran,iz,jdir)%Norm)) &
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)
            end do
          end do
        end do
        deallocate(Atom%Normp)
        nullify(Atom%Normp)

      else

        ! Allocate and initialize out of bound counter
        allocate(outofbound(Atom%ntran))
        outofbound = 0

        ! For each height
        do iz=Rz0,Rz1

          ! Select number of directions
          if (Bstrength(iz).ge.TINYB) then
            znjdir = nrdir
          else
            znjdir = nodir
          end if

          ! For each direction
          do jdir=1,znjdir

            ! For each transition
            do jtran=1,Atom%ntran

              ! No magnetic field
              if (Bstrength(iz).lt.TINYB) then

                ! If storing
                if (Atom%Normp(jtran,iz,jdir)%VRAM) then

                  ! Identify the terms
                  itermf = Atom%fst(jtran)%iterml
                  itermu = Atom%fst(jtran)%itermu

                  ! Run over the J
                  do iU=1,Atom%nJ(itermu)
                    do iL=1,Atom%nJ(itermf)

                      ! Easier to write variable
                      d1 = Atom%Normp(jtran,iz,jdir)%Norm(iL,iU,1,1)

                      ! If the norm is not zero
                      if (d1.gt.TINYN) then

                        ! Check close to 1
                        if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then
                          outofbound(jtran) = outofbound(jtran) + 1
                          if (obadnorm) &
                            call writebadbound(folder, &
                                               Atom%Element,iz,jdir, &
                                               .True.,.False.,jtran, &
                                               iU,iL,1,1,d1)
                        end if

                        Atom%Normp(jtran,iz,jdir)%prof(iL,iU,1,1)%cp &
                           =  dcmplx(dble(Atom%Normp(jtran,iz,jdir)% &
                                      prof(iL,iU,1,1)%cp)/d1, &
                                 dimag(Atom%Normp(jtran,iz,jdir)% &
                                      prof(iL,iU,1,1)%cp))
                      end if

                    end do ! iL
                  end do ! iU

                  ! And deallocate the norm because it is not needed
                  deallocate(Atom%Normp(jtran,iz,jdir)%Norm)

                ! Not storing
                else

                  ! Identify the terms
                  itermf = Atom%fst(jtran)%iterml
                  itermu = Atom%fst(jtran)%itermu

                  ! Run over the J
                  do iU=1,Atom%nJ(itermu)
                    do iL=1,Atom%nJ(itermf)

                      ! Easier to write variable
                      d1 = Atom%Normp(jtran,iz,jdir)%Norm(iL,iU,1,1)

                      ! If the norm is not zero
                      if (d1.gt.TINYN) then

                        ! Check close to 1
                        if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then
                          outofbound(jtran) = outofbound(jtran) + 1
                          if (obadnorm) &
                            call writebadbound(folder, &
                                               Atom%Element,iz,jdir, &
                                              .True.,.False.,jtran, &
                                              iU,iL,1,1,d1)
                        end if

                        Atom%Normp(jtran,iz,jdir)%Norm(iL,iU,1,1) = &
                                                               1d0/d1
                      end if

                    end do ! iL
                  end do ! iU

                end if ! Storing

              ! Yes magnetic field
              else

                ! If storing
                if (Atom%Normp(jtran,iz,jdir)%VRAM) then

                  ! Run over the magnetic components
                  do iMu=1,Atom%nMu(jtran)
                    do iMl=1,Atom%nMl(jtran)
                      do iU=1,Atom%nU(jtran)
                        do iL=1,Atom%nL(jtran)

                          ! Easier to write variable
                          d1 = Atom%Normp(jtran,iz,jdir)% &
                                    Norm(iL,iU,iMl,iMu)

                          ! If the norm is not zero
                          if (d1.gt.TINYN) then

                            ! Check close to 1
                            if (d1.lt.BADNORM.or. &
                                d1.gt.2d0-BADNORM) then
                              outofbound(jtran) = outofbound(jtran)+1
                              if (obadnorm) &
                                call writebadbound(folder, &
                                                   Atom%Element, &
                                                   iz,jdir, &
                                                   .True.,.True., &
                                                   jtran, &
                                                   iMu,iMl,iU,iL,d1)
                            end if

                            Atom%Normp(jtran,iz,jdir)% &
                                 prof(iL,iU,iMl,iMu)%cp = &
                              dcmplx(dble(Atom%Normp(jtran,iz,jdir)% &
                                         prof(iL,iU,iMl,iMu)%cp)/d1, &
                                    dimag(Atom%Normp(jtran,iz,jdir)% &
                                              prof(iL,iU,iMl,iMu)%cp))
                          end if

                        end do ! iL
                      end do ! iU
                    end do ! iMl
                  end do ! iMu

                  ! And deallocate the norm because it is not needed
                  deallocate(Atom%Normp(jtran,iz,jdir)%Norm)

                ! Not storing
                else

                  ! Run over the magnetic components
                  do iMu=1,Atom%nMu(jtran)
                    do iMl=1,Atom%nMl(jtran)
                      do iU=1,Atom%nU(jtran)
                        do iL=1,Atom%nL(jtran)

                          ! Easier to write variable
                          d1 = Atom%Normp(jtran,iz,jdir)% &
                                    Norm(iL,iU,iMl,iMu)

                          ! If the norm is not zero
                          if (d1.gt.TINYN) then

                            ! Check close to 1
                            if (d1.lt.BADNORM.or. &
                                d1.gt.2d0-BADNORM) then
                              outofbound(jtran) = outofbound(jtran)+1
                              if (obadnorm) &
                                call writebadbound(folder, &
                                                   Atom%Element,iz, &
                                                   jdir,.True., &
                                                  .True.,jtran, &
                                                  iMu,iMl,iU,iL,d1)
                            end if

                            Atom%Normp(jtran,iz,jdir)% &
                                Norm(iL,iU,iMl,iMu) = 1d0/d1

                          end if

                        end do ! iL
                      end do ! iU
                    end do ! iMl
                  end do ! iMu

                end if ! Storing

              end if

            end do ! transitions
          end do ! directions
        end do ! heights

        ! Check bad limits
        if (maxval(outofbound).gt.0) then

          ! Check if must speak
          if ((MPID%mpi.and.pid.eq.1).or.(.not.MPID%MPI)) then

            ! For each transition
            do jtran=1,Atom%ntran

                ! Check line is affected
                if (outofbound(jtran).gt.0) then
                  write(umsg,'(A,i4,4A,i6,A)') &
                    ' - Warning: transition ',jtran,' in ', &
                    Atom%Element,' has bad normalization', &
                    ' for the chosen width in ',outofbound(jtran), &
                    ' heights, directions, and components.'
                  call verbose
                end if

            end do ! Transition

          end if ! Must output
        end if ! Any transition had bad normalization

        ! Deallocate
        deallocate(outofbound)

      end if ! Master and MPI

      ! Leave the buffer free
      deallocate(buff1)


      !
      ! Fill file of Voigt profiles
      !
      if (vpfil) then

        ! Only slaves or Master if no MPI
        if (pid.ne.0.or..not.MPID%mpi) then

          ! Allocate buffer
          istep = nfreq*2
          allocate(buff1(istep))

          ! Compute f0size and f1size

          ! For each transition
          do jtran=1,Atom%ntran

            ! If this process have frequencies for this line
            if (Atom%fflag(jtran)%absent) cycle

            Atom%f0size(jtran) = 16d0*(Atom%if0(jtran) - &
                                       Atom%rif0(jtran))
            Atom%f1size(jtran) = 16d0*(Atom%rif1(jtran) - &
                                       Atom%if1(jtran))
          end do ! Transitions

          ! For each direction
          do jdir=1,nrdir

            ! Recover the indexes
            ith = ithv(jdir)
            iph = iphv(jdir)

            ! If emergent
            if (LOS) then

              ct = Geom%L_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = cos(Geom%L_phi(iph))
              sc = sin(Geom%L_phi(iph))

            else

              ct = Geom%V_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = Geom%v_mux(iph)
              sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

            end if

            ! Calculate Doppler shift factor

            vfac = 1d0


            ! For each height
            do iz=1,nz

              ! Thermal part of the Doppler width
              DwT = Atom%cDopp*sqrt(Atmo%T(iz))

              ! Check if relevant
              if (Bstrength(iz).lt.TINYB) then
                ! Skip is not going to be used
                if (jdir.gt.nodir) cycle
              end if

              if (dyn) &
                vfac = 1d0 - atmo%vx(iz)*st*cc - &
                             atmo%vy(iz)*st*sc - &
                             atmo%vz(iz)*ct

              ! For each transition
              do jtran=1,Atom%ntran

                ! If this process have frequencies for this line
                if (Atom%fflag(jtran)%absent) cycle

                ! Output Doppler width
                Dw = Atom%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

                ! Find the term indexes for this transition
                itermf = Atom%fst(jtran)%iterml
                itermu = Atom%fst(jtran)%itermu

                ! Get contributions to damping parameter
                au = Atom%damp(itermu,iz)/Dw
                af = Atom%damp(itermf,iz)/Dw
                auf = Atom%ldamp(jtran,iz)/Dw

                ! Get atomic quantities
                S = Atom%Sval(itermu)

                rLu = Atom%rLval(itermu)
                rJumax = rLu+S
                nMu = nint(2d0*rJumax+1d0)

                rLf = Atom%rLval(itermf)
                rJfmax = rLf+S
                nMf = nint(2d0*rJfmax+1d0)

                ! Get indexes
                if0 = Atom%if0(jtran)
                if1 = Atom%if1(jtran)

                ! Get real indexes
                rif0 = Atom%rif0(jtran)
                rif1 = Atom%rif1(jtran)
                nfreqt = rif1 - rif0 + 1

                ! Prepare jump variable
                loffset = dble(Atom%hvifil) + &
                          Atom%dsize(jdir) + &
                          Atom%zsize(iz)

                ! No magnetic field
                if (Bstrength(iz).lt.TINYB) then

                  ! Prepare jump variable
                  loffset = loffset + Atom%tsize(jtran)

                  ! Initialize component variable
                  icom = 0

                  ! sum over Ju
                  do iU=1,Atom%nJ(itermu)

                    eu = Atom%FSfreq(iU,itermu)/Dw

                    ! Get Ju
                    rJu = Atom%rJval(iU,itermu)

                    ! sum over Jl
                    do mF=1,Atom%nJ(itermf)

                      ! Get Jl
                      rJf = Atom%rJval(mF,itermf)

                      el = Atom%FSfreq(mF,itermf)/Dw

                      ! 6-j
                      f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

                      if (abs(f62).lt.TINYJS) cycle

                      ! Advance component index
                      icom = icom + 1

                      !
                      ! Calculate profile
                      !

                      ! Common quantities
                      Dfreqw = eu - el
                      vfacw = vfac/Dw

                      ! Buffer index
                      i1 = 0

                      ! Inverse normalization
                      d1 = Atom%Normp(jtran,iz,jdir)% &
                                Norm(mF,iU,1,1)

                      ! For each frequency
                      do ifreq=if0,if1

                        ! Compute profile
                        call voigt(Dfreqw - Frec%omega(ifreq)*vfacw, &
                                   au+af+auf,prof)

                        ! Store into buffer
                        buff1(i1+1) = dble(prof)*d1
                        buff1(i1+2) = dimag(prof)

                        ! Advance buffer index
                        i1 = i1 + 2

                      end do ! frequencies

                      ! Open file
                      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                         trim(Atom%vfile), &
                                         MPI_MODE_WRONLY, &
                                         MPI_INFO_NULL, &
                                         funit, ierr)
                      if (ierr.ne.0) goto 1100

                      ! Jump
                      loffsetin = loffset + Atom%f0size(jtran) + &
                                  dble(icom - 1)*16d0*nfreqt
                      do while(loffsetin.gt.offlimit)
                        offset = int(offlimit)
                        call MPI_FILE_SEEK(funit, offset, &
                                           MPI_SEEK_CUR, ierr)
                        loffsetin = loffsetin - offlimit
                      end do
                      offset = int(loffsetin)
                      call MPI_FILE_SEEK(funit, offset, &
                                         MPI_SEEK_CUR, ierr)

                      ! Write
                      call MPI_FILE_WRITE(funit, &
                              buff1(1), i1, &
                              MPI_DOUBLE_PRECISION, &
                              MPI_STATUS_IGNORE, ierr)
                      if (ierr.ne.0) goto 1100

                      ! Close file
                      call MPI_FILE_CLOSE(funit, ierr)

                    end do ! Jf
                  end do ! Ju

                ! Yes magnetic field
                else

                  ! Prepare jump variable
                  loffset = loffset + Atom%tBsize(jtran)

                  ! Initialize component variable
                  icom = 0

                  ! sum over Mu
                  do iMu=1,nMu

                    rMu = -rJumax + dble(iMu-1)

                    ! sum over mu_u
                    do iU=1,Atom%nblk(iMu,itermu)

                      eu = Atom%eval(iU,iMu,itermu,iz)/Dw

                      ! sum over Ml
                      do iMf=1,nMf

                        rMf = -rJfmax + dble(iMf-1)

                        if (nint(abs(rMu-rMf)).gt.1) cycle

                        ! sum over mu_l
                        do mF=1,Atom%nblk(iMf,itermf)

                          el = Atom%eval(mF,iMf,itermf,iz)/Dw

                          ! Advance component index
                          icom = icom + 1

                          !
                          ! Calculate profile
                          !

                          ! Common quantities
                          Dfreqw = eu - el + Atom%Dfreq(jtran)/Dw
                          vfacw = vfac/Dw

                          ! Buffer index
                          i1 = 0

                          ! Inverse normalization
                          d1 = Atom%Normp(jtran,iz,jdir)% &
                                    Norm(mF,iU,iMf,iMu)

                          ! For each frequency
                          do ifreq=if0,if1

                            ! Calculate profile
                            call voigt(Dfreqw - &
                                       Frec%omega(ifreq)*vfacw, &
                                       au+af+auf,prof)

                            ! Store into buffer
                            buff1(i1+1) = dble(prof)*d1
                            buff1(i1+2) = dimag(prof)

                            ! Advance buffer index
                            i1 = i1 + 2

                          end do ! frequencies

                          ! Open file
                          call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                             trim(Atom%vfile), &
                                             MPI_MODE_WRONLY, &
                                             MPI_INFO_NULL, &
                                             funit, ierr)
                          if (ierr.ne.0) goto 1100

                          ! Jump
                          loffsetin = loffset + Atom%f0size(jtran) + &
                                      dble(icom - 1)*16d0*nfreqt

                          do while(loffsetin.gt.offlimit)
                            offset = int(offlimit)
                            call MPI_FILE_SEEK(funit, offset, &
                                               MPI_SEEK_CUR, ierr)
                            loffsetin = loffsetin - offlimit
                          end do
                          offset = int(loffsetin)
                          call MPI_FILE_SEEK(funit, offset, &
                                             MPI_SEEK_CUR, ierr)

                          ! Write
                          call MPI_FILE_WRITE(funit, &
                                  buff1(1), i1, &
                                  MPI_DOUBLE_PRECISION, &
                                  MPI_STATUS_IGNORE, ierr)
                          if (ierr.ne.0) goto 1100

                          ! Close file
                          call MPI_FILE_CLOSE(funit, ierr)

                        end do ! iL
                      end do ! Ml
                    end do ! iU
                  end do ! Mu

                end if ! Magnetic field presence

                ! Deallocate
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)

              end do ! Transitions
            end do ! Directions
          end do ! Heights

          ! Deallocate
          deallocate(Atom%Normp)
          nullify(Atom%Normp)

          ! Allocate a token
          allocate(Atom%Normp(1,1,1))

        end if ! Slave or single CPU
      end if ! Writing Voigt profiles

      ! Control
      call control

      return

1000  umsg = 'Error opening '//trim(Atom%vfile)//' file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing '//trim(Atom%vfile)//' file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine normalize

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes normalization factors for all the Voigt profiles\n
      !!      line(LTElines_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bstrength(dfloat(:)): Magnetic field strength\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!     Frec(Frequency_class): Structure with frequency data
      !!            njdir(integer): Number of directions\n
      !!          ithv(integer(:)): Indexing of polar directions\n
      !!          ithv(integer(:)): Indexing of azimuth directions\n
      !!                            signs\n
      !!               lp(logical): If doing formal solution in this
      !!                            run\n
      !!            ofram(logical): Indicates if out of RAM\n
      !!              LOS(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine getprof_LTE(line,Atmo,Bstrength,Geom,MPID,Frec, &
                             njdir,ithv,iphv,lp,ofram,LOS)

      ! I/O

      type(LTEline_class), intent(inout):: line
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(inout):: Geom
      type(Frequency_class), intent(in):: Frec
      type(MPI_class):: MPID
      logical, intent(in):: lp, LOS
      logical, intent(out):: ofram
      integer, intent(in):: njdir
      integer, dimension(:), intent(in):: ithv,iphv
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: field

      integer:: jdir,iz,jj,nodir,nrdir,znjdir
      integer:: ith,iph,ifreq,if0,if1,iMu,iMf

      double precision:: d1,rMu,rMf,dnubw,auf,Dfreqw
      double precision:: DwT,Dw,vfac,vfacw,ct,st,cc,sc

      complex(kind=8):: prof


      ! Routine name
      urou = 'getprof_LTE'

      ! Check if there is magnetic field
      field = .False.

      ! Initialize
      ofram = .False.

      ! Check magnetic field strength
      do iz=1,nz
        if (Bstrength(iz).ge.TINYB) then
          field = .True.
          exit
        end if
      end do

      ! If LOS and not dynamic, if already iterated, maybe no need to
      ! repeat the normalization
      if (LOS.and..not.dyn.and.lp) then

        ! If no magnetic, copy the file if using it and get out
        if (.not.field) then

          ! Everyone control and return
          call control
          return

        end if ! Non-magnetic
      end if ! LOS, not dynamic, previously normalized

      ! Get real size of direction dimension
      if (dyn) then
        nodir = njdir
        nrdir = njdir
      else
        nodir = 1
        if (field) then
          nrdir = njdir
        else
          nrdir = 1
        end if
      end if


      !
      ! Allocate vector for norm (part I)
      !

      ! Check prof is allocated
      if (associated(line%prof)) then
        do jdir=1,size(line%prof,2)
          do iz=lbound(line%prof,1),ubound(line%prof,1)
            if (allocated(line%prof(iz,jdir)%p)) &
              deallocate(line%prof(iz,jdir)%p)
            if (allocated(line%prof(iz,jdir)%comp)) &
              deallocate(line%prof(iz,jdir)%comp)
          end do
        end do
        deallocate(line%prof)
        nullify(line%prof)
      end if

      ! Structure with the norm for each component
      allocate(line%prof(line%Rz0:Rz1,nrdir))


      !
      ! Allocate the norm array
      !

      ! For each height
      do iz=line%Rz0,Rz1

        ! No magnetic field
        if (Bstrength(iz).lt.TINYB) then

          ! For each direction
          do jdir=1,nodir

            ! Skip unless slave
            if (MPID%mpi.and.pid.eq.0) cycle

            ! Allocate profile itself if storing and it is present
            if (VPRAM.and..not.line%absent) then

              ! Prediction
              d1 = 16d-6*dble(line%if1 - line%if0 + 1)

              ! If no more space
              if (floor(MPID%RAM+d1).gt.RLIM) then

                ! No stored
                line%prof(iz,jdir)%VRAM = .False.
                ofram = .True.

              ! If there is space
              else

                ! Storing
                line%prof(iz,jdir)%VRAM = .True.

                ! Allocate
                allocate(line%prof(iz,jdir)%p(line%if0:line%if1))

                ! Update RAM
                MPID%RAM = MPID%RAM + d1
                MPID%VRAM = MPID%VRAM + d1

              end if ! Space to store

            ! Not storing Voigt
            else

              ! No stored
              line%prof(iz,jdir)%VRAM = .False.

            end if ! Storing

          end do ! directions

        ! Yes magnetic field
        else

          ! For each direction
          do jdir=1,nrdir

            ! Skip if not slave
            if (MPID%mpi.and.pid.eq.0) cycle

            ! Get indexes
            if0 = line%if0
            if1 = line%if1

            !
            ! Count components

            ! Initialize
            jj = 0

            ! sum over Mu
            do iMu=1,line%nMu

              rMu = -line%Ju + dble(iMu-1)

              ! sum over Ml
              do iMf=1,line%nMl

                rMf = -line%Jl + dble(iMf-1)

                if (nint(abs(rMu-rMf)).gt.1) cycle

                jj = jj + 1

              end do ! Ml
            end do ! Mu

            ! Allocate profile itself if storing and it is present
            if (VPRAM.and..not.line%absent) then

              ! Prediction
              d1 = 16d-6*dble(line%if1 - line%if0 + 1)

              ! If no more space
              if (floor(MPID%RAM+d1).gt.RLIM) then

                ! No stored
                line%prof(iz,jdir)%VRAM = .False.
                ofram = .True.

              ! If there is space
              else

                ! Storing
                line%prof(iz,jdir)%VRAM = .True.

                ! Allocate
                allocate(line%prof(iz,jdir)%comp(line%nMl,line%nMu))

                ! sum over Mu
                do iMu=1,line%nMu

                  rMu = -line%Ju + dble(iMu-1)

                  ! sum over Ml
                  do iMf=1,line%nMl

                    rMf = -line%Jl + dble(iMf-1)

                    if (nint(abs(rMu-rMf)).gt.1) cycle

                    ! Allocate
                    allocate(line%prof(iz,jdir)% &
                                  comp(iMf,iMu)%cp(if0:if1))

                  end do ! Ml
                end do ! Mu

                ! Update RAM
                MPID%RAM = MPID%RAM + d1
                MPID%VRAM = MPID%VRAM + d1

              end if ! Space to store

            ! Not storing Voigt
            else

              ! No stored
              line%prof(iz,jdir)%VRAM = .False.

            end if ! Storing

          end do ! directions

        end if ! No magnetic field

      end do ! heights


      !
      ! SLAVE OR SINGLE PROCESSOR
      !
      if (.not.MPID%mpi.or.pid.gt.0) then

        !
        ! Calculate normalization
        !

        ! For each height
        do iz=line%Rz0,Rz1

          ! Thermal part of the Doppler width
          DwT = sqrt(line%cDopp*line%cDopp*Atmo%T(iz) + &
                     Atmo%vmi(iz)*Atmo%vmi(iz))

          ! Select number of directions
          if (Bstrength(iz).ge.TINYB) then
            znjdir = nrdir
          else
            znjdir = nodir
          end if

          ! For each direction
          do jdir=1,znjdir

            ! If not storing, why bother
            if (.not.line%prof(iz,jdir)%VRAM) cycle

            ! Recover the indexes
            ith = ithv(jdir)
            iph = iphv(jdir)

            ! If emergent
            if (LOS) then

              ct = Geom%L_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = cos(Geom%L_phi(iph))
              sc = sin(Geom%L_phi(iph))

            else

              ct = Geom%V_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = Geom%v_mux(iph)
              sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

            end if

            ! Calculate Doppler shift factor

            vfac = 1d0

            if (dyn) &
              vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                           atmo%vz(iz)*ct

            ! Output Doppler width
            Dw = line%Dfreq*DwT

            ! Get indexes
            if0 = line%if0
            if1 = line%if1

            ! Common quantities
            Dfreqw = (line%eu - line%el)/Dw
            vfacw = vfac/Dw
            auf = line%damp(iz)/Dw

            !
            ! No magnetic field
            if (Bstrength(iz).lt.TINYB) then

              !
              ! Calculate profile
              !

              ! For each frequency
              do ifreq=if0,if1

                call voigt(Dfreqw - Frec%omega(ifreq)*vfacw, &
                           auf,prof)

                line%prof(iz,jdir)%p(ifreq) = dble(prof)

              end do ! frequencies

            ! Yes magnetic field
            else

              ! sum over Mu
              do iMu=1,line%nMu

                rMu = -line%Ju + dble(iMu-1)

                ! sum over Ml
                do iMf=1,line%nMl

                  rMf = -line%Jl + dble(iMf-1)

                  if (nint(abs(rMu-rMf)).gt.1) cycle

                  dnubw = B2LK*Bstrength(iz)* &
                          (line%gu*rMu - line%gl*rMf)/Dw

                  !
                  ! Calculate profile
                  !

                  ! For each frequency
                  do ifreq=if0,if1

                    call voigt(Dfreqw + dnubw - &
                               Frec%omega(ifreq)*vfacw, &
                               auf,prof)

                    line%prof(iz,jdir)% &
                         comp(iMf,iMu)%cp(ifreq) = prof

                  end do ! frequencies

                end do ! Ml
              end do ! Mu

            end if ! Magnetic field presence

          end do ! output direction
        end do ! height

      end if ! MPI

      ! Control
      call control

      return

      end subroutine getprof_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Validates an existing file with Voigt profiles\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Bstrength(dfloat(:)): Magnetic field strength\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!     Frec(Frequency_class): Structure with frequency data
      !!            njdir(integer): Number of directions\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!            lnorm(logical): Flag for file correctness
      subroutine validatevoigt(Atom,Bstrength,MPID,Frec, &
                               njdir,Flgsg,lnorm)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(MPI_class), intent(inout):: MPID
      type(Frequency_class), intent(in):: Frec
      type(Fctsg_class):: Flgsg
      double precision, dimension(:), intent(in):: Bstrength
      logical, intent(out):: lnorm
      integer, intent(in):: njdir

      ! Local

      logical:: nonvalid, field

      integer:: i,jtran,itermf,itermu
      integer:: idir,jdir,iz,nodir,nrdir
      integer:: ifreq,if0,if1
      integer:: ncom,ncomNB
      integer:: nMu,nMf
      integer:: iMu,iMf
      integer:: iU,mF
      integer:: nL,nU
      integer:: ios,iaux

      double precision:: d1
      double precision:: rLu,rLf,S
      double precision:: rJumax,rJfmax,rJu,rJf
      double precision:: rMu,rMf,f62
      double precision:: d2,d3,d1T,d2T

      ! Buffers, counter, and sizes for MPI
      integer, dimension(0:nproc-1):: nbf1


      ! Routine name
      urou = 'validatevoigt'

      ! Check if there is magnetic field
      field = .False.

      ! Check magnetic field strength
      do iz=1,nz
        if (Bstrength(iz).ge.TINYB) then
          field = .True.
          exit
        end if
      end do

      ! Get real size of direction dimension
      if (dyn) then
        nodir = njdir
        nrdir = njdir
      else
        nodir = 1
        if (field) then
          nrdir = njdir
        else
          nrdir = 1
        end if
      end if

      ! Check Normp is not allocated
      if (associated(Atom%Normp)) then
        do jdir=1,size(Atom%Normp,3)
          do iz=lbound(Atom%Normp,2),ubound(Atom%Normp,2)
            do jtran=lbound(Atom%Normp,1),ubound(Atom%Normp,1)
              if (allocated(Atom%Normp(jtran,iz,jdir)%prof)) &
                deallocate(Atom%Normp(jtran,iz,jdir)%prof)
              if (allocated(Atom%Normp(jtran,iz,jdir)%Norm)) &
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)
            end do
          end do
        end do
        deallocate(Atom%Normp)
        nullify(Atom%Normp)
      end if

      ! Allocate a token
      allocate(Atom%Normp(1,1,1))

      ! Allocate sizes for subcomponents for each transition
      if (.not.allocated(Atom%nL)) allocate(Atom%nL(Atom%ntran))
      if (.not.allocated(Atom%nU)) allocate(Atom%nU(Atom%ntran))
      if (.not.allocated(Atom%nMl)) allocate(Atom%nMl(Atom%ntran))
      if (.not.allocated(Atom%nMu)) allocate(Atom%nMu(Atom%ntran))


      !
      ! Initialize output
      !
      lnorm = .True.
      nonvalid = .False.


      !
      ! Count indexes
      !

      ! Reset index
      nbf1 = 0

      ! For each transition
      do jtran=1,Atom%ntran

        ! Identify the terms
        itermf = Atom%fst(jtran)%iterml
        itermu = Atom%fst(jtran)%itermu

        ! Get atomic quantities
        S = Atom%Sval(itermu)

        rLu = Atom%rLval(itermu)
        rJumax = rLu+S
        nMu = nint(2d0*rJumax+1d0)

        rLf = Atom%rLval(itermf)
        rJfmax = rLf + S
        nMf = nint(2d0*rJfmax+1d0)

        ! Fill two of the sizes arrays
        Atom%nMu(jtran) = nMu
        Atom%nMl(jtran) = nMf

        ! Initialize counters
        nL = 0
        nU = 0

        ! For each Mu
        do iMu=1,nMu

          rMu = -rJumax + dble(iMu-1)

          ! For each mu_u
          do iU=1,Atom%nblk(iMu,itermu)

            ! Update counter
            if (iU.gt.nU) nu = iMu

            ! For each Mf
            do iMf=1,nMf

              rMf = -rJfmax + dble(iMf-1)

              ! If not allowed, skip
              if (nint(abs(rMu-rMf)).gt.1) cycle

              ! Sum over mu_f
              do mF=1,Atom%nblk(iMf,itermf)

                ! Update counter
                if (mF.gt.nL) nL = mf

              end do ! iL
            end do ! Ml
          end do ! iU
        end do ! Mu

        ! Update variables
        Atom%nU(jtran) = nU
        Atom%nL(jtran) = nL

      end do ! Transitions


      ! Deallocate sizes
      if (allocated(Atom%zsize)) deallocate(Atom%zsize)
      if (allocated(Atom%dsize)) deallocate(Atom%dsize)
      if (allocated(Atom%tsize)) deallocate(Atom%tsize)
      if (allocated(Atom%tBsize)) deallocate(Atom%tBsize)
      if (allocated(Atom%f0size)) deallocate(Atom%f0size)
      if (allocated(Atom%f1size)) deallocate(Atom%f1size)
      if (allocated(Atom%i_Vind)) deallocate(Atom%i_Vind)

      ! Allocate sizes
      allocate(Atom%zsize(nz))
      allocate(Atom%dsize(nrdir))
      allocate(Atom%tsize(Atom%ntran))
      allocate(Atom%tBsize(Atom%ntran))
      allocate(Atom%f0size(Atom%ntran))
      allocate(Atom%f1size(Atom%ntran))

      ! Initialize total blocks of sizes and transition index
      d1T = 0d0
      d2T = 0d0
      Atom%zsize(1) = 0d0
      Atom%dsize(1) = 0d0
      Atom%tsize(1) = 0d0
      Atom%tBsize(1) = 0d0
      Atom%Mncom = 0
      ncom = 0
      ncomNB = 0


      ! Allocate indexing of the order
      allocate(Atom%i_Vind(Atom%ntran))

      ! For each transition
      do jtran=1,Atom%ntran

        ! Get real limits
        if0 = Atom%rif0(jtran)
        if1 = Atom%rif1(jtran)

        ! Get frequency size
        if (Atom%fflag(jtran)%absent) then
          Atom%f0size(jtran) = 0d0
          Atom%f1size(jtran) = 0d0
        else
          Atom%f0size(jtran) = 16d0*(Atom%if0(jtran) - if0)
          Atom%f1size(jtran) = 16d0*(if1 - Atom%if1(jtran))
        end if

        ! Get terms
        itermf = Atom%fst(jtran)%iterml
        itermu = Atom%fst(jtran)%itermu

        ! Get atomic quantities
        S = Atom%Sval(itermu)

        rLu = Atom%rLval(itermu)
        rJumax = rLu+S
        nMu = nint(2d0*rJumax+1d0)

        rLf = Atom%rLval(itermf)
        rJfmax = rLf + S
        nMf = nint(2d0*rJfmax+1d0)

        ! Allocate magnetic indexing
        allocate(Atom%i_Vind(jtran)% &
                      ind(maxval(Atom%nblk(1:nMf,itermf)),nMf, &
                          maxval(Atom%nblk(1:nMu,itermu)),nMu))

        ! Allocate non-magnetic indexing
        allocate(Atom%i_Vind(jtran)% &
                      indNB(Atom%nJ(itermf),Atom%nJ(itermu)))

        ! No magnetic field
        d1 = 0d0

        ! Reset index
        i = 0

        ! sum over Ju
        do iU=1,Atom%nJ(itermu)

          ! Get Ju
          rJu = Atom%rJval(iU,itermu)

          ! sum over Jl
          do mF=1,Atom%nJ(itermf)

            ! Get Jl
            rJf = Atom%rJval(mF,itermf)

            ! 6-j
            f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

            if (abs(f62).lt.TINYJS) cycle

            d1 = d1 + dble(if1 - if0 + 1)

            ! Store and advance index
            i = i + 1
            Atom%i_Vind(jtran)%indNB(mF,iU) = i

          end do ! Final levels
        end do ! Upper levels

        ! Store number of components
        Atom%i_Vind(jtran)%ncomNB = i
        if (i.gt.ncomNB) ncomNB = i

        ! Yes magnetic field
        d2 = 0d0

        ! Reset index
        i = 0

        ! sum over Mu
        do iMu=1,nMu

          rMu = -rJumax + dble(iMu-1)

          ! sum over mu_u
          do iU=1,Atom%nblk(iMu,itermu)

            ! sum over Ml
            do iMf=1,nMf

              rMf = -rJfmax + dble(iMf-1)

              if (nint(abs(rMu-rMf)).gt.1) cycle

              ! sum over mu_l
              do mF=1,Atom%nblk(iMf,itermf)

                d2 = d2 + dble(if1 - if0 + 1)

                ! Store and advance index
                i = i + 1
                Atom%i_Vind(jtran)%ind(mF,iMf,iU,iMu) = i

              end do ! mu_l
            end do ! Ml
          end do ! mu_u
        end do ! Mu

        ! Store number of components
        Atom%i_Vind(jtran)%ncom = i
        if (i.gt.ncomNB) ncom = i

        ! If not the last, add size to next
        if (jtran.lt.Atom%ntran) then
          Atom%tsize(jtran+1) = Atom%tsize(jtran) + d1
          Atom%tBsize(jtran+1) = Atom%tBsize(jtran) + d2
        end if

        ! Add to the total size
        d1T = d1T + d1
        d2T = d2T + d2

        ! deallocate indexes of magnetic components
        if (maxval(Bstrength).lt.TINYB) &
          deallocate(Atom%i_Vind(jtran)%ind)

        ! deallocate indexes of non magnetic components
        if (minval(Bstrength).ge.TINYB) &
          deallocate(Atom%i_Vind(jtran)%indNB)

      end do ! Transitions

      ! If no magnetic field
      if (maxval(Bstrength).lt.TINYB) ncom = 0
      ! If only magnetic field
      if (minval(Bstrength).ge.TINYB) ncomNB = 0
      ! Get maximum
      Atom%Mncom = max(ncom,ncomNB)

      ! Inititlize
      d3 = 0d0

      ! For each height
      do iz=1,nz

        ! If previous height was not magnetic
        if (Bstrength(iz).lt.TINYB) then
          if (iz.lt.nz) Atom%zsize(iz+1) = Atom%zsize(iz) + d1T
          d3 = d3 + d1T
        ! If previous was magnetic
        else
          if (iz.lt.nz) Atom%zsize(iz+1) = Atom%zsize(iz) + d2T
          d3 = d3 + d2T
        end if

      end do ! Heights

      ! For each direction
      do idir=2,nrdir
        Atom%dsize(idir) = Atom%dsize(idir-1) + d3
      end do

      Atom%dsize = Atom%dsize*16d0
      Atom%zsize = Atom%zsize*16d0
      Atom%tsize = Atom%tsize*16d0
      Atom%tBsize = Atom%tBsize*16d0

      ! Size of header
      Atom%hvifil = 4*4 + & ! Dimension integers
                    8*nrdir + & ! Directions sizes
                    8*nz + & ! Height sizes
                    4*2*Atom%ntran + & ! Line limits
                    8*2*Atom%ntran + & ! Line sizes
                    8*nfreq ! Frequencies

      ! If not the master, wait for diagnostic
      if (pid.gt.0) then

        call control
        call MPI_BCAST(lnorm, 1, MPI_LOGICAL, 0, &
                       MPI_COMM_RT, ierr)
        return

      end if


      ! Only the master gets down here
      do while(.True.)

        ! Open files
        open(200, file=trim(Atom%vfile), status='unknown', &
             iostat=ios, err=1000, access='stream', &
             action='read', form='unformatted')

        ! Read dimensions

        ! Angles
        read(200, err=1100) iaux
        if (iaux.ne.nrdir) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Numbers of angles do not coincide:',iaux,nrdir
          call verbose
          exit
        end if

        ! Heights
        read(200, err=1100) iaux
        if (iaux.ne.nz) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Numbers of heights do not coincide:',iaux,nz
          call verbose
          exit
        end if

        ! Frequencies
        read(200, err=1100) iaux
        if (iaux.ne.nfreq) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Numbers of frequencies do not coincide:',iaux,nfreq
          call verbose
          exit
        end if

        ! Transitions
        read(200, err=1100) iaux
        if (iaux.ne.Atom%ntran) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Number of transitions do not coincide:', &
            iaux,Atom%ntran
          call verbose
          exit
        end if

        ! Directional size
        do jdir=1,nrdir
          read(200, err=1100) d1
          if (abs(d1-Atom%dsize(jdir)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
              ' # Directional sizes do not coincide:', &
              d1,Atom%dsize(jdir)
            call verbose
            nonvalid = .True.
            exit
          end if
        end do
        if (nonvalid) exit

        ! Height size
        do iz=1,nz
          read(200, err=1100) d1
          if (abs(d1-Atom%zsize(iz)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
              ' # Height sizes do not coincide:',d1,Atom%zsize(iz)
            nonvalid = .True.
            call verbose
            exit
          end if
        end do
        if (nonvalid) exit

        ! For each transition
        do jtran=1,Atom%ntran

          read(200, err=1100) iaux
          if (iaux.ne.Atom%rif0(jtran)) then
            write(umsg,'(A,1x,i3,1x,i3)') &
              ' # Initial indexes of transition do not coincide:', &
              iaux,Atom%rif0(jtran)
            call verbose
            nonvalid = .True.
            exit
          end if

          read(200, err=1100) iaux
          if (iaux.ne.Atom%rif1(jtran)) then
            write(umsg,'(A,1x,i3,1x,i3)') &
                ' # Final indexes of transition do not coincide:', &
                iaux,Atom%rif1(jtran)
            call verbose
            nonvalid = .True.
            exit
          end if

          read(200, err=1100) d1
          if (abs(d1-Atom%tsize(jtran)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
               ' # Transition sizes without field do not coincide:', &
               d1,Atom%tsize(jtran)
            call verbose
            nonvalid = .True.
            exit
          end if

          read(200, err=1100) d1
          if (abs(d1-Atom%tBsize(jtran)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
                ' # Transition sizes with field do not coincide:', &
                d1,Atom%tBsize(jtran)
            call verbose
            nonvalid = .True.
            exit
          end if

        end do ! Transitions
        if (nonvalid) exit

        ! Write frequency
        do ifreq=1,nfreq
          read(200, err=1100) d1
          if (abs(d1-Frec%omega(ifreq)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
                ' # Frequencies do not coincide:',d1,Frec%omega(ifreq)
            call verbose
            nonvalid = .True.
            exit
          end if
        end do
        if (nonvalid) exit

        ! The master does not need the indexes
        if (pid.eq.0.and.MPID%mpi) deallocate(Atom%i_Vind)

        ! If we are here, file is fine
        lnorm = .False.
        exit

      end do ! Infinite loop

      ! Close file
      close(200, err=1100)

      ! Control
      call control

      ! Share diagnostic
      call MPI_BCAST(lnorm, 1, MPI_LOGICAL, 0, &
                     MPI_COMM_RT, ierr)

      return

1000  umsg = 'Error opening '//trim(Atom%vfile)//' file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1100  umsg = 'Error reading '//trim(Atom%vfile)//' file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine validatevoigt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes normalization factors for the multilevel Voigt
      !! profiles\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!    folder(character(500)): Output folder path\n
      !!            njdir(integer): Number of directions\n
      !!          ithv(integer(:)): Indexing of polar directions\n
      !!          ithv(integer(:)): Indexing of azimuth directions\n
      !!              lio(logical): If doing formal solution in this
      !!                            run\n
      !!            ofram(logical): Indicates if out of RAM\n
      !!              LOS(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine normalizeI(Atom,Atmo,Geom,MPID,Frec,folder, &
                            njdir,ithv,iphv,lio,ofram,LOS)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(inout):: Geom
      type(MPI_class), intent(inout):: MPID
      type(Frequency_class), intent(in):: Frec
      character(len=500), intent(in):: folder
      logical, intent(in):: lio, LOS
      logical, intent(out):: ofram
      integer, intent(in):: njdir
      integer, dimension(:), intent(in):: ithv,iphv

      ! Local

      character(len=4):: record

      logical:: extracomm

      integer:: i,i1,jtran,ktran,jftran,itermf,itermu,iJf,iJu
      integer:: idir,jdir,iz,jj,istep,nodir
      integer:: ith,iph,ifreq,if0,if1
      integer:: ios,lcheckram
      integer, dimension(:), allocatable:: outofbound
      integer, dimension(:,:,:), allocatable:: checkram

      double precision:: prof,d1,d2,W0,W1
      double precision:: el,eu
      double precision:: au,af,auf,atuf,Dfreqw
      double precision:: DwT,Dw,vfac,vfacw
      double precision:: ct,st,cc,sc
      double precision:: loffset

      ! Buffers and sizes for MPI
      logical, dimension(:,:,:), allocatable:: tosend
      integer:: finished
      integer, dimension(5):: id, id1_b
      integer, dimension(0:nproc-1):: nbf1
      integer, dimension(:,:,:), allocatable:: size1
      double precision, dimension(:), allocatable:: buff1

      ! MPI offset type
      integer(kind=MPI_OFFSET_KIND):: offset

      ! Routine name
      urou = 'normalizeI'

      ! Initialize
      ofram = .False.


      ! If LOS and not dynamic, if already iterated, no need to
      ! repeat the normalization
      if (LOS.and..not.dyn.and.lio) then

        ! If using files, just copy to new name
        if (vifil.and.pid.eq.0) then

          ! Orginal
          open(200,file='voigt-I-G-'//trim(Atom%file_label), &
               access='stream',status='old',action='read', &
               iostat=ios)

          ! New
          open(300,file=trim(Atom%vfile),access='stream', &
               status='unknown',action='write',iostat=ios)

          ! Copy records until finished
          do while (.True.)

            ! Read characters
            read(200, end=3000) record

            ! Write the same integer
            write(300) record

            cycle
3000        exit

          end do

          close(200)
          close(300)

        end if ! Copy vfile

        ! Everyone control and return
        call control
        return

      end if ! LOS, not dynamic, previously normalized


      !
      ! Initialize comm flag
      !
      if (nproc.gt.cpulimit) then
        extracomm = .True.
      else
        extracomm = .False.
      end if

      ! Get real size of the direction dimension
      if (dyn) then
        nodir = njdir
      else
        nodir = 1
      end if


      !
      ! If using a file, prepare for it
      !
      if (vifil) then

        ! Deallocate sizes
        if (allocated(Atom%zsize)) deallocate(Atom%zsize)
        if (allocated(Atom%dsize)) deallocate(Atom%dsize)
        if (allocated(Atom%tsize)) deallocate(Atom%tsize)
        if (allocated(Atom%f0size)) deallocate(Atom%f0size)

        ! Allocate sizes
        allocate(Atom%zsize(nz))
        allocate(Atom%dsize(nodir))
        allocate(Atom%tsize(Atom%nftran))
        allocate(Atom%f0size(Atom%nftran))

        ! Initialize total blocks of sizes and transition index
        d1 = 0d0
        ktran = 0
        Atom%zsize(1) = 0d0
        Atom%dsize(1) = 0d0
        Atom%tsize(1) = 0d0

        ! For each transition
        do jtran=1,Atom%ntran

          ! For each FS transition
          do jftran=1,Atom%fst(jtran)%nt

            ! Advance index
            ktran = ktran + 1

            ! Skip first
            if (ktran.eq.1) then
              ! Add last
              d2 = dble(Atom%rif1(Atom%ntran) - &
                        Atom%rif0(Atom%ntran) + 1)
              d1 = d1 + d2
              cycle
            end if

            ! If first of the FS transitions
            if (jftran.eq.1) then

              ! Get indexes
              if0 = Atom%rif0(jtran-1)
              if1 = Atom%rif1(jtran-1)

            ! Not the first FS
            else

              ! Get indexes
              if0 = Atom%rif0(jtran)
              if1 = Atom%rif1(jtran)

            end if ! First FS transition

            ! Size of previous
            d2 = dble(if1 - if0 + 1)
            Atom%tsize(ktran) = Atom%tsize(ktran-1) + d2
            d1 = d1 + d2

          end do ! FS transitions
        end do ! Transitions

        ! Add first height to total
        d2 = d1

        ! For each height
        do iz=2,nz

          Atom%zsize(iz) = Atom%zsize(iz-1) + d1
          d2 = d2 + d1

        end do ! Heights

        ! For each direction
        do idir=2,nodir
          Atom%dsize(idir) = Atom%dsize(idir-1) + d2
        end do

        Atom%dsize = Atom%dsize*8d0
        Atom%zsize = Atom%zsize*8d0
        Atom%tsize = Atom%tsize*8d0

        ! Size of header
        Atom%hvifil = 4*4 + & ! Dimension integers
                      8*nodir + & ! Directions sizes
                      8*nz + & ! Height sizes
                      4*2*Atom%nftran + & ! Line limits
                      8*Atom%nftran + & ! Line sizes
                      8*nfreq ! Frequencies

        !
        ! Only Master writes
        !
        if (pid.eq.0) then

          ! Open files
          open(200, file=trim(Atom%vfile), status='unknown', &
               iostat=ios, err=1000, access='stream', &
               action='write', form='unformatted')

          ! Write dimensions
          write(200, err=1100) nodir
          write(200, err=1100) nz
          write(200, err=1100) nfreq
          write(200, err=1100) Atom%nftran
          write(200, err=1100) Atom%dsize
          write(200, err=1100) Atom%zsize

          ! Initialize index
          ktran = 0

          ! For each transition
          do jtran=1,Atom%ntran

            ! For each FS transition
            do jftran=1,Atom%fst(jtran)%nt

              ! Advance index
              ktran = ktran + 1

              write(200, err=1100) Atom%rif0(jtran)
              write(200, err=1100) Atom%rif1(jtran)
              write(200, err=1100) Atom%tsize(ktran)

            end do ! FS transitions
          end do ! Transitions

          ! Write frequency
          write(200, err=1100) Frec%omega

          ! And close this
          close(200, err=1100)

        end if ! Master

        ! Control
        call control

      end if ! Voigt file


      !
      ! Allocate vector for norm (part I)
      !

      ! Check Normp is not allocated
      if (associated(Atom%Normp)) then
        do jdir=1,size(Atom%Normp,3)
          do iz=lbound(Atom%Normp,2),ubound(Atom%Normp,2)
            do jtran=lbound(Atom%Normp,1),ubound(Atom%Normp,1)
              if (allocated(Atom%Normp(jtran,iz,jdir)%prof)) then
                deallocate(Atom%Normp(jtran,iz,jdir)%prof)
              else if (allocated(Atom%Normp(jtran,iz,jdir)%Norm)) then
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)
              end if
            end do
          end do
        end do
        deallocate(Atom%Normp)
        nullify(Atom%Normp)
      end if

      ! Structure with the norm for each component
      allocate(Atom%Normp(Atom%ntran,Rz0:Rz1,nodir))
      ! Allocate size for MPI
      allocate(size1(Atom%ntran,Rz0:Rz1,nodir))
      ! Allocate bool for MPI
      allocate(tosend(Atom%ntran,Rz0:Rz1,nodir))
      ! Initialize tosend
      tosend = .False.

      ! Allocate the check for the master
      if (pid.eq.0.and.MPID%mpi) then
        allocate(checkram(Atom%ntran,Rz0:Rz1,nodir))
        if (VIRAM) checkram = 1
      end if



      !
      ! Allocate the norm array
      !

      ! For each direction
      do jdir=1,nodir

        ! For each height
        do iz=Rz0,Rz1

          ! For each transition
          do jtran=1,Atom%ntran

            ! Allocate
            allocate(Atom%Normp(jtran,iz,jdir)%Norm( &
                       Atom%fst(jtran)%nt,1,1,1))
            Atom%Normp(jtran,iz,jdir)%Norm = 0d0

            ! Determine size of this array
            size1(jtran,iz,jdir) = Atom%fst(jtran)%nt

            ! Skip unless worker
            if (MPID%mpi.and.pid.eq.0) cycle

            ! Allocate profile itself if storing and it is present
            if (VIRAM.and..not.Atom%fflag(jtran)%absent) then

              ! Predict
              d1 = 8d-6*dble((Atom%if1(jtran) - &
                              Atom%if0(jtran) + 1)* &
                              Atom%fst(jtran)%nt)

              ! If no more space
              if (floor(MPID%RAM+d1).gt.RLIM) then

                ! No stored
                Atom%Normp(jtran,iz,jdir)%VRAM = .False.
                ofram = .True.

              ! If there is space
              else

                ! Storing
                Atom%Normp(jtran,iz,jdir)%VRAM = .True.

                ! Allocate
                allocate(Atom%Normp(jtran,iz,jdir)%prof( &
                         Atom%fst(jtran)%nt,1,1,1))

                ! For each fs transition, allocate profile
                do jj=1,Atom%fst(jtran)%nt

                  ! Allocate
                  allocate(Atom%Normp(jtran,iz,jdir)%prof( &
                              jj,1,1,1)%p( &
                              Atom%if0(jtran):Atom%if1(jtran)))

                end do ! fs transitions

                MPID%RAM = MPID%RAM + d1
                MPID%VRAM = MPID%VRAM + d1

              end if ! Space to store

            ! Not storing Voigt
            else

              ! No stored
              Atom%Normp(jtran,iz,jdir)%VRAM = .False.

              ! Add normalization
              d1 = 8d-6*dble(Atom%fst(jtran)%nt)
              MPID%RAM = MPID%RAM + d1
              MPID%VRAM = MPID%VRAM + d1

            end if

          end do ! For each transition
        end do ! For each direction
      end do ! For each height

      ! Check the maximum size to transfer
      ios = maxval(size1)

      ! Gather the maximum size that each processor is holding
      do while (.True.)
        call MPI_ALLGATHER(ios,1,MPI_INTEGER,nbf1(0),1, &
                           MPI_INTEGER,MPI_COMM_RT,ierr)
        if (ierr.eq.0) exit
      end do

      ! Allocate buffers
      allocate(buff1(nbf1(pid)))


      !
      ! MASTER
      !

      if (MPID%mpi.and.pid.eq.0) then

        !
        ! Initialize finished
        !
        finished = 1


        !
        ! Calculate normalization
        !

        do while (finished.lt.nproc)

          ! Receive the informative package with indexes
          if (extracomm) then
            do while (.True.)
              call MPI_recv(id(1),5,MPI_INTEGER, &
                            MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do
          else
            call MPI_recv(id(1),5,MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
          end if

          ! If ending signal
          if (id(1).lt.0) then
            finished = finished + 1
            cycle
          end if

          if (extracomm) then
            do while (.True.)
              call MPI_SEND(id(1),1,MPI_INTEGER,id(1),id(1), &
                            MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do
          end if

          ! Receive the buffer with the integral
          if (extracomm) then
            do while (.True.)
              call MPI_recv(buff1(1), nbf1(id(1)), &
                            MPI_DOUBLE_PRECISION, id(1), &
                            id(1), MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do
          else
            call MPI_recv(buff1(1), nbf1(id(1)), &
                          MPI_DOUBLE_PRECISION, id(1), &
                          id(1), MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
          end if

          ! Store checkram
          if (VIRAM.and.id(5).lt.0) &
            checkram(id(2),id(3),id(4)) = -1

          ! For each FS transition
          do jj=1,Atom%fst(id(2))%nt

            ! Accumulate the sub-integrals
            Atom%Normp(id(2),id(3),id(4))% &
                 Norm(jj,1,1,1) = &
                    Atom%Normp(id(2),id(3),id(4))% &
                       Norm(jj,1,1,1) + buff1(jj)

          end do ! FS transition
        end do ! Communications to do


      !
      ! SLAVE OR SINGLE PROCESSOR
      !

      else

        !
        ! Calculate normalization
        !

        ! For each direction
        do jdir=1,nodir

          ! Recover the indexes
          ith = ithv(jdir)
          iph = iphv(jdir)

          ! If emergent
          if (LOS) then

            ct = Geom%L_mu(ith)
            st = sqrt(1d0 - ct*ct)
            cc = cos(Geom%L_phi(iph))
            sc = sin(Geom%L_phi(iph))

          else

            ct = Geom%V_mu(ith)
            st = sqrt(1d0 - ct*ct)
            cc = Geom%v_mux(iph)
            sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

          end if

          ! For each height
          do iz=Rz0,Rz1

            ! Thermal part of the Doppler width
            DwT = Atom%cDopp*sqrt(Atmo%T(iz))

            ! Calculate Doppler shift factor

            vfac = 1d0

            if (dyn) &
              vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                           atmo%vz(iz)*ct

            ! For each transition
            do jtran=1,Atom%ntran

              ! If the line is in this process, skip
              if (.not.Atom%fflag(jtran)%absent) then

                ! Flag to send
                tosend(jtran,iz,jdir) = .True.

                ! Find the term indexes for this transition
                do i=1,Atom%nMulti-1
                  do i1=i+1,Atom%nMulti
                    if (Atom%irad(i,i1).eq.jtran) then
                      itermf = i
                      itermu = i1
                    end if
                  end do
                end do

                ! Get contributions to damping parameter
                au = Atom%damp(itermu,iz)
                af = Atom%damp(itermf,iz)
                auf = Atom%ldamp(jtran,iz)

                ! Get indexes
                if0 = Atom%if0(jtran)
                if1 = Atom%if1(jtran)
                ! Get Weights
                W0 = Atom%W0(jtran)
                W1 = Atom%W1(jtran)


                !
                ! Proper normalization
                !

                ! For each FS transition
                do jj=1,Atom%fst(jtran)%nt

                  ! Find the level indexes for this transition
                  do i=1,Atom%nJ(itermf)
                    do i1=1,Atom%nJ(itermu)
                      if (Atom%fst(jtran)%irad(i1,i).eq.jj) then
                        iJf = i
                        iJu = i1
                      end if
                    end do
                  end do

                  ! Get the FS energies
                  eu = Atom%FSfreq(iJu,itermu)
                  el = Atom%FSfreq(iJf,itermf)

                  ! Output Doppler width
                  Dw = (eu - el)*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)


                  !
                  ! Calculate profile
                  !

                  ! Common quantities
                  Dfreqw = (eu - el)/Dw
                  vfacw = vfac/Dw
                  d1 = 1d-5/(sqrt(PI)*Dw)
                  atuf = (au+af+auf)/Dw

                  ! Boundaries

                  ! Lower
                  call voigtI(Dfreqw - Frec%omega(if0)*vfacw, &
                             atuf,prof)
                  Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1) = &
                                                    prof*(W0*d1)

                  ! Store
                  if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                    Atom%Normp(jtran,iz,jdir)%prof(jj,1,1,1)% &
                                                    p(if0) = prof

                  ! For each frequency
                  do ifreq=if0+1,if1-1

                    call voigtI(Dfreqw - Frec%omega(ifreq)*vfacw, &
                               atuf,prof)

                    ! Add to the integral
                    Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1) = &
                    Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1) + &
                                      prof*(Frec%W_freq(ifreq)*d1)

                    ! Store
                    if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                      Atom%Normp(jtran,iz,jdir)%prof(jj,1,1,1)% &
                                                    p(ifreq) = prof

                  end do ! Frequencies

                  ! Upper
                  call voigtI(Dfreqw - Frec%omega(if1)*vfacw, &
                             atuf,prof)
                  Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1) = &
                  Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1) + &
                                                    prof*(W1*d1)

                  ! Store
                  if (Atom%Normp(jtran,iz,jdir)%VRAM) &
                    Atom%Normp(jtran,iz,jdir)%prof(jj,1,1,1)% &
                                                    p(if1) = prof

                end do ! Fine transition

              end if ! Presence test

              if (MPID%mpi.and.tosend(jtran,iz,jdir)) then

                !
                ! Share data with master
                !

                ! Check last send was received
                if (.not.extracomm) then
                  call MPI_WAIT(MPID%request1,MPI_STATUS_IGNORE, &
                                ierr)
                  call MPI_WAIT(MPID%request2,MPI_STATUS_IGNORE, &
                                ierr)
                end if

                ! Send the indexing data
                if (Atom%Normp(jtran,iz,jdir)%VRAM) then
                  id1_b = (/ pid, jtran, iz, jdir, 1 /)
                else
                  id1_b = (/ pid, jtran, iz, jdir,-1 /)
                end if
                if (extracomm) then
                  do while (.True.)
                    call MPI_SEND(id1_b(1),5,MPI_INTEGER,0,0, &
                                  MPI_COMM_RT,ierr)
                    if (ierr.eq.0) exit
                  end do
                else
                  call MPI_ISEND(id1_b(1),5,MPI_INTEGER,0,0, &
                                 MPI_COMM_RT,MPID%request1,ierr)
                end if

                ! Reorder the normalization into the buffer
                buff1(1:size1(jtran,iz,jdir)) = &
                             Atom%Normp(jtran,iz,jdir)%Norm(:,1,1,1)

                if (extracomm) then
                  do while (.True.)
                    call MPI_recv(id1_b(1),1,MPI_INTEGER, &
                                  0, pid, MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)
                    if (ierr.eq.0) exit
                  end do
                end if

                ! Send the actual normalization values
                if (extracomm) then
                  do while (.True.)
                    call MPI_SEND(buff1(1),nbf1(pid), &
                                  MPI_DOUBLE_PRECISION, 0, pid, &
                                  MPI_COMM_RT, ierr)
                    if (ierr.eq.0) exit
                  end do
                else
                  call MPI_ISEND(buff1(1),nbf1(pid), &
                                 MPI_DOUBLE_PRECISION, 0, pid, &
                                 MPI_COMM_RT, &
                                 MPID%request2, ierr)
                end if
              end if ! MPI

            end do ! output transition
          end do ! output direction
        end do ! height

        ! If MPI send finished signal
        if (MPID%mpi) then

          ! Check last send was received
          if (.not.extracomm) &
            call MPI_WAIT(MPID%request1,MPI_STATUS_IGNORE,ierr)

          ! Send the indexing data
          id = (/ -1, -1, -1, -1, -1 /)
          if (extracomm) then
            do while (.True.)
              call MPI_SEND(id(1),5,MPI_INTEGER,0,0, &
                             MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do
          else
            call MPI_ISEND(id(1),5,MPI_INTEGER,0,0, &
                            MPI_COMM_RT,MPID%request1,ierr)
          end if ! Extracommunication
        end if ! MPI
      end if ! Master or slave (or not mpi)

      ! If MPI
      if (MPID%mpi) then

        !
        ! Broadcast the results
        !

        ! Slaves allocate checkram
        if (pid.gt.0.and.VIRAM) then
          ! Allocate checkram
          allocate(checkram(Atom%ntran,Rz0:Rz1,nodir))
        end if

        ! Send norm

        ! For each direction
        do jdir=1,nodir

          ! For each height
          do iz=Rz0,Rz1

            ! For each transition
            do jtran=1,Atom%ntran

              ! Alternative bcast
              if (MPID%altbcast) then

                ! If not master, receive first
                if (pid.ne.0) then

                  ! Receive Norm
                  call MPI_RECV(Atom%Normp(jtran,iz,jdir)% &
                                Norm(1,1,1,1), &
                                size1(jtran,iz,jdir), &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                end if ! No master

                ! For each send
                do istep=1,MPID%nsend

                  ! Send Norm
                  call MPI_ISEND(Atom%Normp(jtran,iz,jdir)% &
                                 Norm(1,1,1,1), &
                                 size1(jtran,iz,jdir), &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,1), ierr)

                end do ! Sends

                ! For each slave
                do istep=1,MPID%nsend

                  ! Wait for everyone to receive the radiation data
                  ! continuing
                  call MPI_WAIT(MPID%requestA(istep,1), &
                                MPI_STATUS_IGNORE,ierr)

                end do ! Sends

              ! Normal bcast
              else

                ! Share Norm
                call MPI_BCAST(Atom%Normp(jtran,iz,jdir)% &
                                    Norm(1,1,1,1), &
                               size1(jtran,iz,jdir), &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

              end if ! Type of bcast

            end do ! transitions
          end do ! directions
        end do ! heights

        ! If storing in RAM
        if (VIRAM) then

          ! Send checkram
          lcheckram = Atom%ntran*Rnz*nodir

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive Norm
              call MPI_RECV(checkram(1,Rz0,1),lcheckram, &
                            MPI_INTEGER,MPID%recv, pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No master

            ! For each send
            do istep=1,MPID%nsend

              ! Send Norm
              call MPI_ISEND(checkram(1,Rz0,1),lcheckram, &
                             MPI_INTEGER,MPID%lsend(istep), &
                             MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,1), ierr)

            end do ! Sends

            ! For each slave
            do istep=1,MPID%nsend

              ! Wait for everyone to receive the radiation data
              ! continuing
              call MPI_WAIT(MPID%requestA(istep,1), &
                            MPI_STATUS_IGNORE,ierr)

            end do ! Sends

          ! Normal bcast
          else

            ! Share Norm
            call MPI_BCAST(checkram(1,Rz0,1),lcheckram, &
                           MPI_INTEGER, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast

          ! Slaves deal with it
          if (pid.gt.0) then

            ! For each direction
            do jdir=1,nodir

              ! For each height
              do iz=Rz0,Rz1

                ! For each transition
                do jtran=1,Atom%ntran

                  ! If was saving but cannot
                  if (Atom%Normp(jtran,iz,jdir)%VRAM.and. &
                      checkram(jtran,iz,jdir).lt.0) then

                    ! Not storing now
                    Atom%Normp(jtran,iz,jdir)%VRAM = .False.

                    ! For each fs transition, allocate profile
                    do jj=1,Atom%fst(jtran)%nt

                      ! Deallocate
                      deallocate(Atom%Normp(jtran,iz,jdir)%prof( &
                                            jj,1,1,1)%p)

                    end do ! fs transitions

                    ! And deallocate prof
                    deallocate(Atom%Normp(jtran,iz,jdir)%prof)

                    ! Remove profile
                    d1 = 8d-6*dble((Atom%if1(jtran) - &
                                    Atom%if0(jtran) + 1)* &
                                    Atom%fst(jtran)%nt)
                    MPID%RAM = MPID%RAM - d1
                    MPID%VRAM = MPID%VRAM - d1

                    ! Add norm
                    d1 = 8d-6*dble(Atom%fst(jtran)%nt)
                    MPID%RAM = MPID%RAM + d1
                    MPID%VRAM = MPID%VRAM + d1

                  end if ! Was string and not now

                end do ! transition
              end do ! height
            end do ! direction

          end if ! Slaves
        end if ! If saving in RAM
      end if ! MPI


      !
      ! Calculate multiplicative factor (instead of division factor)
      !

      ! If MPI, master does not need this
      if (MPID%mpi.and.pid.eq.0) then

        do jdir=1,size(Atom%Normp,3)
          do iz=Rz0,Rz1
            do jtran=1,Atom%ntran
              if (allocated(Atom%Normp(jtran,iz,jdir)%Norm)) &
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)
            end do
          end do
        end do
        deallocate(Atom%Normp)
        nullify(Atom%Normp)

      else

        ! Allocate and initialize out of bound counter
        allocate(outofbound(Atom%nftran))
        outofbound = 0

        ! For each direction
        do jdir=1,nodir

          ! For each height
          do iz=Rz0,Rz1

            ! Initialize transition index
            i1 = 0

            ! For each transition
            do jtran=1,Atom%ntran

              ! If storing
              if (Atom%Normp(jtran,iz,jdir)%VRAM) then

                ! For each FS transition
                do jj=1,Atom%fst(jtran)%nt

                  ! Advance index
                  i1 = i1 + 1

                  ! Easier to write variable
                  d1 = Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1)

                  ! If the norm is not zero
                  if (d1.gt.TINYN) then

                    ! Check close to 1
                    if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then
                      outofbound(i1) = outofbound(i1) + 1
                      if (obadnorm) &
                        call writebadbound(folder, &
                                           Atom%Element,iz,jdir, &
                                           .False.,.False.,jtran, &
                                           jj,1,1,1,d1)
                    end if

                    Atom%Normp(jtran,iz,jdir)%prof(jj,1,1,1)%p = &
                         Atom%Normp(jtran,iz,jdir)%prof(jj,1,1,1)%p/d1

                  end if

                end do ! FS transitions

                ! And deallocate the norm because it is not needed
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)

              ! Not storing
              else

                ! For each FS transition
                do jj=1,Atom%fst(jtran)%nt

                  ! Advance index
                  i1 = i1 + 1

                  ! Easier to write variable
                  d1 = Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1)

                  ! If the norm is not zero
                  if (d1.gt.TINYN) then

                    ! Check close to 1
                    if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then
                      outofbound(i1) = outofbound(i1) + 1
                      if (obadnorm) &
                        call writebadbound(folder, &
                                           Atom%Element,iz,jdir, &
                                           .False.,.False.,jtran, &
                                           jj,1,1,1,d1)
                    end if

                    Atom%Normp(jtran,iz,jdir)%Norm(jj,1,1,1) = 1d0/d1

                  end if

                end do ! FS transitions

              end if ! Storing

            end do ! transitions
          end do ! Directions
        end do ! heights

        ! Check bad limits
        if (maxval(outofbound).gt.0) then

          ! Check if must speak
          if ((MPID%mpi.and.pid.eq.1).or.(.not.MPID%mpi)) then

            ! Initialize index
            i1 = 0

            ! For each transition
            do jtran=1,Atom%ntran

              ! For each FS transition
              do jj=1,Atom%fst(jtran)%nt

                ! Advance index
                i1 = i1 + 1

                ! Check line is affected
                if (outofbound(i1).gt.0) then
                  write(umsg,'(A,i4,",",i4,4A,i4,A)') &
                    ' # Warning: transition ',jtran,jj,' in ', &
                    Atom%Element,' has bad normalization', &
                    ' for the chosen width in ',outofbound(i1), &
                    ' heights X directions.'
                  call verbose
                end if

              end do ! FS transition
            end do ! Transition

          end if ! Must output
        end if ! Any transition had bad normalization

        ! Deallocate
        deallocate(outofbound)

      end if ! Master and MPI

      ! Leave the buffer free
      deallocate(buff1)


      !
      ! Fill file of Voigt profiles
      !
      if (vifil) then

        ! Only slaves or Master if no MPI
        if (pid.ne.0.or..not.MPID%mpi) then

          ! Allocate buffer
          allocate(buff1(nfreq))

          !
          ! Compute frequency size

          ! Initialize index
          ktran = 0

          ! For each transition
          do jtran=1,Atom%ntran

            ! Check if present
            if (Atom%fflag(jtran)%absent) then
              ktran = ktran + Atom%fst(jtran)%nt
              cycle
            end if

            ! Get indexes
            if0 = Atom%if0(jtran)
            if1 = Atom%if1(jtran)

            ! For each FS transition
            do jftran=1,Atom%fst(jtran)%nt

              ! Advance index
              ktran = ktran + 1

              ! Define
              Atom%f0size(ktran) = 8d0*dble(if0 - Atom%rif0(jtran))

            end do

          end do

          ! For each direction
          do jdir=1,nodir

            ! Recover the indexes
            ith = ithv(jdir)
            iph = iphv(jdir)

            ! If emergent
            if (LOS) then

              ct = Geom%L_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = cos(Geom%L_phi(iph))
              sc = sin(Geom%L_phi(iph))

            else

              ct = Geom%V_mu(ith)
              st = sqrt(1d0 - ct*ct)
              cc = Geom%v_mux(iph)
              sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

            end if

            ! For each height
            do iz=1,nz

              ! Thermal part of the Doppler width
              DwT = Atom%cDopp*sqrt(Atmo%T(iz))

              ! Calculate Doppler shift factor

              vfac = 1d0

              if (dyn) &
                vfac = 1d0 - atmo%vx(iz)*st*cc - &
                             atmo%vy(iz)*st*sc - &
                             atmo%vz(iz)*ct

              ! Initialize index
              ktran = 0

              ! For each transition
              do jtran=1,Atom%ntran

                ! Check if present
                if (Atom%fflag(jtran)%absent) then
                  ktran = ktran + Atom%fst(jtran)%nt
                  cycle
                end if

                ! Find the term indexes for this transition
                itermf = Atom%fst(jtran)%iterml
                itermu = Atom%fst(jtran)%itermu

                ! Get contributions to damping parameter
                au = Atom%damp(itermu,iz)
                af = Atom%damp(itermf,iz)
                auf = Atom%ldamp(jtran,iz)

                ! Get indexes
                if0 = Atom%if0(jtran)
                if1 = Atom%if1(jtran)

                ! For each FS transition
                do jftran=1,Atom%fst(jtran)%nt

                  ! Advance index
                  ktran = ktran + 1

                  ! Find the level indexes for this transition
                  do i=1,Atom%nJ(itermf)
                    do i1=1,Atom%nJ(itermu)
                      if (Atom%fst(jtran)%irad(i1,i).eq.jftran) then
                        iJf = i
                        iJu = i1
                      end if
                    end do
                  end do

                  ! Get the FS energies
                  eu = Atom%FSfreq(iJu,itermu)
                  el = Atom%FSfreq(iJf,itermf)

                  ! Output Doppler width
                  Dw = (eu - el)*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

                  ! Open file
                  call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                     trim(Atom%vfile), &
                                     MPI_MODE_WRONLY, &
                                     MPI_INFO_NULL, &
                                     funit, ierr)
                  if (ierr.ne.0) goto 1100


                  ! Jump
                  loffset = dble(Atom%hvifil) + &
                            Atom%dsize(jdir) + &
                            Atom%zsize(iz) + &
                            Atom%tsize(ktran) + &
                            Atom%f0size(ktran)
                  do while(loffset.gt.offlimit)
                    offset = int(offlimit)
                    call MPI_FILE_SEEK(funit, offset, &
                                       MPI_SEEK_CUR, ierr)
                    loffset = loffset - offlimit
                  end do
                  offset = int(loffset)
                  call MPI_FILE_SEEK(funit, offset, &
                                     MPI_SEEK_CUR, ierr)


                  !
                  ! Compute profiles and write
                  !

                  ! Common quantities
                  Dfreqw = (eu - el)/Dw
                  vfacw = vfac/Dw
                  atuf = (au+af+auf)/Dw

                  ! For each frequency
                  do ifreq=if0,if1

                    call voigtI(Dfreqw - Frec%omega(ifreq)*vfacw, &
                                atuf,buff1(ifreq))

                  end do ! Frequencies

                  ! Normalize
                  buff1(if0:if1) = buff1(if0:if1)* &
                          Atom%Normp(jtran,iz,jdir)%Norm(jftran,1,1,1)

                  ! Write
                  istep = if1 - if0 + 1
                  call MPI_FILE_WRITE(funit, buff1(if0), &
                                      istep, MPI_DOUBLE_PRECISION, &
                                      MPI_STATUS_IGNORE, ierr)
                  if (ierr.ne.0) goto 1100

                  ! Close file
                  call MPI_FILE_CLOSE(funit, ierr)

                end do ! Fine transition

                ! Deallocate
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)

              end do ! Transitions
            end do ! Heights
          end do ! Directions

          ! Deallocate
          deallocate(Atom%Normp)
          nullify(Atom%Normp)

          ! Allocate a token
          allocate(Atom%Normp(1,1,1))

        end if ! Slave or single CPU
      end if ! Writing Voigt profiles

      ! Control
      call control

      return

1000  umsg = 'Error opening '//trim(Atom%vfile)//' file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing '//trim(Atom%vfile)//' file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine normalizeI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes profiles for intensity LTE lines\n
      !!      line(LTElines_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            njdir(integer): Number of directions\n
      !!          ithv(integer(:)): Indexing of polar directions\n
      !!          ithv(integer(:)): Indexing of azimuth directions\n
      !!              lio(logical): If doing formal solution in this
      !!                            run\n
      !!            ofram(logical): Indicates if out of RAM\n
      !!              LOS(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine getprofI_LTE(line,Atmo,Geom,MPID,Frec,njdir,ithv, &
                              iphv,lio,ofram,LOS)

      ! I/O

      type(LTEline_class), intent(inout):: line
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(inout):: Geom
      type(MPI_class), intent(inout):: MPID
      type(Frequency_class), intent(in):: Frec
      logical, intent(in):: lio, LOS
      logical, intent(out):: ofram
      integer, intent(in):: njdir
      integer, dimension(:), intent(in):: ithv,iphv

      ! Local

      integer:: jdir,iz,nodir,ith,iph,ifreq,if0,if1

      double precision:: prof,d1,el,eu,atuf,Dfreqw
      double precision:: DwT,Dw,vfac,vfacw,ct,st,cc,sc


      ! Routine name
      urou = 'getprofI_LTE'

      ! Initialize
      ofram = .False.


      ! If LOS and not dynamic, if already iterated, no need to
      ! repeat the normalization
      if (LOS.and..not.dyn.and.lio) then

        ! Everyone control and return
        call control
        return

      end if ! LOS, not dynamic, previously normalized

      ! Get real size of the direction dimension
      if (dyn) then
        nodir = njdir
      else
        nodir = 1
      end if

      !
      ! Allocate vector for norm (part I)
      !

      ! Check prof is allocated
      if (associated(line%prof)) then
        do jdir=1,size(line%prof,2)
          do iz=lbound(line%prof,1),ubound(line%prof,1)
            if (allocated(line%prof(iz,jdir)%p)) &
              deallocate(line%prof(iz,jdir)%p)
            if (allocated(line%prof(iz,jdir)%comp)) &
              deallocate(line%prof(iz,jdir)%comp)
          end do
        end do
        deallocate(line%prof)
        nullify(line%prof)
      end if

      ! Structure with the norm for each component
      allocate(line%prof(line%Rz0:Rz1,nodir))

      !
      ! Allocate the norm array
      !

      ! For each direction
      do jdir=1,nodir

        ! For each height
        do iz=line%Rz0,Rz1

          ! Skip unless worker
          if (MPID%mpi.and.pid.eq.0) cycle

          ! Allocate profile itself if storing and it is present
          if (LVIRAM.and..not.line%absent) then

            ! Predict
            d1 = 8d-6*dble(line%if1 - line%if0 + 1)

            ! If no more space
            if (floor(MPID%RAM+d1).gt.RLIM) then

              ! No stored
              line%prof(iz,jdir)%VRAM = .False.
              ofram = .True.

            ! If there is space
            else

              ! Storing
              line%prof(iz,jdir)%VRAM = .True.

              ! Allocate
              allocate(line%prof(iz,jdir)%p(line%if0:line%if1))

              MPID%RAM = MPID%RAM + d1
              MPID%VRAM = MPID%VRAM + d1

            end if ! Space to store

          ! Not storing Voigt
          else

            ! No stored
            line%prof(iz,jdir)%VRAM = .False.

          end if

        end do ! For each direction
      end do ! For each height


      !
      ! SLAVE OR SINGLE PROCESSOR
      !

      if (pid.gt.0.or..not.MPID%mpi) then

        !
        ! Calculate profiles
        !

        ! For each direction
        do jdir=1,nodir

          ! Recover the indexes
          ith = ithv(jdir)
          iph = iphv(jdir)

          ! If emergent
          if (LOS) then

            ct = Geom%L_mu(ith)
            st = sqrt(1d0 - ct*ct)
            cc = cos(Geom%L_phi(iph))
            sc = sin(Geom%L_phi(iph))

          else

            ct = Geom%V_mu(ith)
            st = sqrt(1d0 - ct*ct)
            cc = Geom%v_mux(iph)
            sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

          end if

          ! For each height
          do iz=line%Rz0,Rz1

            ! If not storing, why bother
            if (.not.line%prof(iz,jdir)%VRAM) cycle

            ! Thermal part of the Doppler width
            DwT = sqrt(line%cDopp*line%cDopp*Atmo%T(iz) + &
                       Atmo%vmi(iz)*Atmo%vmi(iz))

            ! Calculate Doppler shift factor
            vfac = 1d0

            if (dyn) &
              vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                           atmo%vz(iz)*ct

            ! Get indexes
            if0 = line%if0
            if1 = line%if1

            ! Get the FS energies
            eu = line%Eu
            el = line%El

            ! Output Doppler width
            Dw = (eu - el)*DwT

            !
            ! Calculate profile
            !

            ! Common quantities
            Dfreqw = (eu - el)/Dw
            vfacw = vfac/Dw
            atuf = line%damp(iz)/Dw

            ! For each frequency
            do ifreq=if0,if1

              call voigtI(Dfreqw - Frec%omega(ifreq)*vfacw, &
                          atuf,prof)

              ! Store
              line%prof(iz,jdir)%p(ifreq) = prof

            end do ! Frequencies

          end do ! output direction
        end do ! height

      end if ! Worker

      ! Control
      call control

      return

      end subroutine getprofI_LTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Validates an existing file with Voigt profiles\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            njdir(integer): Number of directions\n
      !!            lnorm(logical): Flag for file correctness
      subroutine validatevoigtI(Atom,Frec,njdir,lnorm)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Frequency_class), intent(in):: Frec
      logical, intent(out):: lnorm
      integer, intent(in):: njdir

      ! Local

      logical:: nonvalid
      integer:: jdir,jtran,jftran,ktran,nodir
      integer:: if0,if1,iz,iaux,ifreq,ios
      double precision:: d1,d2

      ! Routine name
      urou = 'validatevoigtI'

      ! Get real direction dimension
      if (dyn) then
        nodir = njdir
      else
        nodir = 1
      end if

      ! Check Normp is not allocated
      if (associated(Atom%Normp)) then
        do jdir=1,size(Atom%Normp,3)
          do iz=lbound(Atom%Normp,2),ubound(Atom%Normp,2)
            do jtran=lbound(Atom%Normp,1),ubound(Atom%Normp,1)
              if (allocated(Atom%Normp(jtran,iz,jdir)%prof)) then
                deallocate(Atom%Normp(jtran,iz,jdir)%prof)
              else if (allocated(Atom%Normp(jtran,iz,jdir)%Norm)) then
                deallocate(Atom%Normp(jtran,iz,jdir)%Norm)
              end if
            end do
          end do
        end do
        deallocate(Atom%Normp)
        nullify(Atom%Normp)
      end if

      ! Allocate a token
      allocate(Atom%Normp(1,1,1))


      !
      ! Initialize output
      !
      lnorm = .True.
      nonvalid = .False.


      ! Deallocate sizes
      if (allocated(Atom%zsize)) deallocate(Atom%zsize)
      if (allocated(Atom%dsize)) deallocate(Atom%dsize)
      if (allocated(Atom%tsize)) deallocate(Atom%tsize)
      if (allocated(Atom%f0size)) deallocate(Atom%f0size)

      ! Allocate sizes
      allocate(Atom%zsize(nz))
      allocate(Atom%dsize(nodir))
      allocate(Atom%tsize(Atom%nftran))
      allocate(Atom%f0size(Atom%nftran))

      ! Initialize total blocks of sizes and transition index
      d1 = 0d0
      ktran = 0
      Atom%zsize(1) = 0d0
      Atom%dsize(1) = 0d0
      Atom%tsize(1) = 0d0

      ! For each transition
      do jtran=1,Atom%ntran

        ! For each FS transition
        do jftran=1,Atom%fst(jtran)%nt

          ! Advance index
          ktran = ktran + 1

          ! Skip first
          if (ktran.eq.1) then
            ! Add last
            d2 = dble(Atom%rif1(Atom%ntran) - &
                      Atom%rif0(Atom%ntran) + 1)
            d1 = d1 + d2

            if (Atom%fflag(jtran)%absent) then
              Atom%f0size(ktran) = 0d0
            else
              Atom%f0size(ktran) = 8d0*dble(Atom%if0(jtran) - &
                                            Atom%rif0(jtran))
            end if

            cycle

          end if

          ! If first of the FS transitions
          if (jftran.eq.1) then

            ! Get indexes
            if0 = Atom%rif0(jtran-1)
            if1 = Atom%rif1(jtran-1)

          ! Not the first FS
          else

            ! Get indexes
            if0 = Atom%rif0(jtran)
            if1 = Atom%rif1(jtran)

          end if ! First FS transition

          ! Size of previous
          d2 = dble(if1 - if0 + 1)
          Atom%tsize(ktran) = Atom%tsize(ktran-1) + d2
          d1 = d1 + d2

          ! Frequency size
          if (Atom%fflag(jtran)%absent) then
            Atom%f0size(ktran) = 0d0
          else
            Atom%f0size(ktran) = 8d0*dble(Atom%if0(jtran) - &
                                          Atom%rif0(jtran))
          end if

        end do ! FS transitions
      end do ! Transitions

      ! Add first height to total
      d2 = d1

      ! For each height
      do iz=2,nz

        Atom%zsize(iz) = Atom%zsize(iz-1) + d1
        d2 = d2 + d1

      end do ! Heights

      ! For each direction
      do jdir=2,nodir
        Atom%dsize(jdir) = Atom%dsize(jdir-1) + d2
      end do

      Atom%dsize = Atom%dsize*8d0
      Atom%zsize = Atom%zsize*8d0
      Atom%tsize = Atom%tsize*8d0

      ! Size of header
      Atom%hvifil = 4*4 + & ! Dimension integers
                    8*nodir + & ! Directions sizes
                    8*nz + & ! Height sizes
                    4*2*Atom%nftran + & ! Line limits
                    8*Atom%nftran + & ! Line sizes
                    8*nfreq ! Frequencies

      ! If not the master, wait for diagnostic
      if (pid.gt.0) then

        call control
        call MPI_BCAST(lnorm, 1, MPI_LOGICAL, 0, &
                       MPI_COMM_RT, ierr)
        return

      end if


      ! Only the master gets down here

      !
      ! Check with the file
      !

      do while(.True.)

        ! Open files
        open(200, file=trim(Atom%vfile), status='unknown', &
             iostat=ios, err=1000, access='stream', &
             action='read', form='unformatted')

        ! Read dimensions

        ! Angles
        read(200, err=1100) iaux
        if (iaux.ne.nodir) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Numbers of angles do not coincide:',iaux,nodir
          call verbose
          exit
        end if

        ! Heights
        read(200, err=1100) iaux
        if (iaux.ne.nz) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Numbers of heights do not coincide:',iaux,nz
          call verbose
          exit
        end if

        ! Frequencies
        read(200, err=1100) iaux
        if (iaux.ne.nfreq) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Numbers of frequencies do not coincide:',iaux,nfreq
          call verbose
          exit
        end if

        ! Transitions
        read(200, err=1100) iaux
        if (iaux.ne.Atom%nftran) then
          write(umsg,'(A,1x,i3,1x,i3)') &
            ' # Number of transitions do not coincide:',iaux, &
            Atom%nftran
          call verbose
          exit
        end if

        ! Directional size
        do jdir=1,nodir
          read(200, err=1100) d1
          if (abs(d1-Atom%dsize(jdir)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
              ' # Directional sizes do not coincide:', &
              d1,Atom%dsize(jdir)
            call verbose
            nonvalid = .True.
            exit
          end if
        end do
        if (nonvalid) exit

        ! Height size
        do iz=1,nz
          read(200, err=1100) d1
          if (abs(d1-Atom%zsize(iz)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
              ' # Height sizes do not coincide:',d1,Atom%zsize(iz)
            nonvalid = .True.
            call verbose
            exit
          end if
        end do
        if (nonvalid) exit

        ! Initialize index
        ktran = 0

        ! For each transition
        do jtran=1,Atom%ntran

          ! For each FS transition
          do jftran=1,Atom%fst(jtran)%nt

            ! Advance index
            ktran = ktran + 1

            read(200, err=1100) iaux
            if (iaux.ne.Atom%rif0(jtran)) then
              write(umsg,'(A,1x,i3,1x,i3)') &
                ' # Initial indexes of transition do not coincide:', &
                iaux,Atom%rif0(jtran)
              call verbose
              nonvalid = .True.
              exit
            end if

            read(200, err=1100) iaux
            if (iaux.ne.Atom%rif1(jtran)) then
              write(umsg,'(A,1x,i3,1x,i3)') &
                  ' # Final indexes of transition do not coincide:', &
                  iaux,Atom%rif1(jtran)
              call verbose
              nonvalid = .True.
              exit
            end if

            read(200, err=1100) d1
            if (abs(d1-Atom%tsize(ktran)).gt.TINYO) then
              write(umsg,'(A,1x,es22.15,1x,es22.15)') &
                  ' # Transition sizes do not coincide:', &
                  d1,Atom%tsize(ktran)
              call verbose
              nonvalid = .True.
              exit
            end if

          end do ! FS transitions
          if (nonvalid) exit
        end do ! Transitions
        if (nonvalid) exit

        ! Write frequency
        do ifreq=1,nfreq
          read(200, err=1100) d1
          if (abs(d1-Frec%omega(ifreq)).gt.TINYO) then
            write(umsg,'(A,1x,es22.15,1x,es22.15)') &
                ' # Frequencies do not coincide:',d1,Frec%omega(ifreq)
            call verbose
            nonvalid = .True.
            exit
          end if
        end do
        if (nonvalid) exit

        ! If we are here, file is fine
        lnorm = .False.
        exit

      end do ! Infinite loop

      ! Close file
      close(200, err=1100)

      ! Control
      call control

      ! Share diagnostic
      call MPI_BCAST(lnorm, 1, MPI_LOGICAL, 0, &
                     MPI_COMM_RT, ierr)

      return

1000  umsg = 'Error opening '//trim(Atom%vfile)//' file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1100  umsg = 'Error reading '//trim(Atom%vfile)//' file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine validatevoigtI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Outputs bad normalization data
      subroutine writebadbound(folder,element,iz,jdir,pol,field, &
                               jtran,a1,a2,a3,a4,d1)

      ! I/O
      character(len=2), intent(in):: element
      character(len=500), intent(in):: folder
      logical, intent(in):: pol,field
      integer, intent(in):: iz,jdir,jtran,a1,a2,a3,a4
      double precision, intent(in):: d1

      ! Local
      logical:: exists
      character(LEN=5):: CPUC

      ! Get ID in character
      write(CPUC,'(I0.5)') gpid

      !
      ! ERROR
      inquire(file=trim(folder)//'/badnorm'//CPUC, exist=exists)
      if(.not.exists)then
        open(800,file=trim(folder)//'/badnorm'//CPUC)
      else
        open(800,file=trim(folder)//'/badnorm'//CPUC, &
             position='append')
      endif

      ! 1D case
      if (run_mode.eq.0) then

        ! If polarized case
        if (pol) then

          ! If there is field
          if (field) then

            write(800,'("Atom",1x,A2,1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Components",'// &
                      '4(1x,i3),1x,"Norm",1x,f18.16)') &
              element,iz,jdir,jtran,a1,a2,a3,a4,d1

          ! If no field
          else

            write(800,'("Atom",1x,A2,1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Components",'// &
                      '2(1x,i3),1x,"Norm",1x,f18.16)') &
              element,iz,jdir,jtran,a1,a2,d1

          end if ! Field presence

        ! Intensity case
        else

            write(800,'("Atom",1x,A2,1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Component",'// &
                      '1x,i3,1x,"Norm",1x,f18.16)') &
              element,iz,jdir,jtran,a1,d1

        end if ! Polarization or intensity

      ! Any other case
      else

        ! If polarized case
        if (pol) then

          ! If there is field
          if (field) then

            write(800,'("Atom",1x,A2,1x,'// &
                      '"LOS",1x,"(",i4,",",i4,")",1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Components",'// &
                      '4(1x,i3),1x,"Norm",1x,f18.16)') &
              element,icoords(1:2),iz,jdir,jtran,a1,a2,a3,a4,d1

          ! If no field
          else

            write(800,'("Atom",1x,A2,1x,'// &
                      '"LOS",1x,"(",i4,",",i4,")",1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Components",'// &
                      '2(1x,i3),1x,"Norm",1x,f18.16)') &
              element,icoords(1:2),iz,jdir,jtran,a1,a2,d1

          end if ! Field presence

        ! Intensity case
        else

            write(800,'("Atom",1x,A2,1x,'// &
                      '"LOS",1x,"(",i4,",",i4,")",1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Component",'// &
                      '1x,i3,1x,"Norm",1x,f18.16)') &
              element,icoords(1:2),iz,jdir,jtran,a1,d1

        end if ! Polarization or intensity
      end if ! Pure 1D or 1.5D (rest)

      !
      ! Close
      close(800)

      ! Return
      return

      end subroutine writebadbound

!#####################################################################
!#####################################################################
!#####################################################################

      !> Just make normalization factors unity\n
      !!  Atom(Atom_class(:)): Structure with the atomic data
      subroutine normalize_cle(Atom)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom

      ! Local

      integer:: ia,jtran,iMu,iU,iMf,mf,itermu,itermf
      integer:: nMu,nMf,nL,nU,d1,d2

      double precision:: S,rLu,rLf,rJumax,rJfmax,rMu,rMf


      ! Routine name
      urou = 'normalize_cle'

      ! For each atom
      do ia=1,nA

        ! Allocate vector for norm (part I)
        allocate(Atom(ia)%Normp(Atom(ia)%ntran,1,1))
        ! Allocate sizes for subcomponents for each transition
        if (.not.allocated(Atom(ia)%nL)) &
          allocate(Atom(ia)%nL(Atom(ia)%ntran))
        if (.not.allocated(Atom(ia)%nU)) &
          allocate(Atom(ia)%nU(Atom(ia)%ntran))
        if (.not.allocated(Atom(ia)%nMl)) &
          allocate(Atom(ia)%nMl(Atom(ia)%ntran))
        if (.not.allocated(Atom(ia)%nMu)) &
          allocate(Atom(ia)%nMu(Atom(ia)%ntran))

        ! For each transition
        do jtran=1,Atom(ia)%ntran

          ! Identify the terms
          itermf = Atom(ia)%fst(jtran)%iterml
          itermu = Atom(ia)%fst(jtran)%itermu

          ! Get atomic quantities
          S = Atom(ia)%Sval(itermu)

          rLu = Atom(ia)%rLval(itermu)
          rJumax = rLu+S
          nMu = nint(2d0*rJumax+1d0)

          rLf = Atom(ia)%rLval(itermf)
          rJfmax = rLf + S
          nMf = nint(2d0*rJfmax+1d0)

          ! Fill two of the sizes arrays
          Atom(ia)%nMu(jtran) = nMu
          Atom(ia)%nMl(jtran) = nMf

          ! Initialize counters
          nL = 0
          nU = 0

          ! For each Mu
          do iMu=1,nMu

            rMu = -rJumax + dble(iMu-1)

            ! For each mu_u
            do iU=1,Atom(ia)%nblk(iMu,itermu)

              ! Update counter
              if (iU.gt.nU) nu = iMu

              ! For each Mf
              do iMf=1,nMf

                rMf = -rJfmax + dble(iMf-1)

                ! If not allowed, skip
                if (nint(abs(rMu-rMf)).gt.1) cycle

                ! Sum over mu_f
                do mF=1,Atom(ia)%nblk(iMf,itermf)

                  ! Update counter
                  if (mF.gt.nL) nL = mf

                end do ! iL
              end do ! Ml
            end do ! iU
          end do ! Mu

          ! Update variables
          Atom(ia)%nU(jtran) = nU
          Atom(ia)%nL(jtran) = nL

          !
          ! Allocate the norm array
          !

          d1 = max(Atom(ia)%nJ(itermf),Atom(ia)%nL(jtran))
          d2 = max(Atom(ia)%nJ(itermu),Atom(ia)%nU(jtran))

          ! Allocate
          allocate(Atom(ia)%Normp(jtran,1,1)%Norm(d1,d2, &
                   Atom(ia)%nMl(jtran),Atom(ia)%nMu(jtran)))
          Atom(ia)%Normp(jtran,1,1)%VRAM = .False.
          Atom(ia)%Normp(jtran,1,1)%Norm = 1d0

        end do ! transitions
      end do ! Atoms

      return

      end subroutine normalize_cle

!#####################################################################
!#####################################################################
!#####################################################################

      end module normalizer_mod
