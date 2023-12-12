      !> Solution reader and output writer
      module iosolution_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Contributors:
!     Hao Li (IAC)
!  Start:
!     04/20/2016
!  Last version:
!     09/25/2023 V3.0.16
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:   V3.0.16 - Updated to term- and transition-wise
!                             K cut limits (TdPA)
!                           - Avoid computing tensor components
!                             with negative Q values (TdPA)
!
!     09/25/2023:   V3.0.15 - Change names for population and
!                             departure coefficient files (TdPA)
!
!     08/07/2023:   V3.0.14 - readsol now returns if the Stokes
!                             parameters could be read (TdPA)
!                           - Use czero to initialize instead of the
!                             explicit 0+j0 (TdPA)
!
!     07/05/2023:   V3.0.13 - Bugfix: There was no writing of Stokes
!                             parameters in the inversion if it was
!                             intensity only (TdPA)
!                           - Bugfix: Forgot to code in the writing
!                             of the contribution function in the
!                             inversion when the height range is
!                             not cropped (TdPA)
!
!     07/03/2023:   V3.0.12 - Added getsol, setsol, setstk, setstkI,
!                             setctr, setctrI, and settau to manage
!                             the input/output of the solution of
!                             the forward problem when running
!                             inversions (TdPA)
!                           - Added an explicit branch for the
!                             inversion when writing the emerging
!                             Stokes parameters (TdPA)
!                           - Added writectrI_inv, writectr_inv, and
!                             writetau_inv to store in file the
!                             contribution function and height of
!                             optical depth equal to one for the
!                             inversion (TdPA)
!                           - Bugfix: Wrong upper limit for Stokes in
!                             readsol when not reading the full
!                             Stokes array in the polarization
!                             problem (TdPA)
!
!     04/25/2023:   V3.0.11 - Bugfix: When saving the intensity
!                             solution file, the definition of AV_int
!                             was the opposite as it should be (TdPA)
!
!     03/21/2023:   V3.0.10 - Bugfix: When saving the intensity
!                             solution file, the branches to store
!                             Stokes or J00C were interchanged (TdPA)
!
!     03/08/2023:    V3.0.9 - When readsol fails to open or read the
!                             indicated input file, it writes the
!                             path in the error message (TdPA)
!
!     02/14/2023:    V3.0.8 - The inversion also calls readsol,
!                             writesol, writesolI, writestk,
!                             writestokI, writectr, writectrI,
!                             writetau, and writeatmo  so it had to
!                             be taken into consideration for
!                             selecting the filename (HL)
!                           - The total population is not updated
!                             in readsol for the inversion mode (HL)
!                           - Readsol accounts for the split in
!                             quadratures (TdPA)
!                           - writesolI uses axiali (TdPA)
!                           - Added GeomI exclusive for use in the
!                             intensity problem (TdPA)
!
!     11/24/2022:    V3.0.7 - Added write_CLEgeom and write_CLE
!                             subroutines (TdPA)
!
!     11/10/2022:    V3.0.6 - Bugfix: The populations of the regions
!                             outside of the RT height range are
!                             not normalized (TdPA)
!                           - Bugfix: There was a departure
!                             coefficient writing that was directed
!                             to the population file unit (TdPA)
!
!     10/26/2022:    V3.0.5 - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.4 - Implemented the height axis restriction
!                             option (TdPA)
!                           - Changed unit numbers in writesol and
!                             writesolI so they are unique (TdPA)
!                           - Bugfix: The offset when writing the
!                             emergent intensity (intensity only case)
!                             had jumps 4 times larger than it should.
!                             It also wrote 16 times more characters
!                             than it should. Fixed both sizes (TdPA)
!
!     07/27/2022:    V3.0.3 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!                           - funit is now global (TdPA)
!
!     07/08/2022:    V3.0.2 - Bugfix: When the output is limited, in
!                             the contribution function the order
!                             height-frequency must be imposed with
!                             an unoptimal additional loop (TdPA)
!                           - Bugfix: Collisions between a term with
!                             itself, or a level with itself, should
!                             not be filtered out in the output (TdPA)
!                           - Bugfix: tau was not being stored in the
!                             atmospheric file when height was the
!                             input (TdPA)
!                           - The atmospheric model is now stored in
!                             double precision (TdPA)
!                           - Changed the name of the solution folder
!                             for the 1.5D case (TdPA)
!
!     06/30/2022:    V3.0.1 - Bugfixes: Quick changes while preparing
!                             the last commit led to mistakes. Solved
!                             some typos, a missing allocation, and
!                             a wrong label in goto (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o readsol, writesol, and writesolI now
!                                choose the solution file path and
!                                filename depending on the run type.
!                              o The first two inputs in writesol and
!                                writesolI have been changed to the
!                                Input structure.
!                              o Writing the rhoKQ, JKQ, Stokesout,
!                                .pop, and .dep files, and even the
!                                solution file is technically optional
!                                now.
!                              o Added returns for when we need to
!                                abort the run.
!                              o Changed the MPI communicator from
!                                MPI_COMM_WORLD to MPI_COMM_RT.
!                              o Added an additional input for
!                                writestk, writestkI, writectr,
!                                writectrI, writetau, writecols,
!                                writedamp, writeback, and writeatmo
!                                so they can handle 1.5D writing.
!                             (TdPA)
!
!     04/07/2022:    V2.0.3 - Added the possibility to write the
!                             frequency dependent radiation field
!                             tensors as part of saving the solution
!                             file (TdPA)
!
!     09/30/2021:    V2.0.2 - Needed to change behaviour when reading
!                             Stokes in an axial problem while keeping
!                             Stokes parameters (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/18/2021:    V2.0.1 - Removed redundant call to control in
!                             writeback (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Removed many control calls (TdPA)
!
!     01/12/2021:   V1.7.14 - Use KSTK to deal with the radiation
!                             field part when reading and writing
!                             the solution file (TdPA)
!
!     09/11/2020:   V1.7.13 - Added RAM counters for radiation field
!                             quantities (TdPA)
!                           - Store departure coefficients in its
!                             own file (TdPA)
!
!     07/10/2020:   V1.7.12 - Bugfix: Error in allocation of Stokes
!                             when reading an angle-averaged solution
!                             with PRD and dynamics (TdPA)
!
!     05/11/2020:   V1.7.11 - Bugfix: In readsol, the master was
!                             sending the wrong amount of elements
!                             for the populations and there was an
!                             extra argument in the isend call for
!                             the density matrix, both when there was
!                             domain decomposition (TdPA)
!                           - The writeback routine now also writes
!                             the frequency vector in the background
!                             quantities file (TdPA)
!
!     04/14/2020:   V1.7.10 - Bugfix: Since 11/13/2019, 1.11.0 version
!                             of rtcoeffiaux_mod, it was necessary to
!                             initialize the Atom(:)%popu when reading
!                             the density matrix Atom(:)%crho (TdPA)
!
!     03/05/2020:    V1.7.9 - Avoided extra space when saving file
!                             with populations for elements with a
!                             single letter (TdPA)
!                           - Changed format of atmos.dat, added
!                             column for atomic H density (TdPA)
!                           - Added routine wAtmo (TdPA)
!
!     12/17/2019:    V1.7.8 - Added continuum opacity to the output
!                             atmosphere (TdPA)
!                           - Fixed header of atmospheric file (TdPA)
!
!     12/10/2019:    V1.7.7 - Completed the writeatmo routine (TdPA)
!
!     11/19/2019:    V1.7.6 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     10/31/2019:    V1.7.5 - J00S need to be allocated even if not
!                             used due to the term correction (TdPA)
!
!     09/26/2019:    V1.7.4 - Added writeatmo (TdPA)
!
!     09/13/2019:    V1.7.3 - Changed condition to write and read
!                             Stokes instead of JKQ, as the angle-
!                             average case with velocities requires
!                             the former now (TdPA)
!
!     08/19/2019:    V1.7.2 - Ignoring the reading of J00S. Commented
!                             the reading and the communication among
!                             CPU (TdPA)
!                           - Bugfix: Wrong message when the number of
!                             polar directions in the run and in the
!                             read solution file were different (TdPA)
!
!     08/16/2019:    V1.7.1 - Bugfix: was using a bad identifier for
!                             MPI without domain decomposition, was
!                             checking MPI%nsend.gt.0 instead of
!                             MPI%mpi (TdPA)
!
!     05/08/2019:    V1.7.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!
!     04/16/2019:    V1.6.2 - Introduced back the normal broadcasting
!                             from mpi (TdPA)
!
!     03/12/2019:    V1.6.1 - Bugfix: Did not change two of the write
!                             statements to the new units (TdPA)
!
!     02/20/2019:    V1.6.0 - New verbosity (TdPA)
!                           - Now uses units in the hundreds (TdPA)
!                           - Now uses specific TINY variables (TdPA)
!
!     01/23/2019:    V1.5.2 - Read JKQC when the solution has
!                             different number of angles if it was
!                             angle-averaged (TdPA)
!
!     11/19/2018:    V1.5.1 - Writedamp prints full damping for each
!                             transition (TdPA)
!                           - Added identificators for collision,
!                             damping, and background files (TdPA)
!
!     11/06/2018:    V1.5.0 - Added writecols, writedamp, and
!                             writeback (TdPA)
!
!     10/01/2018:    V1.4.6 - When storing the intensity solution, it
!                             also saves JoutI, RhooutI, and
!                             StokesoutI (TdPA)
!
!     09/06/2018:    V1.4.5 - Bugfix: When not doing PRD it did not
!                             make sense to check the AV variable. Now
!                             AV at writing and inAV at reading depend
!                             on both booleans (TdPA)
!                           - Possibility to write the intensity
!                             solution with a different name so the
!                             polarization solution does not overwrite
!                             it (TdPA)
!
!     07/30/2018:    V1.4.4 - Bugfix: When reading and applying cut on
!                             K do not try to rotate rhoKQ that does
!                             not exist (TdPA)
!
!     06/01/2018:    V1.4.3 - Do not display warning about the way of
!                             reading Stokes parameters if not going
!                             to read them anyways (TdPA)
!
!     05/31/2018:    V1.4.2 - Missing space in warning about reading
!                             non-axial solution as initialization of
!                             an axial problem (TdPA)
!                           - Bugfix: When the problem was axial and
!                             the solution was not, the reading of
!                             Stokes parameters was wrong, because it
!                             did not read the full azimuth and
!                             assigned the wrong Stokes for different
!                             polar angles (TdPA)
!
!     05/16/2018:    V1.4.1 - Stokes has dimension nPh, instead of
!                             nPh2, in the azimuthal index (TdPA)
!                           - Modified readsol to account for that
!                             index change (TdPA)
!
!     10/24/2017:    V1.4.0 - Writes the population file too (TdPA)
!
!     10/09/2017:    V1.3.1 - Forgot to share the total population for
!                             each atom (TdPA)
!
!     09/27/2017:    V1.3.0 - Implemented send_tree algorithm to
!                             avoid the terrible scaling of the native
!                             mpi_bcast (TdPA)
!
!     09/22/2017:    V1.2.0 - Possibility to limit K (TdPA)
!                           - Bugfix: In writesolI, was not writting
!                             J, was not defined (TdPA)
!
!     09/13/2017:    V1.1.2 - Bugfix: rhonull for population when
!                             reading intensity solution was set to
!                             .True. by comparing with TINY15,
!                             instead with 0 (TdPA)
!
!     09/12/2017:    V1.1.1 - Bugfix: The writing of rhooutI was
!                             wrong (TdPA)
!                           - Not writing logical rhonull, but a proxy
!                             integer (TdPA)
!
!     08/04/2017:    V1.1.0 - Only the master reads, the others just
!                             recieve (TdPA)
!
!     07/06/2017:    V1.0.4 - Bugfix: Because dimensions are not
!                             specified, the Stokes look for the
!                             contribution function must go from 1 to
!                             4, instead of from 0 to 3 (TdPA)
!
!     06/29/2017:    V1.0.3 - Bugfix: When reading from intensity
!                             solution every rhonull was made True
!                             (TdPA)
!
!     06/22/2017:    V1.0.2 - Redone writesolI (TdPA)
!                           - Now the code can read multilevel or
!                             multiterm solutions (TdPA)
!
!     05/26/2017:    V1.0.1 - Bugfix: There was a iz limited by the
!                             CPU limits in readsol, it should take
!                             all the nodes (TdPA)
!                           - Bugfix: JKQC was only being stored if
!                             there was stimulated emission, a copy
!                             and paste bug (TdPA)
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
!  readsol
!    Restores rhoKQ, JKQ and Stokes from files
!
!  getsol
!    Restores rhoKQ, JKQ and Stokes from RAM
!
!  setsol
!    Stores rhoKQ, JKQ and Stokes in RAM
!
!  writesol
!    Stores rhoKQ, JKQ and Stokes in file
!
!  writesolI
!    Stores rhoKQ, JKQ and Stokes from solverI in file
!
!  writestk
!    Stores emergent Stokes in file
!
!  writestkI
!    Stores emergent Stokes from emergenceI in file
!
!  write_CLEgeom
!    Write the geometry of a LOS in CLE in file
!
!  write_CLE
!    Write emergent Stokes and optical depth in CLE in file
!
!  writectr
!    Stores contribution function in file
!
!  writectr_inv
!    Stores contribution function from inversion in file
!
!  writectrI
!    Stores contribution function from emergenceI in file
!
!  writectrI_inv
!    Stores contribution function from intensity inversion in file
!
!  setctr
!    Stores contribution function in RAM
!
!  setctrI
!    Stores contribution function for intensity in RAM
!
!  writetau
!    Stores height where optical depth is equal to one in file
!
!  writetau_inv
!    Stores height where optical depth is equal to one from inversion
!  in file
!
!  settau
!    Stores height where optical depth is equal to one in RAM
!  in file
!
!  writecols
!    Stores inelastic collisions of active atoms in file
!
!  writedamp
!    Stores damping parameter of active atoms in file
!
!  writeback
!    Stores background continuum quantities in file
!
!  writeatmo
!    Stores atmospheric data in ASCII file atmos.dat
!
!  wAtmo
!    Writes an updated atmospheric file ready to be used by HanleRT
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use parameters_mod , only : pi , TINYB , TINYR , cZero , c , me
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Reads a file with a solution calculated in a previous run of
      !! the code.\n
      !!    filename(character(:)): Name of the file to read\n
      !!     GeomI(Geometry_class): Structure with geometry data for
      !!                            the intensity problem\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine readsol(filename,GeomI,Geom,MPID,Flgsg,Bfield,Atom, &
                         read_stokes, &
                         Stokes,JKQ,JKQS,JKQC, &
                         Stokes0,J00,J00S,J00C,J00P)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Geometry_class), intent(in):: Geom,GeomI
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(MPI_class), intent(inout):: MPID
      character(len=500), intent(in):: filename
      logical, intent(out):: read_stokes
      double precision,dimension(:,:,:,:), allocatable:: Stokes0
      double precision,dimension(:,:,:,:,:), allocatable:: Stokes
      double precision,dimension(:,:), allocatable:: J00
      double precision,dimension(:,:), allocatable:: J00S
      double precision,dimension(:,:), allocatable:: J00C
      double precision,dimension(:,:,:), allocatable:: J00P
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      character(len=2):: label
      character(len=8):: scoord

      logical:: inaxial,instm,inAV

      integer:: ios,ia,iz,ifreq,i,ith,iph,iS,itran,jtran,istep
      integer:: it,iJ,iJ1,K,iQ,iR,iR0,iR1
      integer:: ia1,ia2,ia3,ia4,ia5
      integer:: ilabel,inaxial_int,instm_int,inAV_int
      integer:: iproc,rsize,jsize,psize,csize,ssize,nsize
      integer, dimension(:), allocatable:: ibuff

      double precision:: rJ,rJ1,da1,da2,rho0

      complex(kind=8):: integr
      complex(kind=8),dimension(-2:2,0:2,nz):: JKQaux
      complex(kind=8),dimension(-nkx:nkx,nz):: rhoKQaux

      ! Routine name
      urou = 'readsol'

      !
      ! Open file
      !
      if (pid.eq.0) then

        !
        ! If 1D or inversion
        !
        if (run_mode.le.0) then

          ! Open file
          open (200,file=trim(filename), &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='read', form='unformatted')
        !
        ! If 1.5D
        !
        else if (run_mode.eq.1) then

          ! Get LOS index
          write(scoord,'(I0.8)') icoords(3)

          ! Open file
          open (200,file=trim(filename)//'/Solution-'//scoord, &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='read', form='unformatted')

        end if
      end if ! Master


      !
      ! Dimensions, can be use for error handling
      !
      if (pid.eq.0) then
        read (200,err=1100) label
        if (label.eq.'sp') then
          ilabel = 0
        else if (label.eq.'si') then
          ilabel = 1
        else
          ilabel = -1
        end if
      end if

      ! If MPI
      if (MPID%mpi) then

        ! Control
        call control
        if (laborted) return

        ! Alternative bcast
        if (MPID%altbcast) then

          ! If not master, receive first
          if (pid.ne.0) then

            ! Receive label
            call MPI_RECV(ilabel, 1, MPI_INTEGER, &
                          MPID%recv, 1000000+pid, &
                          MPI_COMM_RT, MPI_STATUS_IGNORE, ierr)

          end if ! No Master

          ! For each send
          do istep=1,MPID%nsend

            ! Send label
            call MPI_ISEND(ilabel, 1, MPI_INTEGER, &
                           MPID%lsend(istep), &
                           1000000+MPID%lsend(istep), &
                           MPI_COMM_RT, MPID%requestA(istep,6), &
                           ierr)
          end do ! sends

          ! For each slave to send
          do iproc=1,MPID%nsend

            ! Wait for everyone to receive the radiation data before
            ! continuing
            call MPI_WAIT(MPID%requestA(iproc,6), &
                          MPI_STATUS_IGNORE,ierr)

          end do

        ! Normal bcast
        else

          call MPI_BCAST(ilabel, 1, MPI_INTEGER, 0, MPI_COMM_RT, &
                         ierr)

        end if ! Type of bcast

      end if ! MPI

      ! Important to check the file label
      if(ilabel.ne.0.and.ilabel.ne.1) then
        umsg = 'The specified solution file does '// &
               'not have the correct ID of a solution file'
        call aborted
        return
      end if

      !
      ! Allocate J00 for photoionizations (common for both)
      !
      allocate(J00P(nxphot,2,Rz0:Rz1))
      MPID%RRAM = 8d-6*dble(2*nxphot*Rnz)

      !
      ! Polarization read
      !
      if (ilabel.eq.0) then

        ! Initialize J00P
        J00P = 0d0

        !
        ! Allocations
        !

        ! If we are doing angle averaged, we only need one height,
        ! we allocate two to store the emergence in the quadrature
        if (KSTK) then
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1
          MPID%RRAM = MPID%RRAM + &
                      8d-6*dble(4*nfreq*Geom%nph*Geom%nth*Rnz)
        else
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1
          MPID%RRAM = MPID%RRAM + &
                      8d-6*dble(4*nfreq*Geom%nph*Geom%nth*2)
        end if

        ! JKQ for absorptivity
        allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))

        ! JKQ for stimulated emission
        allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))

        ! JKQ frequency dependent
        allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))

        ! Allocated memory
        MPID%RRAM = MPID%RRAM + &
                    8d-6*dble(Rnz*2*5*3*(2*nxtran + nfreq))

        ! If the master
        if (pid.eq.0) then
          read (200,err=1100) ia1
          read (200,err=1100) ia2
          read (200,err=1100) ia3
          read (200,err=1100) ia4
          read (200,err=1100) ia5
          read (200,err=1100) inaxial_int
          read (200,err=1100) instm_int
          read (200,err=1100) inAV_int
        end if

        ! If doing MPI
        if (MPID%mpi) then

          allocate(ibuff(8))

          ! Control
          call control
          if (laborted) return

          ! If the master
          if (pid.eq.0) then
            ibuff(1) = ia1
            ibuff(2) = ia2
            ibuff(3) = ia3
            ibuff(4) = ia4
            ibuff(5) = ia5
            ibuff(6) = inaxial_int
            ibuff(7) = instm_int
            ibuff(8) = inAV_int
          end if

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive dimensions
              call MPI_RECV(ibuff(1), 8, MPI_INTEGER, &
                            MPID%recv, 2000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No Master

            ! For each send
            do istep=1,MPID%nsend

              ! Send dimensions
              call MPI_ISEND(ibuff(1), 8, MPI_INTEGER, &
                             MPID%lsend(istep), &
                             2000000+MPID%lsend(istep), &
                             MPI_COMM_RT, MPID%requestA(istep,7), &
                             ierr)

            end do ! sends

          ! Normal bcast
          else

            call MPI_BCAST(ibuff(1), 8, MPI_INTEGER, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast

          ! If a slave
          if (pid.ne.0) then
            ia1 = ibuff(1)
            ia2 = ibuff(2)
            ia3 = ibuff(3)
            ia4 = ibuff(4)
            ia5 = ibuff(5)
            inaxial_int = ibuff(6)
            instm_int = ibuff(7)
            inAV_int = ibuff(8)
          end if ! slave

          ! If alternative bcast
          if (MPID%altbcast) then

            ! For each slave to send
            do iproc=1,MPID%nsend

              ! Wait for everyone to receive the radiation data before
              ! continuing
              call MPI_WAIT(MPID%requestA(iproc,7), &
                            MPI_STATUS_IGNORE,ierr)

            end do

          end if ! Alternative bcast
        end if ! MPI

        ! Convert to logical
        if(inaxial_int.eq.1)then
          inaxial = .True.
        else
          inaxial = .False.
        end if

        ! Convert to logical
        if(instm_int.eq.1)then
          instm = .True.
        else
          instm = .False.
        end if

        ! Convert to logical
        if(inAV_int.eq.1)then
          inAV = .True.
        else
          inAV = .False.
        end if

        ! Flag to read stokes parameter
        read_stokes = .True.

        !
        ! Dimension checking
        !

        ! Check height nodes
        if (ia2.ne.nz) then
          umsg = 'Solution file with different number of '// &
                 'heights.'
          call aborted
          return
        end if

        ! Check number of atoms
        if (ia5.ne.nA) then
          umsg = 'Solution file with different number of '// &
                 'atoms.'
          call aborted
          return
        end if

        ! Check number of frequencies, this one does not produce an
        ! abortion
        if (ia1.ne.nfreq) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of frequencies in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes and J^K_Q(nu).'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Check polar nodes, this one does not produce an abortion
        if (ia3.ne.Geom%nTh.and.read_stokes.and..not.inAV) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of polar directions in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes.'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Check azimuthal nodes, this one does not produce an abortion
        if (ia4.ne.Geom%nPh.and.(ia4.ne.1.and.Geom%nPh.ne.1).and. &
            read_stokes.and..not.inAV) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'and is not axial, ignoring Stokes.'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Check azimuthal nodes for AD redistribution, this one does
        ! not produce an abortion
        if (ia4.ne.Geom%nPh2.and.(ia4.ne.1).and. &
            read_stokes.and..not.inAV) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'for eps^(2) and is not axial, ignoring Stokes.'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Warning when reading non axial from AD for axial AD
        if (.not.inaxial.and.axial.and.(.not.(AV.and..not.dyn)).and. &
            .not.inAV.and.read_stokes.and.pid.eq.0) then
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read.'
          call verbose
        end if

        ! Warning when reading non axial from AD for axial AV
        if (.not.inaxial.and.axial.and.(AV.and..not.dyn).and. &
            .not.inAV.and.read_stokes.and.pid.eq.0) then
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read to compute JKQC.'
          call verbose
        end if


        !
        ! Population and RhoKQ
        !

        ! For each atom
        do ia=1,nA

          ! If alternative bcast
          if (MPID%mpi.and.MPID%altbcast) then

            ! For each slave
            do iproc=1,MPID%nsend

              ! Wait for everyone to receive the radiation data before
              ! continuing and reset the buffers
              call MPI_WAIT(MPID%requestA(iproc,1), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,6), &
                            MPI_STATUS_IGNORE,ierr)

            end do ! Processors

          end if ! Alternative Bcast

          ! Only the master reads
          if (pid.eq.0) then

            ! For each height
            do iz=1,nZ

              ! Read population at this height
              read (200,err=1100) da1
              if (run_mode.ge.0) Atom(ia)%n(iz) = da1

            end do ! heights

            ! For each terms
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)

                ! For each level
                do iJ1=1,Atom(ia)%nJ(it)

                  ! Get J'
                  rJ1 = Atom(ia)%rJval(iJ1,it)

                  ! For each K
                  do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                    ! For each Q
                    do iQ=-K,K

                      ! Get rho(J,J')KQ index
                      if (K.le.Atom(ia)%Kcut(it)) &
                        iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! For each height
                      do iz=1,nZ

                        ! Read real and imaginary parts
                        read (200,err=1100) da1,da2

                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        if (K.le.Atom(ia)%Kcut(it)) then

                          Atom(ia)%crho(iR,iz) = dcmplx(da1,da2)

                          ! Auxiliar variable for rotation
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end if

                      end do ! heights

                    end do ! Q

                    !
                    ! Rotate rhoKQ
                    !

                    ! Only if above the limit
                    if (K.le.Atom(ia)%Kcut(it)) then

                      ! For each height in the CPU domain
                      do iz=Rz0,Rz1

                        ! If there is magnetic field
                        if (Bfield%Bstrength(iz).gt.TINYB) then

                          ! Rotate the rhoKQ in the auxiliar variable
                          call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                    Bfield%Btheta(iz), &
                                    Bfield%Bphi(iz),1)

                          ! Store the rotated result in the rhoKQ
                          ! array
                          do iQ=-K,K
                            iR = Atom(ia)%irho(it)% &
                                          Jrho(iJ1,iJ)%kq(iQ,K)
                            Atom(ia)%crho(iR,iz) = rhoKQaux(iQ,iz)
                          end do ! Q

                        end if ! There is B field

                      end do ! heights

                    end if ! K < Kcut

                  end do ! K
                end do ! J'
              end do ! J
            end do ! terms

          end if ! Master


          ! If there are slaves
          if (MPID%mpi) then

            ! Control
            call control
            if (laborted) return

            nsize = nZ

            rsize = Atom(ia)%ndim*RnZ

            ! Alternative bcast
            if (MPID%altbcast) then

              ! If not master, receive first
              if (pid.ne.0) then

                ! Receive n
                call MPI_RECV(Atom(ia)%n(1), nsize, &
                              MPI_DOUBLE_PRECISION, &
                              MPID%recv, 2000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

                ! Receive rhoKQ
                call MPI_RECV(Atom(ia)%crho(1,Rz0), rsize, &
                              MPI_DOUBLE_COMPLEX, &
                              MPID%recv, 3000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No Master

              ! For each send
              do istep=1,MPID%nsend

                ! Send n
                call MPI_ISEND(Atom(ia)%n(1), nsize, &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               2000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,6), &
                               ierr)

                ! Send goout
                call MPI_ISEND(Atom(ia)%crho(1,Rz0), rsize, &
                               MPI_DOUBLE_COMPLEX, &
                               MPID%lsend(istep), &
                               3000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,1), &
                               ierr)

              end do ! sends

            ! Normal bcast
            else

              ! Send n
              call MPI_BCAST(Atom(ia)%n(1), nsize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

              ! Send rho
              call MPI_BCAST(Atom(ia)%crho(1,Rz0), rsize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)

            end if ! Type of bcast
          end if ! MPI

        end do ! atoms

        !
        ! Flag the null rho(J,J')KQ
        !

        ! For each atom
        do ia=1,nA

          ! Initialize to non-null
          Atom(ia)%rhonull = .False.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)

              ! Get the rho00 indexes
              iR0 = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! For each level
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J'
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! Get the rho00 indexes
                iR1 = Atom(ia)%irho(it)%Jrho(iJ1,iJ1)%kq(0,0)

                ! For each K
                do K=nint(abs(rJ-rJ1)), &
                     min(nint(rJ+rJ1),Atom(ia)%Kcut(it))

                  ! For each Q
                  do iQ=-K,K

                    ! Get rho(J,J')KQ index
                    iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                    ! For each height in this domain
                    do iz=Rz0,Rz1

                      ! Get the inverse of rho00
                      rho0 = 1d0/sqrt(abs(Atom(ia)%crho(iR0,iz))* &
                                      abs(Atom(ia)%crho(iR1,iz)))

                      ! If rhoKQ/rho00 is lesser than double precision
                      ! flag null
                      if (abs(Atom(ia)%crho(iR,iz)*rho0).lt.TINYR) &
                        Atom(ia)%rhonull(iR,iz) = .True.

                    end do ! heights

                  end do ! Q
                end do ! K
              end do ! J'
            end do ! J
          end do ! terms
        end do ! atoms


        !
        ! JKQ
        !

        ! Only the master reads
        if (pid.eq.0) then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Atomic shift
              jtran = itran + Atom(ia)%tshift

              ! For each K
              do K=0,2

                ! In K cut
                if (K.le.Atom(ia)%Krad(itran)) then

                  ! For each Q
                  do iQ=-K,K

                    ! For each height
                    do iz=1,nZ

                      ! Read real and imaginary parts of JKQ
                      read (200,err=1100) da1,da2
                      JKQaux(iQ,K,iz) = dcmplx(da1,da2)

                    end do ! heights
                  end do ! Q

                ! Out of cut
                else

                  ! Jump
                  iQ = (2*K+1)*nZ*16
                  call fseek(200,iQ,1)
                  JKQaux(:,K,:) = cZero

                end if

              end do ! K

              !
              ! Rotate JKQ
              !

              ! For each height
              do iz=Rz0,Rz1

                ! Register JKQ in the array
                JKQ(:,:,jtran,iz) = JKQaux(:,:,iz)

                ! If there is magnetic field, rotate
                if (Bfield%Bstrength(iz).gt.TINYB) &
                call fieldB(JKQ(:,:,jtran,iz),1,Flgsg, &
                            Bfield%Btheta(iz),Bfield%Bphi(iz),1)

              end do ! heights
            end do ! transitions
          end do ! atoms

        end if ! Master

        ! If there are slaves
        if (MPID%mpi) then

          ! Control
          call control
          if (laborted) return

          jsize = 15*nxtran*RnZ

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive JKQ
              call MPI_RECV(JKQ(-2,0,1,Rz0), jsize, &
                            MPI_DOUBLE_COMPLEX,  &
                            MPID%recv, 4000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No Master

            ! For each send
            do istep=1,MPID%nsend

              ! Send JKQ
              call MPI_ISEND(JKQ(-2,0,1,Rz0), jsize, &
                             MPI_DOUBLE_COMPLEX, &
                             MPID%lsend(istep), &
                             4000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,2), ierr)

            end do ! Sends

          ! Normal bcast
          else

            ! Share JKQ
            call MPI_BCAST(JKQ(-2,0,1,Rz0), jsize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast
        end if ! MPI

        ! If there is stimulated emission in the input
        if(instm)then

          ! Only the master reads
          if (pid.eq.0) then

            ! For each atom
            do ia=1,nA

              ! For each transition
              do itran=1,Atom(ia)%ntran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tshift

                ! For each K
                do K=0,2

                  ! In K cut
                  if (K.le.Atom(ia)%Krad(itran)) then

                    ! For each Q
                    do iQ=-K,K

                      ! For each height
                      do iz=1,nZ

                        ! Read real and imaginary parts of JKQS
                        read (200,err=1100) da1,da2

                        ! If in this process domain, register
                        if (stm) JKQaux(iQ,K,iz) = dcmplx(da1,da2)


                      end do ! heights
                    end do ! Q

                  ! Out of cut
                  else

                    ! Jump
                    iQ = (2*K+1)*nZ*16
                    call fseek(200,iQ,1)
                    JKQaux(:,K,:) = cZero

                  end if

                end do ! K

                !
                ! Rotate JKQ
                !

                ! If we are currently doing stimulated emission
                if(stm)then

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Register JKQS in the array
                    JKQS(:,:,jtran,iz) = JKQaux(:,:,iz)

                    ! If there is magnetic field, rotate
                    if (Bfield%Bstrength(iz).gt.TINYB) &
                    call fieldB(JKQS(:,:,jtran,iz),1,Flgsg, &
                                Bfield%Btheta(iz),Bfield%Bphi(iz),1)

                  end do ! heights
                end if ! stimulated emission
              end do ! transitions
            end do ! atoms

          end if ! Master

          ! If there are slaves
          if (MPID%mpi) then

            ! Control
            call control
            if (laborted) return

            jsize = 15*nxtran*RnZ

            ! Alternative bcast
            if (MPID%altbcast) then

              ! If not master, receive first
              if (pid.ne.0) then

                ! Receive JKQ
                call MPI_RECV(JKQS(-2,0,1,Rz0), jsize, &
                              MPI_DOUBLE_COMPLEX,  &
                              MPID%recv, 5000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No Master

              ! For each send
              do istep=1,MPID%nsend

                ! Send JKQ
                call MPI_ISEND(JKQS(-2,0,1,Rz0), jsize, &
                               MPI_DOUBLE_COMPLEX, &
                               MPID%lsend(istep), &
                               5000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,3), ierr)

              end do ! Sends

            ! Normal bcast
            else

              call MPI_BCAST(JKQS(-2,0,1,Rz0), jsize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)

            end if ! Type of bcast
          end if ! MPI

        ! No stimulated emission in the input
        else

          if(stm)JKQS = JKQ

        end if ! Stimulated emission in the input


        !
        ! Stokes or JKQC
        !

        ! If we can read Stokes
        if (read_stokes) then

          ! If input is only JKQ
          if (inAV) then

            ! Only the master reads
            if (pid.eq.0) then

              ! For each height
              do iz=1,nz

                ! For each frequency
                do ifreq=1,nfreq

                  ! For each K
                  do K=0,2

                    ! If within K cut
                    if (K.le.Krad) then

                      ! For each Q
                      do iQ=-K,K

                        ! Read real and imaginary parts of JKQS
                        read (200,err=1100) da1,da2
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                        JKQC(iQ,K,ifreq,iz) = dcmplx(da1,da2)

                      end do ! Q

                    ! Out of cut
                    else

                      ! Jump
                      iQ = (2*K+1)*16
                      call fseek(200,iQ,1)
                      JKQC(:,K,ifreq,iz) = cZero

                    end if

                  end do ! K
                end do ! frequency
              end do ! heights

            end if ! Master

            ! If there are slaves
            if (MPID%mpi) then

              ! Control
              call control
              if (laborted) return

              csize = 15*nfreq*RnZ

              ! Alternative bcast
              if (MPID%altbcast) then

                ! If not master, receive first
                if (pid.ne.0) then

                  ! Receive JKQC
                  call MPI_RECV(JKQC(-2,0,1,Rz0), csize, &
                                MPI_DOUBLE_COMPLEX, &
                                MPID%recv, 6000000+pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                end if ! No Master

                ! For each send
                do istep=1,MPID%nsend

                  ! Send JKQC
                  call MPI_ISEND(JKQC(-2,0,1,Rz0), csize, &
                                 MPI_DOUBLE_COMPLEX, &
                                 MPID%lsend(istep), &
                                 6000000+MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,4), ierr)

                end do ! Sends

              ! Normal bcast
              else

                call MPI_BCAST(JKQC(-2,0,1,Rz0), csize, &
                               MPI_DOUBLE_COMPLEX, 0, &
                               MPI_COMM_RT, ierr)

              end if
            end if ! MPI

            ! If we are doing angle dependent, the first iteration
            ! must be angle averaged if the input is angle average
            if (.not.AV) tbAD = .True.

          ! If input is AD
          else

            !
            ! Currently doing AV or not PRD, no need to read Stokes
            ! Therefore, no need to keep Stokes
            !
            if (.not.KSTK.or..not.(dyn.or..not.AV)) then

              ! Master
              if(pid.eq.0)then

                JKQC = dcmplx(.0D0,.0D0)

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! For each Stokes parameter
                        do iS=0,3

                          ! If the input is not axial or is the first
                          ! azimuth
                          if(.not.inaxial.or.iph.eq.1)then

                            ! Read the data point
                            read (200,err=1100) da1

                            if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                            if (.not.axial.or.iph.eq.1) then

                              ! Store
                              if (KSTK) then
                                Stokes(iS,ifreq,iph,ith,iz) = da1
                              else
                                Stokes(iS,ifreq,iph,ith,Rz0) = da1
                              end if

                            end if

                          ! If the input is axial
                          else

                            if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                            ! Store
                            if (KSTK) then
                              Stokes(iS,ifreq,iph,ith,iz) = &
                                            Stokes(iS,ifreq,1,ith,iz)
                            else
                              Stokes(iS,ifreq,iph,ith,1) = &
                                            Stokes(iS,ifreq,1,ith,1)
                            end if

                          end if ! axial input

                        end do ! Stokes parameters

                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! For each K
                        do K=0,Krad

                          ! For each Q
                          do iQ=0,K

                            ! Sum over Stokes parameters of Stokes*TKQ
                            integr = sum(Stokes(:,ifreq,iph,ith,1)* &
                                         Geom%TS(:,iQ,K,iph,ith))

                            ! Add contribution to the JKQC integral
                            JKQC(iQ,K,ifreq,iz) = &
                                               JKQC(iQ,K,ifreq,iz) + &
                                 integr*Geom%W_mu(ith)*Geom%W_mux(iph)

                          end do ! Q
                        end do ! K
                      end do ! frequencies
                    end do ! azimuthal nodes
                  end do ! polar nodes
                end do ! heights

                ! Complete
                do K=1,Krad
                  do iQ=1,K
                    JKQC(-iQ,K,:,Rz0:Rz1) = Flgsg%sg(iQ)* &
                                           conjg(JKQC(iQ,K,:,Rz0:Rz1))
                  end do
                end do

              end if ! Master

              ! If there are slaves
              if (MPID%mpi) then

                ! Control
                call control
                if (laborted) return

                csize = 15*nfreq*RnZ

                ! Alternative bcast
                if (MPID%altbcast) then

                  ! If not master, receive first
                  if (pid.ne.0) then

                    ! Receive JKQC
                    call MPI_RECV(JKQC(-2,0,1,Rz0), csize, &
                                  MPI_DOUBLE_COMPLEX, &
                                  MPID%recv, 6000000+pid, &
                                  MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, &
                                  ierr)

                  end if ! No Master

                  ! For each send
                  do istep=1,MPID%nsend

                    ! Send JKQC
                    call MPI_ISEND(JKQC(-2,0,1,Rz0), csize, &
                                   MPI_DOUBLE_COMPLEX, &
                                   MPID%lsend(istep), &
                                   6000000+MPID%lsend(istep), &
                                   MPI_COMM_RT, &
                                   MPID%requestA(istep,4), &
                                   ierr)

                  end do ! Sends

                ! Normal bcast
                else

                  call MPI_BCAST(JKQC(-2,0,1,Rz0), csize, &
                                 MPI_DOUBLE_COMPLEX, 0, &
                                 MPI_COMM_RT, ierr)

                end if ! Type of bcast
              end if ! MPI

            !
            ! Currently doing AD
            !
            else

              ! Master
              if(pid.eq.0)then

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph=1,Geom%nPh

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! For each Stokes parameter
                        do iS=0,3

                          ! If the input is not axial or it is the
                          ! first azimuth
                          if(.not.inaxial.or.iph.eq.1)then

                            ! Read the data
                            read (200,err=1100) da1

                            if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                            if (.not.axial.or.iph.eq.1) &
                              Stokes(iS,ifreq,iph,ith,iz) = da1

                          ! Input is axial
                          else

                            if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                            Stokes(iS,ifreq,iph,ith,iz) = &
                                            Stokes(iS,ifreq,1,ith,iz)

                          end if ! axial input

                        end do ! Stokes parameters
                      end do ! frequencies
                    end do ! azimuthal nodes
                  end do ! polar nodes
                end do ! heights

              end if ! Master

              ! If there are slaves
              if (MPID%mpi) then

                ! Control
                call control
                if (laborted) return

                ssize = 4*nfreq*Geom%nTh*Geom%nPh*Rnz

                ! Alternative bcast
                if (MPID%altbcast) then

                  ! If not master, receive first
                  if (pid.ne.0) then

                    ! Receive Stokes
                    call MPI_RECV(Stokes(0,1,1,1,Rz0), ssize, &
                                  MPI_DOUBLE_PRECISION, &
                                  MPID%recv, 6000000+pid, &
                                  MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, &
                                  ierr)

                  end if ! No Master

                  ! For each send
                  do istep=1,MPID%nsend

                    ! Send Stokes
                    call MPI_ISEND(Stokes(0,1,1,1,Rz0), ssize, &
                                   MPI_DOUBLE_PRECISION, &
                                   MPID%lsend(istep), &
                                   6000000+MPID%lsend(istep), &
                                   MPI_COMM_RT, &
                                   MPID%requestA(istep,4), ierr)

                  end do ! Sends

                ! Normal bcast
                else

                  call MPI_BCAST(Stokes(0,1,1,1,Rz0), ssize, &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                end if ! Type of bcast
              end if ! MPI

              !
              ! Calculate JKQC from Stokes
              !

              JKQC = dcmplx(.0D0,.0D0)

              ! For each height
              do iz=Rz0,Rz1

                ! For each frequency
                do ifreq=1,nfreq

                  ! For each K
                  do K=0,Krad

                    ! For each Q
                    do iQ=0,K

                      ! For each polar direction
                      do ith=1,Geom%nTh

                        ! For each azimuthal direction
                        do iph=1,Geom%nPh

                          ! Sum over Stokes parameters of Stokes*TKQ
                          integr = sum(Stokes(:,ifreq,iph,ith,iz)* &
                                     Geom%TS(:,iQ,K,iph,ith))

                          ! Add contribution to the JKQC integral
                          JKQC(iQ,K,ifreq,iz) = JKQC(iQ,K,ifreq,iz) &
                               + integr*Geom%W_mu(ith)*Geom%W_mux(iph)

                        end do ! azimuthal nodes
                      end do ! polar nodes
                    end do ! Q
                  end do ! K
                end do ! frequencies
              end do ! heights

              ! Complete
              do K=1,Krad
                do iQ=1,K
                  JKQC(-iQ,K,:,Rz0:Rz1) = Flgsg%sg(iQ)* &
                                          conjg(JKQC(iQ,K,:,Rz0:Rz1))
                end do
              end do

            end if ! Doing AV or AD

            !
            ! Q!=0 JKQC
            !

            ! If the input was axial or we are doing axial
            if(Geom%axial.or.inaxial)then

              ! Kill the Q!=0 components
              JKQC(-2:-1,1:2,:,:) = cZero
              JKQC(1:2,1:2,:,:) = cZero

            endif ! axiality
          end if ! Input AV or AD
        end if ! Can read Stokes


      !
      ! INTENSITY READ
      !
      else if (ilabel.eq.1) then

        !
        ! Allocations
        !

        ! If we are doing angle averaged, we only need one height,
        ! we allocate two to store the emergence in the quadrature
        if (KSTK) then
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1
          MPID%RRAM = MPID%RRAM + &
                      8d-6*dble(nfreq*GeomI%nPh*GeomI%nTh*Rnz)
        else
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1
          MPID%RRAM = MPID%RRAM + &
                      8d-6*dble(nfreq*GeomI%nPh*GeomI%nTh*2)
        end if
        Stokes0 = 0d0

        ! J00 for absorptivity
        allocate(J00(nxt,Rz0:Rz1))

        ! J00 for stimulated emission
        allocate(J00S(nxt,Rz0:Rz1))

        ! J00 frequency dependent
        allocate(J00C(nfreq,Rz0:Rz1))

        ! Allocated memory
        MPID%RRAM = MPID%RRAM + 8d-6*dble(Rnz*(2*nxt + nfreq))

        ! If the master
        if (pid.eq.0) then
          read (200,err=1100) ia1
          read (200,err=1100) ia2
          read (200,err=1100) ia3
          read (200,err=1100) ia4
          read (200,err=1100) ia5
          read (200,err=1100) inaxial_int
          read (200,err=1100) instm_int
          read (200,err=1100) inAV_int
        end if

        ! If doing MPI
        if (MPID%mpi) then

          allocate(ibuff(8))

          ! Control
          call control
          if (laborted) return

          ! If the master
          if (pid.eq.0) then
            ibuff(1) = ia1
            ibuff(2) = ia2
            ibuff(3) = ia3
            ibuff(4) = ia4
            ibuff(5) = ia5
            ibuff(6) = inaxial_int
            ibuff(7) = instm_int
            ibuff(8) = inAV_int
          end if

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive dimenions
              call MPI_RECV(ibuff(1), 8, MPI_INTEGER, &
                            MPID%recv, 2000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No Master

            ! For each send
            do istep=1,MPID%nsend

              ! Send dimensions
              call MPI_ISEND(ibuff(1), 8, MPI_INTEGER, &
                             MPID%lsend(istep), &
                             2000000+MPID%lsend(istep), &
                             MPI_COMM_RT, MPID%requestA(istep,7), &
                             ierr)

            end do ! sends

          ! Normal bcast
          else

            call MPI_BCAST(ibuff(1), 8, MPI_INTEGER, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast

          ! If a slave
          if (pid.ne.0) then
            ia1 = ibuff(1)
            ia2 = ibuff(2)
            ia3 = ibuff(3)
            ia4 = ibuff(4)
            ia5 = ibuff(5)
            inaxial_int = ibuff(6)
            instm_int = ibuff(7)
            inAV_int = ibuff(8)
          end if ! slave

          ! Alternative bcast
          if (MPID%altbcast) then

            ! For each slave to send
            do iproc=1,MPID%nsend

                ! Wait for everyone to receive the radiation data
                ! before continuing
                call MPI_WAIT(MPID%requestA(iproc,7), &
                              MPI_STATUS_IGNORE,ierr)

            end do

          end if ! Alternative bcast
        end if ! MPI

        ! Convert to logical
        if(inaxial_int.eq.1)then
          inaxial = .True.
        else
          inaxial = .False.
        end if

        ! Convert to logical
        if(instm_int.eq.1)then
          instm = .True.
        else
          instm = .False.
        end if

        ! Convert to logical
        if(inAV_int.eq.1)then
          inAV = .True.
        else
          inAV = .False.
        end if

        ! Flag to read stokes parameter
        read_stokes = .True.

        !
        ! Dimension checking
        !

        ! Check height nodes
        if (ia2.ne.nz) then
          umsg = 'Solution file with different number of '// &
                 'heights.'
          call aborted
          return
        end if

        ! Check number of atoms
        if (ia5.ne.nA) then
          umsg = 'Solution file with different number of '// &
                 'atoms.'
          call aborted
          return
        end if

        ! Check number of frequencies, this one does not produce an
        ! abortion
        if (ia1.ne.nfreq) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of frequencies in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes and J^K_Q(nu).'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Check polar nodes, this one does not produce an abortion
        if (ia3.ne.GeomI%nTh.and.read_stokes) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of polar directions in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes.'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Check azimuthal nodes, this one does not produce an abortion
        if (ia4.ne.GeomI%nPh.and.(ia4.ne.1.and.GeomI%nPh.ne.1).and. &
            read_stokes) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'and is not axial, ignoring Stokes.'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Check azimuthal nodes for AD redistribution, this one does
        ! not produce an abortion
        if (ia4.ne.GeomI%nPh2.and.(ia4.ne.1).and. &
            read_stokes) then
          if (pid.eq.0) then
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'for eps^(2) and is not axial, ignoring Stokes.'
            call verbose
          end if
          read_stokes = .False.
        end if

        ! Warning when reading non axial from AD for axial AD
        if (.not.inaxial.and.axiali.and.(.not.(AVI.and..not.dyn)) &
            .and..not.inAV.and.read_stokes.and.pid.eq.0) then
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read.'
          call verbose
        end if

        ! Warning when reading non axial from AD for axial AV
        if (.not.inaxial.and.axiali.and.(AVI.and..not.dyn).and. &
            .not.inAV.and.read_stokes.and.pid.eq.0) then
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read to compute J00C.'
          call verbose
        end if

        !
        ! Population and RhoKQ
        !

        ! For each atom
        do ia=1,nA

          ! If alternative bcast
          if (MPID%mpi.and.MPID%altbcast) then

            ! For each slave
            do iproc=1,MPID%nsend

              ! Wait for everyone to receive the radiation data before
              ! continuing and reset the buffers
              call MPI_WAIT(MPID%requestA(iproc,1), &
                            MPI_STATUS_IGNORE,ierr)
              call MPI_WAIT(MPID%requestA(iproc,6), &
                            MPI_STATUS_IGNORE,ierr)

            end do ! Processors

          end if ! Domain decomposition

          ! Only the master reads
          if (pid.eq.0) then

            ! Initialize pop
            Atom(ia)%crho = cZero

            ! For each height
            do iz=1,nZ

              ! Read population at this height
              read (200,err=1100) da1
              if (run_mode.ge.0) Atom(ia)%n(iz) = da1

            end do ! heights

            ! For each terms
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)

                ! Get component index
                iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                ! For each height
                do iz=1,nZ

                  ! Read real and imaginary parts
                  read (200,err=1100) da1
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                  Atom(ia)%crho(iR,iz) = dcmplx(da1,0d0)

                end do ! heights
              end do ! J
            end do ! terms

          end if ! Master

          ! If there are slaves
          if (MPID%mpi) then

            ! Control
            call control
            if (laborted) return

            nsize = nZ

            rsize = Atom(ia)%ndim*RnZ

            ! Alternative bcast
            if (MPID%altbcast) then

              ! If not master, receive first
              if (pid.ne.0) then

                ! Receive n
                call MPI_RECV(Atom(ia)%n(1), nsize, &
                              MPI_DOUBLE_PRECISION, &
                              MPID%recv, 2000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

                ! Receive rho00
                call MPI_RECV(Atom(ia)%crho(1,Rz0), rsize, &
                              MPI_DOUBLE_COMPLEX, &
                              MPID%recv, 3000000+pid, &
                              MPI_COMM_RT, MPI_STATUS_IGNORE, &
                              ierr)

              end if ! No Master

              ! For each send
              do istep=1,MPID%nsend

                ! Send n
                call MPI_ISEND(Atom(ia)%n(1), nsize, &
                               MPI_DOUBLE_PRECISION, &
                               MPID%lsend(istep), &
                               2000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,6), &
                               ierr)

                ! Send rho00
                call MPI_ISEND(Atom(ia)%crho(1,Rz0), rsize, &
                               MPI_DOUBLE_COMPLEX, &
                               MPID%lsend(istep), &
                               3000000+MPID%lsend(istep), &
                               MPI_COMM_RT, &
                               MPID%requestA(istep,1), &
                               ierr)

              end do ! sends

            ! Normal bcast
            else

              ! Share n
              call MPI_BCAST(Atom(ia)%n(1), nsize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

              ! Share rho00
              call MPI_BCAST(Atom(ia)%crho(1,Rz0), rsize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)

            end if ! bcast type
          end if ! MPI

        end do ! atoms

        !
        ! Flag the null rho(J,J')KQ
        !

        ! For each atom
        do ia=1,nA

          ! Initialize to non-null
          Atom(ia)%rhonull = .True.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)
              rJ = sqrt(2d0*rJ + 1d0)

              ! Get level index
              i = Atom(ia)%irho(it)%irho_ij(iJ)

              ! Get component index
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              do iz=Rz0,Rz1

                ! If rhoKQ/rho00 is lesser than double precision
                ! flag null
                if (abs(Atom(ia)%crho(iR,iz)).gt.0d0) &
                  Atom(ia)%rhonull(iR,iz) = .False.

                ! Define population
                Atom(ia)%popu(i,iz) = dble(Atom(ia)%crho(iR,iz))*rJ

              end do ! heights
            end do ! J
          end do ! terms
        end do ! atoms


        !
        ! J00
        !

        ! Only the master reads
        if (pid.eq.0) then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%nftran

              ! Apply atomic index
              jtran = itran + Atom(ia)%tfshift

              ! For each height
              do iz=1,nZ

                ! Read real and imaginary parts of JKQ
                read (200,err=1100) da1
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                J00(jtran,iz) = da1

              end do ! heights
            end do ! transitions
          end do ! atoms

        end if ! Master

        ! If there are slaves
        if (MPID%mpi) then

          ! Control
          call control
          if (laborted) return

          jsize = nxt*RnZ

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive J00
              call MPI_RECV(J00(1,Rz0), jsize, &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, 4000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No Master

            ! For each send
            do istep=1,MPID%nsend

              ! Send J00
              call MPI_ISEND(J00(1,Rz0), jsize, &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             4000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,2), ierr)

            end do ! Sends

          ! Normal bcast
          else

            call MPI_BCAST(J00(1,Rz0), jsize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, &
                           ierr)

          end if ! Type of bcast
        end if ! MPI

        ! If there is stimulated emission in the input
        if(instm)then

          ! Only the master reads
          if (pid.eq.0) then

            ! For each atom
            do ia=1,nA

              ! For each transition
              do itran=1,Atom(ia)%nftran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tfshift

                ! For each height
                do iz=1,nZ

                  ! Read real and imaginary parts of JKQ
                  read (200,err=1100) da1
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                  if (stm) J00S(jtran,iz) = da1

                end do ! heights
              end do ! transitions
            end do ! atoms

          end if ! Master

         !! If there are slaves
         !if (MPID%mpi) then

         !  ! Control
         !  call control

         !  jsize = nxt*RnZ

         !  ! Alternative bcast
         !  if (MPID%altbcast) then

         !    ! If not master, receive first
         !    if (pid.ne.0) then

         !      ! Receive J00S
         !      call MPI_RECV(J00S(1,Rz0), jsize, &
         !                    MPI_DOUBLE_PRECISION,  &
         !                    MPID%recv, 5000000+pid, &
         !                    MPI_COMM_RT, MPI_STATUS_IGNORE, &
         !                    ierr)

         !    end if ! No Master

         !    ! For each send
         !    do istep=1,MPID%nsend

         !      ! Send J00S
         !      call MPI_ISEND(J00S(1,Rz0), jsize, &
         !                     MPI_DOUBLE_PRECISION, &
         !                     MPID%lsend(istep), &
         !                     5000000+MPID%lsend(istep), &
         !                     MPI_COMM_RT, &
         !                     MPID%requestA(istep,3), ierr)

         !    end do ! Sends

         !  ! Normal bcast
         !  else

         !    call MPI_BCAST(J00S(1,Rz0), jsize, &
         !                   MPI_DOUBLE_PRECISION, 0, &
         !                   MPI_COMM_RT, ierr)

         !  end if ! Type of bcast
         !end if ! MPI

        ! No stimulated emission in the input
       !else

       !  if(stm)J00S = J00

        end if ! Stimulated emission in the input

        !
        ! J00P
        !

        ! Only the master reads
        if (pid.eq.0) then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%nphot

              ! Apply atomic shift
              jtran = itran + Atom(ia)%pshift

              ! For each height
              do iz=1,nZ

                ! Read real and imaginary parts of JKQ
                read (200,err=1100) da1
                read (200,err=1100) da2
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                J00P(jtran,1,iz) = da1
                J00P(jtran,2,iz) = da2

              end do ! heights
            end do ! transitions
          end do ! atoms

        end if ! Master

        ! If there are slaves
        if (MPID%mpi) then

          ! Control
          call control
          if (laborted) return

          psize = nxphot*2*RnZ

          ! Alternative bcast
          if (MPID%altbcast) then

            ! If not master, receive first
            if (pid.ne.0) then

              ! Receive J00 for b-f transitions
              call MPI_RECV(J00P(1,1,Rz0), psize, &
                            MPI_DOUBLE_PRECISION,  &
                            MPID%recv, 6000000+pid, &
                            MPI_COMM_RT, MPI_STATUS_IGNORE, &
                            ierr)

            end if ! No master

            ! For each send
            do istep=1,MPID%nsend

              ! Send J00P
              call MPI_ISEND(J00P(1,1,1), psize, &
                             MPI_DOUBLE_PRECISION, &
                             MPID%lsend(istep), &
                             6000000+MPID%lsend(istep), &
                             MPI_COMM_RT, &
                             MPID%requestA(istep,5), ierr)

            end do ! Sends

          ! Normal bcast
          else

            call MPI_BCAST(J00P(1,1,1), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Type of bcast
        end if ! MPI


        !
        ! Stokes or J00C
        !

        ! If we can read Stokes
        if (read_stokes) then

          ! If input is AV
          if (inAV) then

            ! Only the master reads
            if (pid.eq.0) then

              ! For each height
              do iz=1,nz

                ! For each frequency
                do ifreq=1,nfreq

                  ! Read real and imaginary parts of JKQS
                  read (200,err=1100) da1
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle
                  J00C(ifreq,iz) = da1

                end do ! frequency
              end do ! heights

            end if ! Master

            ! If there are slaves
            if (MPID%mpi) then

              ! Control
              call control
              if (laborted) return

              csize = nfreq*RnZ

              ! Alternative bcast
              if (MPID%altbcast) then

                ! If not master, receive first
                if (pid.ne.0) then

                  ! Receive J00C
                  call MPI_RECV(J00C(1,Rz0), csize, &
                                MPI_DOUBLE_PRECISION,  &
                                MPID%recv, 7000000+pid, &
                                MPI_COMM_RT, MPI_STATUS_IGNORE, &
                                ierr)

                end if ! No master

                ! For each send
                do istep=1,MPID%nsend

                  ! Send J00C
                  call MPI_ISEND(J00C(1,Rz0), csize, &
                                 MPI_DOUBLE_PRECISION, &
                                 MPID%lsend(istep), &
                                 7000000+MPID%lsend(istep), &
                                 MPI_COMM_RT, &
                                 MPID%requestA(istep,4), ierr)

                end do ! Sends

              ! Normal bcast
              else

                call MPI_BCAST(J00C(1,Rz0), csize, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

              end if ! Type of bcast
            end if ! MPI

            ! If we are doing angle dependent, the first iteration
            ! must be angle averaged if the input is angle average
            if (.not.AV) tbAD = .True.

          ! If input is AD
          else

            !
            ! Currently doing AV without velocities, no need to
            ! read Stokes
            !
            if (.not.KSTK.or..not.(dyn.or..not.AV)) then

              ! Master
              if(pid.eq.0)then

                J00C = 0d0

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,GeomI%nTh

                    ! For each azimuthal direction
                    do iph=1,max(GeomI%nPh,ia4)

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! If the input is not axial or is the first
                        ! azimuth
                        if(.not.inaxial.or.iph.eq.1)then

                          ! Read the data point
                          read (200,err=1100) da1

                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          if (.not.axiali.or.iph.eq.1) then

                            ! Store
                            if (KSTK) then
                              Stokes0(ifreq,iph,ith,iz) = da1
                            else
                              Stokes0(ifreq,iph,ith,1) = da1
                            end if

                          end if

                        ! If the input is axial
                        else

                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          ! Store
                          if (KSTK) then
                            Stokes0(ifreq,iph,ith,iz) = &
                                                Stokes0(ifreq,1,ith,1)
                          else
                            Stokes0(ifreq,iph,ith,1) = &
                                                Stokes0(ifreq,1,ith,1)
                          end if

                        end if ! axial input

                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! Add contribution to the JKQC integral
                        J00C(ifreq,iz) = J00C(ifreq,iz) + &
                                        Stokes0(ifreq,iph,ith,1)* &
                                        GeomI%W_mu(ith)* &
                                        GeomI%W_mux(iph)

                      end do ! frequencies
                    end do ! azimuthal nodes
                  end do ! polar nodes
                end do ! heights

              end if ! Master

              ! If there are slaves
              if (MPID%mpi) then

                ! Control
                call control
                if (laborted) return

                csize = nfreq*RnZ

                ! Alternative bcast
                if (MPID%altbcast) then

                  ! If not master, receive first
                  if (pid.ne.0) then

                    ! Receive J00C
                    call MPI_RECV(J00C(1,Rz0), csize, &
                                  MPI_DOUBLE_PRECISION,  &
                                  MPID%recv, 7000000+pid, &
                                  MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)

                  end if ! No master

                  ! For each send
                  do istep=1,MPID%nsend

                    ! Send J00C
                    call MPI_ISEND(J00C(1,Rz0), csize, &
                                   MPI_DOUBLE_PRECISION, &
                                   MPID%lsend(istep), &
                                   7000000+MPID%lsend(istep), &
                                   MPI_COMM_RT, &
                                   MPID%requestA(istep,4), ierr)

                  end do ! Sends

                ! Normal bcast
                else

                  call MPI_BCAST(J00C(1,Rz0), csize, &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                end if ! Type of bcast
              end if ! MPI

            !
            ! Currently doing AD
            !
            else

              ! Master
              if(pid.eq.0)then

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,GeomI%nTh

                    ! For each azimuthal direction
                    do iph=1,max(GeomI%nPh,ia4)

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! If the input is not axial or it is the first
                        ! azimuth
                        if(.not.inaxial.or.iph.eq.1)then

                          ! Read the data
                          read (200,err=1100) da1

                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          if (.not.axiali.or.iph.eq.1) &
                            Stokes0(ifreq,iph,ith,iz) = da1

                        ! Input is axial
                        else

                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          Stokes0(ifreq,iph,ith,iz) = &
                                         Stokes0(ifreq,1,ith,iz)

                        end if ! axial input

                      end do ! frequencies
                    end do ! azimuthal nodes
                  end do ! polar nodes
                end do ! heights

              end if ! Master

              ! If there are slaves
              if (MPID%mpi) then

                ! Control
                call control
                if (laborted) return

                ssize = nfreq*GeomI%nTh*GeomI%nPh*Rnz

                ! Alternative bcast
                if (MPID%altbcast) then

                  ! If not master, receive first
                  if (pid.ne.0) then

                    call MPI_RECV(Stokes0(1,1,1,Rz0), ssize, &
                                  MPI_DOUBLE_PRECISION,  &
                                  MPID%recv, 7000000+pid, &
                                  MPI_COMM_RT, &
                                  MPI_STATUS_IGNORE, ierr)

                  end if ! No master

                  ! For each send
                  do istep=1,MPID%nsend

                    call MPI_ISEND(Stokes0(1,1,1,Rz0), ssize, &
                                   MPI_DOUBLE_PRECISION, &
                                   MPID%lsend(istep), &
                                   7000000+MPID%lsend(istep), &
                                   MPI_COMM_RT, &
                                   MPID%requestA(istep,4), ierr)

                  end do ! Sends

                ! Normal bcast
                else

                  call MPI_BCAST(Stokes0(1,1,1,Rz0), ssize, &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                end if ! Type of bcast
              end if ! MPI

              !
              ! Calculate JKQC from Stokes
              !

              J00C = 0d0

              ! For each height
              do iz=Rz0,Rz1

                ! For each frequency
                do ifreq=1,nfreq

                  ! For each polar direction
                  do ith=1,GeomI%nTh

                    ! For each azimuthal direction
                    do iph=1,GeomI%nPh

                      ! Add contribution to the JKQC integral
                      J00C(ifreq,iz) = J00C(ifreq,iz) + &
                                       Stokes0(ifreq,iph,ith,iz)* &
                                       GeomI%W_mu(ith)* &
                                       GeomI%W_mux(iph)

                    end do ! azimuthal nodes
                  end do ! polar nodes
                end do ! frequencies
              end do ! heights

            end if ! Doing AV or AD
          end if ! Input AV or AD
        end if ! Can read Stokes
      end if ! Type of input

      ! Close unit
      if (pid.eq.0) close (200)

      ! If alternative bcast
      if (MPID%mpi.and.MPID%altbcast) then

        ! For each slave
        do iproc=1,MPID%nsend

          ! Wait for everyone to receive the radiation data before
          ! continuing and reset the buffers
          call MPI_WAIT(MPID%requestA(iproc,1), &
                        MPI_STATUS_IGNORE,ierr)
          call MPI_WAIT(MPID%requestA(iproc,2), &
                        MPI_STATUS_IGNORE,ierr)
          call MPI_WAIT(MPID%requestA(iproc,3), &
                        MPI_STATUS_IGNORE,ierr)
          call MPI_WAIT(MPID%requestA(iproc,4), &
                        MPI_STATUS_IGNORE,ierr)
          if (ilabel.eq.1) &
          call MPI_WAIT(MPID%requestA(iproc,5), &
                        MPI_STATUS_IGNORE,ierr)
          call MPI_WAIT(MPID%requestA(iproc,6), &
                        MPI_STATUS_IGNORE,ierr)

        end do ! Processors

      end if ! Domain decomposition

      ! Control
      call control

      return

1000  umsg = 'Error opening solution file '//trim(filename)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error reading solution file '//trim(filename)
      close(100)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine readsol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the arrays from a previous solution in RAM.\n
      !!    SolF(Solution_F_class): Class with solution data\n
      !!     GeomI(Geometry_class): Structure with geometry data for
      !!                            the intensity problem\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!        intensity(logical): Need to fetch intensity
      !!                            solution
      subroutine getsol(SolF,GeomI,Geom,MPID,Flgsg,Bfield,Atom, &
                        Stokes,JKQ,JKQS,JKQC, &
                        Stokes0,J00,J00S,J00C,J00P, &
                        intensity)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Geometry_class), intent(in):: Geom,GeomI
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(MPI_class), intent(inout):: MPID
      type(Solution_F_class), intent(in):: SolF
      logical, intent(in):: intensity
      double precision,dimension(:,:,:,:), allocatable:: Stokes0
      double precision,dimension(:,:,:,:,:), allocatable:: Stokes
      double precision,dimension(:,:), allocatable:: J00
      double precision,dimension(:,:), allocatable:: J00S
      double precision,dimension(:,:), allocatable:: J00C
      double precision,dimension(:,:,:), allocatable:: J00P
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      integer:: ia,iz,i,itran,it,iJ,iJ1,K,iQ,iR,iR0,iR1,psize

      double precision:: rJ,rJ1,rho0

      complex(kind=8),dimension(-nkx:nkx,nz):: rhoKQaux


      ! Routine name
      urou = 'getsol'

      !
      ! If there is a polarization solution
      !
      if (intensity) then

        !
        ! Allocations
        !

        ! If we are doing angle averaged, we only need one height,
        ! we allocate two to store the emergence in the quadrature
        if (KSTK) then
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1
          MPID%RRAM = 8d-6*dble(nfreq*GeomI%nph*GeomI%nTh*Rnz)
        else
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1
          MPID%RRAM = 8d-6*dble(nfreq*GeomI%nph*GeomI%nTh*2)
        end if

        ! J00 for absorptivity
        allocate(J00(nxt,Rz0:Rz1))

        ! J00 for stimulated emission
        allocate(J00S(nxt,Rz0:Rz1))

        ! J00 frequency dependent
        allocate(J00C(nfreq,Rz0:Rz1))

        ! J00 for photoionizations
        allocate(J00P(nxphot,2,Rz0:Rz1))

        ! Compute allocated memory
        MPID%RRAM = MPID%RRAM + 8d-6*dble(Rnz*(nxphot*2 + &
                                          2*nxt + nfreq))

        ! Master
        if (pid.eq.0) then

          ! If keeping Stokes
          if (KSTK) then

            ! Save
            Stokes0 = dble(SolF%i_StkI_b(:,:,:,Rz0:Rz1))

            ! And share
            if (MPID%mpi) then
            psize = nfreq*GeomI%nph*GeomI%nth*Rnz
            call MPI_BCAST(Stokes0(1,1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)
            end if

          ! Otherwise
          else

            Stokes0 = 0d0

          end if

          ! J00 bar
          J00 = SolF%i_J00_b(:,Rz0:Rz1)

          ! Share
          if (MPID%mpi) then
          psize = nxt*RnZ
          call MPI_BCAST(J00(1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)
          end if

          ! J00 photo
          J00P = SolF%i_J00P_b(:,:,Rz0:Rz1)

          ! Share
          if (MPID%mpi) then
          psize = nxphot*2*Rnz
          call MPI_BCAST(J00P(1,1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)
          end if

          ! J00 freq. dependent
          J00C = SolF%i_J00C_b(:,Rz0:Rz1)

          ! Share
          if (MPID%mpi) then
          psize = nfreq*Rnz
          call MPI_BCAST(J00C(1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)
          end if

        ! Slaves
        else

          ! If keeping Stokes
          if (KSTK) then

            ! Get
            psize = nfreq*GeomI%nph*GeomI%nth*Rnz
            call MPI_BCAST(Stokes0(1,1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Otherwise
          else

            Stokes0 = 0d0

          end if

          ! Get
          psize = nxt*RnZ
          call MPI_BCAST(J00(1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Get
          psize = nxphot*2*Rnz
          call MPI_BCAST(J00P(1,1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Get
          psize = nfreq*Rnz
          call MPI_BCAST(J00C(1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

        end if

        ! Copy in J00S
        J00S = J00

        ! For each atom
        do ia=1,nA

          ! Initialize
          Atom(ia)%crho = cZero

          ! Master
          if (pid.eq.0) then

            ! Copy popu
            Atom(ia)%popu = SolF%i_rhoes_b(ia)%rho

          end if ! Master/slave

          !
          ! And share populations
          !

          ! Share
          if (MPID%mpi) then
          psize = Rnz*Atom(ia)%nlevel
          call MPI_BCAST(Atom(ia)%popu(1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)
          end if

          ! For each height
          do iz=Rz0,Rz1

            ! For each term
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)!,1,-1

                ! Level and KQ
                i = Atom(ia)%irho(it)%irho_ij(iJ)
                iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)
                rJ = sqrt(2d0*rJ + 1d0)

                ! Set up rho
                Atom(ia)%crho(iR,iz) = &
                                     dcmplx(Atom(ia)%popu(i,iz), 0d0)

                ! Define population
                Atom(ia)%popu(i,iz) = Atom(ia)%popu(i,iz)*rJ

              end do
            end do
          end do


          !
          ! Flag the null rho(J,J')KQ
          !

          ! Initialize to non-null
          Atom(ia)%rhonull = .True.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)
              rJ = sqrt(2d0*rJ + 1d0)

              ! Get level index
              i = Atom(ia)%irho(it)%irho_ij(iJ)

              ! Get component index
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              do iz=Rz0,Rz1

                ! If rhoKQ/rho00 is lesser than double precision
                ! flag null
                if (abs(Atom(ia)%crho(iR,iz)).gt.0d0) &
                  Atom(ia)%rhonull(iR,iz) = .False.

              end do ! heights
            end do ! J
          end do ! terms
        end do ! Atoms

      ! Polarization
      else

        ! If we are doing angle averaged, we only need one height,
        ! we allocate two to store the emergence in the quadrature
        if (KSTK) then
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1
          MPID%RRAM = 8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*Rnz)
        else
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1
          MPID%RRAM = 8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*2)
        end if

        ! JKQ for absorptivity
        allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))

        ! JKQ for stimulated emission
        allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))

        ! JKQ frequency dependent
        allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))

        ! J00 for photoionizations
        allocate(J00P(nxphot,2,Rz0:Rz1))

        ! Compute allocated memory
        MPID%RRAM = MPID%RRAM + 8d-6*dble(Rnz*(nxphot*2 + &
                                          2*5*3*(2*nxtran + nfreq)))

        ! Master
        if (pid.eq.0) then

          ! If keeping Stokes
          if (KSTK) then

            ! Save
            Stokes = SolF%i_Stk_b(:,:,:,:,Rz0:Rz1)

            ! And share
            if (MPID%mpi) then
            psize = 4*nfreq*Geom%nph*Geom%nth*Rnz
            call MPI_BCAST(Stokes(0,1,1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)
            end if

          ! Otherwise
          else

            Stokes = 0d0

          end if

          ! JKQ bar
          JKQ = SolF%i_JKQ_b(:,:,:,Rz0:Rz1)

          !
          ! Rotate JKQ
          !

          ! For each height
          do iz=Rz0,Rz1

            ! No field, skip
            if (Bfield%Bstrength(iz).le.TINYB) cycle

            ! For each transition
            do itran=1,nxtran

              ! Rotate
              if (Bfield%Bstrength(iz).gt.TINYB) &
              call fieldB(JKQ(:,:,itran,iz),1,Flgsg, &
                          Bfield%Btheta(iz),Bfield%Bphi(iz),1)

            end do ! Transitions
          end do ! heights

          ! Share
          if (MPID%mpi) then
          psize = 15*nxtran*RnZ
          call MPI_BCAST(JKQ(-2,0,1,Rz0), psize, &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)
          end if

          ! JKQS bar
          if (stm) then

            JKQS = SolF%i_JKQS_b(:,:,:,Rz0:Rz1)

            !
            ! Rotate JKQ
            !

            ! For each height
            do iz=Rz0,Rz1

              ! No field, skip
              if (Bfield%Bstrength(iz).le.TINYB) cycle

              ! For each transition
              do itran=1,nxtran

                ! Rotate
                call fieldB(JKQS(:,:,itran,iz),1,Flgsg, &
                            Bfield%Btheta(iz),Bfield%Bphi(iz),1)

              end do ! Transitions
            end do ! heights

            ! Share
            if (MPID%mpi) then
            psize = 15*nxtran*RnZ
            call MPI_BCAST(JKQS(-2,0,1,Rz0), psize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)
            end if
          end if ! Stimulated emission

          ! JKQ freq. dependent
          JKQC = SolF%i_JKQC_b(:,:,:,Rz0:Rz1)

          ! Share
          if (MPID%mpi) then
          psize = 15*nfreq*Rnz
          call MPI_BCAST(JKQC(-2,0,1,Rz0), psize, &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)
          end if

        ! Slaves
        else

          ! If keeping Stokes
          if (KSTK) then

            ! Get
            psize = 4*nfreq*Geom%nph*Geom%nth*Rnz
            call MPI_BCAST(Stokes(0,1,1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Otherwise
          else

            Stokes = 0d0

          end if

          ! Get
          psize = 15*nxtran*RnZ
          call MPI_BCAST(JKQ(-2,0,1,Rz0), psize, &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

          ! If stimulatted
          if (stm) then

            ! Get
            psize = 15*nxtran*RnZ
            call MPI_BCAST(JKQS(-2,0,1,Rz0), psize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

          end if ! Stimulated emission

          ! Get
          psize = 15*nfreq*Rnz
          call MPI_BCAST(JKQC(-2,0,1,Rz0), psize, &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Master/slave

        ! For each atom
        do ia=1,nA

          ! Initialize
          Atom(ia)%crho = cZero

          ! Master
          if (pid.eq.0) then

            ! For each terms
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)

                ! For each level
                do iJ1=1,Atom(ia)%nJ(it)

                  ! Get J'
                  rJ1 = Atom(ia)%rJval(iJ1,it)

                  ! For each K
                  do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                    ! For each Q
                    do iQ=-K,K

                      ! Get rho(J,J')KQ index
                      if (K.le.Atom(ia)%Kcut(it)) &
                        iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! For each height
                      do iz=1,nZ

                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        if (K.le.Atom(ia)%Kcut(it)) then

                          Atom(ia)%crho(iR,iz) = &
                                        SolF%i_rhoes_b(ia)%crho(iR,iz)

                          ! Auxiliar variable for rotation
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end if

                      end do ! heights
                    end do ! Q

                    !
                    ! Rotate rhoKQ
                    !

                    ! Only if above the limit
                    if (K.le.Atom(ia)%Kcut(it)) then

                      ! For each height in the CPU domain
                      do iz=Rz0,Rz1

                        ! If there is magnetic field
                        if (Bfield%Bstrength(iz).le.TINYB) cycle

                        ! Rotate the rhoKQ in the auxiliar variable
                        call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                  Bfield%Btheta(iz), &
                                  Bfield%Bphi(iz),1)

                        ! Store the rotated result in the rhoKQ
                        ! array
                        do iQ=-K,K
                          iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)
                          Atom(ia)%crho(iR,iz) = rhoKQaux(iQ,iz)
                        end do ! Q
                      end do ! heights

                    end if ! K < Kcut

                  end do ! K
                end do ! J'
              end do ! J
            end do ! terms

          end if ! Master/slave

          ! Share
          if (MPID%mpi) then
            psize = Atom(ia)%ndim*Rnz
            call MPI_BCAST(Atom(ia)%crho(1,Rz0), psize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)
          end if


          !
          ! Flag the null rho(J,J')KQ
          !

          ! Initialize to non-null
          Atom(ia)%rhonull = .False.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)

              ! Get the rho00 indexes
              iR0 = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! For each level
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J'
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! Get the rho00 indexes
                iR1 = Atom(ia)%irho(it)%Jrho(iJ1,iJ1)%kq(0,0)

                ! For each K
                do K=nint(abs(rJ-rJ1)), &
                     min(nint(rJ+rJ1),Atom(ia)%Kcut(it))

                  ! For each Q
                  do iQ=-K,K

                    ! Get rho(J,J')KQ index
                    iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                    ! For each height in this domain
                    do iz=Rz0,Rz1

                      ! Get the inverse of rho00
                      rho0 = 1d0/sqrt(abs(Atom(ia)%crho(iR0,iz))* &
                                      abs(Atom(ia)%crho(iR1,iz)))

                      ! If rhoKQ/rho00 is lesser than double precision
                      ! flag null
                      if (abs(Atom(ia)%crho(iR,iz)*rho0).lt.TINYR) &
                        Atom(ia)%rhonull(iR,iz) = .True.

                    end do ! heights

                  end do ! Q
                end do ! K
              end do ! J'
            end do ! J
          end do ! terms
        end do ! atoms

      end if ! Intensity/polarization initialization

      end subroutine getsol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Store the solution in the SolF structure.\n
      !!    SolF(Solution_F_class): Class with solution data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!  Stokes0(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!        intensity(logical): Need to store just the intensity
      !!                            solution
      subroutine setsol(SolF,Flgsg,Bfield,Atom, &
                        Stokes,JKQ,JKQS,JKQC, &
                        Stokes0,J00,J00S,J00C,J00P, &
                        intensity)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: intensity
      double precision,dimension(:,:,:,:), allocatable:: Stokes0
      double precision,dimension(:,:,:,:,:), allocatable:: Stokes
      double precision,dimension(:,:), allocatable:: J00
      double precision,dimension(:,:), allocatable:: J00S
      double precision,dimension(:,:), allocatable:: J00C
      double precision,dimension(:,:,:), allocatable:: J00P
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8),dimension(:,:,:,:), allocatable:: JKQC

      ! Local

      integer:: ia,it,iJ,iJ1,iR,ii,itran,iz,K,iQ

      double precision:: rJ,rJ1

      complex(kind=8),dimension(-nkx:nkx,nz):: rhoKQaux


      ! Only master
      if (pid.gt.0) return


      !
      ! Intensity
      !
      if (intensity) then

        !
        ! Populations
        !

        ! For each atom
        do ia=1,nA

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level J within the term
            do iJ=1,Atom(ia)%nJ(it)!,1,-1

              ! Get the level and rho00 index
              ii = Atom(ia)%irho(it)%irho_ij(iJ)
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! Save rho00
              SolF%i_rhoes(ia)%rho(ii,:) = dble(Atom(ia)%crho(iR,:))

            end do ! Sublevel
          end do ! Term
        end do ! Atom

        !
        ! Mean radiation field tensors
        !
        SolF%i_J00(:,1:Rz0-1) = 0d0
        SolF%i_J00(:,Rz1+1:nz) = 0d0
        SolF%i_J00(:,Rz0:Rz1) = J00(:,Rz0:Rz1)

        !
        ! Radiation field tensors
        !
        SolF%i_J00C(:,1:Rz0-1) = 0d0
        SolF%i_J00C(:,Rz1+1:nz) = 0d0
        SolF%i_J00C(:,Rz0:Rz1) = J00C(:,Rz0:Rz1)

        !
        ! Photoionization
        !
        SolF%i_J00P(:,:,1:Rz0-1) = 0d0
        SolF%i_J00P(:,:,Rz1+1:nz) = 0d0
        SolF%i_J00P(:,:,Rz0:Rz1) = J00P(:,:,Rz0:Rz1)

        !
        ! Stokes
        !
        if (KSTK) then
          SolF%i_StkI(:,:,:,1:Rz0-1) = 0d0
          SolF%i_StkI(:,:,:,Rz1+1:nz) = 0d0
          SolF%i_StkI(:,:,:,Rz0:Rz1) = Stokes0(:,:,:,Rz0:Rz1)
        end if

      !
      ! Polarization
      !
      else

        !
        ! Density matrix
        !

        ! For each atom
        do ia=1,nA

          ! Initialize
          SolF%i_rhoes(ia)%crho = cZero

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level J within the term
            do iJ=1,Atom(ia)%nJ(it)!,1,-1

              ! Get J value
              rJ = Atom(ia)%rJval(iJ,it)

              ! For each J value within the term
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J value
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! For each K
                do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                  !
                  ! Rotate rhoKQ nto the vertical reference frame
                  if (K.le.Atom(ia)%Kcut(it)) then

                    ! For each height
                    do iz=1,nz

                      ! If out of bounds
                      if (iz.lt.Rz0.or.iz.gt.Rz1) then

                        rhoKQaux(-K:K,iz) = cZero

                      ! In bounds
                      else

                        ! For each Q
                        do iQ=-K,K

                          ! Get the rhoKQ index
                          iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                          ! Store the corresponding rhoKQ into an
                          ! auxiliar variable
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end do ! Q

                        ! If there is non-zero magnetic field, rotate
                        if (Bfield%Bstrength(iz).gt.TINYB) &
                          call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                   -Bfield%Btheta(iz), &
                                   -Bfield%Bphi(iz),-1)

                      end if ! Height bounds

                    end do ! heights

                    ! For each Q
                    do iQ=-K,K

                      ! Get the index
                      iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! Save
                      SolF%i_rhoes(ia)%crho(iR,:) = &
                                                   Atom(ia)%crho(iR,:)

                    end do ! Q

                  end if ! K<=Kcut

                end do ! K
              end do ! J'
            end do ! J
          end do ! Terms
        end do ! Atoms

        !
        ! Mean radiation field tensors
        !

        ! Heights
        do iz=1,nz

          ! Out of limits
          if (iz.lt.Rz0.or.iz.gt.Rz1) then
            SolF%i_JKQ(:,:,:,iz) = cZero
            if (stm) SolF%i_JKQS(:,:,:,iz) = cZero
            cycle
          end if

          ! Store
          SolF%i_JKQ(:,:,:,iz) = JKQ(:,:,:,iz)
          if (stm) SolF%i_JKQS(:,:,:,iz) = JKQS(:,:,:,iz)

          ! If there is a magnetic field
          if (Bfield%Bstrength(iz).gt.TINYB) then

            ! For each transition
            do itran=1,nxtran

              ! Rotate
              call fieldB(SolF%i_JKQ(:,:,itran,iz),1,Flgsg, &
                          -Bfield%Btheta(iz),-Bfield%Bphi(iz),-1)
              if (stm) &
              call fieldB(SolF%i_JKQS(:,:,itran,iz),1,Flgsg, &
                          -Bfield%Btheta(iz),-Bfield%Bphi(iz),-1)

            end do

          end if ! Magnetic field

        end do ! Heights

        !
        ! Radiation field tensors
        !
        SolF%i_JKQC(:,:,:,1:Rz0-1) = cZero
        SolF%i_JKQC(:,:,:,Rz1+1:nz) = cZero
        SolF%i_JKQC(:,:,:,Rz0:Rz1) = JKQC(:,:,:,Rz0:Rz1)

        !
        ! Stokes
        !
        if (KSTK) then
          SolF%i_Stk(:,:,:,:,1:Rz0-1) = 0d0
          SolF%i_Stk(:,:,:,:,Rz1+1:nz) = 0d0
          SolF%i_Stk(:,:,:,:,Rz0:Rz1) = Stokes(:,:,:,:,Rz0:Rz1)
        end if

      end if ! Polarization/intensity

      end subroutine setsol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the computed solution that can be read
      !! by the code as initialization, and separate files for the
      !! radiation field and density tensors.\n
      !!        Input(Input_class): Structure with settings data\n
      !!        suff(character(:)): Suffix for the file names of the
      !!                            radiation field and density matrix
      !!                            files\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!              z(dfloat(:)): Heights array\n
      !! Stokes(dfloat(:,:,:,:,:)): Stokes parameters\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!   JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      subroutine writesol(Input,suff,omega,Geom, &
                          Flgsg,Bfield,Atom,z,Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(in):: Geom
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      character(len=4), intent(in):: suff
      double precision, dimension(:), intent(in):: omega,z
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                     giz0:giz1), intent(in):: Stokes
      complex(kind=8), &
             dimension(-2:2,0:2,nxtran,Rz0:Rz1), intent(in):: JKQ
      complex(kind=8), &
             dimension(-2:2,0:2,nxtran,Rz0:Rz1), intent(in):: JKQS
      complex(kind=8), &
             dimension(-2:2,0:2,nfreq,Rz0:Rz1), intent(in):: JKQC

      ! Local

      integer, parameter:: izero = 0
      integer, parameter:: ione = 1
      double precision, parameter:: dzero = 0d0

      character(len=8):: scoord
      character(len=500):: filename

      logical:: saveJKQnu, saverKQ, saveJKQ, saveP
      logical::  saveD, saveS, saveSol
      logical:: laux

      integer:: ierr,ii,iab,iran
      integer:: ios,ia,iz,ifreq,i,ith,iph,iS,itran,jtran
      integer:: it,iJ,iJ1,K,iQ,iR
      integer:: axial_int,stm_int,AV_int

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: rJ,rJ1,loffset

      complex(kind=8),dimension(-2:2,0:2,nz):: JKQaux
      complex(kind=8),dimension(-nkx:nkx,nz):: rhoKQaux

      !
      ! Slaves just wait
      !
      if (pid.gt.0) then
        call control
        return
      end if

      !
      ! Translate
      !
      filename = Input%folder
      saveSol = Input%keep_sol
      saveP = Input%keep_pop.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveD = Input%keep_dep.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saverKQ = Input%keep_rhoKQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveJKQ = Input%keep_JKQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveS = Input%keep_stokesQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveJKQnu = Input%keep_jkqnu.and.(suff.eq.'NONE'.or. &
                                        run_mode.eq.0)

      ! If not writing anything, come back
      if (.not.(saveSol.or.saveP.or.saveD.or.saverKQ.or.saveJKQ.or. &
                saveS.or.saveJKQnu)) then
        call control
        return
      end if

      ! Routine name
      urou = 'writesol'

      !
      ! Open files
      !

      !
      ! To write the solution
      !
      if (saveSol) then

        ! If 1D and inversion
        if (run_mode.le.0) then

          ! Open solution file
          open (200,file=trim(filename)//'/Solution', &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='write', &
                form='unformatted')

        ! If 1.5D
        else if (run_mode.eq.1) then

          ! Get LOS index
          write(scoord,'(I0.8)') icoords(3)

          ! Open solution file
          open (200,file=trim(filename)// &
                '/Solution-folder/Solution-'// &
                scoord, status='unknown', iostat=ios, err=1000, &
                access='stream', action='write', form='unformatted')

        end if ! 1D vs 1.5D

      end if

      ! If saving JKQ
      if (saveJKQ) then

        ! To write the final JKQ
        if (suff.eq.'NONE') then
          ! If 1D or inversion
          if (run_mode.le.0) then
            open (300,file=trim(filename)//'/Jout', &
                  status='unknown', iostat=ios, err=1002, &
                  access='stream', action='write', &
                  form='unformatted')
          ! 1.5D
          else if (run_mode.eq.1) then
            open (300,file=trim(filename)// &
                  '/Solution-folder/Jout-'//scoord,&
                  status='unknown', iostat=ios, err=1002, &
                  access='stream', action='write', &
                  form='unformatted')
          end if
        ! With suffix
        else
            open (300,file=trim(filename)//'/Jout_'//suff, &
                  status='unknown', iostat=ios, err=1002, &
                  access='stream', action='write', &
                  form='unformatted')
        end if

      end if

      ! If saving rhoKQ
      if (saverKQ) then

        ! To write the final rhoKQ
        if (suff.eq.'NONE') then
          ! If 1D or inversion
          if (run_mode.le.0) then
            open (400,file=trim(filename)//'/Rhoout', &
                  status='unknown', iostat=ios, err=1003, &
                  access='stream', action='write', &
                  form='unformatted')
          ! If 1.5D
          else if (run_mode.eq.1) then
            open (400,file=trim(filename)// &
                  '/Solution-folder/Rhoout-'// &
                  scoord, status='unknown', iostat=ios, err=1003, &
                  access='stream', action='write', &
                  form='unformatted')
          end if
        else
          open (400,file=trim(filename)//'/Rhoout_'//suff, &
                status='unknown', iostat=ios, err=1003, &
                access='stream', action='write', &
                form='unformatted')
        end if
      end if

      !
      ! Convert logicals to integers
      !

      ! Axial symmetry
      if(Geom%axial)then
          axial_int = 1
      else
          axial_int = 0
      end if

      ! Stimulated emission
      if(stm)then
          stm_int = 1
      else
          stm_int = 0
      end if

      ! Angle averaged redistribution function
      if (KSTK) then
          AV_int = 0
      else
          AV_int = 1
      end if

      !
      ! Write headers with dimensions and flags
      !

      ! If saving solution
      if (saveSol) then

        ! Solution file
        write (200,err=1100) 'sp'
        write (200,err=1100) nfreq
        write (200,err=1100) nZ
        write (200,err=1100) Geom%nTh
        write (200,err=1100) Geom%nPh
        write (200,err=1100) nA
        write (200,err=1100) axial_int
        write (200,err=1100) stm_int
        write (200,err=1100) AV_int

      end if

      ! If saving JKQ
      if (saveJKQ) then

        ! JKQ file
        write(300,err=1102) 'bj'
        write(300,err=1102) stm
        write(300,err=1102) nZ
        write(300,err=1102) nA
        write(300,err=1102) nxtran
        write(300,err=1102) z

      end if

      ! If saving rhoKQ
      if (saverKQ) then

        ! rhoKQ file
        write(400,err=1103) 'br'
        write(400,err=1103) nZ
        write(400,err=1103) nA
        write(400,err=1103) z

      end if

      ! If 1.5D, prepare pop y dep buffers
      if (run_mode.eq.1) then
        if (saveP.or.saveD) &
          allocate(buffer(maxval(Input%lim_pop%nbuff)/4))
      end if


      !
      ! Write the data
      !

      ! Only if saving anything
      if (saveSol.or.saveP.or.saveD.or.saverKQ) then

      !
      ! Population and rhoKQ

      ! For each atom
      do ia=1,nA

        !
        ! If saving population or departure coeff.
        !
        if (saveP.or.saveD) then

          !
          ! If 1D
          !
          if (run_mode.eq.0) then

            ! If writing populations
            if (saveP) then

              ! To write the populations
              open (500,file=trim(filename)//'/'// &
                    trim(Atom(ia)%file_label)//'.pop', &
                    status='unknown',iostat=ios, err=1004, &
                    access='stream', action='write', &
                    form='unformatted')

              write(500,err=1104) 'bp'
              write(500,err=1104) nZ
              write(500,err=1104) Atom(ia)%nlevel

            end if

            ! If writing departure c.
            if (saveD) then

              ! To write departure c.
              open (600,file=trim(filename)//'/'// &
                    trim(Atom(ia)%file_label)//'.dep', &
                    status='unknown',iostat=ios, err=1005, &
                    access='stream', action='write', &
                    form='unformatted')

              write(600,err=1105) 'bb'
              write(600,err=1105) nZ
              write(600,err=1105) Atom(ia)%nlevel

            end if

            ! For each height
            do iz=1,nZ

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! For each level
                do it=1,Atom(ia)%nlevel

                  ! Write populations
                  if (saveP) &
                    write (500,err=1104) Atom(ia)%popu(it,iz)
                  ! Write departure
                  if (saveD) &
                    write (600,err=1104) Atom(ia)%popu(it,iz)/ &
                                         Atom(ia)%populte(it,iz)

                end do

              ! In bounds
              else

                ! For each term
                do it=1,Atom(ia)%nMulti

                  ! For each level J within the term
                  do iJ=1,Atom(ia)%nJ(it)!,1,-1

                    ! Get J value and level index
                    rJ = Atom(ia)%rJval(iJ,it)
                    i = Atom(ia)%irho(it)%irho_ij(iJ)
                    iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                    ! Write populations
                    if (saveP) &
                      write (500,err=1104) Atom(ia)%n(iz)* &
                                       sqrt(2d0*rJ+1d0)* &
                                       dble(Atom(ia)%crho(iR,iz))

                    ! Write departure coeff
                    if (saveD) &
                      write (600,err=1105) Atom(ia)%n(iz)* &
                                       sqrt(2d0*rJ+1d0)* &
                                       dble(Atom(ia)%crho(iR,iz))/ &
                                       Atom(ia)%populte(i,iz)
                  end do ! Levels
                end do ! Terms

              end if ! Height bounds

            end do ! Heights

            if (saveP) close(500)
            if (saveD) close(600)

          !
          ! If 1.5D
          !
          else if (run_mode.eq.1) then

            ! Populations
            if (saveP.and.Input%lim_pop%nbuff(ia).gt.0) then

              ! Open file to write the populations
              call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.pop', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
              if (ierr.ne.0) goto 1004

              !
              ! Column offset
              !

              ! Get offset
              loffset = dble(icoords(3)-1)* &
                        dble(Input%lim_pop%nbuff(ia)) + &
                        dble(Input%lim_pop%head_size)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
                if (ierr.ne.0) goto 1014
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1014

              ! Initialize buffer
              ii = 0

              ! If specified
              if (Input%lim_pop%nran.gt.0) then

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Advance buffer
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(i,iz))

                    end do ! Entries

                  ! In bounds
                  else

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Get necessary data
                      it = Atom(ia)%term(i)
                      iJ = Atom(ia)%sublevel(i)
                      rJ = Atom(ia)%rJval(iJ,it)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Advance buffer
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz)))

                    end do ! Ranges to print

                  end if ! Height bounds

                end do ! Heights

              ! Everything
              else

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each level
                    do it=1,Atom(ia)%nlevel

                      ! Write populations
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(it,iz))
                    end do ! Levels

                  ! In bounds
                  else

                    ! For each term
                    do it=1,Atom(ia)%nMulti

                      ! For each level J within the term
                      do iJ=1,Atom(ia)%nJ(it)!,1,-1

                        ! Get J value and level index
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Write populations
                        ii = ii + 1
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                          sqrt(2d0*rJ+1d0)* &
                                          dble(Atom(ia)%crho(iR,iz)))
                      end do ! Levels
                    end do ! Terms

                  end if ! Height bounds

                end do ! Heights

              end if ! Specific or everything

              ! Write buffer
              call MPI_FILE_WRITE(funit,buffer(1), &
                                  Input%lim_pop%nbuff(ia)/4, &
                                  MPI_REAL,MPI_STATUS_IGNORE,ierr)
              if (ierr.ne.0) goto 1304

            end if! Saving populations

            ! Departure coefficients
            if (saveD.and.Input%lim_pop%nbuff(ia).gt.0) then

              ! Open file to write the populations
              call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.dep', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
              if (ierr.ne.0) goto 1005

              !
              ! Column offset
              !

              ! Get offset
              loffset = dble(icoords(3)-1)* &
                        dble(Input%lim_pop%nbuff(ia)) + &
                        dble(Input%lim_pop%head_size)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
                if (ierr.ne.0) goto 1015
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1015

              ! Initialize buffer
              ii = 0

              ! If specified
              if (Input%lim_pop%nran.gt.0) then

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Advance buffer
                      ii = ii +1
                      buffer(ii) = real(Atom(ia)%popu(i,iz)/ &
                                        Atom(ia)%populte(i,iz))

                    end do ! Ranges to print

                  ! In bounds
                  else

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Get necessary data
                      it = Atom(ia)%term(i)
                      iJ = Atom(ia)%sublevel(i)
                      rJ = Atom(ia)%rJval(iJ,it)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Advance buffer
                      ii = ii +1
                      buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz))/ &
                                        Atom(ia)%populte(i,iz))

                    end do ! Ranges to print

                  end if ! Height bounds

                end do ! Height

              ! Everything
              else

                ! Each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each term
                    do it=1,Atom(ia)%nlevel

                      ! Write populations
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(it,iz)/ &
                                        Atom(ia)%populte(i,iz))

                    end do ! Levels

                  ! In bounds
                  else

                    ! For each term
                    do it=1,Atom(ia)%nMulti

                      ! For each level J within the term
                      do iJ=1,Atom(ia)%nJ(it)!,1,-1

                        ! Get J value and level index
                        rJ = Atom(ia)%rJval(iJ,it)
                        i = Atom(ia)%irho(it)%irho_ij(iJ)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Write populations
                        ii = ii + 1
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))/ &
                                         Atom(ia)%populte(i,iz))

                      end do ! Levels
                    end do ! Terms

                  end if ! Height bounds

                end do ! Heights

              end if ! Specific or everything

              ! Write buffer
              call MPI_FILE_WRITE(funit,buffer(1), &
                                  Input%lim_pop%nbuff(ia)/4, &
                                  MPI_REAL,MPI_STATUS_IGNORE,ierr)
              if (ierr.ne.0) goto 1305

            end if! Saving populations

          end if! 1D vs 1.5D
        end if ! Saving populations or departure coefficients

        ! Write the population of the atom into solution
        if (saveSol) write (200,err=1100) Atom(ia)%n

        ! Saving rhoKQ
        if (saverKQ) then

          ! Write the population of the atom into rhoKQ files
          write (400,err=1103) Atom(ia)%n

          ! Write number of terms to rhoKQ file
          write (400,err=1103) Atom(ia)%nMulti

        end if

        ! For each term
        do it=1,Atom(ia)%nMulti

          ! Write number of levels to rhoKQ file
          if (saverKQ) write(400,err=1103) Atom(ia)%nJ(it)

          ! For each level J within the term
          do iJ=1,Atom(ia)%nJ(it)!,1,-1

            ! Get J value
            rJ = Atom(ia)%rJval(iJ,it)

            ! For each J value within the term
            do iJ1=1,Atom(ia)%nJ(it)!,1,-1

              ! Get J value
              rJ1 = Atom(ia)%rJval(iJ1,it)

              ! Write J values in rhoKQ file
              if (saverKQ) then
                write (400,err=1103) nint(2d0*rJ)
                write (400,err=1103) nint(2d0*rJ1)
              end if

              ! For each K
              do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                !
                ! Rotate rhoKQ nto the vertical reference frame
                if (K.le.Atom(ia)%Kcut(it)) then

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      rhoKQaux(-K:K,iz) = cZero

                    ! In bounds
                    else

                      ! For each Q
                      do iQ=-K,K

                        ! Get the rhoKQ index
                        iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                        ! Store the corresponding rhoKQ into an
                        ! auxiliar variable
                        rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                      end do ! Q

                      ! If there is non-zero magnetic field, rotate
                      if (Bfield%Bstrength(iz).gt.TINYB) &
                        call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                 -Bfield%Btheta(iz), &
                                 -Bfield%Bphi(iz),-1)

                    end if ! Height bounds

                  end do ! heights

                  ! For each Q
                  do iQ=-K,K

                    ! Get the index
                    iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                    ! For each height
                    do iz=1,nZ

                      ! Write rhoKQ into rhoKQ file, and the null flag
                      if (saverKQ) then
                        write (400,err=1103) dble(rhoKQaux(iQ,iz))
                        write (400,err=1103) dimag(rhoKQaux(iQ,iz))
                        ! If out of bounds
                        if (iz.lt.Rz0.or.iz.gt.Rz1) then
                          write (400,err=1103) ione
                        ! In bounds
                        else
                          if (Atom(ia)%rhonull(iR,iz)) then
                            write (400,err=1103) ione
                          else
                            write (400,err=1103) izero
                          end if
                        end if ! Height bounds
                      end if

                      ! Write rhoKQ into solution file
                      if (saveSol) &
                        write (200,err=1100) dble(rhoKQaux(iQ,iz)), &
                                             dimag(rhoKQaux(iQ,iz))

                    end do ! heights
                  end do ! Q

                else

                  ! For each Q
                  do iQ=-K,K

                    ! For each height
                    do iz=1,nZ

                      ! Write rhoKQ into rhoKQ file, and the null flag
                      if (saverKQ) then
                        write (400,err=1103) dzero
                        write (400,err=1103) dzero
                        write (400,err=1103) ione
                      end if
                      ! Write rhoKQ into solution file
                      if (saveSol) &
                        write (200,err=1100) dzero,dzero

                    end do ! heights
                  end do ! Q

                end if ! K<=Kcut

              end do ! K
            end do ! J'
          end do ! J
        end do ! Terms
      end do ! Atoms

      end if


      !
      ! JKQ

      ! Only if saving anything
      if (saveSol.or.saveJKQ) then

      ! For each atom
      do ia=1,nA

        ! Write number of transitions in JKQ file
        if (saveJKQ) write(300,err=1102) Atom(ia)%ntran

        ! For each transition
        do itran=1,Atom(ia)%ntran

          ! Apply atomic shift
          jtran = itran + Atom(ia)%tshift

          !
          ! Rotate JKQ

          ! For each height
          do iz=1,nz

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              JKQaux(:,:,iz) = cZero

            ! In bounds
            else

              ! get the jkq into an auxiliar variable
              JKQaux(:,:,iz) = JKQ(:,:,jtran,iz)

              ! If there is a non-zero magnetic field, rotate
              ! into the vertical reference frame
              if (Bfield%Bstrength(iz).gt.TINYB) &
              call fieldB(JKQaux(:,:,iz),1,Flgsg,-Bfield%Btheta(iz), &
                          -Bfield%Bphi(iz),-1)

            end if ! Height bounds

          end do ! heights

          ! For each K
          do K=0,2

            ! For each Q
            do iQ=-K,K

              ! For each height
              do iz=1,nZ

                ! Write the JKQ into the JKQ file
                if (saveJKQ) &
                  write(300,err=1102) dble(JKQaux(iQ,K,iz)), &
                                      dimag(JKQaux(iQ,K,iz))

                ! Write the JKQ into the solution file
                if (saveSol) &
                  write(200,err=1100) dble(JKQaux(iQ,K,iz)), &
                                      dimag(JKQaux(iQ,K,iz))

              end do ! heights
            end do ! Q
          end do ! K
        end do ! transitions
      end do ! Atoms

      ! If there is stimulated emission, write JKQS
      if(stm)then

        ! For each atom
        do ia=1,nA

          ! For each transition
          do itran=1,Atom(ia)%ntran

            ! Apply atomic shift
            jtran = itran + Atom(ia)%tshift

            !
            ! Rotate JKQ

            ! For each height
            do iz=1,nz

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                JKQaux(:,:,iz) = cZero

              ! In bounds
              else

                ! get the jkq into an auxiliar variable
                JKQaux(:,:,iz) = JKQS(:,:,jtran,iz)

                ! If there is a non-zero magnetic field, rotate
                ! into the vertical reference frame
                if (Bfield%Bstrength(iz).gt.TINYB) &
                call fieldB(JKQaux(:,:,iz),1,Flgsg, &
                            -Bfield%Btheta(iz), &
                            -Bfield%Bphi(iz),-1)

              end if ! Height bounds

            end do ! heights

            ! For each K
            do K=0,2

              ! For each Q
              do iQ=-K,K

                ! For each height
                do iz=1,nZ

                  ! Write the JKQS into the JKQ file
                  if (saveJKQ) &
                    write (300,err=1102) dble(JKQaux(iQ,K,iz)), &
                                         dimag(JKQaux(iQ,K,iz))

                  ! Write the JKQS into the solution file
                  if (saveSol) &
                    write (200,err=1100) dble(JKQaux(iQ,K,iz)), &
                                         dimag(JKQaux(iQ,K,iz))

                end do ! heights
              end do ! Q
            end do ! K
          end do ! transitions
        end do ! atoms

      end if ! stimulated emission
      end if ! Saving anything


      !
      ! Radiation field solution
      !

      ! Saving solution
      if (saveSol) then

      !
      ! Keeping full Stokes
      if (KSTK) then

        ! For each height
        do iz=1,nZ

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each degree of freedom
            do ith=1,Geom%nTh
            do iph=1,Geom%nPh
            do ifreq=1,nfreq
            do iS=0,3

              ! Write 0
              write (200,err=1100) dzero

            end do
            end do
            end do
            end do

          ! In bounds
          else

            ! For each polar direction
            do ith=1,Geom%nTh

              ! For each azimuthal direction
              do iph=1,Geom%nPh

                ! For each frequency
                do ifreq=1,nfreq

                  ! For each Stokes parameter
                  do iS=0,3

                    ! Write Stokes parameter into the solution file
                    write (200,err=1100) Stokes(iS,ifreq,iph,ith,iz)

                  end do ! Stokes parameters
                end do ! Frequencies
              end do ! azimuthal directions
            end do ! polar directions

          end if ! Height bounds

        end do ! heights

      ! Not keeping Stokes
      else

        ! For each height
        do iz=1,nZ

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each degree of freedom
            do ifreq=1,nfreq
            do K=0,2
            do iQ=-K,K

              ! Write 0 0
              write (200,err=1100) dzero,dzero

            end do
            end do
            end do

          ! In bounds
          else

            ! For each frequency
            do ifreq=1,nfreq

              ! For each K
              do K=0,2

                ! For each Q
                do iQ=-K,K

                  ! Write the JKQ(k) into the solution file
                  write (200,err=1100) dble(JKQC(iQ,K,ifreq,iz)), &
                                       dimag(JKQC(iQ,K,ifreq,iz))

                end do ! Q
              end do ! K
            end do ! frequencies

          end if ! Height bounds

        end do ! heights

      end if ! AV or AD

      !
      ! Close files
      !
      close (200)

      end if ! Saving solution

      if (saveJKQ) close (300)
      if (saverKQ) close (400)

      ! Check if storing anything else
      if (.not.saveJKQnu.and..not.saveS) then
        call control
        return
      end if


      !
      ! Store the stokes in the quadrature
      !

      ! If there is no suffix
      if (suff.eq.'NONE') then

        ! Open file
        if (saveS) then
          ! If 1D
          if (run_mode.eq.0) then
            open (250,file=trim(filename)//'/Stokesout', &
                  status='unknown', iostat=ios, err=1001, &
                  access='stream', action='write', &
                  form='unformatted')
          ! If 1.5D
          else if (run_mode.eq.1) then
            open (250,file=trim(filename)// &
                  '/Solution-folder/Stokesout-'//scoord, &
                  status='unknown', iostat=ios, err=1001, &
                  access='stream', action='write', &
                  form='unformatted')
          end if
        end if

        ! JKQnu file
        if (saveJKQnu) then
          ! If 1D
          if (run_mode.eq.0) then
            open (350,file=trim(filename)//'/JKQnuout', &
                  status='unknown', iostat=ios, err=1006, &
                  access='stream', action='write', &
                  form='unformatted')
          ! If 1.5D
          else if (run_mode.eq.1) then
            open (350,file=trim(filename)// &
                  '/Solution-folder/JKQnuout-'//scoord, &
                  status='unknown', iostat=ios, err=1006, &
                  access='stream', action='write', &
                  form='unformatted')
          end if
        end if

      ! If there is suffix
      else

        ! Open file
        if (saveS) &
          open (250,file=trim(filename)//'/Stokesout_'//suff, &
                status='unknown', iostat=ios, err=1001, &
                access='stream', action='write', &
                form='unformatted')

        ! JKQnu file
        if (saveJKQnu) &
          open (350,file=trim(filename)//'/JKQnuout_'//suff, &
                status='unknown', iostat=ios, err=1006, &
                access='stream', action='write', &
                form='unformatted')

      end if ! suffix


      !
      ! Write flag and dimensions
      !
      if (saveS) then
        write(250,err=1101) 'bo'
        write(250,err=1101) Nfreq
      end if
      if (saveJKQnu) then
        write(350,err=1106) 'ko'
        write(350,err=1106) nz
        write(350,err=1106) nfreq
      end if

      !
      ! Write data
      !

      ! Frequency axis
      if (saveS) write(250,err=1101) omega
      if (saveJKQnu) write(350,err=1106) omega

      ! Stokes out file
      if (saveS) then

        ! Number of directions
        write(250,err=1101) Geom%nTh/2,Geom%nPh

        ! For each polar direction
        do ith=1,Geom%nTh

          ! Ignore the ones going down
          if (Geom%V_mu(ith).lt.0) cycle

          ! For each azimuthal direction
          do iph=1,Geom%nPh

            ! Write the angles (DEG) of this quadrature direction
            write(250,err=1101) Geom%V_theta(ith)*180D0/pi, &
                                Geom%V_phi(iph)*180D0/pi

            ! Write the emergent Stokes parameters
            write(250,err=1101) transpose(Stokes(:,:,iph,ith,giz0))

          end do ! azimuthal directions
        end do ! polar directions

      end if ! Stokes out file

      ! JKQ(nu) file
      if (saveJKQnu) then

        ! For each height
        do iz=1,nz

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each degree of freedom
            do ifreq=1,nfreq
            do K=0,2
            do iQ=-K,K

              ! Write 0
              write(350,err=1106) dzero,dzero

            end do
            end do
            end do

          ! In bounds
          else

            ! For each frequency
            do ifreq=1,nfreq
              ! For each K multipole
              do K=0,2
                ! For each Q multipole
                do iQ=-K,K

                  ! Write
                  write(350,err=1106) JKQC(iQ,K,ifreq,iz)

                end do ! Q
              end do ! K
            end do ! Frequency

          end if ! Height bounds

        end do ! Height
      end if ! JKQ(nu) file

      !
      ! Close the Stokesout file
      !
      if (saveS) close(250)
      if (saveJKQnu) close(350)

      ! Control
      call control

      return

1000  umsg = 'Error opening solution file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing solution file'
      close(200)
      inquire(unit=300, opened=laux)
      if (laux) close(300)
      inquire(unit=400, opened=laux)
      if (laux) close(400)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening Stokesout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing Stokesout file'
      close(250)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1002  umsg = 'Error opening Jout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1102  umsg = 'Error writing Jout file'
      close(300)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1003  umsg = 'Error opening Rout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1103  umsg = 'Error writing Rout file'
      close(400)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1004  umsg = 'Error opening Population file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1014  umsg = 'Error seeking Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1104  umsg = 'Error writing Population file'
      close(500)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1304  umsg = 'Error writing Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1005  umsg = 'Error opening Departure file'
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1015  umsg = 'Error seeking Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1105  umsg = 'Error writing Departure file'
      close(600)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1305  umsg = 'Error writing Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1006  umsg = 'Error opening JKQnuout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1106  umsg = 'Error writing JKQnuout file'
      close(350)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writesol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the computed solution that can be read
      !! by the code as initialization, and separate files for the
      !! radiation field and density tensors. This is the version
      !! for only intensity.\n
      !!        Input(Input_class): Structure with settings data\n
      !!        suff(character(:)): Suffix for the file names of the
      !!                            radiation field and density matrix
      !!                            files\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!              z(dfloat(:)): Heights array\n
      !!   Stokes(dfloat(:,:,:,:)): Intensity\n
      !!          J00(dfloat(:,:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         J00S(dfloat(:,:)): Mean intensity integrated over
      !!                            emission profile\n
      !!         J00C(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(dfloat(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      !!             keep(logical): Determines the name of the
      !!                            Solution file
      subroutine writesolI(Input,suff,omega,Geom, &
                           Atom,z,Stokes,J00,J00S,J00C,J00P,keep)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(in):: Geom
      character(len=4), intent(in):: suff
      logical, intent(in):: keep
      double precision, dimension(:), intent(in):: omega,z
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        intent(in):: Stokes
      double precision,dimension(nxt,Rz0:Rz1), intent(in):: J00, J00S
      double precision,dimension(nfreq,Rz0:Rz1), intent(in):: J00C
      double precision,dimension(nxphot,2,Rz0:Rz1), intent(in):: J00P

      ! Local

      integer, parameter:: izero = 0
      integer, parameter:: ione = 1
      double precision, parameter:: dzero = 0d0

      character(len=8):: scoord
      character(len=500):: filename

      logical:: saveJ00nu, saverKQ, saveJKQ
      logical:: saveP, saveD, saveS, saveSol
      logical:: laux

      integer:: ierr,ii,iab
      integer:: ios,ia,iz,ifreq,i,ith,iph,itran,jtran
      integer:: K,iQ,iR,it,iJ,iran
      integer:: axial_int,stm_int,AV_int


      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: rJ,loffset


      !
      ! Slaves just wait
      !
      if (pid.gt.0) then
        call control
        return
      end if

      !
      ! Translate
      !
      filename = Input%folder
      saveSol = Input%keep_sol
      saveP = Input%keep_pop.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveD = Input%keep_dep.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saverKQ = Input%keep_rhoKQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveJKQ = Input%keep_JKQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveS = Input%keep_stokesQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveJ00nu = Input%keep_jkqnu.and. &
                  (suff.eq.'NONE'.or.run_mode.eq.0)

      ! If not writing anything, come back
      if (.not.(saveSol.or.saveP.or.saveD.or.saverKQ.or.saveJKQ.or. &
                saveS.or.saveJ00nu)) then
        call control
        return
      end if

      ! Routine name
      urou = 'writesolI'

      !
      ! Open files
      !

      ! Writing solution
      if (saveSol) then

        ! If 1D or inversion
        if (run_mode.le.0) then

          ! To write the solution
          if (keep) then
            open (200,file=trim(filename)//'/SolutionI', &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')
          else
            open (200,file=trim(filename)//'/Solution', &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')
          end if

        ! If 1.5D
        else if (run_mode.eq.1) then

          ! Get LOS index
          write(scoord,'(I0.8)') icoords(3)

          ! To write the solution
          if (keep) then
            open (200,file=trim(filename)// &
                  '/Solution-folder/SolutionI-'//scoord, &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')
          else
            open (200,file=trim(filename)// &
                  '/Solution-folder/Solution-'//scoord, &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')
          end if

        end if ! 1D vs 1.5D
      end if ! Saving solution file

      ! If saving JKQ
      if (saveJKQ) then
        ! To write the final JKQ
        if (suff.eq.'NONE') then
          ! If 1D or inversion
          if (run_mode.le.0) then
            if (keep) then
              open (300,file=trim(filename)//'/JoutI', &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              open (300,file=trim(filename)//'/Jout', &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')
            end if ! Keep intensity sol.
          ! If 1.5D
          else if (run_mode.eq.1) then
            if (keep) then
              open (300,file=trim(filename)// &
                    '/Solution-folder/JoutI-'//scoord, &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              open (300,file=trim(filename)// &
                    '/Solution-folder/Jout-'//scoord, &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')
            end if ! Keep intentisy sol.
          end if ! 1D vs 1.5D
        ! Not final JKQ
        else
          open (300,file=trim(filename)//'/Jout_'//suff, &
                status='unknown', iostat=ios, err=1002, &
                access='stream', action='write', &
                form='unformatted')
        end if ! Final
      end if ! Saving JKQ

      ! If saving rhoKQ
      if (saverKQ) then
        ! To write the final rhoKQ
        if (suff.eq.'NONE') then
          ! If 1D or inversion
          if (run_mode.le.0) then
            if (keep) then
              open (400,file=trim(filename)//'/RhooutI', &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              open (400,file=trim(filename)//'/Rhoout', &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')
            end if ! Keep intensity sol.
          ! If 1.5D
          else if (run_mode.eq.1) then
            if (keep) then
              open (400,file=trim(filename)// &
                    '/Solution-folder/RhooutI-'//scoord, &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              open (400,file=trim(filename)// &
                    '/Solution-folder/Rhoout-'//scoord, &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')
            end if ! Keep intensity sol.
          end if ! 1D vs 1.5D
        ! Not final rhoKQ
        else
          open (400,file=trim(filename)//'/Rhoout_'//suff, &
                status='unknown', iostat=ios, err=1003, &
                access='stream', action='write', &
                form='unformatted')
        end if
      end if ! Saving rhoKQ

      !
      ! Convert logicals to integers
      !

      ! Axial symmetry
      if(Geom%axial)then
          axial_int = 1
      else
          axial_int = 0
      end if

      ! Stimulated emission
      if(stm)then
          stm_int = 1
      else
          stm_int = 0
      end if

      ! Angle averaged redistribution function
      if (KSTK) then
          AV_int = 0
      else
          AV_int = 1
      end if

      !
      ! Write headers with dimensions and flags
      !

      ! If saving solution
      if (saveSol) then

        ! Solution file
        write (200,err=1100) 'si'
        write (200,err=1100) nfreq
        write (200,err=1100) nZ
        write (200,err=1100) Geom%nTh
        write (200,err=1100) Geom%nPh
        write (200,err=1100) nA
        write (200,err=1100) axial_int
        write (200,err=1100) stm_int
        write (200,err=1100) AV_int

      end if

      ! If saving JKQ
      if (saveJKQ) then

        ! JKQ file
        write(300,err=1102) 'bj'
        write(300,err=1102) stm
        write(300,err=1102) nZ
        write(300,err=1102) nA
        write(300,err=1102) nxt
        write(300,err=1102) z

      end if

      ! If saving rhoKQ
      if (saverKQ) then

        ! rhoKQ file
        write(400,err=1103) 'br'
        write(400,err=1103) nZ
        write(400,err=1103) nA
        write(400,err=1103) z

      end if

      ! If 1.5D, prepare pop y dep buffers
      if (run_mode.eq.1) then
        if (saveP.or.saveD) &
          allocate(buffer(maxval(Input%lim_pop%nbuff)/4))
      end if

      !
      ! Write the data
      !

      ! If saving anything
      if (saveSol.or.saverKQ.or.saveP.or.saveD) then

      !
      ! Population and rhoKQ

      ! For each atom
      do ia=1,nA

        !
        ! If saving population or departure coeff.
        !
        if (saveP.or.saveD) then

          !
          ! If 1D
          !
          if (run_mode.eq.0) then

            ! If writing populations
            if (saveP) then

              ! To write the populations
              open (500,file=trim(filename)//'/'// &
                    trim(Atom(ia)%file_label)//'.pop', &
                    status='unknown',iostat=ios, err=1004, &
                    access='stream', action='write', &
                    form='unformatted')

              write(500,err=1104) 'bp'
              write(500,err=1104) nZ
              write(500,err=1104) Atom(ia)%nlevel

            end if

            ! If writing departure c.
            if (saveD) then

              ! To write departure c.
              open (600,file=trim(filename)//'/'// &
                    trim(Atom(ia)%file_label)//'.dep', &
                    status='unknown',iostat=ios, err=1005, &
                    access='stream', action='write', &
                    form='unformatted')

              write(600,err=1105) 'bb'
              write(600,err=1105) nZ
              write(600,err=1105) Atom(ia)%nlevel

            end if

            ! For each height
            do iz=1,nZ

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! For each term
                do i=1,Atom(ia)%nlevel

                  ! Write populations
                  if (saveP) &
                    write (500,err=1104) Atom(ia)%popu(i,iz)

                  ! Write departure coeff
                  if (saveD) &
                    write (600,err=1105) Atom(ia)%popu(i,iz)/ &
                                         Atom(ia)%populte(i,iz)

                end do ! Levels

              ! If in bounds
              else

                ! For each term
                do i=1,Atom(ia)%nlevel

                  ! Get J
                  rJ = Atom(ia)%rJval(Atom(ia)%sublevel(i), &
                       Atom(ia)%term(i))
                  iR = Atom(ia)%irho(Atom(ia)%term(i))% &
                                Jrho(Atom(ia)%sublevel(i), &
                                     Atom(ia)%sublevel(i))%kq(0,0)

                  ! Write populations
                  if (saveP) &
                    write (500,err=1104) Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))

                  ! Write departure coeff
                  if (saveD) &
                    write (600,err=1105) Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))/ &
                                         Atom(ia)%populte(i,iz)

                end do ! Levels

              end if ! Height bounds

            end do ! Heights

            if (saveP) close(500)
            if (saveD) close(600)

          !
          ! If 1.5D
          !
          else if (run_mode.eq.1) then

            ! Populations
            if (saveP.and.Input%lim_pop%nbuff(ia).gt.0) then

              ! Open file to write the populations
              call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.pop', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
              if (ierr.ne.0) goto 1004

              !
              ! Column offset
              !

              ! Get offset
              loffset = dble(icoords(3)-1)* &
                        dble(Input%lim_pop%nbuff(ia)) + &
                        dble(Input%lim_pop%head_size)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
                if (ierr.ne.0) goto 1014
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1014

              ! Initialize buffer
              ii = 0

              ! If specified
              if (Input%lim_pop%nran.gt.0) then

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Advance buffer
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(i,iz))
                    end do ! Ranges to print

                  ! In bounds
                  else

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Get necessary data
                      it = Atom(ia)%term(i)
                      iJ = Atom(ia)%sublevel(i)
                      rJ = Atom(ia)%rJval(iJ,it)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Advance buffer
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz)))
                    end do ! Ranges to print

                  end if ! Height bounds

                end do ! Heights

              ! Everything
              else

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each level
                    do i=1,Atom(ia)%nlevel

                      ! Write populations
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(i,iz))
                    end do ! Levels

                  ! In bounds
                  else

                    ! For each level
                    do i=1,Atom(ia)%nlevel

                      ! Get necessary data
                      it = Atom(ia)%term(i)
                      iJ = Atom(ia)%sublevel(i)
                      rJ = Atom(ia)%rJval(iJ,it)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Write populations
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz)))
                    end do ! Levels

                  end if ! Height bounds

                end do ! Heights

              end if ! Specific or everything

              ! Write buffer
              call MPI_FILE_WRITE(funit,buffer(1), &
                                  Input%lim_pop%nbuff(ia)/4, &
                                  MPI_REAL,MPI_STATUS_IGNORE,ierr)
              if (ierr.ne.0) goto 1304

            end if! Saving populations

            ! Departure coefficients
            if (saveD.and.Input%lim_pop%nbuff(ia).gt.0) then

              ! Open file to write the populations
              call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.dep', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
              if (ierr.ne.0) goto 1005

              !
              ! Column offset
              !

              ! Get offset
              loffset = dble(icoords(3)-1)* &
                        dble(Input%lim_pop%nbuff(ia)) + &
                        dble(Input%lim_pop%head_size)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
                if (ierr.ne.0) goto 1015
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1015

              ! Initialize buffer
              ii = 0

              ! If specified
              if (Input%lim_pop%nran.gt.0) then

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Advance buffer
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(i,iz)/ &
                                        Atom(ia)%populte(i,iz))
                    end do ! Ranges to print

                  ! In bounds
                  else

                    ! For each entry to write
                    do iran=1,Input%lim_pop%nran

                      ! Atom and transition
                      iab = Input%lim_pop%indx(1,iran)
                      if (ia.ne.iab) cycle
                      i = Input%lim_pop%indx(2,iran)

                      ! Get necessary data
                      it = Atom(ia)%term(i)
                      iJ = Atom(ia)%sublevel(i)
                      rJ = Atom(ia)%rJval(iJ,it)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Advance buffer
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%n(ii)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz))/ &
                                        Atom(ia)%populte(i,iz))
                    end do ! Ranges to print

                  end if ! Height bounds

                end do ! Heights

              ! Everything
              else

                ! For each height
                do iz=1,nz

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! For each level
                    do i=1,Atom(ia)%nlevel

                      ! Write populations
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%popu(i,iz)/ &
                                        Atom(ia)%populte(i,iz))
                    end do ! Levels

                  ! In bounds
                  else

                    ! For each level
                    do i=1,Atom(ia)%nlevel

                      ! Get necessary data
                      it = Atom(ia)%term(i)
                      iJ = Atom(ia)%sublevel(i)
                      rJ = Atom(ia)%rJval(iJ,it)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Write populations
                      ii = ii + 1
                      buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz))/ &
                                        Atom(ia)%populte(i,iz))
                    end do ! Levels

                  end if ! Height bounds

                end do ! Heights

              end if ! Specific or everything

              ! Write buffer
              call MPI_FILE_WRITE(funit,buffer(1), &
                                  Input%lim_pop%nbuff(ia)/4, &
                                  MPI_REAL,MPI_STATUS_IGNORE,ierr)
              if (ierr.ne.0) goto 1305

            end if! Saving populations

          end if! 1D vs 1.5D
        end if ! Saving populations or departure coefficients

        ! Write the population of the atom into solution file
        if (saveSol) write (200,err=1100) Atom(ia)%n

        ! Saving rhoKQ
        if (saverKQ) then

          ! Write the population of the atom into rhoKQ file
          write (400,err=1100) Atom(ia)%n

          ! Write number of terms to rhoKQ file
          write (400,err=1103) Atom(ia)%nlevel

        end if

        ! For each term
        do i=1,Atom(ia)%nlevel

          ! Write number of levels to rhoKQ file
          if (saverKQ) write (400,err=1103) 1

          ! Get J value and level index
          iR = Atom(ia)%irho(Atom(ia)%term(i))% &
                        Jrho(Atom(ia)%sublevel(i), &
                             Atom(ia)%sublevel(i))%kq(0,0)

          ! Get J
          rJ = Atom(ia)%rJval(Atom(ia)%sublevel(i),Atom(ia)%term(i))

          ! Write J values in rhoKQ file
          if (saverKQ) then
            write (400,err=1103) nint(2d0*rJ)
            write (400,err=1103) nint(2d0*rJ)
          end if

          ! For each K
          do K=0,nint(2d0*rJ)

            ! For each Q
            do iQ=-K,K

              if (K.lt.1) then

                ! For each height
                do iz=1,nZ

                  ! If out of bounds
                  if (iz.lt.Rz0.or.iz.gt.Rz1) then

                    ! Write rhoKQ into rhoKQ file, and the null flag
                    if (saverKQ) then
                      write (400,err=1103) dzero,dzero
                      write (400,err=1103) ione
                    end if

                    ! Write rhoKQ into solution file
                    if (saveSol) &
                      write (200,err=1100) dzero

                  ! In bounds
                  else

                    ! Write rhoKQ into rhoKQ file, and the null flag
                    if (saverKQ) then
                      write (400,err=1103) dble(Atom(ia)%crho(iR,iz))
                      write (400,err=1103) dzero
                      if (Atom(ia)%rhonull(iR,iz)) then
                        write (400,err=1103) ione
                      else
                        write (400,err=1103) izero
                      end if
                    end if

                    ! Write rhoKQ into solution file
                    if (saveSol) &
                      write (200,err=1100) dble(Atom(ia)%crho(iR,iz))

                  end if ! Height bounds

                end do ! heights

              else

                ! For each height
                do iz=1,nZ

                  ! Write rhoKQ into rhoKQ file, and the null flag
                  if (saverKQ) then
                    write (400,err=1103) dzero
                    write (400,err=1103) dzero
                    write (400,err=1103) ione
                  end if

                end do ! heights

              end if

            end do ! Q
          end do ! K
        end do ! Levels
      end do ! Atoms

      end if ! Saving rho or Solution


      !
      ! J00

      ! If saving anything
      if (saveSol.or.saveJKQ) then

      ! For each atom
      do ia=1,nA

        ! Write number of transitions in JKQ file
        if (saveJKQ) write(300,err=1102) Atom(ia)%nftran

        ! For each transition
        do itran=1,Atom(ia)%nftran

          ! Apply atomic shift
          jtran = itran + Atom(ia)%tfshift

          ! For each height
          do iz=1,nZ

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! Write the JKQ into the JKQ file
              if (saveJKQ) write(300,err=1102) dzero,dzero

              ! Write the JKQ into the solution file
              if (saveSol) &
                write(200,err=1100) dzero

            ! In bounds
            else

              ! Write the JKQ into the JKQ file
              if (saveJKQ) write(300,err=1102) J00(jtran,iz),dzero

              ! Write the JKQ into the solution file
              if (saveSol) &
                write(200,err=1100) J00(jtran,iz)

            end if ! Height bounds

          end do ! heights

          ! For each K
          do K=1,2

            ! For each Q
            do iQ=-K,K

              ! For each height
              do iz=1,nZ

                ! Write the JKQ into the JKQ file
                if (saveJKQ) write(300,err=1102) dzero,dzero

              end do ! heights
            end do ! Q
          end do ! K
        end do ! transitions
      end do ! Atoms

      ! If there is stimulated emission, write JKQS
      if(stm)then

        ! For each atom
        do ia=1,nA

          ! For each transition
          do itran=1,Atom(ia)%nftran

            ! Apply atomic shift
            jtran = itran + Atom(ia)%tfshift

            ! For each height
            do iz=1,nZ

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! Write the JKQ into the JKQ file
                if (saveJKQ) write(300,err=1102) dzero,dzero

                ! Write the JKQ into the solution file
                if (saveSol) &
                  write(200,err=1100) dzero

              ! In bounds
              else

                ! Write the JKQ into the JKQ file
                if (saveJKQ) write(300,err=1102) J00S(jtran,iz),dzero

                ! Write the JKQ into the solution file
                if (saveSol) &
                  write(200,err=1100) J00S(jtran,iz)

              end if ! Height bounds

            end do ! heights

            ! Saving JKQ
            if (saveJKQ) then

              ! For each K
              do K=1,2

                ! For each Q
                do iQ=-K,K

                  ! For each height
                  do iz=1,nZ

                    ! Write the JKQS into the JKQ file
                    write (300,err=1102) dzero,dzero

                  end do ! heights
                end do ! Q
              end do ! K

            end if

          end do ! transitions
        end do ! atoms

      end if ! stimulated emission

      end if ! Saving anything

      !
      ! J00P

      ! Saving solution
      if (saveSol) then

      ! For each atom
      do ia=1,nA

        ! For each transition
        do itran=1,Atom(ia)%nphot

          ! Apply atomic shift
          jtran = itran + Atom(ia)%pshift

          ! For each height
          do iz=1,nZ

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              write(200,err=1100) dzero,dzero

            ! In bounds
            else

              ! Write the J00P into the solution file
              write(200,err=1100) J00P(jtran,1,iz)
              write(200,err=1100) J00P(jtran,2,iz)

            end if ! Height bounds

          end do ! heights

        end do ! transitions
      end do ! Atoms


      !
      ! Radiation field solution
      !

      !
      ! If not angle-averaged
      if (KSTK) then

        ! For each height
        do iz=1,nZ

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each degree of freedom
            do ith=1,Geom%nTh
            do iph=1,Geom%nPh
            do ifreq=1,nfreq

              ! Write 0
              write (200,err=1100) dzero

            end do
            end do
            end do

          ! In bounds
          else

            ! For each polar direction
            do ith=1,Geom%nTh

              ! For each azimuthal direction
              do iph=1,Geom%nPh

                ! For each frequency
                do ifreq=1,nfreq

                  ! Write Stokes parameter into the solution file
                  write (200,err=1100) Stokes(ifreq,iph,ith,iz)

                end do ! Frequencies
              end do ! azimuthal directions
            end do ! polar directions

          end if ! Height bounds

        end do ! heights

      !
      ! If Angle-averaged
      else

        ! For each height
        do iz=1,nZ

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each frequency
            do ifreq=1,nfreq

              ! Write 0
              write (200,err=1100) dzero

            end do ! frequencies

          ! In bounds
          else

            ! For each frequency
            do ifreq=1,nfreq

              ! Write the JKQ(k) into the solution file
              write (200,err=1100) J00C(ifreq,iz)

            end do ! frequencies

          end if ! Height bounds

        end do ! heights

      end if ! AV or AD


      !
      ! Close files
      !
      close (200)

      end if ! Saving Solution

      if (saveJKQ) close (300)
      if (saverKQ) close (400)

      ! Check if storing anything else
      if (.not.saveJ00nu.and..not.saveS) then
        call control
        return
      end if


      !
      ! Store the stokes in the quadrature
      !

      ! If there is no suffix
      if (suff.eq.'NONE') then

        ! Open Stokes file
        if (saveS) then
          ! If 1D
          if (run_mode.eq.0) then
            ! Keeping intensity sol.
            if (keep) then
              ! Open file
              open (250,file=trim(filename)//'/StokesoutI', &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              ! Open file
              open (250,file=trim(filename)//'/Stokesout', &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')
            end if
          ! If 1.5D
          else if (run_mode.eq.1) then
            ! Keeping intensity sol.
            if (keep) then
              ! Open file
              open (250,file=trim(filename)// &
                    '/Solution-folder/StokesoutI-'//scoord, &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              ! Open file
              open (250,file=trim(filename)// &
                    '/Solution-folder/Stokesout-'//scoord, &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')
            end if
          end if
        end if ! Saving Stokes

        ! Open JKQ(nu) file
        if (saveJ00nu) then
          ! If 1D
          if (run_mode.eq.0) then
            ! Keeping intensity solution
            if (keep) then
              open (350,file=trim(filename)//'/JKQnuoutI', &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              open (350,file=trim(filename)//'/JKQnuout', &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')

            end if
          ! If 1.5D
          else if (run_mode.eq.1) then
            ! Keeping intensity solution
            if (keep) then
              open (350,file=trim(filename)// &
                    '/Solution-folder/JKQnuoutI-'//scoord, &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')
            else
              open (350,file=trim(filename)// &
                    '/Solution-folder/JKQnuout-'//scoord, &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')

            end if
          end if
        end if ! Saving J00(nu)

      ! If there is suffix
      else

        ! Open file
        open (250,file=trim(filename)//'/Stokesout_'//suff, &
              status='unknown', iostat=ios, err=1001, &
              access='stream', action='write', &
              form='unformatted')

        ! JKQnu file
        if (saveJ00nu) &
          open (350,file=trim(filename)//'/JKQnuout_'//suff, &
                status='unknown', iostat=ios, err=1006, &
                access='stream', action='write', &
                form='unformatted')

      end if ! suffix


      !
      ! Write flag and dimensions
      !
      if (saveS) then
        write(250,err=1101) 'bo'
        write(250,err=1101) Nfreq
      end if
      if (saveJ00nu) then
        write(350,err=1106) 'ko'
        write(350,err=1106) nz
        write(350,err=1106) nfreq
      end if

      !
      ! Write data
      !

      ! Frequency axis
      if (saveS) write(250,err=1101) omega
      if (saveJ00nu) write(350,err=1106) omega

      ! Saving Stokes
      if (saveS) then

        ! Number of directions
        write(250,err=1101) Geom%nTh/2,Geom%nPh

        ! For each polar direction
        do ith=1,Geom%nTh

          ! Ignore the ones going down
          if (Geom%V_mu(ith).lt.0) cycle

          ! For each azimuthal direction
          do iph=1,Geom%nPh

            ! Write the angles (DEG) of this quadrature direction
            write(250,err=1101) Geom%V_theta(ith)*180D0/pi, &
                                Geom%V_phi(iph)*180D0/pi

            ! Write the emergent Stokes parameters
            write(250,err=1101) Stokes(:,iph,ith,giz0)

            do ifreq=1,nfreq*3
              write(250,err=1101) dzero
            end do

          end do ! azimuthal directions
        end do ! polar directions

      end if ! Saving Stokes

      ! Save J00
      if (saveJ00nu) then

        ! For each height
        do iz=1,nz

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each frequency
            do ifreq=1,nfreq
              write(350,err=1106) dzero
            end do

          ! In bounds
          else

            write(350,err=1106) J00C(:,iz)

          end if ! Height bounds

        end do ! Heights

      end if


      !
      ! Close the Stokesout file
      !
      if (saveS) close(250)
      if (saveJ00nu) close(350)

      ! Control
      call control

      return

1000  umsg = 'Error opening solution file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing solution file'
      close(200)
      inquire(unit=300, opened=laux)
      if (laux) close(300)
      inquire(unit=400, opened=laux)
      if (laux) close(400)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening Stokesout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing Stokesout file'
      close(250)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1002  umsg = 'Error opening Jout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1102  umsg = 'Error writing Jout file'
      close(300)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1003  umsg = 'Error opening Rout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1103  umsg = 'Error writing Rout file'
      close(400)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1004  umsg = 'Error opening Population file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1014  umsg = 'Error seeking Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1104  umsg = 'Error writing Population file'
      close(500)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1304  umsg = 'Error writing Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1005  umsg = 'Error opening Departure file'
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1015  umsg = 'Error seeking Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1105  umsg = 'Error writing Departure file'
      close(600)
      if (saveSol) close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1305  umsg = 'Error writing Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1006  umsg = 'Error opening JKQnuout file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1106  umsg = 'Error writing JKQnuout file'
      close(350)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writesolI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the emerging Stokes parameters.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!       Stokes(dfloat(:,:)): Stokes parameters\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writestk(filename,iph,ith,omega,Geom,Stokes,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega
      double precision,dimension(:,:), intent(in):: Stokes

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,iran,ii,jj,i0,i1,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writestk'


      !
      ! Convert the integers into appropriate length strings
      !
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D or inversion
      !
      if (run_mode.eq.0) then

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/Stokes_'//trim(cth)//'_'// &
             trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream', action='write', form='unformatted')


        !
        ! Write content
        !

        ! Identification
        write(200,err=1100) 'be'

        ! Number of frequencies, angles of LOS
        write(200,err=1100) nfreq,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Stokes vector
        write(200,err=1100) transpose(Stokes)

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Stokes_'//trim(cth)// &
                           '_'//trim(cph),MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each Stokes parameter
          do jj=1,4

            ! For each entry to write
            do iran=1,buff%nran

              ! Atom and transition
              i0 = buff%indx(1,iran)
              i1 = buff%indx(2,iran)
              nn = i1-i0+1

              ! Fill buffer
              buffer(ii+1:ii+nn) = real(reshape(Stokes(jj,i0:i1), &
                                                (/ nn /)))
              ! Advance index
              ii = ii + nn

            end do ! Ranges
          end do ! Stokes parameters

        ! All
        else

          ! Size
          nn = nfreq*4

          ! Fill buffer
          buffer = real(reshape(transpose(Stokes), (/ nn /)))

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      !
      ! If inversion
      !
      else if (run_mode.eq.-1) then

        ! Size
        nn = buff%buffer_size/4

        ! Allocate buffer
        allocate(buffer(nn))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Stokes', &
                           MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Fill buffer
        buffer = real(reshape(transpose(Stokes), (/ nn /)))

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      return

1000  umsg = 'Error opening Stokes file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking StokesI file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing Stokes file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writestk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the emerging intensity.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!         Stokes(dfloat(:)): Stokes parameters\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writestkI(filename,iph,ith,omega,Geom,StokesI,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: StokesI

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,iran,ii,i0,i1,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset
      double precision, dimension(nfreq,0:3):: Stokes


      ! Routine name
      urou = 'writestkI'


      !
      ! Convert the integers into appropriate length strings
      !
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        !
        ! Convert intensity into polarization array
        !
        Stokes(:,1:3) = 0d0
        Stokes(:,0) = StokesI

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/StokesI_'// trim(cth)//'_'// &
             trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream', action='write', form='unformatted')


        !
        ! Write content
        !

        ! Identification
        write(200,err=1100) 'be'

        ! Number of frequencies, angles of LOS
        write(200,err=1100) nfreq,Geom%L_theta(ith)*180d0/pi, &
                          Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Stokes vector
        write(200,err=1100) Stokes

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))
        buffer = 0e0

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/StokesI_'//trim(cth)// &
                           '_'//trim(cph),MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size)/4 + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            buffer(ii+1:ii+nn) = real(reshape(StokesI(i0:i1), &
                                               (/ nn /)))
            ! Advance index
            ii = ii + nn

          end do ! Ranges

        ! All
        else

          ! Fill buffer
          buffer(1:nfreq) = real(StokesI)

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/16, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      !
      ! If inversion
      !
      else if (run_mode.eq.-1) then

        ! Size
        nn = buff%buffer_size/16

        ! Allocate buffer
        allocate(buffer(nn))
        buffer = 0e0

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Stokes', &
                           MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size)/4 + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Fill buffer
        buffer = real(StokesI)

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      return

1000  umsg = 'Error opening StokesI file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking StokesI file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing StokesI file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing StokesI file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writestkI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes the emerging Stokes parameters to RAM.\n
      !!     e_Stokes(dfloat(:,:)): Output stokes parameters\n
      !!       Stokes(dfloat(:,:)): Stokes parameters\n
      !!  buff(IO_helper_class(:)): Info about what to store\n
      !!           lrange(logical): Account for buffer range
      subroutine setstk(e_Stokes,Stokes,buff,lrange)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      logical, intent(in):: lrange
      double precision,dimension(:,:), intent(in):: Stokes
      double precision,dimension(:,:), intent(out):: e_Stokes

      ! Local

      integer:: iran,ii,jj,i0,i1,nn


      ! Only master
      if (pid.gt.0) return

      ! If specified
      if (buff%nran.gt.0.and.lrange) then

        ! For each Stokes parameter
        do jj=1,4

          ! Initialize buffer
          ii = 0

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            e_Stokes(jj,ii+1:ii+nn) = Stokes(jj,i0:i1)

            ! Advance index
            ii = ii + nn

          end do ! Ranges
        end do ! Stokes parameters

      ! All
      else

        ! Copy
        e_Stokes = Stokes

      end if

      return

      end subroutine setstk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes the emerging intensity to RAM.\n
      !!     e_Stokes(dfloat(:,:)): Output stokes parameters\n
      !!       Stokes(dfloat(:,:)): Stokes parameters\n
      !!  buff(IO_helper_class(:)): Info about what to store\n
      !!           lrange(logical): Account for buffer range
      subroutine setstkI(e_Stokes,Stokes,buff,lrange)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      logical, intent(in):: lrange
      double precision,dimension(:), intent(in):: Stokes
      double precision,dimension(:,:), intent(out):: e_Stokes

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Only master
      if (pid.gt.0) return

      ! If specified
      if (buff%nran.gt.0.and.lrange) then

        ! Initialize buffer
        ii = 0

        ! For each entry to write
        do iran=1,buff%nran

          ! Atom and transition
          i0 = buff%indx(1,iran)
          i1 = buff%indx(2,iran)
          nn = i1-i0+1

          ! Fill buffer
          e_Stokes(1,ii+1:ii+nn) = Stokes(i0:i1)

          ! Advance index
          ii = ii + nn

        end do ! Ranges

      ! All
      else

        ! Copy
        e_Stokes(1,:) = Stokes

      end if

      return

      end subroutine setstkI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the emerging Stokes parameters and another
      !! with the optical depth for the CLE mode.\n
      !!    filename(character(:)): Name of the file to write\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!  buff(IO_helper_class(:)): Info about what to store\n
      !!             wtau(logical): If writing tau
      subroutine write_CLEgeom(filename,Atmo,buff,wtau)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      logical, intent(in):: wtau

      ! Local

      integer:: ierr

      integer(kind=MPI_OFFSET_KIND):: offset

      double precision, dimension(2):: buffer


      ! Cartesian does not need to write
      if (Atmo%mode.eq.0) return

      ! If no master, return
      if (pid.gt.0) then
        call control
        return
      end if

      ! Routine name
      urou = 'write_CLEgeom'

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Stokes', &
                         MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Get offset and data
      !

      ! If slab or non-cartesian
      if (Atmo%mode.eq.1.or.Atmo%mode.eq.2) then

        ! Compute offset
        offset = (icoords(3)-1)*16 + buff%head_size

        ! Slab
        if (Atmo%mode.eq.1) then

          buffer(1) = Atmo%z(1)
          buffer(2) = Atmo%ypos

        ! Non-cartesian
        else if (Atmo%mode.eq.2) then

          buffer(1) = Atmo%ypos
          buffer(2) = Atmo%zpos

        end if

      end if ! Atmo mode

      ! Go to offset
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Write buffer
      call MPI_FILE_WRITE(funit, buffer(1), 2, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_STATUS_IGNORE, ierr)
      if (ierr.ne.0) goto 1300

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      !
      ! Tau?
      !
      if (wtau) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Tau', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1100


        ! Go to offset
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1110

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1), 2, &
                            MPI_DOUBLE_PRECISION,&
                            MPI_STATUS_IGNORE, ierr)
        if (ierr.ne.0) goto 1310

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if ! Tau

      ! Control
      call control

      return

1000  umsg = 'Error opening Stokes file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error opening Tau file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1110  umsg = 'Error seeking Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1310  umsg = 'Error writing Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine write_CLEgeom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the emerging Stokes parameters and another
      !! with the optical depth for the CLE mode.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              if0(integer): Initial frequency index\n
      !!              if1(integer): Final frequency index\n
      !!               nf(integer): Number of frequencies\n
      !!       Stokes(dfloat(:,:)): Stokes parameters\n
      !!            tau(dfloat(:)): Optical depth\n
      !!  buff(IO_helper_class(:)): Info about what to store\n
      !!             wtau(logical): If writing tau
      subroutine write_CLE(filename,if0,if1,nf,Stokes,tau, &
                           buff,wtau)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      logical, intent(in):: wtau
      integer, intent(in):: if0,if1,nf
      double precision, dimension(:), intent(in):: tau
      double precision,dimension(:,:), intent(in):: Stokes

      ! Local

      integer:: ierr,is,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'write_CLE'

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Stokes', &
                         MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size+buff%geom_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      !
      ! If splitting frequencies
      !
      if (nproc.gt.1) then

        ! Allocate buffer
        nn = nf
        allocate(buffer(nn))

        ! For each Stokes parameter
        do is=1,4

          ! Initial wavelength shift
          if (if0.gt.1) then
            loffset = dble((if0-1)*4)
            do while(loffset.gt.offlimit)
              offset = int(offlimit)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1010
              loffset = loffset - offlimit
            end do
            offset = int(loffset)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1010
          end if

          ! Store in buffer
          buffer = real(Stokes(is,:))

          ! Write buffer
          call MPI_FILE_WRITE(funit,buffer(1),nn, &
                              MPI_REAL,MPI_STATUS_IGNORE,ierr)
          if (ierr.ne.0) goto 1300

          ! After V, skip
          if (is.eq.4) exit

          ! Final wavelength shift
          if (if1.lt.buff%nn) then
            loffset = dble((buff%nn-if1)*4)
            do while(loffset.gt.offlimit)
              offset = int(offlimit)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1010
              loffset = loffset - offlimit
            end do
            offset = int(loffset)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1010
          end if

        end do ! Stokes parameter

      !
      ! Serial
      !
      else

        ! Allocate buffer
        nn = size(stokes)
        allocate(buffer(nn))

        ! Add to buffer
        buffer = real(reshape(transpose(Stokes),(/nn/)))

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

      end if ! Serial/MPI (freqs)

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      !
      ! Tau?
      !
      if (wtau) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Tau', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1100


        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size/4) + &
                  dble(buff%head_size+buff%geom_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1110
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1110

        ! Allocate buffer
        nn = nf

        ! Store in buffer
        buffer(1:nn) = real(tau)

        !
        ! If splitting frequencies
        !
        if (nproc.gt.1) then

          ! Initial wavelength shift
          if (if0.gt.1) then
            loffset = dble(if0-1)*4
            do while(loffset.gt.offlimit)
              offset = int(offlimit)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1110
              loffset = loffset - offlimit
            end do
            offset = int(loffset)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1110
          end if
        end if ! MPI (freqs)

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1310

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if ! Tau

      return

1000  umsg = 'Error opening Stokes file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error opening Tau file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1110  umsg = 'Error seeking Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1310  umsg = 'Error writing Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine write_CLE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the contribution function.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!              z(dfloat(:)): Heights array\n
      !!      Contr(dfloat(:,:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writectr(filename,iph,ith,omega,Geom,z,Contr,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega,z
      double precision,dimension(:,:,:), intent(in):: Contr

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,ifreq,iran,iz,ii,jj,i0,i1,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectr'

      !
      ! Convert the integers into appropriate length strings
      !
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/Contribution_'//trim(cth)// &
             '_'//trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream',action='write',form='unformatted')


        !
        ! Write content
        !

        ! Identification
        write(200,err=1100) 'bc'

        ! Number of frequencies, number of heights, angles of LOS
        write(200,err=1100) nfreq,nZ,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Height axis
        write(200,err=1100) z

        ! Limited heights
        if (nz.ne.Rnz) then

          ! For each Stokes parameter
          do ios=1,4

            ! For each height less than lower limit, write zero
            do iz=1,Rz0-1
              ! For each frequency
              do ifreq=1,nfreq
                write(200,err=1100) 0d0
              end do ! Frequency
            end do ! Height

            ! Write contribution function
            write(200,err=1100) Contr(ios,:,:)

            ! For each height larger than upper limit, write zero
            do iz=Rz1+1,nz
              ! For each frequency
              do ifreq=1,nfreq
                write(200,err=1100) 0d0
              end do ! Frequency
            end do ! Height

          end do ! For each Stokes parameter

        ! All heights
        else

          ! For each Stokes parameter
          do ios=1,4

            ! Write contribution function, order: is, iz, ifreq
            write(200,err=1100) Contr(ios,:,:)

          end do ! Stokes parameters

        end if ! Limited heights

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Contribution_'// &
                           trim(cth)//'_'//trim(cph), &
                           MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each Stokes parameter
          do jj=1,4

            ! For each height
            do iz=1,nz

              ! Out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! For each entry to write
                do iran=1,buff%nran

                  ! Atom and transition
                  i0 = buff%indx(1,iran)
                  i1 = buff%indx(2,iran)
                  nn = i1-i0+1

                  ! Fill buffer
                  buffer(ii+1:ii+nn) = real(0)

                  ! Advance index
                  ii = ii + nn

                end do ! Ranges

              ! In bounds
              else

                ! For each entry to write
                do iran=1,buff%nran

                  ! Atom and transition
                  i0 = buff%indx(1,iran)
                  i1 = buff%indx(2,iran)
                  nn = i1-i0+1

                  ! Fill buffer
                  buffer(ii+1:ii+nn) = &
                      real(reshape(Contr(jj,i0:i1,iz), (/ nn /)))

                  ! Advance index
                  ii = ii + nn

                end do ! Ranges

              end if ! Height bounds

            end do ! Heights
          end do ! Stokes parameters

        ! All
        else

          ! Not full height range
          if (nz.ne.Rnz) then

            ! For each Stokes parameter
            do jj=1,4

              ! For each height
              do iz=1,nz

                ! Out of bounds
                if (iz.lt.Rz0.or.iz.gt.Rz1) then

                  ! Fill buffer
                  buffer(ii+1:ii+nfreq) = real(0)

                ! In bounds
                else

                  ! Fill buffer
                  buffer(ii+1:ii+nfreq) = &
                      real(reshape(Contr(jj,:,iz), (/ nfreq /)))

                end if ! Height bounds

                ! Advance index
                ii = ii + nfreq

              end do ! Heights
            end do ! Stokes parameters

          ! Full height range
          else

            ! Size
            nn = nfreq*nz

            ! For each Stokes parameter
            do jj=1,4

              ! Fill buffer
              buffer(ii+1:ii+nn) = real(reshape(Contr(jj,:,:), &
                                        (/ nn /)))
              ii = ii + nn
            end do

          end if ! Full height range
        end if ! Full frequency range

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writectr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the contribution function.\n
      !!    filename(character(:)): Name of the file to write\n
      !!        Contr(real(:,:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writectr_inv(filename,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      real, dimension(:,:,:), intent(in):: Contr

      ! Local

      integer:: ierr,iz,ii,jj,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectr_inv'

      ! Size
      nn = buff%buffer_size/4

      ! Allocate buffer
      allocate(buffer(nn))

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(filename)//'/Contribution', &
                         MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Initialize offset buffer
      ii = 0

      ! Not full height range
      if (nz.ne.Rnz) then

        ! For each Stokes parameter
        do jj=1,4

          ! For each height
          do iz=1,nz

            ! Out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! Fill buffer
              buffer(ii+1:ii+buff%nn) = real(0)

            ! In bounds
            else

              ! Fill buffer
              buffer(ii+1:ii+buff%nn) = &
                                reshape(Contr(jj,:,iz), (/ buff%nn /))

            end if ! Height bounds

            ! Advance index
            ii = ii + buff%nn

          end do ! Heights
        end do ! Stokes parameters

      ! Full height range
      else

        ! For each Stokes parameter
        do jj=1,4

          ! Fill buffer
          buffer(ii+1:ii+buff%nn*nz) = &
                              reshape(Contr(jj,:,:), (/ buff%nn*nz /))

          ! Advance index
          ii = ii + buff%nn*nz

        end do ! Stokes parameters

      end if ! Full height range

      ! Write buffer
      call MPI_FILE_WRITE(funit,buffer(1),nn, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writectr_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the contribution function to the
      !! intensity.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!              z(dfloat(:)): Heights array\n
      !!        Contr(dfloat(:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writectrI(filename,iph,ith,omega,Geom,z,Contr,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega,z
      double precision,dimension(:,:), intent(in):: Contr

      ! Local

      character(len=4):: cph,cth

      integer:: ios,iz,ifreq,ierr,iran,ii,i0,i1,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectrI'


      !
      ! Convert the integers into appropriate length strings
      !
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/Contribution_'//trim(cth)// &
             '_'//trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream',action='write',form='unformatted')

        !
        ! Write content
        !

        ! Identification
        write(200,err=1100) 'bc'

        ! Number of frequencies, number of heights, angles of LOS
        write(200,err=1100) nfreq,nZ,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Height axis
        write(200,err=1100) z

        ! Write zeros before lower limit
        do iz=1,Rz0-1
          ! For each frequency
          do ifreq=1,nfreq
            write(200,err=1100) 0d0
          end do ! Frequencies
        end do ! Height below lower limit

        ! Write contribution function, order: is, iz, ifreq
        write(200,err=1100) Contr


        ! Write zeros above upper limit
        do iz=Rz1+1,nZ
          ! For each frequency
          do ifreq=1,nfreq
            write(200,err=1100) 0d0
          end do ! Frequencies
        end do ! Height above upper limit

        ! Fill the other Stokes parameters with 0
        do ios=1,3

          ! For each height
          do iz=1,nZ

            ! For each frequency
            do ifreq=1,nfreq

              write(200,err=1100) 0d0

            end do ! Frequencies
          end do ! Heights
        end do ! Stokes parameters

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))
        buffer = 0e0

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Contribution_'// &
                           trim(cth)//'_'//trim(cph), &
                           MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! Heights
          do iz=1,nz

            ! Out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! For each entry to write
              do iran=1,buff%nran

                ! Atom and transition
                i0 = buff%indx(1,iran)
                i1 = buff%indx(2,iran)
                nn = i1-i0+1

                ! Fill buffer
                buffer(ii+1:ii+nn) = real(0)

                ! Advance index
                ii = ii + nn

              end do ! Ranges

            ! In bounds
            else

              ! For each entry to write
              do iran=1,buff%nran

                ! Atom and transition
                i0 = buff%indx(1,iran)
                i1 = buff%indx(2,iran)
                nn = i1-i0+1

                ! Fill buffer
                buffer(ii+1:ii+nn) = &
                    real(reshape(Contr(i0:i1,iz), (/ nn /)))

                ! Advance index
                ii = ii + nn

              end do ! Ranges

            end if ! Height bounds

          end do ! Heights

        ! All
        else

          ! Limited range
          if (nz.ne.Rnz) then

            ! For each height
            do iz=1,nz

              ! Out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! Fill buffer
                buffer(ii+1:ii+nfreq) = real(0)

              ! In bounds
              else

                ! Fill buffer
                buffer(ii+1:ii+nfreq) = &
                    real(reshape(Contr(:,iz), (/ nfreq /)))

              end if ! Height bounds

              ! Advance index
              ii = ii + nfreq

            end do ! Heights

          ! Full range
          else

            ! Size
            nn = nfreq*nz

            ! Fill buffer
            buffer(ii+1:ii+nn) = real(reshape(Contr, (/ nn /)))

          end if ! Height range
        end if ! Limited ranges

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writectrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the contribution function to the
      !! intensity.\n
      !!    filename(character(:)): Name of the file to write\n
      !!          Contr(real(:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writectrI_inv(filename,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      real, dimension(:,:), intent(in):: Contr

      ! Local

      integer:: iz,ierr,ii,nn,nt

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectrI_inv'


      ! Size
      nt = buff%buffer_size/4
      nn = nt/4

      ! Allocate buffer
      allocate(buffer(nt))
      buffer = 0e0

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(filename)//'/Contribution', &
                         MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Limited range
      if (nz.ne.Rnz) then

        ! Initialize buffer
        ii = 0

        ! For each height
        do iz=1,nz

          ! Out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! Fill buffer
            buffer(ii+1:ii+buff%nn) = real(0)

          ! In bounds
          else

            ! Fill buffer
            buffer(ii+1:ii+buff%nn) = &
                                   reshape(Contr(:,iz), (/ buff%nn /))

          end if ! Height bounds

          ! Advance index
          ii = ii + buff%nn

        end do ! Heights

      ! Full range
      else

        ! Fill buffer
        buffer(1:buff%nn*nz) = reshape(Contr, (/ buff%nn*nz /))

      end if ! Height range

      ! Write buffer
      call MPI_FILE_WRITE(funit,buffer(1),nt, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writectrI_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes the contribution function to RAM.\n
      !!    e_Contr(dfloat(:,:,:)): Output contribution function\n
      !!      Contr(dfloat(:,:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine setctr(e_Contr,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      double precision, dimension(:,:,:), intent(in):: Contr
      real, dimension(:,:,:), intent(out):: e_Contr

      ! Local

      integer:: iran,ii,jj,i0,i1,nn


      ! Only master
      if (pid.gt.0) return

      ! Zero out of bounds
      e_Contr(:,:,1:Rz0-1) = 0.0
      e_Contr(:,:,Rz1+1:nz) = 0.0

      ! If specified
      if (buff%nran.gt.0) then

        ! For each Stokes parameter
        do jj=1,4

          ! Initialize buffer
          ii = 0

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            e_Contr(jj,ii+1:ii+nn,Rz0:Rz1) = &
                                         real(Contr(jj,i0:i1,Rz0:Rz1))

            ! Advance index
            ii = ii + nn

          end do ! Ranges
        end do ! Stokes parameters

      ! All
      else

        ! Copy
        e_Contr(:,:,Rz0:Rz1) = real(Contr(:,:,Rz0:Rz1))

      end if ! Full frequency range

      return

      end subroutine setctr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes the intensity contribution function to RAM.\n
      !!    e_Contr(dfloat(:,:,:)): Output contribution function\n
      !!      Contr(dfloat(:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine setctrI(e_Contr,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      double precision, dimension(:,:), intent(in):: Contr
      real, dimension(:,:,:), intent(out):: e_Contr

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Only master
      if (pid.gt.0) return

      ! Zero out of bounds
      e_Contr(:,:,1:Rz0-1) = 0.0
      e_Contr(:,:,Rz1+1:nz) = 0.0

      ! If specified
      if (buff%nran.gt.0) then

        ! Initialize buffer
        ii = 0

        ! For each entry to write
        do iran=1,buff%nran

          ! Atom and transition
          i0 = buff%indx(1,iran)
          i1 = buff%indx(2,iran)
          nn = i1-i0+1

          ! Fill buffer
          e_Contr(1,ii+1:ii+nn,Rz0:Rz1) = real(Contr(i0:i1,Rz0:Rz1))

          ! Advance index
          ii = ii + nn

        end do ! Ranges

      ! All
      else

        ! Copy
        e_Contr(1,:,Rz0:Rz1) = real(Contr(:,Rz0:Rz1))

      end if ! Full frequency range

      return

      end subroutine setctrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the height where the optical depth is
      !! equal to one.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!            tau(dfloat(:)): Height where the optical depth is
      !!                            equal to one\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writetau(filename,iph,ith,omega,Geom,tau,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega
      double precision,dimension(:), intent(in):: tau

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,iran,ii,i0,i1,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writetau'


      !
      ! Convert the integers into appropriate length strings
      !
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/Tau_'//trim(cth)// &
             '_'//trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream', action='write', form='unformatted')

        !
        ! Write content
        !

        ! Identification
        write(200,err=1100) 'bt'

        ! Number of frequencies, angles of LOS
        write(200,err=1100) nfreq,nZ,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Height of tau=1
        write(200,err=1100) tau

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Tau_'//trim(cth)// &
                           '_'//trim(cph),MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            buffer(ii+1:ii+nn) = real(tau(i0:i1))

            ! Advance index
            ii = ii + nn

          end do ! Ranges

        ! All
        else

          ! Fill buffer
          buffer = real(tau)

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      return

1000  umsg = 'Error opening tau file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing tau file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writetau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes a file with the height where the optical depth is
      !! equal to one.\n
      !!    filename(character(:)): Name of the file to write\n
      !!              tau(real(:)): Height where the optical depth is
      !!                            equal to one\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writetau_inv(filename,tau,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      real, dimension(:), intent(in):: tau

      ! Local

      integer:: ierr,nn

      integer(kind=MPI_OFFSET_KIND):: offset

      double precision:: loffset


      ! Routine name
      urou = 'writetau_inv'

      ! Size
      nn = buff%buffer_size/4

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(filename)//'/Tau', &
                         MPI_MODE_WRONLY, &
                         MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Write buffer
      call MPI_FILE_WRITE(funit,tau(1),nn, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      return

1000  umsg = 'Error opening tau file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writetau_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes to RAM the height where the optical depth is equal to
      !! one.\n
      !!          e_tau(dfloat(:)): Output height where the optical
      !!                            depth is equal to one\n
      !!            tau(dfloat(:)): Height where the optical depth is
      !!                            equal to one\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine settau(e_tau,tau,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      double precision, dimension(:), intent(in):: tau
      real, dimension(:), intent(out):: e_tau

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Only master
      if (pid.gt.0) return

      ! Initialize buffer
      ii = 0

      ! If specified
      if (buff%nran.gt.0) then

        ! For each entry to write
        do iran=1,buff%nran

          ! Atom and transition
          i0 = buff%indx(1,iran)
          i1 = buff%indx(2,iran)
          nn = i1-i0+1

          ! Fill buffer
          e_tau(ii+1:ii+nn) = real(tau(i0:i1))

          ! Advance index
          ii = ii + nn

        end do ! Ranges

      ! All
      else

        ! Fill buffer
        e_tau = real(tau)

      end if

      return

      end subroutine settau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes into file the inelastic collisions for all active
      !! atoms.\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!     folder(character(:)): Path to the output folder\n
      !!  btt(IO_helper_class(:)): Info about what to store\n
      !!  bll(IO_helper_class(:)): Info about what to store
      subroutine writecols(Atom,folder,btt,bll)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(IO_helper_class), intent(in):: btt,bll
      character(len=500), intent(in):: folder

      ! Local

      integer:: ierr,ii,iz,it1,it
      integer:: ia,iterm,iterm1,ilevel,ilevel1,ios

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      !
      ! Slaves just wait
      !
      if (pid.gt.0) then
        call control
        return
      end if

      ! Routine name
      urou = 'writecols'

      !
      ! 1D
      !
      if (run_mode.eq.0) then

        ! Open file to write into
        open (200,file=trim(folder)//'/cols-TT', status='unknown', &
              iostat=ios, err=1000, access='stream', action='write', &
              form='unformatted')
        open (300,file=trim(folder)//'/cols-LL', status='unknown', &
              iostat=ios, err=1001, access='stream', action='write', &
              form='unformatted')

        ! Identification
        write(200,err=1100) 'ct'
        write(300,err=1100) 'cl'

        ! Write number of atoms
        write(200,err=1100) NA
        write(300,err=1101) NA

        ! Write number of heights
        write(200,err=1100) NZ
        write(300,err=1101) NZ

        ! For each atom
        do ia=1,nA

          ! Write number of terms/levels
          write(200,err=1100) Atom(ia)%nmulti
          write(300,err=1101) Atom(ia)%nlevel

          ! For each term pair
          do iterm=1,Atom(ia)%nMulti
            do iterm1=1,Atom(ia)%nMulti

              write(200,err=1100) Atom(ia)%Ccoeff(iterm1,iterm,:)

            end do
          end do

          ! For each level pair
          do ilevel=1,Atom(ia)%nlevel
            do ilevel1=1,Atom(ia)%nlevel

              write(300,err=1101) Atom(ia)%CcoeffJ(ilevel1,ilevel,:)

            end do
          end do

        end do ! Atoms

        ! Close files
        close(200)
        close(300)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        if (btt%buffer_size.gt.bll%buffer_size) then
          allocate(buffer(btt%buffer_size/4))
        else
          allocate(buffer(bll%buffer_size/4))
        end if

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/cols-TT', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(btt%buffer_size) + &
                  dble(btt%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (btt%nran.gt.0) then

          ! For each entry to write
          do ia=1,btt%nran
            ! For each height
            do iz=1,nz
              ! Advance
              ii = ii + 1
              buffer(ii) = real(Atom(btt%indx(1,ia))% &
                                Ccoeff(btt%indx(2,ia), &
                                       btt%indx(3,ia),iz))
            end do ! Height
            ! For each height
            do iz=1,nz
              ! Advance
              ii = ii + 1
              buffer(ii) = real(Atom(btt%indx(1,ia))% &
                                Ccoeff(btt%indx(3,ia), &
                                       btt%indx(2,ia),iz))
            end do ! Height
          end do ! Atom

        ! All
        else

          ! For each atom
          do ia=1,nA
            ! For each pair of terms
            do it=1,Atom(ia)%nMulti
            do it1=1,Atom(ia)%nMulti
              ! For each height
              do iz=1,nz
                ! Advance
                ii = ii + 1
                buffer(ii) = real(Atom(ia)%Ccoeff(it1,it,iz))
              end do ! Height
            end do ! Term
            end do ! Term
          end do ! Atom

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),btt%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/cols-LL', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1001

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(bll%buffer_size) + &
                  dble(bll%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1011
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1011

        ! Initialize buffer
        ii = 0

        ! If specified
        if (bll%nran.gt.0) then

          ! For each entry to write
          do ia=1,bll%nran
            ! For each height
            do iz=1,nz
              ! Advance
              ii = ii + 1
              buffer(ii) = real(Atom(bll%indx(1,ia))% &
                                CcoeffJ(bll%indx(2,ia), &
                                        bll%indx(3,ia),iz))
            end do ! Height
            ! For each height
            do iz=1,nz
              ! Advance
              ii = ii + 1
              buffer(ii) = real(Atom(bll%indx(1,ia))% &
                                CcoeffJ(bll%indx(3,ia), &
                                        bll%indx(2,ia),iz))
            end do ! Height
          end do ! Atom

        ! All
        else

          ! For each atom
          do ia=1,nA
            ! For each pair of terms
            do it=1,Atom(ia)%nlevel
            do it1=1,Atom(ia)%nlevel
              ! For each height
              do iz=1,nz
                ! Advance
                ii = ii + 1
                buffer(ii) = real(Atom(ia)%CcoeffJ(it1,it,iz))
              end do ! Height
            end do ! Term
            end do ! Term
          end do ! Atom

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),bll%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1301

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      ! Control
      call control
      return

1000  umsg = 'Error opening collisions TT file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking collisions TT file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing collisions TT file'
      close(200)
      close(300)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing collisions TT file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening collisions LL file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1011  umsg = 'Error seeking collisions LL file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing collisions LL file'
      close(200)
      close(300)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1301  umsg = 'Error writing collisions LL file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writecols

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes into file the damping parameter for all active
      !! atoms.\n
      !!       Atom(Atom_class): Structure with the atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!   folder(character(:)): Path to the output folder\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writedamp(Atom,Atmo,folder,buff)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder

      ! Local

      integer:: ierr,ii
      integer:: iran,ia,iterml,itermu,itran,i,i1,ios
      double precision, dimension(nz):: Dw, DwT

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      !
      ! Slaves just wait
      !
      if (pid.gt.0) then
        call control
        return
      end if

      ! Routine name
      urou = 'writedamp'


      !
      ! 1D
      !
      if (run_mode.eq.0) then

        ! Open file to write into
        open (200,file=trim(folder)//'/damping', status='unknown', &
              iostat=ios, err=1000, access='stream', action='write', &
              form='unformatted')

        ! Identification
        write(200,err=1100) 'da'

        ! Write number of atoms
        write(200,err=1100) NA

        ! Write number of heights
        write(200,err=1100) NZ

        ! For each atom
        do ia=1,nA

          ! Write number of terms/levels and transitions
          write(200,err=1100) Atom(ia)%ntran

          ! Thermal width
          DwT = Atom(ia)%cDopp*sqrt(Atmo%T)

          ! For each transition
          do itran=1,Atom(ia)%ntran

            ! Doppler width
            Dw = Atom(ia)%Dfreq(itran)*sqrt(DwT*DwT + &
                                            Atmo%vmi*Atmo%vmi)
            Dw = 1d0/Dw

            ! Identify involved terms
            itermu = -1
            do i=1,Atom(ia)%nMulti-1
              do i1=i+1,Atom(ia)%nMulti
                if (Atom(ia)%irad(i,i1).eq.itran) then
                  iterml = i
                  itermu = i1
                  exit
                end if
              end do
              if (itermu.ge.0) exit
            end do

            write(200,err=1100) (Atom(ia)%ldamp(itran,:) + &
                                 Atom(ia)%damp(itermu,:) + &
                                 Atom(ia)%damp(iterml,:))*Dw

          end do

        end do ! Atoms

        ! Close files
        close(200)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/damping', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            ia = buff%indx(1,iran)
            itran = buff%indx(2,iran)

            ! Thermal width
            DwT = Atom(ia)%cDopp*sqrt(Atmo%T)

            ! Doppler width
            Dw = Atom(ia)%Dfreq(itran)* &
                 sqrt(DwT*DwT + Atmo%vmi*Atmo%vmi)
            Dw = 1d0/Dw

            ! Identify involved terms
            itermu = -1
            do i=1,Atom(ia)%nMulti-1
              do i1=i+1,Atom(ia)%nMulti
                if (Atom(ia)%irad(i,i1).eq.itran) then
                  iterml = i
                  itermu = i1
                  exit
                end if
              end do
              if (itermu.ge.0) exit
            end do

            ! Advance buffer
            buffer(ii+1:ii+nz) = real((Atom(ia)%ldamp(itran,:) + &
                                       Atom(ia)%damp(itermu,:) + &
                                       Atom(ia)%damp(iterml,:))*Dw)
            ii = ii + nz

          end do ! Ranges

        ! All
        else

          ! For each atom
          do ia=1,nA

            ! Thermal width
            DwT = Atom(ia)%cDopp*sqrt(Atmo%T)

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Doppler width
              Dw = Atom(ia)%Dfreq(itran)*sqrt(DwT*DwT + &
                                              Atmo%vmi*Atmo%vmi)
              Dw = 1d0/Dw

              ! Identify involved terms
              itermu = -1
              do i=1,Atom(ia)%nMulti-1
                do i1=i+1,Atom(ia)%nMulti
                  if (Atom(ia)%irad(i,i1).eq.itran) then
                    iterml = i
                    itermu = i1
                    exit
                  end if
                end do
                if (itermu.ge.0) exit
              end do

              ! Advance buffer
              buffer(ii+1:ii+nz) = real((Atom(ia)%ldamp(itran,:) + &
                                         Atom(ia)%damp(itermu,:) + &
                                         Atom(ia)%damp(iterml,:))*Dw)
              ii = ii + nz

            end do ! Transitions
          end do ! Atoms

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      ! Control
      call control

      return

1000  umsg = 'Error opening damping file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking damping file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing damping file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing damping file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writedamp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes into file the background continuum quantities.\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!      folder(character(:)): Path to the output folder\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writeback(Cont,omega,folder,MPID,buff)

      ! I/O

      type(Continuum_class):: Cont
      type(MPI_class):: MPID
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder
      double precision, dimension(:), intent(in):: omega

      ! Local

      integer:: ierr,i0,i1,nn,ii,jj,iran
      integer:: iz,nfl,ndl,nxdir,idir
      integer:: iproc,iproc0,bsize,psize,ios
      integer, dimension(:), allocatable:: ndir

      integer(kind=MPI_OFFSET_KIND):: offset

      double precision:: loffset

      double precision, dimension(:), allocatable:: dbuffer

      real, dimension(:), allocatable:: buffer

      ! Routine name
      urou = 'writeback'

      ! Master opens file and writes dimensions
      if (pid.eq.0) then

        !
        ! 1D
        !
        if (run_mode.eq.0) then

          open (200,file=trim(folder)//'/background', &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='write', form='unformatted')

          ! Identification
          write(200,err=1100) 'ba'

          ! Write number of frequencies
          write(200,err=1100) nfreq

          ! Write frequencies
          write(200,err=1100) omega

          ! Write number of heights
          write(200,err=1100) NZ

        end if
      end if

      ! If MPI
      if (MPID%mpi) then

        ! Control
        call control

        ! Allocate number of directions for every cpu
        allocate(ndir(0:nproc-1))
        if (pid.eq.0) Cont%ndir = 0

        ! Master gather number of directions
        call MPI_GATHER(Cont%ndir, 1, MPI_INTEGER, ndir(0), 1, &
                        MPI_INTEGER, 0, MPI_COMM_RT, ierr)

        ! MASTER
        if (pid.eq.0) then

          ! Look for the maximum number of directions
          nxdir = maxval(ndir)

          ! Last processor checked
          iproc0 = 1

          ! Maximum dimension to receive
          bsize = MPID%nxfreq*nz*3*nxdir

          ! Allocate buffer
          allocate(dbuffer(bsize))

          ! Allocate background
          allocate(Cont%c(nfreq,3,nxdir,nz))
          Cont%c = 0d0

          ! For each process
          do iproc=1,nproc-1

            nfl = MPID%nf(iproc)
            ndl = ndir(iproc)
            psize = nz*nfl*ndl*3

            ! Receive buffer
            do while (.True.)
              call MPI_RECV(dbuffer, psize, MPI_DOUBLE_PRECISION, &
                            iproc, iproc, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! If the processor had the maximum number of
            ! directions
            if (ndl.eq.nxdir) then

              Cont%c(MPID%if0(iproc):MPID%if1(iproc),:,:,:) = &
                       reshape(dbuffer(1:psize), (/ nfl,3,ndl,nz /))

            ! If the buffer is smaller than the maximum
            else

              ! We have to replicate the directions, so for each
              ! height
              do iz=1,nz

                ! For each direction
                do idir=1,nxdir

                  Cont%c(MPID%if0(iproc):MPID%if1(iproc),:,idir,iz)= &
                        reshape(dbuffer((iz-1)*(nfl*3)+1:iz*nfl*3), &
                                                        (/ nfl,3 /))

                end do ! directions
              end do ! heights

            end if ! Maximum number of directions

          end do ! Processors

        ! SLAVE
        else

          ! Compute size of data to send
          psize = nz*MPID%nf(pid)*3*Cont%ndir

          ! Send data
          do while (.True.)
            call MPI_SEND(Cont%c(MPID%if0(pid),1,1,1), &
                          psize, MPI_DOUBLE_PRECISION, 0, pid, &
                          MPI_COMM_RT, ierr)
            if (ierr.eq.0) exit
          end do

        end if ! Master or slave

      ! If serial
      else

        ! Save nxdir
        nxdir = Cont%ndir

      end if ! MPI

      ! If Master
      if (pid.eq.0) then

        !
        ! 1D
        !
        if (run_mode.eq.0) then

          ! Write number of directions
          write(200,err=1100) nxdir

          ! Write the full block
          write(200,err=1100) Cont%c

          ! If MPI, deallocate background
          if (MPID%mpi) deallocate(Cont%c)

          ! Close files
          close(200)

        !
        ! 1.5D
        !
        else if (run_mode.eq.1) then

          ! Allocate buffer
          allocate(buffer(buff%buffer_size/4))

          ! Open file
          call MPI_FILE_OPEN(MPI_COMM_SELF, &
                             trim(folder)//'/background', &
                             MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                             ierr)
          if (ierr.ne.0) goto 1000

          !
          ! Column offset
          !

          ! Get offset
          loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                    dble(buff%head_size)
          do while(loffset.gt.offlimit)
            offset = int(offlimit)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1010
            loffset = loffset - offlimit
          end do
          offset = int(loffset)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010

          ! Initialize buffer
          ii = 0

          ! If specified
          if (buff%nran.gt.0) then

            ! For each height
            do iz=1,nz

              ! For each variable
              do jj=1,3

                ! For each entry to write
                do iran=1,buff%nran

                  ! Atom and transition
                  i0 = buff%indx(1,iran)
                  i1 = buff%indx(2,iran)
                  nn = i1-i0+1

                  ! Fill buffer
                  buffer(ii+1:ii+nn) = real(Cont%c(i0:i1,jj,1,iz))

                  ! Advance index
                  ii = ii + nn

                end do ! Ranges
              end do ! Variables
            end do ! Heights

          ! All
          else

            ! Size
            nn = nfreq*nz*3

            ! Fill buffer
            buffer = real(reshape(Cont%c(:,:,1,:), (/ nn /)))

          end if

          ! Write buffer
          call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                              MPI_REAL,MPI_STATUS_IGNORE,ierr)
          if (ierr.ne.0) goto 1300

          ! Close file
          call MPI_FILE_CLOSE(funit, ierr)

        end if ! Run mode
      end if ! Master

      ! Control
      call control

      return

1000  umsg = 'Error opening background file'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking background file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing background file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing background file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writeback

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes into file the atmospheric data.\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Bfield(Bfield_class): Structure with magnetic field data\n
      !!    folder(character(:)): Path to the output folder\n
      !!  buff(IO_helper_class(:)): Info about what to store
      subroutine writeatmo(Atmo,Bfield,folder,buff)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder

      ! Local

      integer:: iz, ios
      integer:: ierr,ii

      integer(kind=MPI_OFFSET_KIND):: offset

      double precision, dimension(:), allocatable:: buffer

      double precision:: loffset

      ! Routine name
      urou = 'writeatmo'

      ! If not master, return
      if (pid.gt.0) then
        call control
        return
      end if

      !
      ! 1D or inversion
      !
      if (run_mode.le.0) then

        ! Master opens file and writes dimensions
        open (200,file=trim(folder)//'/atmos.dat', &
              status='unknown', iostat=ios, err=1000, action='write')

        ! If tau scale
        if (ztau) then

          ! Write header
          write(200,'(A)',err=1100) '! Atmospheric data'
          write(200,'(A)',err=1100) '!             Height [km]'// &
                                     '                tau cont'//&
                                     '        Chi cont [cm^-1]'//&
                                     '         Temperature [K]'//&
                                     ' Gas pressure [dyn/cm^2]'//&
                                     '        Density [g/cm^3]'//&
                                     '   Magnetic strength [G]'//&
                                     '    Magnetic inclination'//&
                                     '        Magnetic azimuth'//&
                                     '       Velocity x [km/s]'//&
                                     '       Velocity y [km/s]'//&
                                     '       Velocity z [km/s]'//&
                                     '  Microturbulence [km/s]'//&
                                     ' e^- pressure [dyn/cm^2]'//&
                                     '     e^- density [cm^-3]'//&
                                     '       H density [cm^-3]'//&
                                     '  atom H density [cm^-3]'//&
                                     '      H- density [cm^-3]'//&
                                     '    HI_0 density [cm^-3]'//&
                                     '    HI_1 density [cm^-3]'//&
                                     '    HI_2 density [cm^-3]'//&
                                     '    HI_3 density [cm^-3]'//&
                                     '    HI_4 density [cm^-3]'//&
                                     '     HII density [cm^-3]'

          ! For each height, output
          do iz=1,nz
            write(200,'(1x,24(2x,es22.15))') &
                                         Atmo%zalt(iz), &
                                         Atmo%z(iz), &
                                         Atmo%chi500(iz), &
                                         Atmo%T(iz), &
                                         Atmo%Pg(iz), &
                                         Atmo%rho(iz), &
                                         Bfield%Bstrength(iz), &
                                         Bfield%Btheta(iz), &
                                         Bfield%Bphi(iz), &
                                         Atmo%vx(iz)*1d6*c, &
                                         Atmo%vy(iz)*1d6*c, &
                                         Atmo%vz(iz)*1d6*c, &
                                         Atmo%vmi(iz)*1d6*c, &
                                         Atmo%Pe(iz), &
                                         Atmo%ne(iz), &
                                         Atmo%nHT(iz), &
                                         Atmo%nHa(iz), &
                                         Atmo%nHm(iz), &
                                         Atmo%nh(iz,1), &
                                         Atmo%nh(iz,2), &
                                         Atmo%nh(iz,3), &
                                         Atmo%nh(iz,4), &
                                         Atmo%nh(iz,5), &
                                         Atmo%nh(iz,6)
          end do

        ! If height scale
        else

          ! Write header
          write(200,'(A)',err=1100) '! Atmospheric data'
          write(200,'(A)',err=1100) '!             Height [km]'// &
                                     '                tau cont'//&
                                     '        Chi cont [cm^-1]'//&
                                     '         Temperature [K]'//&
                                     ' Gas pressure [dyn/cm^2]'//&
                                     '        Density [g/cm^3]'//&
                                     '   Magnetic strength [G]'//&
                                     '    Magnetic inclination'//&
                                     '        Magnetic azimuth'//&
                                     '       Velocity x [km/s]'//&
                                     '       Velocity y [km/s]'//&
                                     '       Velocity z [km/s]'//&
                                     '  Microturbulence [km/s]'//&
                                     ' e^- pressure [dyn/cm^2]'//&
                                     '     e^- density [cm^-3]'//&
                                     '       H density [cm^-3]'//&
                                     '  atom H density [cm^-3]'//&
                                     '      H- density [cm^-3]'//&
                                     '    HI_0 density [cm^-3]'//&
                                     '    HI_1 density [cm^-3]'//&
                                     '    HI_2 density [cm^-3]'//&
                                     '    HI_3 density [cm^-3]'//&
                                     '    HI_4 density [cm^-3]'//&
                                     '     HII density [cm^-3]'

          ! For each height, output
          do iz=1,nz
            write(200,'(1x,24(2x,es22.15))') &
                                         Atmo%z(iz)*1d-5, &
                                         Atmo%zalt(iz), &
                                         Atmo%chi500(iz), &
                                         Atmo%T(iz), &
                                         Atmo%Pg(iz), &
                                         Atmo%rho(iz), &
                                         Bfield%Bstrength(iz), &
                                         Bfield%Btheta(iz), &
                                         Bfield%Bphi(iz), &
                                         Atmo%vx(iz)*1d6*c, &
                                         Atmo%vy(iz)*1d6*c, &
                                         Atmo%vz(iz)*1d6*c, &
                                         Atmo%vmi(iz)*1d6*c, &
                                         Atmo%Pe(iz), &
                                         Atmo%ne(iz), &
                                         Atmo%nHT(iz), &
                                         Atmo%nHa(iz), &
                                         Atmo%nHm(iz), &
                                         Atmo%nh(iz,1), &
                                         Atmo%nh(iz,2), &
                                         Atmo%nh(iz,3), &
                                         Atmo%nh(iz,4), &
                                         Atmo%nh(iz,5), &
                                         Atmo%nh(iz,6)
          end do
        end if

        ! Close files
        close(200)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/atmo.hrt', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize
        ii = 0

        ! If tau scale
        if (ztau) then

          buffer(ii+1:ii+nz) = Atmo%zalt
          ii = ii + nz
          buffer(ii+1:ii+nz) = Atmo%z
          ii = ii + nz

        ! Height scale
        else

          buffer(ii+1:ii+nz) = Atmo%z*1d-5
          ii = ii + nz
          buffer(ii+1:ii+nz) = Atmo%zalt
          ii = ii + nz

        end if

        ! Chi500
        buffer(ii+1:ii+nz) = Atmo%chi500
        ii = ii + nz
        ! T
        buffer(ii+1:ii+nz) = Atmo%T
        ii = ii + nz
        ! Pg
        buffer(ii+1:ii+nz) = Atmo%Pg
        ii = ii + nz
        ! rho
        buffer(ii+1:ii+nz) = Atmo%rho
        ii = ii + nz
        ! Bx
        buffer(ii+1:ii+nz) = Bfield%bstrength* &
                             sin(Bfield%Btheta)* &
                             cos(Bfield%Bphi)
        ii = ii + nz
        ! By
        buffer(ii+1:ii+nz) = Bfield%bstrength* &
                             sin(Bfield%Btheta)* &
                             sin(Bfield%Bphi)
        ii = ii + nz
        ! Bz
        buffer(ii+1:ii+nz) = Bfield%bstrength* &
                             cos(Bfield%Btheta)
        ii = ii + nz
        ! vx
        buffer(ii+1:ii+nz) = Atmo%vx*1d6*c
        ii = ii + nz
        ! vy
        buffer(ii+1:ii+nz) = Atmo%vy*1d6*c
        ii = ii + nz
        ! vz
        buffer(ii+1:ii+nz) = Atmo%vz*1d6*c
        ii = ii + nz
        ! vmi
        buffer(ii+1:ii+nz) = Atmo%vmi*1d6*c
        ii = ii + nz
        ! Pe
        buffer(ii+1:ii+nz) = Atmo%Pe
        ii = ii + nz
        ! ne
        buffer(ii+1:ii+nz) = Atmo%ne
        ii = ii + nz
        ! nHT
        buffer(ii+1:ii+nz) = Atmo%nHT
        ii = ii + nz
        ! nHa
        buffer(ii+1:ii+nz) = Atmo%nHa
        ii = ii + nz
        ! nHm
        buffer(ii+1:ii+nz) = Atmo%nHm
        ii = ii + nz
        ! nH
        buffer(ii+1:ii+nz*6) = reshape(Atmo%nH, (/ nz*6 /))

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if

      ! Control
      call control

      return

1000  umsg = 'Error opening atmospheric file to write'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking atmospheric file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing atmospheric file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing atmospheric file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
      return

      end subroutine writeatmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Writes atmosphere in readable model format\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!       outtypo(integer): What variable to write in the
      !!                         output\n
      !!   folder(character(:)): Path to the output folder\n
      !! filename(character(:)): Path to the original atmospheric
      !!                         file
      subroutine wAtmo(Atmo,outtypo,folder,filename)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      character(len=500), intent(in):: folder,filename
      integer, intent(in):: outtypo

      ! Local

      integer:: iz, ios
      integer:: typo

      double precision:: zfac, vfac

      ! Routine name
      urou = 'wAtmo'

      ! If not master, return
      if (pid.gt.0) then
        call control
        return
      end if

      ! Decide which type of atmosphere we are writing
      if (outtypo.gt.0) then
        typo = outtypo - 1
      else
        typo = Atmo%typo
      end if

      ! Master opens file and writes dimensions
      open (200,file=trim(folder)//'/atmo.atmos', &
            status='unknown', iostat=ios, err=1000, action='write')

      ! Header
      write(200,'(A)',err=1100) '* Written by HANLERT as an '// &
                                'update of '//trim(filename)
      write(200,'(A)',err=1100) '*'
      write(200,'(A)',err=1100) '  HANLERT-model'

      ! Scale
      if (ztau) then
        zfac = 1d0
        write(200,'(A,1x,es22.16)',err=1100) &
          'TAU SCALE',1d2/Atmo%tfreq
      else
        zfac = 1d-5
        write(200,'(A,1x,es22.16)',err=1100) &
          'HEIGHT SCALE',1d2/Atmo%tfreq
      end if

      ! Log g
      write(200,'(A)',err=1100) '*'
      write(200,'(A)',err=1100) '* LG G'
      write(200,'(2x,es22.16)',err=1100) Atmo%logg

      ! Height nodes
      write(200,'(A)',err=1100) '*'
      write(200,'(A)',err=1100) '* NDEP'
      write(200,'(i6)',err=1100) Atmo%nz

      ! Velocity factor
      vfac = c*1d6

      ! First block
      write(200,'(A)',err=1100) '*'

      ! If Standard
      if (typo.eq.0) then

        ! Tau
        if (ztau) then
          write(200,'(A)',err=1100) '*    Optical depth         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else
          write(200,'(A)',err=1100) '*      HEIGHT (km)         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        end if

        ! For every height
        do iz=1,nZ
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%ne(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac
        end do

        ! Hydrogen populations
        write(200,'(A)',err=1100) '*'
        write(200,'(A)',err=1100) '* HYDROGEN POPULATIONS (cm^-3)'
        write(200,'(A)',err=1100) '*     NH(1)      '// &
                                  '      NH(2)      '// &
                                  '      NH(3)      '// &
                                  '      NH(4)      '// &
                                  '      NH(5)      '// &
                                  '       NP        '

        ! For every height
        do iz=1,nZ
          write(200,'(6(1x,es16.10))',err=1100) Atmo%nh(iz,:)
        end do

        ! Helium populations
        if (Atmo%nhe(1,1).gt.-0.1) then

          write(200,'(A)',err=1100) '*'
          write(200,'(A)',err=1100) '* HELIUM POPULATIONS (cm^-3)'

          ! For every height
          do iz=1,nZ
            write(200,'(4(1x,es16.10))',err=1100) Atmo%nhe(iz,:)
          end do

        end if ! Helium populations present


      ! If electron number density
      else if (typo.eq.1) then

        ! Tau
        if (ztau) then
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        end if

        ! For every height
        do iz=1,nZ
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%ne(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac
        end do

        ! Write tag
        write(200,'(A)',err=1100) 'ne'


      ! If electron pressure
      else if (typo.eq.2) then

        ! Tau
        if (ztau) then
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)     Pe (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)     Pe (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        end if

        ! For every height
        do iz=1,nZ
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%Pe(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac
        end do

        write(200,'(A)',err=1100) 'pe'


      ! If electron mass density
      else if (typo.eq.3) then

        ! Tau
        if (ztau) then
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)  dens_e (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)  dens_e (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        end if

        ! For every height
        do iz=1,nZ
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%ne(iz)*me*1d3,Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac
        end do

        write(200,'(A)',err=1100) 'rhoe'


      ! If gas pressure
      else if (typo.eq.4) then

        ! Tau
        if (ztau) then
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)     Pg (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)     Pg (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        end if

        ! For every height
        do iz=1,nZ
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%Pg(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac
        end do

        write(200,'(A)',err=1100) 'pg'


      ! If gas mass density
      else if (typo.eq.5) then

        ! Tau
        if (ztau) then
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)    dens (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)    dens (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        end if

        ! For every height
        do iz=1,nZ
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%rho(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac
        end do

        write(200,'(A)',err=1100) 'rho'

      end if ! Type of output

      ! Close files
      close(200)

      ! Control
      call control

      return

1000  umsg = 'Error opening atmospheric model file to write'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1100  umsg = 'Error writing atmospheric model file'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine wAtmo

!#####################################################################
!#####################################################################
!#####################################################################

      end module iosolution_mod
