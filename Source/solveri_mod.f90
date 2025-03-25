      !> Solution of NLTE problem of the first kind
      module solveri_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     25/03/2025 V4.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     25/03/2025:    V4.0.5 - If there are velocities, it is necessary
!                             to integrate in directions even if there
!                             is no PRD (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    The dimensions are ready for blends, but they are not taken
!  into account.
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
!  solveI
!    Self-consistent solution for the RT problem for intensity
!
!  solveI_predict
!    Calculate the approximate RAM space necessary for the routines
!  solving the self-consistent problem for intensity
!
!  solveI_init
!    Initialize arrays and pointers for the self-consistent solution
!  of the problem for intensity
!
!  solveI_manager
!    Parallel-master task for the self-consistent problem for
!  intensity
!
!  solveI_RT
!    Solve the RTE for intensity for the angular quadrature
!
!  solveI_SEE
!    Advance the populations by solving the SEE and applying NG
!  acceleration if requested
!
!  emergentI
!    Formal solution for the RTE for intensity
!
!  emergentI_predict
!    Calculate the approximate RAM space necessary for the routines
!  solving the formal solution for the RTE for intensity
!
!  emergentI_init
!    Initialize arrays and pointers for the formal solution for the
!  RTE for intensity
!
!  emergentI_manager
!    Parallel-master task for the formal solution for the RTE for
!  intensity
!
!  emergentI_RT
!    Solve the RTE for intensity for a given LOS
!
!  solveJ
!    Self-consistent solution for the problem for continuum intensity
!
!  solveJ_init
!    Initialize arrays and pointers for the self-consistent solution
!  of the problem for continuum intensity
!
!  solveJ_manager
!    Parallel-master task for the self-consistent problem for
!  continuum intensity
!
!  solveJ_RT
!    Solve the RTE for continuum intensity for the angular quadrature
!
!  JKQgen
!    Conversion of radiation field tensors from multi-level intensity
!  to, in general multi-term, polarization
!
!  JKQgen_init
!    Initialize arrays and pointers for the conversion of radiation
!  field tensors from multi-level intensity to, in general multi-term,
!  polarization
!
!  JKQgen_manager
!    Parallel-master task for the conversion of radiation field
!  tensors from multi-level intensity to, in general multi-term,
!  polarization
!
!  JKQgen_RT
!    Calculate JKQ tensors from the solution of the intensity RT
!  problem
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use boundary_mod
      use commons_mod
      use comovingprd_mod
      use fieldb_mod
      use free_mod
      use iosolution_mod
      use jcalci_mod
      use mrc_mod
      use ng_mod
      use normalizer_mod
      use parameters_mod , only : kb, cSaha, fktoJ
      use rtcoeffi_mod
      use rtstepi_mod
      use seei_mod
      use setmpi_mod
      use types_mod

      ! Maximum buffer for NG_int
      double precision, parameter:: maxbuffer_NG = 500d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Self-consistent solution for the RT problem for intensity\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!      Rho_old(Rhoc_class(:)): Structure to store the density
      !!                              matrix of the previous
      !!                              iteration\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!              lload(logical): If a previous solution was
      !!                              read\n
      !!     Stokes(double(:,:,:,:)): Intensity\n
      !!            J00(double(:,:)): Mean intensity integrated over
      !!                              the absorption profile\n
      !!           J00S(double(:,:)): Mean intensity integrated over
      !!                              the emission profile\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence\n
      !!         J00P(double(:,:,:)): Intensity integrals in the
      !!                              photoionization rates
      subroutine solveI(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                        Geom,MPID,Input,lload,Stokes, &
                        J00,J00S,J00C,J00P)
      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Rhoc_class), dimension(:), intent(inout):: Rho_old
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      logical, intent(in):: lload
      double precision, &
             dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             intent(inout):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(inout):: J00C
      double precision, dimension(nxt,Rz0:Rz1), intent(inout):: J00
      double precision, dimension(nxt,Rz0:Rz1), intent(inout):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1), &
                        intent(inout):: J00P

      ! Local

      type(MRC_class):: MRC

      character(LEN=20):: iterS

      logical:: doNG,goout,gooutprd,force_ALI,PRDl,ADD,ADT,RIRAM
      logical:: lp_exu,lALI,NGI

      integer:: ierr,iter,iterr,ia,if0,if1,npz,ntpz,nsend
      integer:: NG_dim,NG_entry
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset
      double precision, dimension(nxb,nxt,Rz0:Rz1):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1):: LambdaP
      double precision, dimension(:), allocatable, target:: LO
      double precision, dimension(:,:), allocatable:: NG_scratch
      double precision, dimension(:,:), allocatable:: J00C_n

      ! Receivers

      double precision, dimension(:), allocatable, target:: Stokes_r
      double precision, dimension(:), allocatable, target:: LambdaL_r
      double precision, dimension(:), allocatable, target:: LambdaP_r
      double precision, dimension(:), allocatable, target:: Prof_r

      ! Senders

      double precision, dimension(:,:,:,:), allocatable:: Stokes_s
      double precision, dimension(:,:,:,:,:), allocatable:: Prof_s
      double precision, dimension(:,:,:,:), allocatable:: rLine_s
      double precision, dimension(:,:,:,:), allocatable:: rPhot_s

      ! Pointers

      double precision, dimension(:,:), pointer:: data1M,data1O
      double precision, dimension(:,:), pointer:: data2O
      double precision, dimension(:), pointer:: rLineO
      double precision, dimension(:), pointer:: rPhotO
      double precision, dimension(:), pointer:: p_exu


#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,-2)
#endif
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P,Input%folder,-2)
#endif

      ! Nullify pointers
      nullify(data1M,data1O,data2O,rLineO,rPhotO,p_exu)

      ! Initialize variables
      call solveI_init(Atom,Frec,Geom,MPID,Input,lload, &
                       NG_dim,NG_entry,NG_scratch,NGI,doNG,goout, &
                       force_ALI,PRDl,ADD,ADT,RIRAM,lp_exu, &
                       if0,if1,npz,ntpz,nsend, &
                       Stokes_r,LambdaL_r,LambdaP_r,Prof_r, &
                       Stokes_s,rLine_s,rPhot_s,Prof_s, &
                       data1M,data1O,data2O,rLineO,rPhotO, &
                       LO,p_exu,J00C_n)

      ! Initialize PRD
      if (PRDl) &
        call initialize_emissI(Atom,Atmo,Geom,Red)

      ! Control
      call control
      if (laborted) goto 2000

      ! If measuring performance
      if (MPID%mpi.and.Input%mpi_perf.and.pid.eq.0) &
        call report_mpi_timeI(Input%folder,Input%ID, &
                             0,0,0,.False.)

      !
      ! Start iterations
      !

      ! For each iteration between the limits specified
      do iter=Input%iteri_min,Input%iteri_max

        ! Flags for physics in Stokes
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics in populations
        if (iter.le.Input%allownphys_pop) then
          if (.not.nphysR) nphysR = .True.
        else
          if (nphysR) nphysR = .False.
        end if

        ! For each atom
        do ia=1,nA

          ! Copy current density matrix
          Rho_old(ia)%crho = Atom(ia)%crho(:,Rz0:Rz1)

        end do ! Atoms

        ! Flag for ALI
        lALI = iter.gt.Input%ALI_delay.or.force_ALI

        ! Internal PRD iterations
        do iterr=1,Input%iteri_prd

          ! Compute second order emissivity
          if (PRD) &
            call comoving_emissI2ord(Atom,Atmo,Geom,Frec,Red, &
                                     Stokes,J00,J00C,0,0, &
                                     Input%PRD_int_mode,.False.)

          ! Master
          if (pid.eq.0) then

            ! Reset new
            J00 = 0d0
            J00P = 0d0
            J00C_n = 0d0
            if (lALI) then
              LambdaL = 0d0
              LambdaP = 0d0
            end if

          end if ! Master

          ! If MPI and master
          if (MPID%mpi.and.pid.eq.0) then

            ! Call manager
            call solveI_manager(Atom,Atmo,Frec,Geom,MPID,Input, &
                                lALI,nsend,iter,iterr,npz,ntpz, &
                                Stokes_r,LambdaL_r, &
                                LambdaP_r,Prof_r,LambdaL,LambdaP, &
                                Stokes,J00C_n,J00,J00S,J00P)

          ! Serial or slave
          else

            ! Solve RTE
            call solveI_RT(Atom,LTElines,Atmo,Cont,Frec,Red,Geom, &
                           MPID,lALI,Input%ALI_photo,lp_exu,if0,if1, &
                           Stokes_s,rLine_s,rPhot_s,Prof_s, &
                           data1M,data1O,data2O,rLineO,rPhotO,LO, &
                           LambdaL,LambdaP,p_exu,Stokes,J00,J00S, &
                           J00C,J00P,J00C_n)

          end if ! Manage or compute

          !
          ! Master
          !
          if (pid.eq.0) then

            ! Calculate MRC for J if PRD
            if (PRD.and.Input%iteri_prd.gt.1) then

              ! Call the routine
              call MRCJ_sb(J00C_n,J00C,MRC)

              ! Convert cm into km
              MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

              ! Global Máster talks
              if (gpid.eq.0) then

                ! If first PRD iteration and second actual
                ! iteration
                if (iterr.eq.1.and.iter.eq.2) then
                  umsg = '         PRD            MRC(J^0_0)'// &
                         ' Freq_index  Wavelength '// &
                         'Height_index Height(km)'
                  call verbose
                end if

                ! Write in stdout
                write(umsg,'(2x,"PRD it:",1x,i3,2x,es20.12,'// &
                           '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                           iterr,MRC%values(2,1),MRC%indexes(1,1), &
                           1d2/Frec%omega(MRC%indexes(1,1)), &
                           MRC%indexes(2,1),MRC%values(1,1)
                call verbose

              end if ! Globalmaster

              ! If pass exit criteria
              if (MRC%values(2,1).le.Input%mrci_r.or. &
                  iterr.eq.Input%iteri_prd) then

                ! Flag out of PRD iteration
                gooutprd = .True.

              ! No pass
              else

                ! Keep going
                gooutprd = .False.

              end if ! Exit criteria

            ! If no PRD or internal iterations
            else

              ! Always go out
              gooutprd = .True.

            end if ! PRD and J iterations

            ! If aborting, obviously go out
            if (laborted) gooutprd = .True.

          end if ! Master

          ! If MPI
          if (MPID%mpi) then

            ! Master, overwrite
            if (pid.eq.0) J00C = J00C_n

            !
            ! Share if we are finished or not (in PRD)
            call MPI_BCAST(gooutprd, 1, MPI_LOGICAL, 0, &
                           MPI_COMM_RT, ierr)

            ! Control
            if (laborted) goto 2000


            !
            ! Share the radiation field information
            !

            ! Share J00
            call MPI_BCAST(J00(1,Rz0), MPID%sizei6(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Share J00S
           !if (stm) &
           !  call MPI_BCAST(J00S(1,Rz0), MPID%sizei6(0), &
           !                 MPI_DOUBLE_PRECISION, 0, &
           !                 MPI_COMM_RT, ierr)

            ! Share J00C
            call MPI_BCAST(J00C(1,Rz0), MPID%sizei7(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Share intensity if doing A-D PRD
            if (PRDl.and.ADD) &
              call MPI_BCAST(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                            MPI_DOUBLE_PRECISION, 0, &
                            MPI_COMM_RT, ierr)

            ! If we are going to do SEE
            if (gooutprd) then

              ! Share J00 for b-f transitions
              call MPI_BCAST(J00P(1,1,Rz0), MPID%sizei3(0), &
                            MPI_DOUBLE_PRECISION, 0, &
                            MPI_COMM_RT, ierr)

              ! ALI
              if (lALI) then

                ! Share Lambda operator for b-b transitions
                call MPI_BCAST(LambdaL(1,1,Rz0), MPID%sizei10(0), &
                              MPI_DOUBLE_PRECISION, 0, &
                              MPI_COMM_RT, ierr)

                ! Share Lambda operator for b-f transitions
                if (Input%ALI_photo) &
                  call MPI_BCAST(LambdaP(1,1,1,Rz0), MPID%sizei2(0), &
                                MPI_DOUBLE_PRECISION, 0, &
                                MPI_COMM_RT, ierr)

              end if ! ALI
            end if ! To do SEE

          ! Serial
          else

            ! Overwrite
            J00C = J00C_n

            ! Control
            if (laborted) goto 2000

          end if

          ! If we are finished, exit
          if (gooutprd) exit

        end do ! Internal PRD iterations

        ! Control
        if (laborted) goto 2000

        ! Advance populations
        call solveI_SEE(Atom,Rho_old,Atmo,Geom,Input,MRC,lALI,PRDl, &
                        ADD,ADT,NGI,doNG,goout,iter,NG_dim,NG_entry, &
                        NG_scratch,LambdaL,LambdaP, &
                        Stokes,J00,J00S,J00C,J00P)

        ! We can switch now to AD if we had AV input
        if (tbAD.and.ADT) then
          AVI = .False.
          tbAD = .False.
          IRAM = RIRAM
        end if

        !
        ! Save partial solution
        !

        ! Check if we want to store partial results
        if(Input%storei.and.mod(iter,Input%storei_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesolI(Input,iterS,Frec%omega,Geom, &
                        !Atom,Atmo%z,Stokes,J00,J00S,J00C,J00P, &
                         Atom,Atmo%z,Stokes,J00,J00,J00C,J00P, &
                         .False.)

        ! Or have a control check
        else

          ! Control check
          call control

        end if ! Store partial solution

        ! Control
        if (laborted) goto 2000

        ! If in the mandatory non-PRD
        if (iter.le.Input%PRD_delay) goout = .False.

        ! Recover the PRD variable and avoid to go out if
        ! started with CRD
        if (iter.eq.Input%PRD_delay) PRD = PRDl

        ! MPI, share goout
        if (MPID%mpi) &
          call MPI_BCAST(goout, 1, MPI_LOGICAL, 0, &
                         MPI_COMM_RT, ierr)
  
        ! If going out, but with no ALI
        if (goout.and..not.lALI.and.Input%ALI_force) then

          ! Don't leave and activate ALI
          goout = .False.
          force_ALI = .True.

        end if ! Going out without ALI but we are forced to do ALI

        ! If we are finished, everyone exits
        if (goout) exit

      end do ! Iterations

      ! If 1.5D, Master save MRC
      if (pid.eq.0.and.run_mode.eq.1.and.Input%keep_MRC) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(Input%folder)//'/MRC', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1001

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1011
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1011

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = 0e0

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1101

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if ! Save MRC

#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,-1)
#endif
#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P,Input%folder,-1)
#endif

      !
      ! Clean pointers
      !
2000  if (associated(data1M).or.associated(data1O)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
        deallocate(rLineO,rPhotO)
        nullify(rLineO,rPhotO)
      end if
      if (.not.MPID%mpi) then
        if (lp_exu) then
          nullify(p_exu)
        else
          deallocate(p_exu)
          nullify(p_exu)
        end if
      end if

      ! Free PRD
      if (PRDl) call free_e2ord(Red)

      return

1001  umsg = 'Error opening MRC file'
      urou = 'solveI'
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
      if (associated(data1M)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
        deallocate(rLineO,rPhotO)
        nullify(rLineO,rPhotO)
      end if
      if (.not.MPID%mpi) then
        if (lp_exu) then
          nullify(p_exu)
        else
          deallocate(p_exu)
          nullify(p_exu)
        end if
      end if
      if (PRD) call free_e2ord(Red)
      return
1011  umsg = 'Error seeking MRC file'
      urou = 'solveI'
      call MPI_FILE_CLOSE(funit, ierr)
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
      if (associated(data1M)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
        deallocate(rLineO,rPhotO)
        nullify(rLineO,rPhotO)
      end if
      if (.not.MPID%mpi) then
        if (lp_exu) then
          nullify(p_exu)
        else
          deallocate(p_exu)
          nullify(p_exu)
        end if
      end if
      if (PRD) call free_e2ord(Red)
      return
1101  umsg = 'Error writing MRC file'
      urou = 'solveI'
      call MPI_FILE_CLOSE(funit, ierr)
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
      if (associated(data1M)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
        deallocate(rLineO,rPhotO)
        nullify(rLineO,rPhotO)
      end if
      if (.not.MPID%mpi) then
        if (lp_exu) then
          nullify(p_exu)
        else
          deallocate(p_exu)
          nullify(p_exu)
        end if
      end if
      if (PRD) call free_e2ord(Red)
      return

      end subroutine solveI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the approximate RAM space necessary for the
      !! routines solving the self-consistent problem for intensity\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!     Input(Input_class): Structure with configuration data
      subroutine solveI_predict(Atom,Atmo,Frec,Red,Geom,MPID,Input)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input

      ! Local

      logical:: AD,ADD,ADT

      integer:: ia,iaux,NG_dim,if0,if1,nth,nph,npz,ntpz,mfreq


      ! Initialize expected RAM
      TRAMc = 0d0

      ! Predict PRD
      if (PRD) call predict_emissI(Atom,Atmo,Geom,Red)

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! Really go back to angle dependent
      ADT = .not.AV.and.AD

      ! NG quantities

      ! If NG acceleration
      if (Input%NGI) then

        ! Initialize NG dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + Rnz*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRD) then

          ! If we need Stokes
          if (ADD) then

            ! Add contirbution
            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            ! Add contribution
            NG_dim = NG_dim + nfreq*RnZ

          end if ! ADD
        end if ! PRD

        ! Check if it requires too much buffer
        if (dble(NG_dim)*8d-6.le.maxbuffer_NG) then

          ! Master
          if (pid.eq.0) then

            ! Add RAM
            TRAMc = TRAMc + 8d-6*dble(NG_dim*(Input%NGI_ord+2))

          ! Slaves
          else

            ! Add RAM
            TRAMc = TRAMc + 8d-6*dble(NG_dim)

          end if ! Master/slave
        end if ! NG
      end if ! NG

      !
      ! Allocations
      !

      ! MPI
      if (MPID%mpi) then

        ! CPU limits
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)
        mfreq = if1 - if0 + 1

        ! Master
        if (pid.eq.0) then

          ! Alternative
          if (MPID%alternI) then

            ! To receive Intensity chunks
            iaux = MPID%nxfreq*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

            ! To receive Lambda operator for b-b transitions
            iaux = MPID%nxtfreqi*nxb*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

            ! To receive Lambda operator for b-f transitions
            if (Input%ALI_photo) then
              iaux = MPID%nxpfreq*nxb*Rnz
              TRAMc = TRAMc + 8d-6*dble(iaux)
            end if

            ! To receive profile information
            iaux = MPID%nxtfreqi*2*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

          ! Normal
          else

            ! Dimensions
            npz = Geom%nph*Rnz
            ntpz = Geom%nth*npz

            ! To receive Intensity chunks
            iaux = MPID%nxfreq*ntpz
            TRAMc = TRAMc + 8d-6*dble(iaux)

            ! To receive Lambda operator for b-b transitions
            iaux = MPID%nxtfreqi*nxb*ntpz
            if (iaux.lt.1) iaux = 1
            TRAMc = TRAMc + 8d-6*dble(iaux)

            ! To receive Lambda operator for b-f transitions
            if (Input%ALI_photo) then
              iaux = MPID%nxpfreq*nxb*ntpz
              if (iaux.lt.1) iaux = 1
              TRAMc = TRAMc + 8d-6*dble(iaux)
            end if

            ! To receive profile information
            iaux = MPID%nxtfreqi*2*ntpz
            if (iaux.lt.1) iaux = 1
            TRAMc = TRAMc + 8d-6*dble(iaux)

          end if

          ! Norm, BLam and BStk
          TRAMc = TRAMc + 8d-6*dble(nxt*Geom%nph*Geom%nth* &
                                    Rnz*(2 + nxb + 2))

          ! New J00C
          TRAMc = TRAMc + 8d-6*dble(Rnz*nfreq)

        ! Slave
        else

          ! Alternative
          if (MPID%alternI) then

            ! Trivial
            nTh = 1
            nPh = 1

          ! Normal
          else

            ! From Geom
            nTh = Geom%nTh
            nPh = Geom%nPh

          end if

          ! To send Intensity chunks
          TRAMc = TRAMc + 8d-6*dble(Rnz*nPh*nTh*mfreq)

          ! To send profile information
          TRAMc = TRAMc + 8d-6*dble(Rnz*nPh*nTh*Frec%ntfreqi*2)

          ! To send Lambda operator for b-b transitions
          TRAMc = TRAMc + 8d-6*dble(Rnz*nPh*nTh*Frec%ntfreqi)

          ! To send Lambda operator for b-f transitions
          if (Input%ALI_photo) &
            TRAMc = TRAMc + 8d-6*dble(Rnz*nPh*nTh*Frec%npfreq)

          ! Common (Master and slave)
          ! Allocate O pointers
          TRAMc = TRAMc + 8d-6*dble(Frec%ntfreqi*3 + Frec%npfreq)

          ! Allocate M and O pointers for RT coeff
          ! Allocate vector for Lambda operator
          TRAMc = TRAMc + 8d-6*dble(mfreq*7)

        end if

      ! Serial
      else

        ! Allocate O pointers
        TRAMc = TRAMc + 8d-6*dble(Frec%ntfreqi*3 + Frec%npfreq)

        ! Allocate M and O pointers for RT coeff
        TRAMc = TRAMc + 8d-6*dble(nfreq*6)

        ! Initialize exponential is not in memory
        if (.not.(PIRAM.and.Frec%pif1.ge.Frec%pif0)) then
          TRAMc = TRAMc + 8d-6
        end if

        ! Vector for Lambda operator
        TRAMc = TRAMc + 8d-6*dble(nfreq)

        ! New J00C
        TRAMc = TRAMc + 8d-6*dble(Rnz*nfreq)

      end if

      return

      end subroutine solveI_predict

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize arrays and pointers for the self-consistent
      !! solution of the problem for intensity\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration data
      !!             lload(logical): If a previous solution was
      !!                             read\n
      !!            NG_dim(integer): Size of NG entry\n
      !!          NG_entry(integer): Index of NG entries\n
      !!    NG_scracth(double(:,:)): Data for NG iteration\n
      !!               NGI(logical): If doing NG acceleration\n
      !!              doNG(logical): If doing NG iteration\n
      !!             goout(logical): Convergence flag\n
      !!         force_ALI(logical): Force ALI in iteration\n
      !!              PRDl(logical): If doing PRD\n
      !!               ADD(logical): If Stokes needed for PRD-AD\n
      !!               ADT(logical): If initializing PRD-AD with
      !!                             PRD-AA\n
      !!             RIRAM(logical): If storing intensity
      !!                             redistribution\n
      !!            lp_exu(logical): If available pre-computed
      !!                             exponentials\n
      !!               if0(integer): Initial frequency index\n
      !!               if1(integer): Final frequency index\n
      !!               npz(integer): Size in z and phi\n
      !!              ntpz(integer): Size in z, theta, and phi\n
      !!             nsend(integer): Expected messages per iteration\n
      !!        Stokes_r(double(:)): Receiver buffer for Stokes\n
      !!       LambdaL_r(double(:)): Receiver buffer for Lambda 
      !!                             lines\n
      !!       LambdaP_r(double(:)): Receiver buffer for Lambda
      !!                             photoionization\n
      !!          Prof_r(double(:)): Receiver buffer for profiles\n
      !!  Stokes_s(double(:,:,:,:)): Sender buffer for Stokes\n
      !!   rLine_s(double(:,:,:,:)): Sender buffer for Lambda lines\n
      !!   rPhot_s(double(:,:,:,:)): Sender buffer for Lambda
      !!                             photoionizations\n
      !!    Prof_s(double(:,:,:,:)): Sender buffer for profiles\n
      !!        data1M(double(:,:)): RT coeff point M\n
      !!        data1O(double(:,:)): RT coeff point O\n
      !!        data2O(double(:,:)): Profiles point O\n
      !!          rLineO(double(:)): Line Lambda\n
      !!          rPhotO(double(:)): Photoionization Lambda\n
      !!              LO(double(:)): Lambda operator\n
      !!           p_exu(double(:)): Pointer to pre-calculated
      !!                             exponential\n
      !!        J00C_n(double(:,:)): New mean intensity with frequency
      !!                             dependence
      subroutine solveI_init(Atom,Frec,Geom,MPID,Input,lload,NG_dim, &
                             NG_entry,NG_scratch,NGI,doNG,goout, &
                             force_ALI,PRDl,ADD,ADT,RIRAM,lp_exu, &
                             if0,if1,npz,ntpz,nsend, &
                             Stokes_r,LambdaL_r,LambdaP_r,Prof_r, &
                             Stokes_s,rLine_s,rPhot_s,Prof_s, &
                             data1M,data1O,data2O,rLineO,rPhotO, &
                             LO,p_exu,J00C_n)
      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      logical, intent(in):: lload
      logical, intent(out):: doNG,goout,force_ALI,NGI
      logical, intent(out):: PRDl,ADD,ADT,RIRAM,lp_exu
      integer, intent(out):: NG_dim,NG_entry,if0,if1,npz,ntpz,nsend
      double precision, dimension(:,:), &
                        allocatable, intent(out):: NG_scratch
      double precision, dimension(:), &
                        allocatable, intent(out):: Stokes_r
      double precision, dimension(:), &
                        allocatable, intent(out):: LambdaL_r
      double precision, dimension(:), &
                        allocatable, intent(out):: LambdaP_r
      double precision, dimension(:), &
                        allocatable, intent(out):: Prof_r
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00C_n
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(out):: Stokes_s
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Prof_s
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(out):: rLine_s
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(out):: rPhot_s
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data2O
      double precision, dimension(:), &
                        pointer, intent(inout):: rLineO
      double precision, dimension(:), &
                        pointer, intent(inout):: rPhotO
      double precision, dimension(:), &
                        allocatable, intent(out):: LO
      double precision, dimension(:), pointer, intent(inout):: p_exu

      ! Local

      logical:: AD,laux,cfile

      integer:: ia,iaux,ios,nth,nph


      ! Initialize converged flag
      goout = .False.

      ! Initialize force ALI
      force_ALI = .False.

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! Really go back to angle dependent
      ADT = .not.AV.and.AD

      ! If storing redistribution
      RIRAM = IRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AVI = .True.
        if (ADT) IRAM = .False.
      end if

      ! Store if it was PRD (first iteration always CRD)
      PRDl = PRD
      if(.not.lload) PRD = .False.

      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      NG_dim = 0
      doNG = .False.
      NGI = Input%NGI

      ! If NG acceleration
      if (NGI) then

        ! Initialize NG dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + Rnz*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRDl) then

          ! If we need Stokes
          if (ADD) then

            ! Add contribution to size
            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            ! Add contribution to size
            NG_dim = NG_dim + nfreq*RnZ

          end if ! ADD
        end if ! PRD

        ! Check if it requires too much buffer
        if (dble(NG_dim)*8d-6.gt.maxbuffer_NG) then

          ! No NG
          NGI = .False.

          ! Master only
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' # The buffer for NG acceleration '// &
                   'is too big. Not doing NG.'
            call verbose

          end if ! Master
        end if ! Too big of a buffer

        ! If finally doing it, allocate
        if (NGI) then

          ! Master
          if (pid.eq.0) then

            ! Allocate space to store NG data
            allocate(NG_scratch(NG_dim, Input%NGI_ord+2))

          ! Slaves
          else

            ! Allocate space to store NG data
            allocate(NG_scratch(NG_dim, 1))

          end if ! Master/slave
        end if ! NG
      end if ! NG

      !
      ! Allocations
      !

      ! MPI
      if (MPID%mpi) then

        ! lp_exu initialize
        lp_exu = .True.

        ! CPU limits
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! Master
        if (pid.eq.0) then

          ! Alternative
          if (MPID%alternI) then

            ! To receive Intensity chunks
            iaux = MPID%nxfreq*Rnz
            allocate(Stokes_r(iaux))

            ! To receive Lambda operator for b-b transitions
            iaux = MPID%nxtfreqi*nxb*Rnz
            allocate(LambdaL_r(iaux))

            ! To receive Lambda operator for b-f transitions
            if (Input%ALI_photo) then
              iaux = MPID%nxpfreq*nxb*Rnz
              allocate(LambdaP_r(iaux))
            end if

            ! To receive profile information
            iaux = MPID%nxtfreqi*2*Rnz
            allocate(Prof_r(iaux))

            ! Number of messages
            nsend = MPID%nnd*Geom%nTh*Geom%nPh

          ! Normal
          else

            ! Dimensions
            npz = Geom%nph*Rnz
            ntpz = Geom%nth*npz

            ! To receive Intensity chunks
            iaux = MPID%nxfreq*ntpz
            allocate(Stokes_r(iaux))

            ! To receive Lambda operator for b-b transitions
            iaux = MPID%nxtfreqi*nxb*ntpz
            if (iaux.lt.1) iaux = 1
            allocate(LambdaL_r(iaux))

            ! To receive Lambda operator for b-f transitions
            if (Input%ALI_photo) then
              iaux = MPID%nxpfreq*nxb*ntpz
              if (iaux.lt.1) iaux = 1
              allocate(LambdaP_r(iaux))
            end if

            ! To receive profile information
            iaux = MPID%nxtfreqi*2*ntpz
            if (iaux.lt.1) iaux = 1
            allocate(Prof_r(iaux))

            ! Number of messages
            nsend = MPID%nnd

          end if

          ! New J00C
          allocate(J00C_n(nfreq,Rz0:Rz1))

        ! Slave
        else

          ! Alternative
          if (MPID%alternI) then

            ! Trivial
            nTh = 1
            nPh = 1

          ! Normal
          else

            ! From Geom
            nTh = Geom%nTh
            nPh = Geom%nPh

          end if

          ! To send Intensity chunks
          allocate(Stokes_s(if0:if1,Rz0:Rz1,nPh,nTh))

          ! To send profile information
          allocate(Prof_s(Frec%ntfreqi,2,Rz0:Rz1,nPh,nTh))

          ! To send Lambda operator for b-b transitions
          allocate(rLine_s(Frec%ntfreqi,Rz0:Rz1,nPh,nTh))

          ! To send Lambda operator for b-f transitions
          if (Input%ALI_photo) &
            allocate(rPhot_s(Frec%npfreq,Rz0:Rz1,nPh,nTh))

          ! Common (Master and slave)
          ! Allocate O pointers
          allocate(data2O(Frec%ntfreqi,2))
          allocate(rLineO(Frec%ntfreqi))
          allocate(rPhotO(Frec%npfreq))

          ! Allocate M and O pointers for RT coeff
          allocate(data1M(MPID%nf(pid),3))
          allocate(data1O(MPID%nf(pid),3))

          ! Allocate vector for Lambda operator
          allocate(LO(if0:if1))

        end if ! Master/slave

      ! Serial
      else

        ! CPU limits
        if0 = 1
        if1 = nfreq

        ! Allocate O pointers
        allocate(data2O(Frec%ntfreqi,2))
        allocate(rLineO(Frec%ntfreqi))
        allocate(rPhotO(Frec%npfreq))

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(nfreq,3))
        allocate(data1O(nfreq,3))

        ! Initialize exponential is not in memory
        if (PIRAM.and.Frec%pif1.ge.Frec%pif0) then
          lp_exu = .True.
        ! No exu allocated or needed
        else
          lp_exu = .False.
          allocate(p_exu(1))
        end if

        ! Vector for Lambda operator
        allocate(LO(nfreq))

        ! New J00
        allocate(J00C_n(nfreq,Rz0:Rz1))

      end if ! MPI/serial

      !
      ! Initialization messages
      !

      ! Global Master
      if (gpid.eq.0) then

        ! Announce we are starting
        umsg = '    Iteration          MRC(rho^0_0) Atom_index '// &
               'Level_index Height_index Height(km)'
        call verbose

      end if ! Global master

      !
      ! Write to file
      !

      ! Only global Máster
      if (gpid.eq.0.and.Input%keep_MRC) then

        ! If appending
        if (Input%appendMRC) then

          ! Check if file exists
          inquire(file=trim(Input%folder)//'/MRCI', exist=laux)

          ! Create if it does not exist
          cfile = .not.laux

        ! Not appending
        else

          ! Create
          cfile = .True.

        end if ! Appending or not the MRC file

        ! If we need to create it
        if (cfile) then

          ! Create file
          open(800, file=trim(Input%folder)//'/MRCI', &
               action='write',iostat=ios,err=1000)

          ! Write header
          write(800,'(A)',err=1100) &
                       '!   Iteration          MRC(rho^0_0) '// &
                       'Atom_index Level_index Height_index '// &
                       'Height(km)'

          ! Close
          close(800)

        end if ! Create file
      end if ! Global Master keeping MRC

      return

1000  umsg = 'Error opening MRCI file'
      urou = 'solveI_init'
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
1100  umsg = 'Error writing MRCI file'
      urou = 'solveI_init'
      close(800)
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
      return

      end subroutine solveI_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Parallel-master task for the self-consistent problem for
      !! intensity\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!             lALI(logical): If doing ALI\n
      !!            nsend(integer): Expected messages per iteration\n
      !!             iter(integer): Iteration number\n
      !!            iterr(integer): Radiation interation number\n
      !!              npz(integer): Size in z and phi\n
      !!             ntpz(integer): Size in z, theta, and phi\n
      !!       Stokes_r(double(:)): Receiver buffer for Stokes\n
      !!      LambdaL_r(double(:)): Receiver buffer for Lambda 
      !!                            lines\n
      !!      LambdaP_r(double(:)): Receiver buffer for Lambda
      !!                            photoionization\n
      !!         Prof_r(double(:)): Receiver buffer for profiles\n
      !!    LambdaL(double(:,:,:)): Lambda operator for lines\n
      !!  LambdaP(double(:,:,:,:)): Lambda operator for
      !!                            photoionizations\n
      !!   Stokes(double(:,:,:,:)): Intensity\n
      !!         J00C(double(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!          J00(double(:,:)): Mean intensity integrated over the
      !!                            absorption profile\n
      !!         J00S(double(:,:)): Mean intensity integrated over the
      !!                            emission profile\n
      !!       J00P(double(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine solveI_manager(Atom,Atmo,Frec,Geom,MPID,Input,lALI, &
                                nsend,iter,iterr,npz,ntpz, &
                                Stokes_r,LambdaL_r, &
                                LambdaP_r,Prof_r,LambdaL,LambdaP, &
                                Stokes,J00C,J00,J00S,J00P)
      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      logical, intent(in):: lALI
      integer, intent(in):: iter,iterr,nsend,npz,ntpz
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: Stokes_r
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: LambdaL_r
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: LambdaP_r
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: Prof_r
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target, intent(out):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(out):: J00C
      double precision, dimension(nxt,Rz0:Rz1), intent(out):: J00
      double precision, dimension(nxt,Rz0:Rz1), intent(out):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1), &
                        intent(out):: J00P
      double precision, dimension(nxb,nxt,Rz0:Rz1), &
                        intent(out):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1), &
                        intent(out):: LambdaP

      ! Local

      logical:: deal

      integer:: id,info_b,if0l,if1l,ip0l,ip1l,nfl,nftl,nfpl
      integer:: itpz,ith,iph,iz,ia,itran,jtran,ftran,fjtran,fftran
      integer, dimension(3):: info_c

      double precision:: WA,daux,dsm
      double precision, &
             dimension(2,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: BStk
      double precision, &
             dimension(nxb,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: BLam
      double precision, &
             dimension(2,nxt,Geom%nph,Geom%nth,Rz0:Rz1):: Norm

      ! Pointers

      double precision, dimension(:), pointer:: p_exu
      double precision, dimension(:,:), pointer:: p_MStk
      double precision, dimension(:,:), pointer:: p_MrLine,p_MrPhot
      double precision, dimension(:,:,:), pointer:: p_MProf


      ! Nullify
      nullify(p_exu,p_MStk,p_MrLine,p_MrPhot,p_MProf)

      ! Initialize
      J00 = 0d0
      J00P = 0d0
      J00C = 0d0
      BStk = 0d0
      BLam = 0d0
      Norm = 0d0
      if (lALI) then
        LambdaL = 0d0
        LambdaP = 0d0
      end if

      !
      ! Alternative MPI
      !
      if (MPID%alternI) then

        ! Expected packages
        do id=1,nsend

          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_c(1),3,MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, &
                          MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Get indexes
          info_b = info_c(1)
          ith = info_c(2)
          iph = info_c(3)

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive intensity
          do while (.True.)
            call MPI_recv(Stokes_r(1), MPID%sizei4(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), MPID%sizei5(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          1+info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! If ALI
          if (lALI) then

            ! Receive Lambda operator for b-b transition
            do while (.True.)
              call MPI_recv(LambdaL_r(1), &
                            MPID%sizei9(info_b), &
                            MPI_DOUBLE_PRECISION, info_b, &
                            2+info_b, &
                            MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! Receive Lambda operator for b-f transition
            if (Input%ALI_photo) then
              do while (.True.)
                call MPI_recv(LambdaP_r(1), &
                              MPID%sizei0(info_b), &
                              MPI_DOUBLE_PRECISION, info_b, &
                              3+info_b, &
                              MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do
            end if

          end if ! ALI

          ! If measuring performance
          if (Input%mpi_perf) &
            call report_mpi_timeI(Input%folder,Input%ID, &
                                  info_b,iter,iterr,.True.)

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          ip0l = Frec%Mpif0(info_b)
          ip1l = Frec%Mpif1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreqi(info_b)
          nfpl = Frec%Mnpfreq(info_b)

          ! Pointers
          p_MStk(if0l:if1l,1:Rnz) => &
                              Stokes_r(1:MPID%sizei4(info_b))
          p_MProf(1:nftl,1:2,1:Rnz) => &
                                Prof_r(1:MPID%sizei5(info_b))
          if (lALI) then
            p_MrLine(1:nftl,1:Rnz) => &
                             LambdaL_r(1:MPID%sizei9(info_b))
            if (Input%ALI_photo) &
              p_MrPhot(1:nfpl,1:Rnz) => &
                             LambdaP_r(1:MPID%sizei0(info_b))
          end if

          ! Get angular weight
          WA = Geom%W_mu(ith)*Geom%W_mux(iph)

          ! Initialize deal variable
          deal = .False.

          ! Each height
          do iz=Rz0,Rz1

            ! Determine where to store intensity
            if (KSTK.or.iz.eq.giz0) &
              Stokes(if0l:if1l,iph,ith,iz) = &
                                          p_MStk(if0l:if1l,iz)
            ! Point to exu values
            if (PIRAM.and.ip1l.ge.ip0l) then
              p_exu => Frec%exu(ip0l:ip1l,iz)
            else
              allocate(p_exu(1))
              deal = .True.
            end if

            ! Calculate frequency integral for b-b quantities
            call FIntI_line(Atom,MPID,Frec%W_freq,info_b, &
                            p_MStk(:,iz),p_MrLine(:,iz), &
                            p_MProf(:,:,iz), &
                            Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz), &
                            BLam(:,:,iph,ith,iz), &
                            lALI)

            ! Calculate rest of integrals
            call FIntI_rest(Atom,MPID,Frec%omega, &
                            Frec%W_freq,ip0l,ip1l, &
                            Atmo%T(iz),info_b,WA, &
                            p_MStk(:,iz),p_MrPhot(:,iz), &
                            J00P(:,:,iz),J00C(:,iz), &
                            LambdaP(:,:,:,iz),lALI, &
                            Input%ALI_photo,p_exu)

            ! Nullify exponential pointer
            if (deal) deallocate(p_exu)
            nullify(p_exu)

          end do ! heights
        end do ! Expected packages

      !
      ! Normal MPI
      !
      else

        ! Expected packages
        do id=1,nsend

          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_b,1,MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive intensity
          do while (.True.)
            call MPI_recv(Stokes_r(1), MPID%sizei4(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), MPID%sizei5(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          1+info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Only if ALI iteration
          if (lALI) then

            ! Receive Lambda operator for b-b transition
            do while (.True.)
              call MPI_recv(LambdaL_r(1), MPID%sizei9(info_b), &
                            MPI_DOUBLE_PRECISION, info_b, &
                            2+info_b, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! Receive Lambda operator for b-f transition
            if (Input%ALI_photo) then
              do while (.True.)
                call MPI_recv(LambdaP_r(1), MPID%sizei0(info_b), &
                              MPI_DOUBLE_PRECISION, info_b, &
                              3+info_b, MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do
            end if
          end if ! ALI iteration

          ! If measuring performance
          if (Input%mpi_perf) &
            call report_mpi_timeI(Input%folder,Input%ID, &
                                  info_b,iter,iterr,.True.)

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          ip0l = Frec%Mpif0(info_b)
          ip1l = Frec%Mpif1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreqi(info_b)
          nfpl = Frec%Mnpfreq(info_b)

          ! Pointers
          p_MStk(if0l:if1l,1:ntpz) => &
                                  Stokes_r(1:MPID%sizei4(info_b))
          p_MProf(1:nftl,1:2,1:ntpz) => &
                                    Prof_r(1:MPID%sizei5(info_b))
          ! Point to actual Lambda data
          if (lALI) then
            p_MrLine(1:nftl,1:ntpz) => &
                                 LambdaL_r(1:MPID%sizei9(info_b))
            if (Input%ALI_photo) &
              p_MrPhot(1:nfpl,1:ntpz) => &
                                 LambdaP_r(1:MPID%sizei0(info_b))
          else
            ! Point to whatever
            p_MrLine(1:1,1:ntpz) => Prof_r(1:ntpz)
            if (Input%ALI_photo) &
              p_MrPhot(1:1,1:ntpz) => Prof_r(1:ntpz)
          end if

          ! Initialize deal
          deal = .False.

          ! If we are aborting, just skip this
          if (laborted) cycle

          ! Compute line quantities
          do itpz=1,ntpz

            ! Get indexes
            ith = (itpz-1)/npz
            iph = (itpz - npz*ith - 1)/Rnz
            iz = itpz - Rnz*iph - npz*ith + Rz0 - 1
            ith = ith + 1
            iph = iph + 1

            ! Determine where to store intensity
            if (KSTK.or.iz.eq.Rz0) &
              Stokes(if0l:if1l,iph,ith,iz) = &
                                            p_MStk(if0l:if1l,itpz)

            ! Calculate frequency integral for b-b quantities
            call FIntI_line(Atom,MPID,Frec%W_freq,info_b, &
                            p_MStk(:,itpz),p_MrLine(:,itpz), &
                            p_MProf(:,:,itpz), &
                            Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz), &
                            BLam(:,:,iph,ith,iz), &
                            lALI)

            ! Point to exu values
            if (PIRAM.and.ip1l.ge.ip0l) then
              p_exu => Frec%exu(ip0l:ip1l,iz)
            else
              allocate(p_exu(1))
              deal = .True.
            end if

            ! Get angular weight
            WA = Geom%W_mu(ith)*Geom%W_mux(iph)

            !
            ! Calculate rest of integrals
            !
            call FIntI_rest(Atom,MPID,Frec%omega, &
                            Frec%W_freq,ip0l,ip1l, &
                            Atmo%T(iz),info_b,WA, &
                            p_MStk(:,itpz),p_MrPhot(:,itpz), &
                            J00P(:,:,iz),J00C(:,iz), &
                            LambdaP(:,:,:,iz),lALI, &
                            Input%ALI_photo,p_exu)


            ! Nullify pointer
            if (deal) deallocate(p_exu)
            nullify(p_exu)

          end do ! Heights and directions
        end do ! Expected packages

      end if ! Type of MPI

      ! Nullify pointers
      nullify(p_MStk,p_MProf,p_MrLine,p_MrPhot)

      !
      ! Apply weights to J00, J00S, and Lambda operator
      ! and normalize
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each polar direction
        do ith=1,Geom%nTh

          ! For each azimuthal direction
          do iph=1,Geom%nph

            ! Get the angular integral weight
            WA = Geom%W_mu(ith)*Geom%W_mux(iph)

            ! For each atom
            do ia=1,nA

              ! For each FS transition
              do ftran=1,Atom(ia)%nftran

                ! Apply shift
                fjtran = ftran + Atom(ia)%tfshift

                ! Get the weight
                if (Norm(1,fjtran,iph,ith,iz).gt.0d0) then

                  ! Inverse norm and angular weight
                  daux = WA/Norm(1,fjtran,iph,ith,iz)

                  ! Integrate angle
                  J00(fjtran,iz) = J00(fjtran,iz) + &
                                   BStk(1,fjtran,iph,ith,iz)* &
                                   daux

                  ! For each transition blended with ftran
                  if (lALI) then

                    ! For each blend
                    do fftran=1,nxb

                      ! Add contribution
                      LambdaL(fftran,fjtran,iz) = &
                               LambdaL(fftran,fjtran,iz) + &
                               BLam(fftran,fjtran,iph,ith,iz)* &
                               daux

                    end do ! Blends

                  end if ! ALI

                end if

               !! If there is stimulated emission
               !if (stm) then

               !  ! Get the weight
               !  if (Norm(2,fjtran,iph,ith,iz).gt.0d0) then

               !    daux = WA/Norm(2,fjtran,iph,ith,iz)

               !    ! Integrate angle
               !    J00S(fjtran,iz) = J00S(fjtran,iz) + &
               !                    BStk(2,fjtran,iph,ith,iz)* &
               !                    daux

               !  end if

               !end if ! Stimulated emission

              end do ! FS transition
            end do ! Atoms
          end do ! azimuthal directions
        end do ! polar directions

        !
        ! Add the Saha factor (ne*Zeta) to J00P
        !

        ! Argument of the exponential
        WA = fktoJ/kb/Atmo%T(iz)

        ! Part that does not depend on the line
        daux = cSaha*Atmo%ne(iz)/(Atmo%T(iz)**(1.5d0))

        ! For each atom
        do ia=1,nA

          ! For each b-f transition
          do itran=1,Atom(ia)%nphot

            ! Apply shift
            jtran = itran + Atom(ia)%pshift

            ! Calculate the multiplicative factor
            dsm = daux*exp(Atom(ia)%phot(itran)%edge*WA)* &
                  Atom(ia)%phot(itran)%glu

            ! Apply it to the emission integral
            J00P(jtran,2,iz) = J00P(jtran,2,iz)*dsm

            ! Apply it to the emission Lambda operator
            if (lALI.and.Input%ALI_photo) &
              LambdaP(:,jtran,2,iz) = LambdaP(:,jtran,2,iz)*dsm

          end do ! b-f transitions
        end do ! atoms
      end do ! heights

      ! Return
      return

      ! Deceive compiler
      WA = J00S(1,1)

      end subroutine solveI_manager

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the RTE for intensity for the angular quadrature\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structure with the LTE line
      !!                              data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!               lALI(logical): If doing ALI\n
      !!               ALIp(logical): If doing ALI in bound-free\n
      !!             lp_exu(logical): If available pre-computed
      !!                              exponentials\n
      !!               if0(integer): Initial frequency index\n
      !!               if1(integer): Final frequency index\n
      !!  Stokes_s(double(:,:,:,:)): Sender buffer for Stokes\n
      !!   rLine_s(double(:,:,:,:)): Sender buffer for Lambda lines\n
      !!   rPhot_s(double(:,:,:,:)): Sender buffer for Lambda
      !!                             photoionizations\n
      !!    Prof_s(double(:,:,:,:)): Sender buffer for profiles\n
      !!        data1M(double(:,:)): RT coeff point M\n
      !!        data1O(double(:,:)): RT coeff point O\n
      !!        data2O(double(:,:)): Profiles point O\n
      !!          rLineO(double(:)): Line Lambda\n
      !!          rPhotO(double(:)): Photoionization Lambda\n
      !!              LO(double(:)): Lambda operator\n
      !!     LambdaL(double(:,:,:)): Lambda operator for lines\n
      !!   LambdaP(double(:,:,:,:)): Lambda operator for
      !!                             photoionizations\n
      !!           p_exu(double(:)): Pointer to pre-calculated
      !!                             exponential\n
      !!    Stokes(double(:,:,:,:)): Intensity\n
      !!           J00(double(:,:)): Mean intensity integrated over
      !!                             the absorption profile\n
      !!          J00S(double(:,:)): Mean intensity integrated over
      !!                             the emission profile\n
      !!          J00C(double(:,:)): Mean intensity with frequency
      !!                             dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates\n
      !!        J00C_n(double(:,:)): New mean intensity with frequency
      !!                             dependence
      subroutine solveI_RT(Atom,LTElines,Atmo,Cont,Frec,Red,Geom, &
                           MPID,lALI,ALIp,lp_exu,if0,if1,Stokes_s, &
                           rLine_s,rPhot_s,Prof_s,data1M,data1O, &
                           data2O,rLineO,rPhotO,LO,LambdaL,LambdaP, &
                           p_exu,Stokes,J00,J00S,J00C,J00P,J00C_n)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: lALI,ALIp,lp_exu
      integer, intent(in):: if0,if1
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: Stokes_s
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Prof_s
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: rLine_s
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: rPhot_s
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data2O
      double precision, dimension(:), pointer, intent(inout):: rLineO
      double precision, dimension(:), pointer, intent(inout):: rPhotO
      double precision, dimension(:), allocatable, &
                        target, intent(inout):: LO
      double precision, dimension(:), pointer, intent(inout):: p_exu
      double precision, &
             dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             intent(inout):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: J00C
      double precision, dimension(nxt,Rz0:Rz1), intent(inout):: J00
      double precision, dimension(nxt,Rz0:Rz1), intent(inout):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1), &
                        intent(inout):: J00P
      double precision, dimension(nxb,nxt,Rz0:Rz1), &
                        intent(inout):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1), &
                        intent(inout):: LambdaP
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C_n

      ! Local

      integer:: ith,iph,jdir,jth,jph,iz,diz,iz0,iz1,m,o,p,nodir,mjdir
      integer:: ia,itran,jtran,ftran,fftran,jftran,ifreq,iil,iip,mfreq
      integer, dimension(3):: info_c

      double precision:: mu_inv,dsm,dsp,ct,st,cc,sc,vfac

      ! Pointers

      double precision, dimension(:,:), pointer:: data1P
      double precision, dimension(:,:), pointer:: data2P
      double precision, dimension(:), pointer:: rLineP
      double precision, dimension(:), pointer:: rPhotP
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP, p_LO


      ! Nullify pointers
      nullify(p_K0M,p_SM,p_StkM,p_K0O,p_SP,p_StkO,p_K0P,p_SP)
      nullify(data1P,data2P,rLineP,rPhotP)

      ! Allocate and initialize in serial
      if (pid.eq.0) then
        J00 = 0d0
        J00S = 0d0
        J00P = 0d0
        J00C_n = 0d0
        if (lALI) then
          LambdaL = 0d0
          if (ALIp) LambdaP = 0d0
        end if
      end if

      ! Frequency size
      mfreq = if1 - if0 + 1

      ! Initialize Doppler shift
      vfac = 1d0

      ! Maximum direction index for norm
      mjdir = ubound(Red%idzao,3)

      !
      ! Ratiation Transfer
      !

      !  For each polar direction
      do ith=1,Geom%nTh

        ! If error
        if (laborted) exit

        ! Calculate inverse of cosine of polar direction
        mu_inv = 1d0/Geom%V_mu(ith)

        ! Determine the direction of propagation for indexes
        diz = -int(sign(1d0, Geom%V_mu(ith)))

        ! Determine the first and last height indexes to run
        ! over
        iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
        iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

        ! Trigonometry for Doppler shift
        if (dyn) then
          ct = Geom%v_mu(ith)
          st = sqrt(1d0 - ct*ct)
        end if

        ! For each azimuthal direction
        do iph=1,Geom%nPh

          ! If error
          if (laborted) exit

          ! Alternative MPI
          if (MPID%alternI) then
            jth = 1
            jph = 1
          ! Normal MPI
          else
            jth = ith
            jph = iph
          end if

          ! Get direction index
          jdir = Geom%i_geom(iph,ith)
          nodir = min(mjdir,Geom%njdir,jdir)


          !
          ! First height
          !

          ! Get Doppler shift
          if (dyn) then

            ! Trigonometry
            cc = Geom%v_mux(iph)
            sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

            ! Amplitude
            vfac = sqrt(Atmo%vx(iz0)*Atmo%vx(iz0) + &
                        Atmo%vy(iz0)*Atmo%vy(iz0) + &
                        Atmo%vz(iz0)*Atmo%vz(iz0))

            ! Doppler shift
            if (vfac.gt.TINYVEL) then

              ! Shift
              vfac = 1d0 - Atmo%vx(iz0)*st*cc - &
                           Atmo%vy(iz0)*st*sc - &
                           Atmo%vz(iz0)*ct

            ! Static
            else

              ! No shift
              vfac = 1d0

            end if
          end if ! If dynamic

          ! If going down, get top boundary
          if(diz.eq.1)then

            ! Call top boundary
            call topI(if0,if1,data1M(:,3))

          ! If going up, get bottom boundary
          else

            ! Call bottom boundary
            call bottomI(Frec%omega,Atmo%T(iz0), &
                         vfac,if0,if1,data1M(:,3))

          endif ! propagation direction

          ! Identify current height
          o = iz0

          ! Calculate radiative coefficients
          call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                        o,jdir,nodir,if0,if1, &
                        J00C(:,o),Cont%ndir, &
                        Cont%c(:,:,:,o), &
                        rLineO(:),rPhotO(:), &
                        data1M(:,1:2),data2O(:,:),.True.)

          ! If error
          if (laborted) exit

          ! If MPI
          if (pid.gt.0) then

            !
            ! Store in buffer
            !

            ! Intensity
            Stokes_s(:,o,jph,jth) = data1M(:,3)

            ! Profiles
            Prof_s(:,:,o,jph,jth) = data2O(:,:)

            ! If ALI iteration
            if (lALI) then

              ! b-b Lambda operator (bottom boundary does not
              ! contribute)
              rLine_s(:,o,jph,jth) = 0d0

              ! b-f Lambda operator (bottom boundary does not
              ! contribute)
              if (ALIp) rPhot_s(:,o,jph,jth) = 0d0

            end if

          ! If serial
          else

            ! Save Stokes
            if (KSTK) Stokes(:,iph,ith,o) = data1M(:,3)

            ! The initial point does not contribute to Lambda
            ! operator
            if (lALI) then
              rLineO = 0d0
              if (ALIp) rPhotO = 0d0
            end if

            ! Point to exu values
            if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

            !
            ! Calculate integrals
            !
            call JcalcI(Atom,Geom,Frec%omega,Frec%W_freq, &
                        Frec%pif0,Frec%pif1, &
                        Atmo%T(o),Atmo%ne(o),iph,ith, &
                        data1M(:,3),rLineO,rPhotO,data2O, &
                        J00(:,o),J00S(:,o),J00P(:,:,o), &
                        J00C_n(:,o),LambdaL(:,:,o), &
                        LambdaP(:,:,:,o),lALI,ALIp,p_exu)

          end if ! MPI/serial

          ! Identify next height
          p = iz0 + diz

          ! Get Doppler shift
          if (dyn) then

            ! Amplitude
            vfac = sqrt(Atmo%vx(p)*Atmo%vx(p) + &
                        Atmo%vy(p)*Atmo%vy(p) + &
                        Atmo%vz(p)*Atmo%vz(p))

            ! Doppler shift
            if (vfac.gt.TINYVEL) then

              ! Shift
              vfac = 1d0 - Atmo%vx(p)*st*cc - &
                           Atmo%vy(p)*st*sc - &
                           Atmo%vz(p)*ct

            ! Static
            else

              ! No shift
              vfac = 1d0

            end if
          end if ! If dynamic

          ! Calculate radiative coefficients
          call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                        p,jdir,nodir,if0,if1, &
                        J00C(:,p),Cont%ndir, &
                        Cont%c(:,:,:,p), &
                        rLineO(:),rPhotO(:), &
                        data1O(:,1:2),data2O(:,:),.True.)

          ! If error
          if (laborted) exit


          !
          ! Intermediate heights
          !

          ! For each height
          do iz=iz0,iz1,diz

            ! We treat the boundaries outside
            if (iz.eq.iz0.or.iz.eq.iz1) cycle

            ! Allocate P pointers
            allocate(data1P(MPID%nf(pid),3))
            allocate(data2P(Frec%ntfreqi,2))
            allocate(rLineP(Frec%ntfreqi))
            allocate(rPhotP(Frec%npfreq))

            ! Identify heights
            m = iz - diz
            o = iz
            p = iz + diz

            ! Calculate distance to previous point
            dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

            ! Calculate distance to the next point
            dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

            ! If tau scale
            if (ztau) then

              ! Get geometric distance
              dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(m))
              dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(p))
            end if

            ! Get Doppler shift
            if (dyn) then

              ! Amplitude
              vfac = sqrt(Atmo%vx(p)*Atmo%vx(p) + &
                          Atmo%vy(p)*Atmo%vy(p) + &
                          Atmo%vz(p)*Atmo%vz(p))

              ! Doppler shift
              if (vfac.gt.TINYVEL) then

                ! Shift
                vfac = 1d0 - Atmo%vx(p)*st*cc - &
                             Atmo%vy(p)*st*sc - &
                             Atmo%vz(p)*ct

              ! Static
              else

                ! No shift
                vfac = 1d0

              end if
            end if ! If dynamic

            ! Calculate radiative coefficients
            call RTCoeffI(Frec,Red,Atom,LTElines,Atmo, &
                          vfac,p,jdir,nodir,if0,if1, &
                          J00C(:,p),Cont%ndir, &
                          Cont%c(:,:,:,p), &
                          rLineP(:),rPhotP(:), &
                          data1P(:,1:2),data2P(:,:),.True.)

            ! Point to the data
            p_K0M  => data1M(:,1)
            p_SM   => data1M(:,2)
            p_StkM => data1M(:,3)
            p_K0O  => data1O(:,1)
            p_SO   => data1O(:,2)
            p_StkO => data1O(:,3)
            p_K0P  => data1P(:,1)
            p_SP   => data1P(:,2)
            p_LO => LO(:)

            ! Apply short characteristics BESSER
            call RTStepI(o,ith,iph,mfreq, &
                         dsm,dsp,p_K0M,p_SM,p_K0O, &
                         p_SO,p_K0P,p_SP,p_StkM, &
                         p_StkO,p_LO,lALI,.True.)

            !
            ! Combine the value of lambda operator with the
            ! transition strength
            !
            if (lALI) then

              ! Initialize indexes for rLine and rPhot
              iil = 0
              iip = 0

              ! For each atom
              do ia=1,nA

                ! For each b-b trantision
                do itran=1,Atom(ia)%ntran

                  ! If this CPU does not have frequencies in
                  ! this line, skip
                  if (Atom(ia)%fflag(itran)%absent) cycle

                  ! For each FS transition
                  do ftran=1,Atom(ia)%fst(itran)%nt

                    ! Get the sequential index of FS transition
                    fftran = Atom(ia)%ifst_ij(ftran,itran)

                    ! Apply shift
                    jftran = fftran + Atom(ia)%tfshift

                    ! For each frequency
                    do ifreq=Atom(ia)%if0(itran), &
                             Atom(ia)%if1(itran)

                      ! Advance and store scaled
                      iil = iil + 1
                      rLineO(iil) = LO(ifreq)*rLineO(iil)

                    end do ! frequency
                  end do ! FS transition
                end do ! b-b transition

                ! If ALI photoionization
                if (ALIp) then

                  ! For each b-f transition
                  do itran=1,Atom(ia)%nphot

                    ! If this CPU does not have frequencies in
                    ! this transition, skip
                    if (Atom(ia)%phot(itran)%absent) cycle

                    ! Apply shift
                    jtran = itran + Atom(ia)%pshift

                    ! For each frequency
                    do ifreq=Atom(ia)%phot(itran)%if0, &
                             Atom(ia)%phot(itran)%if1

                      ! Advance and store scaled
                      iip = iip + 1
                      rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                    end do ! frequency
                  end do ! b-f transition

                end if ! If ALI photoionization

              end do ! atom


              ! If MPI
              if (pid.gt.0) then

                ! b-b Lambda operator
                rLine_s(:,o,jph,jth) = rLineO(:)

                ! Send b-f Lambda operator
                if (ALIp) &
                  rPhot_s(:,o,jph,jth) = rPhotO(:)

              end if ! MPI

            end if ! ALI

            ! If MPI
            if (pid.gt.0) then

              !
              ! Store in buffer
              !

              ! Intensity
              Stokes_s(:,o,jph,jth) = data1O(:,3)

              ! Profiles
              Prof_s(:,:,o,jph,jth) = data2O(:,:)

            ! If serial
            else

              ! Save Stokes
              if (KSTK) Stokes(:,iph,ith,o) = data1O(:,3)

              ! Point to exu values
              if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

              !
              ! Calculate integrals
              !
              call JcalcI(Atom,Geom,Frec%omega,Frec%W_freq, &
                          Frec%pif0,Frec%pif1, &
                          Atmo%T(o),Atmo%ne(o),iph,ith, &
                          data1O(:,3),rLineO,rPhotO,data2O, &
                          J00(:,o),J00S(:,o),J00P(:,:,o), &
                          J00C_n(:,o),LambdaL(:,:,o), &
                          LambdaP(:,:,:,o),lALI,ALIp,p_exu)

            end if

            ! Shift data (O->M, P->O)
            deallocate(data1M,data2O)
            data1M => data1O
            data1O => data1P
            data2O => data2P
            nullify(data1P,data2P)
            if (lALI) then
              deallocate(rLineO,rPhotO)
              rLineO => rLineP
              rPhotO => rPhotP
            else
              deallocate(rLineP,rPhotP)
            end if
            nullify(rLineP,rPhotP)

            ! If error
            if (laborted) exit

          end do ! Intermediate heights

          ! If error
          if (laborted) exit


          !
          ! Last height
          !

          ! Identify heights
          m = iz1 - diz
          o = iz1

          ! Calculate distance to previous point
          dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

          ! If tau scale
          if (ztau) &
            dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                Atmo%chi500(m))

          ! Point to the data
          p_K0M  => data1M(:,1)
          p_SM   => data1M(:,2)
          p_StkM => data1M(:,3)
          p_K0O  => data1O(:,1)
          p_SO   => data1O(:,2)
          p_StkO => data1O(:,3)
          p_LO => LO(:)

          ! Apply short characteristics LINEAR
          call RTStepI(o,ith,iph,MPID%nf(pid), &
                       dsm,dsp,p_K0M,p_SM,p_K0O, &
                       p_SO,p_K0P,p_SP,p_StkM, &
                       p_StkO,p_LO,lALI,.False.)

          ! If error
          if (laborted) exit


          !
          ! Combine the value of lambda operator with the
          ! transition strength
          !

          if (lALI) then

            ! Initialize indexes for rLine and rPhot
            iil = 0
            iip = 0

            ! For each atom
            do ia=1,nA

              ! For each b-b trantision
              do itran=1,Atom(ia)%ntran

                ! If this CPU does not have frequencies in
                ! this line, skip
                if (Atom(ia)%fflag(itran)%absent) cycle

                ! For each FS transition
                do ftran=1,Atom(ia)%fst(itran)%nt

                  ! Get the sequential index of FS transition
                  fftran = Atom(ia)%ifst_ij(ftran,itran)

                  ! Apply atomic shift
                  jftran = fftran + Atom(ia)%tfshift

                  ! For each frequency
                  do ifreq=Atom(ia)%if0(itran), &
                           Atom(ia)%if1(itran)

                    iil = iil + 1
                    rLineO(iil) = LO(ifreq)*rLineO(iil)

                  end do ! frequency
                end do ! FS transition
              end do ! b-b transition

              ! If ALI photoionization
              if (ALIp) then

                ! For each b-f transition
                do itran=1,Atom(ia)%nphot

                  ! If this CPU does not have frequencies in
                  ! this transition, skip
                  if (Atom(ia)%phot(itran)%absent) cycle

                  ! Apply atomic shift
                  jtran = itran + Atom(ia)%pshift

                  ! For each frequency
                  do ifreq=Atom(ia)%phot(itran)%if0, &
                           Atom(ia)%phot(itran)%if1

                    iip = iip + 1
                    rPhotO(iip) = LO(ifreq)*rPhotO(iip)

                  end do ! frequency
                end do ! b-f transition

              end if ! If ALI photoionization

            end do ! atom

            ! MPI
            if (pid.gt.0) then

              ! b-b Lambda operator
              rLine_s(:,o,jph,jth) = rLineO(:)

              ! b-f Lambda operator
              if (ALIp) &
                rPhot_s(:,o,jph,jth) = rPhotO(:)

            end if

          end if ! ALI

          ! MPI
          if (pid.gt.0) then

            !
            ! Store in buffer
            !

            ! Intensity
            Stokes_s(:,o,jph,jth) = data1O(:,3)

            ! Profiles
            Prof_s(:,:,o,jph,jth) = data2O(:,:)

            !
            ! Alternative
            !

            if (MPID%alternI) then

              !
              ! Send to master
              !

              ! Send indexes
              info_c = (/ pid, ith, iph /)
              do while (.True.)
                call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

              ! Send Stokes
              do while (.True.)
                call MPI_SEND(Stokes_s(if0,Rz0,1,1), &
                              MPID%sizei4(pid), &
                              MPI_DOUBLE_PRECISION, 0, pid, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

              ! Send profiles
              do while (.True.)
                call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                              MPID%sizei5(pid), &
                              MPI_DOUBLE_PRECISION, &
                              0, 1+pid, MPI_COMM_RT, &
                              ierr)
                if (ierr.eq.0) exit
              end do

              ! ALI
              if (lALI) then

                ! Send b-b Lambda operator
                do while (.True.)
                  call MPI_SEND(rLine_s(1,Rz0,1,1), &
                                MPID%sizei9(pid), &
                                MPI_DOUBLE_PRECISION, &
                                0, 2+pid, MPI_COMM_RT, &
                                ierr)
                  if (ierr.eq.0) exit
                end do

                ! Send b-f Lambda operator
                if (ALIp) then
                  do while (.True.)
                    call MPI_SEND(rPhot_s(1,Rz0,1,1), &
                                  MPID%sizei0(pid), &
                                  MPI_DOUBLE_PRECISION, &
                                  0, 3+pid, MPI_COMM_RT, &
                                  ierr)
                    if (ierr.eq.0) exit
                  end do
                end if

              end if ! ALI
            end if ! Alternative MPI

          ! Serial
          else

            if (KSTK.or.o.eq.1) &
              Stokes(:,iph,ith,o) = data1O(:,3)

            ! Point to exu values
            if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

            !
            ! Calculate integrals
            !
            call JcalcI(Atom,Geom,Frec%omega,Frec%W_freq, &
                        Frec%pif0,Frec%pif1, &
                        Atmo%T(o),Atmo%ne(o),iph,ith, &
                        data1O(:,3),rLineO,rPhotO,data2O, &
                        J00(:,o),J00S(:,o),J00P(:,:,o), &
                        J00C_n(:,o),LambdaL(:,:,o), &
                        LambdaP(:,:,:,o),lALI,ALIp,p_exu)

          end if ! MPI/serial

        enddo ! azimuthal angles
      enddo ! polar angles

      ! Error
      if (laborted) goto 2000

      ! If normal MPI
      if (pid.gt.0.and..not.MPID%alternI) then

        !
        ! Send to master
        !

        ! Send indexes
        do while (.True.)
          call MPI_SEND(pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                      ierr)
          if (ierr.eq.0) exit
        end do

        ! Send Stokes
        do while (.True.)
          call MPI_SEND(Stokes_s(if0,Rz0,1,1), &
                        MPID%sizei4(pid), &
                        MPI_DOUBLE_PRECISION, 0, pid, &
                        MPI_COMM_RT, ierr)
          if (ierr.eq.0) exit
        end do

        ! Send profiles
        do while (.True.)
          call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                        MPID%sizei5(pid), &
                        MPI_DOUBLE_PRECISION, &
                        0, 1+pid, MPI_COMM_RT, &
                        ierr)
          if (ierr.eq.0) exit
        end do

        ! If ALI
        if (lALI) then

          ! Send b-b Lambda operator
          do while (.True.)
            call MPI_SEND(rLine_s(1,Rz0,1,1), &
                          MPID%sizei9(pid), &
                          MPI_DOUBLE_PRECISION, &
                          0, 2+pid, MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

          ! Send b-f Lambda operator
          if (ALIp) then
            do while (.True.)
              call MPI_SEND(rPhot_s(1,Rz0,1,1), &
                            MPID%sizei0(pid), &
                            MPI_DOUBLE_PRECISION, &
                            0, 3+pid, MPI_COMM_RT, &
                            ierr)
              if (ierr.eq.0) exit
            end do
          end if

        end if ! ALI
      end if ! Normal MPI

      ! Nullify pointers
2000  nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP)
      nullify(p_StkM,p_StkO,p_LO)

      ! Manage errors in MPI
      if (laborted.and.pid.gt.0) then

        ! Aternative MPI
        if (MPID%alternI) then

          ! Send error
          info_c = (/ -pid, ith, iph /)
          do while (.True.)
            call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                          MPI_COMM_RT,ierr)
            if (ierr.eq.0) exit
          end do

        ! Normal MPI
        else

          ! Send error
          do while (.True.)
            call MPI_SEND(-pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

        end if ! Type of MPI
      end if ! Error and doing MPI

      ! Go back
      return

      end subroutine solveI_RT

!#####################################################################
!#####################################################################
!#####################################################################

      !> Advance the populations by solving the SEE and applying NG
      !! acceleration if requested\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!    Rho_old(Rhoc_class(:)): Structure to store the density
      !!                            matrix of the previous
      !!                            iteration\n
      !!          Atmo(Atmo_class): Structure with atmospheric
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!            MRC(MRC_class): Structure with MRC data\n
      !!             lALI(logical): If doing ALI\n
      !!             PRDl(logical): If doing PRD\n
      !!              ADD(logical): If Stokes needed for PRD-AD\n
      !!              ADT(logical): If initializing PRD-AD with
      !!                            PRD-AA\n
      !!              NGI(logical): If doing NG acceleration\n
      !!             doNG(logical): If doing NG iteration\n
      !!            goout(logical): If converged\n
      !!             iter(integer): Current interation\n
      !!           NG_dim(integer): Size of NG entry\n
      !!         NG_entry(integer): Index of NG entries\n
      !!   NG_scracth(double(:,:)): Data for NG iteration\n
      !!    LambdaL(double(:,:,:)): Lambda operator for lines\n
      !!  LambdaP(double(:,:,:,:)): Lambda operator for
      !!                            photoionizations\n
      !!   Stokes(double(:,:,:,:)): Intensity\n
      !!          J00(double(:,:)): Mean intensity integrated over the
      !!                            absorption profile\n
      !!         J00S(double(:,:)): Mean intensity integrated over the
      !!                            emission profile\n
      !!         J00C(double(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!       J00P(double(:,:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine solveI_SEE(Atom,Rho_old,Atmo,Geom,Input,MRC,lALI, &
                            PRDl,ADD,ADT,NGI,doNG,goout,iter,NG_dim, &
                            NG_entry,NG_scratch,LambdaL,LambdaP, &
                            Stokes,J00,J00S,J00C,J00P)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Rhoc_class), dimension(:), intent(in):: Rho_old
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Input_class), intent(in):: Input
      type(MRC_class), intent(inout):: MRC
      logical, intent(in):: lALI,ADD,ADT,NGI,PRDl
      logical, intent(out):: goout
      logical, intent(inout):: doNG
      integer, intent(in):: iter,NG_dim
      integer, intent(inout):: NG_entry
      double precision, dimension(nxb,nxt,Rz0:Rz1), &
                        intent(in):: LambdaL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1), &
                        intent(in):: LambdaP
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target, intent(inout):: Stokes
      double precision, dimension(nxt,Rz0:Rz1), intent(in):: J00
      double precision, dimension(nxt,Rz0:Rz1), intent(in):: J00S
      double precision, dimension(nxphot,2,Rz0:Rz1), intent(in):: J00P
      double precision, dimension(nfreq,Rz0:Rz1), intent(inout):: J00C
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: NG_scratch

      ! Local

      integer:: ia,itran,jtran,fftran,jftran,iz,m,o,p
      integer:: iterm,iJ,ith,iph,ifreq,ing,ios

      double precision:: daux


#ifdef DEBUGJ00
      if (pid.eq.0) call dump_j00(Atom,J00,J00S,J00P, &
                                  Input%folder,iter)
#endif
#ifdef DEBUGLAMBDA
      if (pid.eq.0) call dump_lambda(Atom,LambdaL,LambdaP, &
                                     Input%folder,iter)
#endif

      ! For each atom
      do ia=1,nA

        ! Limiting indexes
        itran = Atom(ia)%tfshift + 1
        jtran = itran + Atom(ia)%nftran - 1
        fftran = Atom(ia)%pshift + 1
        jftran = fftran + Atom(ia)%nphot - 1

        ! For each height
        do iz=Rz0,Rz1

          ! Solve the SEE
#ifdef DEBUGSEE
          call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                   !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                    J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                    LambdaL(:,itran:jtran,iz), &
                    LambdaP(:,fftran:jftran,:,iz),iz,lALI, &
                    Input%ALI_photo,Input%ALI_allow_off,Input)
#else
          call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                   !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                    J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                    LambdaL(:,itran:jtran,iz), &
                    LambdaP(:,fftran:jftran,:,iz),iz,lALI, &
                    Input%ALI_photo,INPUT%ALI_allow_off)
#endif

        end do ! heights
      end do ! atoms
#ifdef DEBUGRHO00
      if (pid.eq.0) call dump_rho00(Atom,Input%folder,iter)
#endif

      ! Control
      call control
      if (laborted) return

      !
      ! NG acceleration
      !

      ! Check if doing NG acceleration
      if(NGI.and.iter.gt.Input%NGI_delay)then

        ! Advance ntry. The master does not need it if
        ! if accelerated the intensity
        NG_entry = NG_entry + 1

        ! If Master
        if (pid.eq.0) then

          ! Initialize index
          o = 0

          ! For each atom
          do ia=1,NA

            ! For each term
            do iterm=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(iterm)

                ! Get level index
                p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                ! For each height
                do iz=Rz0,Rz1

                  o = o + 1

                  ! Store rho
                  NG_scratch(o,NG_entry) = dble(Atom(ia)%crho(p,iz))

                enddo ! Heights
              enddo ! Levels
            enddo ! Term
          enddo ! Atoms

          ! Save last atomic index
          p = o

          ! If PRD
          if (PRDl) then

            ! If ADD
            if (ADD) then

              ! For each height
              do iz=Rz0,Rz1

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! For each azimuthal direction
                  do iph=1,Geom%nPh

                    ! For each frequency
                    do ifreq=1,nfreq

                      ! Advance index
                      o = o + 1

                      ! Store Stokes
                      NG_scratch(o,NG_entry) = &
                                            Stokes(ifreq,iph,ith,iz)

                    end do ! Frequencies
                  end do ! Azimuthal directions
                end do ! Polar directions
              end do ! Heights

            ! If AV
            else

              ! For each height
              do iz=Rz0,Rz1

                ! For each frequency
                do ifreq=1,nfreq

                  ! Advance index
                  o = o + 1

                  ! Store Stokes
                  NG_scratch(o,NG_entry) = J00C(ifreq,iz)

                end do ! Frequencies
              end do ! Heights

            end if ! AV or AD or dyn
          end if ! PRD

          ! Call NG and check if it should be processed
          call NG(NG_dim,p,Input%NGI_ord,NG_scratch,NG_entry,doNG)

        ! Slave
        else

          ! If wrong order
          if (Input%NGI_ord.lt.1.or.Input%NGI_ord.gt.5) then

            ! Do not do
            doNG = .False.

          ! Valid order
          else

            ! Check if Master is in NG step
            if (NG_entry.gt.(Input%NGI_ord+1)) then

              ! Do step
              doNG = .True.

            ! Not a NG step
            else

              ! Do not do
              doNG = .False.

            end if ! NG step
          end if ! order

        end if ! Master of slave

        ! If communication is needed
        if (doNG) then

          ! If Master, send last
          if (pid.eq.0) then
            ing = NG_entry
          ! Slaves send the one they have
          else
            ing = 1
          end if

          ! MPI
          if (nproc.gt.1) then

            ! Share NG iteration
            call MPI_BCAST(NG_scratch(1,ing), NG_dim, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          end if ! MPI

          ! Reconstruct NG data

          ! Initialize index
          o = 0

          ! For each atom
          do ia=1,NA

            ! For each term
            do iterm=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(iterm)

                ! Get level index
                daux = sqrt(2d0*Atom(ia)%rJval(iJ,iterm) + 1d0)
                m = Atom(ia)%irho(iterm)%irho_ij(iJ)
                p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                ! For each height
                do iz=Rz0,Rz1

                  ! Advance index
                  o = o + 1

                  ! Accelerate rho
                  Atom(ia)%crho(p,iz) = dcmplx(NG_scratch(o,ing), 0d0)
                  Atom(ia)%popu(m,iz) = NG_scratch(o,ing)*daux

                enddo ! Heights
              enddo ! Levels
            enddo ! Terms
          enddo ! Atoms

          ! If PRD
          if (PRDl) then

            ! If ADD
            if (ADD) then

              ! For each height
              do iz=Rz0,Rz1

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! For each azimuthal direction
                  do iph=1,Geom%nPh

                    ! For each frequency
                    do ifreq=1,nfreq

                      ! Advance index
                      o = o + 1

                      ! Store Stokes
                      Stokes(ifreq,iph,ith,iz) = NG_scratch(o,ing)

                    end do ! Frequencies
                  end do ! Azimuthal directions
                end do ! Polar directions
              end do ! Heights

            ! If AV
            else

              ! For each height
              do iz=Rz0,Rz1

                ! For each frequency
                do ifreq=1,nfreq

                  ! Advance index
                  o = o + 1

                  ! Store Stokes
                  J00C(ifreq,iz) = NG_scratch(o,ing)

                end do ! Frequencies
              end do ! Heights

            end if ! AV or AD or dyn
          end if ! PRD

          ! Verbose
          if(gpid.eq.0) then
            umsg = 'NG acceleration'
            call verbose
          end if

          ! Reset entry
          NG_entry = 0

        end if ! Communication
      endif ! NG acceleration


      !
      !  Calculate MRC
      !

      ! Only the master does
      if (pid.eq.0) then

        ! Calculate MRC
        call MRCI_sb(Atom,Rho_old,MRC)

        ! Convert cm into km
        MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

        ! Check exit criteria
        if (MRC%values(2,1).le.Input%mrci_i) goout = .True.

        ! Check in first iteration if was a fix population case
        ! but we want PRD
        if (iter.eq.Input%iteri_min.and.PRDl.and. &
          MRC%values(2,1).lt.1d-16) goout = .False.
        ! If going to AD from AA
        if (iter.eq.Input%iteri_min.and.tbAD.and.ADT) &
          goout = .False.
        ! If doing PRD and just had NG, do not leave
        if (NGI.and.PRDl) then
          if (doNG) goout = .False.
        end if

      end if

      ! Only the global Master does
      if (gpid.eq.0) then

        ! Write in stdout
        write(umsg,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                   '2x,f9.3)') &
        iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
        MRC%indexes(2,1),MRC%values(1,1)
        call verbose

        ! File
        if (Input%keep_MRC) then

          ! Write in MRC file
          open(800, file=trim(Input%folder)//'/MRCI', &
               iostat=ios,err=1000,position='append')
          write(800,'(4x,i9,2x,es20.12,2x,i9,3x,i9,2x,i11,'// &
                    '2x,f9.3)', err=1100) &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(2,1),MRC%values(1,1)
          close(800)

        end if ! Output to file
      end if ! Global master

      return

1000  umsg = 'Error opening MRCI file'
      urou = 'solveI_SEE'
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
1100  umsg = 'Error writing MRCI file'
      urou = 'solveI'
      close(800)
      if (nproc.eq.1) then
        call aborted
      else
        call abortedS(umsg,urou,.True.,.True.)
        call control
      end if
      return

      ! Deceive compiler
      daux = J00S(1,1)

      end subroutine solveI_SEE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Formal solution for the RTE for intensity\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!     Stokes(double(:,:,:,:)): Intensity\n
      !!            J00(double(:,:)): Mean intensity integrated over
      !!                              the absorption profile\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence\n
      !!      SolF(Solution_F_class): Structure with the solution of
      !!                              the self-consistent problem and
      !!                              the corresponding emergent
      !!                              profiles, contribution function,
      !!                              and height for optical depth
      !!                              equal to one
      subroutine emergentI(Atom,LTElines,Atmo,Cont,Frec,Red,Geom, &
                           MPID,Input,Stokes,J00,J00C,SolF)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      type(Geometry_class), intent(inout):: Geom
      type(Solution_F_class), intent(inout):: SolF
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target, intent(in):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: J00C
      double precision, dimension(nxt,Rz0:Rz1), intent(in):: J00

      ! Local

      logical:: doRT,doMG,ADD,l1

      integer:: icount,ncount,ith,iph,if0,if1,tau1size

      double precision, dimension(:,:), allocatable:: ContrG
      double precision, dimension(:,:), allocatable, target:: tau
      double precision, dimension(:,:), allocatable, target:: tau1

      ! Receivers

      double precision, dimension(:), allocatable:: Stokes_r
      double precision, dimension(:), allocatable:: Contr_r

      ! Senders

      double precision, dimension(:), allocatable:: Stokes_s
      double precision, dimension(:,:), allocatable:: Contr_s
      double precision, dimension(:,:), allocatable:: tau1_s

      ! Pointers

      double precision, dimension(:), pointer:: etaIM
      double precision, dimension(:,:), pointer:: data1M,data1O

      ! Nullify pointers
      nullify(data1M,data1O,etaIM)

      ! Initialize arrays
      call emergentI_init(Geom,MPID,Input,SolF,ADD,tau1size, &
                          if0,if1,Stokes_r,Contr_r,Stokes_s,Contr_s, &
                          tau1,ContrG,tau,tau1_s,etaIM,data1M, &
                          data1O)

      ! Initialize PRD
      if (PRD) call initialize_emissI(Atom,Atmo,Geom,Red)

      ! Control
      call control
      if (laborted) goto 2000

      ! Signal if not a master
      doRT = pid.gt.0.or.nproc.eq.1
      doMG = .not.doRT

      ! Reset progress counter
      icount = 0

      ! Determine number of directions to do
      ncount = Geom%nThLOS*Geom%nPhLOS

      !
      ! Formal solutions
      !

      ! For each polar LOS direction
      do ith=1,Geom%nThLOS

        ! For each azimuthal LOS direction
        do iph=1,Geom%nPhLOS

          ! Communicate which direction we are doing if global Master
          if (gpid.eq.0) then
            icount = icount + 1
            write(umsg,'(A,i4,A,i4)') &
                       '   Doing direction ',icount,' of ',ncount
            call verbose
          end if

          ! If PRD
          if (PRD) then

            ! If AD, get geometry
            if (.not.AVI) call get_scattering_los(Geom,ith,iph)

            ! Compute emiss2ord
            call comoving_emissI2ord(Atom,Atmo,Geom,Frec,Red, &
                                     Stokes,J00,J00C,ith,iph, &
                                     Input%PRD_int_mode,.True.)
          end if ! PRD

          ! If dynamic, normalize for this LOS
          if (dyn) &
            call normalizeI(Atom,LTElines,Atmo,Geom,MPID,Frec,Red, &
                            Input%folder,l1,ith,iph,.True.)

          ! If manager
          if (doMG) then

            ! Call manager task
            call emergentI_manager(Atmo,Frec,Geom,MPID,Input,SolF, &
                                   ith,iph,Stokes_r,Contr_r,tau1, &
                                   ContrG)

          ! Slave or serial
          else

            ! Call RT
            call emergentI_RT(Atom,LTElines,Atmo,Cont,Frec,Red, &
                              Geom,MPID,Input,SolF,ADD,tau1size, &
                              if0,if1,ith,iph,Contr_s,Stokes_s, &
                              tau1,tau,tau1_s,etaIM,data1M,data1O, &
                              J00C)
          end if

          ! If dynamic, free LOS norms
          if (dyn) call free_norm(Red,.False.)

          ! Failure
          if (laborted) goto 2000

        end do ! Azimuth
      end do ! Polar

      ! Control
      call control

      !
      ! Clean pointers
      !
2000  if (associated(data1M)) then
        deallocate(data1M,data1O)
        nullify(data1M,data1O)
      end if
      if (associated(etaIM)) then
        deallocate(etaIM)
        nullify(etaIM)
      end if

      ! Free PRD
      if (PRD) call free_e2ord(Red)

      end subroutine emergentI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the approximate RAM space necessary for the
      !! routines solving the formal solution for the RTE for
      !! intensity\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!       MPID(MPI_class): Structure with MPI data\n
      !!    Input(Input_class): Structure with configuration data
      subroutine emergentI_predict(Atom,Atmo,Red,Geom,MPID,Input)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Red_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input

      ! Local

      integer:: iaux,mfreq


      ! Initialize
      TRAMc = 0d0

      ! Predict PRD
      if (PRD) call predict_emissI(Atom,Atmo,Geom,Red)

      !
      ! Allocations
      !

      ! MPI
      if (MPID%mpi) then

        ! CPU limits
        mfreq = MPID%if1(pid) - MPID%if0(pid) + 1

        ! Master
        if (pid.eq.0) then

          ! To receive Intensity chunks
          iaux = MPID%nxfreq
          TRAMc = TRAMc + 8d-6*dble(iaux)

          ! If calculating height of tau=1, allocate
          if (Input%out_tau1) &
            TRAMc = TRAMc + 8d-6*dble(2*nfreq)

          ! If calculating contribution function, allocate
          if (Input%out_contr) then
            iaux = MPID%nxfreq*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux + nfreq*Rnz)
          end if

          ! If inverting, need to return the output
          if (run_mode.eq.-1) then

            ! Return Stokes parameters
            SRAMc = SRAMc + 8d-6*dble(nfreq*Geom%nPhLOS* &
                                            Geom%nThLOS)

            ! Return height where optical depth is one
            if (Input%out_tau1) &
              SRAMc = SRAMc + 4d-6*dble(Input%lim_tau%nn* &
                                        Geom%nPhLOS* &
                                        Geom%nThLOS)

            ! Return contribution function
            if (Input%out_contr) &
              SRAMc = SRAMc + 4d-6*dble(Input%lim_ctr%nn* &
                                        Geom%nPhLOS* &
                                        Geom%nThLOS*Rnz)

          end if ! Inversion

        ! Slave
        else

          ! Allocate M and O pointers for RT coeff
          TRAMc = TRAMc + 8d-6*dble(mfreq*6)

          ! To send Intensity chunks
          TRAMc = TRAMc + 8d-6*mfreq

          ! If calculating tau 1 or contribution function, allocate
          if (input%out_tau1.or.input%out_contr) then

            TRAMc = TRAMc + 8d-6*dble(mfreq*(5 + Rnz))

            ! If calculating contribution function, allocate
            if (input%out_contr) &
              TRAMc = TRAMc + 8d-6*dble(mfreq*(1 + Rnz))

          end if ! tau1 output

        end if ! Master/slave

      ! Serial
      else

        ! Allocate M and O pointers for RT coeff
        TRAMc = TRAMc + 8d-6*dble(6*nfreq)

        ! If calculating height of tau=1, allocate
        if (Input%out_tau1.or.Input%out_contr) &
          TRAMc = TRAMc + 8d-6*dble(nfreq*(3 + 2*Rnz))

        ! If inverting, need to return the output
        if (run_mode.eq.-1) then

          ! Returning Stokes output
          SRAMc = SRAMc + 8d-6*dble(nfreq*Geom%nPhLOS* &
                                          Geom%nThLOS)

          ! If output height tau equal 1
          if (Input%out_tau1) &
            SRAMc = SRAMc + 4d-6*dble(Geom%nPhLOS* &
                                      Geom%nThLOS* &
                                      Input%lim_tau%nn)

          ! If output contribution function
          if (Input%out_contr) &
            SRAMc = SRAMc + 4d-6*dble(Geom%nPhLOS* &
                                      Geom%nThLOS* &
                                      Input%lim_ctr%nn*Rnz)
        end if ! Inversion
      end if ! Serial

      end subroutine emergentI_predict

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize arrays and pointers for the formal solution for
      !! the RTE for intensity\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!  SolF(Solution_F_class): Structure with the solution of the
      !!                          self-consistent problem and the
      !!                          corresponding emergent profiles,
      !!                          contribution function, and height
      !!                          for optical depth equal to one\n
      !!            ADD(logical): If Stokes needed for PRD-AD\n
      !!       tau1size(integer): Size sender tau1 package\n
      !!            if0(integer): Initial frequency index\n
      !!            if1(integer): Final frequency index\n
      !!     Stokes_r(double(:)): Receiver Stokes buffer\n
      !!      Contr_r(double(:)): Receiver contribution function
      !!                          buffer\n
      !!     Stokes_s(double(:)): Sender Stokes buffer\n
      !!    Contr_s(double(:,:)): Sender contribution function
      !!                          buffer\n
      !!       tau1(double(:,:)): Data height tau equal to 1\n
      !!     ContrG(double(:,:)): Data contribution function\n
      !!        tau(double(:,:)): Data optical depth\n
      !!     tau1_s(double(:,:)): Height optical depth equal to
      !!                          one sender buffer\n
      !!        etaIM(double(:)): Absorptivity point M\n
      !!     data1M(double(:,:)): RT coeff point M\n
      !!     data1O(double(:,:)): RT coeff point O
      subroutine emergentI_init(Geom,MPID,Input,SolF,ADD,tau1size, &
                                if0,if1,Stokes_r,Contr_r,Stokes_s, &
                                Contr_s,tau1,ContrG,tau,tau1_s, &
                                etaIM,data1M,data1O)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(out):: ADD
      integer, intent(out):: tau1size,if0,if1
      double precision, dimension(:), &
                        allocatable, intent(out):: Stokes_r
      double precision, dimension(:), &
                        allocatable, intent(out):: Contr_r
      double precision, dimension(:), &
                        allocatable, intent(out):: Stokes_s
      double precision, dimension(:,:), &
                        allocatable, intent(out):: Contr_s
      double precision, dimension(:,:), &
                        allocatable, intent(out):: tau1
      double precision, dimension(:,:), &
                        allocatable, intent(out):: ContrG
      double precision, dimension(:,:), &
                        allocatable, intent(out):: tau
      double precision, dimension(:,:), &
                        allocatable, intent(out):: tau1_s
      double precision, dimension(:), pointer, intent(inout):: etaIM
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data1M,data1O

      ! Local

      logical:: AD

      integer:: iaux


      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      !
      ! Allocations
      !

      ! MPI
      if (MPID%mpi) then

        ! CPU limits
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! Master
        if (pid.eq.0) then

          ! To receive Intensity chunks
          iaux = MPID%nxfreq
          allocate(Stokes_r(iaux))

          ! If calculating height of tau=1, allocate
          if (Input%out_tau1) allocate(tau1(2,nfreq))

          ! If calculating contribution function, allocate
          if (Input%out_contr) then
            iaux = MPID%nxfreq*Rnz
            allocate(Contr_r(iaux))
            allocate(ContrG(nfreq,Rz0:Rz1))
          end if

          ! If inverting, need to return the output
          if (run_mode.eq.-1) then
            allocate(SolF%e_Stk(0:0,nfreq,Geom%nPhLOS,Geom%nThLOS))
            if (Input%out_tau1) &
            allocate(SolF%e_tau1(Input%lim_tau%nn,Geom%nPhLOS, &
                                 Geom%nThLOS))
            if (Input%out_contr) &
            allocate(SolF%e_Ctr(0:0,Input%lim_ctr%nn,Rz0:Rz1, &
                                Geom%nPhLOS,Geom%nThLOS))
          end if ! Inversion

        ! Slave
        else

          ! Allocate M and O pointers for RT coeff
          allocate(data1M(MPID%nf(pid),3))
          allocate(data1O(MPID%nf(pid),3))

          ! To send Intensity chunks
          allocate(Stokes_s(MPID%if0(pid):MPID%if1(pid)))

          ! If calculating tau 1 or contribution function, allocate
          if (input%out_tau1.or.input%out_contr) then

            ! Allocate arrays
            allocate(tau(if0:if1,Rz0:Rz1))
            allocate(tau1(2,if0:if1))
            allocate(tau1_s(2,if0:if1))
            tau1size = MPID%nf(pid)*2
            allocate(etaIM(MPID%nf(pid)))

            ! If calculating contribution function, allocate
            if (input%out_contr) &
              allocate(Contr_s(1:MPID%nf(pid),Rz0:Rz1))

          ! No calculating tau or contribution function
          else

            ! Nullify pointer
            nullify(etaIM)

          end if ! tau1 output

        end if ! Master/slave

      ! Serial
      else

        ! CPU limits
        if0 = 1
        if1 = nfreq

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(nfreq,3))
        allocate(data1O(nfreq,3))

        ! If calculating height of tau=1, allocate
        if (Input%out_tau1.or.Input%out_contr) then
          allocate(tau(nfreq,Rz0:Rz1))
          allocate(etaIM(nfreq))
          allocate(tau1(2,nfreq))
          allocate(Contr_s(nfreq,Rz0:Rz1))
        ! Not allocating, nullify
        else
          nullify(etaIM)
        end if

        ! If inverting, need to return the output
        if (run_mode.eq.-1) then
          allocate(SolF%e_Stk(0:0,nfreq,Geom%nPhLOS,Geom%nThLOS))
          if (Input%out_tau1) &
          allocate(SolF%e_tau1(Input%lim_tau%nn, &
                               Geom%nPhLOS,Geom%nThLOS))
          if (Input%out_contr) &
          allocate(SolF%e_Ctr(0:0,Input%lim_ctr%nn,Rz0:Rz1, &
                              Geom%nPhLOS,Geom%nThLOS))
        end if ! Inversion

      end if ! MPI/serial

      end subroutine emergentI_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Parallel-master task for the formal solution for the RTE for
      !! intensity\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Frec(Frequency_class): Structure with frequency data\n
      !!    Geom(Geometry_class): Structure with geometric data\n
      !!         MPID(MPI_class): Structure with MPI data\n
      !!      Input(Input_class): Structure with configuration data\n
      !!  SolF(Solution_F_class): Structure with the solution of the
      !!                          self-consistent problem and the
      !!                          corresponding emergent profiles,
      !!                          contribution function, and height
      !!                          for optical depth equal to one\n
      !!            ith(integer): Index polar direction\n
      !!            iph(integer): Index azimutal direction\n
      !!     Stokes_r(double(:)): Receiver Stokes buffer\n
      !!      Contr_r(double(:)): Receiver contribution function
      !!                          buffer\n
      !!       tau1(double(:,:)): Data height tau equal to 1\n
      !!     ContrG(double(:,:)): Data contribution function
      subroutine emergentI_manager(Atmo,Frec,Geom,MPID,Input,SolF, &
                                   ith,iph,Stokes_r,Contr_r,tau1, &
                                   ContrG)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      type(Solution_F_class), intent(inout):: SolF
      integer, intent(in):: iph,ith
      double precision, dimension(:), &
                        allocatable, intent(inout):: Stokes_r
      double precision, dimension(:), &
                        allocatable, intent(inout):: Contr_r
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: tau1
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: ContrG

      ! Local

      integer:: id,rpid,ierr,tau1size,sshift,iz
      integer, dimension(2):: info_b

      double precision, dimension(nfreq):: Stokes_out


      ! If calculating height of tau=1
      if (Input%out_tau1) then

        ! For each range
        do id=1,MPID%nnd

          ! Receive ID
          call MPI_recv(rpid,1,MPI_INTEGER, &
                        MPI_ANY_SOURCE, 2, &
                        MPI_COMM_RT, MPI_STATUS_IGNORE, &
                        ierr)

          ! Get package size
          tau1size = MPID%nf(rpid)*2

          ! Receive height where optical depth is one
          call MPI_recv(tau1(1,MPID%if0(rpid)), &
                        tau1size, &
                        MPI_DOUBLE_PRECISION, rpid, &
                        3+rpid, MPI_COMM_RT, &
                        MPI_STATUS_IGNORE, ierr)

        end do ! Ranges

        ! Inversion
        if (run_mode.eq.-1) then

          ! Keep tau1
          call settau(SolF%e_tau1(:,iph,ith),tau1(2,:), &
                      Input%lim_tau)

        ! Synthesis
        else

          ! Store the height where tau=1
          call writetau(Input%folder,iph,ith,Frec%omega,Geom, &
                        tau1(2,:),Input%lim_tau)
          call control
          if (laborted) return

        end if ! Inversion
      end if ! calculate tau=1

      ! For each frequency domain
      do id=1,MPID%nnd

        !
        ! Receive data from a slave
        !

        ! Receive indexing data
        call MPI_recv(info_b,2,MPI_INTEGER, &
                      MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                      MPI_STATUS_IGNORE, ierr)

        ! If calculating contribution function
        if (Input%out_contr) then

          ! Receive contribution function data
          call MPI_recv(Contr_r(1), MPID%sizei14(info_b(1)), &
                        MPI_DOUBLE_PRECISION, info_b(1), &
                        1+info_b(1), MPI_COMM_RT, &
                        MPI_STATUS_IGNORE, ierr)

          ! Reset shift in index
          sshift = 0

          ! For each height
          do iz=Rz0,Rz1

            ! Rearrange the contribution function
            ContrG(MPID%if0(info_b(1)):MPID%if1(info_b(1)), &
                   iz) = Contr_r(sshift+1: &
                                 sshift+MPID%nf(info_b(1)))

            ! Update the shift in the buffer
            sshift = sshift + MPID%nf(info_b(1))

          end do ! heights

        end if ! If contribution function

        ! If it is not a boundary, we are not receiving
        ! anything else
        if (info_b(2).lt.0) cycle

        ! Receive intensity
        call MPI_recv(Stokes_r(1), MPID%nf(info_b(1)), &
                      MPI_DOUBLE_PRECISION, info_b(1), &
                      info_b(1), MPI_COMM_RT, &
                      MPI_STATUS_IGNORE, ierr)

        ! Rearrange the intensity
        Stokes_out(MPID%if0(info_b(1)):MPID%if1(info_b(1))) = &
                                 Stokes_r(1:MPID%nf(info_b(1)))

      end do ! Frequency domains

      !
      ! Output
      !

      ! If inverting
      if (run_mode.eq.-1) then

        ! Keep Stokes
        call setstkI(SolF%e_Stk(:,:,iph,ith),Stokes_out, &
                     Input%lim_stk,.False.)

        ! Keep contribution function
        if (Input%out_contr) &
          call setctrI(SolF%e_Ctr(:,:,:,iph,ith),ContrG, &
                       Input%lim_ctr)

      ! Synthesis, to file
      else

        ! Write stokes
        call writestkI(Input%folder,iph,ith,Frec%omega,Geom, &
                       Stokes_out,Input%lim_stk)
        if (laborted) return

        ! Write contribution function
        if (Input%out_contr) then
          call writectrI(Input%folder,iph,ith,Frec%omega,Geom, &
                         Atmo%z,ContrG,Input%lim_ctr)
          if (laborted) return
        end if

      end if ! Inversion or synthesis

      ! If stdout terminal, remove bar and say completed
      if (gpid.eq.0) then
        umsg = '   Completed'
        call verbose
      end if

      end subroutine emergentI_manager

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the RTE for intensity for a given LOS\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!      SolF(Solution_F_class): Structure with the solution of
      !!                              the self-consistent problem and
      !!                              the corresponding emergent
      !!                              profiles, contribution function,
      !!                              and height for optical depth
      !!                              equal to one
      !!                ADD(logical): If Stokes needed for PRD-AD\n
      !!           tau1size(integer): Size sender tau1 package\n
      !!                if0(integer): Initial frequency index\n
      !!                if1(integer): Final frequency index\n
      !!                ith(integer): Index polar direction\n
      !!                iph(integer): Index azimutal direction\n
      !!        Contr_s(double(:,:)): Sender contribution function
      !!                              buffer\n
      !!         Stokes_s(double(:)): Sender Stokes buffer\n
      !!           tau1(double(:,:)): Data height tau equal to 1\n
      !!            tau(double(:,:)): Data optical depth\n
      !!         tau1_s(double(:,:)): Height optical depth equal to
      !!                              one sender buffer\n
      !!            etaIM(double(:)): Absorptivity point M\n
      !!         data1M(double(:,:)): RT coeff point M\n
      !!         data1O(double(:,:)): RT coeff point O\n
      !!     Stokes(double(:,:,:,:)): Intensity\n
      !!            J00(double(:,:)): Mean intensity integrated over
      !!                              the absorption profile\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence
      subroutine emergentI_RT(Atom,LTElines,Atmo,Cont,Frec,Red, &
                              Geom,MPID,Input,SolF,ADD,tau1size, &
                              if0,if1,ith,iph,Contr_s,Stokes_s, &
                              tau1,tau,tau1_s,etaIM,data1M,data1O, &
                              J00C)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Geometry_class), intent(inout):: Geom
      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(in):: Input
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: ADD
      integer, intent(in):: tau1size,ith,iph,if0,if1
      double precision, dimension(:), &
                        allocatable, intent(inout):: Stokes_s
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: Contr_s
      double precision, dimension(:,:), &
                        allocatable, target, intent(inout):: tau
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: tau1
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: tau1_s
      double precision, dimension(:), &
                        pointer, intent(inout):: etaIM
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: J00C

      ! Local

      integer:: iz,iz0,iz1,diz,m,o,p,op,mfreq,jdir
      integer, dimension(2):: info_b

      double precision:: mu_inv,dsm,dsp,dzm,dzp,ct,st,cc,sc,vfac

      ! Dummies
      double precision, dimension(1):: ad1
      double precision, dimension(1,1):: ad2

      ! Pointers
      double precision, dimension(:), pointer:: p_K0M, p_SM, p_StkM
      double precision, dimension(:), pointer:: p_K0O, p_SO, p_StkO
      double precision, dimension(:), pointer:: p_K0P, p_SP
      double precision, dimension(:), pointer:: tauM
      double precision, dimension(:), pointer:: etaIO
      double precision, dimension(:,:), pointer:: data1P


      ! Nullify pointers
      nullify(p_K0M,p_SM,p_StkM,p_K0O,p_SP,p_StkO,p_K0P,p_SP)
      nullify(tauM,etaIO,data1P)

      ! Get direction index
      jdir = Geom%i_geom(iph,ith)

      ! Calculate inverse of cosine of polar direction
      mu_inv = 1d0/Geom%L_mu(ith)

      ! Determine the direction of propagation for indexes
      diz = -int(sign(1d0, Geom%L_mu(ith)))

      ! Determine the first and last height indexes to run over
      iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
      iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

      ! Initialize index of Stokes
      op = 1

      ! Number of frequencies
      mfreq = if1 - if0 + 1

      ! Trigonometry for Doppler shift
      if (dyn) then
        ct = Geom%L_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = cos(Geom%L_phi(iph))
        sc = sin(Geom%L_phi(iph))
      else
        vfac = 1d0
      end if


      !
      ! If calculating height of tau=1
      !
      if (Input%out_tau1.or.Input%out_contr) then

        ! Reset cummulative quantities
        tau = 0d0
        tau1(1,:) = 0
        tau1(2,:) = Atmo%z(Rz0)
        tauM => tau(if0:if1,Rz0)

        ! First height

        ! Top boundary
        o = Rz0

        ! Get Doppler shift
        if (dyn) then

          ! Amplitude
          vfac = sqrt(Atmo%vx(o)*Atmo%vx(o) + &
                      Atmo%vy(o)*Atmo%vy(o) + &
                      Atmo%vz(o)*Atmo%vz(o))

          ! Doppler shift
          if (vfac.gt.TINYVEL) then

            ! Shift
            vfac = 1d0 - Atmo%vx(o)*st*cc - &
                         Atmo%vy(o)*st*sc - &
                         Atmo%vz(o)*ct

          ! Static
          else

            ! No shift
            vfac = 1d0

          end if
        end if ! If dynamic

        ! Calculate absorptivity
        call RTAbsI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                    o,jdir,if0,if1,Cont%ndir, &
                    Cont%c(:,:,:,o),etaIM)

        ! Middle heights
        do iz=Rz0+1,Rz1

          ! Allocate O pointer
          allocate(etaIO(MPID%nf(pid)))

          ! Determine indexes
          m = iz - 1
          o = iz

          ! Calculate distance to previous point
          dsm = abs(Atmo%z(o) - Atmo%z(m))*mu_inv

          ! If tau scale
          if (ztau) &
            dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                Atmo%chi500(m))

          ! Get Doppler shift
          if (dyn) then

            ! Amplitude
            vfac = sqrt(Atmo%vx(o)*Atmo%vx(o) + &
                        Atmo%vy(o)*Atmo%vy(o) + &
                        Atmo%vz(o)*Atmo%vz(o))

            ! Doppler shift
            if (vfac.gt.TINYVEL) then

              ! Shift
              vfac = 1d0 - Atmo%vx(o)*st*cc - &
                           Atmo%vy(o)*st*sc - &
                           Atmo%vz(o)*ct

            ! Static
            else

              ! No shift
              vfac = 1d0

            end if
          end if ! If dynamic

          ! Calculate opacity current point
          call RTAbsI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                      o,jdir,if0,if1,Cont%ndir, &
                      Cont%c(:,:,:,o),etaIO)

          ! Accumulate tau
          call RTtauI(dsm,mfreq,Atmo%z(m),Atmo%z(o), &
                      etaIM,etaIO,tauM,tau(if0:if1,o), &
                      tau1(:,if0:if1))


          ! Shift the opacity value and previous tau
          deallocate(etaIM)
          etaIM => etaIO
          nullify(etaIO)
          tauM => tau(if0:if1,o)

        end do ! Heights

        ! Free tauM
        nullify(tauM)

        !
        ! Send to master the tau=1 if needed
        if (pid.gt.0.and.Input%out_tau1) then

          ! Wait for last send to finish
          call MPI_WAIT(MPID%request7,MPI_STATUS_IGNORE,ierr)
          call MPI_WAIT(MPID%request8,MPI_STATUS_IGNORE,ierr)

          ! Send indexes
          call MPI_ISEND(pid,1,MPI_INTEGER,0,2,MPI_COMM_RT, &
                         MPID%request7,ierr)

          ! Send tau1
          tau1_s = tau1
          call MPI_ISEND(tau1_s(1,if0), tau1size, &
                         MPI_DOUBLE_PRECISION, &
                         0,3+pid,MPI_COMM_RT, &
                         MPID%request8,ierr)

          ! If synthesis, control
          if (run_mode.ne.-1) then
            call control
            if (laborted) goto 2000
          end if

        ! Serial outputting z_tau=1
        else if (Input%out_tau1) then

          ! Inversion
          if (run_mode.eq.-1) then

            ! Keep tau1
            call settau(SolF%e_tau1(:,iph,ith),tau1(2,:), &
                        Input%lim_tau)

          ! Synthesis
          else

            ! Write to file
            call writetau(Input%folder,iph,ith,Frec%omega,Geom, &
                          tau1(2,:),Input%lim_tau)

          end if ! Inversion
        end if ! Slave or serial master
      end if ! Calculate tau=1 for output or tau for contribution

      !
      ! Actual intensity emergence
      !

      !
      ! First height
      !

      ! Get Doppler shift
      if (dyn) then

        ! Amplitude
        vfac = sqrt(Atmo%vx(iz0)*Atmo%vx(iz0) + &
                    Atmo%vy(iz0)*Atmo%vy(iz0) + &
                    Atmo%vz(iz0)*Atmo%vz(iz0))

        ! Doppler shift
        if (vfac.gt.TINYVEL) then

          ! Shift
          vfac = 1d0 - Atmo%vx(iz0)*st*cc - &
                       Atmo%vy(iz0)*st*sc - &
                       Atmo%vz(iz0)*ct

        ! Static
        else

          ! No shift
          vfac = 1d0

        end if
      end if ! If dynamic

      ! If going down, get top boundary
      if(diz.eq.1)then

        ! Call top boundary
        call topI(if0,if1,data1M(:,3))

      ! If going up, get bottom boundary
      else

        ! Call bottom boundary
        call bottomI(Frec%omega,Atmo%T(iz0),vfac,if0,if1, &
                     data1M(:,3))

      endif ! propagation direction

      ! Identify current height
      o = iz0

      ! Index for Stokes
      if (PRD.and.ADD) op = o

      ! Calculate radiative coefficients
      call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                    o,jdir,1,if0,if1,J00C(:,o), &
                    Cont%ndir,Cont%c(:,:,:,o), &
                    ad1,ad1,data1M(:,1:2), &
                    ad2,.False.)

      ! If calculating contribution function
      if (Input%out_contr) then

        ! Salve, wait till last communication was received
        if (pid.gt.0) &
          call MPI_WAIT(MPID%request5,MPI_STATUS_IGNORE,ierr)

        ! Contribution at first point is zero
        Contr_s(:,o) = 0d0

      end if ! computing contribution function

      ! Identify next height
      p = iz0 + diz

      ! Index for Stokes
      if (PRD.and.ADD) op = p

      ! Get Doppler shift
      if (dyn) then

        ! Amplitude
        vfac = sqrt(Atmo%vx(p)*Atmo%vx(p) + &
                    Atmo%vy(p)*Atmo%vy(p) + &
                    Atmo%vz(p)*Atmo%vz(p))

        ! Doppler shift
        if (vfac.gt.TINYVEL) then

          ! Shift
          vfac = 1d0 - Atmo%vx(p)*st*cc - &
                       Atmo%vy(p)*st*sc - &
                       Atmo%vz(p)*ct

        ! Static
        else

          ! No shift
          vfac = 1d0

        end if
      end if ! If dynamic

      ! Calculate radiative coefficients
      call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                    p,jdir,1,if0,if1,J00C(:,p), &
                    Cont%ndir,Cont%c(:,:,:,p), &
                    ad1,ad1,data1O(:,1:2), &
                    ad2,.False.)

      !
      ! Intermediate heights
      !

      ! For each height this CPU has assigned
      do iz=iz0,iz1,diz

        ! We treat the boundaries outside
        if (iz.eq.iz0.or.iz.eq.iz1) cycle

        ! Allocate P pointers
        allocate(data1P(MPID%nf(pid),3))

        ! Identify heights
        m = iz - diz
        o = iz
        p = iz + diz

        ! Calculate distance to previous point
        dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

        ! Calculate distance to the next point
        dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

        ! If tau scale
        if (ztau) then

          ! Get geometrical distance
          dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                              Atmo%chi500(m))
          dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                              Atmo%chi500(p))
        end if

        ! Index for Stokes
        if (PRD.and.ADD) op = p

        ! Get Doppler shift
        if (dyn) then

          ! Amplitude
          vfac = sqrt(Atmo%vx(p)*Atmo%vx(p) + &
                      Atmo%vy(p)*Atmo%vy(p) + &
                      Atmo%vz(p)*Atmo%vz(p))

          ! Doppler shift
          if (vfac.gt.TINYVEL) then

            ! Shift
            vfac = 1d0 - Atmo%vx(p)*st*cc - &
                         Atmo%vy(p)*st*sc - &
                         Atmo%vz(p)*ct

          ! Static
          else

            ! No shift
            vfac = 1d0

          end if
        end if ! If dynamic

        ! Calculate radiative coefficients
        call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                      p,jdir,1,if0,if1,J00C(:,p), &
                      Cont%ndir,Cont%c(:,:,:,p), &
                      ad1,ad1,data1P(:,1:2), &
                      ad2,.False.)

        ! Point to the data
        p_K0M  => data1M(:,1)
        p_SM   => data1M(:,2)
        p_StkM => data1M(:,3)
        p_K0O  => data1O(:,1)
        p_SO   => data1O(:,2)
        p_StkO => data1O(:,3)
        p_K0P  => data1P(:,1)
        p_SP   => data1P(:,2)

        ! Apply short characteristics BESSER
        call RTStepI(o,ith,iph,mfreq,dsm,dsp,p_K0M,p_SM,p_K0O, &
                     p_SO,p_K0P,p_SP,p_StkM, &
                     p_StkO,ad1,.False.,.True.)

        ! If calculating contribution function
        if (Input%out_contr) then

          ! Calculate vertical distance between points
          if (ztau) then
            dzm = dsm/mu_inv
            dzp = dsp/mu_inv
          else
            dzm = Atmo%z(o) - Atmo%z(m)
            dzp = Atmo%z(p) - Atmo%z(o)
          end if

          ! Compute contribution
          call RTContrI(dsm,dsp,mfreq,dzm,dzp, &
                        p_K0M,p_K0O,p_SO,p_K0P, &
                        tau(if0:if1,o),Contr_s(:,o), &
                        .True.)

        end if ! calculating contribution

        ! Shift data (O->M, P->O)
        deallocate(data1M)
        data1M => data1O
        data1O => data1P
        nullify(data1P)

      end do ! Intermediate heights

      !
      ! Last height
      !

      ! Identify heights
      m = iz1 - diz
      o = iz1

      ! Calculate distance to previous point
      dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

      ! If tau scale, get geometrical distance
      if (ztau) &
        dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                            Atmo%chi500(m))

      ! Point to the data
      p_K0M  => data1M(:,1)
      p_SM   => data1M(:,2)
      p_StkM => data1M(:,3)
      p_K0O  => data1O(:,1)
      p_SO   => data1O(:,2)
      p_StkO => data1O(:,3)

      ! Apply short characteristics LINEAR
      call RTStepI(o,ith,iph,mfreq, &
                   dsm,dsp,p_K0M,p_SM,p_K0O, &
                   p_SO,p_K0P,p_SP,p_StkM, &
                   p_StkO,ad1,.False.,.False.)

      ! If calculating contribution function
      if (Input%out_contr) then

        ! Calculate vertical distance between points
        if (ztau) then
          dzm = dsm/mu_inv
        else
          dzm = Atmo%z(o) - Atmo%z(m)
        end if

        ! Calculate contribution function
        call RTContrI(dsm,dsp,mfreq,dzm,dzp,p_K0M, &
                      p_K0O,p_SO,p_K0P,tau(if0:if1,o), &
                      Contr_s(:,o),.False.)

      end if ! calculating contribution function

      ! Serial
      if (pid.eq.0) then

        ! If calculating contribution function
        if (Input%out_contr) then

          ! If inverting
          if (run_mode.eq.-1) then

            ! Keep contribution function
            call setctrI(SolF%e_Ctr(:,:,:,iph,ith),Contr_s, &
                         Input%lim_ctr)

          ! Synthesis
          else

            ! Write contribution function
            call writectrI(Input%folder,iph,ith,Frec%omega,Geom, &
                           Atmo%z,Contr_s,Input%lim_ctr)
          end if
        end if ! calculating contribution function

        ! If inverting
        if (run_mode.eq.-1) then

          ! Keep Stokes
          call setstkI(SolF%e_Stk(:,:,iph,ith),p_StkO, &
                       Input%lim_stk,.False.)

        ! Synthesis
        else

          ! Write stokes
          call writestkI(Input%folder,iph,ith,Frec%omega,Geom, &
                         p_StkO,Input%lim_stk)

        end if

        ! Communicate we finished this direction
        if (gpid.eq.0) then
          umsg = '   Completed'
          call verbose
        end if

      ! Slave
      else

        ! Wait till last communication was received
        call MPI_WAIT(MPID%request3,MPI_STATUS_IGNORE,ierr)
        call MPI_WAIT(MPID%request4,MPI_STATUS_IGNORE,ierr)

        !
        ! Send to master
        !

        ! Send indexes
        info_b = (/ pid, 1 /)
        call MPI_ISEND(info_b(1),2,MPI_INTEGER, &
                       0,0,MPI_COMM_RT,MPID%request3, &
                       ierr)

        ! If calculating contribution function
        if (Input%out_contr) then

          ! Send contribution function
          call MPI_ISEND(contr_s(1,Rz0), &
                         MPID%sizei14(pid), &
                         MPI_DOUBLE_PRECISION, &
                         0, 1+pid, MPI_COMM_RT, &
                         MPID%request5, ierr)
        end if

        ! Send intensity
        Stokes_s = data1O(:,3)
        call MPI_ISEND(Stokes_s(if0), MPID%nf(pid), &
                       MPI_DOUBLE_PRECISION, 0, pid, &
                       MPI_COMM_RT, MPID%request4, ierr)

      end if ! Serial or slave

      ! Free pointers
2000  nullify(p_K0M,p_SM,p_StkM,p_K0O,p_SP,p_StkO,p_K0P,p_SP)

      end subroutine emergentI_RT

!#####################################################################
!#####################################################################
!#####################################################################

      !> Self-consistent solution for the problem for continuum
      !! intensity\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!     Stokes(double(:,:,:,:)): Intensity\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence
      subroutine solveJ(Atmo,Atom,LTElines,Cont,Frec,Red,Geom,MPID, &
                        Input,Stokes,J00C)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        target, intent(inout):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(inout):: J00C

      ! Local

      type(MRC_class):: MRC

      logical:: ADD,PRDl,goout

      integer:: if0,if1,iter,ntpz,npz,nsend

      double precision, dimension(:,:), allocatable:: J00C_n

      ! Receivers

      double precision, dimension(:), allocatable:: Stokes_r

      ! Senders

      double precision, dimension(:,:,:,:), allocatable:: Stokes_s

      ! Pointers

      double precision, dimension(:,:), pointer:: data1M,data1O


      ! Initialization
      call solveJ_init(Geom,MPID,Input,if0,if1,npz,ntpz,nsend,PRDl, &
                       ADD,Stokes_r,Stokes_s,J00C_n,data1M,data1O)

      ! Control
      call control
      if (laborted) goto 2000

      !
      ! Start iterations
      !

      ! For each iteration between the limits specified
      do iter=1,Input%iter_j

        ! Master, reset
        if (pid.eq.0) J00C_n = 0d0

        ! If MPI and master
        if (MPID%mpi.and.pid.eq.0) then

          ! Call manager
          call solveJ_manager(Geom,MPID,npz,ntpz,nsend, &
                              Stokes_r,Stokes,J00C_n)

        ! Serial or slave
        else

          ! Call formal solution
          call solveJ_RT(Atom,Atmo,LTElines,Cont,Frec,Red,Geom,MPID, &
                         Input,if0,if1,data1M,data1O, &
                         Stokes_s,Stokes,J00C,J00C_n)

        end if ! Manage or compute

        ! Control
        call control
        if (laborted) goto 2000

        ! Master compute MRC
        if (pid.eq.0) then

          !
          ! Calculate MRC for J
          !

          ! Call the routine
          call MRCJ_sb(J00C_n,J00C,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Write in stdout
          if (gpid.eq.0) then
            write(umsg,'(3x,"J it:",1x,i3,2x,es20.12,'// &
                       '2x,i9,2x,f10.4,4x,i9,2x,f9.3)') &
                       iter,MRC%values(2,1),MRC%indexes(1,1), &
                       1d2/Frec%omega(MRC%indexes(1,1)), &
                       MRC%indexes(2,1),MRC%values(1,1)
            call verbose
          end if

          ! Check convergence
          goout = MRC%values(2,1).lt.Input%mrcj

        end if ! Master

        ! Master overwrite
        if (pid.eq.0) J00C = J00C_n

        ! If MPI
        if (MPID%mpi) then

          !
          ! Share the radiation information
          !

          ! Share J00C
          call MPI_BCAST(J00C(1,Rz0), MPID%sizei7(0), &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

          ! Share intensity if doing A-D PRD
          if (PRD.and.ADD) &
            call MPI_BCAST(Stokes(1,1,1,giz0), MPID%sizei8(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Finished
          call MPI_BCAST(goout,1,MPI_LOGICAL,0, &
                         MPI_COMM_RT, ierr)

        end if ! MPI

        ! Finish?
        if (goout) exit

      end do ! Iterations

      ! B-B, restore PRD
      if (Input%init_J_bb) PRD = PRDl

      !
      ! Clean slave pointers
      !
2000  if (pid.gt.0.or..not.MPID%mpi) then
        deallocate(data1M,data1O)
        nullify(data1M,data1O)
      end if

      return

      end subroutine solveJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize arrays and pointers for the self-consistent
      !! solution of the problem for continuum intensity\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                              data\n
      !!               if0(integer): Initial frequency index\n
      !!               if1(integer): Final frequency index\n
      !!               npz(integer): Size in z and phi\n
      !!              ntpz(integer): Size in z, theta, and phi\n
      !!             nsend(integer): Expected messages per iteration\n
      !!              PRDl(logical): If originally PRD\n
      !!               ADD(logical): If storing Stokes for PRD-AD\n
      !!        Stokes_r(double(:)): Receiver buffer\n
      !!  Stokes_s(double(:,:,:,:)): Sender buffer\n
      !!        J00C_n(double(:,:)): New value of J00C\n
      !!        data1M(double(:,:)): RT coeff point M\n
      !!        data1O(double(:,:)): RT coeff point O\n
      subroutine solveJ_init(Geom,MPID,Input,if0,if1,npz,ntpz,nsend, &
                             PRDl,ADD,Stokes_r,Stokes_s,J00C_n, &
                             data1M,data1O)
      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      logical, intent(out):: PRDl,ADD
      integer, intent(out):: if0,if1,npz,ntpz,nsend
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, dimension(:), &
                        allocatable, intent(out):: Stokes_r
      double precision, dimension(:,:), allocatable, &
                        intent(out):: J00C_n
      double precision, dimension(:,:,:,:), allocatable, &
                        intent(out):: Stokes_s

      ! Local

      logical:: AD

      integer:: iaux,nth,nph


      ! Announce we are starting
      if (gpid.eq.0) then
        umsg = '   Iteration            MRC(J^0_0) Freq_index  '// &
               'Wavelength Height_index Height(km)'
        call verbose
      end if

      ! Initialize angle depended flag
      AD = .not.AVI
      ADD = AD.or.dyn

      ! B-B
      if (Input%init_J_bb) then

        ! Store if it was PRD and say it is not PRD now
        PRDl = PRD
        PRD = .False.

      end if

      ! MPI
      if (MPID%mpi) then

        ! CPU limits
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! Master
        if (pid.eq.0) then

          ! Alternative
          if (MPID%alternJ) then

            ! To receive Intensity chunks
            iaux = MPID%nxfreq*Rnz
            nsend = Geom%nTh*Geom%nPh*MPID%nnd

          ! Normal
          else

            ! Dimensions
            npz = Geom%nph*Rnz
            ntpz = Geom%nth*npz
            nsend = MPID%nnd

            ! To receive Intensity chunks
            iaux = MPID%nxfreq*Geom%nTh*Geom%nPh*Rnz

          end if

          ! Allocate receiver
          allocate(Stokes_r(iaux))

          ! Allocate new J00C
          allocate(J00C_n(nfreq,Rz0:Rz1))

        ! Slave
        else

          ! Alternative
          if (MPID%alternJ) then

            ! Trivial
            nTh = 1
            nPh = 1

          ! Normal
          else

            ! From Geom
            nTh = Geom%nTh
            nPh = Geom%nPh

          end if ! Alternative/normal MPI

          ! To send Intensity chunks
          allocate(Stokes_s(if0:if1,Rz0:Rz1,nPh,nTh))

          ! Allocate M and O pointers for RT coeff
          allocate(data1M(MPID%nf(pid),3))
          allocate(data1O(MPID%nf(pid),3))

        end if ! Master/slave

      ! Serial
      else

        ! CPU limits
        if0 = 1
        if1 = nfreq

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(nfreq,3))
        allocate(data1O(nfreq,3))

        ! Allocate new Stokes and J00C
        allocate(Stokes_s(nfreq,Geom%nPh,Geom%nTh,giz0:giz1))
        allocate(J00C_n(nfreq,Rz0:Rz1))

      end if ! MPI/serial

      end subroutine solveJ_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Parallel-master task for the self-consistent problem for
      !! continuum intensity\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!          MPID(MPI_class): Structure with MPI data\n
      !!             npz(integer): Size in z and phi\n
      !!            ntpz(integer): Size in z, theta, and phi\n
      !!           nsend(integer): Expected messages per iteration\n
      !!      Stokes_r(double(:)): Receiver buffer\n
      !!  Stokes(double(:,:,:,:)): Intensity\n
      !!        J00C(double(:,:)): Mean intensity with frequency
      !!                           dependence
      subroutine solveJ_manager(Geom,MPID,npz,ntpz,nsend, &
                                Stokes_r,Stokes,J00C)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      integer, intent(in):: npz,ntpz,nsend
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: Stokes_r
      double precision, &
             dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             target, intent(out):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(out):: J00C

      ! Local

      integer:: iz,ith,iph,itz,itpz,id,info_b,ierr,if0l,if1l,nfl
      integer, dimension(3):: info_c

      double precision:: WA

      ! Pointer

      double precision, dimension(:,:), pointer:: p_MStk


      ! Nullify
      nullify(p_MStk)

      !
      ! Alternative MPI
      !
      if (MPID%alternJ) then

        ! For expected messages
        do id=1,nsend

          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_c(1),3, &
                          MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Get indexes
          info_b = info_c(1)
          ith = info_c(2)
          iph = info_c(3)

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive intensity
          do while (.True.)
            call MPI_recv(Stokes_r(1), &
                          MPID%sizei4b(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          nfl = MPID%nf(info_b)

          ! Pointers
          p_MStk(if0l:if1l,Rz0:Rz1) => &
                               Stokes_r(1:MPID%sizei4b(info_b))

          ! Get angular weight
          WA = Geom%W_mu(ith)*Geom%W_mux(iph)

          ! Each height
          do iz=Rz0,Rz1

            ! Determine where to store intensity
            if (KSTK.or.iz.eq.Rz0) &
              Stokes(if0l:if1l,iph,ith,iz) = &
                                            p_MStk(if0l:if1l,iz)

            !
            ! Calculate frequency integral
            !
            call FIntJ(MPID,WA,info_b,p_MStk(:,iz),J00C(:,iz))

          end do ! heights
        end do ! Expected messages

      ! Normal MPI
      else

        ! For expected messages
        do id=1,nsend

          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_b,1, &
                          MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive intensity
          do while (.True.)
            call MPI_recv(Stokes_r(1), &
                          MPID%sizei4b(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          nfl = MPID%nf(info_b)

          ! Pointer
          p_MStk(if0l:if1l,1:ntpz) => &
                                   Stokes_r(1:MPID%sizei4b(info_b))

          ! For each height
          do iz=Rz0,Rz1

            ! For each polar direction
            do ith=1,Geom%nth

              ! Partial indexing
              itz = iz + (ith-1)*npz - Rz0 + 1

              ! For each azimuth
              do iph=1,Geom%nph

                ! Get frequency weight
                WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                ! Get running index
                itpz = itz + Rnz*(iph-1)

                ! Determine where to store intensity
                if (KSTK.or.iz.eq.Rz0) &
                  Stokes(if0l:if1l,iph,ith,iz) = &
                                              p_MStk(if0l:if1l,itpz)

                !
                ! Calculate frequency integral
                !
                call FIntJ(MPID,WA,info_b,p_MStk(:,itpz),J00C(:,iz))


              enddo ! azimuthal directions
            enddo ! polar directions
          enddo ! Heights
        end do ! Expected messages

      end if ! Type of MPI

      ! Nullify pointer
      nullify(p_MStk)

      end subroutine solveJ_manager

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the RTE for continuum intensity for the angular
      !! quadrature\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!                if0(integer): Initial frequency index\n
      !!                if1(integer): Final frequency index\n
      !!         data1M(double(:,:)): RT coeff point M\n
      !!         data1O(double(:,:)): RT coeff point O\n
      !!   Stokes_s(double(:,:,:,:)): Sender buffer\n
      !!     Stokes(double(:,:,:,:)): Intensity\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence\n
      !!         J00C_n(double(:,:)): Mean intensity with frequency
      !!                              dependence (new)
      subroutine solveJ_RT(Atom,Atmo,LTElines,Cont,Frec,Red,Geom, &
                           MPID,Input,if0,if1,data1M,data1O, &
                           Stokes_s,Stokes,J00C,J00C_n)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Continuum_class), intent(in):: Cont
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      integer, intent(in):: if0,if1
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, &
             dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             target, intent(inout):: Stokes
      double precision, dimension(nfreq,Rz0:Rz1), intent(inout):: J00C
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C_n
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: Stokes_s

      ! Local

      integer:: ith,iph,jdir,diz,iz,iz0,iz1,m,o,p,op,mfreq,nodir,mjdir
      integer,dimension(3):: info_c

      double precision:: mu_inv,dsm,dsp,ct,st,cc,sc,vfac
      double precision, dimension(1):: daux1

      ! Dummy

      double precision, dimension(:), allocatable:: ad1
      double precision, dimension(:,:), allocatable:: ad2
      double precision, dimension(:,:,:), allocatable:: ad3

      ! Pointers

      double precision, dimension(:), pointer:: p_K0M,p_SM,p_StkM
      double precision, dimension(:), pointer:: p_K0O,p_SO,p_StkO
      double precision, dimension(:), pointer:: p_K0P,p_SP
      double precision, dimension(:,:), pointer:: data1P


      ! Nullify pointers
      nullify(p_K0M,p_SM,p_StkM,p_K0O,p_SP,p_StkO,p_K0P,p_SP)
      nullify(data1P)

      ! Number of frequencies
      mfreq = if1 - if0 + 1

      ! Initialize op
      op = 1

      ! Initialize Doppler shift
      vfac = 1d0

      ! Maximum direction index for norm
      mjdir = ubound(Red%idzao,3)

      ! If including bound-bound
      if (Input%init_J_bb) then
        allocate(ad1(1),ad2(1,1),ad3(1,1,1))
        ad1 = 0d0
        ad2 = 0d0
        ad3 = 0d0
      end if

      !
      ! Ratiation Transfer
      !

      !  For each polar direction
      do ith=1,Geom%nTh

        ! Calculate inverse of cosine of polar direction
        mu_inv = 1d0/Geom%V_mu(ith)

        ! Determine the direction of propagation for indexes
        diz = -int(sign(1d0, Geom%V_mu(ith)))


        ! Determine the first and last height indexes to run
        ! over

        iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
        iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

        ! Trigonometry for Doppler shift
        if (dyn) then
          ct = Geom%v_mu(ith)
          st = sqrt(1d0 - ct*ct)
        end if

        ! For each azimuthal direction
        do iph=1,Geom%nPh

          ! Get direction index
          jdir = Geom%i_geom(iph,ith)
          nodir = min(mjdir,Geom%njdir,jdir)


          !
          ! First height
          !

          ! Get Doppler shift
          if (dyn) then

            ! Trigonometry
            cc = Geom%v_mux(iph)
            sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

            ! Amplitude
            vfac = sqrt(Atmo%vx(iz0)*Atmo%vx(iz0) + &
                        Atmo%vy(iz0)*Atmo%vy(iz0) + &
                        Atmo%vz(iz0)*Atmo%vz(iz0))

            ! Doppler shift
            if (vfac.gt.TINYVEL) then

              ! Shift
              vfac = 1d0 - Atmo%vx(iz0)*st*cc - &
                           Atmo%vy(iz0)*st*sc - &
                           Atmo%vz(iz0)*ct

            ! Static
            else

              ! No shift
              vfac = 1d0

            end if
          end if ! Dynamic

          ! If going down, get top boundary
          if(diz.eq.1)then

            ! Call top boundary
            call topI(if0,if1,data1M(:,3))

          ! If going up, get bottom boundary
          else

            ! Call bottom boundary
            call bottomI(Frec%omega,Atmo%T(iz0), &
                         vfac,if0,if1,data1M(:,3))

          endif ! propagation direction

          ! Identify current height
          o = iz0

          ! Calculate radiative coefficients
          if (Input%init_J_bb) then
            call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                          o,jdir,nodir,if0,if1, &
                          J00C(:,o),Cont%ndir, &
                          Cont%c(:,:,:,o), &
                          ad1,ad1,data1M(:,1:2),ad2,.False.)
          else
            call RTCoeffJ(jdir,if0,if1,J00C(:,o), &
                          Cont%ndir,Cont%c(:,:,:,o), &
                          data1M(:,1:2))
          end if

          ! If MPI
          if (pid.gt.0) then

            ! Alternative
            if (MPID%alternJ) then

              ! Store in buffer
              Stokes_s(:,o,1,1) = data1M(:,3)

            ! Normal
            else

              ! Stokes in buffer
              Stokes_s(:,o,iph,ith) = data1M(:,3)

            end if

          ! If serial
          else

            ! Save Stokes
            if (KSTK) Stokes_s(:,iph,ith,o) = data1M(:,3)

            ! Calculate integrals
            call JcalcJ(Geom,iph,ith,data1M(:,3),J00C_n(:,o))

          end if ! MPI/serial

          ! Identify next height
          p = iz0 + diz

          ! Get Doppler shift
          if (dyn) then

            ! Amplitude
            vfac = sqrt(Atmo%vx(p)*Atmo%vx(p) + &
                        Atmo%vy(p)*Atmo%vy(p) + &
                        Atmo%vz(p)*Atmo%vz(p))

            ! Doppler shift
            if (vfac.gt.TINYVEL) then

              ! Shift
              vfac = 1d0 - Atmo%vx(p)*st*cc - &
                           Atmo%vy(p)*st*sc - &
                           Atmo%vz(p)*ct

            ! Static
            else

              ! No shift
              vfac = 1d0

            end if
          end if ! If dynamic

          ! Calculate radiative coefficients
          if (Input%init_J_bb) then
            call RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac, &
                          p,jdir,nodir,if0,if1, &
                          J00C(:,p),Cont%ndir, &
                          Cont%c(:,:,:,p), &
                          ad1,ad1,data1O(:,1:2),ad2,.False.)
          else
            call RTCoeffJ(jdir,if0,if1,J00C(:,p), &
                          Cont%ndir,Cont%c(:,:,:,p), &
                          data1O(:,1:2))
          end if

          !
          ! Intermediate heights
          !

          ! For each height this CPU has assigned
          do iz=iz0,iz1,diz

            ! We treat the boundaries outside
            if(iz.eq.iz0.or.iz.eq.iz1)cycle

            ! Allocate P pointers
            allocate(data1P(MPID%nf(pid),3))

            ! Identify heights
            m = iz - diz
            o = iz
            p = iz + diz

            ! Calculate distance to previous point
            dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

            ! Calculate distance to the next point
            dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

            ! If tau scale
            if (ztau) then
              dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(m))
              dsp = abs(dsp)*2d0/(Atmo%chi500(o) + &
                                  Atmo%chi500(p))
            end if

            ! Get Doppler shift
            if (dyn) then

              ! Amplitude
              vfac = sqrt(Atmo%vx(p)*Atmo%vx(p) + &
                          Atmo%vy(p)*Atmo%vy(p) + &
                          Atmo%vz(p)*Atmo%vz(p))

              ! Doppler shift
              if (vfac.gt.TINYVEL) then

                ! Shift
                vfac = 1d0 - Atmo%vx(p)*st*cc - &
                             Atmo%vy(p)*st*sc - &
                             Atmo%vz(p)*ct

              ! Static
              else

                ! No shift
                vfac = 1d0

              end if
            end if ! If dynamic

            ! Calculate radiative coefficients
            if (Input%init_J_bb) then
              call RTCoeffI(Frec,Red,Atom,LTElines,Atmo, &
                            vfac,p,jdir,nodir,if0,if1, &
                            J00C(:,p),Cont%ndir, &
                            Cont%c(:,:,:,p), &
                            ad1,ad1,data1P(:,1:2),ad2,.False.)
            else
              call RTCoeffJ(jdir,if0,if1,J00C(:,p), &
                            Cont%ndir,Cont%c(:,:,:,p), &
                            data1P(:,1:2))
            end if

            ! Point to the data
            p_K0M  => data1M(:,1)
            p_SM   => data1M(:,2)
            p_StkM => data1M(:,3)
            p_K0O  => data1O(:,1)
            p_SO   => data1O(:,2)
            p_StkO => data1O(:,3)
            p_K0P  => data1P(:,1)
            p_SP   => data1P(:,2)

            ! Apply short characteristics BESSER
            call RTStepI(o,ith,iph,mfreq, &
                         dsm,dsp,p_K0M,p_SM,p_K0O, &
                         p_SO,p_K0P,p_SP,p_StkM, &
                         p_StkO,daux1,.False.,.True.)

            ! If MPI
            if (pid.gt.0) then

              ! Alternative
              if (MPID%alternJ) then

                ! Store in buffer
                Stokes_s(:,o,1,1) = data1O(:,3)

              ! Normal
              else

                ! Stokes in buffer
                Stokes_s(:,o,iph,ith) = data1O(:,3)

              end if

            ! If serial
            else

              ! Save Stokes
              if (KSTK) Stokes_s(:,iph,ith,o) = data1O(:,3)

              ! Calculate integrals
              call JcalcJ(Geom,iph,ith,data1O(:,3),J00C_n(:,o))

            end if ! MPI/serial

            ! Shift data (O->M, P->O)
            deallocate(data1M)
            data1M => data1O
            data1O => data1P
            nullify(data1P)

          end do ! Intermediate heights


          !
          ! Last height
          !

          ! Identify heights
          m = iz1 - diz
          o = iz1

          ! Calculate distance to previous point
          dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

          ! If tau scale
          if (ztau) &
            dsm = abs(dsm)*2d0/(Atmo%chi500(o) + &
                                Atmo%chi500(m))

          ! Point to the data
          p_K0M  => data1M(:,1)
          p_SM   => data1M(:,2)
          p_StkM => data1M(:,3)
          p_K0O  => data1O(:,1)
          p_SO   => data1O(:,2)
          p_StkO => data1O(:,3)

          ! Apply short characteristics LINEAR
          call RTStepI(o,ith,iph,mfreq, &
                       dsm,dsp,p_K0M,p_SM,p_K0O, &
                       p_SO,p_K0P,p_SP,p_StkM, &
                       p_StkO,daux1,.False.,.False.)

          ! If MPI
          if (pid.gt.0) then

            ! Alternative
            if (MPID%alternJ) then

              !
              ! Send to master if error
              !

              ! If had an error
              if (laborted) then

                ! Send error
                info_c = (/ -pid, ith, iph /)
                do while (.True.)
                  call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                                MPI_COMM_RT,ierr)
                  if (ierr.eq.0) exit
                end do

                cycle

              end if

              ! Store in buffer
              Stokes_s(:,o,1,1) = data1O(:,3)

              !
              ! Send to master
              !

              ! Send indexes
              info_c = (/ pid, ith, iph /)
              do while (.True.)
                call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                              MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do

              ! Send Stokes
              do while (.True.)
                call MPI_SEND(Stokes_s(if0,Rz0,1,1), &
                              MPID%sizei4b(pid), &
                              MPI_DOUBLE_PRECISION, 0, pid, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

            ! Normal
            else

              ! Stokes in buffer
              Stokes_s(:,o,iph,ith) = data1O(:,3)

            end if

          ! If serial
          else

            ! Save Stokes
            if (KSTK) Stokes_s(:,iph,ith,o) = data1O(:,3)

            ! Calculate integrals
            call JcalcJ(Geom,iph,ith,data1O(:,3),J00C_n(:,o))

          end if ! MPI/serial

        enddo ! azimuthal angles
      enddo ! polar angles

      ! If MPI and normal
      if (pid.gt.0.and..not.MPID%alternJ) then

        !
        ! Send to master
        !

        ! If had an error
        if (laborted) then

          ! Send error
          do while (.True.)
            call MPI_SEND(-pid,1,MPI_INTEGER,0,0, &
                          MPI_COMM_RT,ierr)
            if (ierr.eq.0) exit
          end do

        ! No problems
        else

          ! Send indexes
          do while (.True.)
            call MPI_SEND(pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

          ! Send Stokes
          do while (.True.)
            call MPI_SEND(Stokes_s(if0,Rz0,1,1), &
                          MPID%sizei4b(pid), &
                          MPI_DOUBLE_PRECISION, 0, pid, &
                          MPI_COMM_RT, ierr)
            if (ierr.eq.0) exit
          end do

        end if ! Errors
      end if ! If normal MPI

      ! Serial and storing Stokes
      if (pid.eq.0.and.KSTK) then

        ! Overwrite
        Stokes = Stokes_s

      end if

      ! Free pointers
      nullify(p_K0M,p_K0O,p_K0P,p_SM,p_SO,p_SP,p_StkM,p_StkO)

      ! Return
      return

      end subroutine solveJ_RT

!#####################################################################
!#####################################################################
!#####################################################################

      !> Conversion of radiation field tensors from multi-level
      !! intensity to, in general multi-term, polarization\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!     Rho_old(Rhoc_class(:)): Structure to store the density
      !!                             matrix of the previous
      !!                             iteration\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!             Red(Red_class): Structure with redistribution
      !!                             input frequency data,
      !!                             redistribution function data,
      !!                             and profile or normalization
      !!                             data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!             Pcorr(logical): If the mean intensities are to be
      !!                             corrected\n
      !!       Bfield(Bfield_blass): Structure with magnetic field
      !!                             data\n
      !!              rnPh(integer): Allocation size for Stokes\n
      !!   Stokes0(double(:,:,:,:)): Intensity\n
      !!           J00(double(:,:)): Mean intensity integrated over
      !!                             the absorption profile\n
      !!          J00S(double(:,:)): Mean intensity integrated over
      !!                             the emission profile\n
      !!          J00C(double(:,:)): Mean intensity with frequency
      !!                             dependence\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates
      subroutine JKQgen(Atom,Rho_old,Atmo,Frec,Red,Geom,MPID,Input, &
                        Flgsg,Pcorr,Bfield,rnPh,Stokes0,J00,J00S, &
                        J00C,Stokes,JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Rhoc_class), dimension(:), intent(inout):: Rho_old
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Input_class), intent(in):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: Pcorr
      integer, intent(in):: rnPh
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: Stokes0
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(inout):: J00P
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Stokes
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQC

      ! Local

      type(MRC_class):: MRC

      logical:: ADD

      integer:: if0,if1,ifreq,iph,nth,nph,npz,ntpz,nsend,ierr
      integer:: iz,ia,itran,jtran,ftran,fftran,jftran
      integer, dimension(:), allocatable:: psize

      double precision, dimension(:,:), allocatable:: LambdaL
      double precision, dimension(:,:,:), allocatable:: LambdaP
      double precision, dimension(:,:,:,:,:), allocatable:: Norm
      double precision, dimension(:,:,:,:,:), allocatable:: BStk

      ! Buffers
      ! Receivers
      double precision, dimension(:), allocatable:: Prof_r
      ! Senders
      double precision, dimension(:,:,:,:,:), allocatable:: Prof_s

      ! Pointers
      double precision, dimension(:,:), pointer:: data2O


      ! Nullify pointers
      nullify(data2O)

      ! Initialize
      call JKQgen_init(Atom,Rho_old,Frec,Geom,MPID,ADD,if0,if1, &
                       nth,nph,npz,ntpz,nsend,psize,Norm,BStk, &
                       LambdaL,LambdaP,JKQ,JKQS,Prof_r,Prof_s,data2O)

      ! Control
      call control
      if (laborted) goto 2000

      ! If MPI and master
      if (MPID%mpi.and.pid.eq.0) then

        ! Call manager
        call JKQgen_manager(Atom,Frec,Geom,MPID,ADD,nth,nph,npz, &
                            ntpz,nsend,psize,Norm,BStk, &
                            Prof_r,Stokes0,J00C,JKQ,JKQS)

      ! Serial or slave
      else

        ! Get profile data
        call JKQgen_RT(Atom,Atmo,Frec,Red,Geom,MPID,Flgsg,Bfield, &
                       ADD,nth,nph,if0,if1,psize,data2O,Prof_s, &
                       Stokes0,J00C,JKQ,JKQS)

      end if

      ! Control
      call control
      if (laborted) goto 2000

      ! If MPI
      if (MPID%mpi) then

        ! Share JKQ
        call MPI_BCAST(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                       MPI_DOUBLE_COMPLEX, 0, &
                       MPI_COMM_RT, ierr)

        ! Share JKQS if stimulated emission
        if (stm) &
          call MPI_BCAST(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)
      end if

      !
      ! If Polarization initial correction
      if (Pcorr) then

        !
        ! Make the J00 and J00S flat
        !

        ! For each height
        do iz=Rz0,Rz1

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Apply shift
              jtran = itran + Atom(ia)%tshift

              ! For each FS transition within this transition
              do ftran=1,Atom(ia)%fst(itran)%nt

                ! Get the global transition index
                fftran = Atom(ia)%ifst_ij(ftran,itran)

                ! Apply shift
                jftran = fftran + Atom(ia)%tfshift

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00(jftran,iz) = dble(JKQ(0,0,jtran,iz))

                ! Each transition gets the same contribution of
                ! the term-term J00, as an initial condition
                J00S(jftran,iz) = dble(JKQS(0,0,jtran,iz))

              end do ! fs transition
            end do ! term transition
          end do ! atom
        end do ! height


        !
        ! Solve SEE
        !

        ! For each atom
        do ia=1,nA

          ! Limiting indexes
          itran = Atom(ia)%tfshift + 1
          jtran = itran + Atom(ia)%nftran - 1
          fftran = Atom(ia)%pshift + 1
          jftran = fftran + Atom(ia)%nphot - 1

          ! For each height
          do iz=Rz0,Rz1

            ! Solve the SEE
#ifdef DEBUGSEE
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,.False.,.True.,Input)
#else
            call SEEI(Atom(ia),Rho_old(ia),J00(itran:jtran,iz), &
                      J00(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                     !J00S(itran:jtran,iz),J00P(fftran:jftran,:,iz), &
                      LambdaL,LambdaP,iz,.False.,.False.,.True.)
#endif

          end do ! heights
        end do ! atoms

        !
        !  Calculate MRC
        !

        ! Only the master does
        if (pid.eq.0) then

          ! Calculate MRC
          call MRC_sb(Atom,Rho_old,Input%anisotropy_only,MRC)

          ! Convert cm into km
          MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5

          ! Write in stdout
          if (gpid.eq.0) then
            write(umsg,'(A,1x,es22.12)') ' - The conversion to '// &
                       'terms have changed the populations a '// &
                       'maximum of ',MRC%values(2,1)
            call verbose
          end if

        end if
      end if ! P correction

      !
      ! Stokes
      !

      ! Allocate the extra Stokes parameters (memory counter in hanle)
      if (KSTK) then
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz1))
      else
        allocate(Stokes(0:3,nfreq,rnPh,Geom%nTh,Rz0:Rz0+1))
      end if

      ! Initialize to 0
      Stokes = 0d0

      ! Get the intensity from the intensity array
      if (KSTK) then
        if (axiali.and..not.axial) then
          do iph=1,Geom%nPh
            Stokes(0,:,iph,:,:) = Stokes0(:,1,:,:)
          end do
        else
          Stokes(0,:,:,:,:) = Stokes0
        end if
      end if

      ! Deallocate the intensity array (memory discounted in hanle)
      deallocate(Stokes0)

      ! Deallocate J00 and J00S (memory discounted in hanle)
      deallocate(J00,J00S)


      !
      ! JKQ frequency dependent
      !

      ! Allocate JKQC (memory counted in hanle)
      allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
      JKQC = cZero

      ! Get J00C into the complex array
      do iz=Rz0,Rz1
        do ifreq=1,nfreq

          JKQC(0,0,ifreq,iz) = dcmplx(J00C(ifreq,iz),0d0)

        end do ! frequencies
      end do ! heights

      ! Deallocate J00C (memory discounted in hanle)
      deallocate(J00C)

      ! Control
      call control

2000  if (associated(data2O)) then
        deallocate(data2O)
        nullify(data2O)
      end if

      ! Put back the AD
      if (tbAD) AV = .False.

      end subroutine JKQgen

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize arrays and pointers for the conversion of
      !! radiation field tensors from multi-level intensity to, in
      !! general multi-term, polarization\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!    Rho_old(Rhoc_class(:)): Structure to store the density
      !!                            matrix of the previous iteration\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!              ADD(logical): If Stokes needed for PRD-AD\n
      !!              if0(integer): Initial frequency index\n
      !!              if1(integer): Final frequency index\n
      !!              nth(integer): Number of polar directions to do\n
      !!              nph(integer): Number of azimuth to do\n
      !!              npz(integer): Size in z and phi\n
      !!             ntpz(integer): Size in z, theta, and phi\n
      !!            nsend(integer): Number of MPI communications\n
      !!         psize(integer(:)): Package size for MPI\n
      !!    LambdaL(double(:,:,:)): Lambda operator for lines\n
      !!  LambdaP(double(:,:,:,:)): Lambda operator for
      !!                            photoionizations\n
      !!   Norm(double(:,:,:,:,:)): Profile integral\n
      !!   BStk(double(:,:,:,:,:)): Integrand for radiation field
      !!                            tensors\n
      !!     JKQ(dcomplx(:,:,:,:)): Radiation field tensors integrated
      !!                            over the absorption profile\n
      !!    JKQS(dcomplx(:,:,:,:)): Radiation field tensors integrated
      !!                            over the emission profile\n
      !!         Prof_r(double(:)): Receiver buffer for profiles\n
      !!   Prof_s(double(:,:,:,:)): Sender buffer for profiles\n
      !!       data2O(double(:,:)): Profiles point O
      subroutine JKQgen_init(Atom,Rho_old,Frec,Geom,MPID,ADD,if0, &
                             if1,nth,nph,npz,ntpz,nsend,psize,Norm, &
                             BStk,LambdaL,LambdaP,JKQ,JKQS,Prof_r, &
                             Prof_s,data2O)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Rhoc_class), dimension(:), intent(inout):: Rho_old
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      logical, intent(out):: ADD
      integer, intent(out):: if0,if1,nth,nph,npz,ntpz,nsend
      integer, dimension(:), allocatable, intent(out):: psize
      double precision, dimension(:), &
                        allocatable, intent(out):: Prof_r
      double precision, dimension(:,:), &
                        allocatable, intent(out):: LambdaL
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: LambdaP
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Prof_s
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Norm
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: BStk
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQS
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data2O

      ! Local

      logical:: AD

      integer:: int1,int2,ia,iaux,iproc


      ! Trick to have AV input for AD calculation
      if (tbAD) AV = .True.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! Select the angular limits
      if ((PRD.and.ADD).or.dyn) then
        nth = Geom%nTh
        nph = Geom%nPh
      else
        nth = 1
        nph = 1
      end if

      ! Initialize lambda operator to 0
      int1 = 1
      int2 = 1
      do ia=1,nA
        if (Atom(ia)%nftran.gt.int1) int1 = Atom(ia)%nftran
        if (Atom(ia)%nphot.gt.int2) int2 = Atom(ia)%nphot
      end do
      allocate(LambdaL(nxb,int1))
      LambdaL = 0d0
      allocate(LambdaP(nxb,int2,2))
      LambdaP = 0d0

      ! MPI
      if (MPID%mpi) then

        ! Allocate
        allocate(psize(0:nproc-1))

        ! CPU limits
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! Alternative
        if (MPID%alternJgen) then

          ! Calculate size of transmition package
          do iproc=0,nproc-1
            psize(iproc) = 2*Frec%ntfreq*Rnz
          end do

          ! Master
          if (pid.eq.0) then

            ! Norm
            allocate(Norm(2,nxtran,nph,nth,Rz0:Rz1))

            ! Cumulative JKQ
            allocate(BStk(2,nxtran,nph,nth,Rz0:Rz1))

            ! Zero out
            Norm = 0d0
            BStk = 0d0

            ! To receive profile information
            iaux = MPID%nxtfreq*2*Rnz
            allocate(Prof_r(iaux))

            ! Number of communications
            nsend = nTh*nPh*MPID%nnd

          ! Slave
          else

            ! Allocate O pointers
            allocate(data2O(Frec%ntfreq,2))

            ! To send profile information
            allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1,1,1))

          end if ! Master/slave

        ! Normal
        else

          ! Calculate size of transmition package
          do iproc=0,nproc-1
            psize(iproc) = 2*Frec%Mntfreq(iproc)*nph*nth*Rnz
          end do

          ! Master
          if (pid.eq.0) then

            ! Dimensions
            npz = nPh*Rnz
            ntpz = nTh*npz

            ! Norm
            allocate(Norm(2,nxtran,nph,nth,Rz0:Rz1))

            ! Cumulative JKQ
            allocate(BStk(2,nxtran,nph,nth,Rz0:Rz1))

            ! To receive profile information
            iaux = MPID%nxtfreq*2*nTh*nPh*Rnz
            allocate(Prof_r(iaux))

            ! Zero out
            Norm = 0d0
            BStk = 0d0

            ! Number of communications
            nsend = MPID%nnd

          ! Slave
          else

            ! Allocate O pointers
            allocate(data2O(Frec%ntfreq,2))

            ! To send profile information
            allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1,nph,nth))

          end if ! Master/slave

        end if

      ! Serial
      else

        ! CPU limits
        if0 = 1
        if1 = nfreq

        ! Allocate O pointers
        allocate(data2O(Frec%ntfreq,2))

      end if

      !
      ! Allocate Radiation field
      !
      ! Allocate JKQ and JKQS and initialize to 0 (memory counted in
      ! hanle)
      allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
      allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))

      ! Master zero out
      if (pid.eq.0) then
        JKQ = cZero
        JKQS = cZero
      end if

      ! For each atom
      do ia=1,nA
        Rho_old(ia)%crho = Atom(ia)%crho(:,Rz0:Rz1)
      end do

      end subroutine JKQgen_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Parallel-master task for the conversion of radiation field
      !! tensors from multi-level intensity to, in general multi-term,
      !! polarization\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!              ADD(logical): If Stokes needed for PRD-AD\n
      !!              nth(integer): Number of polar directions to do\n
      !!              nph(integer): Number of azimuth to do\n
      !!              npz(integer): Size in z and phi\n
      !!             ntpz(integer): Size in z, theta, and phi\n
      !!            nsend(integer): Number of MPI communications\n
      !!         psize(integer(:)): Package size for MPI\n
      !!   Norm(double(:,:,:,:,:)): Profile integral\n
      !!   BStk(double(:,:,:,:,:)): Integrand for radiation field
      !!                            tensors\n
      !!         Prof_r(double(:)): Receiver buffer for profiles\n
      !!  Stokes0(double(:,:,:,:)): Intensity\n
      !!         J00C(double(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over the absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over the emission profile
      subroutine JKQgen_manager(Atom,Frec,Geom,MPID,ADD,nth,nph,npz, &
                                ntpz,nsend,psize,Norm,BStk, &
                                Prof_r,Stokes0,J00C,JKQ,JKQS)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: ADD
      integer, intent(in):: nth,nph,npz,ntpz,nsend
      integer, dimension(:), allocatable, intent(in):: psize
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: Prof_r
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Norm
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: BStk
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: Stokes0
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQS

      ! Local

      integer:: ith,iph,iz,if0l,if1l,nfl,nftl,id,ierr,info_b
      integer:: ia,itran,jtran,itpz
      integer, dimension(3):: info_c

      double precision:: WA,daux
      double precision, dimension(:,:,:), pointer:: p_MProf


      ! Nullify
      nullify(p_MProf)

      ! Reset radiation field variables
      Norm = 0d0
      BStk = 0d0

      !
      ! Alternative MPI
      !
      if (MPID%alternJgen) then

        ! For each communication
        do id=1,nsend

          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_c(1), 3, MPI_INTEGER, &
                          MPI_ANY_SOURCE, 0, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Get indexes
          info_b = info_c(1)
          ith = info_c(2)
          iph = info_c(3)

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), psize(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreq(info_b)

          ! Pointers
          p_MProf(1:nftl,1:2,Rz0:Rz1) => Prof_r(1:psize(info_b))

          ! Each height
          do iz=Rz0,Rz1

            !
            ! Calculate frequency integral
            !
            if ((PRD.and.ADD).or.dyn) then
              if (axiali) then
                call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                            Stokes0(:,1,ith,iz), &
                            p_MProf(:,:,iz),Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz))
              else
                call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                            Stokes0(:,iph,ith,iz), &
                            p_MProf(:,:,iz),Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz))
              end if
            else
              call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                          J00C(:,iz),p_MProf(:,:,iz), &
                          Norm(:,:,iph,ith,iz), &
                          BStk(:,:,iph,ith,iz))

            end if

          end do ! heights

        end do ! Communications

      !
      ! Normal MPI
      !
      else

        ! For each communication
        do id=1,nsend

          !
          ! Receive data from a slave
          !

          ! Receive indexing data
          do while (.True.)
            call MPI_recv(info_b,1, MPI_INTEGER, MPI_ANY_SOURCE, &
                          0, MPI_COMM_RT, MPI_STATUS_IGNORE, &
                          ierr)
            if (ierr.eq.0) exit
          end do

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), psize(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreq(info_b)

          ! Pointers
          p_MProf(1:nftl,1:2,1:ntpz) => Prof_r(1:psize(info_b))

          ! Compute line quantities
          do itpz=1,ntpz

            ! Get indexes
            ith = (itpz-1)/npz
            iph = (itpz - npz*ith - 1)/Rnz
            iz = itpz - Rnz*iph - npz*ith + Rz0 - 1
            ith = ith + 1
            iph = iph + 1

            !
            ! Calculate frequency integral
            !
            if ((PRD.and.ADD).or.dyn) then
              if (axiali) then
                call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                            Stokes0(:,1,ith,iz), &
                            p_MProf(:,:,itpz),Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz))
              else
                call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                            Stokes0(:,iph,ith,iz), &
                            p_MProf(:,:,itpz),Norm(:,:,iph,ith,iz), &
                            BStk(:,:,iph,ith,iz))
              end if
            else
              call FJgInt(Atom,MPID,Frec%W_freq,info_b, &
                          J00C(:,iz),p_MProf(:,:,itpz), &
                          Norm(:,:,iph,ith,iz), &
                          BStk(:,:,iph,ith,iz))

            end if

          end do ! heights/directions
        end do ! Communications

      end if ! Type of MPI

      ! Nullify pointers
      nullify(p_MProf)

      !
      ! Apply weights to JKQ, JKQS, and normalize if
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each polar direction
        do ith=1,nTh

          ! For each azimuthal direction
          do iph=1,nph

            ! Get the angular integral weight
            if ((PRD.and.ADD).or.dyn) then
              WA = Geom%W_mu(ith)*Geom%W_mux(iph)
            else
              WA = 1d0
            end if

            ! For each atom
            do ia=1,nA

              ! For each FS transition
              do itran=1,Atom(ia)%ntran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tshift

                ! Get the weight
                if (Norm(1,jtran,iph,ith,iz).gt.0d0) then

                  ! Inverse norm and angular weight
                  daux = WA/Norm(1,jtran,iph,ith,iz)

                  ! Integrate angle
                  JKQ(0,0,jtran,iz) = JKQ(0,0,jtran,iz) + &
                        dcmplx(BStk(1,jtran,iph,ith,iz)*daux, 0d0)

                end if

                ! If there is stimulated emission
                if (stm) then

                  ! Get the weight
                  if (Norm(2,jtran,iph,ith,iz).gt.0d0) then

                    ! Inverse norm and angular weight
                    daux = WA/Norm(2,jtran,iph,ith,iz)

                    ! Integrate angle
                    JKQS(0,0,jtran,iz) = JKQS(0,0,jtran,iz) + &
                        dcmplx(BStk(2,jtran,iph,ith,iz)*daux, 0d0)

                  end if ! Non-zero weight
                end if ! Stimulated emission

              end do ! transitions
            end do ! atoms
          end do ! azimuthal directions
        end do ! polar directions
      end do ! heights

      end subroutine JKQgen_manager

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate JKQ tensors from the solution of the intensity RT
      !! problem\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            input frequency data,
      !!                            redistribution function data,
      !!                            and profile or normalization
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!      Bfield(Bfield_blass): Structure with magnetic field
      !!                            data\n
      !!              ADD(logical): If Stokes needed for PRD-AD\n
      !!              nth(integer): Number of polar directions to do\n
      !!              nph(integer): Number of azimuth to do\n
      !!              if0(integer): Initial frequency index\n
      !!              if1(integer): Final frequency index\n
      !!         psize(integer(:)): Package size for MPI\n
      !!       data2O(double(:,:)): Profiles point O\n
      !!   Prof_s(double(:,:,:,:)): Sender buffer for profiles\n
      !!  Stokes0(double(:,:,:,:)): Intensity\n
      !!         J00C(double(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!    JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over the absorption profile\n
      !!   JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over the emission profile
      subroutine JKQgen_RT(Atom,Atmo,Frec,Red,Geom,MPID,Flgsg, &
                           Bfield,ADD,nth,nph,if0,if1,psize,data2O, &
                           Prof_s,Stokes0,J00C,JKQ,JKQS)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Fctsg_class), intent(inout):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: ADD
      integer, intent(in):: nth,nph,if0,if1
      integer, dimension(:), allocatable, intent(in):: psize
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(in):: Stokes0
      double precision, dimension(:,:), allocatable, intent(in):: J00C
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Prof_s
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQS
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data2O

      ! Local

      integer:: ith,iph,jdir,iz,ierr
      integer, dimension(3):: info_c

      double precision:: WA,ct,st,cc,sc,vfac

      complex(kind=8), dimension(:,:,:), pointer:: TSo,TKQo


      ! Nullify pointers
      nullify(TSo,TKQo)

      ! Initialize Doppler shift
      vfac = 1d0

      ! For each polar direction
      do ith=1,nTh

        ! Trigonometry for Doppler shift
        if (dyn) then
          ct = Geom%v_mu(ith)
          st = sqrt(1d0 - ct*ct)
        end if

        ! For each azimuthal direction
        do iph=1,nPh

          ! If serial (master)
          if (pid.eq.0) then

            ! Select the angular limits
            if ((PRD.and.ADD).or.dyn) then
              WA = Geom%W_mu(ith)*Geom%W_mux(iph)
            else
              WA = 1d0
            end if
          end if

          ! Trigonometry Doppler shift
          if (dyn) then
            cc = Geom%v_mux(iph)
            sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)
          end if

          ! Get direction index
          jdir = Geom%i_geom(iph,ith)

          ! Point TKQ_S
          TSo => Geom%TS(:,:,:,jdir)

          ! For each height this CPU has assigned
          do iz=Rz0,Rz1

            ! Select correct TKQout
            if (Bfield%Bstrength(iz).gt.TINYB) then
              TKQo => Geom%TB(:,:,:,jdir,iz)
            else
              TKQo => TSo
            end if

            ! Get Doppler shift
            if (dyn) then

              ! Amplitude
              vfac = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                          Atmo%vy(iz)*Atmo%vy(iz) + &
                          Atmo%vz(iz)*Atmo%vz(iz))

              ! Doppler shift
              if (vfac.gt.TINYVEL) then

                ! Shift
                vfac = 1d0 - Atmo%vx(iz)*st*cc - &
                             Atmo%vy(iz)*st*sc - &
                             Atmo%vz(iz)*ct

              ! Static
              else

                ! No shift
                vfac = 1d0

              end if
            end if ! If dynamic

            ! RT coefficients
            call Termprof(Frec,Red,Atom,Atmo,vfac,Flgsg,Geom, &
                          Bfield,iz,jdir,if0,if1,TSo,TKQo, &
                          data2O)

            ! Serial (master)
            if (pid.eq.0) then

              ! If PRD and angle-dependent or dynamic
              if ((PRD.and.ADD).or.dyn) then

                ! If intensity axially symmetric
                if (axiali) then

                  ! Get contribution
                  call Jgen(Atom,Frec%W_freq,WA, &
                            Stokes0(:,1,ith,iz),data2O, &
                            JKQ(:,:,:,iz),JKQS(:,:,:,iz))

                ! Intensity not axially symmetric
                else

                  ! Get contribution
                  call Jgen(Atom,Frec%W_freq,WA, &
                            Stokes0(:,iph,ith,iz),data2O, &
                            JKQ(:,:,:,iz),JKQS(:,:,:,iz))

                end if ! Intensity axial symmetry

              ! No PRD or static
              else

                ! Get contribution
                call Jgen(Atom,Frec%W_freq,WA,J00C(:,iz),data2O, &
                          JKQ(:,:,:,iz),JKQS(:,:,:,iz))

              end if ! PRD and dynamic

            ! Alternative MPI (slave)
            else if (MPID%alternJgen) then

              ! Store in buffer
              Prof_s(:,:,iz,1,1) = data2O(:,:)

            ! Normal MPI (slave)
            else

              ! Store in buffer
              Prof_s(:,:,iz,iph,ith) = data2O(:,:)

            end if ! Serial/type MPI

          end do ! Heights

          ! Alternative MPI (slave)
          if (pid.gt.0.and.MPID%alternJgen) then

            !
            ! Send to master if error
            !

            ! If had an error
            if (laborted) then

              ! Send error
              info_c = (/ -pid, ith, iph /)
              do while (.True.)
                call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                              MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do

              cycle

            end if

            !
            ! Send to master
            !

            ! Send indexes
            info_c = (/ pid, ith, iph /)
            do while (.True.)
              call MPI_SEND(info_c(1),3,MPI_INTEGER,0,0, &
                            MPI_COMM_RT,ierr)
              if (ierr.eq.0) exit
            end do

            ! Send profiles
            do while (.True.)
              call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                            psize(pid), MPI_DOUBLE_PRECISION, &
                            0, pid, MPI_COMM_RT, ierr)
              if (ierr.eq.0) exit
            end do

          end if ! Alternative MPI

        end do ! Azimuths
      end do ! Polar

      ! If normal MPI (slave)
      if (pid.gt.0.and..not.MPID%alternJgen) then

        !
        ! Send to master
        !

        ! If had an error
        if (laborted) then

          ! Send indexes
          do while (.True.)
            call MPI_SEND(-pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

        ! No problems
        else

          ! Send indexes
          do while (.True.)
            call MPI_SEND(pid,1,MPI_INTEGER,0,0,MPI_COMM_RT, &
                          ierr)
            if (ierr.eq.0) exit
          end do

          ! Send profiles
          do while (.True.)
            call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                          psize(pid), MPI_DOUBLE_PRECISION, &
                          0, pid, MPI_COMM_RT, ierr)
            if (ierr.eq.0) exit
          end do
        end if ! Errors
      end if ! Normal MPI

      ! Free
      nullify(TSo,TKQo)

      end subroutine JKQgen_RT

!#####################################################################
!#####################################################################
!#####################################################################

      end module solveri_mod
