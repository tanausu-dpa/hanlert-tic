      !> Solves the NLTE problem of the second kind
      module solver_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     20/02/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     20/02/2025:    V4.0.1 - Added argument to the MRCJKQ_sb and
!                             MRC_sb calls (TdPA)
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
!  solve
!    Self-consistent solution for the RT problem
!
!  solve_predict
!    Calculate the approximate RAM necessary for the routines solving
!  the self-consistent problem
!
!  solve_init
!    Initialize arrays and pointers for the self-consistent solution
!
!  solve_manager
!    Parallel-master task for the self-consistent problem
!
!  solve_RT
!    Solve the RTE for the angular quadrature
!
!  solve_SEE
!    Advance the density matrices by solving the SEE and applying NG
!  acceleration if requested
!
!  emergent
!    Formal solution for the RTE
!
!  emergent_predict
!    Calculate the approximate RAM space necessary for the routines
!  solving the formal solution for the RTE
!
!  emergent_init:
!    Initialize arrays and pointers for the formal solution for the
!  RTE
!
!  emergent_manager
!    Parallel-master task for the formal solution for the RTE
!
!  emergent_RT
!    Solve the RTE for a given LOS
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use boundary_mod
      use commons_mod
      use comovingprd_mod
      use free_mod
      use gauss_mod , only : setTKQLOS
      use iosolution_mod
      use jcalc_mod
      use mrc_mod
      use ng_mod
      use normalizer_mod
      use parameters_mod , only : B2L , cZero , kb, cSaha , fktoJ
      use rtcoeff_mod
      use rtstep_mod
      use rtstepi_mod
      use see_mod
      use setmpi_mod
      use types_mod

      ! Maximum buffer for NG_int
      double precision, parameter:: maxbuffer_NG = 500d0

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Self-consistent solution for the RT problem\n
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
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!   JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the
      !!                              radiation tensors\n
      !!   Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!      JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                              integrated over the absorption
      !!                              profile\n
      !!     JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                              integrated over the emission
      !!                              profile\n
      !!     JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                              frequency dependence
      subroutine solve(Atom,LTElines,Rho_old,Atmo,Cont,Frec,Red, &
                       Bfield,Geom,MPID,Input,Flgsg,JKQ_asym, &
                       Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Rhoc_class), dimension(:), intent(inout):: Rho_old
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Fctsg_class), intent(inout):: Flgsg
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      double precision, &
             dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             target, intent(inout):: Stokes
      complex(kind=8), dimension(:,:,:), intent(in):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(inout):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(inout):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(inout):: JKQC

      ! Local

      type(MRC_class):: MRC

      character(LEN=20):: iterS

      logical:: doNG,goout,gooutprd,ADD,RPRAM,lp_exu,NGP

      integer:: NG_dim,NG_entry,if0,if1,npz,ntpz,nsend,iter,iterr,ia

      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(3):: bMRC

      double precision:: loffset
      double precision, dimension(nxphot,2,Rz0:Rz1):: J00P
      double precision, dimension(:,:), allocatable:: NG_scratch

      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC_n
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC_old

      ! Receivers

      double precision, dimension(:), allocatable, target:: Stokes_r
      double precision, dimension(:), allocatable, target:: Prof_r

      ! Senders

      double precision, dimension(:,:,:,:,:), allocatable:: Stokes_s
      double precision, dimension(:,:,:,:,:), allocatable:: Prof_s

      ! Pointers

      double precision, dimension(:), pointer:: p_exu
      double precision, dimension(:,:), pointer:: data2O
      double precision, dimension(:,:,:), pointer:: data1M,data1O


#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,-2)
#endif
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,-2)
#endif

      ! Nullify pointers
      nullify(data1M,data1O,data2O,p_exu)

      ! Initialize
      call solve_init(Atom,Frec,Geom,MPID,Input,NG_dim,NG_entry, &
                      NG_scratch,NGP,doNG,goout,ADD,RPRAM,lp_exu, &
                      if0,if1,npz,ntpz,nsend,Stokes_r,Prof_r, &
                      Stokes_s,Prof_s,data1M,data1O,data2O,p_exu, &
                      JKQC_old,JKQC_n)

      ! Initialize PRD
      if (PRD) call initialize_emiss(Atom,Geom,Red)

      ! Control
      call control
      if (laborted) goto 2000

      ! If MPI and measuring performance
      if (MPID%mpi.and.Input%mpi_perf.and.pid.eq.0) &
        call report_mpi_time(Input%folder,Input%ID, &
                             0,0,.False.)

      !
      ! Start iterations
      !

      ! For each iteration between the limits specified
      do iter=Input%iter_min,Input%iter_max

        ! Flags for physics in Stokes
        if (iter.le.Input%allownphys_stk) then
          if (.not.nphysS) nphysS = .True.
        else
          if (nphysS) nphysS = .False.
        end if

        ! Flags for physics in density matrix
        if (iter.le.Input%allownphys_rho) then
          if (.not.nphysR) nphysR = .True.
        else
          if (nphysR) nphysR = .False.
        end if

        ! For each atom
        do ia=1,nA

          ! Copy density matrix
          Rho_old(ia)%crho = Atom(ia)%crho

        end do

        ! Internal PRD iterations
        do iterr=1,Input%iter_prd

          ! If PRD iterations
          if (pid.eq.0.and.PRD.and.Input%iter_prd.gt.1) &
            JKQC_old = JKQC(0:2,:,:,:)

          ! Compute second order emissivity
          if (PRD) &
            call comoving_emiss2ord(Atom,Atmo,Geom,Frec,Red, &
                                    Flgsg,Bfield,Stokes,JKQ_asym, &
                                    JKQ,JKQC,0,0,Input%PRD_int_mode, &
                                    .False.)

          ! If MPI and master
          if (MPID%mpi.and.pid.eq.0) then

            ! Call manager
            call solve_manager(Atom,Atmo,Frec,Geom,Bfield,MPID, &
                               Input,Flgsg,nsend,iter,npz,ntpz, &
                               Stokes_r,Prof_r,JKQ_asym,Stokes, &
                               JKQ,JKQS,JKQC,J00P)

          ! Serial or slave
          else
  
            ! Solve RTE
            call solve_RT(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                          Geom,MPID,Input,Flgsg,lp_exu,if0,if1, &
                          Stokes_s,Prof_s,data1M,data1O,data2O, &
                          p_exu,JKQ_asym,Stokes,JKQ,JKQS,JKQC, &
                          JKQC_n,J00P)

          end if ! Manage or compute

          ! Control
          call control
          if (laborted) goto 2000

          !
          ! Master
          !
          if (pid.eq.0) then

            ! Calculate MRC for J if PRD
            if (PRD.and.Input%iter_prd.gt.1) then

              ! Call the routine
              call MRCJKQ_sb(JKQC,JKQC_old,Input%anisotropy_only,MRC)

              ! Convert cm into km
              MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5
              MRC%values(1,2) = Atmo%z(MRC%indexes(2,2))*1d-5

              ! Global Máster talks
              if (gpid.eq.0) then

                ! If first PRD iteration and second actual
                ! iteration
                if (iterr.eq.1.and.iter.eq.1) then
                  umsg = '         PRD            MRC(J^0_0)'// &
                         ' Freq_index  Wavelength '// &
                         'Height_index Height(km)'// &
                         '  MRC(J^K_Q)'// &
                         ' Freq_index  Wavelength '// &
                         'Height_index Height(km)  K  Q'
                  call verbose
                end if

                ! Write in stdout
                write(umsg,'(2x,"PRD it:",1x,i3,2x,es20.12,'// &
                           '2x,i9,2x,f10.4,4x,i9,2x,f9.3,'// &
                           '2x,es20.12,2x,i9,2x,f10.4,4x,i9,'// &
                           '2x,f9.3,1x,i2,1x,i2)') &
                           iterr,MRC%values(2,1),MRC%indexes(1,1), &
                           1d2/Frec%omega(MRC%indexes(1,1)), &
                           MRC%indexes(2,1),MRC%values(1,1), &
                           MRC%values(2,2),MRC%indexes(1,2), &
                           1d2/Frec%omega(MRC%indexes(1,2)), &
                           MRC%indexes(2,2),MRC%values(1,2), &
                           MRC%indexes(3,2),MRC%indexes(4,2)
                call verbose

              end if ! Globalmaster

              ! If pass exit criteria
              if ((MRC%values(2,1).le.Input%mrc_r.and. &
                   MRC%values(2,2).le.Input%mrc_p_r).or. &
                  iterr.eq.Input%iter_prd) then

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

            !
            ! Share the radiation information
            !

            ! Share if we are finished or not (in PRD)
            call MPI_BCAST(gooutprd, 1, MPI_LOGICAL, 0, &
                           MPI_COMM_RT, ierr)

            ! Share JKQ
            call MPI_BCAST(JKQ(-2,0,1,Rz0), MPID%size6(0), &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

            ! Share JKQS if stimulated emission
            if (stm) &
            call MPI_BCAST(JKQS(-2,0,1,Rz0), MPID%size6(0), &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

            ! Share JKQC
            call MPI_BCAST(JKQC(-2,0,1,Rz0), MPID%size7(0), &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

            ! Share Stokes if doing A-D PRD
            if (PRD.and.ADD) &
            call MPI_BCAST(Stokes(0,1,1,1,giz0), MPID%size8(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Share J00 for b-f transitions
            call MPI_BCAST(J00P(1,1,Rz0), MPID%sizei3(0), &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)
          end if ! MPI

          ! If we are finished, exit
          if (gooutprd) exit

        end do ! Internal PRD iterations

        ! Control
        if (laborted) goto 2000

        ! Advance density matrix
        call solve_SEE(Atom,Rho_old,Atmo,Bfield,Geom,Input,Flgsg, &
                       MRC,NG_dim,NG_entry,NG_scratch,NGP,doNG, &
                       goout,ADD,iter,Stokes,JKQ,JKQS,JKQC,J00P)

        ! We can swith now to AD if we had AV input
        if (tbAD) then
          AV = .False.
          tbAD = .False.
          PRAM = RPRAM
          ! Force at least two iterations
          goout = .False.
        end if

        !
        ! Save partial solution
        !

        ! Check if we want to store partial results
        if(Input%store.and.mod(iter,Input%store_step).eq.0)then

          ! Only the master writes
          write(iterS,'(i0.4)') iter
          call writesol(Input,iterS,Frec%omega,Geom,Flgsg, &
                        Bfield,Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)

        end if

        ! If MPI, share if finished
        if (MPID%mpi) then

          ! Share goout
          call MPI_BCAST(goout, 1, MPI_LOGICAL, 0, &
                         MPI_COMM_RT, ierr)

        end if

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
        if (ierr.ne.0) goto 1000

        !
        ! Jump columns
        !

        ! Get offset
        loffset = 11d0 + dble(icoords(3)-1)*12d0
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Store in buffer
        bMRC(1) = real(iter)
        bMRC(2) = real(MRC%values(2,1))
        bMRC(3) = real(MRC%values(2,2))

        ! Write
        call MPI_FILE_WRITE(funit,bMRC(1),3, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if ! 1.5D synthesis and master keeping MRC file

      ! Control
      call control

      !
      ! Clean pointers
      !
2000  if (associated(data1O)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
      end if
      if (.not.lp_exu) deallocate(p_exu)
      nullify(p_exu)

      ! Free PRD
      if (PRD) call free_e2ord(Red)


#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,-1)
#endif
#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,-1)
#endif

      ! Done
      return

1000  umsg = 'Error opening MRC file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      if (associated(data1O)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
      end if
      if (.not.lp_exu) deallocate(p_exu)
      nullify(p_exu)
      if (PRD) call free_e2ord(Red)
      return
1010  umsg = 'Error seeking MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      if (associated(data1O)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
      end if
      if (.not.lp_exu) deallocate(p_exu)
      nullify(p_exu)
      if (PRD) call free_e2ord(Red)
      return
1300  umsg = 'Error writing MRC file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      if (associated(data1O)) then
        deallocate(data1M,data1O,data2O)
        nullify(data1M,data1O,data2O)
      end if
      if (.not.lp_exu) deallocate(p_exu)
      nullify(p_exu)
      if (PRD) call free_e2ord(Red)
      return

      end subroutine solve

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the approximate RAM necessary for the routines
      !! solving the self-consistent problem\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!     Input(Input_class): Structure with configuration data
      subroutine solve_predict(Atom,Frec,Red,Geom,MPID,Input)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input


      ! Local

      logical:: AD,ADD

      integer:: ia,iaux,NG_dim,if0,if1,mfreq


      ! Initialize
      TRAMc = 0d0

      ! Predict PRD
      if (PRD) call predict_emiss(Atom,Geom,Red)

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! NG quantities

      ! If NG acceleration
      if (Input%NG) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + Rnz*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRD) then

          ! If we need Stokes
          if (ADD) then

            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            NG_dim = NG_dim + nfreq*RnZ

          end if !ADD
        end if ! PRD

        ! If it requires too much buffer
        if (dble(NG_dim)*8d-6.le.maxbuffer_NG) then

          ! Master
          if (pid.eq.0) then
            TRAMc = TRAMc + 8d-6*dble(NG_dim*(Input%NG_ord+2))
          ! Slave
          else
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
        mfreq = if1-if0+1

        ! Master
        if (pid.eq.0) then

          ! Alternative
          if (MPID%alternP) then

            ! To receive Intensity chunks
            iaux = 4*MPID%nxfreq*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

            ! To receive profile information
            iaux = MPID%nxtfreq*2*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

          ! Normal
          else

            ! To receive Intensity chunks
            iaux = 4*MPID%nxfreq*Geom%nTh*Geom%nPh*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

            ! To receive profile information
            iaux = MPID%nxtfreq*2*Geom%nTh*Geom%nPh*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)

          end if ! Type of MPI

          ! Norm
          TRAMc = TRAMc + 8d-6*dble(2*nxtran*Geom%nph*Geom%nth*Rnz)

          ! BStk
          TRAMc = TRAMc + 16d-6*dble(15*2*nxtran*Geom%nph* &
                                                 Geom%nth*Rnz)

        ! Slave
        else

          ! M and O pointers for RT coeff
          TRAMc = TRAMc + 8d-6*dble(2*4*mfreq*6)

          ! O pointers for profile
          TRAMc = TRAMc + 8d-6*dble(Frec%ntfreq*2)

          ! Alternative MPI
          if (MPID%alternP) then

            ! To send Intensity chunks
            TRAMc = TRAMc + 8d-6*dble(4*mfreq*Rnz)

            ! To send profile information
            TRAMc = TRAMc + 8d-6*dble(Frec%ntfreq*2*Rnz)

          ! Normal MPI
          else

            ! To send Intensity chunks
            TRAMc = TRAMc + 8d-6*dble(4*mfreq*Rnz* &
                                      Geom%nTh*Geom%nPh)

            ! To send profile information
            TRAMc = TRAMc + 8d-6*dble(Frec%ntfreq*2*Rnz* &
                                      Geom%nTh*Geom%nPh)

          end if ! Type of MPI
        end if ! Master/slave

      ! Serial
      else

        ! Common (Master and slave)
        ! O pointers
        TRAMc = TRAMc + 8d-6*dble(Frec%ntfreq*2)

        ! Allocate M and O pointers for RT coeff
        TRAMc = TRAMc + 8d-6*dble(2*4*nfreq*6)

        ! Initialize exponential is not in memory
        if (.not.(PIRAM.and.Frec%pif1.ge.Frec%pif0)) &
          TRAMc = TRAMc + 8d-6

        ! New JKQC
        TRAMc = TRAMc + 16d-6*dble(9*Rnz*nfreq)

      end if ! MPI/serial

      ! Master with PRD iterations, JKQC previous iteration
      if (pid.eq.0.and.PRD.and.Input%iter_prd.gt.1) &
        TRAMc = TRAMc + 16d-6*dble(9*nfreq*Rnz)

      end subroutine solve_predict

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize arrays and pointers for the self-consistent
      !! solution\n
      !!          Atom(Atom_class(:)): Structures with atomic data\n
      !!        Frec(Frequency_class): Structure with frequency data\n
      !!         Geom(Geometry_class): Structure with geometric data\n
      !!              MPID(MPI_class): Structure with MPI data\n
      !!           Input(Input_class): Structure with configuration
      !!                               data\n
      !!              NG_dim(integer): Size of NG entry\n
      !!            NG_entry(integer): Index of NG entries\n
      !!      NG_scracth(double(:,:)): Data for NG iteration\n
      !!                 NGP(logical): If doing NG acceleration\n
      !!                doNG(logical): If doing NG iteration\n
      !!               goout(logical): Convergence flag\n
      !!                 ADD(logical): If Stokes needed for PRD-AD\n
      !!               RPRAM(logical): If originally saving redis.
      !!                               functions\n
      !!              lp_exu(logical): If available pre-computed
      !!                               exponentials\n
      !!                 if0(integer): Initial frequency index\n
      !!                 if1(integer): Final frequency index\n
      !!                 npz(integer): Size in z and phi\n
      !!                ntpz(integer): Size in z, theta, and phi\n
      !!               nsend(integer): Expected messages per
      !!                               iteration\n
      !!          Stokes_r(double(:)): Receiver buffer for Stokes\n
      !!            Prof_r(double(:)): Receiver buffer for profiles\n
      !!  Stokes_s(double(:,:,:,:,:)): Sender buffer for Stokes\n
      !!    Prof_s(double(:,:,:,:,:)): Sender buffer for profiles\n
      !!        data1M(double(:,:,:)): RT coeff point M\n
      !!        data1O(double(:,:,:)): RT coeff point O\n
      !!          data2O(double(:,:)): Profiles point O\n
      !!             p_exu(double(:)): Pointer to pre-calculated
      !!                               exponential\n
      !!  JKQC_old(dcomplex(:,:,:,:)): Old radiation field tensors
      !!                               with frequency dependence
      !!    JKQC_n(dcomplex(:,:,:,:)): New radiation field tensors
      !!                               with frequency dependence
      subroutine solve_init(Atom,Frec,Geom,MPID,Input,NG_dim, &
                            NG_entry,NG_scratch,NGP,doNG,goout,ADD, &
                            RPRAM,lp_exu,if0,if1,npz,ntpz,nsend, &
                            Stokes_r,Prof_r,Stokes_s,Prof_s, &
                            data1M,data1O,data2O,p_exu, &
                            JKQC_old,JKQC_n)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      logical, intent(out):: NGP,doNG,goout,ADD,RPRAM,lp_exu
      integer, intent(out):: NG_dim,NG_entry,if0,if1,npz,ntpz,nsend
      double precision, dimension(:,:), &
                        allocatable, intent(out):: NG_scratch
      double precision, dimension(:), &
                        allocatable, intent(out):: Stokes_r
      double precision, dimension(:), &
                        allocatable, intent(out):: Prof_r
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Stokes_s
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Prof_s
      double precision, dimension(:,:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data2O
      double precision, dimension(:), &
                        pointer, intent(inout):: p_exu
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQC_old
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQC_n

      ! Local

      logical:: AD,laux,cfile

      integer:: ia,iaux,ios


      ! Initialize converged flag
      goout = .False.

      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! If storing redistribution
      RPRAM = PRAM

      ! Trick to have AV input for AD calculation
      if (tbAD) then
        AV = .True.
        PRAM = .False.
      end if

      !
      ! Initialize NG quantities
      !

      ! Initialize entry index
      NG_entry = 0
      doNG = .False.
      NGP = Input%NG

      ! If NG acceleration
      if (Input%NG) then

        ! Initialize NG rho00 dimension
        NG_dim = 0

        ! For each atom
        do ia=1,nA
          NG_dim = NG_dim + Rnz*Atom(ia)%nlevel
        end do

        ! If doing PRD
        if (PRD) then

          ! If we need Stokes
          if (ADD) then

            NG_dim = NG_dim + nfreq*Geom%nPh*Geom%nTh*(giz1-giz0+1)

          ! If we need J00C
          else

            NG_dim = NG_dim + nfreq*RnZ

          end if !ADD
        end if ! PRD

        ! If it requires too much buffer
        if (dble(NG_dim)*8d-6.gt.maxbuffer_NG) then

          ! No NG acceleration
          NGP = .False.

          ! Issue warning
          umsg = ' # The buffer for NG acceleration '// &
                 'is too big. Not doing NG.'
          call verbose

        end if ! NG space too large

        ! If finally doing it, allocate
        if (NGP) then

          ! Master
          if (pid.eq.0) then
            allocate(NG_scratch(NG_dim, Input%NG_ord+2))
          else
            allocate(NG_scratch(NG_dim, 1))
          end if

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
          if (MPID%alternP) then

            ! MPI packages
            nsend = MPID%nnd*Geom%nTh*Geom%nPh

            ! To receive Intensity chunks
            iaux = 4*MPID%nxfreq*Rnz
            allocate(Stokes_r(iaux))

            ! To receive profile information
            iaux = MPID%nxtfreq*2*Rnz
            allocate(Prof_r(iaux))

          ! Normal
          else

            ! MPI packages
            nsend = MPID%nnd

            ! Dimensions
            npz = Geom%nph*Rnz
            ntpz = Geom%nth*npz

            ! To receive Intensity chunks
            iaux = 4*MPID%nxfreq*Geom%nTh*Geom%nPh*Rnz
            allocate(Stokes_r(iaux))

            ! To receive profile information
            iaux = MPID%nxtfreq*2*Geom%nTh*Geom%nPh*Rnz
            allocate(Prof_r(iaux))

          end if ! Type of MPI

        ! Slave
        else

          ! Allocate M and O pointers for RT coeff
          allocate(data1M(4,MPID%nf(pid),6))
          allocate(data1O(4,MPID%nf(pid),6))

          ! Allocate O pointers
          allocate(data2O(Frec%ntfreq,2))

          ! Alternative MPI
          if (MPID%alternP) then

            ! To send Intensity chunks
            allocate(Stokes_s(0:3,if0:if1,Rz0:Rz1,1,1))

            ! To send profile information
            allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1,1,1))

          ! Normal MPI
          else

            ! To send Intensity chunks
            allocate(Stokes_s(0:3,if0:if1,Rz0:Rz1,Geom%nPh,Geom%nTh))

            ! To send profile information
            allocate(Prof_s(Frec%ntfreq,2,Rz0:Rz1,Geom%nPh,Geom%nTh))

          end if ! Type of MPI
        end if ! Master/slave

      ! Serial
      else

        ! CPU limits
        if0 = 1
        if1 = nfreq

        ! Common (Master and slave)
        ! Allocate O pointers
        allocate(data2O(Frec%ntfreq,2))

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(4,nfreq,6))
        allocate(data1O(4,nfreq,6))

        ! Initialize exponential is not in memory
        if (PIRAM.and.Frec%pif1.ge.Frec%pif0) then
          lp_exu = .True.
        ! No exu allocated or needed
        else
          lp_exu = .False.
          allocate(p_exu(1))
        end if

        ! Allocate new JKQ
        allocate(JKQC_n(0:2,0:2,nfreq,Rz0:Rz1))

      end if ! MPI/serial

      ! Master with PRD iterations, JKQC previous iteration
      if (pid.eq.0.and.PRD.and.Input%iter_prd.gt.1) &
        allocate(JKQC_old(0:2,0:2,nfreq,Rz0:Rz1))


      !
      ! Initialization messages
      !

      ! Global Master
      if (gpid.eq.0) then

        ! Announce we are starting
        umsg = '   Iteration          MRC(rho^0_0) Atom_index '// &
               'Term_index    2J  Height_index Height(km) |   '// &
               '       MRC(rho^K_Q) Atom_index Term_index '// &
               "   2J   2J' Height_index Height(km)  K  Q"
        call verbose
      end if

      ! Open the file to store MRC
      if (gpid.eq.0.and.Input%keep_MRC) then

        ! If appending
        if (Input%appendMRC) then

          ! If does not exist, create it
          inquire(file=trim(Input%folder)//'/MRC', exist=laux)

          ! Create file is does not exist
          cfile = .not.laux

        ! Not appending
        else

          ! Create file
          cfile = .True.

        end if

        ! Creating MRC file
        if (cfile) then

          ! Open new file
          open(800, file=trim(Input%folder)//'/MRC', &
               action='write',iostat=ios,err=1000)

          ! Write header
          write(800,'(A)',err=1100) &
                       '!  Iteration          MRC(rho^0_0) '// &
                       'Atom_index Term_index    2J  '// &
                       'Height_index Height(km) |   '// &
                       '       MRC(rho^K_Q) Atom_index '// &
                       'Term_index    2J   2J'// &
                       "' Height_index Height(km)  K  Q"

          ! Close file
          close(800)

        end if ! Creating file
      end if ! Global master keeping MRC

      return

1000  umsg = 'Error opening MRC file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing MRC file'
      close(800)
      return

      end subroutine solve_init

!#####################################################################
!#####################################################################
!#####################################################################

      !>    Parallel-master task for the self-consistent problem\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!           Atmo(Atmo_class): Structure with atmospheric
      !!                             data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!             nsend(integer): Expected messages per iteration\n
      !!              iter(integer): Iteration number\n
      !!               npz(integer): Size in z and phi\n
      !!              ntpz(integer): Size in theta and phi\n
      !!        Stokes_r(double(:)): Receiver buffer for Stokes\n
      !!          Prof_r(double(:)): Receiver buffer for profiles\n
      !!  JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                             tensors\n
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
      subroutine solve_manager(Atom,Atmo,Frec,Geom,Bfield,MPID, &
                               Input,Flgsg,nsend,iter,npz,ntpz, &
                               Stokes_r,Prof_r,JKQ_asym,Stokes, &
                               JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      integer, intent(in):: nsend,iter,npz,ntpz
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: Stokes_r
      double precision, dimension(:), &
                        allocatable, target, intent(inout):: Prof_r
      double precision, &
             dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             target, intent(out):: Stokes
      double precision, dimension(nxphot,2,Rz0:Rz1), &
                        intent(out):: J00P
      complex(kind=8), dimension(:,:,:), intent(in):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(out):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(out):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(out):: JKQC

      ! Local

      logical:: deal

      integer:: ith,iph,jdir,iz,itpz,if0l,if1l,ip0l,ip1l,nfl,nftl
      integer:: ia,itran,jtran,K,iQ,ierr,id,info_b,if0p
      integer, dimension(3):: info_c

      double precision:: WA,daux,daux2
      double precision, &
             dimension(2,nxtran,Geom%njdir,Rz0:Rz1):: Norm
      complex(kind=8), &
         dimension(0:2,0:2,2,nxtran,Geom%njdir,Rz0:Rz1):: BStk

      ! Pointers

      double precision, dimension(:), pointer:: p_exu
      double precision, dimension(:,:,:), pointer:: p_MStk
      double precision, dimension(:,:,:), pointer:: p_MProf
      complex(kind=8), dimension(:,:,:), pointer:: TSo,TKQo


      ! Nullify pointers
      nullify(p_exu,p_MStk,p_MProf,TSo,TKQo)

      ! Reset radiation field variables
      JKQ = cZero
      JKQS = cZero
      J00P = 0d0
      Norm = 0d0
      BStk = cZero
      JKQC = cZero

      !
      ! Alternative MPI
      !
      if (MPID%alternP) then

        ! MPI packages
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
          jdir = Geom%i_geom(iph,ith)

          ! Point TKQ_S
          TSo => Geom%TSo(:,:,:,jdir)

          ! Flag error
          if (info_b.lt.0) laborted = .True.

          ! Continue?
          if (info_b.lt.0) cycle

          ! Receive intensity
          do while (.True.)
            call MPI_recv(Stokes_r(1), MPID%size4(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), MPID%size5(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          1+info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! If measuring performance
          if (Input%mpi_perf) &
            call report_mpi_time(Input%folder,Input%ID, &
                                 info_b,iter,.True.)

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          ip0l = Frec%Mpif0(info_b)
          ip1l = Frec%Mpif1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreq(info_b)
          if0p = MPID%if0(info_b) - 1

          ! Pointers
          p_MStk(0:3,if0l:if1l,Rz0:Rz1) => &
                                 Stokes_r(1:MPID%size4(info_b))
          p_MProf(1:nftl,1:2,Rz0:Rz1) => &
                                   Prof_r(1:MPID%size5(info_b))

          ! Do not process if leaving
          if (laborted) cycle

          ! Get angular weight
          WA = Geom%W_mu(ith)*Geom%W_mux(iph)

          ! Initialize deal flag
          deal = .False.

          ! Each height
          do iz=Rz0,Rz1

            ! Determine where to store intensity
            if (KSTK.or.iz.eq.Rz0) &
              Stokes(0:3,if0l:if1l,iph,ith,iz) = &
                                        p_MStk(0:3,if0l:if1l,iz)
            ! Point to exu values
            if (PIRAM.and.ip1l.ge.ip0l) then
              p_exu => Frec%exu(ip0l:ip1l,iz)
            else
              allocate(p_exu(1))
              deal = .True.
            end if

            ! Select correct TKQout
            if (Bfield%Bstrength(iz).gt.TINYB) then
              TKQo => Geom%TBo(:,:,:,jdir,iz)
            else
              TKQo => TSo
            end if

            ! Calculate frequency integral for b-b quantities
            call FInt_line(Atom,Frec%W_freq, &
                           Frec%Mlif0(info_b), &
                           Frec%Mlif1(info_b), &
                           if0p,info_b,p_MStk(:,:,iz), &
                           p_MProf(:,:,iz),TKQo, &
                           Norm(:,:,jdir,iz), &
                           BStk(:,:,:,:,jdir,iz))

            !
            ! Calculate rest of integrals
            !
            call FInt_rest(Atom,Frec%omega, &
                           Frec%W_freq,ip0l,ip1l,if0l,if1l,if0p, &
                           Atmo%T(iz),info_b,WA, &
                           p_MStk(:,:,iz),TSo,J00P(:,:,iz), &
                           JKQC(:,:,:,iz),p_exu)

            ! Nullify pointer
            if (deal) deallocate(p_exu)
            nullify(p_exu)

          end do ! heights
        end do ! MPI packages

      !
      ! Normal MPI
      !
      else

        ! Each MPI communication
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
            call MPI_recv(Stokes_r(1), MPID%size4(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! Receive profile
          do while (.True.)
            call MPI_recv(Prof_r(1), MPID%size5(info_b), &
                          MPI_DOUBLE_PRECISION, info_b, &
                          1+info_b, MPI_COMM_RT, &
                          MPI_STATUS_IGNORE, ierr)
            if (ierr.eq.0) exit
          end do

          ! If measuring performance
          if (Input%mpi_perf) &
            call report_mpi_time(Input%folder,Input%ID, &
                                 info_b,iter,.True.)

          ! Shorter variables
          if0l = MPID%if0(info_b)
          if1l = MPID%if1(info_b)
          ip0l = Frec%Mpif0(info_b)
          ip1l = Frec%Mpif1(info_b)
          nfl = MPID%nf(info_b)
          nftl = Frec%Mntfreq(info_b)
          if0p = MPID%if0(info_b) - 1

          ! Pointers
          p_MStk(0:3,if0l:if1l,1:ntpz) => &
                                     Stokes_r(1:MPID%size4(info_b))
          p_MProf(1:nftl,1:2,1:ntpz) => Prof_r(1:MPID%size5(info_b))

          ! Do not process if leaving
          if (laborted) cycle

          ! Initialize deal flag
          deal = .False.

          ! Compute line quantities
          do itpz=1,ntpz

            ! Get indexes
            ith = (itpz-1)/npz
            iph = (itpz - npz*ith - 1)/Rnz
            iz = itpz - Rnz*iph - npz*ith + Rz0 - 1
            ith = ith + 1
            iph = iph + 1
            jdir = Geom%i_geom(iph,ith)

            ! Point TKQ_S
            TSo => Geom%TSo(:,:,:,jdir)

            ! Select correct TKQout
            if (Bfield%Bstrength(iz).gt.TINYB) then
              TKQo => Geom%TBo(:,:,:,jdir,iz)
            else
              TKQo => TSo
            end if

            ! Determine where to store intensity
            if (KSTK.or.iz.eq.Rz0) &
              Stokes(0:3,if0l:if1l,iph,ith,iz) = &
                                          p_MStk(0:3,if0l:if1l,itpz)

            ! Calculate frequency integral for b-b quantities
            call FInt_line(Atom,Frec%W_freq, &
                           Frec%Mlif0(info_b),Frec%Mlif1(info_b), &
                           if0p,info_b,p_MStk(:,:,itpz), &
                           p_MProf(:,:,itpz),TKQo, &
                           Norm(:,:,jdir,iz), &
                           BStk(:,:,:,:,jdir,iz))

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
            call FInt_rest(Atom,Frec%omega, &
                           Frec%W_freq,ip0l,ip1l,if0l,if1l,if0p, &
                           Atmo%T(iz),info_b,WA, &
                           p_MStk(:,:,itpz),TSo,J00P(:,:,iz), &
                           JKQC(:,:,:,iz),p_exu)

            ! Nullify pointer
            if (deal) deallocate(p_exu)
            nullify(p_exu)

          end do ! Directions and heights
        end do ! frequency domains

      end if ! Type of MPI

      ! Nullify pointers
      nullify(p_MStk,p_MProf)

      !
      ! Apply weights to JKQ, JKQS, and normalize
      !

      ! For each height
      do iz=Rz0,Rz1

        ! Initialize
        jdir = 0

        ! For each polar direction
        do ith=1,Geom%nTh

          ! For each azimuthal direction
          do iph=1,Geom%nph

            ! Advance direction index
            jdir = jdir + 1

            ! Get the angular integral weight
            WA = Geom%W_mu(ith)*Geom%W_mux(iph)

            ! For each atom
            do ia=1,nA

              ! For each FS transition
              do itran=1,Atom(ia)%ntran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tshift

                ! Get the weight
                if (Norm(1,jtran,jdir,iz).gt.0d0) then

                  ! Inverse norm and angular weight
                  daux = WA/Norm(1,jtran,jdir,iz)

                  ! Integrate angle
                  JKQ(0:2,:,jtran,iz) = JKQ(0:2,:,jtran,iz) + &
                                 BStk(:,:,1,jtran,jdir,iz)*daux

                end if

                ! If there is stimulated emission
                if (stm) then

                  ! Get the weight
                  if (Norm(2,jtran,jdir,iz).gt.0d0) then

                    ! Inverse norm and angular weight
                    daux = WA/Norm(2,jtran,jdir,iz)

                    ! Integrate angle
                    JKQS(0:2,:,jtran,iz) = &
                                 JKQS(0:2,:,jtran,iz) + &
                                 BStk(:,:,2,jtran,jdir,iz)*daux

                  end if

                end if ! Stimulated emission
              end do ! transitions
            end do ! atoms
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

            ! Apply atomic shift
            jtran = itran + Atom(ia)%pshift

            ! Calculate the multiplicative factor
            daux2 = daux*exp(Atom(ia)%phot(itran)%edge*WA)* &
                    Atom(ia)%phot(itran)%glu

            ! Apply it to the emission integral
            J00P(jtran,2,iz) = J00P(jtran,2,iz)*daux2

          end do ! b-f transitions
        end do ! atoms
      end do ! heights

      ! If not axial, need to complete JKQ
      if (.not.axial) then

        ! Continuum
        do K=1,Krad
          do iQ=1,K
            JKQC(-iQ,K,:,:) = Flgsg%sg(iQ)*conjg(JKQC(iQ,K,:,:))
          end do
        end do

        ! For each atom
        do ia=1,nA

          ! For each FS transition
          do itran=1,Atom(ia)%ntran

            ! Apply atomic shift
            jtran = itran + Atom(ia)%tshift

            ! K
            do K=1,Atom(ia)%Krad(itran)

              ! Q
              do iQ=1,K

                ! Conjugate
                JKQ(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                     conjg(JKQ(iQ,K,jtran,:))

                ! Stimulated, conjugate
                if (stm) &
                  JKQS(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                        conjg(JKQS(iQ,K,jtran,:))

              end do ! Q
            end do ! K
          end do ! Transition
        end do ! Atom
      end if ! Not axial

      !
      ! Add the Ad-hoc asymmetries
      !
      if (Input%nasym.gt.0) &
        call addJKQasym(Bfield,Flgsg,JKQ_asym,JKQ,JKQS,JKQC)

      ! Clean pointers
      nullify(TSo,TKQo)

      end subroutine solve_manager

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the RTE for the angular quadrature\n
      !!          Atom(Atom_class(:)): Structures with atomic data\n
      !!   LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!             Atmo(Atmo_class): Structure with atmospheric
      !!                               data\n
      !!        Cont(Continuum_class): Structure with background
      !!                               opacity data\n
      !!        Frec(Frequency_class): Structure with frequency data\n
      !!               Red(Red_class): Structure with redistribution
      !!                               input frequency data,
      !!                               redistribution function data,
      !!                               and profile or normalization
      !!                               data\n
      !!         Bfield(Bfield_class): Structure with magnetic field
      !!                               data\n
      !!         Geom(Geometry_class): Structure with geometric data\n
      !!              MPID(MPI_class): Structure with MPI data\n
      !!           Input(Input_class): Structure with configuration
      !!                               data\n
      !!           Flgsg(Fctsg_class): Structure with factorials,
      !!                               signs, and J-symbols\n
      !!              lp_exu(logical): If available pre-computed
      !!                               exponentials\n
      !!                 if0(integer): Initial frequency index\n
      !!                 if1(integer): Final frequency index\n
      !!  Stokes_s(double(:,:,:,:,:)): Sender buffer for Stokes\n
      !!    Prof_s(double(:,:,:,:,:)): Sender buffer for profiles\n
      !!        data1M(double(:,:,:)): RT coeff point M\n
      !!        data1O(double(:,:,:)): RT coeff point O\n
      !!          data2O(double(:,:)): Profiles point O\n
      !!             p_exu(double(:)): Pointer to pre-calculated
      !!                               exponential\n
      !!    JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the
      !!                               radiation tensors\n
      !!    Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!       JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                               integrated over the absorption
      !!                               profile\n
      !!      JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                               integrated over the emission
      !!                               profile\n
      !!      JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                               frequency dependence\n
      !!    JKQC_n(dcomplex(:,:,:,:)): New radiation field tensors
      !!                               with frequency dependence\n
      !!          J00P(double(:,:,:)): Intensity integrals in the
      !!                               photoionization rates
      subroutine solve_RT(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                          Geom,MPID,Input,Flgsg,lp_exu,if0,if1, &
                          Stokes_s,Prof_s,data1M,data1O,data2O, &
                          p_exu,JKQ_asym,Stokes,JKQ,JKQS,JKQC, &
                          JKQC_n,J00P)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Bfield_class), intent(in):: Bfield
      type(Geometry_class), intent(in):: Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Fctsg_class), intent(inout):: Flgsg
      logical, intent(in):: lp_exu
      integer, intent(in):: if0,if1
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Stokes_s
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Prof_s
      double precision, &
          dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
          intent(out):: Stokes
      double precision, dimension(nxphot,2,Rz0:Rz1), &
                       intent(out):: J00P
      complex(kind=8), dimension(:,:,:):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(out):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(out):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(inout):: JKQC
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQC_n
      double precision, dimension(:,:,:), &
                        pointer, intent(inout):: data1M,data1O
      double precision, dimension(:,:), &
                        pointer, intent(inout):: data2O
      double precision, dimension(:), &
                        pointer, intent(inout):: p_exu

      ! Local

      integer:: ith,iph,jth,jph,jdir,mfreq,iz,iz0,iz1,diz,m,o,p
      integer:: ia,itran,jtran,K,iQ
      integer, dimension(3):: info_c

      double precision:: vfac,mu_inv,dsm,dsp,ct,st,cc,sc,WA

      ! Pointers

      double precision, dimension(:,:), pointer:: data2P
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP
      double precision, dimension(:,:,:), pointer:: data1P
      complex(kind=8), dimension(:,:,:), pointer:: TSo,TKQo


      ! Nullify pointers
      nullify(p_K0M,p_K1M,p_K2M,p_SM,p_StkM)
      nullify(p_K0O,p_K1O,p_K2O,p_SO,p_StkO)
      nullify(p_K0P,p_SP)
      nullify(data1P,data2P)
      nullify(TSo,TKQo)

      ! Initialize in serial
      if (pid.eq.0) then

        ! Reset radiation field variables
        Stokes = 0d0
        JKQ = cZero
        JKQS = cZero
        JKQC_n = cZero
        J00P = 0d0

      end if

      ! Frequency size
      mfreq = if1 - if0 + 1

      ! Initialize Doppler shift
      vfac = 1d0

      !
      ! Radiation transfer
      !

      ! For each polar direction
      do ith=1,Geom%nTh

        ! Error
        if (laborted) exit

        ! Calculate inverse of cosine of polar direction
        mu_inv = 1d0/Geom%V_mu(ith)

        ! Determine the direction of propagation for indexes
        diz = -int(sign(1d0, Geom%V_mu(ith)))

        ! Determine the first and last height indexes to run over
        iz0 = ((Rz0 - Rz1)*diz + Rz1 + Rz0)/2
        iz1 = ((Rz1 - Rz0)*diz + Rz1 + Rz0)/2

        ! Trigonometry for Doppler shift
        if (dyn) then
          ct = Geom%V_mu(ith)
          st = sqrt(1d0 - ct*ct)
        end if

        ! For each azimuthal direction
        do iph=1,Geom%nPh

          ! Error
          if (laborted) exit

          ! Get direction index
          jdir = Geom%i_geom(iph,ith)

          ! Point TKQ_S
          TSo => Geom%TSo(:,:,:,jdir)

          ! Alternative MPI
          if (MPID%alternP) then
            jth = 1
            jph = 1
          else
            jth = ith
            jph = iph
          end if

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
            call top(if0,if1,data1M(:,:,6))

          ! If going up, get bottom boundary
          else

            ! Call bottom boundary
            call bottom(Frec%omega,Atmo%T(iz0),vfac, &
                        if0,if1,data1M(:,:,6))

          endif ! propagation direction

          ! Identify current height
          o = iz0

          ! Select correct TKQout
          if (Bfield%Bstrength(o).gt.TINYB) then
            TKQo => Geom%TBo(:,:,:,jdir,o)
          else
            TKQo => TSo
          end if

          ! Calculate radiative coefficients
          call RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                       Geom,vfac,o,jdir,if0,if1, &
                       JKQC(:,:,:,o),Cont%ndir, &
                       Cont%c(:,:,:,o),Bfield,TSo,TKQo, &
                       data1M(:,:,1:5),data2O,.True.)

          ! Error
          if (laborted) exit

          ! If MPI
          if (pid.gt.0) then

            !
            ! Store in buffer
            !

            ! Stokes
            Stokes_s(:,:,o,jph,jth) = data1M(:,:,6)

            ! Profiles
            Prof_s(:,:,o,jph,jth) = data2O

          ! Serial
          else

            ! Store Stokes
            if (KSTK) Stokes(:,:,iph,ith,o) = data1M(:,:,6)

            ! Point to exu values
            if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

            ! Angular weight
            WA = Geom%W_mu(ith)*Geom%W_mux(iph)

            !
            ! Calculate integrals
            !
            call Jcalc(Atom,Frec%omega,Frec%W_freq, &
                       Frec%lif0,Frec%lif1, &
                       Frec%pif0,Frec%pif1, &
                       Atmo%T(o),Atmo%ne(o),WA, &
                       data1M(:,:,6),data2O,TSo,TKQo, &
                       JKQ(:,:,:,o),JKQS(:,:,:,o),J00P(:,:,o), &
                       JKQC_n(:,:,:,o),p_exu)

          end if

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

          ! Select correct TKQout
          if (Bfield%Bstrength(p).gt.TINYB) then
            TKQo => Geom%TBo(:,:,:,jdir,p)
          else
            TKQo => TSo
          end if

          ! Calculate radiative coefficients
          call RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                       Geom,vfac,p,jdir,if0,if1, &
                       JKQC(:,:,:,p),Cont%ndir, &
                       Cont%c(:,:,:,p),Bfield,TSo,TKQo, &
                       data1O(:,:,1:5),data2O,.True.)

          ! Error
          if (laborted) exit

          !
          ! Intermediate heights
          !

          ! For each height this CPU has assigned
          do iz=iz0,iz1,diz

            ! We treat the boundaries outside
            if(iz.eq.iz0.or.iz.eq.iz1)cycle

            ! Allocate P pointers
            allocate(data1P(4,MPID%nf(pid),6))
            allocate(data2P(Frec%ntfreq,2))

            ! Identify heights
            m = iz - diz
            o = iz
            p = iz + diz

            ! Calculate distance to previous point
            dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

            ! Caculate quantities of the next point
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

            ! Select correct TKQout
            if (Bfield%Bstrength(p).gt.TINYB) then
              TKQo => Geom%TBo(:,:,:,jdir,p)
            else
              TKQo => TSo
            end if

            ! RT coefficients
            call RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                         Geom,vfac,p,jdir,if0,if1, &
                         JKQC(:,:,:,p),Cont%ndir, &
                         Cont%c(:,:,:,p),Bfield,TSo,TKQo, &
                         data1P(:,:,1:5),data2P,.True.)

            ! Point to the data
            p_K0M  => data1M(:,:,1)
            p_K1M  => data1M(:,:,2)
            p_K2M  => data1M(:,:,3)
            p_SM   => data1M(:,:,5)
            p_StkM => data1M(:,:,6)
            p_K0O  => data1O(:,:,1)
            p_K1O  => data1O(:,:,2)
            p_K2O  => data1O(:,:,3)
            p_SO   => data1O(:,:,5)
            p_StkO => data1O(:,:,6)
            p_K0P  => data1P(:,:,1)
            p_SP   => data1P(:,:,5)

            ! Apply short characteristics BESSER
            call RTStep(o,ith,iph,mfreq, &
                        dsm,dsp,p_K0M,p_K1M,p_K2M, &
                        p_SM,p_K0O,p_K1O,p_K2O, &
                        p_SO,p_K0P,p_SP,p_StkM, &
                        p_StkO,.True.)

            ! If MPI
            if (pid.gt.0) then

              !
              ! Store in buffer
              !

              ! Stokes
              Stokes_s(:,:,o,jph,jth) = data1O(:,:,6)

              ! Profiles
              Prof_s(:,:,o,jph,jth) = data2O

            ! Serial
            else

              ! Store Stokes
              if (KSTK) Stokes(:,:,iph,ith,o) = data1O(:,:,6)

              ! Point to exu values
              if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

              ! Select correct TKQout
              if (Bfield%Bstrength(o).gt.TINYB) then
                TKQo => Geom%TBo(:,:,:,jdir,o)
              else
                TKQo => TSo
              end if

              ! Calculate contribution to JKQ
              call Jcalc(Atom,Frec%omega,Frec%W_freq, &
                         Frec%lif0,Frec%lif1, &
                         Frec%pif0,Frec%pif1, &
                         Atmo%T(o),Atmo%ne(o),WA, &
                         data1O(:,:,6),data2O,TSo,TKQo, &
                         JKQ(:,:,:,o),JKQS(:,:,:,o), &
                         J00P(:,:,o),JKQC_n(:,:,:,o),p_exu)

            end if ! MPI/serial

            ! Shift data (O->M, P->O)
            deallocate(data1M,data2O)
            data1M => data1O
            data1O => data1P
            data2O => data2P
            nullify(data1P,data2P)

            ! Error
            if (laborted) exit

          end do ! Intermediate heights

          ! Error
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
          p_K0M  => data1M(:,:,1)
          p_K1M  => data1M(:,:,2)
          p_K2M  => data1M(:,:,3)
          p_SM   => data1M(:,:,5)
          p_StkM => data1M(:,:,6)
          p_K0O  => data1O(:,:,1)
          p_K1O  => data1O(:,:,2)
          p_K2O  => data1O(:,:,3)
          p_SO   => data1O(:,:,5)
          p_StkO => data1O(:,:,6)

          ! Apply short characteristics LINEAR
          call RTStep(o,ith,iph,mfreq, &
                      dsm,dsp,p_K0M,p_K1M,p_K2M, &
                      p_SM,p_K0O,p_K1O,p_K2O, &
                      p_SO,p_K0P,p_SP,p_StkM, &
                      p_StkO,.False.)

          ! Error
          if (laborted) exit

          ! If MPI
          if (pid.gt.0) then

            !
            ! Store in buffer
            !

            ! Stokes
            Stokes_s(:,:,o,jph,jth) = data1O(:,:,6)

            ! Profiles
            Prof_s(:,:,o,jph,jth) = data2O

            !
            ! Alternative
            !
            if (MPID%alternP) then

              !
              ! Send to master if error
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
                call MPI_SEND(Stokes_s(0,if0,Rz0,1,1), &
                              MPID%size4(pid), &
                              MPI_DOUBLE_PRECISION, 0, pid, &
                              MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

              ! Send profiles
              do while (.True.)
                call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                              MPID%size5(pid), &
                              MPI_DOUBLE_PRECISION, &
                              0, 1+pid, MPI_COMM_RT, ierr)
                if (ierr.eq.0) exit
              end do

            end if ! Alternative MPI

          ! Serial
          else

            ! Store Stokes
            if (KSTK) Stokes(:,:,iph,ith,o) = data1O(:,:,6)

            ! Point to exu values
            if (lp_exu) p_exu => Frec%exu(Frec%pif0:Frec%pif1,o)

            ! Select correct TKQout
            if (Bfield%Bstrength(o).gt.TINYB) then
              TKQo => Geom%TBo(:,:,:,jdir,o)
            else
              TKQo => TSo
            end if

            !
            ! Calculate integrals
            !
            call Jcalc(Atom,Frec%omega,Frec%W_freq, &
                       Frec%lif0,Frec%lif1, &
                       Frec%pif0,Frec%pif1, &
                       Atmo%T(o),Atmo%ne(o),WA, &
                       data1O(:,:,6),data2O,TSo,TKQo, &
                       JKQ(:,:,:,o),JKQS(:,:,:,o),J00P(:,:,o), &
                       JKQC_n(:,:,:,o),p_exu)

          end if ! MPI/serial
        end do ! Azimuthal angles
      end do ! Polar angles

      ! Error, jump to end
      if (laborted) goto 2000

      ! If MPI not alternative
      if (pid.gt.0.and..not.MPID%alternP) then

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
          call MPI_SEND(Stokes_s(0,if0,Rz0,1,1), &
                        MPID%size4(pid), &
                        MPI_DOUBLE_PRECISION, 0, pid, &
                        MPI_COMM_RT, ierr)
          if (ierr.eq.0) exit
        end do

        ! Send profiles
        do while (.True.)
          call MPI_SEND(Prof_s(1,1,Rz0,1,1), &
                        MPID%size5(pid), MPI_DOUBLE_PRECISION, &
                        0, 1+pid, MPI_COMM_RT, ierr)
          if (ierr.eq.0) exit
        end do

      ! Serial
      else if (pid.eq.0) then

        ! Update variables
        JKQC(0:2,:,:,:) = JKQC_n

        ! Complete JKQ
        if (.not.axial) then

          ! Continuum
          do K=1,Krad
            do iQ=1,K
              JKQC(-iQ,K,:,:) = Flgsg%sg(iQ)*conjg(JKQC(iQ,K,:,:))
            end do
          end do

          ! For each atom
          do ia=1,nA

            ! For each FS transition
            do itran=1,Atom(ia)%ntran

              ! Apply atomic shift
              jtran = itran + Atom(ia)%tshift

              ! K
              do K=1,Atom(ia)%Krad(itran)

                ! Q
                do iQ=1,K

                  ! Conjugate
                  JKQ(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                       conjg(JKQ(iQ,K,jtran,:))

                  ! Stimulated, conjugate
                  if (stm) &
                    JKQS(-iQ,K,jtran,:) = Flgsg%sg(iQ)* &
                                         conjg(JKQS(iQ,K,jtran,:))

                end do ! Q
              end do ! K
            end do ! Transition
          end do ! Atom
        end if ! Not axial

        !
        ! Add the Ad-hoc asymmetries
        !
        if (Input%nasym.gt.0) &
          call addJKQasym(Bfield,Flgsg,JKQ_asym,JKQ,JKQS,JKQC)

      end if ! Normal MPI/serial

      ! Clean pointers
2000  nullify(p_K0M,p_K1M,p_K2M,p_SM,p_StkM)
      nullify(p_K0O,p_K1O,p_K2O,p_SO,p_StkO)
      nullify(p_K0P,p_SP)
      nullify(TSo,TKQo)

      ! If aborting MPI, deal with messages
      if (laborted.and.pid.gt.0) then

        ! Alternative MPI
        if (MPID%alternP) then

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
      end if ! Error while doing MPI

      end subroutine solve_RT

!#####################################################################
!#####################################################################
!#####################################################################

      !> Advance the density matrices by solving the SEE and applying
      !! NG acceleration if requested\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!     Rho_old(Rhoc_class(:)): Structure to store the density
      !!                             matrix of the previous
      !!                             iteration\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!             MRC(MRC_class): Structure with the Maximum
      !!                             Relative Change data\n
      !!            NG_dim(integer): Size of NG entry\n
      !!          NG_entry(integer): Index of NG entries\n
      !!    NG_scracth(double(:,:)): Data for NG iteration\n
      !!              doNG(logical): If doing NG iteration\n
      !!             goout(logical): If converged\n
      !!               ADD(logical): If Stokes needed for PRD-AD\n
      !!              iter(integer): Current iteration\n
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
      subroutine solve_SEE(Atom,Rho_old,Atmo,Bfield,Geom,Input, &
                           Flgsg,MRC,NG_dim,NG_entry,NG_scratch, &
                           NGP,doNG,goout,ADD,iter,Stokes, &
                           JKQ,JKQS,JKQC,J00P)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Rhoc_class), dimension(:), intent(in):: Rho_old
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Geometry_class), intent(in):: Geom
      type(Fctsg_class), intent(inout):: Flgsg
      type(Input_class), intent(in):: Input
      type(MRC_class), intent(inout):: MRC
      logical, intent(inout):: doNG
      logical, intent(in):: ADD,NGP
      logical, intent(out):: goout
      integer, intent(in):: iter,NG_dim
      integer, intent(inout):: NG_entry
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: NG_scratch
      double precision, dimension(nxphot,2,Rz0:Rz1), intent(in):: J00P
      double precision, &
          dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
          intent(inout):: Stokes
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(in):: JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(inout):: JKQC

      ! Local

      integer:: ia,itran,jtran,iphot,jphot,iz,ith,iph,ifreq
      integer:: o,iterm,iJ,p,ing,ios

      double precision:: larmor


#ifdef DEBUGJKQ
      if (pid.eq.0) call dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS, &
                                  Input%folder,iter)
#endif

      ! For each atom
      do ia=1,nA

        ! Limiting indexes
        itran = Atom(ia)%tshift + 1
        jtran = itran + Atom(ia)%ntran - 1
        iphot = Atom(ia)%pshift + 1
        jphot = iphot + Atom(ia)%nphot - 1

        ! For each height
        do iz=Rz0,Rz1

          ! Set magnetic data
          larmor = B2L*Bfield%Bstrength(iz)

          ! Solve the SEE
          call SEE(Atom(ia),JKQ(:,:,itran:jtran,iz), &
                   JKQS(:,:,itran:jtran,iz), &
                   J00P(iphot:jphot,:,iz),larmor,Flgsg,iz)

        end do ! heights
      end do ! atoms
#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,iter)
#endif

      ! Control
      call control
      if (laborted) goto 2000


      !
      ! NG acceleration
      !

      ! Check if doing NG acceleration
      if(NGP.and.iter.gt.Input%NG_delay)then

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
          if (PRD) then

            ! If ADD
            if (ADD) then

              ! For each height
              do iz=giz0,giz1

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
                                          Stokes(0,ifreq,iph,ith,iz)

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
                  NG_scratch(o,NG_entry) = dble(JKQC(0,0,ifreq,iz))

                end do ! Frequencies
              end do ! Heights

            end if ! AV or AD
          end if ! PRD

          ! Call NG and check if it should be processed
          call NG(NG_dim,p,Input%NGI_ord,NG_scratch,NG_entry,doNG)

        ! Slave
        else

          ! If wrong order
          if (Input%NG_ord.lt.1.or.Input%NG_ord.gt.5) then

            ! Do not do
            doNG = .False.

          ! Valid order
          else

            ! Check if Master is in NG step
            if (NG_entry.gt.(Input%NG_ord+1)) then

              ! Do step
              doNG = .True.

            ! Not a NG step
            else

              ! Not a NG step
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

          ! Share NG iteration
          ! If MPI
          if (nproc.gt.1) &
            call MPI_BCAST(NG_scratch(1,ing), NG_dim, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

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
                p = Atom(ia)%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

                ! For each height
                do iz=Rz0,Rz1

                  ! Advance index
                  o = o + 1

                  ! Accelerate rho
                  Atom(ia)%crho(p,iz) = &
                                      dcmplx(NG_scratch(o,ing), 0d0)
                enddo ! Heights
              enddo ! Levels
            enddo ! Terms
          enddo ! Atoms

          ! If PRD
          if (PRD) then

            ! If ADD
            if (ADD) then

              ! For each height
              do iz=giz0,giz1

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! For each azimuthal direction
                  do iph=1,Geom%nPh

                    ! For each frequency
                    do ifreq=1,nfreq

                      ! Advance index
                      o = o + 1

                      ! Store Stokes
                      Stokes(0,ifreq,iph,ith,iz) = NG_scratch(o,ing)

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
                  JKQC(0,0,ifreq,iz) = &
                                      dcmplx(NG_scratch(o,ing), 0d0)
                end do ! Frequencies
              end do ! Heights

            end if ! AV or AD
          end if ! PRD

          ! Global master
          if (gpid.eq.0) then

            ! Verbose
            umsg = 'NG acceleration'
            call verbose

          end if ! Global master

          ! Reset entry index
          NG_entry = 0

        end if ! Communication
      endif ! NG acceleration

      !
      !  Calculate MRC
      !

      ! Only the master does
      if (pid.eq.0) then

        ! Calculate MRC
        call MRC_sb(Atom,Rho_old,Input%anisotropy_only,MRC)

        ! Convert cm into km
        MRC%values(1,1) = Atmo%z(MRC%indexes(2,1))*1d-5
        MRC%values(1,2) = Atmo%z(MRC%indexes(2,2))*1d-5

        ! Check exit criteria
        if (MRC%values(2,1).le.Input%mrc_i.and. &
            MRC%values(2,2).le.Input%mrc_p) &
          goout = .True.

        ! If doing PRD and just had NG, do not leave
        if (NGP.and.PRD) then
          if (doNG) goout = .False.
        end if

      end if

      ! Only the global Master
      if (gpid.eq.0) then

        ! Write in stdout
        write(umsg, &
        '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
        '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,1x,i4,3x,i11,'// &
        '2x,f9.3,1x,i2,1x,i2)') &
        iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
        MRC%indexes(5,1),MRC%indexes(2,1), &
        MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
        MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
        MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
        MRC%indexes(6,2)
        call verbose

        ! If file
        if (Input%keep_MRC) then

          ! Open old file
          open(800, file=trim(Input%folder)//'/MRC', &
               iostat=ios,err=1000,position='append')

          ! Write in MRC file
          write(800, &
          '(3x,i9,2x,es20.12,2x,i9,2x,i9,2x,i4,3x,i11,'// &
          '2x,f9.3,4x,es20.12,2x,i9,2x,i9,2x,i4,2x,i4,2x,i11,'// &
          '2x,f9.3,1x,i2,1x,i2)', err=1100) &
          iter,MRC%values(2,1),MRC%indexes(1,1),MRC%indexes(3,1), &
          MRC%indexes(5,1),MRC%indexes(2,1), &
          MRC%values(1,1),MRC%values(2,2),MRC%indexes(1,2), &
          MRC%indexes(3,2),MRC%indexes(4,2),MRC%indexes(5,2), &
          MRC%indexes(2,2),MRC%values(1,2),MRC%indexes(6,1), &
          MRC%indexes(6,2)

          ! Close file
          close(800)

        end if ! If MRC in file
      end if ! Global master

      ! Return
2000  return

1000  umsg = 'Error opening MRC file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing MRC file'
      close(800)
      return

      end subroutine solve_SEE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Formal solution for the RTE\n
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
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!   JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the
      !!                              radiation tensors\n
      !!   Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!      JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                              integrated over the absorption
      !!                              profile\n
      !!     JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                              frequency dependence\n
      !!      SolF(Solution_F_class): Structure with the solution of
      !!                              the self-consistent problem and
      !!                              the corresponding emergent
      !!                              profiles, contribution function,
      !!                              and height for optical depth
      !!                              equal to one
      subroutine emergent(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                          Geom,MPID,Input,Flgsg,JKQ_asym,Stokes, &
                          JKQ,JKQC,SolF)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Fctsg_class), intent(inout):: Flgsg
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      type(Geometry_class), intent(inout):: Geom
      type(Bfield_class), intent(in):: Bfield
      type(Solution_F_class), intent(inout):: SolF
      double precision, &
             dimension(0:3,nfreq,Geom%nPh,Geom%nTh,giz0:giz1), &
             target, intent(in):: Stokes
      complex(kind=8), dimension(:,:,:), intent(in):: JKQ_asym
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(in):: JKQC

      ! Local

      logical:: ADD,doRT,doMG,l1

      integer:: tau1size,if0,if1,ith,iph,icount,ncount

      double precision, dimension(:,:), allocatable, target:: tau1
      double precision, dimension(:,:), allocatable, target:: tau
      double precision, dimension(:,:,:), allocatable:: ContrG

      ! Receivers

      double precision, dimension(:), allocatable:: Stokes_r
      double precision, dimension(:), allocatable:: Contr_r

      ! Senders

      double precision, dimension(:,:), allocatable:: Stokes_s
      double precision, dimension(:,:), allocatable:: tau1_s
      double precision, dimension(:,:,:), allocatable:: Contr_s

      ! Pointers

      double precision, dimension(:), pointer:: etaIM
      double precision, dimension(:,:,:), pointer:: data1M,data1O


      ! Nullify
      nullify(data1M,data1O,etaIM)

      ! Initialize arrays
      call emergent_init(Geom,MPID,Input,SolF,ADD,tau1size, &
                         if0,if1,Stokes_r,Contr_r,Stokes_s, &
                         Contr_s,tau1,ContrG,tau,tau1_s, &
                         etaIM,data1M,data1O)

      ! Predict PRD
      if (PRD) call initialize_emiss(Atom,Geom,Red)

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

      ! For each LOS polar direction
      do ith=1,Geom%nThLOS

        ! For each LOS azimuthal direction
        do iph=1,Geom%nPhLOS

          ! Communicate which direction we are doing if global Master
          if (gpid.eq.0) then
            icount = icount + 1
            write(umsg,'(A,i4,A,i4)') &
                       '   Doing direction ',icount,' of ',ncount
            call verbose

          end if

          ! Set geometrical tensors
          call setTKQLOS(Geom,Flgsg,Bfield,ith,iph)

          ! If dynamic, normalize profiles
          if (dyn) &
            call normalize(Atom,LTElines,Atmo,Bfield%Bstrength,Geom, &
                           MPID,Frec,Red,Flgsg,Input%folder, &
                           l1,ith,iph,.True.)

          ! If PRD
          if (PRD) then

              ! If AD, get geometry
              if (.not.AV) call get_scattering_los(Geom,ith,iph)

              ! Compute emiss2ord
              call comoving_emiss2ord(Atom,Atmo,Geom,Frec,Red,Flgsg, &
                                      Bfield,Stokes,JKQ_asym, &
                                      JKQ,JKQC,ith,iph, &
                                      Input%PRD_int_mode,.True.)
          end if

          ! If manager
          if (doMG) then

            ! Call manager tark
            call emergent_manager(Atmo,Frec,Geom,MPID,Input,SolF, &
                                  ith,iph,Stokes_r,Contr_r,tau1, &
                                  ContrG)

          ! Slave or serial
          else

            ! Call RT
            call emergent_RT(Atom,LTElines,Atmo,Cont,Frec,Red, &
                             Bfield,Geom,MPID,Input,Flgsg,SolF, &
                             ADD,tau1size,if0,if1,ith,iph,Stokes_s, &
                             Contr_s,tau1,tau,tau1_s,etaIM, &
                             data1M,data1O,JKQC)

          end if

          ! Free LOS geometrical tensors
          deallocate(Geom%TSL,Geom%TBL)
          nullify(Geom%TSL,Geom%TBL)

          ! If dynamic, free LOS norms
          if (dyn) call free_norm(Red,.False.)

          ! Failure
          call control
          if (laborted) goto 2000

        end do ! Azimuth
      end do ! Polar

      ! Clean
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

      return

      end subroutine emergent

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the approximate RAM space necessary for the
      !! routines solving the formal solution for the RTE\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!       MPID(MPI_class): Structure with MPI data\n
      !!    Input(Input_class): Structure with configuration data
      subroutine emergent_predict(Atom,Red,Geom,MPID,Input)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Red_class), intent(in):: Red
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input

      ! Local

      integer:: iaux,mfreq


      ! Initialize
      TRAMc = 0d0

      ! Predict PRD
      if (PRD) call predict_emiss(Atom,Geom,Red)

      ! If MPI
      if (MPID%mpi) then

        ! CPU sizes
        mfreq = MPID%if1(pid) - MPID%if0(pid) + 1

        ! Master
        if (pid.eq.0) then

          ! To receive Intensity chunks
          iaux = 4*MPID%nxfreq
          TRAMc = TRAMc + 8d-6*dble(iaux)

          ! If calculating height of tau=1
          if (Input%out_tau1) &
            TRAMc = TRAMc + 8d-6*dble(2*nfreq)

          ! If calculating contribution function
          if (Input%out_contr) then
            iaux = 4*MPID%nxfreq*Rnz
            TRAMc = TRAMc + 8d-6*dble(iaux)
            TRAMc = TRAMc + 8d-6*dble(4*nfreq*Rnz)
          end if

          ! If inverting, need to return the output
          if (run_mode.eq.-1) then

            ! Stokes
            SRAMc = SRAMc + 8d-6*dble(4*nfreq*Geom%nPhLOS* &
                                              Geom%nThLOS)

            ! If keeping height for tau equal 1
            if (Input%out_tau1) &
              SRAMc = SRAMc + 8d-6*dble(4*Geom%nPhLOS* &
                                          Geom%nThLOS* &
                                          Input%lim_tau%nn)

            ! If keeping contribution
            if (Input%out_contr) &
              SRAMc = SRAMc + 8d-6*dble(4*Geom%nPhLOS* &
                                          Geom%nThLOS*Rnz* &
                                          Input%lim_ctr%nn)
          end if ! Inversion

        ! Slave
        else

          ! M and O pointers for RT coeff
          TRAMc = TRAMc + 8d-6*dble(6*mfreq*6)

          ! To send Stokes chunks
          TRAMc = TRAMc + 8d-6*dble(4*mfreq)

          ! If calculating tau 1 or contribution function
          if (input%out_tau1.or.input%out_contr) then

            ! Taus
            TRAMc = TRAMc + 8d-6*dble(mfreq*(Rnz+5))

            ! If calculating contribution function
            if (input%out_contr) &
            TRAMc = TRAMc + 8d-6*dble(mfreq*4*(Rnz+1))

          end if ! tau1 output
        end if ! Master/slave

      ! Serial
      else

        ! M and O pointers for RT coeff
        TRAMc = TRAMc + 8d-6*dble(6*nfreq*6)

        ! If calculating height of tau=1
        if (Input%out_tau1.or.Input%out_contr) &
          TRAMc = TRAMc + 8d-6*dble(nfreq*(Rnz*4+4))

        ! If inverting, need to return the output
        if (run_mode.eq.-1) then

          ! Stokes
          SRAMc = SRAMc + 8d-6*dble(4*nfreq*Geom%nPhLOS* &
                                    Geom%nThLOS)

          ! If keeping height for tau equal 1
          if (Input%out_tau1) &
            SRAMc = SRAMc + 8d-6*dble(4*Geom%nPhLOS* &
                                      Geom%nThLOS* &
                                      Input%lim_tau%nn)

          ! If keeping contribution
          if (Input%out_contr) &
            SRAMc = SRAMc + 8d-6*dble(4*Geom%nPhLOS* &
                                      Geom%nThLOS*Rnz* &
                                      Input%lim_ctr%nn)
        end if ! Inversion
      end if ! MPI/serial

      end subroutine emergent_predict

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize arrays and pointers for the formal solution for
      !! the RTE\n
      !!    Geom(Geometry_class): Structure with geometric data\n
      !!         MPID(MPI_class): Structure with MPI data\n
      !!      Input(Input_class): Structure with configuration data\n
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
      !!   Stokes_s(double(:,:)): Sender Stokes buffer\n
      !!  Contr_s(double(:,:,:)): Sender contribution function
      !!                          buffer\n
      !!       tau1(double(:,:)): Data height tau equal to 1\n
      !!   ContrG(double(:,:,:)): Data contribution function\n
      !!        tau(double(:,:)): Data optical depth\n
      !!     tau1_s(double(:,:)): Height optical depth equal to one
      !!                          sender buffer\n
      !!        etaIM(double(:)): Absorptivity point M\n
      !!     data1M(double(:,:)): RT coeff point M\n
      !!     data1O(double(:,:)): RT coeff point O
      subroutine emergent_init(Geom,MPID,Input,SolF,ADD,tau1size, &
                               if0,if1,Stokes_r,Contr_r,Stokes_s, &
                               Contr_s,tau1,ContrG,tau,tau1_s, &
                               etaIM,data1M,data1O)
      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(out):: ADD
      integer, intent(out):: if0,if1,tau1size
      double precision, dimension(:,:), &
                        allocatable, intent(out):: tau1
      double precision, dimension(:,:), &
                        allocatable, intent(out):: tau
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: ContrG
      double precision, dimension(:), &
                        allocatable, intent(out):: Stokes_r
      double precision, dimension(:), &
                        allocatable, intent(out):: Contr_r
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: Contr_s
      double precision, dimension(:,:), &
                        allocatable, intent(out):: Stokes_s
      double precision, dimension(:,:), &
                        allocatable, intent(out):: tau1_s
      double precision, dimension(:), pointer, intent(inout):: etaIM
      double precision, dimension(:,:,:), pointer, &
                        intent(inout):: data1M,data1O

      ! Local

      logical:: AD

      integer:: iaux


      ! Initialize angle depended flag
      AD = .not.AV
      ADD = AD.or.dyn

      ! If MPI
      if (MPID%mpi) then

        ! CPU limits
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! Master
        if (pid.eq.0) then

          ! To receive Intensity chunks
          iaux = 4*MPID%nxfreq
          allocate(Stokes_r(iaux))

          ! If calculating height of tau=1, allocate
          if (Input%out_tau1) allocate(tau1(2,nfreq))

          ! If calculating contribution function, allocate
          if (Input%out_contr) then
            iaux = 4*MPID%nxfreq*Rnz
            allocate(Contr_r(iaux))
            allocate(ContrG(0:3,nfreq,Rz0:Rz1))
          end if

          ! If inverting, need to return the output
          if (run_mode.eq.-1) then
            allocate(SolF%e_Stk(0:3,nfreq,Geom%nPhLOS,Geom%nThLOS))
            if (Input%out_tau1) &
            allocate(SolF%e_tau1(Input%lim_tau%nn, &
                                 Geom%nPhLOS,Geom%nThLOS))
            if (Input%out_contr) &
            allocate(SolF%e_Ctr(0:3,Input%lim_ctr%nn,Rz0:Rz1, &
                                Geom%nPhLOS,Geom%nThLOS))
          end if ! Inversion

        ! Slave
        else

          ! Allocate M and O pointers for RT coeff
          allocate(data1M(4,MPID%nf(pid),6))
          allocate(data1O(4,MPID%nf(pid),6))

          ! To send Stokes chunks
          allocate(Stokes_s(0:3,if0:if1))

          ! If calculating tau 1 or contribution function, allocate
          if (input%out_tau1.or.input%out_contr) then

            allocate(tau(if0:if1,Rz0:Rz1))
            allocate(tau1(2,if0:if1))
            allocate(tau1_s(2,if0:if1))
            tau1size = MPID%nf(pid)*2
            allocate(etaIM(if0:if1))

            ! If calculating contribution function, allocate
            if (input%out_contr) then

              allocate(Contr_s(0:3,1:MPID%nf(pid),Rz0:Rz1))

            end if ! contribution output

          ! No tau1 or contribution function
          else

            ! Nullify
            nullify(etaIM)

          end if ! tau1 output

        end if ! Master/slave

      ! Serial
      else

        ! CPU limits
        if0 = 1
        if1 = nfreq

        ! Allocate M and O pointers for RT coeff
        allocate(data1M(4,nfreq,6))
        allocate(data1O(4,nfreq,6))

        ! If calculating height of tau=1, allocate
        if (Input%out_tau1.or.Input%out_contr) then

          allocate(tau(nfreq,Rz0:Rz1))
          allocate(etaIM(nfreq))
          allocate(tau1(2,nfreq))
          allocate(Contr_s(0:3,nfreq,Rz0:Rz1))

        ! No tau1 or contribution function
        else

          ! Nullify
          nullify(etaIM)

        end if ! tau1 output

        ! If inverting, need to return the output
        if (run_mode.eq.-1) then

          ! Stokes
          allocate(SolF%e_Stk(0:3,nfreq,Geom%nPhLOS,Geom%nThLOS))

          ! If storing height for optical depth equal to 1
          if (Input%out_tau1) &
            allocate(SolF%e_tau1(Input%lim_tau%nn,Geom%nPhLOS, &
                                 Geom%nThLOS))

          ! If storing contribution function
          if (Input%out_contr) &
            allocate(SolF%e_Ctr(0:3,Input%lim_ctr%nn,Rz0:Rz1, &
                                Geom%nPhLOS,Geom%nThLOS))

        end if ! Inversion
      end if ! MPI/serial

      end subroutine emergent_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Parallel-master task for the formal solution for the RTE\n
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
      !!   ContrG(double(:,:,:)): Data contribution function\n
      subroutine emergent_manager(Atmo,Frec,Geom,MPID,Input,SolF, &
                                  ith,iph,Stokes_r,Contr_r,tau1, &
                                  ContrG)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Input_class), intent(in):: Input
      type(Solution_F_class), intent(inout):: SolF
      integer, intent(in):: ith,iph
      double precision, dimension(:,:,:), &
                        allocatable, intent(inout):: ContrG
      double precision, dimension(:), &
                        allocatable, intent(inout):: Stokes_r
      double precision, dimension(:), &
                        allocatable, intent(inout):: Contr_r
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: tau1

      ! Local

      integer:: id,rpid,ierr,tau1size,sshift,iz,ifreq
      integer, dimension(2):: info_b

      double precision, dimension(0:3,nfreq):: Stokes_out


      ! If calculating height of tau=1
      if (Input%out_tau1) then

        ! For each frequency group
        do id=1,MPID%nnd

          ! Receive ID
          call MPI_recv(rpid,1,MPI_INTEGER, &
                        MPI_ANY_SOURCE, 2, &
                        MPI_COMM_RT, MPI_STATUS_IGNORE, &
                        ierr)

          ! Size of tau package
          tau1size = MPID%nf(rpid)*2

          ! Receive tau
          call MPI_recv(tau1(1,MPID%if0(rpid)), &
                        tau1size, &
                        MPI_DOUBLE_PRECISION, rpid, &
                        3+rpid, MPI_COMM_RT, &
                        MPI_STATUS_IGNORE, ierr)

        end do ! Frequency ranges

        ! If inversion
        if (run_mode.eq.-1) then

          ! Keep tau1
          call settau(SolF%e_tau1(:,iph,ith),tau1(2,:), &
                      Input%lim_tau)

        ! If synthesis
        else

          ! Store the height where tau=1
          call writetau(Input%folder,iph,ith,Frec%omega,Geom, &
                        tau1(2,:),Input%lim_tau)
          call control
          if (laborted) return

        end if ! Inversion/synthesis
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
          call MPI_recv(Contr_r(1), MPID%size3(info_b(1)), &
                        MPI_DOUBLE_PRECISION, info_b(1), &
                        1+info_b(1), MPI_COMM_RT, &
                        MPI_STATUS_IGNORE, ierr)

          ! Reset shift in index
          sshift = 0

          ! For each height
          do iz=Rz0,Rz1

            ! Rearrange the contribution function
            do ifreq=0,MPID%nf(info_b(1))-1

              ! Save contribution function
              ContrG(:,MPID%if0(info_b(1))+ifreq,iz) = &
                                  Contr_r(sshift+4*ifreq+1: &
                                          sshift+4*(ifreq+1))

            end do ! frequencies

            ! Update the shift in the buffer
            sshift = sshift + 4*MPID%nf(info_b(1))

          end do ! heights

        end if ! if contribution function

        ! If it is not a boundary, we are not receiving
        ! anything else
        if (info_b(2).lt.0) cycle

        ! Receive Stokes
        call MPI_recv(Stokes_r(1), MPID%size10(info_b(1)), &
                      MPI_DOUBLE_PRECISION, info_b(1), &
                      info_b(1), MPI_COMM_RT, &
                      MPI_STATUS_IGNORE, ierr)

        ! Rearrange the Stokes
        do ifreq=0,MPID%nf(info_b(1))-1

          ! Save Stokes
          Stokes_out(:,MPID%if0(info_b(1))+ifreq) = &
                              Stokes_r(4*ifreq+1:4*(ifreq+1))

        end do ! frequencies
      end do ! Frequenct blocks

      ! If inverting
      if (run_mode.eq.-1) then

        ! Keep Stokes
        call setstk(SolF%e_Stk(:,:,iph,ith),Stokes_out, &
                    Input%lim_stk,.False.)

        ! Keep contribution function
        if (Input%out_contr) &
          call setctr(SolF%e_Ctr(:,:,:,iph,ith),ContrG, &
                      Input%lim_ctr)

      ! Synthesis
      else

        ! Write stokes
        call writestk(Input%folder,iph,ith,Frec%omega,Geom, &
                        Stokes_out,Input%lim_stk)
        if (laborted) return

        ! If writing contribution function
        if (Input%out_contr) then

          ! Write contribution function
          call writectr(Input%folder,iph,ith,Frec%omega,Geom, &
                        Atmo%z,ContrG,Input%lim_ctr)
          if (laborted) return

        end if ! Write contribution function
      end if ! Inversion/synthesis

      ! Say completed
      if (gpid.eq.0) then
        umsg = '   Completed'
        call verbose
      end if

      end subroutine emergent_manager

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the RTE for a given LOS\n
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
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
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
      !!       Stokes_s(double(:,:)): Sender Stokes buffer\n
      !!      Contr_s(double(:,:,:)): Sender contribution function
      !!                              buffer\n
      !!           tau1(double(:,:)): Data height tau equal to 1\n
      !!            tau(double(:,:)): Data optical depth\n
      !!         tau1_s(double(:,:)): Height optical depth equal to
      !!                              one sender buffer\n
      !!            etaIM(double(:)): Absorptivity point M\n
      !!         data1M(double(:,:)): RT coeff point M\n
      !!         data1O(double(:,:)): RT coeff point O
      !!   JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the
      !!                              radiation tensors\n
      !!   Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!      JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                              integrated over the absorption
      !!                              profile\n
      !!     JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                              frequency dependence
      subroutine emergent_RT(Atom,LTElines,Atmo,Cont,Frec,Red, &
                             Bfield,Geom,MPID,Input,Flgsg,SolF, &
                             ADD,tau1size,if0,if1,ith,iph,Stokes_s, &
                             Contr_s,tau1,tau,tau1_s,etaIM, &
                             data1M,data1O,JKQC)
      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(in):: Cont
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Bfield_class), intent(in):: Bfield
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(in):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: ADD
      integer, intent(in):: if0,if1,ith,iph,tau1size
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: tau1
      double precision, dimension(:,:), &
                        allocatable, intent(inout), target:: tau
      double precision, dimension(:,:,:), &
                        allocatable, intent(inout):: Contr_s
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: Stokes_s
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: tau1_s
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(in):: JKQC
      double precision, dimension(:), pointer, intent(inout):: etaIM
      double precision, dimension(:,:,:), &
                        pointer, intent(inout):: data1M,data1O

      ! Local

      integer:: ierr,jdir,iz,iz0,iz1,diz,m,o,p,op,mfreq
      integer, dimension(2):: info_b

      double precision:: mu_inv,ct,st,cc,sc,vfac,dsm,dsp,dzm,dzp

      ! Dummy

      double precision, dimension(1,1):: ad2

      ! Pointers

      double precision, dimension(:), pointer:: tauM
      double precision, dimension(:), pointer:: etaIO
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M
      double precision, dimension(:,:), pointer:: p_K2M
      double precision, dimension(:,:), pointer:: p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O
      double precision, dimension(:,:), pointer:: p_K2O
      double precision, dimension(:,:), pointer:: p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP
      double precision, dimension(:,:,:), pointer:: data1P
      complex(kind=8), dimension(:,:,:), pointer:: TSo,TKQo


      ! Nullify pointers
      nullify(p_K0M,p_K1M,p_K2M,p_SM,p_StkM)
      nullify(p_K0O,p_K1O,p_K2O,p_SO,p_StkO)
      nullify(p_K0P,p_SP,etaIO,tauM,data1P)
      nullify(TSo,TKQo)

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
      mfreq = if1-if0+1

      ! Point to TKQ_S
      TSo => Geom%TSL(:,:,:,1)

      ! Trigonometry for Doppler shift
      if (dyn) then
        ct = Geom%L_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = cos(Geom%L_phi(iph))
        sc = sin(Geom%L_phi(iph))
      else
        vfac = 1d0
      end if

      ! If slave
      if (pid.gt.0) then

        ! Wait till last communication was received
        call MPI_WAIT(MPID%request3,MPI_STATUS_IGNORE,ierr)
        call MPI_WAIT(MPID%request4,MPI_STATUS_IGNORE,ierr)
        call MPI_WAIT(MPID%request5,MPI_STATUS_IGNORE,ierr)

      end if

      !
      ! If calculating height of tau=1
      !
      if (Input%out_tau1.or.Input%out_contr) then

        !
        ! First height
        !

        ! Reset cummulative quantities
        tau = 0d0
        tau1(1,:) = 0
        tau1(2,:) = Atmo%z(Rz0)
        tauM => tau(if0:if1,Rz0)

        ! Top boundary
        o = Rz0

        ! Select correct TKQout
        if (Bfield%Bstrength(o).gt.TINYB) then
          TKQo => Geom%TBL(:,:,:,1,o)
        else
          TKQo => TSo
        end if

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
        call RTAbs(Frec,Red,Atom,LTElines,Atmo,vfac,Flgsg, &
                   o,jdir,if0,if1,Cont%ndir, &
                   Cont%c(:,:,:,o),Bfield,TKQo,etaIM)

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

          ! Select correct TKQout
          if (Bfield%Bstrength(o).gt.TINYB) then
            TKQo => Geom%TBL(:,:,:,1,o)
          else
            TKQo => TSo
          end if

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
          call RTAbs(Frec,Red,Atom,LTElines,Atmo,vfac,Flgsg, &
                     o,jdir,if0,if1,Cont%ndir, &
                     Cont%c(:,:,:,o),Bfield,TKQo,etaIO)

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

        ! Nullify pointer
        nullify(tauM)

        ! If outputting the height of tau=1
        if (input%out_tau1) then

          ! Master (serial)
          if (pid.eq.0) then

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

          ! Slave (MPI)
          else

            !
            ! Send to master the tau=1
            !

            ! Wait for last send to finish
            call MPI_WAIT(MPID%request7,MPI_STATUS_IGNORE,ierr)
            call MPI_WAIT(MPID%request8,MPI_STATUS_IGNORE,ierr)

            ! Send indexes
            call MPI_ISEND(pid,1,MPI_INTEGER, &
                           0,2,MPI_COMM_RT, &
                           MPID%request7,ierr)

            ! Send tau data
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

          end if ! Serial/MPI
        end if ! Output tau

      end if ! calculate tau=1

      !
      ! Actual emergence
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
        call top(if0,if1,data1M(:,:,6))

      ! If going up, get bottom boundary
      else

        ! Call bottom boundary
        call bottom(Frec%omega,Atmo%T(iz0),vfac, &
                    if0,if1,data1M(:,:,6))

      endif ! propagation direction

      ! Identify current height
      o = iz0

      ! Select correct TKQout
      if (Bfield%Bstrength(o).gt.TINYB) then
        TKQo => Geom%TBL(:,:,:,1,o)
      else
        TKQo => TSo
      end if

      ! Index for Stokes
      if (PRD.and.ADD) op = o

      ! Calculate radiative coefficients
      call RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                   Geom,vfac,o,jdir,if0,if1, &
                   JKQC(:,:,:,o),Cont%ndir, &
                   Cont%c(:,:,:,o),Bfield,TSo,TKQo, &
                   data1M(:,:,1:5),ad2,.False.)

      ! If calculating contribution function
      if (Input%out_contr) then

        ! Salve, wait till last communication was received
        if (pid.gt.0) &
          call MPI_WAIT(MPID%request5,MPI_STATUS_IGNORE,ierr)
          
        ! Store in buffer. The first does not contribute
        Contr_s(:,:,o) = 0d0

      end if

      ! Identify next height
      p = iz0 + diz

      ! Select correct TKQout
      if (Bfield%Bstrength(p).gt.TINYB) then
        TKQo => Geom%TBL(:,:,:,1,p)
      else
        TKQo => TSo
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
      call RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                   Geom,vfac,p,jdir,if0,if1, &
                   JKQC(:,:,:,p),Cont%ndir, &
                   Cont%c(:,:,:,p),Bfield,TSo,TKQo, &
                   data1O(:,:,1:5),ad2,.False.)

      !
      ! Intermediate heights
      !

      ! For each height this CPU has assigned
      do iz=iz0,iz1,diz

        ! We treat the boundaries outside
        if (iz.eq.iz0.or.iz.eq.iz1) cycle

        ! Allocate P pointers
        allocate(data1P(4,MPID%nf(pid),6))

        ! Identify heights
        m = iz - diz
        o = iz
        p = iz + diz

        ! Calculate distance to previous point
        dsm = (Atmo%z(o) - Atmo%z(m))*mu_inv

        ! Caculate quantities of the next point
        dsp = (Atmo%z(p) - Atmo%z(o))*mu_inv

        ! If tau scale
        if (ztau) then

          ! Get geometrical distance
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

        ! Index for Stokes
        if (PRD.and.ADD) op = p

        ! Select correct TKQout
        if (Bfield%Bstrength(p).gt.TINYB) then
          TKQo => Geom%TBL(:,:,:,1,p)
        else
          TKQo => TSo
        end if

        ! RT coefficients
        call RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                     Geom,vfac,p,jdir,if0,if1, &
                     JKQC(:,:,:,p),Cont%ndir, &
                     Cont%c(:,:,:,p),Bfield,TSo,TKQo, &
                     data1P(:,:,1:5),ad2,.False.)

        ! Point to the data
        p_K0M  => data1M(:,:,1)
        p_K1M  => data1M(:,:,2)
        p_K2M  => data1M(:,:,3)
        p_SM   => data1M(:,:,5)
        p_StkM => data1M(:,:,6)
        p_K0O  => data1O(:,:,1)
        p_K1O  => data1O(:,:,2)
        p_K2O  => data1O(:,:,3)
        p_SO   => data1O(:,:,5)
        p_StkO => data1O(:,:,6)
        p_K0P  => data1P(:,:,1)
        p_SP   => data1P(:,:,5)

        ! Apply short characteristics BESSER
        call RTStep(o,ith,iph,mfreq, &
                    dsm,dsp,p_K0M,p_K1M,p_K2M, &
                    p_SM,p_K0O,p_K1O,p_K2O, &
                    p_SO,p_K0P,p_SP,p_StkM, &
                    p_StkO,.True.)

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
          call RTContr(mfreq,dsm,dsp,dzm,dzp,p_K0M, &
                       p_K0O,p_K1O,p_K2O,p_SO,p_K0P, &
                       p_StkO,tau(if0:if1,o),Contr_s(:,:,o),.True.)

        end if ! calculating contribution function

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
      p_K0M  => data1M(:,:,1)
      p_K1M  => data1M(:,:,2)
      p_K2M  => data1M(:,:,3)
      p_SM   => data1M(:,:,5)
      p_StkM => data1M(:,:,6)
      p_K0O  => data1O(:,:,1)
      p_K1O  => data1O(:,:,2)
      p_K2O  => data1O(:,:,3)
      p_SO   => data1O(:,:,5)
      p_StkO => data1O(:,:,6)

      ! Apply short characteristics LINEAR
      call RTStep(o,ith,iph,mfreq, &
                  dsm,dsp,p_K0M,p_K1M,p_K2M, &
                  p_SM,p_K0O,p_K1O,p_K2O, &
                  p_SO,p_K0P,p_SP,p_StkM, &
                  p_StkO,.False.)

      ! If calculating contribution function
      if (Input%out_contr) then

        ! Calculate vertical distance between points
        if (ztau) then
          dzm = dsm/mu_inv
        else
          dzm = Atmo%z(o) - Atmo%z(m)
        end if

        ! Calculate contribution function
        call RTContr(mfreq,dsm,dsp,dzm,dzp,p_K0M, &
                     p_K0O,p_K1O,p_K2O,p_SO,p_K0P, &
                     p_StkO,tau(if0:if1,o),Contr_s(:,:,o),.False.)

      end if ! calculating contribution function

      ! Serial (Master)
      if (pid.eq.0) then

        ! If calculating contribution function
        if (Input%out_contr) then

          ! Inversion
          if (run_mode.eq.-1) then

            ! Keep contribution function
            call setctr(SolF%e_Ctr(:,:,:,iph,ith),Contr_s, &
                        Input%lim_ctr)

          ! Synthesis
          else

            ! Write contribution function
            call writectr(Input%folder,iph,ith,Frec%omega,Geom, &
                          Atmo%z,Contr_s,Input%lim_ctr)

          end if ! Inversion
        end if ! calculating contribution function

        ! If inverting
        if (run_mode.eq.-1) then

          ! Keep Stokes
          call setstk(SolF%e_Stk(:,:,iph,ith),p_StkO, &
                      Input%lim_stk,.False.)

        ! Synthesis
        else

          ! Write stokes
          call writestk(Input%folder,iph,ith,Frec%omega, &
                        Geom,p_StkO,Input%lim_stk)

        end if ! Inverting

        ! Communicate we finished this direction
        if (gpid.eq.0) then
          umsg = '   Completed'
          call verbose
        end if

      ! MPI (slave)
      else

        !
        ! Send to master
        !

        ! Send indexes
        info_b = (/ pid , 1 /)
        call MPI_ISEND(info_b(1),2,MPI_INTEGER, &
                       0,0, MPI_COMM_RT,MPID%request3, &
                       ierr)

        ! If calculating contribution function
        if (Input%out_contr) then

          ! Send contribution function
          call MPI_ISEND(contr_s(0,1,Rz0), &
                         MPID%size3(pid), MPI_DOUBLE_PRECISION, &
                         0, 1+pid, MPI_COMM_RT, &
                         MPID%request5, ierr)
        end if

        ! Send intensity
        Stokes_s = data1O(:,:,6)
        call MPI_ISEND(Stokes_s(0,MPID%if0(pid)), &
                       MPID%size10(pid), MPI_DOUBLE_PRECISION, &
                       0, pid, MPI_COMM_RT, &
                       MPID%request4, ierr)

      end if

      ! Clean pointers
2000  nullify(p_K0M,p_K1M,p_K2M,p_SM,p_StkM)
      nullify(p_K0O,p_K1O,p_K2O,p_SO,p_StkO)
      nullify(p_K0P,p_SP,TSo,TKQo)

      end subroutine emergent_RT

!#####################################################################
!#####################################################################
!#####################################################################

      end module solver_mod
