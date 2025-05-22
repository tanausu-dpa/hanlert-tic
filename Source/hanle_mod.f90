      !> Flow control of the code
      module hanle_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Contributors:
!     Hao Li (IAC/NSSCC)
!  Start:
!     18/04/2017
!  Last version:
!     15/05/2025 V4.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/05/2025:    V4.0.5 - Generalized declarations of Atom, Atomb,
!                             Mol, and Rho_old to allow for empty
!                             arrays for any of them (TdPA)
!                           - Consider that the integrated radiation
!                             field tensors may not be needed (TdPA)
!                           - Add the existence of active atoms to the
!                             requirements to perform iterations at
!                             all (TdPA)
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
!  hanle
!    Solve the NLTE self-consistent problem and calculate the emergent
!  profiles
!
!  hanle_setup
!    Prepare the background quantities and possible restrictions in
!  the range of heights to consider in the self-consistent NLTE
!  problem. Initialize thermal part of photoionization rates,
!  photoionization auxiliar quantities, and density matrices
!
!  hanle_reback
!    Recalculate the background continuum quantities
!
!  hanle_init
!    Initialize radiation and density matrices
!
!  hanle_intensity
!    Solve the self-consistent NLTE problem for intensity or/and
!  calculate the emergent intensity profiles
!
!  hanle_polarization
!    Solve the self-consistent NLTE problem or the second kind or/and
!  calculate the emergent Stokes profiles
!
!  prepare_buffers
!    Allocate the arrays to store the solution of the forward problem
!  in the inversion mode and configure the synthesis mode for the
!  inversion code, depending on if the synthesis is a trial or a
!  response function
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use background_mod
      use commons_mod
      use diagon_mod
      use free_mod
      use gauss_mod
      use getztau_mod
      use initialize_mod
      use initphotoion_mod
      use initpopu_mod
      use iosolution_mod
      use normalizer_mod
      use omegabuild_mod
      use parameters_mod , only : TINYB
      use setmpi_mod
      use solver_mod
      use solveri_mod
      use strength_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the NLTE self-consistent problem and calculate the
      !! emergent profiles\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!        Atomb(Atom_class(:)): Structures with atomic data for
      !!                              background atoms\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!           Mol(Mol_class(:)): Structures with molecular data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!       GeomI(Geometry_class): Structure with geometric data
      !!                              for the intensity problem\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!          fudge(fudge_class): Structure with fudge data\n
      !!        kurucz(kurucz_class): Structure with Kurucz line
      !!                              data\n
      !!            JKQin(double(:)): Data with ad-hoc JKQ tensors\n
      !!      SolF(Solution_F_class): Structure with the solution of
      !!                              the self-consistent problem and
      !!                              the corresponding emergent
      !!                              profiles, contribution function,
      !!                              and height for optical depth
      !!                              equal to one\n
      !!              lload(logical): If reading an existing solution
      !!                              file\n
      !!                lio(logical): If solving the self-consistent
      !!                              NLTE problem for intensity\n
      !!                lie(logical): If calculating the emergent
      !!                              intensity\n
      !!                 lp(logical): If solving the self-consistent
      !!                              NLTE problem of the 2nd kind\n
      !!                lpe(logical): If calculating the emergent
      !!                              Stokes parameters\n
      !!               free(logical): If allowed to free some input
      !!                              data
      subroutine hanle(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input, &
                       GeomI,Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                       JKQin,SolF,lload,lio,lie,lp,lpe,free)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTElines
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(inout):: Frec
      type(Geometry_class), intent(inout):: GeomI,Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: lload,lio,lie,lp,lpe,free
      double precision, dimension(:), allocatable, intent(in):: JKQin

      ! Local

      type(Bfield_class):: Bfield0
      type(Continuum_class):: Cont
      type(Rhoc_class), allocatable, dimension(:):: Rho_old

      logical:: rlimw = .True.
      logical:: csize
      logical:: rback,rdyn,raxial
      logical:: l1, l2, ofram
      logical:: rVIRAM,rVPRAM,rWIRAM,rWPRAM

      integer:: iph,rnPh

      double precision:: rVPhi,rVmux,rVmuy,rWmux
      double precision, dimension(:,:,:,:), allocatable:: StokesI
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes
      double precision, dimension(:,:), allocatable:: J00
      double precision, dimension(:,:), allocatable:: J00S
      double precision, dimension(:,:), allocatable:: J00C
      double precision, dimension(:,:,:), allocatable:: J00P

      complex(kind=8), dimension(:,:,:), allocatable:: JKQ_asym
      complex(kind=8), dimension(:,:,:), allocatable:: JKQ_asym_fake
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC
#ifdef DEBUGSYN
      ! Shift ID so everything is verbosed
      ! Not necessary if 1D because the global master is here
      if (run_mode.ne.0) gpid = gpid - 1
#endif

      ! Memory count for locally allocated types
      MRAMc = MRAMc + 1d-6*sizeof(Bfield0)
      MRAMc = MRAMc + 1d-6*sizeof(Cont)
      if (nA.gt.0) then
        allocate(Rho_old(nA))
        MRAMc = MRAMc + 1d-6*sizeof(Rho_old)
      end if
      if (allocated(JKQ_asym)) RRAMc = 1d-6*(sizeof(JKQ_asym))

      ! Original RAM storage flags
      rVIRAM = VIRAM
      rVPRAM = VPRAM
      rWIRAM = IRAM
      rWPRAM = PRAM

      ! Original dynamic flag
      rdyn = dyn

      ! If forcing the problem to be static for intensity, but it
      ! is dynamic
      if (dyn.and.Input%static_int) then

        ! Trick the problem pointing the velocity to the array
        ! of zeros
        dyn = .False.
        Atmo%vxa => Atmo%vx
        Atmo%vya => Atmo%vy
        Atmo%vza => Atmo%vz
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros
        Atmo%vz => Atmo%zeros

      ! If axial intensity and not polarization
      else if (axiali.and..not.axial) then

        ! Cheat the velocities by pointing to the array of zeros
        Atmo%vxa => Atmo%vx
        Atmo%vya => Atmo%vy
        Atmo%vx => Atmo%zeros
        Atmo%vy => Atmo%zeros
        nullify(Atmo%vza)

      ! Otherwise, the alternative velocity pointers are null
      else

        ! Initialize
        nullify(Atmo%vxa)
        nullify(Atmo%vya)
        nullify(Atmo%vza)

      end if

      ! Call initialization
      call hanle_setup(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input, &
                       GeomI,Geom,Bfield,Frec,Flgsg,fudge, &
                       kurucz,Cont,Rho_old,free,rback)

      ! Control
      if (laborted) goto 1000


      !
      ! Check if we need to reevaluate sizes
      !

      ! If the vertical dimension has been restricted
      if (Rnz.ne.nz) then

        ! Shorten variables for code flow
        l1 = .not.lload
        l2 = (lio.and..not.lie).or.(lload.and.(lp.or.lpe))

        ! If doing MPI, recalculate the message sizes
        if (MPID%mpi) &
          call setmpi_sizes(MPID,GeomI,Geom,Frec,lio,lp,l1,l2, &
                            Input%ALI_photo,.True.)

      end if ! Restricting the vertical dimension

      ! Control
      if (laborted) goto 1000

      !
      ! Initialize solution
      !
      call hanle_init(Atom,Atmo,Input,GeomI,Geom,Bfield, &
                      Frec,Flgsg,SolF,lload,lio,lie,lp,lpe, &
                      Stokes,JKQ,JKQS,JKQC, &
                      StokesI,J00,J00S,J00C,J00P)

      ! Control
      if (laborted) goto 1000

#ifdef DEBUGATMO
      ! Dump current model atmosphere in ASCII file
      if (pid.eq.0) call dump_atmo(Atmo,Bfield,Input%folder,1)
#endif

      !
      ! Self-consistent NLTE problem for intensity
      !

      !
      ! Solve the NLTE self-consistent problem for intensity
      ! or/and get emergent intensity
      if (lio.or.lie) &
        call hanle_intensity(Atom,LTElines,Atmo,MPID,Input, &
                             GeomI,Bfield,Frec,Flgsg,SolF, &
                             Cont,Rho_old, &
                             StokesI,J00,J00S,J00C,J00P, &
                             lload,lio,lie,lp.or.lpe,rlimw,ofram)

      ! Control
      if (laborted) goto 1000

      !
      ! Finished self-consistent NLTE problem for intensity
      !

      ! If updating the model
      if (Input%update_atmos.ge.0) then

        ! If doing polarization
        if (pid.eq.0.and.(lp.or.lpe)) then

          ! Message
          umsg = 'The option update_atmos is not compatible '// &
                 'with computing the polarization. The code '// &
                 'is stopping after updating the atmosphere'
          call verbose

        end if ! Master and doing polarization

        ! Leave hanle already
        goto 1000

      ! Not updating the model
      else

        ! If chose to recalculate electron number density
        if (Input%redo_ne.eq.1.or.Input%redo_ne.eq.11) then

          ! Master writes
          if (pid.eq.0) then

            ! Message
            umsg = ' # Warning: You chose to redo electrons '// &
                   'at the end of iterations, but not '// &
                   'to update the atmosphere. Electrons are '// &
                   'not recomputed.'
            call verbose

          end if ! Master
        end if ! Recalculating electron number density
      end if ! Updating or not the model

      ! If not in inversion mode
      if (run_mode.ne.-1) then

        !
        ! Deallocate atmospheric quantities if still in RAM

        ! Gas pressure
        if (allocated(Atmo%Pg)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pg)
          deallocate(Atmo%Pg)
        end if

        ! Electron pressure
        if (allocated(Atmo%Pe)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%Pe)
          deallocate(Atmo%Pe)
        end if

        ! Mass density
        if (allocated(Atmo%rho)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%rho)
          deallocate(Atmo%rho)
        end if
      end if ! Not inversion

      ! If permission to free model memory and no need to redo
      ! background opacity later
      if (free.and..not.rback) then

        ! If there are background atoms
        if (allocated(Atomb)) then

          ! Free from RAM
          nAb = 0
          call free_atom_full(Atomb)

        end if ! There are background atoms
      end if ! Permission to free and not redoing continuum opacity

      ! If forcing the problem to be static for intensity, but it
      ! is dynamic
      if (rdyn.and.Input%static_int) then

        ! Return the problem to its original state
        dyn = rdyn

        ! Point the velocities back to original data
        Atmo%vx => Atmo%vxa
        Atmo%vy => Atmo%vya
        Atmo%vz => Atmo%vza

        ! Nullify auxiliar pointers
        nullify(Atmo%vxa,Atmo%vya,Atmo%vza)

      ! If axial intensity and not polarization, un-cheat
      ! the velocities
      else if (axiali.and..not.axial) then

        ! Point the velocities back to original data
        Atmo%vx => Atmo%vxa
        Atmo%vy => Atmo%vya

        ! Nullify auxiliar pointers
        nullify(Atmo%vxa,Atmo%vya)

      end if ! Different velocity between intensity and polarizaiton

      !
      ! Self-consistent NLTE problem of the second-kind
      !

      ! If doing polarization calculations
      if (lp.or.lpe) then

        ! Initialize extra asymmetry radiation field tensors
        call initialize_asym(Input,Flgsg,JKQin,JKQ_asym)

        !
        ! Update radiation RAM if coming from intensity
        !

        ! If did intensity iterations and not emergence or if
        ! we loaded a file which did not contain polarization
        if ((lio.and..not.lie).or. &
            (lload.and..not.allocated(Stokes))) then

          ! This is what is going to be kept in radiation
          ! memory
          RRAMc = 1d-6*sizeof(JKQ_asym)

          !
          ! Pre-compute amount of RAM to fill with radiation later
          !

          ! If need Stokes
          if (KSTK) then

            ! Full size
            RRAMc = RRAMc + 8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*Rnz)

          ! No need of Stokes
          else

            ! Only two heights
            RRAMc = RRAMc + 8d-6*dble(4*nfreq*Geom%nPh*Geom%nTh*2)

          end if ! Keeping Stokes

          ! Size in radiation field tensors
          RRAMc = RRAMc + 16d-6*dble(Rnz*15*nfreq)

          ! If Atoms
          if (nA.gt.0) then

            ! Atomic radiation field tensors
            RRAMc = RRAMc + 1d-6*sizeof(J00P) + &
                            16d-6*dble(Rnz*15*2*nxtran)

          end if ! Atoms
        end if ! Recalculate RAM for radiation

        ! If we need to redo the background, do it
        if (rback) &
          call hanle_reback(Atom,Atomb,Mol,Atmo,MPID,Input,Geom, &
                            Frec,Flgsg,fudge,kurucz,Cont,free)

        !
        ! If there is a magnetic field, and doing a two-step
        ! calculation
        if (maxval(Bfield%Bstrength).gt.TINYB.and. &
            Input%two_step_pol) then

          ! Set zero field structure
          allocate(Bfield0%Bstrength(nz))
          allocate(Bfield0%Btheta(nz))
          allocate(Bfield0%Bphi(nz))
          Bfield0%Bstrength = 0d0
          Bfield0%Btheta = 0d0
          Bfield0%Bphi = 0d0

          ! Memory count
          MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bstrength)
          MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Btheta)
          MRAMc = MRAMc + 1d-6*sizeof(Bfield0%Bphi)

          ! Get current status of geometry
          raxial = axial
          rnPh = Geom%nPh
          rVPhi = Geom%V_phi(1)
          rVmux = Geom%V_mux(1)
          rVmuy = Geom%V_muy(1)
          rWmux = Geom%W_mux(1)
          csize = .False.

          ! Master verbose
          if(gpid.eq.0) then
            umsg = ' - Solving the polarized but '// &
                   'non-magnetic problem'
            call verbose
          end if

          ! If intensity was static or axial
          if (Input%static_int.or.axiali) then

            ! Check axial velocity
            if (maxval(abs(Atmo%vx)).le.0d0.and. &
                maxval(abs(Atmo%vy)).le.0d0) then

              ! And fake it
              axial = .True.
              Geom%nPh = 1
              Geom%V_phi(1) = 0d0
              Geom%V_mux(1) = 1d0
              Geom%V_muy(1) = 1d0
              Geom%W_mux(1) = 1d0

              ! Fake sizes for messages
              csize = .True.
              MPID%size4 = MPID%size4/rnPh
              MPID%size5 = MPID%size5/rnPh

            end if ! No horizontal velocity
          end if ! Intensity was axial or static

          ! If the polarization problem can be axial
          if (axial) then

            ! Call polarization solution WITHOUT field or
            ! emergence with fake and dummy JKQ
            call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                    Geom,Bfield0,Frec,Flgsg,SolF, &
                                    Cont,Rho_old, &
                                    StokesI,J00,J00S,J00C,J00P, &
                                    Stokes,JKQ,JKQS,JKQC, &
                                    JKQ_asym_fake,rnPh,.False., &
                                    lload,lio,lp,.False., &
                                    rlimw,ofram)

            ! If will not be axial
            if (.not.raxial) then

              ! Copy Stokes for the rest of azimuth
              do iph=2,rnPh
                Stokes(:,:,iph,:,:) = Stokes(:,:,1,:,:)
              end do

            end if ! Will not be axial

          ! If the problem cannot be solve as axial anyways
          else

            ! Call polarization solution WITHOUT field
            call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                    Geom,Bfield0,Frec,Flgsg,SolF, &
                                    Cont,Rho_old, &
                                    StokesI,J00,J00S,J00C,J00P, &
                                    Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                    rnPh,.False., &
                                    lload,lio,lp,.False., &
                                    rlimw,ofram)
          end if

          ! Restore the geometry variables
          axial = raxial
          Geom%nPh = rnPh
          Geom%V_phi(1) = rVPhi
          Geom%V_mux(1) = rVmux
          Geom%V_muy(1) = rVmuy
          Geom%W_mux(1) = rWmux

          ! If changed the sizes, retore them
          if (csize) then
            MPID%size4 = MPID%size4*rnPh
            MPID%size5 = MPID%size5*rnPh
          end if

          ! Solve the actual polarization problem
          call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                  Geom,Bfield,Frec,Flgsg,SolF, &
                                  Cont,Rho_old, &
                                  StokesI,J00,J00S,J00C,J00P, &
                                  Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                  Geom%nPh,.True., &
                                  lload,.False.,lp,lpe, &
                                  rlimw,ofram)

        ! Normal run
        else

          ! Solve the polarization problem
          call hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                  Geom,Bfield,Frec,Flgsg,SolF, &
                                  Cont,Rho_old, &
                                  StokesI,J00,J00S,J00C,J00P, &
                                  Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                  Geom%nPh,.True., &
                                  lload,lio,lp,lpe,rlimw,ofram)

        end if ! Zero field first step solution
      end if ! Polarization

      !
      ! Finished self-consistent NLTE problem of the second-kind
      !


      !
      ! Clean memory
      !

      ! Clean local variables
1000  call free_hanle(Atom,Cont,Geom,Frec,Bfield0,Rho_old,JKQ_asym)
#ifdef DEBUGSYN
      ! Recover ID before leaving
      if (run_mode.ne.0) gpid = gpid + 1
#endif

      ! Restore RAM flags
      VIRAM = rVIRAM
      VPRAM = rVPRAM
      IRAM = rWIRAM
      PRAM = rWPRAM

      ! Restore dyn just in case we got here due to error
      dyn = rdyn

      ! Memory count
      if (allocated(JKQ_asym)) then
        RRAMc = 1d-6*(sizeof(JKQ_asym))
      else
        RRAMc = 0d0
      end if

      return

      end subroutine hanle

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the background quantities and possible restrictions
      !! in the range of heights to consider in the self-consistent
      !! NLTE problem. Initialize thermal part of photoionization
      !! rates, photoionization auxiliar quantities, and density
      !! matrices\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!        Atomb(Atom_class(:)): Structures with atomic data for
      !!                              background atoms\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!           Mol(Mol_class(:)): Structures with molecular data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!       GeomI(Geometry_class): Structure with geometric data
      !!                              for the intensity problem\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!          fudge(fudge_class): Structure with fudge data\n
      !!        kurucz(kurucz_class): Structure with Kurucz line
      !!                              data\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!      Rho_old(Rhoc_class(:)): Structure to store the density
      !!                              matrix of the previous
      !!                              iteration\n
      !!               free(logical): If allowed to free some input
      !!                              data
      !!              rback(logical): Indicate if the background needs
      !!                              to be recalculated later
      subroutine hanle_setup(Atom,Atomb,LTElines,Mol,Atmo,MPID, &
                             Input,GeomI,Geom,Bfield,Frec,Flgsg, &
                             fudge,kurucz,Cont,Rho_old, &
                             free,rback)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(inout):: LTElines
      type(Mol_class), dimension(:), allocatable, intent(inout):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Continuum_class), intent(inout):: Cont
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(inout):: Frec
      type(Geometry_class), intent(in):: GeomI,Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Rhoc_class), dimension(:), &
                        allocatable, intent(inout):: Rho_old
      logical, intent(in):: free
      logical, intent(out):: rback

      ! Local

      integer:: ia,ios


      ! Routine name
      urou = 'hanle_setup'

      ! Control
      if (laborted) return

      ! Allocate space for chi scale
      if (.not.allocated(Atmo%chi500)) then
        allocate(Atmo%chi500(nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%chi500)
      end if

      ! Calculate background continuum quantities
      call background(Atom,Atomb,Mol,Atmo,fudge,kurucz, &
                      Input,Frec%omega,Cont,GeomI,MPID,Flgsg)

      ! Memory count
      if (allocated(Cont%c)) then
        BRAMc = 1d-6*sizeof(Cont%c)
      else
        BRAMc = 0d0
      end if

      ! Control
      if (laborted) return

      ! Global master verbosity
      if (gpid.eq.0) then
        umsg = ' - Background quantities calculated'
        call verbose
      end if


      !
      ! Do we need to repeat the background calculation?
      !

      ! Initialize
      rback = .False.

      ! If bound-bound transitions in background
      if (Input%addbb) then

        ! If slave or not mpi
        if (.not.MPID%mpi.or.pid.gt.0) then

          ! If size of directions larger than 1
          if (size(Cont%c,3).gt.1) then

            ! If geometry changes between quadratures
            ! then we repeat
            rback = Geom%nPh.ne.GeomI%nPh.or. &
                    Geom%nTh.ne.GeomI%nTh

          end if ! More than one direction in continuum (b-b present)
        end if ! Slave or not MPI

        ! If MPI, check if anyone found the need to repeat the
        ! continuum calculation
        if (MPID%mpi) &
          call MPI_ALLREDUCE(MPI_IN_PLACE,rback,1,MPI_LOGICAL, &
                             MPI_LOR,MPI_COMM_RT,ios)

      end if ! B-B can be in background

      ! Calculate continuum opacity at reference frequency
      call chi_freq(Atom,Atomb,Mol,Atmo,fudge,Input, &
                    Atmo%tfreq,Atmo%chi500,1,nz,MPID%mpi)

      ! If master doing MPI
      if (pid.eq.0.and.MPID%mpi) then

        ! Remove space for chi scale
        MRAMc = MRAMc - 1d-6*sizeof(Atmo%chi500)
        deallocate(Atmo%chi500)

      end if ! Master and MPI

      ! Control
      if (laborted) return

      !
      ! Compute missing height or tau, whatever is not the
      ! input axis
      call getztau(Atmo,.True.)


      !
      ! Store atmosphere in ASCII
      if (Input%keep_atmo) call writeatmo(Atmo,Bfield, &
                                          Input%folder, &
                                          Input%lim_atmo)
      ! Control
      if (laborted) return

      ! If chi500 is allocated
      if (allocated(Atmo%chi500)) then

        ! If the master and doing MPI and no use for it later
        ! to update or write the model atmosphere
        if (MPID%mpi.and.pid.eq.0.and. &
            (Input%update_atmos.lt.0.and..not.Input%keep_atmo)) then

          ! Free RAM
          MRAMc = MRAMc - 1d-6*sizeof(Atmo%chi500)
          deallocate(Atmo%chi500)

        ! Otherwise
        else

          ! If not using the tau scale and not updating the atmosphere
          if (.not.ztau.and.Input%update_atmos.lt.0.and. &
              .not.Input%keep_atmo) then

            ! Free RAM
            MRAMc = MRAMc - 1d-6*sizeof(Atmo%chi500)
            deallocate(Atmo%chi500)

          end if ! Using the tau scale nor updating the atmosphere
        end if ! If Master and MPI
      end if ! Allocated chi500 variable

      ! If not recalculating electron number density or updating
      ! the model atmosphere
      if (Input%redo_ne.le.0.and. &
          (Input%update_atmos.lt.0.and.Atmo%typo.eq.0)) then

        ! If allocated gas pressure
        if (allocated(Atmo%Pg)) then

          ! If not an inversion (which needs this)
          if (run_mode.ne.-1) then

            ! Clean not to be used atmospheric variables
            MRAMc = MRAMc - 1d-6*(sizeof(Atmo%Pg))
            MRAMc = MRAMc - 1d-6*(sizeof(Atmo%Pe))
            MRAMc = MRAMc - 1d-6*(sizeof(Atmo%rho))
            deallocate(Atmo%Pg)
            deallocate(Atmo%Pe)
            deallocate(Atmo%rho)

          end if ! Not inversion
        end if ! Gas pressure is allocated

        ! Can free 1D quantities and no need to re-compute continuum
        ! opacity later
        if (free.and..not.rback) then

          ! Free background atoms
          nAb = 0
          call free_atom_full(Atomb)

        end if ! Can free 1D quantities and not recalculating back
      end if ! Not recalculating electron number density

      ! Deallocate background atoms and molecules
      if (free.and..not.rback) then
        nM = 0
        call free_mol_full(Mol)
      end if

      ! Control
      call control

      ! Control
      if (laborted) return

      ! If keeping background, call writer
      if (Input%keep_back) call writeback(Cont,Frec%omega, &
                                          Input%folder,MPID, &
                                          Input%lim_back)
      ! Control
      if (laborted) return


      !
      ! Restrict heights
      call restrict_zaxis(Atmo,Input)


      !
      ! If there are LTE lines, check if they need to be restricted
      ! in height
      if (Input%nLTE.gt.0) &
        call restrict_LTE_lines(Atmo,LTElines)


      ! For each active atom
      do ia=1,nA

        ! Compute thermal part of the photoionization rates in SEE
        call setphotoTEI(Atom(ia),Frec,Atmo%T,Atmo%ne,free)

        ! And initialize  the density matrix
        call Initcrho_old(Atom(ia),Rho_old(ia))

      end do ! Active atoms

      ! Control
      if (laborted) return

      ! Master verbosity
      if (gpid.eq.0) then
        umsg = ' - Initialized photoionization quantities (thermal)'
        call verbose
      end if

      !
      ! Store quantities for epsIphoto
      !

      ! If storing photoionization auxiliar variables in RAM
      if (PIRAM) then

        ! Call the initializer of these variables
        call ramphoto(Atom,Frec,Atmo%T)

      ! If not storing photoionization auxiliar variables in RAM
      else

        ! If in serial mode
        if (.not.MPID%mpi) then

          ! Allocate a dummy index to avoid memory 'leaks'
          allocate(Frec%exu(1,1))
          PRAMc = PRAMc + 1d-6*sizeof(Frec%exu)

        ! If in MPI
        else

          ! Just point to null
          nullify(Frec%exu)

        end if ! MPI
      end if ! Storing or not photoionization auxiliar quantities

      end subroutine hanle_setup

!#####################################################################
!#####################################################################
!#####################################################################

      !> Recalculate the background continuum quantities\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!   Atomb(Atom_class(:)): Structures with atomic data for
      !!                         background atoms\n
      !!      Mol(Mol_class(:)): Structures with molecular data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!     Input(Input_class): Structure with configuration data\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                         J-symbols\n
      !!     fudge(fudge_class): Structure with fudge data\n
      !!   kurucz(kurucz_class): Structure with Kurucz line data\n
      !!  Cont(Continuum_class): Structure with background opacity
      !!                         data\n
      !!          free(logical): If allowed to free some input data
      subroutine hanle_reback(Atom,Atomb,Mol,Atmo,MPID,Input,Geom, &
                              Frec,Flgsg,fudge,kurucz,Cont,free)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), allocatable, intent(inout):: Mol
      type(Atmo_class), intent(in):: Atmo
      type(Continuum_class), intent(inout):: Cont
      type(Fctsg_class), intent(in):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: free

      ! Clean current continuum
      if (allocated(Cont%c)) deallocate(Cont%c)

      ! Do it again
      call background(Atom,Atomb,Mol,Atmo,fudge,kurucz, &
                      Input,Frec%omega,Cont,Geom,MPID,Flgsg)

      ! Memory count
      if (allocated(Cont%c)) then
        BRAMc = 1d-6*sizeof(Cont%c)
      else
        BRAMc = 0d0
      end if

      ! Control
      if (laborted) return

      ! Master verbose
      if (gpid.eq.0) then
        umsg = ' - Background quantities re-calculated'
        call verbose
      end if

      ! If we can free the background atoms and molecules here
      if (free) then

        ! Free background atoms
        nAb = 0
        call free_atom_full(Atomb)

        ! Free molecules
        nM = 0
        call free_mol_full(Mol)

      end if ! We can free background atoms and molecules here

      end subroutine hanle_reback

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize radiation and density matrices\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!           Atmo(Atmo_class): Structure with atmospheric
      !!                             data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with geometric data
      !!                             for the intensity problem\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         Flgsg(Fctsg_class): Structure with factorials,
      !!                             signs, and J-symbols\n
      !!     SolF(Solution_F_class): Structure with the solution of
      !!                             the self-consistent problem and
      !!                             the corresponding emergent
      !!                             profiles, contribution function,
      !!                             and height for optical depth
      !!                             equal to one\n
      !!             lload(logical): If reading an existing solution
      !!                             file\n
      !!               lio(logical): If solving the self-consistent
      !!                             NLTE problem for intensity\n
      !!               lie(logical): If calculating the emergent
      !!                             intensity\n
      !!                lp(logical): If solving the self-consistent
      !!                             NLTE problem of the 2nd kind\n
      !!               lpe(logical): If calculating the emergent
      !!                             Stokes parameters\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence\n
      !!   StokesI(double(:,:,:,:)): Intensity\n
      !!           J00(double(:,:)): Mean intensity integrated over
      !!                             the absorption profile\n
      !!          J00S(double(:,:)): Mean intensity integrated over
      !!                             the emission profile\n
      !!          J00C(double(:,:)): Mean intensity with frequency
      !!                             dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates
      subroutine hanle_init(Atom,Atmo,Input,GeomI,Geom,Bfield, &
                            Frec,Flgsg,SolF,lload,lio,lie,lp,lpe, &
                            Stokes,JKQ,JKQS,JKQC, &
                            StokesI,J00,J00S,J00C,J00P)
      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(in):: Flgsg
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: GeomI,Geom
      type(Input_class), intent(in):: Input
      type(Solution_F_class), intent(in):: SolF
      logical, intent(in):: lload,lio,lie,lp,lpe
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(out):: StokesI
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Stokes
      double precision, dimension(:,:), allocatable, intent(out):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: J00P
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQC

      ! Local

      logical:: read_stokes

      integer:: ia


      !
      ! If to load a previous solution
      if (lload) then

        ! Doing inversions
        if (run_mode.eq.-1) then

          ! Get solution from RAM
          call getsol(SolF,GeomI,Geom,Flgsg,Bfield,Atom, &
                      Stokes,JKQ,JKQS,JKQC, &
                      StokesI,J00,J00S,J00C,J00P,lio.or.lie)

          ! Master notify
          if (gpid.eq.0) then
            umsg = ' - Copied from solution'
            call verbose
          end if

        ! Doing synthesis
        else

          ! Master notify
          if (gpid.eq.0) then
            umsg = ' - Reading solution'
            call verbose
          end if

          ! Get solution from file
          call readsol(Input%solution,GeomI,Geom,Flgsg,Bfield, &
                       Atom,read_stokes, &
                       Stokes,JKQ,JKQS,JKQC, &
                       StokesI,J00,J00S,J00C,J00P)

          ! Control
          if (laborted) return

          ! If loaded an intensity solution but not iterating
          ! polarization before computing emergence
          if (lload.and.lpe.and..not.lp.and. &
              .not.lio.and..not.allocated(JKQ)) then

            ! Warning of unsuitability of the solution file for this
            ! specific run
            umsg = ' # Warning: Loading an intensity solution '// &
                   'to compute polarized emergent profiles '// &
                   'without iterating'
            call verbose

          end if ! Unsuitable Solution file for this run

          ! If loaded a solution for a polarization problem but we
          ! need to iterate intensity
          if (lio.and..not.allocated(J00)) then

            ! This is a critical error
            umsg = 'Cannot read multiterm solution as '// &
                   'initialization for multilevel'
            call aborted
            return

          end if ! Wrong solution file for this run

          ! If could not read Stokes and frequency dependent radiation
          ! field tensors
          if (.not.read_stokes) then

            ! If read polarization
            if (allocated(JKQC)) then

              ! Initialize radiation field
              call initialize_failread(Frec,Atmo,Stokes,JKQC)

            ! If read intensity
            else

              ! Initialize radiation field
              call initializeI_failread(Frec,Atmo,StokesI,J00C)

            end if ! Type of solution
          end if ! Could read Stokes
        end if ! Inversion/synthesis


      !
      ! Not reading a previous solution
      !
      else

        ! If going through the intensity problem
        if (lio.or.lie) then

          ! Initialize unpolarized radiation field
          call initializeI(Frec,GeomI,Atmo, &
                           StokesI,J00,J00S,J00C,J00P,Input%mode)

        ! If not going through the intensity problem, but through
        ! the polarized one
        else if (lp.or.lpe) then

          ! Initialize polarized radiation field
          call initialize(Frec,Geom,Atmo, &
                          Stokes,JKQ,JKQS,JKQC,J00P,Input%mode)

        end if ! Intensity/polarization problem

        ! Control
        if (laborted) return

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Radiation Field Initialized'
          call verbose
        end if

        ! For every atom
        do ia=1,nA

          ! Normalize the populations
          call correctpop(Atom(ia),0)

        end do ! Active atoms

        ! Control
        if (laborted) return

      end if

      end subroutine hanle_init

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the self-consistent NLTE problem for intensity or/and
      !! calculate the emergent intensity profiles\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!       GeomI(Geometry_class): Structure with geometric data
      !!                              for the intensity problem\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!      SolF(Solution_F_class): Structure with the solution of
      !!                              the self-consistent problem and
      !!                              the corresponding emergent
      !!                              profiles, contribution function,
      !!                              and height for optical depth
      !!                              equal to one\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!      Rho_old(Rhoc_class(:)): Structure to store the density
      !!                              matrix of the previous
      !!                              iteration\n
      !!    StokesI(double(:,:,:,:)): Intensity\n
      !!            J00(double(:,:)): Mean intensity integrated over
      !!                              the absorption profile\n
      !!           J00S(double(:,:)): Mean intensity integrated over
      !!                              the emission profile\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence\n
      !!         J00P(double(:,:,:)): Intensity integrals in the
      !!                              photoionization rates
      !!              lload(logical): If reading an existing solution
      !!                              file\n
      !!                lio(logical): If solving the self-consistent
      !!                              NLTE problem for intensity\n
      !!                lie(logical): If calculating the emergent
      !!                              intensity\n
      !!                 lp(logical): If solving the self-consistent
      !!                              NLTE problem of the 2nd kind\n
      !!              rlimw(logical): If we can write a message about
      !!                              reaching the RAM limit\n
      !!              ofram(logical): If we have reached the RAM limit
      subroutine hanle_intensity(Atom,LTElines,Atmo,MPID,Input, &
                                 GeomI,Bfield,Frec,Flgsg,SolF, &
                                 Cont,Rho_old, &
                                 StokesI,J00,J00S,J00C,J00P, &
                                 lload,lio,lie,lp,rlimw,ofram)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(inout):: Flgsg
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(inout):: GeomI
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      type(Solution_F_class), intent(inout):: SolF
      type(Continuum_class), intent(in):: Cont
      type(Rhoc_class), dimension(:), &
                        allocatable, intent(inout):: Rho_old
      logical, intent(in):: lload,lio,lie,lp
      logical, intent(inout):: rlimw,ofram
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: StokesI
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(inout):: J00P

      ! Local

      type(Red_class):: Red

      logical:: liter,literJ

      double precision, dimension(1):: ad1
      double precision, dimension(:,:,:,:,:), allocatable:: Stokes

      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQ
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQS
      complex(kind=8), dimension(:,:,:,:), allocatable:: JKQC


      ! Initialize redistribution structure pointers
      nullify(Red%dzao,Red%rzao,Red%pzao,Red%zao)

      ! Initialize memory used in input frequency data
      FRAMc = FRAMc + 1d-6*sizeof(Red)

      ! Check if doing any iteration in J or in populations
      literJ = Input%iter_J.gt.0
      liter = lio.and.(Input%iteri_max.ge.Input%iteri_min).and.nA.gt.0

      !
      ! Preliminars
      !

      ! Doing absolutely anything
      if (liter.or.literJ.or.lie) then

        !
        ! Set number of output directions
        !

        ! If iterating or using lines in J iteration
        if (liter.or.(literJ.and.Input%init_J_bb)) then

          ! Number of output directions is quadrature
          GeomI%njdir = GeomI%nTh*GeomI%nPh

          ! Index geometry
          call geom_index(GeomI,.False.)

        ! Not using quadrature for lines
        else

          ! Number of output directions given by LOS
          GeomI%njdir = 1

          ! Index geometry
          call geom_index(GeomI,.True.)

        end if

        ! If doing PRD and doing any intensity formal solution
        if (Input%iter_ord.eq.2.and.(lie.or.liter)) then

          ! Index redistribution data arrays, allocate for
          ! input frequencies, and estimate minimum memory
          ! for PRD normalization
          call index_red(Atom,Red,ad1,.False.)

          ! Call input frequency axis
          call omegabuildinI(Frec,Red,Atom,Atmo,Input,MPID,ofram)
          if (laborted) goto 1000

          ! Master verbose
          if (gpid.eq.0.and.nA.gt.0) then
            umsg = ' - Input frequency axis initialized (intensity)'
            call verbose
          end if ! Master
        end if ! PRD and any intensity formal solution

        ! If solving NLTE proble, predict RAM necessary for RT
        if (liter) &
          call solveI_predict(Atom,Atmo,Frec,Red,GeomI,MPID,Input)

        ! If we need actual line data
        if (liter.or.(literJ.and.Input%init_J_bb).or.lie) then

          ! Index normalization array and estimate minimum
          ! size for normalization
          call index_norm(Atom,LTElines,Atmo,Bfield,GeomI,Red,.False.)

          ! If doing PRD and doing any intensity formal solution
          if (Input%iter_ord.eq.2.and.(lie.or.liter)) then

            ! Call normalization for 1st order profiles in PRD
            call normalizeI_PRD(Atom,Atmo,Frec,Red,ofram)
            if (laborted) goto 1000

            ! Master verbose
            if (gpid.eq.0.and.nA.gt.0) then
              umsg = ' - Normalized 1st order profiles for '// &
                     'PRD (intensity)'
              call verbose
            end if

            ! If storing profiles
            if (IRAM.or.VIRAM) then

              ! If CPU went beyond ram limit
              if (ofram.and.rlimw) then

                ! Verbose
                write(umsg,'(A,i3,A)') ' # Processor',pid, &
                                       ' reached RAM limit with '// &
                                       'input frequency and PRD '// &
                                       'flat allocations.'
                call verbose

                ! Do not notify anymore in this CPU
                rlimw = .False.

              end if ! Went beyond limit and can notify
            end if ! Storing profiles
          end if ! doing PRD and calculating lines
        end if ! Doing any kind of calculation here
      end if ! Intensity formal solution

      ! If iterating or non-dynamic
      if (.not.dyn.or.liter.or.(literJ.and.Input%init_J_bb)) then

        ! Normalize first order profiles
        call normalization(Atom,LTElines,Atmo,Atmo%zeros,GeomI, &
                           Frec,Red,Input,Flgsg,MPID,rlimw,.False.)
        if (laborted) goto 1000

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Normalized/precalculated profiles (intensity)'
          call verbose
        end if

      end if ! Iterating or non dynamic

      ! If not loading and doing J iterations and doing
      ! intensity
      if (.not.lload.and.literJ.and.lio) then

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Iterating radiation field'
          call verbose
        end if

        ! Solve
        call solveJ(Atmo,Atom,LTElines,Cont,Frec, &
                    Red,GeomI,MPID,Input,StokesI,J00C)
        if (laborted) goto 1000

      end if ! Need to do J iterations

      !
      ! If formal solution or emergence and PRD
      !
      if (Input%iter_ord.eq.2.and.(liter.or.lie)) then

        ! Allocate Warr space
        call allocate_WarrI(Atom,Red,GeomI,ofram)

        ! If CPU went above ram and can send messages
        if (IRAM.and.ofram.and.rlimw) then

          ! Verbose
          write(umsg,'(A,i3,A)') ' # Processor',pid, &
                                 ' reached RAM limit with '// &
                                 'redistribution allocations.'
          call verbose

          ! This CPU cannot send it more times
          rlimw = .False.

        end if ! CPU beyond limits
      end if ! PRD with any intensity calculation

      !
      ! Intensity formal solution
      !

      ! If iterating
      if (liter) then

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Starting intensity iterations'
          call verbose
        end if

        ! If we have to, report RAM use
        if (Input%RAMreport) call RAMreport(Input%folder,0,1)

        ! If measuring performance, do it
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

        ! Solve the NLTE problem of the first kind
        call solveI(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                    Red,GeomI,MPID,Input,lload,StokesI, &
                    J00,J00S,J00C,J00P)
        if (laborted) goto 1000

        ! If measuring performance, report now
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

        ! If intensity was AA and pol is AD, treat as AD
        if (PRD.and..not.AV.and.AVI) tbAD = .False.

      end if ! Iterating

      !
      ! Intensity emergence
      !

      ! If calculating emergent intensity
      if (lie) then

        ! Write the solution file
        if (gpid.eq.0) then
          umsg = ' - Saving intensity solution'
          call verbose
        end if

        ! If inversion
        if (run_mode.eq.-1) then

          ! Set solution if we need to keep it
          if (SolF%keep_solution) &
            call setsol(SolF,Flgsg,Bfield, &
                        Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                        J00,J00S,J00C,J00P,.True.)

        ! If synthesis
        else

          ! Write to file
          call writesolI(Input,'NONE',Frec%omega,GeomI,Atom, &
                         Atmo%z,StokesI,J00,J00S,J00C,J00P,.False.)

        end if ! Running mode

        ! Control
        if (laborted) goto 1000

        ! Index geometry
        call geom_index(GeomI,.True.)

        ! Set number of output directions (LOS are treated one at
        ! a time)
        GeomI%njdir = 1

        ! If dynamic
        if (dyn) then

          ! No point in storing Voigt profiles
          VIRAM = .False.

          ! Free existing normalization data
          call free_norm(Red,.True.)

          ! Index normalization array and estimate minimum
          ! size for normalization
          call index_norm(Atom,LTElines,Atmo,Bfield,GeomI,Red,.False.)

        end if

        ! If angle-dependent
        if (.not.AVI) then

          ! No point in storing redistribution functions
          IRAM = .False.

          ! Free redistribution
          call free_warr(Red)

        end if

        ! Predict RAM necessary to solve the RTE
        call emergentI_predict(Atom,Atmo,Red,GeomI,MPID,Input)

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Emergent intensity'
          call verbose
        end if

        ! If we have to, report RAM use
        if (Input%RAMreport) call RAMreport(Input%folder,0,0)

        ! If measuring performance, report now
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

        ! Perform formal solution
        call emergentI(Atom,LTElines,Atmo,Cont,Frec,Red,GeomI, &
                       MPID,Input,StokesI,J00,J00C,SolF)
        if (laborted) goto 1000

        ! If measuring performance, report now
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

      ! If we are doing polarization but we want to keep the
      ! intensity solution file, and we have a solution to store,
      ! in synthesis mode
      elseif (Input%keepIsol.and.lio.and.lp.and.run_mode.ge.0) then

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Saving intensity multi-level solution'
          call verbose
        end if

        ! Write to file
        call writesolI(Input,'NONE',Frec%omega,GeomI,Atom, &
                       Atmo%z,StokesI,J00,J00S,J00C,J00P,.True.)
        if (laborted) goto 1000

      end if ! Emergence of storing solution

      ! If no calling polarization in synthesis
      if (run_mode.ge.0.and.lio.and..not.lp.and..not.lie) then

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Saving solution'
          call verbose
        end if

        ! Write to file as final
        call writesolI(Input,'NONE',Frec%omega,GeomI,Atom, &
                       Atmo%z,StokesI,J00,J00S,J00C,J00P,.False.)

      end if ! If we need to write the solution here afterall

      ! Clean Red structure
1000  call free_red(Red)
      call free_local_geom(GeomI)

      ! Prediction counter
      TRAMc = 0d0

      return

      end subroutine hanle_intensity

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the self-consistent NLTE problem or the second kind
      !! or/and calculate the emergent Stokes profiles\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!      SolF(Solution_F_class): Structure with the solution of
      !!                              the self-consistent problem and
      !!                              the corresponding emergent
      !!                              profiles, contribution function,
      !!                              and height for optical depth
      !!                              equal to one\n
      !!       Cont(Continuum_class): Structure with background
      !!                              opacity data\n
      !!      Rho_old(Rhoc_class(:)): Structure to store the density
      !!                              matrix of the previous
      !!                              iteration\n
      !!    StokesI(double(:,:,:,:)): Intensity\n
      !!            J00(double(:,:)): Mean intensity integrated over
      !!                              the absorption profile\n
      !!           J00S(double(:,:)): Mean intensity integrated over
      !!                              the emission profile\n
      !!           J00C(double(:,:)): Mean intensity with frequency
      !!                              dependence\n
      !!         J00P(double(:,:,:)): Intensity integrals in the
      !!                              photoionization rates\n
      !!   Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!      JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                              integrated over the absorption
      !!                              profile\n
      !!     JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                              integrated over the emission
      !!                              profile\n
      !!     JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                              frequency dependence\n
      !!   JKQ_asym(dcomplex(:,:,:)): Extra asymmetry for the
      !!                              radiation field tensors\n
      !!               rnPh(integer): Allocation size for Stokes in
      !!                              the azimuth dimension\n
      !!             saving(logical): If generating a solution file\n
      !!              lload(logical): If reading an existing solution
      !!                              file\n
      !!                lio(logical): If solving the self-consistent
      !!                              NLTE problem for intensity\n
      !!                lie(logical): If calculating the emergent
      !!                              intensity\n
      !!                 lp(logical): If solving the self-consistent
      !!                              NLTE problem of the 2nd kind\n
      !!                lpe(logical): If calculating the emergent
      !!                              Stokes parameters\n
      !!              rlimw(logical): If we can write a message about
      !!                              reaching the RAM limit\n
      !!              ofram(logical): If we have reached the RAM limit
      subroutine hanle_polarization(Atom,LTElines,Atmo,MPID,Input, &
                                    Geom,Bfield,Frec,Flgsg,SolF, &
                                    Cont,Rho_old, &
                                    StokesI,J00,J00S,J00C,J00P, &
                                    Stokes,JKQ,JKQS,JKQC,JKQ_asym, &
                                    rnPh,saving, &
                                    lload,lio,lp,lpe,rlimw,ofram)
          
      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(inout):: Flgsg
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(inout):: Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      type(Solution_F_class), intent(inout):: SolF
      type(Continuum_class), intent(in):: Cont
      type(Rhoc_class), dimension(:), &
                        allocatable, intent(inout):: Rho_old
      logical, intent(in):: saving,lload,lio,lp,lpe
      logical, intent(inout):: rlimw,ofram
      integer, intent(in):: rnPh
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: StokesI
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(inout):: J00P
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Stokes
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQC
      complex(kind=8), dimension(:,:,:), &
                       allocatable, intent(in):: JKQ_asym

      ! Local

      type(Red_class):: Red

      logical:: l1,liter,lcorr,field

      integer:: ia


      ! Initialize redistribution pointers
      nullify(Red%dzao,Red%rzao,Red%pzao,Red%zao)

      ! Initialize memory in Red structure
      FRAMc = FRAMc + 1d-6*sizeof(Red)

      ! Check if iterating
      liter = lp.and.nA.gt.0.and.(Input%iter_max.ge.Input%iter_min)

      ! Check if correcting JKQ
      lcorr = (lio.or.(lload.and..not.allocated(Stokes)))

      ! Check if there is a magnetic field
      field = maxval(Bfield%Bstrength).gt.TINYB

      ! Check RT axiality
      RTaxial = axial.and. &
                maxval(Bfield%Bstrength).le.TINYB

      !
      ! Prepare geometrical tensors
      !

      ! If PRD, or iterating, or correcting JKQ, get B geometrical
      ! tensors in quadrature
      if (field.and.(lp.or.lcorr.or.Input%iter_ord.eq.2)) &
        call setTB(Geom,Flgsg,Bfield)

      !
      ! Diagonalize Hamiltonian
      !

      ! If there is magnetic field
      if (field) then

        ! For each atom
        do ia=1,nA

          ! Diagonalize
          call diagon(Atom(ia),Bfield,Input%zeeman_mode,Flgsg)

          ! Compute transition strength in energy representation
          call strength_ev(Atom(ia),Bfield)

        end do ! Atoms

        ! Master verbose
        if(gpid.eq.0.and.nA.gt.0) then
          umsg = ' - Hamiltonian diagonalized'
          call verbose
          umsg = ' - Dipole strengths in energy eigenbasis calculated'
          call verbose
        end if

      ! No magnetic field
      else

        ! For each atom
        do ia=1,nA

          ! Initialize dimensions coming from diagonalization even if
          ! trivial
          call diagon_B0(Atom(ia))

        end do ! Atoms

      end if ! Magnetic field presence

      !
      ! Normalizations
      !

      ! If iterating or doing corrections, normalize
      if (liter.or.lcorr.or.lpe) then

        !
        ! Set number of output directions
        !

        ! Iterating or correcting
        if (liter.or.lcorr) then

          ! Number of output directions is quadrature
          Geom%njdir = Geom%nTh*Geom%nPh

          ! Index geometry
          call geom_index(Geom,.False.)

        ! Only emergence
        else

          ! Number of output directions given by LOS
          Geom%njdir = 1

          ! Index geometry
          call geom_index(Geom,.True.)

        end if

        ! If doing PRD and calculating Stokes profiles
        if (Input%iter_ord.eq.2.and.(lpe.or.liter)) then

          ! Index redistribution data arrays, allocate for
          ! input frequencies, and estimate minimum memory
          ! for PRD normalization
          call index_red(Atom,Red,Bfield%Bstrength,.True.)

          ! Get input frequency axis
          call omegabuildin(Frec,Red,Atom,Atmo, &
                            Bfield%Bstrength,Input,MPID, &
                            ofram)
          if (laborted) goto 1000

          ! Master verbose
          if (gpid.eq.0.and.nA.gt.0) then
            umsg = ' - Input frequency axis initialized'
            call verbose
          end if ! Master
        end if ! PRD

        ! Predict RAM necessary for RT
        if (liter) &
          call solve_predict(Atom,Frec,Red,Geom,MPID,Input)

        ! Index normalization array and estimate minimum
        ! size for normalization
        call index_norm(Atom,LTElines,Atmo,Bfield,Geom,Red,.True.)

        ! If doing PRD and calculating Stokes profiles
        if (Input%iter_ord.eq.2.and.(lpe.or.liter)) then

          ! Call normalization for 1st order profiles in PRD
          call normalize_PRD(Atom,Atmo,Bfield%Bstrength, &
                             Frec,Red,ofram)
          if (laborted) goto 1000

          ! Master verbose
          if (gpid.eq.0.and.nA.gt.0) then
            umsg = ' - Normalized 1st order profiles for PRD'
            call verbose
          end if

          ! If storing profiles
          if (PRAM.or.VPRAM) then

            ! If CPU went beyond ram limit
            if (ofram.and.rlimw) then

              ! Verbose
              write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                  ' reached RAM limit with redistribution ', &
                  'allocation'
              call verbose

              ! Do not notify anymore in this CPU
              rlimw = .False.

            end if ! Went beyond limit and can notify
          end if ! Storing profiles
        end if ! PRD
      end if ! Iterating or doing corrections

      ! If iterating or non dynamic
      if (.not.dyn.or.(liter.or.lcorr)) then

        ! Normalize first order profiles
        call normalization(Atom,LTElines,Atmo,Bfield%Bstrength,Geom, &
                           Frec,Red,Input,Flgsg,MPID,rlimw,.True.)
        if (laborted) goto 1000

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Normalized/precalculated profiles'
          call verbose
        end if

      end if ! Iterating or non dynamic

#ifdef DEBUGRHOKQ
      if (pid.eq.0) call dump_rho(Atom,Input%folder,-3)
#endif
      ! If we need to correct the JKQ to multi-term from the
      ! multi-level solution
      if (lcorr) then

        ! Save the NCHLT variable and set it temporally to false
        l1 = NCHLT
        if (NCHLT) NCHLT = .False.
        ! Correct the JKQ tensors
        call JKQgen(Atom,Rho_old,Atmo,Frec,Red,Geom, &
                    MPID,Input,Flgsg,Input%Pcorr,Bfield, &
                    rnPh,StokesI,J00,J00S,J00C, &
                    Stokes,JKQ,JKQS,JKQC,J00P)
        if (laborted) goto 1000

        ! Only if synthesis
        if (run_mode.ne.-1) then

          ! Master verbose
          if (gpid.eq.0) then
            umsg = ' - Saving intensity solution'
            call verbose
          end if

          ! Write the solution file
          call writesol(Input,'NONE',Frec%omega,Geom,Flgsg,Bfield, &
                        Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)

        end if ! Synthesis

        ! Restore NHCLT variable
        NCHLT = l1

      end if ! Correcting JKQ tensors

      ! If doing PRD
      if (Input%iter_ord.eq.2) then

        ! If non-coherent lower term, check the critical field
        ! condition
        if (NCHLT) call check_nchlt(Atom,JKQ,Bfield)

      end if ! PRD

      ! Any PRD polarization calculation
      if (Input%iter_ord.eq.2.and.(lp.or.lpe)) then

        ! Allocate Warr data
        call allocate_Warr(Atom,Red,Geom,Bfield%Bstrength,ofram)

        ! If CPU went beyond ram limit
        if (PRAM.and.ofram.and.rlimw) then

          ! Verbose
          write(umsg,'(A,i3,A)') ' # Processor',pid, &
                                 ' reached RAM limit with '// &
                                 'redistribution allocations.'
          call verbose

          ! This CPU cannot send it more times
          rlimw = .False.

        end if ! CPU beyond limits
      end if ! PRD with any polarization calculation

      !
      ! Formal Solution
      !

      ! If solving NLTE problem of the second kind
      if (lp) then

        ! If actually iterating
        if (liter) then

          ! Master verbose
          if (gpid.eq.0) then
            umsg = ' - Starting iterations'
            call verbose
          end if

          ! If we have to report RAM use, do it now
          if (Input%RAMreport) call RAMreport(Input%folder,1,1)

          ! If measuring performance, report now
          if (Input%g_perf.and.pid.eq.0) &
            call report_time(Input%folder,Input%ID,.True.)

          ! Solve the NLTE problem of the second kind
          call solve(Atom,LTElines,Rho_old,Atmo,Cont,Frec, &
                     Red,Bfield,Geom,MPID,Input,Flgsg, &
                     JKQ_asym,Stokes,JKQ,JKQS,JKQC)
          if (laborted) goto 1000

          ! If measuring performance, report it now
          if (Input%g_perf.and.pid.eq.0) &
            call report_time(Input%folder,Input%ID,.True.)

          ! Master verbose
          if (gpid.eq.0) then
            umsg = ' - Saving solution'
            call verbose
          end if

          ! If saving the solution
          if (saving) then

            ! If inversion
            if (run_mode.eq.-1) then

              ! Keep the solution if needed
              if (SolF%keep_solution) &
                call setsol(SolF,Flgsg,Bfield, &
                            Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                            J00,J00S,J00C,J00P,.False.)

            ! Synthesis
            else

              ! Write the solution to file
              call writesol(Input,'NONE',Frec%omega,Geom,Flgsg, &
                            Bfield,Atom,Atmo%z,Stokes,JKQ,JKQS,JKQC)

            end if ! Inversion/synthesis
          end if ! Saving solution

        ! No iterations, but inversion, keeping solution, and allowed
        ! to do so
        else if (run_mode.eq.-1.and.SolF%keep_solution.and. &
                 saving) then

          ! Set solution
          call setsol(SolF,Flgsg,Bfield, &
                      Atom,Stokes,JKQ,JKQS,JKQC,StokesI, &
                      J00,J00S,J00C,J00P,.False.)

        end if ! There are actual iterations
      endif ! Formal solution for polarization


      !
      ! Stokes emergence
      !

      ! Need to generate emergent Stokes profiles
      if (lpe) then

        ! Index geometry
        call geom_index(Geom,.True.)

        ! Set number of output directions
        Geom%njdir = 1

        ! If dynamic
        if (dyn) then

          ! No point in storing Voigt profiles
          VPRAM = .False.

          ! Free existing normalization data
          call free_norm(Red,.True.)

          ! Index normalization array and estimate minimum
          ! size for normalization
          call index_norm(Atom,LTElines,Atmo,Bfield,Geom,Red,.True.)

        end if

        ! If angle-dependent
        if (.not.AV) then

          ! No point in storing redistribution functions
          PRAM = .False.

          ! Free redistribution
          call free_warr(Red)

        end if

        ! Predict RAM necessary to solve the RTE
        call emergent_predict(Atom,Red,Geom,MPID,Input)

        ! Master verbose
        if (gpid.eq.0) then
          umsg = ' - Emergent Stokes'
          call verbose
        end if

        ! If we have to report RAM use, do it now
        if (Input%RAMreport) call RAMreport(Input%folder,1,0)

        ! If measuring performance, report it now
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

        ! Perform formal solution
        call emergent(Atom,LTElines,Atmo,Cont,Frec,Red,Bfield, &
                      Geom,MPID,Input,Flgsg,JKQ_asym,Stokes, &
                      JKQ,JKQC,SolF)
        if (laborted) goto 1000

        ! If measuring performance, report it now
        if (Input%g_perf.and.pid.eq.0) &
          call report_time(Input%folder,Input%ID,.True.)

      end if ! Energence

      ! Control errors
      call control

      !
      ! Clean memory
      !

      ! Clean Red structures and geometrical tensors
1000  call free_red(Red)
      call free_local_geom(Geom)

      ! Prediction counter
      TRAMc = 0d0

      return

      end subroutine hanle_polarization

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare solution structure for the inversion problem\n
      !> Allocate the arrays to store the solution of the forward
      !! problem in the inversion mode and configure the synthesis
      !! mode for the inversion code, depending on if the synthesis is
      !! a trial or a response function\n
      !!  SolF(Solution_F_class): Structure with the solution of the
      !!                          self-consistent problem and the
      !!                          corresponding emergent profiles,
      !!                          contribution function, and height
      !!                          for optical depth equal to one\n
      !!      Input(Input_class): Structure with configuration data\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!   GeomI(Geometry_class): Structure with geometric data for
      !!                          the intensity problem\n
      !!    Geom(Geometry_class): Structure with geometric data\n
      !!             RF(logical): If solving the NLTE problem to
      !!                          calculate the response function
      subroutine prepare_buffers(SolF,Input,Atom,GeomI,Geom,RF)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Geometry_class), intent(in):: GeomI, Geom
      type(Input_class), intent(inout):: Input
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: RF

      ! Local

      integer:: ia


      !
      ! Configure run
      !

      ! If run if to compute a response function
      if (RF) then

        ! If initializing from existing density matrix
        if (Input%popuinit) then

          ! Read and calculate mode
          Input%mode = 'B'

        ! If doing from scratch
        else

          ! Write mode
          Input%mode = 'W'

        end if ! If using previous solution

        ! Do not keep the full solution
        SolF%keep_solution = .False.

      ! If the run is a trial
      else

        ! Start from scratch
        Input%mode = 'W'

        ! And keep the full solution
        SolF%keep_solution = .True.

      end if ! Computing response function or a trial

      ! Only master continues from here
      if (pid.gt.0) return

      !
      ! Prepare buffers
      !

      ! If not initialized the space for solutions
      if (SolF%no_initialized) then

          ! If we need to initialize intensity
        if (Input%Type_Inversion.eq.0.or. &
            (Input%Type_Inversion.eq.3.and.Input%force.eq.'I').or. &
            (Input%Type_Inversion.eq.4.and.Input%force.eq.'I')) then

          ! If keeping Stokes parameters
          if (KSTK) then

            ! Allocate
            allocate(SolF%i_StkI(nfreq,GeomI%nPh,GeomI%nTh,nz))
            allocate(SolF%i_StkI_b(nfreq,GeomI%nPh,GeomI%nTh,nz))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_StkI_t(nfreq,GeomI%nPh,GeomI%nTh,nz))

            ! Count memory
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_StkI)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_StkI_b)
            if (Input%LM_Method.eq.1) &
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_StkI_t)
          end if

          ! Allocate
          allocate(SolF%i_J00C(nfreq,nz))
          allocate(SolF%i_J00C_b(nfreq,nz))
          if (Input%LM_Method.eq.1) then
            allocate(SolF%i_J00C_t(nfreq,nz))
          end if

          ! Count memory
          SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00C)
          SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00C_b)
          if (Input%LM_Method.eq.1) then
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00C_t)
          end if

          ! If active atoms
          if (nA.gt.0) then

            ! Allocate
            allocate(SolF%i_J00(nxt,nz))
            allocate(SolF%i_J00P(nxphot,2,nz))
            allocate(SolF%i_J00_b(nxt,nz))
            allocate(SolF%i_J00P_b(nxphot,2,nz))
            if (Input%LM_Method.eq.1) then
              allocate(SolF%i_J00_t(nxt,nz))
              allocate(SolF%i_J00P_t(nxphot,2,nz))
            end if

            ! Count memory
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00P)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00_b)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00P_b)
            if (Input%LM_Method.eq.1) then
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00_t)
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_J00P_t)
            end if

            ! Allocate rhoes
            allocate(SolF%i_rhoes(na))
            allocate(SolF%i_rhoes_b(na))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_rhoes_t(na))

            ! Count memory
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_b)
            if (Input%LM_Method.eq.1) &
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_t)

            ! For each atom
            do ia=1,nA

              ! Allocate
              allocate(SolF%i_rhoes(ia)%rho(Atom(ia)%nlevel,nz))
              allocate(SolF%i_rhoes_b(ia)%rho(Atom(ia)%nlevel,nz))
              if (Input%LM_Method.eq.1) &
                allocate(SolF%i_rhoes_t(ia)%rho(Atom(ia)%nlevel,nz))

              ! Count memory
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes(ia)%rho)
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_b(ia)%rho)
              if (Input%LM_Method.eq.1) &
                SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_t(ia)%rho)

            end do ! Atoms

          end if ! If active atoms
        end if ! Intensity

        ! If doing polarization
        if (Input%Type_inversion.eq.1.or. &
            Input%Type_inversion.eq.2.or. &
            (Input%Type_Inversion.eq.3.and.Input%force.eq.'N').or. &
            (Input%Type_Inversion.eq.4.and.Input%force.eq.'N')) then

          ! If keeping Stokes paramaters
          if (KSTK) then

            ! Allocate
            allocate(SolF%i_Stk(0:3,nfreq,Geom%nPh,Geom%nTh,nz))
            allocate(SolF%i_Stk_b(0:3,nfreq,Geom%nPh,Geom%nTh,nz))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_Stk_t(0:3,nfreq,Geom%nPh,Geom%nTh,nz))

            ! Count memory
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_Stk)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_Stk_b)
            if (Input%LM_Method.eq.1) &
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_Stk_t)

          end if ! Keeping Stokes parameters

          ! Allocate
          allocate(SolF%i_JKQC(-2:2,0:2,nfreq,nz))
          allocate(SolF%i_JKQC_b(-2:2,0:2,nfreq,nz))
          if (Input%LM_Method.eq.1) then
            allocate(SolF%i_JKQC_t(-2:2,0:2,nfreq,nz))
          end if

          ! Count memory
          SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQC)
          SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQC_b)
          if (Input%LM_Method.eq.1) then
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQC_t)
          end if

          ! If active atoms
          if (nA.gt.0) then

            ! Allocate
            allocate(SolF%i_JKQ(-2:2,0:2,nxtran,nz))
            allocate(SolF%i_JKQS(-2:2,0:2,nxtran,nz))
            allocate(SolF%i_JKQ_b(-2:2,0:2,nxtran,nz))
            allocate(SolF%i_JKQS_b(-2:2,0:2,nxtran,nz))
            if (Input%LM_Method.eq.1) then
              allocate(SolF%i_JKQ_t(-2:2,0:2,nxtran,nz))
              allocate(SolF%i_JKQS_t(-2:2,0:2,nxtran,nz))
            end if

            ! Count memory
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQ)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQS)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQ_b)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQS_b)
            if (Input%LM_Method.eq.1) then
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQ_t)
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_JKQS_t)
            end if

            ! Allocate rhoes
            allocate(SolF%i_rhoes(na))
            allocate(SolF%i_rhoes_b(na))
            if (Input%LM_Method.eq.1) &
              allocate(SolF%i_rhoes_t(na))

            ! Count memory
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes)
            SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_b)
            if (Input%LM_Method.eq.1) &
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_t)

            ! For each atom
            do ia=1,nA

              ! Allocate
              allocate(SolF%i_rhoes(ia)%crho(Atom(ia)%ndim,nz))
              allocate(SolF%i_rhoes_b(ia)%crho(Atom(ia)%ndim,nz))
              if (Input%LM_Method.eq.1) &
                allocate(SolF%i_rhoes_t(ia)%crho(Atom(ia)%ndim,nz))

              ! Count memory
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes(ia)%crho)
              SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_b(ia)%crho)
              if (Input%LM_Method.eq.1) &
                SRAMc = SRAMc + 1d-6*sizeof(SolF%i_rhoes_t(ia)%crho)

            end do ! Atoms

          end if ! Active atoms
        end if ! Polarization

        ! Flag initialized
        SolF%no_initialized = .False.

      end if ! Initialize

      !
      ! Free space for emergence outputs
      !

      ! Stokes
      if (allocated(SolF%e_Stk)) then
        SRAMc = SRAMc - 1d-6*sizeof(SolF%e_Stk)
        deallocate(SolF%e_Stk)
      end if

      ! Contribution functions
      if (allocated(SolF%e_Ctr)) then
        SRAMc = SRAMc - 1d-6*sizeof(SolF%e_Ctr)
        deallocate(SolF%e_Ctr)
      end if

      ! Heights for optical depth equal to one
      if (allocated(SolF%e_tau1)) then
        SRAMc = SRAMc - 1d-6*sizeof(SolF%e_tau1)
        deallocate(SolF%e_tau1)
      end if

      end subroutine prepare_buffers

!#####################################################################
!#####################################################################
!#####################################################################

      end module hanle_mod
