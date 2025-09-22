      !> Pre-calculation of emiss2ord integrals
      module comovingprd_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     08/10/2024
!  Last version:
!     22/09/2025 V4.1.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     22/09/2025:    V4.1.2 - Bugfix: Solved the deadlock issue (TdPA)
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
!  predict_emiss
!    Predict the amount of memory space that will be needed by the
!  second order emissivity
!
!  predict_emissI
!    Predict the amount of memory space that will be needed by the
!  second order emissivity for intensity
!
!  initialize_emiss
!    Allocate space needed for the second order emissivity
!
!  initialize_emissI
!    Allocate space needed for the second order emissivity for
!  intensity
!
!  comoving_emiss2ord
!    Compute second order emissivity
!
!  comoving_emissI2ord
!    Compute second order emissivity for intensity
!
!  share_emiss
!    Share emissivity data in one height of the model atmosphere
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use inter_mod
      use parameters_mod , only : c , TINYB , TINYVEL
      use rtcoeffiaux_mod
      use rtcoeffaux_mod
      use types_mod

      ! Deadlock preventers
      integer, parameter:: W_minimal = 4
      integer, parameter:: W_nominal = 8
      double precision, parameter:: mem_limit = 50d6

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Predict the amount of memory space that will be needed by the
      !! second order emissivity\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data
      subroutine predict_emiss(Atom,Geom,Red)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Red_class), intent(in):: Red

      ! Local

      integer:: ia,jtran,indx,ndir,nf,nf2,lnz

      double precision:: MTRAMc,lTRAMc,lTRAMbc


      ! If no atoms, leave
      if (nA.lt.1) return

      ! Master is not going to use these, so leave if master in MPI
      if (pid.eq.0.and.nproc.gt.1) return

      ! If no PRD, how did you even get here
      if (.not.PRD) return

      ! Number of output directions
      ndir = Geom%njdir

      ! Number of heights
      lnz = Rz1_PRD - Rz0 + 1

      ! Initialize maximum
      MTRAMc = 0d0

      ! For each atom
      do ia=1,nA

        ! For each transition
        do jtran=1,Atom(ia)%ntran

          ! If no PRD, skip
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Get zao index
          indx = Red%izao(jtran,ia,Rz0)

          ! Initialize size
          nf = 0
          nf2 = 0

          ! Get size
          if (Red%ao(indx)%tgf1.ge.Red%ao(indx)%tgf0) &
            nf = Red%ao(indx)%tgf1 - Red%ao(indx)%tgf0 + 1
          if (Red%ao(indx)%gf1.ge.Red%ao(indx)%gf0) &
            nf2 = Red%ao(indx)%gf1 - Red%ao(indx)%gf0 + 1

          ! Get size for sender
          lTRAMc = 8d-6*dble(ndir*4*Red%ao(indx)%nn(pid))

          ! If MPI, add receiver
          if (nproc.gt.1) &
            lTRAMc = 8d-6*dble(ndir*4*lnz*nf)

          ! Check auxiliar top
          lTRAMbc = lTRAMc + 8d-6*dble(nf*8*ndir)
          if (lTRAMbc.gt.MTRAMc) MTRAMc = lTRAMbc

          ! Check auxiliar bottom
          lTRAMbc = lTRAMc + 8d-6*dble(3*nf + nf2)
          if (lTRAMbc.gt.MTRAMc) MTRAMc = lTRAMbc

          ! Size for collection
          ! The 4 assumes that will be using splines to
          ! interpolate
          if (nf.gt.0) then
            lTRAMc = 8d-6*dble(nf*(lnz*ndir*8 + 4))
            if (lTRAMc.gt.MTRAMc) MTRAMc = lTRAMc
          end if

          ! If this CPU does not have frequencies
          ! in this line, skip
          if (Atom(ia)%fflag(jtran)%absent) cycle

          ! Skip if transition is outside limits
          if (nf2.le.0) cycle

          ! Size for storage
          TRAMc = TRAMc + 8d-6*dble(ndir*4*nf2*lnz)

        end do ! Transition
      end do ! Atom

      ! Add the maximum working package
      TRAMc = TRAMc + MTRAMc

      end subroutine predict_emiss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Predict the amount of memory space that will be needed by the
      !! second order emissivity for intensity\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data
      subroutine predict_emissI(Atom,Geom,Red)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Red_class), intent(in):: Red

      ! Local

      integer:: ia,jtran,indx,fjtran,ffjtran
      integer:: ndir,nf,nf2,lnz

      double precision:: MTRAMc,lTRAMc,lTRAMbc


      ! If no atoms, leave
      if (nA.lt.1) return

      ! Master is not going to use these, so leave if master MPI
      if (pid.eq.0.and.nproc.gt.1) return

      ! If no PRD, how did you even get here
      if (.not.PRD) return

      ! Number of output directions
      ndir = Geom%njdir

      ! Number of heights
      lnz = Rz1_PRD - Rz0 + 1

      ! Initialize maximum
      MTRAMc = 0d0

      ! If dynamic
      if (dyn) then

        ! Number of output directions in geometry
        ndir = Geom%njdir

      ! No velocity
      else

        ! If angle-averaged
        if (AVI) then

          ! Isotropic
          ndir = 1

        ! If angle-dependent
        else

          ! In general non-isotropic
          ndir = Geom%njdir

        end if ! AA/AD
      end if ! Dynamic

      ! For each atom
      do ia=1,nA

        ! For each transition
        do jtran=1,Atom(ia)%ntran

          ! If no PRD, skip
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

            ! Get zao index
            indx = Red%izao(ffjtran,ia,Rz0)

            ! Initialize size
            nf = 0
            nf2 = 0

            ! Get size
            if (Red%ao(indx)%tgf1.ge.Red%ao(indx)%tgf0) &
              nf = Red%ao(indx)%tgf1 - Red%ao(indx)%tgf0 + 1
            if (Red%ao(indx)%gf1.ge.Red%ao(indx)%gf0) &
              nf2 = Red%ao(indx)%gf1 - Red%ao(indx)%gf0 + 1

            ! Get size for sender
            lTRAMc = 8d-6*dble(ndir*2*Red%ao(indx)%nn(pid))

            ! If MPI, add receiver
            if (nproc.gt.1) &
              lTRAMc = 8d-6*dble(ndir*2*lnz*nf)

            ! Check auxiliar top
            lTRAMbc = lTRAMc + 8d-6*dble(nf*(1 + 2*ndir))
            if (lTRAMbc.gt.MTRAMc) MTRAMc = lTRAMbc

            ! Check auxiliar bottom
            lTRAMbc = lTRAMc + 8d-6*dble(3*nf + nf2)
            if (lTRAMbc.gt.MTRAMc) MTRAMc = lTRAMbc

            ! If this CPU does not have frequencies
            ! in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Skip if transition is outside limits
            if (nf2.le.0) cycle

            ! Size for storage
            TRAMc = TRAMc + 8d-6*dble(ndir*2*nf2*lnz)

          end do ! Transition level
        end do ! Transition term
      end do ! Atom

      ! Add the maximum working package
      TRAMc = TRAMc + MTRAMc

      end subroutine predict_emissI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate space needed for the second order emissivity\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data
      subroutine initialize_emiss(Atom,Geom,Red)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Red_class), intent(inout):: Red

      ! Local

      integer:: ia,jtran,iz,indx,jndx,if0l2,if1l2,ndir


      ! If no atoms, leave
      if (nA.lt.1) return

      ! Master is not going to use these, so leave if master MPI
      if (pid.eq.0.and.nproc.gt.1) return

      ! Number of output directions
      ndir = Geom%njdir

      ! For each height
      do iz=Rz0,Rz1_PRD

        ! For each atom
        do ia=1,nA

          ! For each transition
          do jtran=1,Atom(ia)%ntran

            ! If no PRD, skip
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! If this CPU does not have frequencies
            ! in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Get zao index
            jndx = Red%izao(jtran,ia,Rz0)
            indx = Red%izao(jtran,ia,iz)

            ! Get indexes
            if0l2 = Red%ao(jndx)%gf0
            if1l2 = Red%ao(jndx)%gf1

            ! Skip if outside range
            if (if0l2.gt.if1l2) cycle

            ! Allocate emiss
            allocate(Red%zao(indx)%eps20(ndir,if0l2:if1l2))
            allocate(Red%zao(indx)%eps21(ndir,if0l2:if1l2))
            allocate(Red%zao(indx)%eps22(ndir,if0l2:if1l2))
            allocate(Red%zao(indx)%eps23(ndir,if0l2:if1l2))

          end do ! Transition
        end do ! Atom
      end do ! Height

      end subroutine initialize_emiss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate space needed for the second order emissivity for
      !! intensity\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data
      subroutine initialize_emissI(Atom,Atmo,Geom,Red)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Red_class), intent(inout):: Red

      ! Local

      logical:: lvel

      integer:: ia,jtran,iz,indx,jndx,fjtran,ffjtran,if0l2,if1l2,ndir

      double precision:: vel


      ! If no atoms, leave
      if (nA.lt.1) return

      ! Master is not going to use these, so leave if master MPI
      if (pid.eq.0.and.nproc.gt.1) return

      ! For each height
      do iz=Rz0,Rz1_PRD

        ! Velocity
        if (dyn) then

          ! Compute maximum velocity
          vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                     Atmo%vy(iz)*Atmo%vy(iz) + &
                     Atmo%vz(iz)*Atmo%vz(iz))

          ! Check if large enough
          lvel = vel.gt.TINYVEL

          ! Large enough velocity
          if (lvel) then

            ! Number of output directions in geometry
            ndir = Geom%njdir

          ! No velocity
          else

            ! If angle-averaged
            if (AVI) then

              ! Isotropic
              ndir = 1

            ! If angle-dependent
            else

              ! In general non-isotropic
              ndir = Geom%njdir

            end if ! AA/AD
          end if ! Velocity magnitude

        ! No velocity
        else

          ! Trivial
          vel = 0d0
          lvel = .False.

          ! If angle-averaged
          if (AVI) then

            ! Isotropic
            ndir = 1

          ! If angle-dependent
          else

            ! In general non-isotropic
            ndir = Geom%njdir

          end if ! AA/AD
        end if ! Dynamic

        ! For each atom
        do ia=1,nA

          ! For each transition
          do jtran=1,Atom(ia)%ntran

            ! If no PRD, skip
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! If this CPU does not have frequencies
            ! in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! For each FS transition
            do fjtran=1,Atom(ia)%fst(jtran)%nt

              ! Get the sequential index of this FS transition
              ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

              ! Get zao index
              jndx = Red%izao(ffjtran,ia,Rz0)
              indx = Red%izao(ffjtran,ia,iz)

              ! Get indexes
              if0l2 = Red%ao(jndx)%gf0
              if1l2 = Red%ao(jndx)%gf1

              ! Skip if outside range
              if (if0l2.gt.if1l2) cycle

              ! Allocate emiss
              allocate(Red%zao(indx)%eps20(ndir,if0l2:if1l2))
              allocate(Red%zao(indx)%rpf(ndir,if0l2:if1l2))

            end do ! Transition level
          end do ! Transition term
        end do ! Atom
      end do ! Height

      end subroutine initialize_emissI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute second order emissivity for intensity\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            input frequency data,
      !!                            redistribution function data, and
      !!                            profile or normalization data\n
      !!        Flgsg(Fctsg_class): Structure with factorials,
      !!                            signs, and J-symbols\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !! Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!      JKQa(dcmplex(:,:,:)): Extra asymmetry for the radiation
      !!                            field tensors\n
      !!     JKQ(dcmplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over the absorption profile\n
      !!    JKQC(dcmplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence
      !!         int_mode(integer): Interpolation mode\n
      !!              lth(integer): Current index for polar direction
      !!                            if doing los formal solution\n
      !!              lph(integer): Current index for azimuth
      !!                            direction if doing los formal
      !!                            solution\n
      !!              los(logical): If performing last formal
      !!                            solution
      subroutine comoving_emiss2ord(Atom,Atmo,Geom,Frec,Red, &
                                    Flgsg,Bfield,Stokes,JKQa, &
                                    JKQ,JKQC,lth,lph,int_mode,los)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: los
      integer, intent(in):: lth,lph,int_mode
      double precision, dimension(0:3,nfreq,Geom%nPh, &
                                  Geom%nTh,giz0:giz1), &
                        intent(in):: Stokes
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(in):: JKQC

      ! Local

#ifndef oldmpi
      type(MPI_request), dimension(:), allocatable:: requests
#endif

      logical:: lfield,lvel,first

      integer:: iz0,iz1,Mif0,Mif1,jndx,jj,kk,ll,iran
      integer:: iz,ia,jtran,itermu,itermf,indx,t0,t1,i0,i1
      integer:: if0l2,if1l2,if0tl2,if1tl2,if0Il2,if1Il2,ifreq,nf,nfl
      integer:: ndir,idir,ierr,iph,ith,nth,nph,lnz
#ifndef oldmpi
      integer:: ntest,W_max,W_load,lgiz
      integer, dimension(:), allocatable:: done_index
#endif
      integer, dimension(:), allocatable:: nsend
      integer, dimension(:,:), allocatable:: nf_s,disp

      double precision:: vel,DwT,Dw,ct,st,cc,sc,vfac
      double precision, dimension(:), allocatable:: omega
      double precision, dimension(:,:), allocatable:: splin
      double precision, dimension(:,:), allocatable:: eps0
      double precision, dimension(:,:), allocatable:: eps1
      double precision, dimension(:,:), allocatable:: eps2
      double precision, dimension(:,:), allocatable:: eps3
      double precision, dimension(:,:,:), &
                        allocatable, target:: eps20123_s

      ! Pointers
      type(Reda_class), pointer:: p_fed
      type(Redb_class), pointer:: p_red
      type(Redb2_class), pointer:: p_rwarr
      type(Prof_class), pointer:: p_Norm
      double precision, dimension(:,:), pointer:: eps20
      double precision, dimension(:,:), pointer:: eps21
      double precision, dimension(:,:), pointer:: eps22
      double precision, dimension(:,:), pointer:: eps23
      double precision, dimension(:,:,:), pointer:: eps20123_r
      complex(kind=8), dimension(:,:,:,:), pointer:: TKQo
      type(Redb2_class), target:: p_dummy


      ! Nullify pointers
      nullify(p_fed,p_red,p_rwarr,p_Norm,eps20123_r,TKQo)

      ! If MPI
      if (nproc.gt.1) then
#ifndef oldmpi
        ! Initialize requests
        allocate(requests(Rz0:Rz1_PRD))
        allocate(done_index(Rz0:Rz1_PRD))
        ntest = 0
        do iz=Rz0,Rz1_PRD
          requests(iz) = MPI_REQUEST_NULL
          done_index(iz) = 0
        end do
        ! Allocate frequency size for output
        allocate(nsend(Rz0:Rz1_PRD))
        allocate(nf_s(0:nproc-1,Rz0:Rz1_PRD))
        allocate(disp(0:nproc-1,Rz0:Rz1_PRD))
#else
        ! Allocate frequency size for output
        allocate(nsend(1)
        allocate(nf_s(0:nproc-1,1))
        allocate(disp(0:nproc-1,1))
#endif
      end if ! MPI

      ! If no atoms, leave
      if (nA.lt.1) return

      ! Dummy
      if (.not.PRAM) p_rwarr => p_dummy

      ! Number of heights
      lnz = Rz1_PRD - Rz0 + 1

      ! Directions
      ndir = Geom%njdir

      ! Number of directions
      if (los) then
        nth = 1
        nph = 1
      else
        nth = Geom%nTh
        nph = Geom%nPh
      end if

      ! For each atom
      do ia=1,nA

        ! For each output transition
        do jtran=1,Atom(ia)%ntran

          ! Skip if no PRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          itermf = Atom(ia)%fst(jtran)%iterml

          ! Transition ranges
          t0 = Atom(ia)%tshift + 1
          t1 = t0 + Atom(ia)%ntran - 1

          ! Get index
          jndx = Red%izao(jtran,ia,Rz0)

          ! Point to relevant variable
          p_fed => Red%ao(jndx)

          ! Get limits for this transition
          if0l2 = p_fed%gf0
          if1l2 = p_fed%gf1
          if0tl2 = p_fed%tgf0
          if1tl2 = p_fed%tgf1

          ! If no output, skip
          if ((if1tl2-if0tl2).lt.0) cycle

          ! Initialize first communication
          first = .True.

          ! Allocate emissivity to collect
          allocate(eps0(ndir,if0tl2:if1tl2))
          allocate(eps1(ndir,if0tl2:if1tl2))
          allocate(eps2(ndir,if0tl2:if1tl2))
          allocate(eps3(ndir,if0tl2:if1tl2))
          allocate(eps20(ndir,if0tl2:if1tl2))
          allocate(eps21(ndir,if0tl2:if1tl2))
          allocate(eps22(ndir,if0tl2:if1tl2))
          allocate(eps23(ndir,if0tl2:if1tl2))

          ! Number of frequencies
          nf = if1tl2 - if0tl2 + 1

          ! If doing MPI
          if (nproc.gt.1) then

            ! Allocate emissivity to receive
            allocate(eps20123_r(ndir,4,nf*lnz))
#ifndef oldmpi
            ! Initialize requests
            ntest = 0
            do iz=Rz0,Rz1_PRD
              requests(iz) = MPI_REQUEST_NULL
              done_index(iz) = 0
            end do

            ! Calculate maximum working load
            W_max = min(W_nominal,max(W_minimal, &
                                    nint(mem_limit/ &
                        dble(nf*sizeof(eps20123_r(:,:,1))))))

            ! Initialize last height index sent and load
            lgiz = Rz0 - 1
            W_load = 0

            ! Initialize MPI data
            nsend = 0
            nf_s = 0
            disp = 0
#endif
          end if ! MPI

          ! Has work to do
          if (p_fed%nn(pid).gt.0) then

            ! Get height limits
            iz0 = Rz0 + (p_fed%Mi0(pid)-1)/p_fed%nfreq
            iz1 = Rz0 + (p_fed%Mi1(pid)-1)/p_fed%nfreq

            ! Allocate sender
            allocate(eps20123_s(ndir,4,p_fed%nn(pid)))

            ! Initialize rolling indexes
            jj = 0

          ! Does not have work
          else

            ! Ensure no looping
            iz0 = 0
            iz1 = -1

          end if ! Sanity check

          ! Initialize index
          ll = 0

          ! For each height
          do iz=Rz0,Rz1_PRD

            ! If within limits
            if (iz.ge.iz0.and.iz.le.iz1) then

              ! Magnetic
              if (Bfield%Bstrength(iz).gt.TINYB) then

                ! Flag
                lfield = .True.

                ! If Line of sight
                if (LOS) then

                  ! Point to geometrical tensor
                  TKQo => Geom%TBL(:,:,:,:,iz)

                ! Quadrature
                else

                  ! Point to geometrical tensor
                  TKQo => Geom%TBo(:,:,:,:,iz)

                end if

              ! No magnetic field
              else

                ! Flag
                lfield = .False.

                ! If Line of sight
                if (LOS) then

                  ! Point to geometrical tensor
                  TKQo => Geom%TSL

                ! Quadrature
                else

                  ! Point to geometrical tensor
                  TKQo => Geom%TSo

                end if

              end if ! Magnetic field

              ! Velocity
              if (dyn) then

                ! Compute maximum velocity
                vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                           Atmo%vy(iz)*Atmo%vy(iz) + &
                           Atmo%vz(iz)*Atmo%vz(iz))

                ! Check if large enough
                lvel = vel.gt.TINYVEL

              ! No velocity
              else

                ! Trivial
                vel = 0d0
                lvel = .False.

              end if

              ! Thermal part of the Doppler width
              DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

              ! Add the microt. to Doppler width
              Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                              Atmo%vmi(iz)**2d0)

              ! Get zao index
              indx = Red%izao(jtran,ia,iz)

              ! Get limits for this transition
              if0Il2 = Red%zao(indx)%Igf0
              if1Il2 = Red%zao(indx)%Igf1

              ! Get left frequency indexes
              if (iz.eq.iz0) then
                Mif0 = p_fed%Mi0(pid) - (iz0-Rz0)*p_fed%nfreq
              else
                Mif0 = 1
              end if

              ! Get right frequency index
              if (iz.eq.iz1) then
                Mif1 = p_fed%Mi1(pid) - (iz1-Rz0)*p_fed%nfreq
              else
                Mif1 = p_fed%nfreq
              end if

              ! Initialize
              eps0 = 0d0
              eps1 = 0d0
              eps2 = 0d0
              eps3 = 0d0
              eps20 = 0d0
              eps21 = 0d0
              eps22 = 0d0
              eps23 = 0d0

              ! Emissivity for local calculation
              if ((if1Il2-if0Il2).ge.0) then

                ! Pointers to zao, pzao, and rzao
                p_red => Red%zao(indx)
                p_Norm => Red%pzao(indx)
                if (PRAM) p_rwarr => Red%rzao(indx)

                ! If magnetic
                if (lfield) then

                  ! Call first order emissivity
                  call emiss(Atom(ia),TKQo,Frec%omega,Flgsg,jtran, &
                             itermu,itermf,iz,if0Il2,if1Il2, &
                             ndir,p_Norm,Dw, &
                             eps0(:,if0Il2:if1Il2), &
                             eps1(:,if0Il2:if1Il2), &
                             eps2(:,if0Il2:if1Il2), &
                             eps3(:,if0Il2:if1Il2))

                  ! Call second order emissivity
                  call emiss2ord(Atom(ia),Geom,Atmo%vx(iz), &
                                 Atmo%vy(iz),Atmo%vz(iz),lvel, &
                                 Frec%omega, &
                                 p_fed,p_red,p_rwarr,p_Norm, &
                                 Flgsg,jtran, &
                                 itermu,itermf,iz, &
                                 if0Il2,if1Il2,Mif0,Mif1, &
                                 DwT,Dw,Bfield,Atmo%vmi(iz), &
                                 TKQo,Stokes(:,:,:,:,iz), &
                                 JKQa,JKQ(:,:,:,iz),JKQC(:,:,:,iz), &
                                 eps20(:,if0Il2:if1Il2), &
                                 eps21(:,if0Il2:if1Il2), &
                                 eps22(:,if0Il2:if1Il2), &
                                 eps23(:,if0Il2:if1Il2))

                ! If non-magnetic
                else

                  ! Call first order emissivity
                  call emissNB(Atom(ia),TKQo,Frec%omega,Flgsg,jtran, &
                               itermu,itermf,iz,if0Il2,if1Il2, &
                               ndir,p_Norm,Dw, &
                               eps0(:,if0Il2:if1Il2), &
                               eps1(:,if0Il2:if1Il2), &
                               eps2(:,if0Il2:if1Il2), &
                               eps3(:,if0Il2:if1Il2))

                  ! Call second order emissivity
                  call emiss2ordNB(Atom(ia),Geom,Atmo%vx(iz), &
                                   Atmo%vy(iz),Atmo%vz(iz),lvel, &
                                   Frec%omega, &
                                   p_fed,p_red,p_rwarr,p_Norm, &
                                   Flgsg,jtran, &
                                   itermu,itermf,iz, &
                                   if0Il2,if1Il2,Mif0,Mif1, &
                                   DwT,Dw,Atmo%vmi(iz), &
                                   TKQo,Stokes(:,:,:,:,iz), &
                                   JKQa,JKQ(:,:,:,iz), &
                                   JKQC(:,:,:,iz), &
                                   eps20(:,if0Il2:if1Il2), &
                                   eps21(:,if0Il2:if1Il2), &
                                   eps22(:,if0Il2:if1Il2), &
                                   eps23(:,if0Il2:if1Il2))

                end if ! Magnetic-field

                ! Combine to get total emissivity
                eps20(:,if0Il2:if1Il2) = eps20(:,if0Il2:if1Il2) + &
                                         eps0(:,if0Il2:if1Il2)
                eps21(:,if0Il2:if1Il2) = eps21(:,if0Il2:if1Il2) + &
                                         eps1(:,if0Il2:if1Il2)
                eps22(:,if0Il2:if1Il2) = eps22(:,if0Il2:if1Il2) + &
                                         eps2(:,if0Il2:if1Il2)
                eps23(:,if0Il2:if1Il2) = eps23(:,if0Il2:if1Il2) + &
                                         eps3(:,if0Il2:if1Il2)

                ! Second rolling index
                kk = 0

                ! Save
                do iran=1,p_fed%nran
                  do ifreq=p_fed%if0(iran),p_fed%if1(iran)

                    ! Advance index
                    kk = kk + 1

                    ! Manage MPI
                    if (kk.lt.Mif0) cycle
                    if (kk.gt.Mif1) exit

                    ! Advance index
                    jj = jj + 1

                    ! Save
                    eps20123_s(:,1,jj) = eps20(:,ifreq)
                    eps20123_s(:,2,jj) = eps21(:,ifreq)
                    eps20123_s(:,3,jj) = eps22(:,ifreq)
                    eps20123_s(:,4,jj) = eps23(:,ifreq)

                  end do
                end do

              end if ! Local calculation
            end if ! Within own height limits

            ! If MPI
            if (nproc.gt.1) then
#ifndef oldmpi
              ! If not first call
              if (.not.first) then

                ! Check finished
                call MPI_TESTSOME(size(done_index),requests, &
                                  ntest,done_index, &
                                  MPI_STATUSES_IGNORE,ierr)

                ! Update load
                W_load = W_load - ntest

              end if

              ! Not the first anymore
              first = .False.

              ! Send while we can afford
              do while (W_load.lt.W_max.and.lgiz.lt.iz)

                ! Advance
                W_load = W_load + 1
                lgiz = lgiz + 1

                ! Share
                call share_emiss(p_fed,lgiz,ll,nf,ndir*4, &
                                 nsend(lgiz),nf_s(:,lgiz), &
                                 disp(:,lgiz),requests(lgiz), &
                                 eps20123_r,eps20123_s)

                ! Check finished
                call MPI_TESTSOME(size(done_index),requests, &
                                  ntest,done_index, &
                                  MPI_STATUSES_IGNORE,ierr)

                ! Update load
                W_load = W_load - ntest

              end do ! Send jobs till load is full or we are done
#else
              call share_emiss(p_fed,iz,ll,nf,ndir*4,nsend(1), &
                               nf_s(:,1),disp(:,1), &
                               eps20123_r,eps20123_s)
#endif
            end if ! MPI

          end do ! Heights

          ! Free
          deallocate(eps0,eps1,eps2,eps3)
          deallocate(eps20,eps21,eps22,eps23)
          nullify(eps20,eps21,eps22,eps23)

          ! If serial, just point
          if (nproc.eq.1) eps20123_r => eps20123_s

          !
          ! Get what is needed
          !

          ! For each height
          do iz=Rz0,Rz1_PRD

            ! Get zao index
            indx = Red%izao(jtran,ia,iz)

#ifndef oldmpi
            ! If we are not done with sharing
            if (lgiz.lt.Rz1_PRD) then

              ! Check finished
              call MPI_TESTSOME(size(done_index),requests, &
                                ntest,done_index, &
                                MPI_STATUSES_IGNORE,ierr)

              ! Update load
              W_load = W_load - ntest

              ! Send while we can afford
              do while (W_load.lt.W_max.and.lgiz.lt.Rz1_PRD)

                ! Advance
                W_load = W_load + 1
                lgiz = lgiz + 1

                ! Share
                call share_emiss(p_fed,lgiz,ll,nf,ndir*4, &
                                 nsend(lgiz),nf_s(:,lgiz), &
                                 disp(:,lgiz),requests(lgiz), &
                                 eps20123_r,eps20123_s)

              end do ! Send jobs till load is full or we are done

              ! If we actually need this data now
              do while (lgiz.lt.iz)

                ! Check finished
                call MPI_TESTSOME(size(done_index),requests, &
                                  ntest,done_index, &
                                  MPI_STATUSES_IGNORE,ierr)

                ! Update load
                W_load = W_load - ntest

                ! Send while we can afford
                do while (W_load.lt.W_max.and.lgiz.lt.Rz1_PRD)

                  ! Advance
                  W_load = W_load + 1
                  lgiz = lgiz + 1

                  ! Share
                  call share_emiss(p_fed,lgiz,ll,nf,ndir*4, &
                                   nsend(lgiz),nf_s(:,lgiz), &
                                   disp(:,lgiz),requests(lgiz), &
                                   eps20123_r,eps20123_s)

                end do ! Send jobs till load is full or we are done

              end do ! We need this data, now

            end if ! While we still need to communicate
#endif
            ! If the CPU allocated eps20
            if (allocated(Red%zao(indx)%eps20)) then

              ! Magnetic
              if (Bfield%Bstrength(iz).gt.TINYB) then

                ! Flag
                lfield = .True.

                ! If Line of sight
                if (LOS) then

                  ! Point to geometrical tensor
                  TKQo => Geom%TBL(:,:,:,:,iz)

                ! Quadrature
                else

                  ! Point to geometrical tensor
                  TKQo => Geom%TBo(:,:,:,:,iz)

                end if

              ! No magnetic field
              else

                ! Flag
                lfield = .False.

                ! If Line of sight
                if (LOS) then

                  ! Point to geometrical tensor
                  TKQo => Geom%TSL

                ! Quadrature
                else

                  ! Point to geometrical tensor
                  TKQo => Geom%TSo

                end if

              end if ! Magnetic field

              ! Velocity
              if (dyn) then

                ! Compute maximum velocity
                vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                           Atmo%vy(iz)*Atmo%vy(iz) + &
                           Atmo%vz(iz)*Atmo%vz(iz))

                ! Check if large enough
                lvel = vel.gt.TINYVEL

              ! No velocity
              else

                ! Trivial
                vel = 0d0
                lvel = .False.

              end if ! Velocity

              ! Indexes to point to
              i0 = 1 + (iz-Rz0)*nf
              i1 = i0 + nf - 1
#ifndef oldmpi
              ! If MPI and request was not completed
              if (nproc.gt.1.and. &
                  requests(iz).ne.MPI_REQUEST_NULL) then

                ! Check request
                call MPI_WAIT(requests(iz),MPI_STATUS_IGNORE,ierr)
                W_load = W_load - 1

              end if
#endif
              ! Point to relevant data
              eps20 => eps20123_r(:,1,i0:i1)
              eps21 => eps20123_r(:,2,i0:i1)
              eps22 => eps20123_r(:,3,i0:i1)
              eps23 => eps20123_r(:,4,i0:i1)

              ! If dynamic, we need to interpolate
              if (lvel) then

                ! Number of frequencies
                nfl = if1l2 - if0l2 + 1

                ! Get memory space
                allocate(omega(if0l2:if1l2))
                if (int_mode.eq.1) allocate(splin(nf,3))

                ! Initialize direction
                idir = 0

                ! For every polar direction
                do ith=1,nTh

                  ! LOS
                  if (los) then

                    ! Get cos polar
                    ct = Geom%L_mu(lth)

                  ! Quadrature
                  else

                    ! Get cos polar
                    ct = Geom%V_mu(ith)

                  end if ! LOS/quadrature

                  ! Get sin polar
                  st = sqrt(1d0 - ct*ct)

                  ! For every azimuth
                  do iph=1,nPh

                    ! LOS
                    if (los) then

                      ! Get cos and sin azimuth
                      cc = cos(Geom%L_phi(lph))
                      sc = sin(Geom%L_phi(lph))

                    ! Quadrature
                    else

                      ! Get cos and sin azimuth
                      cc = Geom%v_mux(iph)
                      sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

                    end if ! LOS/quadrature

                    ! Advance index
                    idir = idir + 1

                    ! Velocity factor
                    vfac = 1d0 - atmo%vx(iz)*st*cc - &
                                 atmo%vy(iz)*st*sc - &
                                 atmo%vz(iz)*ct

                    ! Get new omega
                    omega = Frec%omega(if0l2:if1l2)*vfac

                    !
                    ! Splines
                    !
                    if (int_mode.eq.1) then

                      ! Interpolate with splines
                      call spline(Frec%omega(if0tl2:if1tl2), &
                                  eps20(idir,:), &
                                  splin(:,1), &
                                  splin(:,2), &
                                  splin(:,3), nf)
                      do ifreq=if0l2,if1l2
                        Red%zao(indx)%eps20(idir,ifreq) = &
                                  ispline(omega(ifreq), &
                                          Frec%omega(if0tl2:if1tl2), &
                                          eps20(idir,:), &
                                          splin(:,1),splin(:,2), &
                                          splin(:,3),nf)
                      end do

                      ! Interpolate with splines
                      call spline(Frec%omega(if0tl2:if1tl2), &
                                  eps21(idir,:), &
                                  splin(:,1), &
                                  splin(:,2), &
                                  splin(:,3), nf)
                      do ifreq=if0l2,if1l2
                        Red%zao(indx)%eps21(idir,ifreq) = &
                                  ispline(omega(ifreq), &
                                          Frec%omega(if0tl2:if1tl2), &
                                          eps21(idir,:), &
                                          splin(:,1),splin(:,2), &
                                          splin(:,3),nf)
                      end do

                      ! Interpolate with splines
                      call spline(Frec%omega(if0tl2:if1tl2), &
                                  eps22(idir,:), &
                                  splin(:,1), &
                                  splin(:,2), &
                                  splin(:,3), nf)
                      do ifreq=if0l2,if1l2
                        Red%zao(indx)%eps22(idir,ifreq) = &
                                  ispline(omega(ifreq), &
                                          Frec%omega(if0tl2:if1tl2), &
                                          eps22(idir,:), &
                                          splin(:,1),splin(:,2), &
                                          splin(:,3),nf)
                      end do

                      ! Interpolate with splines
                      call spline(Frec%omega(if0tl2:if1tl2), &
                                  eps23(idir,:), &
                                  splin(:,1), &
                                  splin(:,2), &
                                  splin(:,3), nf)
                      do ifreq=if0l2,if1l2
                        Red%zao(indx)%eps23(idir,ifreq) = &
                                  ispline(omega(ifreq), &
                                          Frec%omega(if0tl2:if1tl2), &
                                          eps23(idir,:), &
                                          splin(:,1),splin(:,2), &
                                          splin(:,3),nf)
                      end do

                    !
                    ! Linear interpolation
                    !
                    else if (int_mode.eq.0) then

                      ! Interpolate linear
                      call Intpol_Lin(Frec%omega(if0tl2:if1tl2), &
                                      eps20(idir,:),nf, &
                                      omega(if0l2:if1l2), &
                                      Red%zao(indx)% &
                                          eps20(idir,if0l2:if1l2), &
                                      nfl)

                      ! Interpolate linear
                      call Intpol_Lin(Frec%omega(if0tl2:if1tl2), &
                                      eps21(idir,:),nf, &
                                      omega(if0l2:if1l2), &
                                      Red%zao(indx)% &
                                          eps21(idir,if0l2:if1l2), &
                                      nfl)

                      ! Interpolate linear
                      call Intpol_Lin(Frec%omega(if0tl2:if1tl2), &
                                      eps22(idir,:),nf, &
                                      omega(if0l2:if1l2), &
                                      Red%zao(indx)% &
                                          eps22(idir,if0l2:if1l2), &
                                      nfl)

                      ! Interpolate linear
                      call Intpol_Lin(Frec%omega(if0tl2:if1tl2), &
                                      eps22(idir,:),nf, &
                                      omega(if0l2:if1l2), &
                                      Red%zao(indx)% &
                                          eps22(idir,if0l2:if1l2), &
                                      nfl)

                    end if ! Type of interpolation

                  end do ! Azimuth angle
                end do ! Polar angle

                ! Free 
                deallocate(omega)
                if (int_mode.eq.1) deallocate(splin)

              ! If static
              else

                ! Limits
                i0 = if0l2 - if0tl2 + 1
                i1 = if1l2 - if0tl2 + 1

                ! Just fetch
                Red%zao(indx)%eps20 = eps20(:,i0:i1)
                Red%zao(indx)%eps21 = eps21(:,i0:i1)
                Red%zao(indx)%eps22 = eps22(:,i0:i1)
                Red%zao(indx)%eps23 = eps23(:,i0:i1)

              end if ! Dynamic
#ifndef oldmpi
            ! No data but MPI and request not completed
            else if (nproc.gt.1.and. &
                     requests(iz).ne.MPI_REQUEST_NULL) then

              ! Check reception
              call MPI_WAIT(requests(iz),MPI_STATUS_IGNORE,ierr)
              W_load = W_load - 1
#endif
            end if ! If need to take anything

          end do ! Heights

          ! Free
          nullify(eps20,eps21,eps22,eps23)
          if (nproc.gt.1) deallocate(eps20123_r)
          nullify(eps20123_r)
          deallocate(eps20123_s)

          ! Error
          if (laborted) exit

        end do ! Output transition

        ! Error
        if (laborted) exit

      end do ! Atoms

      ! Nullify pointers
      nullify(p_fed,p_red,p_rwarr,p_Norm,eps20123_r,TKQo)

      return

      end subroutine comoving_emiss2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute second order emissivity for intensity\n
      !!      Atom(Atom_class(:)): Structures with atomic data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!           Red(Red_class): Structure with redistribution
      !!                           input frequency data,
      !!                           redistribution function data, and
      !!                           profile or normalization data\n
      !!  Stokes(double(:,:,:,:)): Stokes parameters\n
      !!         JKQ(double(:,:)): Mean intensity integrated over the
      !!                           absorption profile\n
      !!        JKQC(double(:,:)): Mean intensity with frequency
      !!                           dependence
      !!        int_mode(integer): Interpolation mode\n
      !!             lth(integer): Current index for polar direction
      !!                           if doing los formal solution\n
      !!             lph(integer): Current index for azimuth direction
      !!                           if doing los formal solution\n
      !!             los(logical): If performing last formal
      !!                           solution
      subroutine comoving_emissI2ord(Atom,Atmo,Geom,Frec,Red, &
                                     Stokes,JKQ,JKQC,lth,lph, &
                                     int_mode,los)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Geometry_class), intent(in):: Geom
      logical, intent(in):: los
      integer, intent(in):: lth,lph,int_mode
      double precision, dimension(nfreq,Geom%nPh, &
                                  Geom%nTh,giz0:giz1), &
                        intent(in):: Stokes
      double precision, dimension(nxt,Rz0:Rz1), intent(in):: JKQ
      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: JKQC

      ! Local

#ifndef oldmpi
      type(MPI_request), dimension(:), allocatable:: requests
#endif

      logical:: lvel,first

      integer:: iz0,iz1,Mif0,Mif1,jndx,jj,kk,ll,iran
      integer:: iz,ia,jtran,itermu,itermf,indx,jdir,t0,t1,i0,i1
      integer:: if0l2,if1l2,if0tl2,if1tl2,if0Il2,if1Il2,ifreq,nf,nfl
      integer:: if0jl2,if1jl2,ffjtran,ffktran,ffltran,fjtran,lnz
      integer:: iJf,iJu,ndir,njdir,idir,ierr,iph,ith,nth,nph
#ifndef oldmpi
      integer:: ntest,W_max,W_load,lgiz
      integer, dimension(:), allocatable:: done_index
#endif
      integer, dimension(:), allocatable:: nsend
      integer, dimension(:,:), allocatable:: nf_s,disp

      double precision:: vel,DwT,Dw,ct,st,cc,sc,vfac,Dfreq
      double precision, dimension(:), allocatable:: omega
      double precision, dimension(:), allocatable:: eps0
      double precision, dimension(:,:), allocatable:: splin
      double precision, dimension(:,:,:), &
                        allocatable, target:: eps20rpf_s

      ! Pointers
      type(Reda_class), pointer:: p_fed
      type(Redb_class), pointer:: p_red
      type(Redb2_class), pointer:: p_rwarr
      type(Prof_class), pointer:: p_Norm
      double precision, dimension(:,:), pointer:: eps20,rpf
      double precision, dimension(:,:,:), pointer:: eps20rpf_r
      type(Redb2_class), target:: p_dummy


      ! Nullify pointers
      nullify(p_fed,p_red,p_rwarr,p_Norm,eps20rpf_r)

      ! If MPI
      if (nproc.gt.1) then
#ifndef oldmpi
        ! Initialize requests
        allocate(requests(Rz0:Rz1_PRD))
        allocate(done_index(Rz0:Rz1_PRD))
        ntest = 0
        do iz=Rz0,Rz1_PRD
          requests(iz) = MPI_REQUEST_NULL
          done_index(iz) = 0
        end do
        ! Allocate frequency size for output
        allocate(nsend(Rz0:Rz1_PRD))
        allocate(nf_s(0:nproc-1,Rz0:Rz1_PRD))
        allocate(disp(0:nproc-1,Rz0:Rz1_PRD))
#else
        ! Allocate frequency size for output
        allocate(nsend(1)
        allocate(nf_s(0:nproc-1,1))
        allocate(disp(0:nproc-1,1))
#endif
      end if ! MPI

      ! If no atoms, leave
      if (nA.lt.1) return

      ! Dummy
      if (.not.IRAM) p_rwarr => p_dummy

      ! Number of heights
      lnz = Rz1_PRD - Rz0 + 1

      ! If dynamic
      if (dyn) then

        ! Number of output directions in geometry
        njdir = Geom%njdir

      ! No velocity
      else

        ! If angle-averaged
        if (AVI) then

          ! Isotropic
          njdir = 1

        ! If angle-dependent
        else

          ! In general non-isotropic
          njdir = Geom%njdir

        end if ! AA/AD
      end if ! Dynamic

      ! Number of directions
      if (los) then
        nth = 1
        nph = 1
      else
        nth = Geom%nTh
        nph = Geom%nPh
      end if

      ! For each atom
      do ia=1,nA

        ! For each output transition
        do jtran=1,Atom(ia)%ntran

          ! Skip if no PRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          itermf = Atom(ia)%fst(jtran)%iterml

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Idenfity involved levels
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJf = Atom(ia)%fst(jtran)%ilevell(fjtran)

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
            ffktran = ffjtran + Atom(ia)%tfshift

            ! Get frequency of FS transition
            Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                    Atom(ia)%FSfreq(iJf,itermf)

            ! Real index in trano
            ffltran = Atom(ia)%itrano(ffjtran)

            ! Get izao index
            jndx = Red%izao(ffjtran,ia,Rz0)

            ! Point to relevant variable
            p_fed => Red%ao(jndx)

            ! Get limits for this transition
            if0l2 = p_fed%gf0
            if1l2 = p_fed%gf1
            if0tl2 = p_fed%tgf0
            if1tl2 = p_fed%tgf1

            ! If no output, skip
            if ((if1tl2-if0tl2).lt.0) cycle

            ! Initialize first communication
            first = .True.

            ! Allocate emissivity to collect
            allocate(eps0(if0tl2:if1tl2))

            ! Number of frequencies
            nf = if1tl2 - if0tl2 + 1

            ! If doing MPI
            if (nproc.gt.1) then

              ! Allocate emissivity to receive
              allocate(eps20rpf_r(njdir,2,nf*lnz))
#ifndef oldmpi
              ! Initialize requests
              ntest = 0
              do iz=Rz0,Rz1_PRD
                requests(iz) = MPI_REQUEST_NULL
                done_index(iz) = 0
              end do

              ! Calculate maximum working load
              W_max = min(W_nominal,max(W_minimal, &
                                      nint(mem_limit/ &
                          dble(nf*sizeof(eps20rpf_r(:,:,1))))))

              ! Initialize last height index sent and load
              lgiz = Rz0 - 1
              W_load = 0

              ! Initialize MPI data
              nsend = 0
              nf_s = 0
              disp = 0
#endif
            end if ! MPI

            ! Has work to do
            if (p_fed%nn(pid).gt.0) then

              ! Get height limits
              iz0 = Rz0 + (p_fed%Mi0(pid)-1)/p_fed%nfreq
              iz1 = Rz0 + (p_fed%Mi1(pid)-1)/p_fed%nfreq

              ! Allocate sender
              allocate(eps20rpf_s(njdir,2,p_fed%nn(pid)))

              ! Initialize rolling index
              jj = 0

            ! Does not have work
            else

              ! Ensure no looping
              iz0 = 0
              iz1 = -1

            end if ! Sanity check

            ! Initialize index
            ll = 0

            ! For each height
            do iz=Rz0,Rz1_PRD

              ! If within limits
              if (iz.ge.iz0.and.iz.le.iz1) then

                ! Velocity
                if (dyn) then

                  ! Compute maximum velocity
                  vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                             Atmo%vy(iz)*Atmo%vy(iz) + &
                             Atmo%vz(iz)*Atmo%vz(iz))

                  ! Check velocity magnitude
                  lvel = vel.gt.TINYVEL

                  ! If large enough velocity
                  if (lvel) then

                    ! Directions from geometry structure
                    ndir = Geom%njdir

                  ! Small velocity
                  else

                    ! If angle-averaged
                    if (AVI) then

                      ! Isotropic
                      ndir = 1

                    ! If angle-dependent
                    else

                      ! In general non-isotropic
                      ndir = Geom%njdir

                    end if ! AA/AD
                  end if ! Large velocity

                ! No velocity
                else

                  ! Trivial
                  vel = 0d0
                  lvel = .False.

                  ! If angle-averaged
                  if (AVI) then

                    ! Isotropic
                    ndir = 1

                  ! If angle-dependent
                  else

                    ! In general non-isotropic
                    ndir = Geom%njdir

                  end if ! AA/AD
                end if ! Dynamic

                ! Thermal part of the Doppler width
                DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

                ! Add the microt. to Doppler width
                Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

                ! Allocate emissivity to collect
                allocate(eps20(ndir,if0tl2:if1tl2))
                allocate(rpf(ndir,if0tl2:if1tl2))

                ! Get izao index
                indx = Red%izao(ffjtran,ia,iz)

                ! Get limits for this transition
                if0Il2 = Red%zao(indx)%Igf0
                if1Il2 = Red%zao(indx)%Igf1
                if0jl2 = Red%zao(indx)%ggf0
                if1jl2 = Red%zao(indx)%ggf1

                ! Get left frequency indexes
                if (iz.eq.iz0) then
                  Mif0 = p_fed%Mi0(pid) - (iz0-Rz0)*p_fed%nfreq
                else
                  Mif0 = 1
                end if

                ! Get right frequency index
                if (iz.eq.iz1) then
                  Mif1 = p_fed%Mi1(pid) - (iz1-Rz0)*p_fed%nfreq
                else
                  Mif1 = p_fed%nfreq
                end if

                ! Transition ranges
                t0 = Atom(ia)%tfshift + 1
                t1 = t0 + Atom(ia)%nftran - 1

                ! Initialize
                eps0 = 0d0
                eps20 = 0d0
                rpf = 0d0

                ! Emissivity for local calculation
                if ((if1Il2-if0Il2).ge.0) then

                  ! Pointers to zao, pzao, and rzao
                  p_red => Red%zao(indx)
                  p_Norm => Red%pzao(indx)
                  if (IRAM) p_rwarr => Red%rzao(indx)

                  ! First order
                  call emissI(Atom(ia),Frec%omega,jtran, &
                              itermu,itermf,iJu,iJf,iz, &
                              if0Il2,if1Il2, &
                              p_Norm,Dw,eps0(if0Il2:if1Il2))

                  ! Second order
                  call emissI2ord(Atom(ia),Geom,Atmo%vx(iz), &
                                  Atmo%vy(iz),Atmo%vz(iz),lvel, &
                                  Frec%omega, &
                                  p_fed,p_red,p_rwarr,p_Norm, &
                                  ndir,jtran,fjtran,itermu,itermf, &
                                  iJu,iJf,iz, &
                                  if0Il2,if1Il2,Mif0,Mif1,DwT,Dw, &
                                  Atmo%vmi(iz), &
                                  Stokes(:,:,:,iz), &
                                  JKQ(t0:t1,iz), &
                                  JKQC(if0jl2:if1jl2,iz), &
                                  eps20(:,if0Il2:if1Il2), &
                                  rpf(:,if0Il2:if1Il2))

                  ! Combine
                  do idir=1,ndir
                    eps20(idir,if0Il2:if1Il2) = &
                                             eps0(if0Il2:if1Il2) + &
                                             eps20(idir,if0Il2:if1Il2)
                  end do

                  ! Second rolling index
                  kk = 0

                  ! Save
                  do iran=1,p_fed%nran
                    do ifreq=p_fed%if0(iran),p_fed%if1(iran)

                      ! Advance index
                      kk = kk + 1

                      ! Manage MPI
                      if (kk.lt.Mif0) cycle
                      if (kk.gt.Mif1) exit

                      ! Advance index
                      jj = jj + 1

                      ! Save
                      eps20rpf_s(1:ndir,1,jj) = eps20(:,ifreq)
                      eps20rpf_s(1:ndir,2,jj) = rpf(:,ifreq)

                    end do
                  end do

                end if ! Local calculation

                ! Free
                deallocate(eps20,rpf)
                nullify(eps20,rpf)

              end if ! Within own height limits

              ! If MPI
              if (nproc.gt.1) then
#ifndef oldmpi
                ! If not first call
                if (.not.first) then

                  ! Check finished
                  call MPI_TESTSOME(size(done_index),requests, &
                                    ntest,done_index, &
                                    MPI_STATUSES_IGNORE,ierr)

                  ! Update load
                  W_load = W_load - ntest

                end if

                ! Not the first anymore
                first = .False.

                ! Send while we can afford
                do while (W_load.lt.W_max.and.lgiz.lt.iz)

                  ! Advance
                  W_load = W_load + 1
                  lgiz = lgiz + 1

                  ! Share
                  call share_emiss(p_fed,lgiz,ll,nf,njdir*2, &
                                   nsend(lgiz),nf_s(:,lgiz), &
                                   disp(:,lgiz),requests(lgiz), &
                                   eps20rpf_r,eps20rpf_s)

                  ! Check finished
                  call MPI_TESTSOME(size(done_index),requests, &
                                    ntest,done_index, &
                                    MPI_STATUSES_IGNORE,ierr)

                  ! Update load
                  W_load = W_load - ntest

                end do ! Send jobs till load is full or we are done
#else
                call share_emiss(p_fed,iz,ll,nf,njdir*2,nsend(1), &
                                 nf_s(:,1),disp(:,1), &
                                 eps20rpf_r,eps20rpf_s)
#endif
              end if ! MPI

            end do ! Heights

            ! Free
            deallocate(eps0)

            ! If serial, just point
            if (nproc.eq.1) eps20rpf_r => eps20rpf_s


            !
            ! Get what is needed
            !

            ! For each height
            do iz=Rz0,Rz1_PRD

              ! Get index
              indx = Red%izao(ffjtran,ia,iz)
#ifndef oldmpi
              ! If we are not done with sharing
              if (lgiz.lt.Rz1_PRD) then

                ! Check finished
                call MPI_TESTSOME(size(done_index),requests, &
                                  ntest,done_index, &
                                  MPI_STATUSES_IGNORE,ierr)

                ! Update load
                W_load = W_load - ntest

                ! Send while we can afford
                do while (W_load.lt.W_max.and.lgiz.lt.Rz1_PRD)

                  ! Advance
                  W_load = W_load + 1
                  lgiz = lgiz + 1

                  ! Share
                  call share_emiss(p_fed,lgiz,ll,nf,njdir*2, &
                                   nsend(lgiz),nf_s(:,lgiz), &
                                   disp(:,lgiz),requests(lgiz), &
                                   eps20rpf_r,eps20rpf_s)

                end do ! Send jobs till load is full or we are done

                ! If we actually need this data now
                do while (lgiz.lt.iz)

                  ! Check finished
                  call MPI_TESTSOME(size(done_index),requests, &
                                    ntest,done_index, &
                                    MPI_STATUSES_IGNORE,ierr)

                  ! Update load
                  W_load = W_load - ntest

                  ! Send while we can afford
                  do while (W_load.lt.W_max.and.lgiz.lt.Rz1_PRD)

                    ! Advance
                    W_load = W_load + 1
                    lgiz = lgiz + 1

                    ! Share
                    call share_emiss(p_fed,lgiz,ll,nf,njdir*2, &
                                     nsend(lgiz),nf_s(:,lgiz), &
                                     disp(:,lgiz),requests(lgiz), &
                                     eps20rpf_r,eps20rpf_s)

                  end do ! Send jobs till load is full or we are done

                end do ! We need this data, now

              end if ! While we still need to communicate
#endif
              ! If emissivity is allocated
              if (allocated(Red%zao(indx)%eps20)) then

                ! Velocity
                if (dyn) then

                  ! Compute maximum velocity
                  vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                             Atmo%vy(iz)*Atmo%vy(iz) + &
                             Atmo%vz(iz)*Atmo%vz(iz))

                  ! Check velocity magnitude
                  lvel = vel.gt.TINYVEL

                  ! If large enough velocity
                  if (lvel) then

                    ! Directions from geometry structure
                    ndir = Geom%njdir

                  ! Small velocity
                  else

                    ! If angle-averaged
                    if (AVI) then

                      ! Isotropic
                      ndir = 1

                    ! If angle-dependent
                    else

                      ! In general non-isotropic
                      ndir = Geom%njdir

                    end if ! AA/AD
                  end if ! Large velocity

                ! No velocity
                else

                  ! Trivial
                  vel = 0d0
                  lvel = .False.

                  ! If angle-averaged
                  if (AVI) then

                    ! Isotropic
                    ndir = 1

                  ! If angle-dependent
                  else

                    ! In general non-isotropic
                    ndir = Geom%njdir

                  end if ! AA/AD
                end if ! Dynamic

                ! Indentify pointing indexes
                i0 = 1 + (iz-Rz0)*nf
                i1 = i0 + nf - 1
#ifndef oldmpi
                ! If MPI and request not completed yet
                if (nproc.gt.1.and. &
                    requests(iz).ne.MPI_REQUEST_NULL) then

                  ! Check request
                  call MPI_WAIT(requests(iz),MPI_STATUS_IGNORE,ierr)
                  W_load = W_load - 1

                end if
#endif
                ! Point to relevant data
                eps20 => eps20rpf_r(:,1,i0:i1)
                rpf => eps20rpf_r(:,2,i0:i1)

                ! If dynamic, we need to interpolate
                if (lvel) then

                  ! Number of frequencies
                  nfl = if1l2 - if0l2 + 1

                  ! Get memory space
                  allocate(omega(if0l2:if1l2))
                  if (int_mode.eq.1) allocate(splin(nf,3))

                  ! Initialize direction
                  idir = 0

                  ! For every polar direction
                  do ith=1,nTh

                    ! LOS
                    if (los) then

                      ! Get cos polar
                      ct = Geom%L_mu(lth)

                    ! Quadrature
                    else

                      ! Get cos polar
                      ct = Geom%V_mu(ith)

                    end if ! LOS/quadrature

                    ! Get sin polar
                    st = sqrt(1d0 - ct*ct)

                    ! For every azimuth
                    do iph=1,nPh

                      ! LOS
                      if (los) then

                        ! Get cos and sin azimuth
                        cc = cos(Geom%L_phi(lph))
                        sc = sin(Geom%L_phi(lph))

                      ! Quadrature
                      else

                        ! Get cos and sin azimuth
                        cc = Geom%v_mux(iph)
                        sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

                      end if ! LOS/quadrature

                      ! Advance index
                      idir = idir + 1

                      ! If angle-averaged
                      if (AVI) then

                        ! Corrected index
                        jdir = 1

                      ! If angle-dependent
                      else

                        ! Copy index
                        jdir = idir

                      end if ! AA/AD

                      ! Velocity factor
                      vfac = 1d0 - atmo%vx(iz)*st*cc - &
                                   atmo%vy(iz)*st*sc - &
                                   atmo%vz(iz)*ct

                      ! Get new omega
                      omega = Frec%omega(if0l2:if1l2)*vfac

                      !
                      ! Splines
                      !
                      if (int_mode.eq.1) then

                        ! Interpolate with splines
                        call spline(Frec%omega(if0tl2:if1tl2), &
                                    eps20(jdir,:), &
                                    splin(:,1), &
                                    splin(:,2), &
                                    splin(:,3), nf)
                        do ifreq=if0l2,if1l2
                          Red%zao(indx)%eps20(idir,ifreq) = &
                                  ispline(omega(ifreq), &
                                          Frec%omega(if0tl2:if1tl2), &
                                          eps20(jdir,:), &
                                          splin(:,1),splin(:,2), &
                                          splin(:,3),nf)
                        end do

                        ! Interpolate with splines
                        call spline(Frec%omega(if0tl2:if1tl2), &
                                    rpf(jdir,:), &
                                    splin(:,1), &
                                    splin(:,2), &
                                    splin(:,3), nf)
                        do ifreq=if0l2,if1l2
                          Red%zao(indx)%rpf(idir,ifreq) = &
                                  ispline(omega(ifreq), &
                                          Frec%omega(if0tl2:if1tl2), &
                                          rpf(jdir,:), &
                                          splin(:,1),splin(:,2), &
                                          splin(:,3),nf)
                        end do

                      !
                      ! Linear
                      !
                      else if (int_mode.eq.0) then

                        ! Interpolate linear
                        call Intpol_Lin(Frec%omega(if0tl2:if1tl2), &
                                        eps20(jdir,:),nf, &
                                        omega(if0l2:if1l2), &
                                        Red%zao(indx)% &
                                            eps20(idir,if0l2:if1l2), &
                                        nfl)

                        ! Interpolate linear
                        call Intpol_Lin(Frec%omega(if0tl2:if1tl2), &
                                        rpf(jdir,:),nf, &
                                        omega(if0l2:if1l2), &
                                        Red%zao(indx)% &
                                            rpf(idir,if0l2:if1l2), &
                                        nfl)

                      end if ! Type of interpolation

                    end do ! Azimuth
                  end do ! Polar

                  ! Free
                  deallocate(omega)
                  if (int_mode.eq.1) deallocate(splin)

                ! If static
                else

                  ! Limits
                  i0 = if0l2 - if0tl2 + 1
                  i1 = if1l2 - if0tl2 + 1

                  ! Just fetch
                  Red%zao(indx)%eps20 = eps20(1:ndir,i0:i1)
                  Red%zao(indx)%rpf = rpf(1:ndir,i0:i1)

                end if ! Dynamic
#ifndef oldmpi
              ! No data but MPI and request not completed
              else if (nproc.gt.1.and. &
                       requests(iz).ne.MPI_REQUEST_NULL) then

                ! Check reception
                call MPI_WAIT(requests(iz),MPI_STATUS_IGNORE,ierr)
                W_load = W_load - 1
#endif
              end if ! If need to take anything

            end do ! Heights

            ! Free
            nullify(eps20,rpf)
            if (nproc.gt.1) deallocate(eps20rpf_r)
            nullify(eps20rpf_r)
            deallocate(eps20rpf_s)

          end do ! Output transition FS
        end do ! Output transition term
      end do ! Atoms

      ! Nullify pointers
      nullify(p_fed,p_red,p_rwarr,p_Norm,eps20rpf_r)

      end subroutine comoving_emissI2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Share emissivity data in one height of the model atmosphere\n
      !!       p_fed(Reda_class): Structure with redistribution output
      !!                          frequency data\n
      !!             iz(integer): Height index to share\n
      !!             ll(integer): Counter for contiguous array\n
      !!             nf(integer): Number of frequencies for CPU\n
      !!             nn(integer): Direction and Stokes factor\n
      !!          nsend(integer): Amount of data to send\n
      !!        nf_s(integer(:)): Amount of data to receive\n
      !!        disp(integer(:)): Displacement in receiver buffer\n
      !!   request(MPI_requests): Request for non-blocking
      !!                          communication\n
      !!  datum_r(double(:,:,:)): Receiver buffer\n
      !!  datum_s(double(:,:,:)): Sender buffer
#ifndef oldmpi
      subroutine share_emiss(p_fed,iz,ll,nf,nn,nsend,nf_s,disp, &
                             request,datum_r,datum_s)
#else
      subroutine share_emiss(p_fed,iz,ll,nf,nn,nsend,nf_s,disp, &
                             datum_r,datum_s)
#endif

      ! I/O

      type(Reda_class), intent(in):: p_fed
#ifndef oldmpi
      type(MPI_request), intent(inout):: request
#endif
      integer, intent(in):: iz,nf,nn
      integer, intent(inout):: ll,nsend
      integer, dimension(0:nproc-1), intent(inout):: nf_s,disp
      double precision, dimension(:,:,:), &
                        allocatable, intent(in):: datum_s
      double precision, dimension(:,:,:), &
                        pointer, intent(inout):: datum_r

      ! Local

      integer:: i0,iproc,jz0,jz1,Mif0,Mif1,ierr


      ! Initialize
      nsend = 0
      nf_s = 0
      disp = 0
      i0 = 1

      ! For each CPU
      do iproc=0,nproc-1

        ! Get height limits
        jz0 = Rz0 + (p_fed%Mi0(iproc)-1)/p_fed%nfreq
        jz1 = Rz0 + (p_fed%Mi1(iproc)-1)/p_fed%nfreq

        ! In range
        if (iz.ge.jz0.and.iz.le.jz1) then

          ! Get left frequency indexes
          if (iz.eq.jz0) then
            Mif0 = p_fed%Mi0(iproc) - &
                   (jz0-Rz0)*p_fed%nfreq
          else
            Mif0 = 1
          end if

          ! Get right frequency index
          if (iz.eq.jz1) then
            Mif1 = p_fed%Mi1(iproc) - &
                   (jz1-Rz0)*p_fed%nfreq
          else
            Mif1 = p_fed%nfreq
          end if

          ! Save receiving
          nf_s(iproc) = Mif1 - Mif0 + 1
          disp(iproc) = Mif0 - 1

          ! If this is the CPU
          if (pid.eq.iproc) then

            ! Save number of sends
            nsend = nf_s(pid)

            ! Advance index
            ll = ll + 1

            ! Save initial point
            i0 = ll

            ! Advance rest
            ll = ll + Mif1 - Mif0

          end if ! This is the CPU
        end if ! In height range

      end do ! CPUs

      ! Scale
      nsend = nsend*nn
      nf_s = nf_s*nn
      disp = disp*nn

#ifdef oldmpi
      ! Share
      call MPI_ALLGATHERV(datum_s(1,1,i0),nsend, &
                          MPI_DOUBLE_PRECISION, &
                          datum_r(1,1,1+(iz-Rz0)*nf), &
                          nf_s,disp, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_COMM_RT,ierr)
#else
      ! Share
      call MPI_IALLGATHERV(datum_s(1,1,i0),nsend, &
                           MPI_DOUBLE_PRECISION, &
                           datum_r(1,1,1+(iz-Rz0)*nf), &
                           nf_s,disp, &
                           MPI_DOUBLE_PRECISION, &
                           MPI_COMM_RT,request,ierr)
#endif

      end subroutine share_emiss

!#####################################################################
!#####################################################################
!#####################################################################

      end module comovingprd_mod
