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
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - First version (TdPA)
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

      integer:: ia,jtran,iz,indx,if0l2,if1l2,ndir,nf


      ! Master is not going to use these, so leave if master in MPI
      if (pid.eq.0.and.nproc.gt.1) return

      ! If no PRD, how did you even get here
      if (.not.PRD) return

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

            ! Get zao index
            indx = Red%izao(jtran,ia,iz)

            ! Transition ranges
            if0l2 = Red%zao(indx)%tgf0
            if1l2 = Red%zao(indx)%tgf1
            nf = if1l2 - if0l2 + 1

            ! Size for collection
            ! The 4 assumes that will be using splines to
            ! interpolate
            if (nf.gt.0) TRAMc = TRAMc + 8d-6*dble(nf*(ndir*8 + 4))

            ! If this CPU does not have frequencies
            ! in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Get indexes for frequency
            if0l2 = Red%zao(indx)%gf0
            if1l2 = Red%zao(indx)%gf1
            nf = if1l2 - if0l2 + 1

            ! Skip if transition is outside limits
            if (nf.le.0) cycle

            ! Size for storage
            TRAMc = TRAMc + 8d-6*dble(ndir*4*nf)

          end do ! Transition
        end do ! Atom
      end do ! Height

      end subroutine predict_emiss

!#####################################################################
!#####################################################################
!#####################################################################

      !> Predict the amount of memory space that will be needed by the
      !! second order emissivity for intensity\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data
      subroutine predict_emissI(Atom,Atmo,Geom,Red)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Red_class), intent(in):: Red

      ! Local

      logical:: lvel

      integer:: ia,jtran,iz,indx,fjtran,ffjtran,if0l2,if1l2,ndir,nf

      double precision:: vel


      ! Master is not going to use these, so leave if master MPI
      if (pid.eq.0.and.nproc.gt.1) return

      ! If no PRD, how did you even get here
      if (.not.PRD) return

      ! For each height
      do iz=Rz0,Rz1_PRD

        ! If dynamic
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

            ! For each FS transition
            do fjtran=1,Atom(ia)%fst(jtran)%nt

              ! Get the sequential index of this FS transition
              ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

              ! Get zao index
              indx = Red%izao(ffjtran,ia,iz)

              ! Get indexes
              if0l2 = Red%zao(indx)%tgf0
              if1l2 = Red%zao(indx)%tgf1
              nf = if1l2 - if0l2 + 1

              ! Size for collection
              ! The 5 assumes that will be using splines to
              ! interpolate
              if (nf.gt.0) TRAMc = TRAMc + 8d-6*dble(nf*(ndir*2 + 5))

              ! If this CPU does not have frequencies
              ! in this line, skip
              if (Atom(ia)%fflag(jtran)%absent) cycle

              ! Get indexes
              if0l2 = Red%zao(indx)%gf0
              if1l2 = Red%zao(indx)%gf1
              nf = if1l2 - if0l2 + 1

              ! Skip if transition is outside limits
              if (nf.le.0) cycle

              ! Size for storage
              TRAMc = TRAMc + 8d-6*dble(ndir*2*nf)

            end do ! Transition level
          end do ! Transition term
        end do ! Atom
      end do ! Height

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

      integer:: ia,jtran,iz,indx,if0l2,if1l2,ndir

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
            indx = Red%izao(jtran,ia,iz)

            ! Get indexes
            if0l2 = Red%zao(indx)%gf0
            if1l2 = Red%zao(indx)%gf1

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

      integer:: ia,jtran,iz,indx,fjtran,ffjtran,if0l2,if1l2,ndir

      double precision:: vel


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
              indx = Red%izao(ffjtran,ia,iz)

              ! Get indexes
              if0l2 = Red%zao(indx)%gf0
              if1l2 = Red%zao(indx)%gf1

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

      logical:: lfield, lvel

      integer:: iz,ia,jtran,itermu,itermf,indx,t0,t1
      integer:: if0l2,if1l2,if0tl2,if1tl2,if0Il2,if1Il2,ifreq,nf,nfl
      integer:: ndir,idir,ierr,iph,ith,nth,nph

      double precision:: vel,DwT,Dw,ct,st,cc,sc,vfac
      double precision, dimension(:), allocatable:: omega
      double precision, dimension(:,:), allocatable:: splin
      double precision, dimension(:,:), allocatable:: eps0
      double precision, dimension(:,:), allocatable:: eps1
      double precision, dimension(:,:), allocatable:: eps2
      double precision, dimension(:,:), allocatable:: eps3
      double precision, dimension(:,:), allocatable:: eps20
      double precision, dimension(:,:), allocatable:: eps21
      double precision, dimension(:,:), allocatable:: eps22
      double precision, dimension(:,:), allocatable:: eps23

      ! Pointers
      type(Redb_class), pointer:: p_red
      type(Redb2_class), pointer:: p_rwarr
      type(Prof_class), pointer:: p_Norm
      complex(kind=8), dimension(:,:,:,:), pointer:: TKQo
      type(Redb2_class), target:: p_dummy


      ! Nullify pointers
      nullify(p_red,p_rwarr,p_Norm,TKQo)

      ! Dummy
      if (.not.PRAM) p_rwarr => p_dummy

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

      ! For each height
      do iz=Rz0,Rz1_PRD

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

        ! For each atom
        do ia=1,nA

          ! Thermal part of the Doppler width
          DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

          ! For each output transition
          do jtran=1,Atom(ia)%ntran

            ! Skip if no PRD
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! Get zao index
            indx = Red%izao(jtran,ia,iz)

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            itermf = Atom(ia)%fst(jtran)%iterml

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

            ! Get limits for this transition
            if0l2 = Red%zao(indx)%gf0
            if1l2 = Red%zao(indx)%gf1
            if0tl2 = Red%zao(indx)%tgf0
            if1tl2 = Red%zao(indx)%tgf1
            if0Il2 = Red%zao(indx)%Igf0
            if1Il2 = Red%zao(indx)%Igf1

            ! Transition ranges
            t0 = Atom(ia)%tshift + 1
            t1 = t0 + Atom(ia)%ntran - 1

            ! If no output, skip
            if ((if1tl2-if0tl2).lt.0) cycle

            ! Allocate emissivity to collect
            allocate(eps0(ndir,if0tl2:if1tl2))
            allocate(eps1(ndir,if0tl2:if1tl2))
            allocate(eps2(ndir,if0tl2:if1tl2))
            allocate(eps3(ndir,if0tl2:if1tl2))
            allocate(eps20(ndir,if0tl2:if1tl2))
            allocate(eps21(ndir,if0tl2:if1tl2))
            allocate(eps22(ndir,if0tl2:if1tl2))
            allocate(eps23(ndir,if0tl2:if1tl2))

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
                           Geom%njdir,p_Norm,Dw, &
                           eps0(:,if0Il2:if1Il2), &
                           eps1(:,if0Il2:if1Il2), &
                           eps2(:,if0Il2:if1Il2), &
                           eps3(:,if0Il2:if1Il2))

                ! Call second order emissivity
                call emiss2ord(Atom(ia),Geom,Atmo%vx(iz), &
                               Atmo%vy(iz),Atmo%vz(iz),lvel, &
                               Frec%omega,p_red,p_rwarr,p_Norm, &
                               Flgsg,jtran, &
                               itermu,itermf,iz,if0Il2,if1Il2, &
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
                             Geom%njdir,p_Norm,Dw, &
                             eps0(:,if0Il2:if1Il2), &
                             eps1(:,if0Il2:if1Il2), &
                             eps2(:,if0Il2:if1Il2), &
                             eps3(:,if0Il2:if1Il2))

                ! Call second order emissivity
                call emiss2ordNB(Atom(ia),Geom,Atmo%vx(iz), &
                                 Atmo%vy(iz),Atmo%vz(iz),lvel, &
                                 Frec%omega,p_red,p_rwarr,p_Norm, &
                                 Flgsg,jtran, &
                                 itermu,itermf,iz,if0Il2,if1Il2, &
                                 DwT,Dw,Atmo%vmi(iz), &
                                 TKQo,Stokes(:,:,:,:,iz), &
                                 JKQa,JKQ(:,:,:,iz),JKQC(:,:,:,iz), &
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

            end if ! Local size


            ! If MPI
            if (nproc.gt.1) then

              ! Share
              call MPI_ALLGATHERV(MPI_IN_PLACE,0, &
                                  MPI_DATATYPE_NULL, &
                                  eps20,Red%zao(indx)%nf*ndir, &
                                  ndir*(Red%zao(indx)%Mif0 - 1), &
                                  MPI_DOUBLE_PRECISION, &
                                  MPI_COMM_RT, ierr)
              call MPI_ALLGATHERV(MPI_IN_PLACE,0, &
                                  MPI_DATATYPE_NULL, &
                                  eps21,Red%zao(indx)%nf*ndir, &
                                  ndir*(Red%zao(indx)%Mif0 - 1), &
                                  MPI_DOUBLE_PRECISION, &
                                  MPI_COMM_RT, ierr)
              call MPI_ALLGATHERV(MPI_IN_PLACE,0, &
                                  MPI_DATATYPE_NULL, &
                                  eps22,Red%zao(indx)%nf*ndir, &
                                  ndir*(Red%zao(indx)%Mif0 - 1), &
                                  MPI_DOUBLE_PRECISION, &
                                  MPI_COMM_RT, ierr)
              call MPI_ALLGATHERV(MPI_IN_PLACE,0, &
                                  MPI_DATATYPE_NULL, &
                                  eps23,Red%zao(indx)%nf*ndir, &
                                  ndir*(Red%zao(indx)%Mif0 - 1), &
                                  MPI_DOUBLE_PRECISION, &
                                  MPI_COMM_RT, ierr)
            end if ! MPI


            !
            ! Get what is needed
            !

            ! If the CPU allocated eps20
            if (allocated(Red%zao(indx)%eps20)) then

              ! If dynamic, we need to interpolate
              if (lvel) then

                ! Number of frequencies
                nf = if1tl2 - if0tl2 + 1
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

                ! Just fetch
                Red%zao(indx)%eps20 = eps20(:,if0l2:if1l2)
                Red%zao(indx)%eps21 = eps21(:,if0l2:if1l2)
                Red%zao(indx)%eps22 = eps22(:,if0l2:if1l2)
                Red%zao(indx)%eps23 = eps23(:,if0l2:if1l2)

              end if ! Dynamic
            end if ! If need to take anything

            ! Deallocate emissivity to collect
            deallocate(eps0,eps1,eps2,eps3)
            deallocate(eps20,eps21,eps22,eps23)

          end do ! Output transition
        end do ! Atoms
      end do ! Heights

      ! Nullify pointers
      nullify(p_red,p_rwarr,p_Norm,TKQo)

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

      logical:: lvel

      integer:: iz,ia,jtran,itermu,itermf,indx,jdir,t0,t1
      integer:: if0l2,if1l2,if0tl2,if1tl2,if0Il2,if1Il2,ifreq,nf,nfl
      integer:: if0jl2,if1jl2,ffjtran,ffktran,ffltran,fjtran
      integer:: iJf,iJu,ndir,idir,ierr,iph,ith,nth,nph

      double precision:: vel,DwT,Dw,ct,st,cc,sc,vfac,Dfreq
      double precision, dimension(:), allocatable:: omega
      double precision, dimension(:), allocatable:: eps0
      double precision, dimension(:,:), allocatable:: splin
      double precision, dimension(:,:), allocatable:: eps20,rpf

      ! Pointers
      type(Redb_class), pointer:: p_red
      type(Redb2_class), pointer:: p_rwarr
      type(Prof_class), pointer:: p_Norm
      type(Redb2_class), target:: p_dummy


      ! Nullify pointers
      nullify(p_red,p_rwarr,p_Norm)

      ! Dummy
      if (.not.IRAM) p_rwarr => p_dummy

      ! Number of directions
      if (los) then
        nth = 1
        nph = 1
      else
        nth = Geom%nTh
        nph = Geom%nPh
      end if

      ! For each height
      do iz=Rz0,Rz1_PRD

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

        ! For each atom
        do ia=1,nA

          ! Thermal part of the Doppler width
          DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

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

              ! Get izao index
              indx = Red%izao(ffjtran,ia,iz)

              ! Get frequency of FS transition
              Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                      Atom(ia)%FSfreq(iJf,itermf)

              ! Add the microt. to Doppler width
              Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

              ! Real index in trano
              ffltran = Atom(ia)%itrano(ffjtran)

              ! Get limits for this transition
              if0l2 = Red%zao(indx)%gf0
              if1l2 = Red%zao(indx)%gf1
              if0tl2 = Red%zao(indx)%tgf0
              if1tl2 = Red%zao(indx)%tgf1
              if0Il2 = Red%zao(indx)%Igf0
              if1Il2 = Red%zao(indx)%Igf1
              if0jl2 = Red%zao(indx)%ggf0
              if1jl2 = Red%zao(indx)%ggf1

              ! Transition ranges
              t0 = Atom(ia)%tfshift + 1
              t1 = t0 + Atom(ia)%nftran - 1

              ! If no output, skip
              if ((if1tl2-if0tl2).lt.0) cycle

              ! Allocate emissivity to collect
              allocate(eps0(if0tl2:if1tl2))
              allocate(eps20(ndir,if0tl2:if1tl2))
              allocate(rpf(ndir,if0tl2:if1tl2))

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
                            itermu,itermf,iJu,iJf,iz,if0Il2,if1Il2, &
                            p_Norm,Dw,eps0(if0Il2:if1Il2))

                ! Second order
                call emissI2ord(Atom(ia),Geom,Atmo%vx(iz), &
                                Atmo%vy(iz),Atmo%vz(iz),lvel, &
                                Frec%omega,p_red,p_rwarr,p_Norm, &
                                ndir,jtran,fjtran,itermu,itermf, &
                                iJu,iJf,iz, &
                                if0Il2,if1Il2,DwT,Dw, &
                                Atmo%vmi(iz), &
                                Stokes(:,:,:,iz), &
                                JKQ(t0:t1,iz), &
                                JKQC(if0jl2:if1jl2,iz), &
                                eps20(:,if0Il2:if1Il2), &
                                rpf(:,if0Il2:if1Il2))

                ! Combine
                do idir=1,ndir
                  eps20(idir,if0Il2:if1Il2) = eps0(if0Il2:if1Il2) + &
                                             eps20(idir,if0Il2:if1Il2)
                end do

              end if ! Local calculation

              ! If MPI
              if (nproc.gt.1) then

                ! Share
                call MPI_ALLGATHERV(MPI_IN_PLACE,0, &
                                    MPI_DATATYPE_NULL, &
                                    eps20,Red%zao(indx)%nf*ndir, &
                                    ndir*(Red%zao(indx)%Mif0 - 1), &
                                    MPI_DOUBLE_PRECISION, &
                                    MPI_COMM_RT, ierr)
                call MPI_ALLGATHERV(MPI_IN_PLACE,0, &
                                    MPI_DATATYPE_NULL, &
                                    rpf,Red%zao(indx)%nf*ndir, &
                                    ndir*(Red%zao(indx)%Mif0 - 1), &
                                    MPI_DOUBLE_PRECISION, &
                                    MPI_COMM_RT, ierr)
              end if ! MPI

              !
              ! Get what is needed
              !

              ! If emissivity is allocated
              if (allocated(Red%zao(indx)%eps20)) then

                ! If dynamic, we need to interpolate
                if (lvel) then

                  ! Number of frequencies
                  nf = if1tl2 - if0tl2 + 1
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

                  ! Just fetch
                  Red%zao(indx)%eps20 = eps20(:,if0l2:if1l2)
                  Red%zao(indx)%rpf = rpf(:,if0l2:if1l2)

                end if ! Dynamic
              end if ! If need to take anything

              ! Deallocate emissivity to collect
              deallocate(eps0,eps20,rpf)

            end do ! Output transition FS
          end do ! Output transition term
        end do ! Atoms
      end do ! Heights

      ! Nullify pointers
      nullify(p_red,p_rwarr,p_Norm)

      end subroutine comoving_emissI2ord

!#####################################################################
!#####################################################################
!#####################################################################

      end module comovingprd_mod
